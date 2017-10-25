/************************************************************************
File mmps_prism_intf.c - implementation of parallel interface routines for
                        meshes of prismatic elements

Contains definitions of interface routines:
  mmpr_init_mesh - to initialize the parallel mesh data structure

  mmpr_el_owner - to get an owning processor identifier of an element
  mmpr_el_set_owner - to set an owning processor identifier of an element
  mmpr_fa_owner - to get an owning processor identifier of a face
  mmpr_fa_set_owner - to set an owning processor identifier of a face
  mmpr_ed_owner - to get an owning processor identifier of an edge
  mmpr_ed_set_owner - to set an owning processor identifier of an edge
  mmpr_ve_owner - to get an owning processor identifier of a vertex
  mmpr_ve_set_owner - to set an owning processor identifier of a vertex
  mmpr_el_id_at_owner - to get an local identifier of an element
  mmpr_el_set_id_at_owner - to set an local identifier of an element
  mmpr_fa_id_at_owner - to get an local identifier of a face
  mmpr_fa_set_id_at_owner - to set an local identifier of a face
  mmpr_ed_id_at_owner - to get an local identifier of an edge
  mmpr_ed_set_id_at_owner - to set an local identifier of an edge
  mmr_ve_id_at_owner - to get an local identifier of a vertex
  mmpr_ve_set_id_at_owner - to set an local identifier of a vertex

  mmpr_init_ref - to initialize the parallel process of refinement
  mmpr_refine - to refine an element or the whole mesh in parallel
  mmpr_derefine - to derefine an element or the whole mesh in parallel
  mmpr_final_ref - to finalize the process of parallel refinement
  mmpr_free_mesh - to free space allocated for parallel mesh data structure
------------------------------
History:
	09.2012 - Kazimierz Michalik, parmetis version
	2015 - Kazimierz Michalik, fpcm version
*************************************************************************/

#include "mmph_adapter.h"

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include <boost/foreach.hpp>
#include<algorithm>
#include<stack>
#include<set>

/* interface for the mesh manipulation module */
#include "mmh_intf.h"

#include "ddh_intf.h"

/* interface for the paralle mesh manipulation module */
#include "mmph_intf.h"

#include "uth_intf.h"
#include "pch_intf.h"
#include "uth_io_results.h"
#include "CompositeTransfererWithOwnership.hpp"
#include "CompositeTransferData.hpp"


std::vector<mmpt_mesh>  mmpv_meshes;        /* array of meshes */


// Cleaning functions for all meshes.
void mmpr_atexit_clear_meshes();


mmpt_mesh::mmpt_mesh()
: mesh_id(-1)
  {
        // Register cleaning function.
        atexit(mmpr_atexit_clear_meshes);

        assert(pcr_is_parallel_initialized() == 1);
        transferer = new CompositeTransfererWithOwnership(pcv_my_proc_id, *this);
        ownership.setLocalOwnerID(pcv_my_proc_id);
  }

mmpt_mesh::mmpt_mesh(const mmpt_mesh & other)
    : mesh_id(other.mesh_id)
{
    transferer =new CompositeTransfererWithOwnership(pcv_my_proc_id, *this);
    ownership = other.ownership;
    history.clear();
}

mmpt_mesh::~mmpt_mesh()
{
    clear();
}

void mmpt_mesh::init()
{
    int bufid =-1;
    if(pcv_my_proc_id == PCC_MASTER_PROC_ID) {
        const int n_nodes = mmr_get_nr_node(mesh_id);
        double* coords=new double[3*n_nodes];

        double* it = coords;
        int nno=0;
        while(0!= (nno=mmr_get_next_node_all(mesh_id,nno))) {
            mmr_node_coor(mesh_id,nno,it);
            it+=3;
        }

        glob_ids.createMeshBase(coords,n_nodes);

        delete [] coords; coords=NULL;

        // sending mesh base to all processes
        bufid = pcr_send_buffer_open(MMPC_MSG_ID_MESH_BASE,PCC_DEFAULT_BUFFER_SIZE);
        pcr_buffer_pack_double(MMPC_MSG_ID_MESH_BASE,bufid,3,glob_ids.origin);
        pcr_buffer_pack_double(MMPC_MSG_ID_MESH_BASE,bufid,3,glob_ids.d);
        pcr_buffer_pack_int(MMPC_MSG_ID_MESH_BASE,bufid,3,glob_ids.binary_oper);
        pcr_buffer_pack_char(MMPC_MSG_ID_MESH_BASE,bufid,1,(char*) &glob_ids.in_line_digits);
        pcr_buffer_pack_char(MMPC_MSG_ID_MESH_BASE,bufid,1,(char*) &glob_ids.line_digits);
        pcr_buffer_pack_char(MMPC_MSG_ID_MESH_BASE,bufid,1,(char*) &glob_ids.shift);
    }

    bufid = pcr_buffer_bcast(MMPC_MSG_ID_MESH_BASE,bufid,PCC_MASTER_PROC_ID);

    if(pcv_my_proc_id != PCC_MASTER_PROC_ID) {
        pcr_buffer_unpack_double(MMPC_MSG_ID_MESH_BASE,bufid,3,glob_ids.origin);
        pcr_buffer_unpack_double(MMPC_MSG_ID_MESH_BASE,bufid,3,glob_ids.d);
        pcr_buffer_unpack_int(MMPC_MSG_ID_MESH_BASE,bufid,3,glob_ids.binary_oper);
        pcr_buffer_unpack_char(MMPC_MSG_ID_MESH_BASE,bufid,1,(char*) &glob_ids.in_line_digits);
        pcr_buffer_unpack_char(MMPC_MSG_ID_MESH_BASE,bufid,1,(char*) &glob_ids.line_digits);
        pcr_buffer_unpack_char(MMPC_MSG_ID_MESH_BASE,bufid,1,(char*) &glob_ids.shift);
        pcr_recv_buffer_close(MMPC_MSG_ID_MESH_BASE,bufid);
    }
}

bool mmpt_mesh::check() const
{
    return ownership.check();
}

const ApplyChanges& mmpt_mesh::apply(const TransferOrder &o, const TransferResult &r)
{
    history.push_back(new ApplyChanges());
    ApplyChanges& ap = * *history.rbegin();
    // Return value:
    std::vector<int>&    vertices_old2newIds = ap.vertices_old2newIds;

    // Sender side
    if(o.Source_id == pcv_my_proc_id) {
        if(o.t_policy == TRANSFER_MOVE) {
            mmr_init_ref(mesh_id); // To make sure that mmr_del_elem will work.
            for(int i=0; i < o.Elements.size(); ++i) {
                mmr_del_elem(mesh_id,o.Elements[i]);
                ownership.RemoveID<MMC_ELEMENT>(o.Elements[i]);
            }

            mmr_is_ready_for_proj_dof_ref(mesh_id);
            mmr_final_ref(mesh_id);

            TransferResult::mmpt_faceSideMap::const_iterator it(r.new_sub_bdn_faces.begin());
            for(; it != r.new_sub_bdn_faces.end(); ++it) {
                switch(it->second) {
                case    FaceSide0:
                case    FaceSide1:
                        mmr_fa_set_sub_bnd(mesh_id,it->first,it->second);
                        break;
                case FaceSideBoth:
                    mmr_del_face(mesh_id,it->first);
                    ownership.RemoveID<MMC_FACE>(it->first);
                    break;

                }
            }

        }
//        else { // TRANSFER_COPY
//        }
    }
    else { // Reciever side

        mmr_init_read(mesh_id,
                    mmr_get_max_node_id(mesh_id)+r.new_vertices_coords.size()/3,
                    mmr_get_max_edge_id(mesh_id)+o.Elements.size()*MMC_MAXELEDGES,
                    mmr_get_max_face_id(mesh_id)+o.Elements.size()*MMC_MAXELFAC,
                    mmr_get_max_elem_id(mesh_id)+o.Elements.size());

        mmr_init_ref(mesh_id);

        // MODYFING THE MESH
        // Resizing mesh struct to hold new enitities.
        const std::vector<int> & vertices = r.new_vertices_transfered_ids;
        const std::vector<double>& vts_coords = r.new_vertices_coords;
        const int n_copied_elems = r.new_elems.size();
        mmr_init_read(mesh_id,
                    mmr_get_max_node_id(mesh_id)+vertices.size(),
                    mmr_get_max_edge_id(mesh_id)+n_copied_elems*MMC_MAXELEDGES,
                    mmr_get_max_face_id(mesh_id)+n_copied_elems*MMC_MAXELFAC,
                    mmr_get_max_elem_id(mesh_id)+n_copied_elems);


        vertices_old2newIds.resize( *std::max_element(vertices.begin(),vertices.end()) + 1, -1 );
        std::vector<int>::const_iterator v_it(vertices.begin()),v_end(vertices.end());
        std::vector<double>::const_iterator v_coor_it(vts_coords.begin());

        for(;v_it != v_end; ++ v_it) {
        EntOwn   vtx_ownership;
            int loc_id = EntOwnerships::UNDEFINED;

            // 1. Perhaps vertex owner is not from Source_id, and we have it under some other glob_id?
            std::vector<int>::const_iterator it =
                std::find( r.ownerships2vtx_trans_id.begin(),
                           r.ownerships2vtx_trans_id.end(),
                           *v_it);

            if(it != r.ownerships2vtx_trans_id.end()) {
                const int glob_id_pos = it - r.ownerships2vtx_trans_id.begin();
                vtx_ownership = r.transfer_vertices_ownerships[glob_id_pos];
                loc_id = ownership.GetLocIdForOwnership<MMC_NODE>(vtx_ownership);
            }

            // if vertex is not found we can add it
            if(loc_id == EntOwnerships::UNDEFINED) {
                vertices_old2newIds[*v_it] = mmr_add_node(mesh_id, MMC_AUTO_GENERATE_ID, &(*v_coor_it) );
                //mf_log_info("VE recieved as %d was added as %d.",
//                            *v_it,
//                            vertices_old2newIds[*v_it]);
            }
            else {
                vertices_old2newIds[*v_it] = loc_id;
                //mf_log_info("VE recieved as (%d,%d) is local %d.",
//                            vtx_glob_id.owner, vtx_glob_id.id_at_owner,
//                            loc_id);
            }
            v_coor_it += 3;
        }

        ap.element_pos2newIds.resize(r.new_elems.size(),0);
        for(int e=0; e < r.new_elems.size();  ++e) {
            static const int VERTS_PER_TYPE[8]={0,0,2,3,4,6,8,4};
            const int n_of_vts = VERTS_PER_TYPE[r.new_elems[e].type];
            int elem_new_verts_ids[MMC_MAXELVNO]={0};

            for(int v=0; v < n_of_vts; ++v) {
                mf_check_debug(vertices_old2newIds[r.new_elems[e].vertices[v]] != 0,
                               "Incorrect vertex number(%d)!",
                               vertices_old2newIds[r.new_elems[e].vertices[v]]);
                elem_new_verts_ids[v] = vertices_old2newIds[r.new_elems[e].vertices[v]];
                mf_check_debug(elem_new_verts_ids[v] != 0,
                               "Incorrect mapping from external vertex id(%d) to local id(%d)",
                               vertices_old2newIds[r.new_elems[e].vertices[v]],elem_new_verts_ids[v]);
            }

            ap.element_pos2newIds[e]
                    = mmr_add_elem(mesh_id,MMC_AUTO_GENERATE_ID,
                                 r.new_elems[e].type,
                                 elem_new_verts_ids, NULL,
                                 r.new_elems[e].groupID );

            int faces[MMC_MAXELFAC+1]={0};
            mmr_el_faces(mesh_id, ap.element_pos2newIds[e] ,faces,NULL);
            for(int f=1; f <= faces[0]; ++f) {
                mmr_fa_set_bc(mesh_id,faces[f],r.new_elems[e].BCs[f-1]);
            }

        }

        mmr_is_ready_for_proj_dof_ref(mesh_id);
        mmr_final_ref(mesh_id);

    }


    mmr_test_mesh(mesh_id);
    r.wasApplied = true;
    return ap;
}

const ApplyChanges& mmpt_mesh::apply(const TransferOrder &o, const TransferWithOwnershipResult &r)
{
    const ApplyChanges& ap = apply(o,static_cast<const TransferResult&>(r));

    //mf_debug("Entering apply(ownership)");
    if(o.t_policy == TRANSFER_COPY) {
        // const int id_offset = mmr_get_max_elem_id(mesh_id)+1;

        EntOwn glob_id(o.Source_id, EntOwnerships::UNDEFINED);
        for(int i=0; i < r.new_elems.size(); ++i) {
            glob_id.id_at_owner = r.new_elems[i].id_at_prev_owner;
            ownership.AddOwnership<MMC_ELEMENT>(ap.element_pos2newIds[i],glob_id);
            neighbours.insert( ownership.GetOwner(glob_id) );
            mf_log_info("apply(ownership): EL at %d is (%d,%d)",
                        ap.element_pos2newIds[i],
                        ownership.GetOwner(glob_id),
                        ownership.GetId_at_owner(glob_id));
        }

        for(int i=0; i < o.Elements.size(); ++i) {
            ownership.setBoundarable<MMC_ELEMENT>(o.Elements[i]);
        }

    }

  
    assert(r.transfer_vertices_ownerships.size() <= r.new_vertices_transfered_ids.size());
    for(int v=0; v < r.transfer_vertices_ownerships.size(); ++v) {
        const EntOwn& glob_id = r.transfer_vertices_ownerships[v];
        if(ownership.GetOwner(glob_id)!= pcv_my_proc_id) { // real foregin id
            if(EntOwnerships::UNDEFINED == ownership.GetLocIdForOwnership<MMC_NODE>(glob_id) ) {
                ownership.AddOwnership<MMC_NODE>(ap.vertices_old2newIds[ r.ownerships2vtx_trans_id[v] ],
                                         glob_id);
                neighbours.insert( ownership.GetOwner(glob_id));
                //mf_log_info("apply(ownership): VE at %d is (%d,%d)",
//                            ap.vertices_old2newIds[ r.glob_id2vtx_trans_id[v] ],
//                            glob_id.owner,
//                            glob_id.id_at_owner);
            }
            else {
                //mf_log_info("not_apply(ownership): VE at %d is (%d,%d)",
//                            ap.vertices_old2newIds[ r.glob_id2vtx_trans_id[v] ],
//                            glob_id.owner,
//                            glob_id.id_at_owner);
            }
        }
        else { // it's ours node!
            mf_check( ownership.GetId_at_owner(glob_id)==  ap.vertices_old2newIds[ r.ownerships2vtx_trans_id[v] ],
                    "Process internal node (%d) with wrong glob id(%d,%d)!",
                    ap.vertices_old2newIds[ r.ownerships2vtx_trans_id[v] ],
                    ownership.GetOwner(glob_id),
                    ownership.GetId_at_owner(glob_id));
        }
    }
    return ap;
}

const ApplyChanges& mmpt_mesh::apply(const TransferOrder &o,
                                                const CompositeTransferWithOwnershipResult &r)
{
    const ApplyChanges& ap_prev = apply(o,static_cast<const TransferWithOwnershipResult&>(r) );
    history.push_back(new ApplyChanges(ap_prev));
    ApplyChanges& ap = * *history.rbegin();

    if(! r.transfer_orders.empty() ) {
        // copying recived transfer directions
        ap.transfer_directions.insert(ap.transfer_directions.begin(),
                                      r.transfer_orders.begin(),r.transfer_orders.end());
        // rewriting relative element position as transfer element to correct local id
        for(int o=0; o < ap.transfer_directions.size(); ++o) {
            for(int e=0; e < ap.transfer_directions[o].Elements.size(); ++e) {
                int & el_id = ap.transfer_directions[o].Elements[e];    // id from previous process

                for(int i=0; i < r.new_elems.size(); ++i) {
                    if(r.new_elems[i].id_at_prev_owner == el_id) {
                        //mf_log_info("Rewriting element id to transfer from %d to %d.", el_id, ap.element_pos2newIds.at(i) );
                        el_id = ap.element_pos2newIds.at(i);    //
                        break;
                    }
                }
            }
        }
    }
    return ap;
}

void mmpt_mesh::findContestedExternalVerticesOfElements(const std::vector<int> & element_ids, std::vector<int> & externalVertices) const
{
    // NOTE: ONLY FOR VERTICES!!!

    externalVertices.clear();

    unordered_map<int,int> face2times;
    int faces[MMC_MAXELFAC+1]={0};

    //  check all faces belonging to elements to send:
    for(int i=0; i < element_ids.size(); ++i) {
        mmr_el_faces(mesh_id,element_ids[i],faces,NULL);
        for(int f=1; f <= faces[0]; ++f) {
            ++face2times[faces[f]];
        }
    }

    // if face is twice - it's internal
    // else it is subdomain or bc face

    unordered_map<int,int>::const_iterator it=face2times.begin();
    std::vector<int> newSubdomainBndFaces;
    for(; it != face2times.end(); ++it) {
        mf_check_debug(it->second >= 1, "Wrong number of face negihs(%d)!",it->second);
        mf_check_debug(it->second <= 2, "Wrong number of face negihs(%d)!",it->second);
        if(it->second == 1) {
            // check whether face is not outside (BC) face
            int neigs[2]={0};
            mmr_fa_neig(mesh_id, it->first, neigs, NULL,NULL,NULL,NULL,NULL);
            if(neigs[1] != MMC_BOUNDARY) { // if second neigbour is not at boundary
                newSubdomainBndFaces.push_back(it->first);
            }
        }
    }

    int     nodes[MMC_MAXFAVNO+1]={0};
    for(int f=0; f < newSubdomainBndFaces.size(); ++f) {
        mmr_fa_node_coor(mesh_id,newSubdomainBndFaces[f],nodes,NULL);

        for(int v=1; v <= nodes[0]; ++v) {
            externalVertices.push_back(nodes[v]);
        }
    }

    std::sort(externalVertices.begin(),externalVertices.end());
    externalVertices.erase(std::unique(externalVertices.begin(),externalVertices.end()),
                           externalVertices.end());

}

void mmpt_mesh::findContestedVerticesAssigmentToProcs(const std::vector<TransferOrder> & orders, // in:
                                          std::map<int,EntOwn>& contested_vrts_glob_ids) const  // out:
{
    contested_vrts_glob_ids.clear();
    std::vector<int>    contested_vertices;
    std::map<int,int>       procs_max_vtx_id;

    // Step 1. Find all vertices which stays with this process (current owner).
    // findContestedExternalVerticesOfElements( ALL orders elements superposition)
    std::vector<int>    all_elements;
    for(int o=0; o < orders.size(); ++o) {
        const TransferOrder& order = orders[o];
        if(order.Source_id == pcv_my_proc_id) {
            all_elements.insert(all_elements.begin(),
                                order.Elements.begin(),order.Elements.end());
        }
        procs_max_vtx_id[order.Destination_id] = 0; // Assuming empty procs!!!
    }
    std::sort(all_elements.begin(),all_elements.end());
    all_elements.erase( std::unique(all_elements.begin(), all_elements.end()),
                        all_elements.end()  );

    findContestedExternalVerticesOfElements(all_elements,contested_vertices);

    EntOwn globID_tmp(pcv_my_proc_id, EntOwnerships::UNDEFINED);
    for(int v=0; v < contested_vertices.size();++v) {
        ownership.SetId_at_owner(globID_tmp,contested_vertices[v]);
        contested_vrts_glob_ids[contested_vertices[v]] = globID_tmp;
    }
    contested_vertices.clear();
    // Step 2. Find all contested vertices and setup ownership
    // with greedy strategy (lowest takes ownership).

    for(int o=0; o < orders.size(); ++o) {
        const TransferOrder& order = orders[o];

        if(order.Source_id == pcv_my_proc_id) {

            findContestedExternalVerticesOfElements(order.Elements,contested_vertices);

            for(int i=0; i < contested_vertices.size(); ++i) {
                const int & vertexId = contested_vertices[i];
                if(contested_vrts_glob_ids.find(vertexId) == contested_vrts_glob_ids.end()) {
                    EntOwn globID_tmp(ownership.GetId_at_owner<MMC_NODE>(vertexId),
                                         ownership.GetOwner<MMC_NODE>(vertexId) );
                    if( ownership.GetOwner(globID_tmp) == pcv_my_proc_id ) { // other case was done in step 1.
                        ownership.Set(globID_tmp,order.Destination_id,EntOwnerships::UNDEFINED);
                    }
                    contested_vrts_glob_ids[vertexId] = globID_tmp;
//                    mf_log_info("Setting contested vtx glob_id %d,%d",
//                                globID_tmp.owner, globID_tmp.id_at_owner );

                }

                if(order.Destination_id < ownership.GetOwner(contested_vrts_glob_ids[vertexId])) {
                    // change ownership
//                    mf_log_info("findContestedVerticesAssigmentToProcs: Changing ownership (%d,%d) to (%d,%d)",
//                                contested_vrts_glob_ids[vertexId].owner,
//                                contested_vrts_glob_ids[vertexId].id_at_owner,
//                                order.Destination_id,
//                                EntitiesOwnership::UNDEFINED);

                    ownership.Set(contested_vrts_glob_ids[vertexId],order.Destination_id,EntOwnerships::UNDEFINED);
                }
            }
        }
    }
}

void mmpt_mesh::clear() {
    mesh_id = -1;
    for(int i=0; i < history.size(); ++i) {
        delete history[i];
    }
    history.clear();
    delete transferer;
    transferer = NULL;
}


template<>
void mmpt_mesh::GetOffspring<MMC_ELEMENT>(std::vector<EntOwn> & Parents,
                               std::vector<LID> & Offspring,
                               std::vector<GID> & Offspring_gids) const
{
    int sons[MMC_MAXELSONS+1]={0};

    ownership.getBoundarablesGlobIDs<MMC_ELEMENT>(Parents, EntOwnerships::INNER_BOUNDARABLE);

    Offspring.reserve(Parents.capacity()*MMC_MAXELSONS);

    for(int i=0; i < Parents.size(); ++i)  {
        EntOwn &parent =Parents[i];
        assert(ownership.GetOwner(parent) == pcv_my_proc_id);

        mmr_el_fam(mesh_id,ownership.GetId_at_owner(parent),sons, NULL);
        Offspring.insert(Offspring.end(),sons+1, sons+1+sons[0]);

        // we have to also get information about new vertices
        // created during element refinement
        // that is:  mid-face nodes (for QUAD faces)
        // and mid-edge nodes

        int faces[MMC_MAXELFAC+1]={0};
        mmr_el_faces(mesh_id,ownership.GetId_at_owner(parent), faces, NULL);
        for(int f=1; f <= faces[0]; ++f) {
            if(mmr_fa_type(mesh_id,faces[f]) == MMC_QUAD) {
                bool mid_node_ours=false;
                int neigs[2]={0};
                mmr_fa_neig(mesh_id,faces[f],neigs,NULL,NULL,NULL,NULL,NULL);

                if(neigs[1] == MMC_BOUNDARY) {
                    mid_node_ours=true;
                }
                else {
                    neigs[0] = ownership.GetOwner<MMC_ELEMENT>(neigs[0]);
                    neigs[1] = ownership.GetOwner<MMC_ELEMENT>(neigs[1]);

                    if(std::min(neigs[0],neigs[1]) == pcv_my_proc_id) {
                        mid_node_ours = true;
                    }
                }

                int mid_node=EntOwnerships::UNDEFINED;
                if(mid_node_ours) {
                    mmr_fa_fam(mesh_id,faces[f],NULL,&mid_node);
                }

                Offspring.push_back( mid_node );
                Offspring_gids.push_back( GlobID_Node(mid_node) );
            }
        }

        int edges[MMC_MAXELEDGES+1]={0};
        mmr_el_edges(mesh_id,parent.id_at_owner,edges);
        for(int e=1; e <= edges[0]; ++e) {
            int edge_mid_node = -1;
            mmr_edge_sons(mesh_id,edges[e],NULL,&edge_mid_node);
            mf_check_debug(edge_mid_node != -1, "Incorrect id for edge mid-node!");
            Offspring.push_back( edge_mid_node );
            Offspring_gids.push_back( GlobID_Node(edge_mid_node) );
        }

    }
}

template<>
void mmpt_mesh::SetOffspring<MMC_ELEMENT>(const   std::vector<EntOwn> & Parents,
                                          const   std::vector<int> & Offspring,
                                          const   std::vector<GID> & Offspring_gids)

{
    int offspring_loc_ids[MMC_MAXELSONS+1]={0};


    const int* it_offspring = & Offspring[0];

    for(int p=0; p < Parents.size(); ++p) {
        const int parent_loc_id = ownership.GetLocIdForOwnership<MMC_ELEMENT>(Parents[p]);

        mmr_el_fam(mesh_id,parent_loc_id,offspring_loc_ids,NULL);
        EntOwn tmp_glob_id(ownership.GetOwner(Parents[p]),EntOwnerships::UNDEFINED);
        for(int o=1; o <= offspring_loc_ids[0]; ++o, ++it_offspring) {
            ownership.SetId_at_owner(tmp_glob_id,*it_offspring);
            ownership.AddOwnership<MMC_ELEMENT>(offspring_loc_ids[o],tmp_glob_id);
        }

        // we have to also get information about new vertices
        // created during element refinement
        // that is:  mid-face nodes (for QUAD faces)
        // and mid-edge nodes
        int faces[MMC_MAXELFAC+1]={0};
        mmr_el_faces(mesh_id,parent_loc_id, faces, NULL);
        for(int f=1; f <= faces[0]; ++f) {
            if(mmr_fa_type(mesh_id,faces[f]) == MMC_QUAD) {
                if(*it_offspring != EntOwnerships::UNDEFINED) {
                    int mid_node=-1;
                    mmr_fa_fam(mesh_id,faces[f],NULL,&mid_node);
                    ownership.SetId_at_owner(tmp_glob_id,*it_offspring);
                    ownership.AddOwnership<MMC_NODE>(mid_node,tmp_glob_id);
                }

                ++it_offspring;
            }
        }

        int edges[MMC_MAXELEDGES+1]={0};
        mmr_el_edges(mesh_id,parent_loc_id,edges);
        for(int e=1; e <= edges[0]; ++e, ++it_offspring) {
            int edge_mid_node = -1;
            mmr_edge_sons(mesh_id,edges[e],NULL,&edge_mid_node);
            ownership.SetId_at_owner(tmp_glob_id,*it_offspring);
            ownership.AddOwnership<MMC_NODE>(edge_mid_node,tmp_glob_id);
        }
    }
}


/*------------------------------------------------------------
  mmpr_atexit_clear_meshes - to clear the parallel meshes data structures
------------------------------------------------------------*/
void mmpr_atexit_clear_meshes() {
    for(int i=0; i < mmpv_meshes.size(); ++i) {
        mmpv_meshes.at(i).clear();
    }
    mmpv_meshes.clear();
}


#ifdef __cplusplus
extern "C"{
#endif


int       mmpv_partition_method = MMPC_METIS_PART_METHOD;   /*  method that is used to generate partitions */

/*------------------------------------------------------------
  mmpr_init_mesh - to initialize the parallel mesh data structure
------------------------------------------------------------*/
int mmpr_init_mesh(  /* returns: >0 - Mesh ID, <0 - error code */
  int Control,
  int Mesh_id,
  int Nr_proc,
  int My_proc_id
  )
{
 utr_io_result_cumulative_timer_start(RESULT_TIME_MMPM);
/*++++++++++++++++ executable statements ++++++++++++++++*/
   int rc = Mesh_id;
  // mmpv_nr_sub = Nr_proc;
  mf_check(pcv_my_proc_id == My_proc_id, "Wrong process id!");

  mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

  assert(pmesh.mesh_id == Mesh_id);

//  pmesh.d_mesh.Clear();

  mf_check(pmesh.ownership.getLocalOwnerID() == pcv_my_proc_id, "Wrong local owner id");

  pmesh.init();

    if(pcv_nr_proc >= 2) {
       const int rc = mmpr_create_subdomains(Mesh_id,MMC_INIT_GEN_LEVEL);

      if(rc < 0) {
            mf_fatal_err("mmpr_create_subdomains failed.");
      }


    }
    else {
        mf_log_warn("Running parallel mesh in not parallel env: ignoring subdomains creation");
        rc = -1;
    }

    utr_io_result_cumulative_timer_stop(RESULT_TIME_MMPM);

    pmesh.check();

  return rc;
}



/*---------------------------------------------------------
  mmpr_create_subdomains - to decompose the mesh and create subdomains
---------------------------------------------------------*/
int mmpr_create_subdomains( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,      /* in: approximation field ID  */
  int Control        /* in: generation level as basis for decomposition */
  )
{
    assert(pcr_is_parallel_initialized() == 1);

#ifdef DEBUG_MMM
     printf("mmpr_create_subdomains: starting for Mesh_id %d, Control %d\n",Mesh_id,Control);
#endif

     int rc=0;

    mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    if(pcv_nr_proc < 2) {
            rc = -1;
            return rc;
    }

    if(pcr_is_this_master()) {

        std::vector<int>    el_ids,sdm_sizes;
        el_ids.resize(mmr_get_nr_elem(Mesh_id));
        sdm_sizes.resize(pcv_nr_proc);

        std::vector<int> overlap_sizes, overlap_elems;
        overlap_sizes.resize(pcv_nr_proc);
        overlap_elems.resize(pcv_nr_proc*mmr_get_nr_elem(Mesh_id)); // here was error with memory

        int * part  = NULL;

        rc = ddr_create_subdomains_scheme(Mesh_id, DDC_DEFAULT, pcv_nr_proc,
                              sdm_sizes.data(), el_ids.data(),
                              overlap_sizes.data(),  overlap_elems.data(), & part);


        if(rc != pcv_nr_proc) {
            mf_log_err("ddr_create_subdomains failed.");
            return rc=-2;
        }


        // send core elements and overlap to procs
        // mmpr_find_subdomain_boundary_entities(pmesh,part);

        std::vector<TransferOrder> transfer_orders(pcv_nr_proc+1);


        /// Core Elements transfer orders
        int* p_el_ids = el_ids.data() + sdm_sizes[0];   // move element pointer after MASTER_PROC elements
        //int* p_overlap_elems = overlap_elems.data();
        for(int p=PCC_MASTER_PROC_ID+1; p <= pcv_nr_proc; ++p) {

            TransferOrder el_order(pcv_my_proc_id,p,TRANSFER_MOVE,p_el_ids, sdm_sizes[p-1]);

            transfer_orders[p]= el_order;

            p_el_ids+=sdm_sizes[p-1];
        }

        /// Overlap Elements transfer directions
        int it_ovl_el=0;
         for(int destinationId =PCC_MASTER_PROC_ID; destinationId <= pcv_nr_proc; ++destinationId) {
            for(int el=0; el < overlap_sizes[destinationId-1]; ++el, ++it_ovl_el) {

                // if so - expand it;
                // part is numbered from 0, but overlap_elems from 1, so "-1"
                // source ID have to be numbered from 1, but part is from 0, so "+1"
                const int sourceId = part[overlap_elems[it_ovl_el]-1] + 1;
                assert(sourceId != destinationId);

                TransferOrder   *ovl_order=NULL;
                // find if order from part[overlap_elems[it_ovl_el] to destination already exist,
                for(int i=0; i < transfer_orders[sourceId].relative_transfer_orders.size() && ovl_order==NULL; ++i) {
                    if(transfer_orders[sourceId].relative_transfer_orders[i].Destination_id == destinationId) {
                        ovl_order = & transfer_orders[sourceId].relative_transfer_orders[i];
                    }
                }

                // no such order found, creating new one
                if(ovl_order == NULL) {

                   transfer_orders[sourceId].relative_transfer_orders.push_back(
                                        TransferOrder(sourceId,destinationId,TRANSFER_COPY));
                   ovl_order = & *transfer_orders[sourceId].relative_transfer_orders.rbegin();
                }

                ovl_order->Elements.push_back( overlap_elems[it_ovl_el] ); // +1 because ovl_el ID = ovl_el position + 1

                /// Also, we must note destination process to receive additional transfers.
                // But only once.
                std::vector<TransferOrder>::const_iterator it =
                        std::find(transfer_orders[destinationId].relative_transfer_orders.begin(),
                                  transfer_orders[destinationId].relative_transfer_orders.end(),
                                  TransferOrder(sourceId,destinationId,TRANSFER_COPY) );
                if( it ==  transfer_orders[destinationId].relative_transfer_orders.end() ) {
                    transfer_orders[destinationId].relative_transfer_orders.push_back(TransferOrder(sourceId,destinationId,TRANSFER_COPY));
                }
            }
         }





         const std::vector<TransferOrder> transfers_orders2 = transfer_orders[PCC_MASTER_PROC_ID].relative_transfer_orders;
         transfer_orders.erase(transfer_orders.begin(),transfer_orders.begin()+2);
         std::vector<const TransferResult*> transfer_results;

         pmesh.transferer->doTransfer(transfer_orders,transfer_results);

         for(int i=0; i < transfer_results.size(); ++i) {
             if(transfer_results[i] != NULL) {
                pmesh.apply(transfer_orders[i],
                            static_cast<const CompositeTransferWithOwnershipResult &>(*transfer_results[i]) );
            }
         }

        /// Transfering overlap from MASTER proc to other procs
        std::vector<const TransferResult*> transfer_results2;

        pmesh.transferer->doTransfer(transfers_orders2,transfer_results2);

        for(int i=0; i < transfers_orders2.size(); ++i ) {
            if( transfer_results2[i] != NULL ) {
                pmesh.apply(transfers_orders2[i],
                            static_cast<const TransferWithOwnershipResult &>(*transfer_results2[i]) );
            }
        }

        transfer_orders.clear();
        transfer_results.clear();
        transfer_results2.clear();

    }
    else { // other procs than master
        // recieve core elements and overlap
        TransferOrder recv_el_order(PCC_MASTER_PROC_ID,pcv_my_proc_id,TRANSFER_MOVE );
        const CompositeTransferWithOwnershipResult& r_el = static_cast<const CompositeTransferWithOwnershipResult &>(pmesh.transferer->doTransfer(recv_el_order));
        const ApplyChanges& apply_result = pmesh.apply(recv_el_order, r_el);

        for(int i=0; i < apply_result.transfer_directions.size(); ++i) {
                const TransferWithOwnershipResult& r
                = static_cast<const TransferWithOwnershipResult& >(pmesh.transferer->doTransfer(apply_result.transfer_directions[i]));
                pmesh.apply(apply_result.transfer_directions[i],r);
        }

    }

    return rc;
}


/*------------------------------------------------------------
  mmpr_balance_load - to balance load between processes using ParMetis
------------------------------------------------------------*/
int mmpr_balance_load(
  int Mesh_id  /* in: mesh ID */
  )
{
    utr_io_result_cumulative_timer_start(RESULT_TIME_MMPM);
    // This method is needed only if more than 1 process is working.
    int rc=METIS_OK;
    if(pcr_nr_proc() > 1) {
        assert(pcr_is_parallel_initialized() == 1);


#ifdef DEBUG_MMM
     printf("mmpr_balance_load: starting for Mesh_id %d\n",Mesh_id);
#endif


    // ddr_balance_subdomains(Mesh_id,);

    }
    utr_io_result_cumulative_timer_stop(RESULT_TIME_MMPM);
    return rc;
}

/**-----------------------------------------------------------
mmpr_update_ref_list - to update list of refined elements due to irregularity
                      constraints using inter-processor communication
------------------------------------------------------------*/
int mmpr_update_ref_list(
                   /** returns: >=0 - success code, <0 - error code */
  const int Mesh_id,     /** in: mesh ID */
  int *Nr_ref,     /** in/out: number of refined elements */
  int **List_ref   /** in/out: list of refined elements */
  )
  {
    if(pcv_nr_proc < 2) {
            int rc = -1;
            return rc;
    }

    int& nr_ref = *Nr_ref;
    int* &list_ref = *List_ref;

    std::sort(list_ref, list_ref+nr_ref);
    nr_ref = std::unique(list_ref, list_ref+nr_ref) - list_ref;

    mf_check(*Nr_ref == nr_ref, "Found duplicates %d in ref_list!", *Nr_ref - nr_ref);

    for(int i=0; i < nr_ref; ++i) {
        mf_log_info("list_ref: %d", list_ref[i]);
    }

    mmpt_mesh& pmesh = *mmpr_select_mesh(Mesh_id);

    // Check if any element concerns us at all.

    std::vector<const EntOwn*> elem_to_handle;

    int n_bound_el = pmesh.ownership.getBoundaryOwnerships<MMC_ELEMENT>(list_ref,list_ref+nr_ref, elem_to_handle);

    std::vector<EntOwn> adapt_glob_ids;
    if(n_bound_el > 0) {
        adapt_glob_ids.reserve(n_bound_el);

        int i=0;
        while(adapt_glob_ids.size() < n_bound_el) {
            if(elem_to_handle[i] != &EntOwnerships::OWN_NULL) {
                if(elem_to_handle[i] != &EntOwnerships::OWN_SELF) {
                adapt_glob_ids.push_back(*elem_to_handle[i]);
                    mf_log_info("updateing ref list with (%d,%d)",
                                elem_to_handle[i]->owner, elem_to_handle[i]->id_at_owner);
                }
                else {
                    adapt_glob_ids.push_back(EntOwn(pcv_my_proc_id,list_ref[i]));
                    mf_log_info("updateing ref list with (%d,%d)",
                               pcv_my_proc_id, list_ref[i]);
                }
            }
            ++i;
        }
    }

    for(mmpt_mesh::NEIG_SET::const_iterator it = pmesh.neighbours.begin();
        it != pmesh.neighbours.end(); ++it)
    {

        const int bufid = pcr_send_buffer_open(MMPC_ADAPT_REF_MSG_ID,PCC_DEFAULT_BUFFER_SIZE);
        pcr_buffer_pack_int(MMPC_ADAPT_REF_MSG_ID,bufid,1,&n_bound_el);
        if(n_bound_el > 0) {
            pcr_buffer_pack_int(MMPC_ADAPT_REF_MSG_ID,bufid,
                                2*n_bound_el,reinterpret_cast<const int*>(adapt_glob_ids.data()));
        }

        pcr_buffer_send(MMPC_ADAPT_REF_MSG_ID,bufid,*it);
    }

    mf_check(nr_ref <= mmr_get_nr_elem(Mesh_id), "More adaptables then elements!");

    ////// RECIEVER SIDE
    std::set<int> new_loc_ids;
    for(mmpt_mesh::NEIG_SET::const_iterator it = pmesh.neighbours.begin();
        it != pmesh.neighbours.end(); ++it)
    {
        int recieved=0;
        const int bufid = pcr_buffer_receive(MMPC_ADAPT_REF_MSG_ID,PCC_ANY_PROC,PCC_DEFAULT_BUFFER_SIZE);
        pcr_buffer_unpack_int(MMPC_ADAPT_REF_MSG_ID,bufid,1,& recieved);

        mf_log_info("Recieving %d new elemts to refine", recieved);

        if(recieved > 0) {
            adapt_glob_ids.resize(recieved);
            pcr_buffer_unpack_int(MMPC_ADAPT_REF_MSG_ID,bufid,
                                  2*recieved,reinterpret_cast<int*>(adapt_glob_ids.data()));

            for(int i=0; i < adapt_glob_ids.size(); ++i) {
                int loc_id = pmesh.ownership.GetLocIdForOwnership<MMC_ELEMENT>(adapt_glob_ids[i]);
                 mf_log_info("(%d,%d) recognize at pos %d",
                             pmesh.ownership.GetOwner(adapt_glob_ids[i]),
                             pmesh.ownership.GetId_at_owner(adapt_glob_ids[i]),
                             loc_id);
                if(loc_id != EntOwnerships::UNDEFINED) {
                    // so there is something that we have...
                    new_loc_ids.insert(loc_id);
                    pmesh.ownership.setBoundarable<MMC_ELEMENT>(loc_id);
                }
            }

        }
    }

    if(new_loc_ids.size() > 0) {
        // create new (bigger) array
        int new_nr_ref = nr_ref + new_loc_ids.size();
        int* new_ref_list = (int*) calloc(new_nr_ref, sizeof(int));
        // combine the two ranges....
        // NOTE: list_ref already sorted!
        int* new_ref_list_end = std::set_union(new_loc_ids.begin(), new_loc_ids.end(),
                                               list_ref, list_ref+nr_ref,
                                               new_ref_list);

        new_nr_ref = new_ref_list_end - new_ref_list;

        // free old data
        free(list_ref);
        // switch to new data
        nr_ref = new_nr_ref;
        list_ref = new_ref_list;
        new_ref_list = NULL;

        }
    mf_check(nr_ref <= mmr_get_nr_elem(Mesh_id), "More adaptables then elements!");

      return 1;
  }


/**-----------------------------------------------------------
  mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready
                                  for dofs projection
------------------------------------------------------------*/
int mmpr_is_ready_for_proj_dof_ref(
  int Mesh_id
)
{

    if(mmr_is_ready_for_proj_dof_ref(Mesh_id) < 0) {
        mf_fatal_err("Mesh not ready for dofs projection!");
    }

    if(pcv_nr_proc < 2) {
            int rc = -1;
            return rc;
    }

    mmpt_mesh& pmesh = *mmpr_select_mesh(Mesh_id);

    // Synchronize ID of newly created element and vertices.
    std::vector<EntOwn> el_parents_ownerships;
    std::vector<LID>  el_offspring_loc_ids;
    std::vector<GID>    el_offspring_gids;

    pmesh.GetOffspring<MMC_ELEMENT>(el_parents_ownerships, el_offspring_loc_ids,el_offspring_gids);
    mf_log_info("mmpr_is_ready_for_proj_dof_ref: got offspring ids (%d;%d)",
                el_parents_ownerships.size(),el_offspring_loc_ids.size());

    // Send to neighbours
    mmpt_mesh::NEIG_SET::const_iterator it=pmesh.neighbours.begin();
    for(; it != pmesh.neighbours.end(); ++it) {

        const int bufid = pcr_send_buffer_open(MMPC_ADAPT_PROJ_READY_MSG_ID,PCC_DEFAULT_BUFFER_SIZE);
        int n_par = el_parents_ownerships.size();
        pcr_buffer_pack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid,1,&n_par);
        if(n_par > 0) {
            pcr_buffer_pack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid,2*el_parents_ownerships.size(),
                                reinterpret_cast<const int*>(el_parents_ownerships.data()));
            int n_offspring = el_offspring_loc_ids.size();
            pcr_buffer_pack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid,1,&n_offspring);
            pcr_buffer_pack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid,n_offspring,
                                    el_offspring_loc_ids.data());
            pcr_buffer_pack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid,n_offspring,
                                el_offspring_gids.data());
        }
        pcr_buffer_send(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid,*it);

    }


    // Recieving
    const int n_neig = pmesh.neighbours.size();
    for(int p=0; p < n_neig; ++p) {
        const int bufid = pcr_buffer_receive(MMPC_ADAPT_PROJ_READY_MSG_ID,PCC_ANY_PROC,PCC_DEFAULT_BUFFER_SIZE);
        int n_par=0;
        pcr_buffer_unpack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid, 1, &n_par);
        if(n_par > 0) {
            el_parents_ownerships.resize(n_par);
            pcr_buffer_unpack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid, 2*n_par, reinterpret_cast<int*>(& el_parents_ownerships[0] ));
            int n_offspring=0;
            pcr_buffer_unpack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid, 1, &n_offspring);
            el_offspring_loc_ids.resize(n_offspring);
            pcr_buffer_unpack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid, n_offspring,& el_offspring_loc_ids[0]);
            el_offspring_gids.resize(n_offspring);
            pcr_buffer_unpack_int(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid, n_offspring, el_offspring_gids.data());
        }
        pcr_recv_buffer_close(MMPC_ADAPT_PROJ_READY_MSG_ID,bufid);

        pmesh.SetOffspring<MMC_ELEMENT>(el_parents_ownerships,el_offspring_loc_ids, el_offspring_gids);
        mf_log_info("mmpr_is_ready_for_proj_dof_ref: recieved offspring ids (%d;%d)",
                    el_parents_ownerships.size(),el_offspring_loc_ids.size());
    }


    return 1;
}


/*---------------------------------------------------------
  mmpr_select_mesh - to select the proper mesh
---------------------------------------------------------*/
mmpt_mesh* mmpr_select_mesh( /* returns pointer to the chosen mesh */
			   /* to avoid errors if input is not valid */
			   /* it returns the pointer to the current mesh */
  int Mesh_id    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{

    assert(Mesh_id > 0);

    if(static_cast<unsigned int>(Mesh_id) >= mmpv_meshes.size()) {
        mmpv_meshes.resize(Mesh_id+1);
        mmpv_meshes[Mesh_id].mesh_id=Mesh_id;
        mmpv_meshes[Mesh_id].ownership.setMeshId(Mesh_id);

    }
    assert(static_cast<unsigned int>(Mesh_id) <= mmpv_meshes.size());
    assert(Mesh_id == mmpv_meshes[Mesh_id].mesh_id);
    return & mmpv_meshes[Mesh_id];
}


/*--------------------------------------------------------------------------
  mmpr_el_owner - to get an owning processor identifier of an element
---------------------------------------------------------------------------*/
int mmpr_el_owner(/* returns: >0 -an owning processor identifier of an element */
	         /*          <0 - error code */
  int Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El         /* in: global (within the subdomain) element ID */
  )
{
    assert(Mesh_id > 0);
    assert(El > 0);
    assert(static_cast<unsigned int>(El) <= mmr_get_max_elem_id(Mesh_id));

    const mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    return pmesh.ownership.GetOwner<MMC_ELEMENT>(El);
}

/*--------------------------------------------------------------------------
  mmpr_fa_owner - to get an owning processor identifier of a face
---------------------------------------------------------------------------*/
int mmpr_fa_owner(  /* returns: >0 - an owning processor identifier of a face */
	          /*          <0 - error code */
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa        /* in: global (within the subdomain)  ID */
  )
{
    assert(Mesh_id > 0);
    assert(Fa > 0);
    assert(static_cast<unsigned int>(Fa) <= mmr_get_max_face_id(Mesh_id));

    const mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    return pmesh.ownership.GetOwner<MMC_FACE>(Fa);
}

/*--------------------------------------------------------------------------
  mmpr_ed_owner - to get an owning processor identifier of an edge
---------------------------------------------------------------------------*/
int mmpr_ed_owner(  /* returns: >0 an owning processor identifier of an edge */
	          /*          <0 - error code */
  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed          /* in: global (within the subdomain) edge ID */
  )
{
    assert(Mesh_id > 0);
    assert(Ed > 0);
    assert(static_cast<unsigned int>(Ed) <= mmr_get_max_elem_id(Mesh_id));

    const mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    return pmesh.ownership.GetOwner<MMC_EDGE>(Ed);
}


/*--------------------------------------------------------------------------
  mmpr_ve_owner - to get an owning processor identifier of a vertex
---------------------------------------------------------------------------*/
int mmpr_ve_owner( /* returns: >0 - an owning processor identifier of a vertex */
	         /*          <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve        /* in: global (within the subdomain) vertex ID */
  )
{
    assert(Mesh_id > 0);
    assert(Ve > 0);
    assert(static_cast<unsigned int>(Ve) <= mmr_get_max_node_id(Mesh_id));
    const mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    return pmesh.ownership.GetOwner<MMC_NODE>(Ve);
}


/*--------------------------------------------------------------------------
  mmpr_el_id_at_owner - to get an local identifier of an element
---------------------------------------------------------------------------*/
int mmpr_el_id_at_owner( /* returns: >0 - an local identifier of an element */
	         /*          <0 - error code */
  int Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El         /* in: global (within the subdomain) element ID */
  )
{
    assert(Mesh_id > 0);
    assert(El > 0);
    assert(static_cast<unsigned int>(El) <= mmr_get_max_elem_id(Mesh_id));

    const mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    return pmesh.ownership.GetId_at_owner<MMC_ELEMENT>(El);
}


///*--------------------------------------------------------------------------
//  mmpr_fa_id_at_owner - to get an local identifier of a face
//---------------------------------------------------------------------------*/
//int mmpr_fa_id_at_owner(  /* returns: >0 - an local identifier of a face */
//	          /*          <0 - error code */
//  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
//  int Fa        /* in: global (within the subdomain)  ID */
//  )
//{
//    assert(Mesh_id > 0);
//    assert(Fa > 0);
//    assert(static_cast<unsigned int>(Fa) < mmr_get_max_face_id(Mesh_id));
//
//    const mmpt_owners_map & omap = mmpr_select_mesh(Mesh_id)->ovl_faces;
//    const mmpt_owners_map::const_iterator it = omap.find(Fa);
//    return  it == omap.end() ? Fa : it->second.id_at_owner;
//}


///*--------------------------------------------------------------------------
//  mmpr_ed_id_at_owner - to get an local identifier of an edge
//---------------------------------------------------------------------------*/
//int mmpr_ed_id_at_owner(  /* returns: >0 an local identifier of an edge */
//	          /*          <0 - error code */
//  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
//  int Ed          /* in: global (within the subdomain) edge ID */
//  )
//{
//    assert(Mesh_id > 0);
//    assert(Ed > 0);
//    assert(static_cast<unsigned int>(Ed) < mmr_get_max_edge_id(Mesh_id));
//
//    const mmpt_owners_map & omap = mmpr_select_mesh(Mesh_id)->ovl_edges;
//    const mmpt_owners_map::const_iterator it = omap.find(Ed);
//    return  it == omap.end() ? Ed : it->second.id_at_owner;
//}


/*--------------------------------------------------------------------------
  mmr_ve_id_at_owner - to get an local identifier of a vertex
---------------------------------------------------------------------------*/
int mmpr_ve_id_at_owner( /* returns: >0 - an local identifier of a vertex */
	         /*          <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve        /* in: global (within the subdomain) vertex ID */
  )
{
    assert(Mesh_id > 0);
    assert(Ve > 0);
    assert(static_cast<unsigned int>(Ve) < mmr_get_max_edge_id(Mesh_id));

    const mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    return pmesh.ownership.GetId_at_owner<MMC_NODE>(Ve);
}


/*------------------------------------------------------------
  mmpr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
int mmpr_init_ref(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{

/* local variables */
  mmpt_mesh&  pmesh=*mmpr_select_mesh(Mesh_id);

///*++++++++++++++++ executable statements ++++++++++++++++*/


  return(0);
}


/*------------------------------------------------------------
  mmpr_final_ref - to finalize the process of parallel refinement
------------------------------------------------------------*/
int mmpr_final_ref(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{



  return(0);
}

/**-----------------------------------------------------------
  mmpr_n_neighbouring_procs - to get number of neighburing subdomains assigned to other procs.
  ///
  /// Domain 1 (D1) is neighbouring to the domain 2 (D2)
  /// only when they have at least one common mesh entity.
  /// <=> (D1 n D2) != {0}
  ///

------------------------------------------------------------*/
int mmpr_n_neighbouring_procs( // returns number of procs connected with caller proc
        const int Mesh_id //in: mesh id
        )
{
    return mmpr_select_mesh(Mesh_id)->neighbours.size();
}

/**-----------------------------------------------------------
  mmpr_interproc_connectivity - to get detailed info about neighbouring procs
------------------------------------------------------------*/
int mmpr_neighbouring_procs( // returns number of procs connected with caller proc
        const int Mesh_id,  //in: mesh id
        int * Neigh_proc_ids //out:[mmpr_ipc_n_connected_procs]
                                   // connected negibouring procs ids
        )
{
    const mmpt_mesh& pmesh = *mmpr_select_mesh(Mesh_id);

    mmpt_mesh::NEIG_SET::const_iterator it = pmesh.neighbours.begin();
    for(int i=0;it != pmesh.neighbours.end();++it) {
        Neigh_proc_ids[i++] = *it;
    }
    return pmesh.neighbours.size();
}


/**-----------------------------------------------------------
  mmpr_check_mesh - to check the parallel mesh data structure
------------------------------------------------------------*/
int mmpr_check_mesh(  /** returns: >0 - Mesh ID, <0 - error code */
  int Mesh_id
  )
{
    return mmpr_select_mesh(Mesh_id)->check() ? 1 : 0;
}


#ifdef __cplusplus
}
#endif
