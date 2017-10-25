#include <mpi.h>
#include <algorithm>

#include "Transferer.hpp"
#include "mmph_mesh.h"
#include "mmh_intf.h"
#include "pch_intf.h"
#include "uth_log.h"

namespace mmpt
{



bool    operator==(const TransferOrder& first, const TransferOrder& second)
{
    return first.Source_id == second.Source_id
            && first.Destination_id == second.Destination_id
            && first.t_policy == second.t_policy
            && first.Elements == second.Elements;
}


TransferOrder::TransferOrder(const int Sourceproc_id, const int Destproc_id,
                                 const TRANSFER_POLICY Policy,
                             const int *elements, const int n_elements)
                       : Source_id(Sourceproc_id),
                         Destination_id(Destproc_id),
                         t_policy(Policy),
                         Elements(elements, elements+n_elements)
{
}

Transferer::Transferer(int my_id, const mmpt_mesh &parallel_mesh) :
my_subdomain_id(my_id), pmesh(parallel_mesh), last_applaied_transfer(-1),
  current_orders(NULL)
{

}

Transferer::~Transferer()
{
    for(int i=0; i < transfers_results.size(); ++i) {
        delete transfers_results[i];
    }
    transfers_results.clear();
}

void Transferer::doMassTransferInit(const std::vector<TransferOrder>& orders)
{
    current_orders = & orders;
}

void Transferer::doMassTransferClear(const std::vector<TransferOrder> &orders, std::vector<const TransferResult*> &results)
{
        current_orders = NULL;
}

void Transferer::doTransfer(const std::vector<TransferOrder> &orders, std::vector<const TransferResult*> &results)
{

    // update transfer orders to each other
    doMassTransferInit(orders);

    // realize transfers
    for(int i=0; i < orders.size(); ++i) {
        if( !(orders[i].Elements.empty() && orders[i].Source_id==pcv_my_proc_id) ) {
            results.push_back(&doTransfer(orders[i]));
        }
        else {
            results.push_back(NULL);
        }
    }
    doMassTransferClear(orders,results);
}


const TransferResult &Transferer::doTransfer(const TransferOrder &order)
{
    assert(order.Source_id > 0 || order.Source_id==PCC_ANY_PROC);
    assert(order.Source_id <= pcr_nr_proc() || order.Source_id==PCC_ANY_PROC);
    assert(order.Destination_id > 0);
    assert(order.Destination_id <= pcr_nr_proc());
    assert(order.Destination_id != order.Source_id);

    if( !transfers_results.empty() ) {
        if( transfers_results.back()->wasApplied ) {
            last_applaied_transfer = transfers_results.size();
        }
    }

    transfers_history.push_back(order);
    TransferResult& r = * new TransferResult();
    transfers_results.push_back(&r);

    if(pcr_my_proc_id() == order.Source_id) {

        doSend(order,r);

    }
    else if(pcr_my_proc_id() == order.Destination_id) {

        doRecv(order,r);
    }

    return r;
}

void Transferer::doSend(const TransferOrder &order, TransferResult &r) const
{
    mf_log_info("Realizing transfer of element %d to %d", order.Elements[0], order.Destination_id);

    using std::vector;

    const int N_transfer_elems = order.Elements.size();
    int n_transfered_mesh_entities=0;
    int pcr_rc=0;
    int msg=0;

    assert(order.Elements.data() != NULL);
    assert(N_transfer_elems > 0);

//    mf_check(pmesh.sub_bnd_vtx_flags.size() == mmr_get_max_node_id(pmesh.mesh_id)+1, "Subdomain vertices map has incorrect size %d (should be %d). Data are not initialized?",
//             pmesh.sub_bnd_vtx_flags.size(),
//             mmr_get_max_node_id(pmesh.mesh_id)+1);
//        assert(pmesh.boundary_faces.size() >= 4);
//        assert(pmesh.boundary_faces.size() < mmr_get_nr_face(pmesh.mesh_id));


    const int buffer = pcr_send_buffer_open(MMPC_MSG_ID_COPYING,PCC_DEFAULT_BUFFER_SIZE);

    // sending N_copied_elems
    // >>>>>>>>>>>>> Send 1.
    //pcr_rc = pcr_send_int(Dest_proc_id,0,1,& N_copied_elems);
    pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,1,&N_transfer_elems);

    //mf_debug("ddr_copy_full_elems: 1 sended %d",pcr_rc);

    if(N_transfer_elems > 0) {
        // collecting what will be sended
        vector<int> elem_verts,elem_bcs, elem_groups;
        vector<int> & vertices = r.new_vertices_transfered_ids;
        elem_groups.reserve(N_transfer_elems);
        vertices.reserve(N_transfer_elems*4); // vertices used by elements to send
        elem_verts.reserve( vertices.capacity() );
        elem_bcs.reserve(N_transfer_elems*MMC_MAXELFAC);
        int nodes[MMC_MAXELVNO+1]={0};
        for(int e=0; e < N_transfer_elems; ++e) {
            int el_id = order.Elements[e];
            //mf_debug("Coping el %d",el_id);
            elem_groups.push_back(mmr_el_groupID(pmesh.mesh_id,el_id));
            mmr_el_node_coor(pmesh.mesh_id,el_id,nodes,NULL);
            assert(nodes[0]>0);
            elem_verts.insert( elem_verts.end(), nodes, nodes+1+nodes[0] ); // no of nodes + nodes ids
            vertices.insert( vertices.end(), nodes+1, nodes+1+nodes[0] );  // only nodes ids
            int faces[MMC_MAXELFAC+1]={0};
            mmr_el_faces(pmesh.mesh_id,el_id,faces,NULL);
            elem_bcs.push_back(faces[0]);

            // !!! MODYFING THE MESH
            for(int f=1; f <= faces[0]; ++f) {
                elem_bcs.push_back( mmr_fa_bc(pmesh.mesh_id,faces[f]) );

                if(order.t_policy == TRANSFER_MOVE) {
//                        mmr_init_ref(pmesh.mesh_id); // To make sure that mmr_del_elem will work.
                    // set face neig to MMC_SUB_BND
                    int neig[2]={0};
                    mmr_fa_neig(pmesh.mesh_id,faces[f],neig,NULL,NULL,NULL,NULL,NULL);

                    FaceSubdomainSide side_affected = FaceSideNone;


                    if(neig[0]==el_id) {
                        side_affected == FaceSide0;
                    }
                    else if(neig[1]==el_id){
                        side_affected == FaceSide1;
                    }
                    else {
                        mf_log_err("Wrong face neigs!");
                    }

                    if(side_affected != FaceSideNone) {
                        TransferResult::mmpt_faceSideMap::iterator it
                                = r.new_sub_bdn_faces.find(faces[f]);

                        if(it == r.new_sub_bdn_faces.end()) {
                            r.new_sub_bdn_faces[faces[f]] = side_affected;
                        }
                        else {

                        }

                    }

                }

            }
        }



        // remove duplicate entries
        std::sort(vertices.begin(),vertices.end());
        vertices.resize( std::distance( vertices.begin(), std::unique(vertices.begin(),vertices.end()) ) );


//          This really shouldn't be done in the transferer! At least in the mmpd_..
//          TODO: move to mmpd_
//        // sort-out boundary vertices at the end of collection
//        std::sort(vertices.begin(),vertices.end(), [&](int a, int b){
//            mf_log_info("Checking %d and %d", a,b);
//            return pmesh.sub_bnd_vtx_flags.at(a)< pmesh.sub_bnd_vtx_flags.at(b);
//        }
//        );
//
//        vector<char> vtx_subdomain_flags;
//        for( const int& ve : vertices) {
//            vtx_subdomain_flags.push_back(pmesh.sub_bnd_vtx_flags[ve]);
//        }

        // sending number of vertices
        int vertices_size = vertices.size();
        // >>>>>>>>>>>>> Send 2.
        //pcr_rc=pcr_send_int(Dest_proc_id,0,1, & vertices_size );
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,1, & vertices_size );
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 2 sended %d",pcr_rc);
        // sending vertices identyfiers
        // >>>>>>>>>>>>> Send 3.
        //pcr_rc=pcr_send_int(Dest_proc_id,0, vertices.size(), vertices.data() );
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,vertices.size(), vertices.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 3 sended %d",pcr_rc);

        vector<double> vts_coords;
        vts_coords.resize(vertices.size()*3);

        vector<int>::const_iterator v_it(vertices.begin()), v_end(vertices.end() );
        vector<double>::iterator v_coor_it(vts_coords.begin());
        for(;v_it != v_end; ++v_it) {
            mmr_node_coor(pmesh.mesh_id,*v_it, & (*v_coor_it) );
            v_coor_it += 3;
        }

        // sending verts coordinates
        // >>>>>>>>>>>>> Send 4.
        //pcr_rc=pcr_send_double(Dest_proc_id,0,vts_coords.size(), vts_coords.data() );
        pcr_rc = pcr_buffer_pack_double(MMPC_MSG_ID_COPYING,buffer,vts_coords.size(), vts_coords.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 4 sended %d",pcr_rc);

        // sending el_ids_at_owner
        // >>>>>>>>>>>>> Send 5.
        //pcr_send_int(Source_proc_id,0,N_copied_elems,Copied_elem_ids);
        pcr_rc=pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,N_transfer_elems,order.Elements.data());
        //mf_debug("ddr_transfer_full_elems: 5 sended %d",pcr_rc);

        // copying BC info
        // >>>>>>>>>>>>> Send 6.
        int elem_bcs_size = elem_bcs.size();
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,1,& elem_bcs_size);
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 8 sended %d",pcr_rc);
        // >>>>>>>>>>>>> Send 7.
        //pcr_rc=pcr_send_int(Dest_proc_id,0,elem_verts.size(), elem_verts.data() );
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,elem_bcs.size(), elem_bcs.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems:  9 sended %d",pcr_rc);

        // sending number of elem vertices
        int elem_verts_size = elem_verts.size();
        // >>>>>>>>>>>>> Send 8.
        //pcr_rc=pcr_send_int(Dest_proc_id,0,1, & elem_verts_size );
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,1,& elem_verts_size);
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 6 sended %d",pcr_rc);
        // sending elem vertices
        // >>>>>>>>>>>>> Send 9.
        //pcr_rc=pcr_send_int(Dest_proc_id,0,elem_verts.size(), elem_verts.data() );
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,elem_verts.size(), elem_verts.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 7 sended %d",pcr_rc);

        // elem_groups
        pcr_rc = pcr_buffer_pack_int(MMPC_MSG_ID_COPYING,buffer,elem_groups.size(), elem_groups.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: elem groups sended %d",pcr_rc);

        // >>>>>>>>>>>>> Send 10.

        n_transfered_mesh_entities = N_transfer_elems + vertices.size();

        pcr_rc = pcr_buffer_send(MMPC_MSG_ID_COPYING,buffer,order.Destination_id);
        mf_check(pcr_rc == MPI_SUCCESS, "ddr_transfer_full_elems: sending buffer!");


    }

}

void Transferer::doRecv(const TransferOrder &order, TransferResult &r) const
{
    using std::vector;

    const int N_transfer_elems = order.Elements.size();
    int n_transfered_mesh_entities=0;
    int pcr_rc=0;
    int msg=0;

    // recieving data
    // recieving n_copied_elems
    const int buffer = pcr_buffer_receive(MMPC_MSG_ID_COPYING,order.Source_id,PCC_DEFAULT_BUFFER_SIZE);
    int source_proc_id = order.Source_id;
    if(source_proc_id == PCC_ANY_PROC) {
        source_proc_id = pcr_buffer_source_id(buffer);
    }
    int n_copied_elems=0;
    // >>>>>>>>>>>>> Recv 1.
    //pcr_rc = pcr_receive_int(Source_proc_id,0,1,& n_copied_elems);
    pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,1,& n_copied_elems);
    assert(pcr_rc >= 0);
    //mf_debug("ddr_transfer_full_elems: 1 recived %d",pcr_rc);

    if(n_copied_elems > 0) {
        // recieving number of vertices
        int n_vertices=0;
        // >>>>>>>>>>>>> Recv 2.
        //pcr_rc=pcr_receive_int(Source_proc_id,0,1,& n_vertices);
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,1,& n_vertices);
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 2 recived %d",pcr_rc);
        // recieving vertices
        vector<int>& vertices = r.new_vertices_transfered_ids;
        vertices.resize(n_vertices,0);
        // >>>>>>>>>>>>> Recv 3.
        //pcr_rc=pcr_receive_int(Source_proc_id,0,n_vertices, vertices.data() );
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,n_vertices, vertices.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 3 recived %d",pcr_rc);

        // recieving vertices coords
        vector<double>& vts_coords = r.new_vertices_coords;
        vts_coords.resize(n_vertices*3, 0.0);
        // >>>>>>>>>>>>> Recv 4.
        //pcr_rc=pcr_receive_double(Source_proc_id,0,vts_coords.size(),vts_coords.data());
        pcr_rc = pcr_buffer_unpack_double(MMPC_MSG_ID_COPYING,buffer,vts_coords.size(),vts_coords.data());
        assert(pcr_rc == MPI_SUCCESS);
        //mf_debug("ddr_transfer_full_elems: 4 recived %d",pcr_rc);

        // recieving elem ids at owner
        vector<int> elem_ids_at_owner;
        elem_ids_at_owner.resize(n_copied_elems,0);
        // >>>>>>>>>>>>> Recv 5.
        //pcr_receive_int(Source_proc_id,0,n_copied_elems, elem_ids_at_owner.data() );
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,n_copied_elems, elem_ids_at_owner.data());
        //mf_debug("ddr_transfer_full_elems: 5 recived %d",pcr_rc);

        // recv BC info!
        int n_el_bcs=0;
        // >>>>>>>>>>>>> Recv 6.
        //pcr_rc=pcr_receive_int(Source_proc_id,0,1,& n_el_verts);
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,1,& n_el_bcs);
        assert(pcr_rc >= 0);
        assert( n_el_bcs > 0);
        //mf_debug("ddr_transfer_full_elems: 8 recived %d",pcr_rc);

        vector<int>  elem_bcs(n_el_bcs,0);
        // >>>>>>>>>>>>> Recv 7.
        //pcr_rc=pcr_receive_int(Source_proc_id,0,n_el_verts, elem_verts.data() );
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,n_el_bcs, elem_bcs.data() );
        assert(pcr_rc >= 0);
        //mf_debug("ddr_transfer_full_elems: 9 recived %d",pcr_rc);

        // recieving number of elems verts
        int n_el_verts=0;
        // >>>>>>>>>>>>> Recv 8.
        //pcr_rc=pcr_receive_int(Source_proc_id,0,1,& n_el_verts);
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,1,& n_el_verts);
        assert(pcr_rc >= 0);
        assert( n_el_verts > 0);
        //mf_debug("ddr_transfer_full_elems: 6 recived %d",pcr_rc);

        vector<int>  elem_verts(n_el_verts,0);
        // >>>>>>>>>>>>> Recv 9.
        //pcr_rc=pcr_receive_int(Source_proc_id,0,n_el_verts, elem_verts.data() );
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,n_el_verts, elem_verts.data() );
        assert(pcr_rc >= 0);
        //mf_debug("ddr_transfer_full_elems: 7 recived %d",pcr_rc);

        // elem groups
        std::vector<int>    elem_groups(n_copied_elems,0);
        pcr_rc = pcr_buffer_unpack_int(MMPC_MSG_ID_COPYING,buffer,n_copied_elems, elem_groups.data() );
        assert(pcr_rc >= 0);
        //mf_debug("ddr_transfer_full_elems: 7 recived %d",pcr_rc);

//        // >>>>>>>>>>>>> Recv 10.


        for(int i=0; i < n_vertices; ++i) {
            mf_check_debug(source_proc_id > 0,"Source proc number invalid (%d)", source_proc_id);
            mf_check_debug(source_proc_id != pcv_my_proc_id,"Source proc number invalid (%d)", source_proc_id);
        }

        int n_elems_new=0;
        vector<int>::const_iterator bc_it= elem_bcs.begin();
        vector<int>::const_iterator groups_it= elem_groups.begin();
        vector<int>::const_iterator it = elem_verts.begin(), end = elem_verts.end(), it_elem_ids_at_owner(elem_ids_at_owner.begin());
        for(; it != end; ++it_elem_ids_at_owner, ++groups_it) {
            const int n_of_vts = *it;
            ++it;
            mf_check_debug(*it != 0,"Wrong vertex id");

            const int el_type_from_n_vts[]={0,0,0,0, MMC_TETRA,0, MMC_PRISM};

            ++n_elems_new;

            TransferElement tel(el_type_from_n_vts[n_of_vts],
                                &(*it),
                                & *(bc_it+1),
                                *it_elem_ids_at_owner,
                                *groups_it); //group number

             it+=n_of_vts;

            r.new_elems.push_back(tel);

            bc_it += el_type_from_n_vts[n_of_vts]==MMC_TETRA ? 5 : 6;

        }

        mfp_check_debug(n_copied_elems == n_elems_new, "No. of re-created elements does not match no. of recived elems.!");


        n_transfered_mesh_entities = n_copied_elems + vertices.size();

    }
    pcr_recv_buffer_close(MMPC_MSG_ID_COPYING,buffer);
}


}// namespace
