/************************************************************************
File dds_metis_adapter_intf.c - implementation of interface routines for 
								Metis and ParMetis - mesh partitioning libraries 

Contains definitions of interface routines: 
	
  mar_partition_mesh - to to decompose the mesh and create subdomains by Metis library
  mar_adaptive_repartition - to balance the work load of a graph that corresponds to an adaptively refined mesh
  mar_refine_partition - to improve the quality of an existing a k-way partitioning
  mar_get_parts - get partition vector for the elements of the mesh. Indexes - active element number. values - procesors ids
  mar_get_parts_with_overlap - get partition vectors for the elements of the mesh with overlap info

  mar_set_metis_options - to specify metis behaviour
  mar_set_parmetis_options - to specify parmetis behaviour
  mar_initialize_work - to initialize the internal data structures
  mar_end_work - free memory for internal data structures
  mar_mesh_to_graph - fill lookup_table. Collapse elements into nodes.
  mar_mesh_to_CSR - fill xadj and adjncy data structures (in CSR format)
  mar_mesh_to_distributed_CSR - fill xadj and adjncy data structures (in Distributed CSR format)
  mar_gather_part_local - gather local partitions and merge them to global
  mar_scatter_part_local - satter local partitions from proc 0 to all other
  mar_check_mpi_flags - check is parallel enviroment variables
  mar_minimalize_communication - reorder newly created partion vector to maximalize match to old partiton vector
  mar_set_proc_ids_to_objects - distribute partitions ids to mesh objects
  mar_create_overlap - create overlap, fill data structure, set internal
  mar_find_el_conn_by_nodes - recurent procedure that find neighbours of org_el by his node ids
  mar_set_subdomains_boundaries - to find and set domain boundaries (faces set)

------------------------------  			
History:
    10.2012 - Kazimierz Michali, incorporating into ModFem
	08.2012 - Kamil Wachala, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include<algorithm>
#include<cassert>
#include<metis.h>
#include<parmetis.h>
#include<numeric>
#include<vector>
#include<deque>
//#include<boost/foreach.hpp>
#include <sstream>
#include<iostream>

/* internal header file for the metis adapter module */
#include "mmh_intf.h"
#include "uth_intf.h"
#include "pch_intf.h"
#include "ddh_metis_adapter.h"
#include "../mmph_adapter.h"
#include "uth_log.h" //#include "dbg.h"

#ifdef __cplusplus
extern "C"{
#endif


/* GLOBAL VARIABLES */
//mmpt_data pmesh.data;


// /* Metis */
// idx_t options[METIS_NOPTIONS]; /* allow to fine-tune and modify various aspects of the internal algorithms used by METIS */
// idx_t* xadj = nullptr; /* input CSR format. Indicates on amount of adjected verticles*/
// idx_t* adjncy = nullptr; /* input CSR format. Indicates on adjected verticles */
// idx_t* adjwgt = nullptr; /* weights of the edges */

// /* Parmetis */
// idx_t* vtxdist = nullptr; /* describes how the vertices of the graph are distributed among the processors */
// idx_t** part_local = nullptr; /* upon successful completion stores the local partition vector of the mesh (ids for each element) */
// idx_t** xadj_local = nullptr; /* local input CSR format. Indicates on amount of adjected verticles*/
// idx_t** adjncy_local = nullptr; /* local input CSR format. Indicates on adjected verticles */
// idx_t** adjwgt_local = nullptr; /* local weights of the edges */
// real_t* tpwgts = nullptr; /* fraction of vertex weight that should be distributed to each sub-domain for each balance constraint */
// real_t* ubvec = nullptr; /* specify the imbalance tolerance for each vertex weight, with 1 being perfect balance and nparts being perfect imbalance. A value of 1.05 recommended */
// idx_t wgt_flag = 0,*wgtflag = &wgt_flag; /* indicate if the graph is weighted. 0 - No weights, 1 - weights on edges */
// idx_t num_flag = 0, *numflag = &num_flag; /* numbering scheme that is used for the vtxdist, xadj, adjncy, and part. 0 C-style numbering */
// real_t i_tr = 0, *itr = &i_tr; /* ratio of inter-processor communication time compared to data redistribution time. Value of 1000.0 is recommended*/

// MPI_Comm comm_ = MPI_COMM_WORLD, *comm = &comm_; /* pointer to the MPI communicator of the processes that call PARMETIS */

// /* Metis and Parmetis */
// idx_t* part = nullptr; /* upon successful completion stores the partition vector of the mesh (ids for each element) */
// idx_t n_parts = 0, *nparts = &n_parts; /* The number of parts to partition the mesh */
// idx_t obj_val = 0, *objval = &obj_val; /* after complition: Stores the edge-cut or the total communication volume of the partitioning solution  */
// idx_t n_con = 0, *ncon = &n_con; /* the number of balancing constraints */
// idx_t met_parmet_result = 0; /* the value returned by metis or parmetis routines */

// /* Internal */
// // replaced with mmpv_mmpv_my_proc_id  //int mmpv_my_proc_id = 0; /* currently processor id that execude code */
// // replaced with mmpv_nr_sub  //int mmpv_nr_sub = 0; /* total amount of anavabile processors */





///*------------------------------------------------------------
//  mmpr_recreate_owner_tables - to initialize the internal data structures
//------------------------------------------------------------*/
//void mmpr_recreate_owner_tables(int Mesh_id)
//{
//    mmpt_mesh& pmesh = *mmpr_select_mesh(Mesh_id);
//    pmesh.elems.resize(mmr_get_max_elem_max(Mesh_id)+1);
//    pmesh.faces.resize(mmr_get_max_face_max(Mesh_id)+1);
//    pmesh.edges.resize(mmr_get_max_edge_max(Mesh_id)+1);
//    pmesh.nodes.resize(mmr_get_max_node_max(Mesh_id)+1);

//  /* allocate space for nodes' structures */
//    assert(pcr_my_proc_id()!=1 || !pmesh.elems.empty()); // there is no mesh without elements... (at least for master proc)

//   // at least master proc must have mesh
//    assert(pcr_my_proc_id()!=1 || !pmesh.elems.empty());
//    assert(pcr_my_proc_id()!=1 || !pmesh.faces.empty());
//    assert(pcr_my_proc_id()!=1 || !pmesh.edges.empty());
//    assert(pcr_my_proc_id()!=1 || !pmesh.nodes.empty());

//    // Setting up initial ownership info.

//    int my_proc_id = pcr_my_proc_id();
//    for (int Nel=0; (Nel=mmr_get_next_elem_all(Mesh_id,Nel));) {
//        //int rc = mmpr_el_set_owner(Mesh_id,Nel,my_proc_id);
//        //rc = mmpr_el_set_id_at_owner(Mesh_id,Nel,Nel);
//        pmesh.elems[Nel].owner = MMPC_MY_OWNERSHIP;
//        pmesh.elems[Nel].id_at_owner = Nel;
//    }
//    for (int Nfa=0; (Nfa=mmr_get_next_face_all(Mesh_id,Nfa));) {
////        int rc = mmpr_fa_set_owner(Mesh_id,Nfa,my_proc_id);
////        assert(rc >= 0);
////        rc = mmpr_fa_set_id_at_owner(Mesh_id,Nfa,Nfa);
////        assert(rc >= 0);
//        pmesh.faces[Nfa].owner = MMPC_MY_OWNERSHIP;
//        pmesh.faces[Nfa].id_at_owner = Nfa;

//    }
//    for (int Ned=0; (Ned=mmr_get_next_edge_all(Mesh_id,Ned));) {
////        int rc = mmpr_ed_set_owner(Mesh_id,Ned,my_proc_id);
////        assert(rc >= 0);
////        rc = mmpr_ed_set_id_at_owner(Mesh_id,Ned,Ned);
////        assert(rc >= 0);
//        pmesh.edges[Ned].owner = MMPC_MY_OWNERSHIP;
//        pmesh.edges[Ned].id_at_owner = Ned;
//    }
//    for (int Nno=0; (Nno=mmr_get_next_node_all(Mesh_id,Nno));) {
////        int rc = mmpr_ve_set_owner(Mesh_id,Nno,my_proc_id);
////        assert(rc >= 0);
////        rc = mmpr_ve_set_id_at_owner(Mesh_id,Nno,Nno);
////        assert(rc >= 0);
//        pmesh.nodes[Nno].owner = MMPC_MY_OWNERSHIP;
//        pmesh.nodes[Nno].id_at_owner = Nno;
//    }
//}


/* ROUTINES */

void mmpr_print_ownership(const int Mesh_id) {
    int n=mmr_get_nr_node(Mesh_id);
    int i=0;
    std::cout << "No. of nodes: " << n << "\n";
    while((i=mmr_get_next_node_all(Mesh_id,i)) != 0) {
        std::cout << "Node " << i << " owner: " << mmpr_ve_owner(Mesh_id,i)
                  << " id_at_owner: " << mmpr_ve_id_at_owner(Mesh_id,i) << "\n";
    }

    n=mmr_get_nr_elem(Mesh_id);
    i=0;
    std::cout << "No. of elems: " << n << "\n";
    while((i=mmr_get_next_elem_all(Mesh_id,i)) != 0) {
        std::cout << "Elem " << i << " owner: " << mmpr_el_owner(Mesh_id,i)
                  << " id_at_owner: " << mmpr_el_id_at_owner(Mesh_id,i) << "\n";
    }

}

 // Filling VtxDist table.
int mmpr_PCSR_update_vtxdist(const int Mesh_id,idx_t * Vtx_dist_table)
{
    assert(Vtx_dist_table != nullptr);
    mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);
    std::vector<int> vtx_counts(pcv_nr_proc,0);

    pmesh.data.ne = mmr_get_nr_elem(Mesh_id);

    pcr_allgather_int( & pmesh.data.ne , 1, vtx_counts.data(), 1 );
    std::fill_n(Vtx_dist_table,pcv_nr_proc+1,0);
    // because VtxDist[0]=0
    std::partial_sum(vtx_counts.begin(),vtx_counts.end(),Vtx_dist_table+1);
    assert(Vtx_dist_table[0]==0);
    return 0;
}

int mmpr_PCSR_distribute_elems_to_procs(const int Mesh_id,
                                  const int Source_proc_id)
{

    //ddr_print_ownership(Mesh_id);

    mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);
    ddt_PCSR &Pcsr = pmesh.data.mp.PCSR;
    int n_distributed_elems=0;
    std::vector<int> new_glob_pos;

    if(pcr_my_proc_id() == Source_proc_id) {
        mf_debug("ddr_distribute_elems_to_procs: in sedner proc");
        assert(Pcsr.empty() == false);


        // mark boundary faces
        pmesh.bnd_vtx_flags.resize(mmr_get_nr_node(Mesh_id)+1,0);
        int fa=0;
        while( 0 != (fa=mmr_get_next_face_all(Mesh_id,fa)) ) {
            int neighs[2]={0},neigsSubdomains[2]={0};
            mmr_fa_neig(Mesh_id,fa,neighs,NULL,NULL,NULL,NULL,NULL);


            if(neighs[1] != 0) {//ignore external Boundary Conditions
                if(Pcsr.part_.at(neighs[0]-1) != Pcsr.part_.at(neighs[1]-1)) {
                    pmesh.boundary_faces.push_back(fa);
                    int nodes[5]={0};
                    mmr_fa_node_coor(Mesh_id,fa,nodes,NULL);
                    for(int i=1; i <= nodes[0];++i) {
                        if(pmesh.bnd_vtx_flags[nodes[i]]!=1) {
                            pmesh.bnd_vtx_flags[nodes[i]]=1;
                            mfp_debug("Vertex %d is at boundary.",nodes[i]);
                        }
                   }
                }
            }
        }

        pmesh.boundary_faces.erase(std::unique(pmesh.boundary_faces.begin(),pmesh.boundary_faces.end()), pmesh.boundary_faces.end());
        mfp_debug("\nNo. of boundary faces=%d", pmesh.boundary_faces.size());

        // Splitting initial PCSR into per-proc-PCSRs
        // Waring: element indexes in per-proc-PCSRs will be different!
        typedef ddt_CSR::CsrNodeInternal cNode;
        std::vector<ddt_PCSR> procs_pcsr(pcr_nr_proc());
        std::for_each(Pcsr.begin(),Pcsr.end(), [&](cNode & n) {
            // split into subvectors
            procs_pcsr[*n.ppart].push_back(n);
        });



        new_glob_pos.resize(Pcsr.size(),0);
        int new_pos=Pcsr.vtxDist()[pcv_my_proc_id-1];
        for(int p=0; p < procs_pcsr.size();++p) {
            // Creating mapping array: new_glob_ids[old_id]=new_id
            if(procs_pcsr[p].size() > 0) {
                std::for_each(procs_pcsr[p].begin(),procs_pcsr[p].end(),
                              [&](const cNode & n) {
                                new_glob_pos[(*n.pel_id)-1]=new_pos;
                                ++new_pos;
                              });
            }

        }

        for(int p=0; p < procs_pcsr.size();++p) {
            // Establish correct identyfiers in PCSR adjncy tables.
            for(idx_t & neigh_pos : procs_pcsr[p].adjncy_) {
                neigh_pos = new_glob_pos[neigh_pos]; // +1 for id, not pos
            }

            // Transfering elements (using pre-created lists)
            if(p+1 != Source_proc_id) {
                mmpr_transfer_full_elems(Mesh_id,Source_proc_id,p+1,
                                        procs_pcsr[p].size(), procs_pcsr[p].el_loc_id_.data(),MMPE_MOVE );
                n_distributed_elems+=procs_pcsr[p].size();

                procs_pcsr[p].resize_n_proc(pcv_nr_proc);

                std::vector<idx_t> buf;
                procs_pcsr[p].write(buf);

                //mf_debug("ss: %s",ss.str().data());
                const int n_nums = buf.size();
                mf_debug("ddr_distribute_elems_to_procs: sending archived PCSR with size %d",n_nums);
                const int buffer = pcr_send_buffer_open(MMPC_DISTIRB_MSG_ID,PCC_DEFAULT_BUFFER_SIZE);

                int rc = pcr_buffer_pack_int(MMPC_DISTIRB_MSG_ID,buffer,1,&n_nums);
                assert(rc == MPI_SUCCESS);
                rc = pcr_buffer_pack_int(MMPC_DISTIRB_MSG_ID,buffer,n_nums, buf.data());
                assert(rc == MPI_SUCCESS);

                rc = pcr_buffer_send(MMPC_DISTIRB_MSG_ID,buffer,p+1);
                assert(rc == MPI_SUCCESS);
            }
        }
        assert(new_pos == Pcsr.size());

        // Removing transfered elements from this proc structures.
        Pcsr = std::move(procs_pcsr[pcv_my_proc_id-1]);
    }
    else {
        mf_debug("ddr_distribute_elems_to_procs: in recv proc");
        mf_check_debug(Source_proc_id < pcv_my_proc_id,"Wrong source proc number in distrib (%d)", Source_proc_id);
        mf_check_debug(Source_proc_id > 0,"Wrong source proc number in distrib (%d)", Source_proc_id);
        mmpr_transfer_full_elems( Mesh_id, Source_proc_id, pcr_my_proc_id(), 0, nullptr, MMPE_MOVE );

        int n_nums=0;
        const int buffer = pcr_buffer_receive(MMPC_DISTIRB_MSG_ID, Source_proc_id , PCC_DEFAULT_BUFFER_SIZE);
        int rc = pcr_buffer_unpack_int(MMPC_DISTIRB_MSG_ID, buffer,1,&n_nums);
        //int rc=pcr_receive_int(Source_proc_id,MMPC_DISTIRB_MSG_ID,1,&n_bytes);
        assert(rc == MPI_SUCCESS);
        mf_debug("ddr_distribute_elems_to_procs: reciving PCSR archive size: %d",n_nums);
        assert(n_nums > 0);

        std::vector<idx_t> v;
        v.resize(n_nums, 0);
        rc = pcr_buffer_unpack_int(MMPC_DISTIRB_MSG_ID,buffer,n_nums, const_cast<idx_t*>(v.data()));
        assert(rc == MPI_SUCCESS);

        Pcsr.clear();
        Pcsr.read(v);
        assert( !Pcsr.empty());

        // now in this proc graph nodes represeting elements have new ids
        for(int e=0; e < Pcsr.size(); ++e) {
            Pcsr.el_loc_id_[e]=e+1;
        }

        pcr_recv_buffer_close(MMPC_DISTIRB_MSG_ID,buffer);
    }

    mmr_test_mesh(Mesh_id);

    //ddr_print_ownership(Mesh_id);

    return n_distributed_elems;
}





/*------------------------------------------------------------
  ddr_create_inter_subdomain_connectivity - to initially create ipc
------------------------------------------------------------*/
int mmpr_PCSR_create_IPC(// return >=0 - ok, <0 - error
                         const int Mesh_id // IN: mesh id
                         )
{
    mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);
    int N_elems = pmesh.data.ne;
    assert(N_elems >= 0);
    const idx_t* Xadj = pmesh.data.mp.PCSR.xadj();
    //const idx_t* Adjwgt = pmesh.data.mp.PCSR.adjwgt();
    const idx_t* Adjncy = pmesh.data.mp.PCSR.adjncy();
    const idx_t* Adj_neig = pmesh.data.mp.PCSR.adj_neig_no();
    const idx_t* Vwgt = pmesh.data.mp.PCSR.vwgt();
    const idx_t* VtxDist = pmesh.data.mp.PCSR.vtxDist();
    const idx_t* El_id = pmesh.data.mp.PCSR.el_id();
    mmpt_data::ipc_type & ipc = pmesh.data.ipc;
    mmpt_data::ipc_links_type & ipc_links = pmesh.data.ipc_links;

    int retval=0;

    // Iterate throu adjency table, and cout neighbours from other subdomains.
    const idx_t*  it_xadj=Xadj;
    const idx_t*  it_neig_el=Adjncy;
    const idx_t*  it_vwgt=Vwgt;
    const idx_t* found_in_vtxdist=nullptr; // offset used to identity
    const idx_t* end_vtxdist = VtxDist+pcr_nr_proc()+1;
    const idx_t* min_el_id=VtxDist+pcr_my_proc_id()-1;
    const idx_t* max_el_id=VtxDist+pcr_my_proc_id();

    // Additional auxiliary table.
    std::vector<int> elem_seq;
    elem_seq.reserve(N_elems);
    int nel=0;
    while((nel=mmr_get_next_act_elem(Mesh_id,nel)) != 0) {
        elem_seq.push_back(nel);
    }

   // std::cout << pmesh.data.mp.PCSR;

    assert(elem_seq.size() == N_elems);

    for(int e=0; e < N_elems; ++e, ++it_xadj, ++it_vwgt) {
        for(const idx_t* it_adjency_proc_end = & Adjncy[*(it_xadj+1)];
        it_neig_el != it_adjency_proc_end; ++ it_neig_el) {
            // Compare neigbouring element id with VtxDist,
            // to find wich proc owns this element.
            // Proc with smaller id?
            found_in_vtxdist = nullptr;


            if((*it_neig_el) < (*min_el_id) ) {
                found_in_vtxdist = std::upper_bound( VtxDist+1, min_el_id, (*it_neig_el) ) -1;
                assert(found_in_vtxdist  < min_el_id);
                assert(found_in_vtxdist  >= VtxDist);
            } // or proc with greater id?
            else if((*it_neig_el) > (*max_el_id)) {
                found_in_vtxdist = std::upper_bound( max_el_id, end_vtxdist, (*it_neig_el) ) -1;
                assert(found_in_vtxdist  >= max_el_id);
                assert(found_in_vtxdist  < end_vtxdist);
            }

            if(found_in_vtxdist != nullptr) {
                // Compute owning_proc_id from relative distance from begining of VtxDist.
                const int owning_proc_id = std::distance(VtxDist, found_in_vtxdist )+1; // +1 to make id not rank
                assert(owning_proc_id > 0);

                mf_debug("Neigh %d owner %d",(*it_neig_el)+1,owning_proc_id);

                // Update inter process connectivity (including weight)
                ipc[owning_proc_id]+=Vwgt[e];
                assert(ipc.at(owning_proc_id) > 0);

                // Save information which local face points to the other side of subdomain.
                // Store info about which proc it is and what is local id of element on the other side.
                int faces[MMC_MAXELFAC+1]={0};
                mmr_el_faces(Mesh_id, El_id[e], faces, nullptr);

                //int el_neigs[MMC_MAXELFAC+1] = {0};
                //mmr_el_eq_neig(Mesh_id,El_id[e],el_neigs, nullptr);

                int nth_face = std::distance(& Adjncy[*it_xadj], it_neig_el);
                assert(nth_face >= 0);
                nth_face = (Adj_neig+(*it_xadj))[ nth_face ];
                // const int nth_face = std::distance(el_neigs+1, std::find(el_neigs+1,el_neigs+1+el_neigs[0],(*it_neig_el)+1) );
                assert(nth_face <= faces[0]);
                assert(nth_face >= 0);
                const int local_face_id = faces[ 1+ nth_face ];

                assert(local_face_id > 0);
                assert(local_face_id <= mmr_get_max_face_id(Mesh_id));


                // Checking whether result is correct
#ifndef NDEBUG
                int neig[2]={0};
                mmr_fa_neig(Mesh_id,local_face_id,neig, nullptr, nullptr, nullptr, nullptr, nullptr) ;
                //assert(neig[1] == 0);
                //assert(neig[0] == El_id[e]);
#endif
                assert(std::distance(& Adjncy[*it_xadj], it_neig_el) < MMC_MAXELFAC);
                assert(*it_neig_el >= *found_in_vtxdist);

                ipc_links[local_face_id].owner = owning_proc_id;
                // !!!!!!!!!!!!!!!!!!!! TODO: here is an error !!!!!!!!!!!!!
                ipc_links[local_face_id].id_at_owner = (*it_neig_el - VtxDist[owning_proc_id-1]) + 1; //+1 to make proper id [1..n] not [0..n-1]

                mf_check_debug(ipc_links[local_face_id].owner > 0,"Owner proc id cannot be <0 (now %d)",ipc_links[local_face_id].owner);
                mf_check_debug(ipc_links[local_face_id].id_at_owner > 0,"Id at owner cannot be <0 (now %d)",ipc_links[local_face_id].id_at_owner);

                mf_debug("Cross boundary link from el=%d face=%d to proc=%d el=%d",
                            El_id[e],local_face_id,owning_proc_id,ipc_links[local_face_id].id_at_owner);
            }
        }
    }

    mf_debug("No. of cross boundary links=%d",ipc_links.size());

    return retval;
}



/*------------------------------------------------------------
  ddr_expand_overlap - to expand overlap regions by one element layer
------------------------------------------------------------*/
// NOTE: if overlap does not exist at all, it will be created as 1 element overlap.
int mmpr_IPC_expand_overlap(const int Mesh_id)
{
    mmpt_mesh & pmesh = *mmpr_select_mesh(Mesh_id);

    assert(pmesh.data.ipc.empty() == false);
    assert(pmesh.data.ipc_links.empty() == false);

    int retval=0;
    // We base on inter process connecitivty data (and assume that it is correct)
    // For each neigbouring proces, we seek list of elements, that we want to copy TO us.
    // But we also know, which elements this process want to copy FROM us.
    // NOTE: this 1-element overlap is based on the rule that smaller-id-proc
    // does't hava a copy of overlap elems - it just own them as very own.
    // So, proc have a copy of elements only FROM smaller-id-proc, and own elements
    // that are copied from it by the bigger-id-proc.
    // E.g. First proc will never copy anything from any proc, and no proc will copy from last-id-proc.

    std::map< int, std::vector<int> > elem_to_send_to_procs;

    for(const mmpt_data::ipc_links_type::value_type & it_ipc_links : pmesh.data.ipc_links) {
        if(it_ipc_links.second.owner > pcv_my_proc_id) {
            int neigs[2]={0};

            mmr_fa_eq_neig(Mesh_id,it_ipc_links.first,neigs,nullptr,nullptr);
            assert(neigs[0] != 0);
            elem_to_send_to_procs[it_ipc_links.second.owner].push_back(neigs[0]);
        }
    }

    for(const mmpt_data::ipc_type::value_type & it_ipc : pmesh.data.ipc) {

        if(it_ipc.first < pcr_my_proc_id()) { // smaller ids = bigger priority
//            // First send.
//            ddr_copy_full_elems(Mesh_id, pcr_my_proc_id(), it_ipc.first,
//                                elems_to_send.size(),
//                                elems_to_send.data());
            // Then recieve.
            mmpr_transfer_full_elems(Mesh_id, MPI_ANY_SOURCE, pcr_my_proc_id(),
                                0,nullptr, MMPE_COPY);
        }
        else { // it_ipc.first > pcr_my_proc_id() // bigger ids = smaller priority
            std::vector<int> & elems_to_send = elem_to_send_to_procs.at(it_ipc.first);
            // Duplicate entires are possible, so remove them.
            std::sort(elems_to_send.begin(),elems_to_send.end());
            elems_to_send.resize( std::distance( elems_to_send.begin(), std::unique( elems_to_send.begin(), elems_to_send.end()) ) );

            //            // First recieve.
//            ddr_copy_full_elems(Mesh_id, it_ipc.first, pcr_my_proc_id(),
//                                0,
//                                nullptr);
            // Then send.
            mmpr_transfer_full_elems(Mesh_id, pcr_my_proc_id(), it_ipc.first,
                                elems_to_send.size(),
                                elems_to_send.data(), MMPE_COPY);
        }
        mmr_test_mesh(Mesh_id);

        // After that, we have to update inter process connecitivty data.
        // Remove faces no longer at boundary from ipc data.
        // Put overlap faces into ipc data.
        // (validate that, they have MMC_SUB_BND flags set as neigs)
        // const mmpt_PCSR& csr = pmesh.data.mp.PCSR;
    }

    mmr_test_mesh(Mesh_id);

    return retval;
}



/* INTERNAL ROUTINES (TREATED AS PRIVATE) */




///*------------------------------------------------------------
//  mmpr_end_work - free memory for internal data structures
//------------------------------------------------------------*/
//void mmpr_end_work()
//{
//	/*++++++++++++++++ executable statements ++++++++++++++++*/

////  mmpr_free_CRS(& pmesh.data.mp.PCSR);
  
////	UTM_SAFE_FREE_PTR(pmesh.data.mp.vtxdist);
////	UTM_SAFE_FREE_PTR(pmesh.data.mp.tpwgts);
////	UTM_SAFE_FREE_PTR(pmesh.data.mp.part_local);
//	//UTM_SAFE_FREE_PTR(pmesh.data.mp.adjncy_local);
//	//UTM_SAFE_FREE_PTR(pmesh.data.mp.adjwgt_local);
//	//UTM_SAFE_FREE_PTR(pmesh.data.mp.xadj_local);
//}




///*------------------------------------------------------------
//  mmpr_create_overlap - create overlap, fill data structure, set internal
//------------------------------------------------------------*/
//void mmpr_create_overlap(const int MeshID)
//{
//    /*
//        NOTE:
//                when creating overlap the rule is:
//                    the subdoamin with bigger ID contains border elements of smaller id subdomains.
//    */
	

//}



#ifdef __cplusplus
}
#endif
