#include <algorithm>
#include <vector>

#include <mpi.h>
#include <metis.h>
#include <parmetis.h>


#include "ddh_intf.h"
// #include "ddd_metis_adapter/ddh_metis_adapter.h"
#include "uth_log.h"
#include "PCsr.hpp"
#include "uth_io_results.h"

typedef CSR<idx_t> ddt_CSR;
typedef PCSR<idx_t> ddt_PCSR;

typedef struct
{
  idx_t options[METIS_NOPTIONS]; /* allow to fine-tune and modify various aspects of the internal algorithms used by METIS */
  ddt_PCSR  PCSR; // one graph for each sub-domain

  // idx_t* vtxdist; /* describes how the vertices of the graph are distributed among the processors */
  idx_t* vsize;
//  idx_t** part_local ; /* upon successful completion stores the local partition vector of the mesh (ids for each element) */
  std::vector<real_t> tpwgts ; /* fraction of vertex weight that should be distributed to each sub-domain for each balance constraint */
  real_t ubvec ; /* specify the imbalance tolerance for each vertex weight, with 1 being perfect balance and nparts being perfect imbalance. A value of 1.05 recommended */
  idx_t wgtflag ; /* indicate if the graph is weighted. 0 - No weights, 1 - weights on edges */
  idx_t numflag ; /* numbering scheme that is used for the vtxdist, xadj, adjncy, and part. 0 C-style numbering */
  real_t itr ;  /* ratio of inter-processor communication time compared to data redistribution time. Value of 1000.0 is recommended*/

  MPI_Comm comm; /* pointer to the MPI communicator of the processes that call PARMETIS */

  /* Metis and Parmetis */
  std::vector<idx_t> part; /* upon successful completion stores the partition vector of the mesh (ids for each element) */
  idx_t nparts ; /* The number of parts to partition the mesh */
  idx_t objval ; /* after complition: Stores the edge-cut or the total communication volume of the partitioning solution  */
  idx_t ncon ; /* the number of balancing constraints */
  idx_t met_parmet_result ; /* the value returned by metis or parmetis routines */
  idx_t edgecut;

  void clear() {
     UTM_SAFE_FREE_PTR(vsize);
     tpwgts.clear();
     PCSR.clear();
     part.clear();
  }

} ddt_partitioning_data;


std::vector<ddt_partitioning_data> ddv_meshes_part_data;

int ddr_check_CSR(const int Mesh_id, ddt_CSR & Csr)
{
//    mmpt_CSR &Csr = mmpr_select_mesh(Mesh_id)->data.mp.PCSR;

    mf_check(Csr.size() == mmr_get_nr_elem(Mesh_id), "Incorrect size(=%d) of CSR (should be %d)",Csr.size(),mmr_get_nr_elem(Mesh_id));

    int nel=0;
    int neig[6]={0};

    const idx_t* it_xadj = Csr.xadj();
    const idx_t* it_adjncy = Csr.adjncy();


    while( (nel = mmr_get_next_act_elem(Mesh_id,nel)) != 0 ) {
        mmr_el_eq_neig(Mesh_id,nel,neig,NULL);

        // NOTE: in pcsr bc_cond neigbors (=0) are omitted, therefore pcsr neighbors table can by shorter
        const int n_neigs = *(it_xadj+1) - *(it_xadj);
        assert(n_neigs == neig[0] - std::count(neig+1,neig+1+neig[0],0));

        for(int i=0; i < neig[0]; ++i) {
            if(neig[1+i] != 0) {
                assert(neig[1+i] == ((*it_adjncy) + 1));
                ++it_adjncy;
            }
        }

        ++it_xadj;
    }
    return 0;
}

//--------------------------------------------------------
// ddr_mesh_to_CRS_graph - to save initial mesh (generation=level=1) as graph in compressed storage format
//-------------------------------------------------------
/*
  All of the graph partitioning and sparse matrix ordering routines in METIS take as input the adjacency structure of the
graph and the weights of the vertices and edges (if any).
  The adjacency structure of the graph is stored using the compressed storage format (CSR). The CSR format is a
     widely used scheme for storing sparse graphs. In this format the adjacency structure of a graph with n vertices and
                                                                                          m edges is represented using two arrays xadj and adjncy. The xadj array is of size n + 1 whereas the adjncy
                                                                                          array is of size 2m (this is because for each edge between vertices v and u we actually store both (v; u) and (u; v)).
          The adjacency structure of the graph is stored as follows. Assuming that vertex numbering starts from 0 (C style),
          then the adjacency list of vertex i is stored in array adjncy starting at index xadj[i] and ending at (but not
                                                                                                                 including) index xadj[i+1] (i.e., adjncy[xadj[i]] through and including adjncy[xadj[i+1]-1]). That
          is, for each vertex i, its adjacency list is stored in consecutive locations in the array adjncy, and the array xadj
  is used to point to where it begins and where it ends.
  The weights of the vertices (if any) are stored in an additional array called vwgt. If ncon is the number of weights
  associated with each vertex, the array vwgt contains n ncon elements (recall that n is the number of vertices). The
  weights of the ith vertex are stored in ncon consecutive entries starting at location vwgt[i ncon]. Note that if
  each vertex has only a single weight, then vwgt will contain n elements, and vwgt[i] will store the weight of the
  ith vertex. The vertex-weights must be integers greater or equal to zero. If all the vertices of the graph have the same
  weight (i.e., the graph is unweighted), then the vwgt can be set to NULL.
*/
//-------------------------------------------------------
bool ddr_is_greater0(const int i) { return i > 0; }
int ddr_mesh_to_CRSGraph(// return >=0 - number of vertices in CSR graph, <0 - error
                         const int Mesh_id, // IN: mesh to compress
                         ddt_CSR & Csr
                         )
{
//    mmpt_CSR &Csr = mmpr_select_mesh(Mesh_id)->data.mp.PCSR;
    Csr.clear();
    Csr.reserve(mmr_get_nr_elem(Mesh_id));

    std::vector<int> sons(mmr_get_max_gen(Mesh_id)*MMC_MAXELSONS, 0);
    // loop over all elements, choose only initial elems
    ddt_CSR::CsrNode gh_nd;
    while( (gh_nd.el_id  = mmr_get_next_elem_all( Mesh_id, gh_nd.el_id )) > 0 ) {
        assert(gh_nd.el_id > 0);
        if( mmr_el_gen( Mesh_id, gh_nd.el_id) == MMC_INIT_GEN_LEVEL ) {
            // Determine how many sub-elements are in this initial elem.
            int father = mmr_el_fam_all( Mesh_id, gh_nd.el_id, sons.data() );
            assert(father == MMC_NO_FATH);

            // Weight = number of sub-elements (created during refinement)
            gh_nd.vwgt = sons[0]+1;

            // Gather neigbours info and count real neighbours (filter out boundary conditions flags).
            int neigs[MMC_MAXELFAC+1] = {0};
            int neigs_no[MMC_MAXELFAC+1] = {0,1,2,3,4,5};
            int rc=mmr_el_eq_neig( Mesh_id, gh_nd.el_id, neigs, NULL);
            assert(rc >= 1); // success
            // below: +1 is because (neighs[0] == no. of all neigs), not neig id!
//            gh_nd.n_neighs = std::distance(
//                        gh_nd.neighs,
//                        std::copy_if(neigs+1, neigs+1+neigs[0], gh_nd.neighs, ddr_is_greater0)
//                    );

            int *first=neigs+1,
                    *last=neigs+1+neigs[0],
                    *result=gh_nd.neighs;
            while (first!=last) {
              if (ddr_is_greater0(*first)) {
                *result = *first;
                ++result;
              }
              ++first;
            }
            gh_nd.n_neighs = std::distance(gh_nd.neighs,result);




            if(gh_nd.n_neighs != neigs[0]) {
                for(int i=1,i2=0; i <= neigs[0]; ++i) {
                    if(neigs[i] != 0) { //bc_cond
                        gh_nd.neig_no[i2] = neigs_no[i-1];
                        ++i2;
                    }
                }
            }
            else {
                std::copy(neigs_no,neigs_no+gh_nd.n_neighs,gh_nd.neig_no);
            }

            // Because we number elements from 1, and storing from 0
            // so neighs must be decresed by 1.
//            std::transform(gh_nd.neighs, gh_nd.neighs+gh_nd.n_neighs, gh_nd.neighs ,[](ddt_CSR::Tind & n){ return --n;});

            first = gh_nd.neighs;
            last = gh_nd.neighs+gh_nd.n_neighs;
            result = gh_nd.neighs;

            while (first != last) {
              *result = --(*first);
              ++result; ++first;
            }

            Csr.push_back(gh_nd);
            //std::cout << Csr;
            assert(Csr.size() <= mmr_get_nr_elem(Mesh_id));
        }//!if
    }//!while

    ddr_check_CSR(Mesh_id, Csr);

    return Csr.size();
}


/*------------------------------------------------------------
  ddr_set_metis_options - to specify metis behaviour
------------------------------------------------------------*/
void ddr_set_metis_options(ddt_partitioning_data & data)
{
    data.ncon = 1; // one constraint. Default value
    data.objval = 0; // out parameter. Reset.
    data.itr = 1000; // default value
    data.wgtflag = 1; // weights on the edges only

    /*++++++++++++++++ executable statements ++++++++++++++++*/

    METIS_SetDefaultOptions(data.options); // reset to defalult options for metis

    data.options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    /*
        Specifies the type of objective. Possible values are:
        METIS_OBJTYPE_CUT Edge-cut minimization.
        METIS_OBJTYPE_VOL Total communication volume minimization.
    */

    data.options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    /*
        Specifies the matching scheme to be used during coarsening. Possible values are:
        METIS_CTYPE_RM Random matching.
        METIS_CTYPE_SHEM Sorted heavy-edge matching.
    */

    data.options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW;
    /*
        Determines the algorithm used during initial partitioning. Possible values are:
        METIS_IPTYPE_GROW Grows a bisection using a greedy strategy.
        METIS_IPTYPE_RANDOM Computes a bisection at random followed by a refinement.
        METIS_IPTYPE_EDGE Derives a separator from an edge cut.
        METIS_IPTYPE_NODE Grow a bisection using a greedy node-based strategy.
    */

    data.options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
    /*
        Determines the algorithm used for refinement. Possible values are:
        METIS_RTYPE_FM FM-based cut refinement.
        METIS_RTYPE_GREEDY Greedy-based cut and volume refinement.
        METIS_RTYPE_SEP2SIDED Two-sided node FM refinement.
        METIS_RTYPE_SEP1SIDED One-sided node FM refinement.
    */

    data.options[METIS_OPTION_MINCONN] = 1;
    /*
     *Specifies that the partitioning routines should try to minimize the maximum degree of the subdomain graph,
    i.e., the graph in which each partition is a node, and edges connect subdomains with a shared interface.
    0 Does not explicitly minimize the maximum connectivity.
    1 Explicitly minimize the maximum connectivity.
    */

    data.options[METIS_OPTION_CONTIG] = 1;
    /*
     *Specifies that the partitioning routines should try to produce partitions that are contiguous. Note that if the
    input graph is not connected this option is ignored.
    0 Does not force contiguous partitions.
    1 Forces contiguous partitions.
    */

    data.options[METIS_OPTION_DBGLVL] = METIS_DBG_TIME; //METIS_DBG_INFO | METIS_DBG_MOVEINFO;
    /*
    Specifies the amount of progress/debugging information will be printed during the execution of the algo-
    rithms. The default value is 0 (no debugging/progress information). A non-zero value can be supplied that
    is obtained by a bit-wise OR of the following values.
    METIS_DBG_INFO (1) Prints various diagnostic messages.
    METIS_DBG_TIME (2) Performs timing analysis.
    METIS_DBG_COARSEN (4) Displays various statistics during coarsening.
    METIS_DBG_REFINE (8) Displays various statistics during refinement.
    METIS_DBG_IPART (16) Displays various statistics during initial partitioning.
    METIS_DBG_MOVEINFO (32) Displays detailed information about vertex moves during refine-
                           ment.
    METIS_DBG_SEPINFO (64) Displays information about vertex separators.
    METIS_DBG_CONNINFO (128) Displays information related to the minimization of subdomain
                            connectivity.
    METIS_DBG_CONTIGINFO (256) Displays information related to the elimination of connected com-
                              ponents.
    Note that the numeric values are provided for use with the -dbglvl option of M ETIS’ stand-alone pro-
    grams. For the API routines it is sufficient to OR the above constants.
    */
}

/*------------------------------------------------------------
  mmpr_set_parmetis_options - to specify parmetis behaviour
------------------------------------------------------------*/
void ddr_set_parmetis_options(ddt_partitioning_data & data)
{
    /*++++++++++++++++ executable statements ++++++++++++++++*/
    ddr_set_metis_options(data);

    data.ubvec = 1.05f;

    data.options[0] = 1; // disable default options
    data.options[1] = 0; // random number seed
    data.options[2] = PARMETIS_PSR_COUPLED; // # of partitions == # of processors ==  # of processes
}


////////////////////////// interface /////////////////////////////////

#ifdef __cplusplus
extern "C" {
#endif

int ddr_create_subdomains_scheme(const int Mesh_id,
              ddt_mesh_handling Mesh_handling,
              int N_subdomains,
              int * N_subdomains_elems,
              int * Subdomains_elems,
              int * Overlap_sizes,
              int * Subdomains_elems_overlap,
              int **Part_ptr)   // numbering from 0
{

    utr_io_result_cumulative_timer_start(RESULT_TIME_DDM);

    mf_check_mem(N_subdomains_elems);
    mf_check_mem(Subdomains_elems);


    if(ddv_meshes_part_data.size() <= Mesh_id) {
        ddv_meshes_part_data.resize(Mesh_id+1);
    }

    ddt_partitioning_data & data = ddv_meshes_part_data.at(Mesh_id);

    ddr_mesh_to_CRSGraph(Mesh_id,data.PCSR);

    ddr_set_parmetis_options(data);

    switch(Mesh_handling) {
        case DDC_GLOBAL_MESH:
        mf_log_err("Global mesh not yet supported!");


//        const int n_procs = pcr_nr_proc();
//        std::vector<int> vtcs_distrib(n_procs,0);
//        pcr_allgather_int( & pmesh.data.ne , 1, vtcs_distrib.data(), 1 );

//        // There is a possibilty that mesh has poor initial vertex distribution
//        // (there are some procs with 0 verteices at begining).
//        // Unfortunetly ParMETIS requires vertices to be pre-distributed.
//        // So we have to create some basic distribution to
//        // assign each proc at least one vertex.
//        // So we are looking for procs with no vertices.
//        // When found, previous proc mesh is partitioned.
//        //
//        // We're begining from 1, becouse first process always have mesh.
//        for( int p_rank=1 ; p_rank < pmesh.data.mp.nparts ; ++p_rank) {
//            if(vtcs_distrib[p_rank] == 0) {
//                // so p_id-th proc is first proc with no vertices
//                // last proc with vertices
//                const int proc_w_mesh_id=p_rank; // previous proc

//                if(proc_w_mesh_id == pcr_my_proc_id()) {
//                    // use METIS to sequentialy partition mesh
//                    int n_novert_procs = n_procs-p_rank+1;
//                   ddr_create_subdomains(Mesh_id,DDC_LOCAL_MESH,n_novert_procs,NULL,NULL,NULL);
//                }
//            }
//        }


//        data.met_parmet_result = ParMETIS_V3_PartKway(
//                    data.PCSR.vtxDist(),
//                    data.PCSR.xadj(),
//                    data.PCSR.adjncy(),
//                    data.PCSR.vwgt(),
//                    data.PCSR.adjwgt(),
//                    & data.wgtflag,
//                    & data.numflag,
//                    & data.ncon,
//                    & data.nparts,
//                    data.tpwgts.data(),
//                    & data.ubvec,
//                    data.options,
//                    & data.objval,
//                    data.PCSR.part(),
//                    & data.comm
//                    );

        break;

    case DDC_LOCAL_MESH:
    case DDC_DEFAULT:
        int n_e = data.PCSR.size();
        idx_t* part = data.PCSR.part();
        // for naive problems
        if(N_subdomains == n_e) {
            for(int n=0; n < n_e; ++n) {
                part[n] = n;
            }
        }
        else {
            mf_debug("Lunching METIS.");
            data.met_parmet_result = METIS_PartGraphKway(& n_e,
                                    & data.ncon,
                                    data.PCSR.xadj(),
                                    data.PCSR.adjncy(),
                                    data.PCSR.vwgt(),
                                    data.PCSR.vwgt(),//pmesh.data.mp.vsize,
                                    NULL, // all edges in graph have same weight
                                    & N_subdomains,
                                    data.tpwgts.data(),
                                    & data.ubvec,
                                    data.options,
                                    & data.objval,
                                    part);
        }
        break;
    }

    // Splitting initial PCSR into per-proc-PCSRs
    // Waring: element indexes in per-proc-PCSRs will be different!
    mf_debug("Splitting PCRS");
    typedef ddt_CSR::CsrNodeInternal cNode;
    std::vector<ddt_PCSR> procs_pcsr(N_subdomains);
//    std::for_each(data.PCSR.begin(),data.PCSR.end(), [&](cNode & n) {
//        // split into subvectors
//        procs_pcsr[*n.ppart].push_back(n);
//    });

    for(ddt_PCSR::iterator it = data.PCSR.begin();
        it != data.PCSR.end();
        ++it) {
         procs_pcsr[*(*it).ppart].push_back(*it);
    }

    mf_debug("Rewriting PCRS output.");
    int *current = Subdomains_elems;
    for(int p=0; p < procs_pcsr.size();++p) {
        N_subdomains_elems[p] = procs_pcsr[p].size();
        if(N_subdomains_elems[p] < 1) {
            mf_fatal_err("Subdomain %d has no elements assigned!", p);
        }
        memcpy(current, procs_pcsr[p].el_id(), sizeof(int) * N_subdomains_elems[p]);
        current += N_subdomains_elems[p];
    }

    if((Subdomains_elems_overlap != NULL) && (Overlap_sizes!=NULL) ) {
       int * p_cur_overlap = Subdomains_elems_overlap;

        std::vector<int>  overlap;
        for(int p=0; p < N_subdomains; ++p) {
            data.PCSR.getOverlapElems(overlap,p);

            Overlap_sizes[p]=overlap.size();
            memcpy( p_cur_overlap, overlap.data(), sizeof(int) * overlap.size() );
            p_cur_overlap += overlap.size();
        }

    }

    if(Part_ptr != NULL) {
        if(*Part_ptr == NULL) {
            *Part_ptr = data.PCSR.part();
        }
        else {
            mf_log_err("Expected NULL pointer not referencing to memory adress.");
        }
    }

    utr_io_result_cumulative_timer_stop(RESULT_TIME_DDM);
    return N_subdomains;
}

int ddr_balance_subdomains(const int Mesh_id,
                     int * N_subdomains_elems,
                     int *Subdomains_elems,
                     int ** Subdomains_elems_overlap)
{
    ddt_partitioning_data & data = ddv_meshes_part_data.at(Mesh_id);

    // Once created structures are still valid.
    // This is because we create them for initial elements,
    // and that is NOT CHANGING during computations.
    // The only tables that needs updating
    // are vwgt(vertice weight) and vsize (v. redistribution cost) table.
    // For each initial element,
    // its' vsize[i] have to correspond with number
    // of DOFs for this initial element and all its children.
    // We assume that vsize = vwgt.
    //assert(mmpv_data.mp.vsize != NULL );
    int rc = ParMETIS_V3_AdaptiveRepart(
                data.PCSR.vtxDist(),
                data.PCSR.xadj(),
                data.PCSR.adjncy(),
                data.PCSR.vwgt(),
                data.PCSR.vwgt(),
                data.PCSR.adjwgt(),
                & data.wgtflag,
                & data.numflag,
                & data.ncon,
                & data.nparts,
                data.tpwgts.data(),
                & data.ubvec,
                & data.itr,
                data.options,
                & data.edgecut,
                data.part.data(),
                & data.comm);

    return rc;
}

#ifdef __cplusplus
}
#endif

