/************************************************************************
File ddh_metis_adapter.h - domain decomposition adapter for external mesh 
						   partitioning libraries Metis and ParMetis       

Contains:
  - constants    
  - data types 
  - function headers     
------------------------------  			
History:        
        08.2012 - Kamil Wachala, initial version
		03.2013 - Kazimierz Michalik, incorporating into ModFEM
*************************************************************************/

#ifndef _ddd_metis_adapter_
#define _ddd_metis_adapter_

#include<map>

/* mesh manipulation interface */
#include "mmh_intf.h"

/* parallel mesh manipulation interface */
#include "mmph_intf.h"

// uses
#include "../mmph_adapter.h"

#ifdef __cplusplus
extern "C"{
#endif


/*** LOCAL PROCEDURES  ***/
//--------------------------------------------------------
// ddr_mesh_to_CRS_graph - to save initial mesh (generation=level=1) as graph in compressed storage format
//-------------------------------------------------------
int ddr_mesh_to_CRSGraph(// return >=0 - number of vertices in CSR graph, <0 - error
                         const int Mesh_id);

int ddr_check_CSR(const int Mesh_id) ;

int ddr_PCSR_create_initial_distribution(
        const int Mesh_id);

int mmpr_PCSR_distribute_elems_to_procs(const int Mesh_id,
        const int Source_proc_id);









/**-----------------------------------------------------------
  ddr_transfer_copy_elem - to copy (original remains) element from one process (mesh) to another process (mesh)
------------------------------------------------------------*/
//int ddr_copy_full_elems(
//        const int Mesh_id,
//        const int Source_proc_id,
//        const int Dest_proc_id,
//        const int N_copied_elems,
//        const int *Copied_elem_ids);

/**-----------------------------------------------------------
  ddr_create_inter_subdomain_connectivity - to create and initalize inter subdomain data
------------------------------------------------------------*/
/// Inter subdomain connectivity provides info about
/// 1. general connections info between subdomains
/// 2. detailed info about which [element? face?] is on the other side of subdomain boundary
/// NOTE: this works only for initial elements (father == MMC_NO_FATHER)
int mmpr_PCSR_create_IPC(// return >=0 - number of vertices in CSR graph for local proc, <0 - error
         const int Mesh_id);
/**-----------------------------------------------------------
  mmpr_partition_mesh - to to decompose the mesh and create subdomains by Metis (SFC - Parmetis) library
------------------------------------------------------------*/
int ddr_PCSR_partition_mesh( /* returns: 1 - success, <0 - error code */
  int Mesh_id, /* ID of the current mesh */
  int Part_amount, /* The number of parts to partition the mesh  */
  int Parition_tool /* Used partition tool. Can be: MAC_USE_KWAY_GRAPH_PART_TOOL or MAC_USE_RB_GRAPH_PART_TOOL */
  );

/**-----------------------------------------------------------
  mmpr_adaptive_repartition - to balance the work load of a graph that corresponds to an adaptively refined mesh
------------------------------------------------------------*/
int mmpr_adaptive_repartition( /* returns: 1 - success, <0 - error code */
  int Mesh_id /* ID of the current mesh */
  );

 /**-----------------------------------------------------------
  mmpr_improve_partitioning - to improve the quality of an existing a k-way partitioning
------------------------------------------------------------*/
int ddr_PCSR_improve_partitioning( /* returns: 1 - success, <0 - error code */
  int Mesh_id /* ID of the current mesh */
  );

/**-----------------------------------------------------------
  mmpr_set_metis_options - to specify metis behaviour
------------------------------------------------------------*/
void mmpr_set_metis_options(const int Mesh_id, idx_t options[]);

/**-----------------------------------------------------------
  mmpr_set_parmetis_options - to specify parmetis behaviour
------------------------------------------------------------*/
void mmpr_set_parmetis_options(const int Mesh_id, idx_t options[]);

/**-----------------------------------------------------------
  mmpr_initialize_work - to initialize the internal data structures
------------------------------------------------------------*/
void mmpr_init_data(
  int Mesh_id, /* ID of the current mesh */
  int Part_amount /* The number of parts to partition the mesh  */
  );

/**-----------------------------------------------------------
  mmpr_end_work - free memory for internal data structures
------------------------------------------------------------*/
void mmpr_end_work();

/**-----------------------------------------------------------
  mmpr_create_overlap - create overlap, fill data structure, set internal
------------------------------------------------------------*/
void mmpr_create_overlap(const int MeshID);

/**-----------------------------------------------------------
  ddr_update_overlap - to update one element overlap layer
------------------------------------------------------------*/
// NOTE: if overlap does not exist at all, it will be created as 1 element overlap.
int mmpr_IPC_update(const int Mesh_id,
                   const int* ElemsChanged,
                   const int nElems
                       );
/**-----------------------------------------------------------
  ddr_expand_overlap - to expand overlap regions by one element layer
------------------------------------------------------------*/
// NOTE: if overlap does not exist at all, it will be created as 1 element overlap.
int mmpr_IPC_expand_overlap(const int Mesh_id
                       );


#ifdef __cplusplus
}
#endif


  
#endif
