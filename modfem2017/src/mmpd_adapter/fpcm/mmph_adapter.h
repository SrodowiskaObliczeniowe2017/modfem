#ifndef _mmph_adapter_
#define _mmph_adapter_

#include <metis.h>
#include <parmetis.h>
#include <unordered_map>
#include <vector>
#include <set>

#include "uth_intf.h"

#include "mmph_ipid.h"

#include "Transferer.hpp"

using namespace mmpt;

#ifdef __cplusplus
extern "C"{
#endif

/*** CONSTANTS ***/

#define MAC_USE_KWAY_GRAPH_PART_TOOL   0   /* partitioning tool. Better quality, much time to compute */
#define MAC_USE_RB_GRAPH_PART_TOOL   1   /* partitioning tool. Worse quality, less time to compute */

#define MAC_MAX_ELEMENT_FACES   6   /* max supported facas that that element may have */
#define MAC_MAX_ELEMENT_EGDES   12  /* max supported edges that that element may have. */
#define MAC_MAX_ELEMENT_NODES   8  /* max supported edges that that element may have. */
#define MAC_MAX_NODES_ON_FACE   4  /* max supported edges that that element may have. */
#define MAC_MAX_NODES_COOR_ON_FACE   12  /* max supported edges that that element may have. */

#define MAC_MAX_ELEMENTS_CONNECTED_BY_NODE   384  /* max supported edges that that element may have. 6 faces x 3 genlev (8 el) ^ 2 */
#define MAC_MAX_PARTITIONS   128  /* max supported edges that that element may have. */

#define MAC_PARALLEL_MESH_CORASING   1  /* use OpenMP to speedup lookup_table building */

#define MAC_PRISM_OPTIMIZATION  1  /* use weights during decomposition that indicates amount of data to exchange betweene each face in prism elements */

#define MMC_MAX_NUM_MESH   10   //* maximal number of meshes

#define MMPC_CLASSIC_PART_METHOD   0   //* use classic (old) partition method
#define MMPC_METIS_PART_METHOD   1   //* use Metis and ParMetis libraries to partition a mesh



/*** GLOBAL VARIABLES for the whole module ***/

extern int       mmpv_partition_method;   /*  method that is used to generate partitions */
extern std::vector<mmpt_mesh>  mmpv_meshes;        /* array of meshes */

/*** LOCAL PROCEDURES  ***/

/**--------------------------------------------------------
  mmpr_select_mesh - to select the proper mesh
---------------------------------------------------------*/
mmpt_mesh* mmpr_select_mesh( /* returns pointer to the chosen mesh */
               /* to avoid errors if input is not valid */
               /* it returns the pointer to the current mesh */
  int Mesh_id    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );



/**-----------------------------------------------------------
  mmpr_transfer_full_elems - to move element from one process (mesh) to another process (mesh)
------------------------------------------------------------*/
int mmpr_transfer_full_elems(mmpt_mesh & pmesh,
        const int Source_proc_id,
        const int Dest_proc_id,
        const int N_transfer_elems,
        const int * Transfer_elem_ids,
        const TRANSFER_POLICY T_policy);

#ifdef __cplusplus
}
#endif

#endif
