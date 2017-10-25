/************************************************************************
File apph_dg_prism.h - internal information for parallel overlap for 
   approximation module for dg approximation and prismatic elements       

Contains:
  - constants
  - data types 
  - global variables (for the whole module)
  - function headers                     

------------------------------  			
History:        
        08.2011 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _apph_dg_prism_
#define _apph_dg_prism_

#include "mmh_intf.h"

#include "mmph_intf.h"

#include "aph_intf.h"

#include "apph_intf.h"


#define APPC_MAX_NUM_SUB MMPC_MAX_NUM_SUB // maximal number of subdomains
#define APPC_MAX_NUM_TABLES 10 // several subsystems

/* 
 * data structure to handle exchange of dofs between processors
 * for multi-problem multi-grid solvers 
 */

/* definition of type appt_levels - data structure for mesh levels */
typedef struct {

/* parameters and arrays for exchanging data on interface dof structures */
/* 
   1. data received, updated and send - dof structures owned by the processor,
                                        updated also by others
*/
  //int numexch1[APPC_MAX_NUM_SUB], *dofentexch1[APPC_MAX_NUM_SUB];
  //int *posexch1[APPC_MAX_NUM_SUB], *nrdofexch1[APPC_MAX_NUM_SUB];
/* 
   2. data send and received (after update) - dof structures owned by others,
                                              updated also by the processor
*/
  //int numexch2[APPC_MAX_NUM_SUB], *dofentexch2[APPC_MAX_NUM_SUB]; 
  //int *posexch2[APPC_MAX_NUM_SUB], *nrdofexch2[APPC_MAX_NUM_SUB];
/* 
   3. data updated and send - dof structures owned by the processor, 
                              not updated by others (sometimes equivalent to
                              being on the boundary of others' subdoamins)
*/
  int numsend[APPC_MAX_NUM_SUB], *dofentsend[APPC_MAX_NUM_SUB];
  int *possend[APPC_MAX_NUM_SUB], *nrdofsend[APPC_MAX_NUM_SUB]; 
/* 
   4. data received (after update) - dof structures owned by others, 
                                     not updated by the processor,
				     (sometimes equivalent to the list 
				     of the subdomain's boundary dof structures)
*/
  int numrecv[APPC_MAX_NUM_SUB], *dofentrecv[APPC_MAX_NUM_SUB];
  int *posrecv[APPC_MAX_NUM_SUB], *nrdofrecv[APPC_MAX_NUM_SUB]; 

  int nrdof_int;

} appt_levels;


typedef struct {
  int field_id;
  int nreq_offset; // we select certain DOFs (nreq_offset..nreq_offset+nreq)
  int nreq;
  int nr_levels;           /* number of levels */
  appt_levels *level;       /* array of data structures for levels */
} appt_exchange_tables;

/* GLOBAL VARIABLES */
extern int appv_nr_proc; // number of process(or)s = number of subdomains
extern int appv_my_proc_id; // executing process(or)s ID
extern int appv_nr_exchange_tables;
extern int appv_cur_exchange_table_id;
extern appt_exchange_tables appv_exchange_table[APPC_MAX_NUM_TABLES]; 

#endif
