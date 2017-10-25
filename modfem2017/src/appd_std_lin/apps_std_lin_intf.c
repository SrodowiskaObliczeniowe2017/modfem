/************************************************************************
File apps_std_lin_intf.c - implementation of the interface apph_intf.h,
                            an overlay for parallel MPI execution of
                            approximation module routines - in this case
                            for standard linear approximation (apd_std_lin)

Procedures:
  appr_get_ent_owner - to return owning process(or) ID
  appr_get_ent_id_at_owner

  appr_init_exchange_tables - to initialize data structure related to exchange 
                             of dofs
  appr_create_exchange_tables - to create lists of dof structures
                              exchanged between pairs of processors.
  appr_exchange_dofs - to exchange dofs between processors
  appr_free_exchange_tables - to free data structure related to exchange 
                             of dofs

  appr_sol_vec_norm - to compute a norm of global vector in parallel
  appr_sol_sc_prod - to compute a scalar product of two global vectors 
  appr_get_nr_dofs_owned - utility to calculate the number of dofs owned at
                          the finest level

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
        08.2015 - KB, new structure for multi-subsystem, multigrid solvers
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#include "mmh_intf.h"

#include "mmph_intf.h"

#include "aph_intf.h"

#include "apph_intf.h"

#include "./apph_std_lin.h"

/* interface for parallel communication modules */
#include "pch_intf.h"

#include "lin_alg_intf.h"

#include "uth_log.h"
#include "uth_io_results.h"
#include "uth_system.h"

/* GLOBAL VARIABLES */
int appv_nr_proc=-1; // number of process(or)s = number of subdomains
int appv_my_proc_id=-1; // executing process(or)s ID
int appv_nr_exchange_tables = 0;
int appv_cur_exchange_table_id = -1;
appt_exchange_tables appv_exchange_table[APPC_MAX_NUM_TABLES]; 

/*---------------------------------------------------------
appr_get_ent_owner - to return owning process(or) ID
---------------------------------------------------------*/
extern int appr_get_ent_owner( 
                   /* returns: >=0 -success code - owner ID, <0 -error code */
  int Field_id,    /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{


#ifdef DEBUG_APM
  if(Ent_type!=APC_VERTEX){
    printf("DOF entity type != APC_VERTEX in apps_std_lin_intf! Exiting!\n");
    exit(-1);
  }
#endif

  int mesh_id;
  mesh_id=apr_get_mesh_id(Field_id);
  // for std_lin DOFs are associated with vertices only
  mmpr_ve_owner(mesh_id,Ent_id);
  return(0);
}

/*---------------------------------------------------------
appr_get_ent_id_at_owner
---------------------------------------------------------*/
extern int appr_get_ent_id_at_owner( 
        /* returns: >=0 -success code - entity local ID, <0 -error code */
  int Field_id,    /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{

#ifdef DEBUG_APM
  if(Ent_type!=APC_VERTEX){
    printf("DOF entity type != APC_VERTEX in apps_std_lin_intf! Exiting!\n");
    exit(-1);
  }
#endif

  int mesh_id;
  mesh_id=apr_get_mesh_id(Field_id);
  // for std_lin DOFs are associated with vertices only
  mmpr_ve_id_at_owner(mesh_id,Ent_id);
  return(0);
}


/*---------------------------------------------------------
  appr_init_exchange_tables - to initialize data structure related to exchange 
                             of dofs
---------------------------------------------------------*/
int appr_init_exchange_tables( /* returns: ID of created tables */
  int Nr_proc,      /* number of process(or)s */
  int My_proc_id,   /* executing process(or)s ID */
  int Field_id,
  int Nreq_offset,
  int Nreq,
  int Nr_levels    /* in: number of levels for each field */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  appv_nr_proc = Nr_proc;
  appv_my_proc_id = My_proc_id;
  
  /* proceed only if necessary */
  if(appv_nr_proc==1) return(0);
   
  utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );

  appv_nr_exchange_tables++;
  appv_cur_exchange_table_id = appv_nr_exchange_tables-1;

  if(appv_nr_exchange_tables >= APPC_MAX_NUM_TABLES){
    printf("Increase APPC_MAX_NUM_TABLES in apph_dg_prism.h! Exiting!\n");
    exit(-1);
  }


  int table_id = appv_cur_exchange_table_id;
  appv_exchange_table[table_id].field_id = Field_id;
  appv_exchange_table[table_id].nreq_offset = Nreq_offset;
  appv_exchange_table[table_id].nreq = Nreq;
  appv_exchange_table[table_id].nr_levels = Nr_levels;
  appv_exchange_table[table_id].level = 
      (appt_levels *)calloc(Nr_levels,sizeof(appt_levels));
  
/*kbw
  printf("Initialized exchange tables for field %d\n", Field_id);
  printf("index %d, nreq_offset %d, nreq %d, nr_levels %d\n", 
	 table_id, Nreq_offset, Nreq, Nr_levels); 
/*kew*/
  utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM );
  return(table_id);
}


/*---------------------------------------------------------
  appr_create_exchange_tables - to create lists of dofs
                              exchanged between pairs of processors.
  REMARK. a simplified setting is used, where dof structures are modified
    (updated) only by their owners and data send to other processors.
    Thus both exchange1 and exchange2 groups of arrays are empty and
    only send and receive arrays are filled. In the routine
    appr_create_exchange_tables each processor constructs lists of dof
    structures needed for solver (dofrecv group) and sends them to owning 
    processors. It also receives from other processors requests for owned 
    dof structures needed by other processors (dofsend group).
---------------------------------------------------------*/
int appr_create_exchange_tables( 
                   /* returns: >=0 -success code, <0 -error code */
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /* in: level ID */
  int Nr_dof_ent,  /* in: number of DOF entities in the level */
  int* L_dof_ent_type,/* in: list of DOF entities associated with DOF structs */
  int* L_dof_ent_id,  /* in: list of DOF entities associated with DOF structs */
  int* L_struct_nrdof,    /* in: list of nrdofs for each dof struct */
  int* L_struct_posg,     /* in: list of positions within the global */
                      /*     vector of dofs for each dof struct */
  int* L_elem_to_struct,/*in: list of DOF structs associated with DOF entities */
  int* L_face_to_struct,/*in: list of DOF structs associated with DOF entities */
  int* L_edge_to_struct,/*in: list of DOF structs associated with DOF entities */
  int* L_vert_to_struct /*in: list of DOF structs associated with DOF entities */
  )
{

  appt_levels* level_p; /* pointer to level structure */

  int i,j,k, iaux, istr, istr_int, field_id, mesh_id, message_id, buffer_id;
  int nr_proc, my_proc_id, nreq_offset, nreq;
  int isub, iproc, nno, ipid=0, struct_id, owner_id, no_id_local, no_id_owner;
  double daux;

/*++++++++++++++++ executable statements ++++++++++++++++*/


  /* proceed only if necessary */
  if(appv_nr_proc==1) return(0);

  utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );

  /* rewrite to local variables */
  nr_proc = appv_nr_proc;
  my_proc_id = appv_my_proc_id;
  appv_cur_exchange_table_id = Table_index;

/*kbw
printf("Processor (subdomain) %d entering create_exchange_tables %d for level %d\n", 
       my_proc_id, Table_index, Level_id);
/*kew*/

  if(Level_id >= appv_exchange_table[Table_index].nr_levels){
    printf("Wrong Level_id %d in appr_create_exchange_tables! Exiting!\n", Level_id);
    exit(-1);
  }

  /* pointer to particular level structure */
  level_p = &appv_exchange_table[Table_index].level[Level_id];

  /* associate the suitable mesh with approximation field */ 
  field_id = appv_exchange_table[Table_index].field_id;
  nreq_offset = appv_exchange_table[Table_index].nreq_offset;
  nreq = appv_exchange_table[Table_index].nreq; 
  mesh_id = apr_get_mesh_id(field_id);

  /* initialize */
  for(isub=0;isub<nr_proc;isub++){
    level_p->numrecv[isub]=0;
  }


  /* for all dof structures (for STD_LIN just vertices=nodes) */
  for(istr=0;istr<Nr_dof_ent;istr++){

    nno=L_dof_ent_id[istr];
 
    owner_id = mmpr_ve_owner(mesh_id, nno);
    if(owner_id != my_proc_id){

      /* set info on the dof structure */
      no_id_local = nno;
      no_id_owner = mmpr_ve_id_at_owner(mesh_id, nno);

      level_p->numrecv[owner_id-1]++;

    }
  }

  for(isub=0;isub<nr_proc;isub++){
    if(level_p->numrecv[isub]>0){
      level_p->dofentrecv[isub] = (int *)calloc(level_p->numrecv[isub],sizeof(int));
      level_p->posrecv[isub] = (int *)calloc(level_p->numrecv[isub],sizeof(int));
      level_p->nrdofrecv[isub] = (int *)calloc(level_p->numrecv[isub],sizeof(int));
    }
  }

  /* again initialize */
  for(isub=0;isub<nr_proc;isub++){
    level_p->numrecv[isub]=0;
  }

  /* for all dof structures (for STD_LIN just vertices=nodes) */
  istr_int=-1;
  level_p->nrdof_int=0;
  for(istr=0;istr<Nr_dof_ent;istr++){

    nno=L_dof_ent_id[istr];
 
    owner_id = mmpr_ve_owner(mesh_id, nno);
    if(owner_id != my_proc_id){

      /* set info on the dof structure */
      no_id_local = nno;
      no_id_owner = mmpr_ve_id_at_owner(mesh_id, nno);

/*kbw
      printf("struct %d, element %d, ipid %d, owner %d, id_owner %d, nrdof %d, posg %d\n",
	     istr, nno, ipid, owner_id, no_id_owner, 
	     L_struct_nrdof[istr], L_struct_posg[istr]);
/*kew*/

      struct_id = level_p->numrecv[owner_id-1];
      level_p->dofentrecv[owner_id-1][struct_id]=no_id_owner;
      level_p->posrecv[owner_id-1][struct_id]=L_struct_posg[istr];
      level_p->nrdofrecv[owner_id-1][struct_id]=L_struct_nrdof[istr];
      level_p->numrecv[owner_id-1]++;

    }
    else{

      /* check whether internal dofs are at the beginning */
      istr_int++;
      if(istr_int==istr) level_p->nrdof_int+=L_struct_nrdof[istr];
      else level_p->nrdof_int = 0;

/*kbw
  printf("Processor (subdomain) %d, level %d, node %d (%d), nr_dofs %d s\n",
	 my_proc_id, Level_id, istr, nno, L_struct_nrdof[istr]);
  printf(" added to %d consecutive internal dofs\n", level_p->nrdof_int);
/*kew*/

    }

  }

/*kbw
  printf("Processor (subdomain) %d, level %d, %d consecutive internal dofs\n",
	 my_proc_id, Level_id, level_p->nrdof_int);
/*kew*/

  //  message_id = PCC_CR_EX;
  message_id = 1234+Level_id;

  for(isub=0;isub<nr_proc;isub++){

    // omitting sendig to self
    if(isub != appv_my_proc_id-1) {

        buffer_id = pcr_send_buffer_open(message_id,0);
        pcr_buffer_pack_int(message_id, buffer_id, 1, &my_proc_id);
        pcr_buffer_pack_int(message_id, buffer_id, 1, &level_p->numrecv[isub]);

        if(level_p->numrecv[isub]>0){

          pcr_buffer_pack_int(message_id, buffer_id,
                  level_p->numrecv[isub], level_p->dofentrecv[isub]);
          pcr_buffer_pack_int(message_id, buffer_id,
                  level_p->numrecv[isub], level_p->nrdofrecv[isub]);

    /*kbw
          printf("Processor (subdomain) %d, level %d\n",my_proc_id, Level_id);
          printf("sending to processor %d: %d elements (DOF entities)\n",
             isub+1, level_p->numrecv[isub] );
          for(istr=0;istr<level_p->numrecv[isub];istr++){
        printf("vertex(node) ID %d , posglob %d, nrdof %d\n",
               level_p->dofentrecv[isub][istr],
               level_p->posrecv[isub][istr],
               level_p->nrdofrecv[isub][istr]);
          }
    /*kew*/

        }

        pcr_buffer_send(message_id, buffer_id, isub+1);
    }

  }

  for(iproc=0;iproc<nr_proc;iproc++){

  // omitting sendig to self
  if(iproc != appv_my_proc_id-1) {


    //buffer_id = pcr_buffer_receive(message_id, iproc, 0);
    // iproc version more robust, PCC_ANY_PROC may require the use of proper message ID
    buffer_id = pcr_buffer_receive(message_id, PCC_ANY_PROC, 0);

#ifdef DEBUG_PCM

    printf("Asynchronous messaging (PCC_ANY_PROC) - may lead to incorrect assignement of messages\n");

#endif

    pcr_buffer_unpack_int(message_id, buffer_id, 1, &i);
    isub=i-1;
    pcr_buffer_unpack_int(message_id, buffer_id, 1, &level_p->numsend[isub]);

    if(level_p->numsend[isub]>0){

      level_p->dofentsend[isub] = (int *)malloc(level_p->numsend[isub]*sizeof(int));
      level_p->possend[isub] = (int *)malloc(level_p->numsend[isub]*sizeof(int));
      level_p->nrdofsend[isub] = (int *)malloc(level_p->numsend[isub]*sizeof(int));

      pcr_buffer_unpack_int(message_id, buffer_id, 
			    level_p->numsend[isub], level_p->dofentsend[isub]);
      pcr_buffer_unpack_int(message_id, buffer_id, 
			    level_p->numsend[isub], level_p->nrdofsend[isub]);

      for(istr=0;istr<level_p->numsend[isub];istr++){
	struct_id = L_vert_to_struct[level_p->dofentsend[isub][istr]];
	level_p->possend[isub][istr]=L_struct_posg[struct_id];
      }
	

    }

    pcr_recv_buffer_close(message_id, buffer_id);

/*kbw
      printf("Processor (subdomain) %d, level %d, %d consecutive internal dofs\n",
	     my_proc_id, Level_id, level_p->nrdof_int);
      printf("received from  processor %d: %d elements (DOF entities)\n",
	     isub+1, level_p->numsend[isub] );
      for(istr=0;istr<level_p->numsend[isub];istr++){
	printf("vertex(node) ID %d, posglob %d, nrdof %d\n",
	       level_p->dofentsend[isub][istr],
	       level_p->possend[isub][istr],
	       level_p->nrdofsend[isub][istr]);
      }
/*kew*/
  }

  }

/*kbw
printf("Processor (subdomain) %d leaving create_exchange_tables\n", my_proc_id);
/*kew*/
  utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM);
  return(0);

}

/*---------------------------------------------------------
  appr_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
int appr_exchange_dofs(
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /* in: level ID */
  double* Vec_dofs  /* in: vector of dofs to be exchanged */
  )
{

  appt_levels* level_p; /* pointer to level structure */

  int i, isub, ibl, iproc, buffer_id, message_id; 

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  printf(
    "Processor (subdomain) %d entering exchange_dofs for table %d, level %d\n", 
    appv_my_proc_id, Table_index, Level_id);
/*kew*/

  /* proceed only if necessary */
  if(appv_nr_proc==1) return(0);

  utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );

#ifdef DEBUG
  if(Level_id >= appv_exchange_table[Table_index].nr_levels){
    printf("Wrong Level_id %d in appr_exchange_dofs! Exiting!\n", Level_id);
    exit(-1);
  }
#endif

  appv_cur_exchange_table_id = Table_index;
  /* pointer to particular level structure */
  level_p = &appv_exchange_table[Table_index].level[Level_id];

/*kbw
  printf("Processor (subdomain) %d entering exchange_dofs\n", appv_my_proc_id);
  printf("Input vector, %d internal DOFs:\n", level_p->nrdof_int);
  for(i=0;i<level_p->nrdof_int;i++){
  printf("%20.15lf",Vec_dofs[i]);
  }
  //  for(i=0;i<12;i++){
//  printf("%20.15lf",Vec_dofs[i]);
//  }
//  printf("\n");
//  for(i=1920;i<1920+12;i++){
//  printf("%20.15lf",Vec_dofs[i]);
//  }
//  printf("\n");
//  for(i=level_p->nrdof_int;i<level_p->nrdof_int+12;i++){
//  printf("%20.15lf",Vec_dofs[i]);
//  }
  printf("\n");
/*kew*/

  //  message_id = PCC_EX_DOFS_1;
  message_id = 3134+Level_id;

  for(isub=0;isub<appv_nr_proc;isub++){
    
    if(isub != appv_my_proc_id-1) {
        if(level_p->numsend[isub]>0){

          buffer_id = pcr_send_buffer_open(message_id, 0);

          pcr_buffer_pack_int(message_id, buffer_id, 1, &appv_my_proc_id);

          for(ibl=0;ibl<level_p->numsend[isub];ibl++){
        pcr_buffer_pack_double(message_id, buffer_id,
                       level_p->nrdofsend[isub][ibl],
                       &Vec_dofs[level_p->possend[isub][ibl]]);
          }

          pcr_buffer_send(message_id, buffer_id, isub+1);

    /*kbw
          printf("Processor (subdomain) %d, %d consecutive internal dofs\n",
             appv_my_proc_id, level_p->nrdof_int);
          printf("sending to  processor %d: %d elements (DOF entities)\n",
             isub+1, level_p->numsend[isub] );
    /*kbw
          for(ibl=0;ibl<level_p->numsend[isub];ibl++){
        printf("element ID %d, posglob %d, nrdof %d, dofs:\n",
               level_p->dofentsend[isub][ibl],level_p->possend[isub][ibl],level_p->nrdofsend[isub][ibl]);
        for(i=level_p->possend[isub][ibl];i<level_p->possend[isub][ibl]+level_p->nrdofsend[isub][ibl];i++){
          printf("%20.15lf",Vec_dofs[i]);
        }
        printf("\n");
          }
    /*kew*/

        } /* end if sending dofs to a given subdoamin */
    }

  } /* end loop over subdomains */

  for(iproc=0;iproc<appv_nr_proc;iproc++){

    if(iproc != appv_my_proc_id-1) {
        if(level_p->numrecv[iproc]>0){

	  //buffer_id = pcr_buffer_receive(message_id, iproc, 0);
	  // iproc version more robust, PCC_ANY_PROC may require the use of proper message ID
	  buffer_id = pcr_buffer_receive(message_id, PCC_ANY_PROC, 0);
	  
#ifdef DEBUG_PCM
	  
	  printf("Asynchronous messaging (PCC_ANY_PROC) - may lead to incorrect assignement of messages\n");
	  
#endif
	  
          pcr_buffer_unpack_int(message_id, buffer_id, 1, &i);
          isub=i-1;
	  
          for(ibl=0;ibl<level_p->numrecv[isub];ibl++){
	    
	    pcr_buffer_unpack_double(message_id, buffer_id,
				     level_p->nrdofrecv[isub][ibl],
				     &Vec_dofs[level_p->posrecv[isub][ibl]]);
	    
          }
	  
          pcr_recv_buffer_close(message_id, buffer_id);

    /*kbw
          printf("Processor (subdomain) %d\n",appv_my_proc_id);
          printf("received from processor %d: %d elements (DOF entities)\n",
             isub+1, level_p->numrecv[isub] );
    /*kbw
          for(ibl=0;ibl<level_p->numrecv[isub];ibl++){
        printf("element ID %d, posglob %d, nrdof %d\n",
               level_p->dofentrecv[isub][ibl],level_p->posrecv[isub][ibl],level_p->nrdofrecv[isub][ibl]);
        for(i=level_p->posrecv[isub][ibl];i<level_p->posrecv[isub][ibl]+level_p->nrdofrecv[isub][ibl];i++){
          printf("%20.15lf",Vec_dofs[i]);
        }
        printf("\n");
          }
    /*kew*/

        } /* end if receiving dofs */
    }
  } /* end loop over processors */

/*kbw
printf("Processor (subdomain) %d leaving exchange_dofs\n", appv_my_proc_id);
/*kew*/

  utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM );
  return(1);

}

/*---------------------------------------------------------
  appr_sol_vec_norm - to compute a norm of global vector in parallel
---------------------------------------------------------*/
double appr_sol_vec_norm( /* returns: L2 norm of global Vector */
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /* in: level ID */
  int Nrdof,            /* in: number of vector components */
  double* Vector        /* in: local part of global Vector */
  )
{

  appt_levels* level_p; /* pointer to level structure */

  int IONE=1;

  double daux, vec_norm;
  int buffer_id, message_id, iproc;

/*++++++++++++++++ executable statements ++++++++++++++++*/
 utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );
/*kbw
printf("Processor (subdomain) %d entering sol_vec_norm\n", appv_my_proc_id);
printf("Table_index %d, Level_id %d, Nrdof %d\n", Table_index, Level_id, Nrdof);
/*kew*/

  /* for sequential execution */
  if(appv_nr_proc==1) vec_norm = dnrm2_(&Nrdof, Vector, &IONE);
  else{

    appv_cur_exchange_table_id = Table_index;

#ifdef DEBUG
    if(Level_id >= appv_exchange_table[Table_index].nr_levels){
      printf("Wrong Level_id %d in appr_sol_vec_norm! Exiting!\n", Level_id);
      exit(-1);
    }
#endif

    /* pointer to particular level structure */
    level_p = &appv_exchange_table[Table_index].level[Level_id];

/* 
 * if dofs are ordered so that initial nrdof_int dofs correspond 
 * to internal nodes we can use BLAS 
 */
    if(level_p->nrdof_int>0){

       daux = ddot_(&level_p->nrdof_int, Vector, &IONE, Vector, &IONE);

/*kbw
      printf("Processor (subdomain) %d\n",appv_my_proc_id);
      printf("internal dofs %d, local vector norm %lf\n",
	     level_p->nrdof_int, daux);
/*kew*/

    }
    else{

      printf("Modify appr_vec_norm for not consecutive DOFs!\n");
      exit(-1);

    }

    pcr_allreduce_sum_double(1, &daux, &vec_norm);

/*kbw
      printf("Processor (subdomain) %d\n",appv_my_proc_id);
      printf("computed global vector norm %lf\n", vec_norm );
/*kew*/

    vec_norm=sqrt(vec_norm);

  }

/*kbw
      printf("Nr_proc %d, my_proc_id %d, computed total norm %lf\n",
	     appv_nr_proc, appv_my_proc_id, vec_norm);
/*kew*/

/*kbw
printf("Processor (subdomain) %d leaving sol_vec_norm\n", appv_my_proc_id);
/*kew*/
  utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM );
  return(vec_norm);
}


/*---------------------------------------------------------
  appr_sol_sc_prod - to compute a scalar product of two global vectors 
---------------------------------------------------------*/
double appr_sol_sc_prod( /* retruns: scalar product of Vector1 and Vector2 */
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /* in: level ID */
  int Nrdof,           /* in: number of vector components */
  double* Vector1,     /* in: local part of global Vector */
  double* Vector2      /* in: local part of global Vector */
  )
{

  int IONE=1;

  appt_levels* level_p; /* pointer to level structure */

  double daux, sc_prod;
  int buffer_id, iproc;

/*++++++++++++++++ executable statements ++++++++++++++++*/
 utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );
/*kbw
printf("Processor (subdomain) %d entering sol_sc_prod\n", appv_my_proc_id);
/*kew*/

  /* for sequential execution */
  if(appv_nr_proc==1) sc_prod = ddot_(&Nrdof, Vector1, &IONE, Vector2, &IONE);
  else{
    
#ifdef DEBUG
    if(Level_id >= appv_exchange_table[Table_index].nr_levels){
      printf("Wrong Level_id %d in appr_sol_vec_norm! Exiting!\n", Level_id);
      exit(-1);
    }
#endif

    appv_cur_exchange_table_id = Table_index;
    /* pointer to particular level structure */
    level_p = &appv_exchange_table[Table_index].level[Level_id];
    
    /* 
     * if dofs are ordered so that initial nrdof_int dofs correspond 
     * to internal nodes we can use BLAS 
     */
    if(level_p->nrdof_int>0){
      
      daux = ddot_(&level_p->nrdof_int, Vector1, &IONE, Vector2, &IONE);
      
    }
    else{
      
      printf("Internal DOFs not at the beginning in sol_sc_prod. Exiting!\n");
      exit(0);

    }
        
/*kbw
  printf("Subdomain %d: computed local scalar product %lf for %d DOFS\n",
	 appv_my_proc_id, daux, level_p->nrdof_int);
/*kew*/

    pcr_allreduce_sum_double(1, &daux, &sc_prod);
    
  }
  
/*kbw
      printf("Nr_proc %d, my_proc_id %d, computed total scalar product %lf\n",
	     appv_nr_proc, appv_my_proc_id, sc_prod);
/*kew*/
  
/*kbw
printf("Processor (subdomain) %d leaving sol_sc_prod\n", appv_my_proc_id);
/*kew*/
  utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM );
  return(sc_prod);
}


/**--------------------------------------------------------
  appr_get_nr_dofs_owned - utility to calculate the number of dofs owned at
                          the finest level
---------------------------------------------------------*/
extern int appr_get_nr_dofs_owned( 
  int Table_index    /* in: exchange table index  */
  )
{
    int rv=0;
    utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );
    appv_cur_exchange_table_id = Table_index;
    //mf_log_info("Table_index=%d",Table_index);
    if(pcv_nr_proc > 1) {
      int level_id = appv_exchange_table[Table_index].nr_levels-1;

      appt_levels* level_p = &appv_exchange_table[Table_index].level[level_id];
      //mf_log_info("Table_index=%d, level_id=%d, level_p->nrdof_int=%d",
      //           Table_index,level_id,level_p->nrdof_int);

      rv = (level_p->nrdof_int);
    }
    else {
        rv = apr_get_nrdofs_glob(appv_exchange_table[Table_index].field_id);
    }

    utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM );
    return rv;

}

/*---------------------------------------------------------
  appr_free_exchange_tables - to free data structures related to exchange 
                             of dofs
---------------------------------------------------------*/
int appr_free_exchange_tables( 
		  /* returns: >=0 -success code, <0 -error code */
  int Table_index    /* in: exchange table index  */
  )
{

  int ilev, isub, nr_sub;
  appt_levels* level_p; /* pointer to level structure */

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* proceed only if necessary */
  if(appv_nr_proc==1) return(0);
  
  utr_io_result_cumulative_timer_start(RESULT_TIME_APPM );

  nr_sub = appv_nr_proc;
  appv_cur_exchange_table_id = Table_index;
  
  /* free data structures to manage exchange of dofs */
  for(ilev=0; ilev<appv_exchange_table[Table_index].nr_levels;ilev++){
    
    level_p = &appv_exchange_table[Table_index].level[ilev];
    
    for(isub=0;isub<nr_sub;isub++){
      
      if(level_p->numsend[isub]>0){
	free(level_p->dofentsend[isub]);
	free(level_p->possend[isub]);
	free(level_p->nrdofsend[isub]);
      }
      
      if(level_p->numrecv[isub]>0){
	free(level_p->dofentrecv[isub]);
	free(level_p->posrecv[isub]);
	free(level_p->nrdofrecv[isub]);
      }
      
    }
    
  }
  
  free(appv_exchange_table[Table_index].level);
    
  utr_io_result_cumulative_timer_stop(RESULT_TIME_APPM );

  return(0);
}



