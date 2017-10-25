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

#include "GhostBlockExchanger.h"

int* my_posglobs;
int my_posglob_offset;

void init_posglobs(int* posglobs, int posglob_offset){
	my_posglobs = posglobs;
	my_posglob_offset = posglob_offset;
}
#ifdef PARALLEL
int exchange_global_row_indices(
                   /* returns: >=0 -success code, <0 -error code */
  int Nr_dof_ent,  /* in: number of DOF entities in the level */
  int* L_dof_ent_id,  /* in: list of DOF entities associated with DOF structs */
  int* L_struct_nrdof,    /* in: list of nrdofs for each dof struct */
  int* L_struct_posg,     /* in: list of positions within the global */
                      /*     vector of dofs for each dof struct */
  int* L_thing_to_struct, /*in: list of DOF structs associated with DOF entities */
  int field_id
  )
{
	int i,istr,ibl,iproc,isub;

	int appv_nr_proc=-1;
	int appv_my_proc_id=-1;

	int *posglobs = my_posglobs;
	int posglob_offset = my_posglob_offset;

  MPI_Comm_size(MPI_COMM_WORLD, &appv_nr_proc);


  /* proceed only if necessary */
  if(appv_nr_proc < 2)
	  return 0;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int numrecv[appv_nr_proc];
  for(i = 0; i<appv_nr_proc; i++)
	  numrecv[i] = 0;
  int my_proc_id = rank + 1;
  appv_my_proc_id = my_proc_id;

  printf("my proc id %d\n", my_proc_id);

  int mesh_id = apr_get_mesh_id(field_id);

  /* for all dof structures (for STD_LIN just vertices=nodes) */
  for(istr=0;istr<Nr_dof_ent;istr++){

    int nno=L_dof_ent_id[istr];

    int owner_id = mmpr_ve_owner(mesh_id, nno);
    if(owner_id != my_proc_id){

      /* set info on the dof structure */
      int no_id_local = nno;
      int no_id_owner = mmpr_ve_id_at_owner(mesh_id, nno);

      numrecv[owner_id-1]++;
//      printf("foreign struct my id %d owner id %d %d %d\n", no_id_local, no_id_owner, owner_id, my_proc_id);

    }
  }

  int** dofentrecv = (int **)malloc(appv_nr_proc*sizeof(int*));
  int** posrecv = (int **)malloc(appv_nr_proc*sizeof(int*));
  int** nrdofrecv = (int **)malloc(appv_nr_proc*sizeof(int*));

  for(isub=0;isub<appv_nr_proc;isub++){
    if(numrecv[isub]>0){
      printf("receiving %d structs\n", numrecv[isub]);
      dofentrecv[isub] = (int *)calloc(numrecv[isub],sizeof(int));
      posrecv[isub] = (int *)calloc(numrecv[isub],sizeof(int));
      nrdofrecv[isub] = (int *)calloc(numrecv[isub],sizeof(int));
    }
  }

  /* again initialize */
  for(isub=0;isub<appv_nr_proc;isub++){
    numrecv[isub]=0;
  }

  /* for all dof structures (for STD_LIN just vertices=nodes) */
  int istr_int=-1;

  for(istr=0;istr<Nr_dof_ent;istr++){

    int nno=L_dof_ent_id[istr];

    int owner_id = mmpr_ve_owner(mesh_id, nno);
    if(owner_id != my_proc_id){

      /* set info on the dof structure */
      int no_id_local = nno;
      int no_id_owner = mmpr_ve_id_at_owner(mesh_id, nno);

/*kbw
      printf("struct %d, element %d, ipid %d, owner %d, id_owner %d, nrdof %d, posg %d\n",
	     istr, nno, ipid, owner_id, no_id_owner,
	     L_struct_nrdof[istr], L_struct_posg[istr]);
/*kew*/

      int struct_id = numrecv[owner_id-1];
      dofentrecv[owner_id-1][struct_id]=no_id_owner;
      posrecv[owner_id-1][struct_id]=L_struct_posg[istr];
      nrdofrecv[owner_id-1][struct_id]=L_struct_nrdof[istr];
      numrecv[owner_id-1]++;

    }
    else{
//    	printf("ERROR!\n");
      /* check whether internal dofs are at the beginning */
//      istr_int++;
//      if(istr_int==istr) level_p->nrdof_int+=L_struct_nrdof[istr];
//      else level_p->nrdof_int = 0;

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
  int message_id = 1234+60;

  for(isub=0;isub<appv_nr_proc;isub++){

    // omitting sendig to self
    if(isub != my_proc_id-1) {

        int buffer_id = pcr_send_buffer_open(message_id,0);
        pcr_buffer_pack_int(message_id, buffer_id, 1, &my_proc_id);
        pcr_buffer_pack_int(message_id, buffer_id, 1, &numrecv[isub]);

        if(numrecv[isub]>0){

          pcr_buffer_pack_int(message_id, buffer_id,
                  numrecv[isub], dofentrecv[isub]);
          pcr_buffer_pack_int(message_id, buffer_id,
                  numrecv[isub], nrdofrecv[isub]);

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

  int numsend[appv_nr_proc];
  int** dofentsend = (int** )malloc(appv_nr_proc*sizeof(int*));
  int** possend = (int** )malloc(appv_nr_proc*sizeof(int*));
  int** nrdofsend = (int** )malloc(appv_nr_proc*sizeof(int*));

  for(iproc=0;iproc<appv_nr_proc;iproc++){

  // omitting sendig to self
  if(iproc != appv_my_proc_id-1) {


    //buffer_id = pcr_buffer_receive(message_id, iproc, 0);
    // iproc version more robust, PCC_ANY_PROC may require the use of proper message ID
    int buffer_id = pcr_buffer_receive(message_id, PCC_ANY_PROC, 0);

#ifdef DEBUG_PCM

    printf("Asynchronous messaging (PCC_ANY_PROC) - may lead to incorrect assignement of messages\n");

#endif

    pcr_buffer_unpack_int(message_id, buffer_id, 1, &i);
    int isub=i-1;
    pcr_buffer_unpack_int(message_id, buffer_id, 1, &(numsend[isub]));

    if(numsend[isub]>0){

      dofentsend[isub] = (int *)malloc(numsend[isub]*sizeof(int));
      possend[isub] = (int *)malloc(numsend[isub]*sizeof(int));
      nrdofsend[isub] = (int *)malloc(numsend[isub]*sizeof(int));

      pcr_buffer_unpack_int(message_id, buffer_id,
			    numsend[isub], dofentsend[isub]);
      pcr_buffer_unpack_int(message_id, buffer_id,
			    numsend[isub], nrdofsend[isub]);

      for(istr=0;istr<numsend[isub];istr++){
    	  int struct_id = L_thing_to_struct[dofentsend[isub][istr]];
    	  possend[isub][istr]=L_struct_posg[struct_id] + posglob_offset;
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




  //  message_id = PCC_EX_DOFS_1;
  message_id = 3134+60;

  for(isub=0;isub<appv_nr_proc;isub++){

    if(isub != appv_my_proc_id-1) {
        if(numsend[isub]>0){

          int buffer_id = pcr_send_buffer_open(message_id, 0);

          pcr_buffer_pack_int(message_id, buffer_id, 1, &appv_my_proc_id);

          for(ibl=0;ibl<numsend[isub];ibl++){
        	  pcr_buffer_pack_int(message_id, buffer_id, 1, &(possend[isub][ibl]));
//                       nrdofsend[isub][ibl],
//                       &Vec_dofs[level_p->possend[isub][ibl]]);
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
        if(numrecv[iproc]>0){

	  //buffer_id = pcr_buffer_receive(message_id, iproc, 0);
	  // iproc version more robust, PCC_ANY_PROC may require the use of proper message ID
	  int buffer_id = pcr_buffer_receive(message_id, PCC_ANY_PROC, 0);

#ifdef DEBUG_PCM

	  printf("Asynchronous messaging (PCC_ANY_PROC) - may lead to incorrect assignement of messages\n");

#endif

          pcr_buffer_unpack_int(message_id, buffer_id, 1, &i);
          int isub=i-1;

          for(ibl=0;ibl<numrecv[isub];ibl++){

        	  pcr_buffer_unpack_int(message_id, buffer_id, 1, &posglobs[posrecv[isub][ibl]]);
//        	  printf("Assigning into position %d glob %d\n", posrecv[isub][ibl], posglobs[posrecv[isub][ibl]]);

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



  return(0);

}
#endif
