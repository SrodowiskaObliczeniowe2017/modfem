/*************************************************************
File contains procedures:
 lar_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
  lar_initialize_SM_and_LV - to initialize stiffness matrix and/or load vector
  lar_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local
                           stiffness matrix and load vector
  lar_free_SM_and_LV - to free space for a block structure

------------------------------
History:
	11.2015 - Damian Goik
*************************************************************************/

#include "uth_log.h"
#include "las_petsc_intf.hpp"



struct AMGSolverData amg_solver_data;
int first_ghost_block = -1;
int posglob_offset = 0;

void amg_project_solution_to_level(int level_id)
{
	amg_solver_data.solver->MoveSolutionToLevel(level_id);
}

void amg_project_solution_from_level(int level_id)
{
	amg_solver_data.solver->MoveSolutionFromLevel(level_id);
}

double amg_get_global_diff(double diff){
  double global_diff;
  MPI_Allreduce(&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
  return global_diff;
}

/*---------------------------------------------------------
  lar_petsc_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_petsc_allocate_SM_and_LV( // returns: matrix index in itv_matrices array
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int Nrdof_glob,  /* in: total number of DOFs */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
	)
{
	  amg_solver_data.posglobs = (int*)malloc(sizeof(int) * Nrblocks);
	  memcpy(amg_solver_data.posglobs, Posglob,sizeof(int)* Nrblocks);
//	  for(int i = 0; i<10; i++)
//		  printf("%d\n", Posglob[i]);
	  int block_index;
	  int max_nonzeros_in_row_diagonal = 0;
	  int max_nonzeros_in_row_off_diagonal = 0;
	  int ghost_block_cnt = 0;
	  //TODO: make a better preallocation
	  for(block_index=0;block_index<Nrblocks;block_index++)
	  {
	      if(Nroffbl[block_index]<=0)
	      {
	    	  if(first_ghost_block < 0){
	    		  first_ghost_block = amg_solver_data.posglobs[block_index];
	    		  printf("FIRST GHOST ASSIGNED %d %d\n",first_ghost_block, amg_solver_data.posglobs[block_index]);
	    	  }
	    	  ghost_block_cnt++;
	      }
	  }

	  for(block_index=0;block_index<Nrblocks;block_index++)
	  {
	      int diagonal_nonzeros = Nrdofbl[block_index];
	      int off_diagonal_nonzeros = 0;
	      if(Nroffbl[block_index]>0)
	      {
	    	  int i;
	    	  for(i=0;i<Nroffbl[block_index];i++)
	    	  {
	    		  int intersecting_block_id = L_offbl[block_index][i];
	    		  if(amg_solver_data.posglobs[intersecting_block_id] >= first_ghost_block && first_ghost_block != -1){
	    			  off_diagonal_nonzeros += Nrdofbl[intersecting_block_id];
	    		  }
	    		  else{
					  diagonal_nonzeros += Nrdofbl[intersecting_block_id];
	    		  }
	    	  }
	      }
	      max_nonzeros_in_row_diagonal = std::max(max_nonzeros_in_row_diagonal, diagonal_nonzeros);
	      max_nonzeros_in_row_off_diagonal = std::max(max_nonzeros_in_row_off_diagonal, off_diagonal_nonzeros);
	  }

	  printf("NR of ghost blocks %d %d %d\n",ghost_block_cnt, max_nonzeros_in_row_diagonal, max_nonzeros_in_row_off_diagonal);
	  Nrdof_glob -= ghost_block_cnt;
	  int global_nrdof_glob;
	  int local_nrdof = Nrdof_glob;

	  PetscBool initialized;
	  PetscInitialized(&initialized);
	  if(!initialized)
		  PetscInitialize(NULL,NULL,NULL,NULL);

	  MPI_Allreduce(&Nrdof_glob, &global_nrdof_glob, 1,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  Nrdof_glob = global_nrdof_glob;
	  amg_solver_data.amgInitData = amg_algorithm_data;

	  // ERROR in 3.7.4 - need 2 arguments - KB
	  //	  PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);
	  PetscErrorCode ierr;
	  ierr = MatCreate(PETSC_COMM_WORLD,&(amg_solver_data.A));CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetSizes(amg_solver_data.A,local_nrdof,local_nrdof,Nrdof_glob,Nrdof_glob);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetFromOptions(amg_solver_data.A);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  //MatSetOption(amg_solver_data.A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE);

	  ierr = VecCreate(PETSC_COMM_WORLD,&(amg_solver_data.b));CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = VecSetSizes(amg_solver_data.b,local_nrdof,Nrdof_glob);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = VecSetFromOptions(amg_solver_data.b);CHKERRABORT(PETSC_COMM_WORLD, ierr);



	  ierr = MatMPIAIJSetPreallocation(amg_solver_data.A,max_nonzeros_in_row_diagonal,NULL,max_nonzeros_in_row_off_diagonal,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSeqAIJSetPreallocation(amg_solver_data.A,max_nonzeros_in_row_diagonal,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetUp(amg_solver_data.A); CHKERRABORT(PETSC_COMM_WORLD, ierr);

	  printf("CURR MATRIX ID\n");

		int size;
		MPI_Comm_size(PETSC_COMM_WORLD, &size);
		const PetscInt* ranges = new PetscInt[size+1];
		MatGetOwnershipRanges(amg_solver_data.A,&ranges);
		int rank;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		posglob_offset = ranges[rank];
		init_posglobs(amg_solver_data.posglobs, posglob_offset);
		for(int i = 0; i<Nrblocks; i++){
			if(Nroffbl[i] > 0)
				amg_solver_data.posglobs[i] += posglob_offset;
			else
				amg_solver_data.posglobs[i] = -1;
		}
		printf("FIRST POSGLOB %d \n", posglob_offset);
	  return 0;
	  //return(itv_cur_matrix_id);
}



int lar_petsc_initialize_SM_and_LV(int Matrix_id, int Comp_type)
{
	//for the petsc no need to split allocation and initialization, though part of the allocation call may be places here
	MatZeroEntries(amg_solver_data.A);
	VecSet(amg_solver_data.b, 0.0);
	delete amg_solver_data.solver;
	amg_solver_data.solver = NULL;
	return 0;
}

//#include "mmh_intf.h"
void writeNodeCoordinates()
{
//	int el_nodes[MMC_MAXELVNO+1];
//	double coor[MMC_MAXELVNO*3];
//	int elem_id = mmr_get_next_act_elem(1,0);
//	//int ent_nr = mmr_get_nr_elem(1);
//	FILE* f = fopen("mesh.dat","w+");
//	while(true)
//	{
//		if(elem_id <= 0)
//			break;
//		mmr_el_node_coor(1, elem_id, el_nodes, coor);
//		for(int j = 0; j<el_nodes[0]; j++)
//		{
//			fprintf(f,"%d %lf %lf %lf\n",amg_solver_data.posglobs[el_nodes[j+1]],coor[3*j],coor[3*j+1],coor[3*j+2]);
//		}
//		elem_id = mmr_get_next_act_elem(1,elem_id);
//	}
//	fclose(f);
}

/*------------------------------------------------------------
  lar_petsc_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local
                           stiffness matrix and load vector
------------------------------------------------------------*/
int lar_petsc_assemble_SM_and_LV(
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
                         /*   ITC_SOLVE - solve the system */
                         /*   ITC_RESOLVE - resolve for the new rhs vector */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
  )
{
//	printf("%d\n",first_ghost_block);

	PetscErrorCode ierr;
	int i;
	int j;
	int all_blocks_height = 0;
	for(i = 0; i<Nr_dof_bl; i++)
		all_blocks_height += L_bl_nrdof[i];

	int posglob_column_buffer[1000];

	if(all_blocks_height > 1000){
		printf("POSGLOB BUFFER MAX IS 1000");
		CHKERRABORT(PETSC_COMM_WORLD, 17);
	}

	int column_index = 0;
	for(i = 0; i<Nr_dof_bl; i++){
		int posglob = amg_solver_data.posglobs[L_bl_id[i]];
		for(j = 0; j<L_bl_nrdof[i]; j++){
			posglob_column_buffer[column_index++] = posglob + j;
		}
	}

	//transpose
	for(i = 0; i<all_blocks_height; i++){
		for(j = i+1; j<all_blocks_height; j++){
			double tmp = Stiff_mat[i*all_blocks_height + j];
			Stiff_mat[i*all_blocks_height + j] = Stiff_mat[j*all_blocks_height + i];
			Stiff_mat[j*all_blocks_height + i] = tmp;
		}
	}

	int row_index = 0;
	for(i = 0; i<Nr_dof_bl; i++){
//		int posglob = amg_solver_data.posglobs[L_bl_id[i]];
		if(first_ghost_block == -1 || (amg_solver_data.posglobs[L_bl_id[i]]< first_ghost_block + posglob_offset
				&& amg_solver_data.posglobs[L_bl_id[i]] >= posglob_offset)){

			ierr = MatSetValues(amg_solver_data.A,L_bl_nrdof[i],posglob_column_buffer + row_index,
					column_index,posglob_column_buffer,Stiff_mat + all_blocks_height*row_index,ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);

			ierr = VecSetValues(amg_solver_data.b,L_bl_nrdof[i],posglob_column_buffer + row_index,
					&(Rhs_vect[row_index]),ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);
		}

		row_index += L_bl_nrdof[i];
	}



//	int row_offset = 0;
//	for(i = 0; i<Nr_dof_bl; i++)
//	{
//		int column_offset = 0;
//		for(j = 0; j<Nr_dof_bl; j++)
//		{
//
//			//printf("Przeciecie bloku %d z %d\n",L_bl_id[i], L_bl_id[j]);
//			int k;
//			int l;
//
//			int j_index = amg_solver_data.posglobs[L_bl_id[j]];
//
//			if(first_ghost_block == -1 || (amg_solver_data.posglobs[L_bl_id[i]]< first_ghost_block + posglob_offset
//					&& amg_solver_data.posglobs[L_bl_id[i]] >= posglob_offset)){
//				for(k = 0; k<L_bl_nrdof[i]; k++)
//				{
//					int i_index = amg_solver_data.posglobs[L_bl_id[i]];
//					for(l = 0; l<L_bl_nrdof[j]; l++)
//					{
//						double value = Stiff_mat[(column_offset+k)*all_blocks_height + row_offset +l];
//						//printf("(%d %d) v %lf\n",i_index,j_index,value);
//						//test[i_index][j_index] += value;
//						ierr = MatSetValues(amg_solver_data.A,1,&i_index,1,&j_index,&value,ADD_VALUES);
//						if(ierr != 0){
//							printf("XXXX %d %d\n", i_index, j_index);
//						}
//						CHKERRABORT(PETSC_COMM_WORLD, ierr);
//						i_index++;
//					}
//					j_index++;
//					//printf("\n");
//				}
//			}
//
//			column_offset += L_bl_nrdof[j];
//		}
//		int r;
//		for(r = 0; r<L_bl_nrdof[i]; r++)
//		{
////			printf("FIRST %d\n", first_ghost_block);
//			//printf("rhs %d wartosc %lf\n",r + amg_solver_data.posglobs[L_bl_id[i]],Rhs_vect[r+row_offset]);
//			int row = r + amg_solver_data.posglobs[L_bl_id[i]];
//			if(first_ghost_block == -1 || (amg_solver_data.posglobs[L_bl_id[i]]< first_ghost_block + posglob_offset
//					&& amg_solver_data.posglobs[L_bl_id[i]] >= posglob_offset)){
//				//printf("size: %d\n",size);
//	//			if(row > 0)
//	//				printf("row %d\n",row);
//
//				ierr = VecSetValues(amg_solver_data.b,1,&row,&(Rhs_vect[r+row_offset]),ADD_VALUES);
////				ierr = VecSetValues(amg_solver_data.b,1,&row,&one,ADD_VALUES);
//
//				CHKERRABORT(PETSC_COMM_WORLD, ierr);
//			}
//		}
//		row_offset += L_bl_nrdof[i];
//	}
//	printf("END OF ASSEMBLING\n");
  return(0);
}


int lar_petsc_fill_preconditioner(int Matrix_id)
{



	MPI_Barrier(PETSC_COMM_WORLD);
	//reuse this call to create AMG solver structure + levels
	PetscErrorCode ierr;
	ierr = MatAssemblyBegin(amg_solver_data.A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyEnd(amg_solver_data.A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyBegin(amg_solver_data.b); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(amg_solver_data.b); CHKERRABORT(PETSC_COMM_WORLD, ierr);

	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Optmat123parallel.data", &viewer);
	MatView(amg_solver_data.A, viewer);
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Optvec123parallel.data", &viewer);
	VecView(amg_solver_data.b, viewer);
	PetscViewerDestroy(&viewer);
	MPI_Barrier(PETSC_COMM_WORLD);
	CHKERRABORT(PETSC_COMM_WORLD, 17);

	//writeNodeCoordinates();
	int Ndof;
	ierr = MatGetSize(amg_solver_data.A,&Ndof,NULL);
	ErrorEvaluator* error_evaluator = new ResidualBasedErrorEvaluator(amg_solver_data.A, amg_solver_data.b, 1.0e-7, 1.0e-8, Ndof, -1);
	amg_solver_data.solver = new AMGSolverStructure(amg_solver_data.A, amg_solver_data.b, Ndof, error_evaluator, false);


	clock_t begin = clock();
	amg_solver_data.solver->CreateAMGSolverLevels(amg_solver_data.amgInitData.coarsening_algorithm,
			amg_solver_data.amgInitData.interpolation_algorithm,
			amg_solver_data.amgInitData.strength_threshold,
			amg_solver_data.amgInitData.levels_number);
	//amg_solver_data.solver->ExtractCoarseRowsNumbers();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Total time of level creation: %lf\n",elapsed_secs);

	return 0;
}

int lar_petsc_free_SM_and_LV(int Matrix_id)
{
	if(amg_solver_data.solver != NULL)
	{
		printf("FREE AMG \n");
		MPI_Barrier(PETSC_COMM_WORLD);
		delete amg_solver_data.solver;
		MatDestroy(&(amg_solver_data.A));
		VecDestroy(&(amg_solver_data.b));
		free(amg_solver_data.posglobs);
		amg_solver_data.solver = NULL;
	}
	return 0;
}

void lar_petsc_perform_BJ_or_GS_iterations(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* in,out: vector of unknowns updated */
		/* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{
	time_t start = time(NULL);
	clock_t begin = clock();
	amg_solver_data.solver->RunVCycle(Matrix_id, V);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Total time of solving: %lf %.2f\n",elapsed_secs, (double)(time(NULL) - start));
	/*
	row 0: (0, 2.22222e+06)  (1, 1.11111e+06)  (2, 555555)  (3, 555555)  (4, 277778)  (5, 277778)
	row 1: (0, 1.11111e+06)  (1, 2.22222e+06)  (2, 277778)  (3, 277778)  (4, 555555)  (5, 555555)
	row 2: (0, 555555)  (1, 277778)  (2, 2.22222e+06)  (3, 0.0833333)  (4, 1.11111e+06)  (5, -0.0833333)  (6, 277778)  (7, 555555)
	row 3: (0, 555555)  (1, 277778)  (2, 0.0833333)  (3, 2.22222e+06)  (4, -0.0833333)  (5, 1.11111e+06)  (6, 277778)  (7, 555555)
	row 4: (0, 277778)  (1, 555555)  (2, 1.11111e+06)  (3, -0.0833333)  (4, 2.22222e+06)  (5, 0.0833333)  (6, 555555)  (7, 277778)
	row 5: (0, 277778)  (1, 555555)  (2, -0.0833333)  (3, 1.11111e+06)  (4, 0.0833333)  (5, 2.22222e+06)  (6, 555555)  (7, 277778)
	row 6: (2, 277778)  (3, 277778)  (4, 555555)  (5, 555555)  (6, 2.22222e+06)  (7, 1.11111e+06)
	row 7: (2, 555555)  (3, 555555)  (4, 277778)  (5, 277778)  (6, 1.11111e+06)  (7, 2.22222e+06)
0.976263
-0.647741
-0.945741
-0.427846
0.266126
0.423095
0.976263
-0.945741
-0.647741
-0.427846
0.423095
0.266126


0.9762629342
-0.6477411784
-0.9457409368
-0.4278455596
0.2661264401
0.4230951309
0.9762629342
-0.9457409368
-0.6477411784
-0.4278455596
0.4230951309
0.2661264401
	*/

	//VecView(result, PETSC_VIEWER_STDOUT_WORLD);
	return;
}
