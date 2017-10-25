#include "AMGSolverStructure.hpp"

//FILE* utv_log_out = stdout;

AMGSolverStructure::AMGSolverStructure(Mat matrix, Vec rhs, int size, ErrorEvaluator* error_evaluator,
		bool create_coarsening_scheme)
{
	this->matrix = matrix;
	this->rhs = rhs;
	this->size = size;
	this->error_evaluator = error_evaluator;
	this->create_coarsening_scheme = create_coarsening_scheme;
	this->smoothers = NULL;
}

AMGSolverStructure::~AMGSolverStructure()
{
	EndVCycle();
	for(int i = 1; i<levelsNumber; i++){
		MatDestroy(&amg_level_matrices[i]);
		MatDestroy(&amg_interpolation_matrices[i-1]);
		VecDestroy(&amg_level_rhs[i]);
		VecDestroy(&amg_level_x[i]);
		VecDestroy(&amg_level_residuals[i-1]);
	}
	delete[] amg_level_matrices;
	delete[] amg_interpolation_matrices;
	delete[] amg_level_rhs;
	delete[] amg_level_x;
	delete[] amg_level_residuals;
	delete[] amg_level_sizes;

	for(int i = 0; i<levelsNumber - 1; i++)
		delete cf_splitters[i];
	delete[] cf_splitters;

	if(create_coarsening_scheme)
		for(int i = 0; i<levelsNumber - 1; i++)
			delete coarseRowsNumbers[i];
	delete[] coarseRowsNumbers;

	delete error_evaluator;
}

RSCFSplitter* AMGSolverStructure::GetSplitter(int levelNumber, double strength_threshold, int coarsening_algorithm, int interpolation_algorithm){

	InterpolationStrategy* strategy;
	if(interpolation_algorithm == 1){
		strategy = new GeneralInterpolation();
	}
	else{
		strategy = new DirectInterpolation();
	}

	if(coarsening_algorithm == 1){
		 return new GeneralRSCFSplitter(amg_level_matrices[levelNumber-1],strength_threshold, strategy);
	}
	else{
		return new RSCFSplitter(amg_level_matrices[levelNumber-1],strength_threshold, strategy);
	}
}

void AMGSolverStructure::CreateAMGSolverLevels(int coarsening_algorithm, int interpolation_algorithm, double strength_threshold, int levels_number)
{
	levelsNumber = levels_number;
	InitLevels();
	printf("level size %d\n",amg_level_sizes[0]);
	double matrix_multiplication_time = 0;
//	MatView(amg_level_matrices[0], PETSC_VIEWER_STDOUT_WORLD);
	for(int levelNumber = 1; levelNumber<GetLevelsNumber(); levelNumber++)
	{
		cf_splitters[levelNumber - 1] = GetSplitter(levelNumber, strength_threshold, coarsening_algorithm, interpolation_algorithm);
		cf_splitters[levelNumber - 1]->MakeCFSplitting();
		if(create_coarsening_scheme)
			coarseRowsNumbers[levelNumber - 1] = cf_splitters[levelNumber - 1]->GetNumbersOfCoarseRows();
		amg_interpolation_matrices[levelNumber - 1] = cf_splitters[levelNumber - 1]->GetMatrixFromCoarseToFine();
		clock_t begin = clock();
		MatGetSize(amg_interpolation_matrices[levelNumber - 1],NULL,&(amg_level_sizes[levelNumber]));
		Mat intermediate;

		PetscErrorCode ierr = MatMatMult(amg_level_matrices[levelNumber - 1],amg_interpolation_matrices[levelNumber - 1],
				MAT_INITIAL_MATRIX,PETSC_DEFAULT,&intermediate);
		CHKERRABORT(PETSC_COMM_WORLD, ierr);
		MatTransposeMatMult(amg_interpolation_matrices[levelNumber - 1],intermediate,
				MAT_INITIAL_MATRIX,PETSC_DEFAULT,&(amg_level_matrices[levelNumber]));
//		MatView(amg_level_matrices[levelNumber], PETSC_VIEWER_STDOUT_WORLD);
//		MatView(amg_interpolation_matrices[levelNumber - 1], PETSC_VIEWER_STDOUT_WORLD);
		MatInfo info;
		MatGetInfo(amg_level_matrices[levelNumber],MAT_GLOBAL_SUM,&info);
		printf("%lf\n", info.nz_used/ amg_level_sizes[levelNumber] / amg_level_sizes[levelNumber]);
	    VecCreate(PETSC_COMM_WORLD,&(amg_level_x[levelNumber]));
	    int localRowNumbers;
	    MatGetLocalSize(amg_level_matrices[levelNumber], &localRowNumbers, NULL);
		VecSetSizes(amg_level_x[levelNumber],localRowNumbers,amg_level_sizes[levelNumber]);
		VecSetFromOptions(amg_level_x[levelNumber]);
		VecSet(amg_level_x[levelNumber],0.0);
		VecDuplicate(amg_level_x[levelNumber-1], &(amg_level_residuals[levelNumber-1]));
		VecDuplicate(amg_level_x[levelNumber], &(amg_level_rhs[levelNumber]));
		clock_t end = clock();
		matrix_multiplication_time += double(end - begin) / CLOCKS_PER_SEC;
		printf("level size %d\n",amg_level_sizes[levelNumber]);
		//printf("AFTER MULTIPLICATION\n");
	}
	//printf("MATRIX MULTIPLICATION TIME: %lf\n", matrix_multiplication_time);
}

void AMGSolverStructure::ExtractCoarseRowsNumbers()
{

//	for(int i = 0 ; i<levelsNumber-1; i++)
//	{
//		printf("%d\n",i);
//		for (std::list<int>::iterator it=coarseRowsNumbers[i]->begin() ; it != coarseRowsNumbers[i]->end(); ++it)
//		{
//			printf("%d ", *it);
//		}
//		printf("\n");
//	}

	for(int i = 1; i<levelsNumber-1; i++)
	{
		std::list<int>* fine = coarseRowsNumbers[i-1];
		std::list<int>* coarse = coarseRowsNumbers[i];

	    std::list<int>::iterator fine_it = fine->begin();
	    int previous_advance = 0;
	    int advance = 0;
		for (std::list<int>::iterator it=coarse->begin() ; it != coarse->end();)
		{
			previous_advance = advance;
			advance = *it;
			std::advance(fine_it, advance - previous_advance);
			int coraseRowNumber = *fine_it;
			it = coarse->erase(it);
			coarse->insert(it, coraseRowNumber);
		}
	}
	FILE* f = fopen("mesh.dat","a+");
	for(int i = 0 ; i<levelsNumber-1; i++)
	{
		fprintf(f,"LEVEL\n");
		for (std::list<int>::iterator it=coarseRowsNumbers[i]->begin() ; it != coarseRowsNumbers[i]->end(); ++it)
		{
			fprintf(f,"%d ", *it);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

void AMGSolverStructure::InitVCycle()
{
	EndVCycle();
	smoothers = new Smoother*[levelsNumber];
	for(int i = 0; i<levelsNumber - 1; i++)
	{
		smoothers[i] = new GaussSidelMethod();
		smoothers[i]->PreSmoothing(amg_level_matrices[i],amg_level_rhs[i]);
	}
	smoothers[levelsNumber - 1] = new DirectSolver();
//	smoothers[levelsNumber - 1] = new GaussSidelMethod();
	smoothers[levelsNumber - 1]->PreSmoothing(amg_level_matrices[levelsNumber - 1],amg_level_rhs[levelsNumber - 1]);
}

void AMGSolverStructure::EndVCycle()
{
	if(smoothers != NULL)
		for(int i = 0; i<levelsNumber; i++)
		{
			delete smoothers[i];
		}

	delete[] smoothers;
	smoothers = NULL;
}

Vec AMGSolverStructure::RunVCycle()
{
	InitVCycle();
	int v_cycle_iteration_nr = 0;
	do
	{
		SolveFromLevel(0);
	}
	while(!error_evaluator->Stop(amg_level_x[0], ++v_cycle_iteration_nr));
	EndVCycle();

	Vec result;
	VecDuplicate(amg_level_x[0], &result);
	VecCopy(amg_level_x[0], result);
	return result;
}

//TODO: work around for seq
void sum_vectors(Vec v1, Vec v2){
	int ss;
	VecGetSize(v1,&ss);
	int* tmp_a = (int*)malloc(sizeof(int)*ss);
	for(int i = 0; i<ss; i++)
		tmp_a[i] = i;
	double* tmp_v;
	VecGetArray(v2, &tmp_v);
	VecSetValues(v1,ss,tmp_a,tmp_v,ADD_VALUES);
	VecRestoreArray(v2, &tmp_v);
	free(tmp_a);
	VecAssemblyBegin(v1);
	VecAssemblyEnd(v1);
}

void AMGSolverStructure::SolveFromLevel(int level)
{
	clock_t begin = clock();
	int presmoothing_count = 4;
	int postsmoothing_count = 4;

	if(level == levelsNumber -1)
	{
		SolveFinalLevel();
		return;
	}

	//Mat level_matrix = amg_level_matrices[level];
	smoothers[level]->SetRhs(amg_level_rhs[level]);
	if(level != 0)
		VecSet(amg_level_x[level],0.0);

	//VecView(amg_level_rhs[level], PETSC_VIEWER_STDOUT_WORLD);
	//pre smoothing
	for(int i = 0; i<presmoothing_count; i++)
		smoothers[level]->Smooth(amg_level_x[level]);
	//mf_log_info("%s %d %s %d","Level:  ", level, "Presmoothing iterations: ", presmoothing_count);


	MatResidual(amg_level_matrices[level],amg_level_rhs[level],amg_level_x[level],amg_level_residuals[level]);
	MatMultTranspose(amg_interpolation_matrices[level], amg_level_residuals[level], amg_level_rhs[level+1]);
	SolveFromLevel(level+1);

	//reuse level residual vector so that intermediate multiplication result does not have to be stored in another vec
	MatMult(amg_interpolation_matrices[level],amg_level_x[level+1],amg_level_residuals[level]);
//	VecAXPY(amg_level_x[level],1.0,amg_level_residuals[level]);
	sum_vectors(amg_level_x[level], amg_level_residuals[level]);


	//post smoothing
	for(int i = 0; i<postsmoothing_count; i++)
		smoothers[level]->Smooth(amg_level_x[level]);
	//mf_log_info("%s %d %s %d","Level:  ", level, "Postsmoothing iterations: ", postsmoothing_count);

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	printf("Level %d time: %lf\n",level,elapsed_secs);
}

void AMGSolverStructure::SolveFinalLevel()
{
	int final_level = levelsNumber - 1;
	if(final_level > 0)
		VecSet(amg_level_x[final_level],0.0);
	smoothers[final_level]->SetRhs(amg_level_rhs[final_level]);
	//for(int i = 0; i<8 ;i++)
		smoothers[final_level]->Smooth(amg_level_x[final_level]);
}

void AMGSolverStructure::RunVCycle(int levelId, double* X)
{
	//modfem uses reverse level numbers: highest number in modfem is related to the fines grid, highest numer in amg is related to the coarsest grid
	levelId = GetLevelsNumber() - levelId - 1;
	//(lav matrix id)

	if(smoothers == NULL)
		InitVCycle();

	SmoothLevel(levelId);

	if(levelId == 0)
	{
		PetscScalar* v;
		VecGetArray(amg_level_x[levelId], &v);
		int localRowNumbers;
		MatGetLocalSize(matrix, &localRowNumbers, NULL);
		for(int i = 0; i< localRowNumbers; i++)
		{
			X[i] = v[i];
		}

		VecRestoreArray(amg_level_x[levelId], &v);
	}
}

void AMGSolverStructure::MoveSolutionFromLevel(int level)
{
	level = GetLevelsNumber() - level - 1;
	MatResidual(amg_level_matrices[level],amg_level_rhs[level],amg_level_x[level],amg_level_residuals[level]);
	MatMultTranspose(amg_interpolation_matrices[level], amg_level_residuals[level], amg_level_rhs[level+1]);
	VecSet(amg_level_x[level+1],0.0);
}

void AMGSolverStructure::MoveSolutionToLevel(int level)
{
	level = GetLevelsNumber() - level - 1;
	PetscErrorCode ierr = MatMult(amg_interpolation_matrices[level],amg_level_x[level+1],amg_level_residuals[level]);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAXPY(amg_level_x[level],1.0,amg_level_residuals[level]); CHKERRABORT(PETSC_COMM_WORLD, ierr);
}


void AMGSolverStructure::SmoothLevel(int level)
{
	if(level == levelsNumber -1)
	{
		SolveFinalLevel();
		return;
	}
	smoothers[level]->SetRhs(amg_level_rhs[level]);
	smoothers[level]->Smooth(amg_level_x[level]);
}

int AMGSolverStructure::GetLevelsNumber()
{
	return levelsNumber;
}


void AMGSolverStructure::InitLevels()
{
	amg_level_matrices = new Mat[levelsNumber];
	amg_interpolation_matrices = new Mat[levelsNumber - 1];
	amg_level_rhs = new Vec[levelsNumber];
	amg_level_x = new Vec[levelsNumber];
	amg_level_residuals = new Vec[levelsNumber-1];
	amg_level_sizes = new int[levelsNumber];
	coarseRowsNumbers = new std::list<int>*[levelsNumber-1];

	amg_level_matrices[0] = matrix;
	amg_level_rhs[0] = rhs;
	amg_level_sizes[0] = size;

	VecCreate(PETSC_COMM_WORLD,&(amg_level_x[0]));
    int localRowNumbers;
    MatGetLocalSize(matrix, &localRowNumbers, NULL);
	VecSetSizes(amg_level_x[0],localRowNumbers,size);
	VecSetFromOptions(amg_level_x[0]);
	//VecDuplicate(this->rhs,&(this->amg_level_x[0]));
	VecSet(this->amg_level_x[0],0.0);
	if(levelsNumber > 1)
		VecDuplicate(amg_level_x[0], &(amg_level_residuals[0]));
	cf_splitters = (CFSplitter**) new RSCFSplitter*[levelsNumber - 1];
//	mf_log_info("%s %d","Numebr of AMG levels: ", levelsNumber);
}

void AMGSolverStructure::ChangeSystem(Vec rhs, int max_iterations)
{
	this->rhs = rhs;
	this->amg_level_rhs[0] = rhs;
	VecScale(amg_level_rhs[0], -1.0);
	VecSet(amg_level_x[0],0.0);
	delete error_evaluator;
	error_evaluator = new ResidualBasedErrorEvaluator(matrix, rhs, 1.0e-10, 1.0e-10, size, max_iterations);//new ResidualBasedErrorEvaluator(matrix, rhs, 1.0e-7, 1.0e-8, size);
}

void AMGSolverStructure::CreateNextLevel(int levelNumber)
{

}
