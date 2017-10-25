
#ifndef AMGSOLVERSTRUCTURE_H
#define AMGSOLVERSTRUCTURE_H

#include <petscksp.h>
#include <petscmat.h>
#include "uth_log.h"
#include "../amg/RSCFSplitter.hpp"
#include "../amg/GeneralRSCFSplitter.hpp"
#include "../amg/GaussSidelMethod.hpp"
#include "../amg/DirectSolver.hpp"
#include "../amg/Interpolation/DirectInterpolation.hpp"
#include "../amg/Interpolation/GeneralInterpolation.hpp"
#include "ErrorEvaluator.hpp"
#include "KnownSolutionErrorEvaluator.hpp"
#include "ResidualBasedErrorEvaluator.hpp"


class AMGSolverStructure
{
	public:

		AMGSolverStructure(Mat matrix, Vec rhs, int size, ErrorEvaluator* error_evaluator, bool create_coarsening_scheme);
		~AMGSolverStructure();
		void CreateAMGSolverLevels(int coarsening_algorithm, int interpolation_algorithm, double strength_threshold, int levels_number);
		Vec RunVCycle();
		void RunVCycle(int levelId, double* X);
		void InitVCycle();
		void EndVCycle();
		void ExtractCoarseRowsNumbers();
		void ChangeSystem(Vec rhs, int max_iterations);
		void MoveSolutionFromLevel(int level);
		void MoveSolutionToLevel(int level);

	private:
		RSCFSplitter* GetSplitter(int levelNumber, double strength_threshold, int coarsening_algorithm, int interpolation_algorithm);
		void SmoothLevel(int level);
		void SolveFromLevel(int level);
		void SolveFinalLevel();


		Mat matrix;
		Vec rhs;
		int size;
		bool create_coarsening_scheme;


		Mat* amg_level_matrices;
		Mat* amg_interpolation_matrices;
		Vec* amg_level_rhs;
		Vec* amg_level_x;
		Vec* amg_level_residuals;
		int* amg_level_sizes;
		int levelsNumber;

		CFSplitter** cf_splitters;
		Smoother** smoothers;
		ErrorEvaluator* error_evaluator;
		std::list<int> **coarseRowsNumbers;

		int GetLevelsNumber();
		void InitLevels();
		void CreateNextLevel(int levelNumber);

};

#endif
