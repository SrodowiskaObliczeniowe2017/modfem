/*
 * RSCFSplitter.h
 *
 *  Created on: Jun 6, 2015
 *      Author: damian
 */

#ifndef RSCFSPLITTER_H_
#define RSCFSPLITTER_H_

#include "CFSplitter.hpp"
#include "AMGRows.hpp"
#include "RowDataCollector.hpp"
#include <petscmat.h>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <map>
#include "Interpolation/DirectInterpolation.hpp"
#include "Interpolation/AdvancedRSCFInterpolation.hpp"
#include "uth_log.h"

//Ruge Stuben CF Splitter
class RSCFSplitter : public CFSplitter
{
	public:
		RSCFSplitter(Mat mat, double strength_threshold, InterpolationStrategy* strategy);
		virtual ~RSCFSplitter();
		virtual void MakeCFSplitting();
		void GetSetsWithInfluencedNumebr(std::map<int, int> *c_set, std::map<int, int> *f_set);
		Mat GetMatrixFromCoarseToFine();
		virtual std::list<int>* GetNumbersOfCoarseRows();

	protected:
		virtual PetscScalar GetExtremalValue(const PetscScalar* values, PetscInt size, PetscInt row_index,
				const PetscInt* columns, int first_row_in_range, int range_end);
		virtual bool IsStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
				PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);

	private:
		InterpolationStrategy* strategy;
		double strength_threshold;
		struct row_info* row_info_array;
		struct influenced_info* influenced_info_array;
		PetscInt first_row_in_range;
		PetscInt range_end;
		void ReplaceAfterInfluencedNumberChange(struct row_info* row_info, bool switch_root);
		struct row_info* getHeapParent(int index);
		void replaceRowInfos(struct row_info* row_info1, struct row_info* row_info2);
		void makeHeap(int size);
		void PrepareFromCoarseToFineMat(Mat* from_coarse_to_fine_mat, int rows_number, int columns_number);


};

#endif /* RSCFSPLITTER_H_ */
