#include "SchurComplement.hpp"


void create_s_complement(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat* S, Vec bvp){
	  PetscErrorCode ierr;

	    Mat         AdB;
	    Vec         diag;

	    ierr = VecDuplicate(bvp,&diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_LUMP) {

	    ierr = MatGetRowSum(Avv,diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    const PetscScalar* values;
	    const PetscInt* columns;
	    PetscInt columns_number;
	    int size;
	    MatGetSize(Avv,&size,NULL);
	    for(int i = 0; i<size; i++){
	    	MatGetRow(Avv,i, &columns_number,&columns,&values);
	    	double v = 0;
	    	for(int j = 0; j<columns_number; j++){
	    		if(values[j] < 0)
	    			v -= 2*values[j];
	    	}
	    	VecSetValues(diag, 1, &i, &v, ADD_VALUES);
	    	MatRestoreRow(Avv, i, &columns_number, &columns, &values);
	    }


//	    } else if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_DIAG) {
//	      ierr = MatGetDiagonal(A00,diag);CHKERRQ(ierr);
//	    } else SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Unknown MatSchurComplementAinvType: %D", ainvtype);
	    ierr = VecReciprocal(diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    VecView(diag, PETSC_VIEWER_STDOUT_WORLD);
	    ierr = MatDuplicate(Avp,MAT_COPY_VALUES,&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr = MatDiagonalScale(AdB,diag,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr = VecDestroy(&diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    /* Cannot really reuse Spmat in MatMatMult() because of MatAYPX() -->
	         MatAXPY() --> MatHeaderReplace() --> MatDestroy_XXX_MatMatMult()  */
	    ierr     = MatDestroy(S);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr     = MatMatMult(Apv,AdB,MAT_INITIAL_MATRIX,PETSC_DEFAULT,S);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	    /* TODO: when can we pass SAME_NONZERO_PATTERN? */
	    ierr     = MatAYPX(*S,-1,App,DIFFERENT_NONZERO_PATTERN);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr     = MatDestroy(&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);

}

void create_s_complement_Vmm(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat* S, Mat Vmm, Vec bvp){
	  PetscErrorCode ierr;

	    Mat         AdB;
	    Vec         diag;

	    ierr = VecDuplicate(bvp,&diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_LUMP) {

	    MatGetDiagonal(Vmm,diag);


//	    } else if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_DIAG) {
//	      ierr = MatGetDiagonal(A00,diag);CHKERRQ(ierr);
//	    } else SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Unknown MatSchurComplementAinvType: %D", ainvtype);
	    ierr = VecReciprocal(diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    VecView(diag, PETSC_VIEWER_STDOUT_WORLD);
	    ierr = MatDuplicate(Avp,MAT_COPY_VALUES,&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr = MatDiagonalScale(AdB,diag,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr = VecDestroy(&diag);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    /* Cannot really reuse Spmat in MatMatMult() because of MatAYPX() -->
	         MatAXPY() --> MatHeaderReplace() --> MatDestroy_XXX_MatMatMult()  */
	    ierr     = MatDestroy(S);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr     = MatMatMult(Apv,AdB,MAT_INITIAL_MATRIX,PETSC_DEFAULT,S);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	    /* TODO: when can we pass SAME_NONZERO_PATTERN? */
	    ierr     = MatAYPX(*S,-1,App,DIFFERENT_NONZERO_PATTERN);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr     = MatDestroy(&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);

}

Mat create_s_complement_spai(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat restricted_matrix, Mat* S){
	  PetscErrorCode ierr;
//	  	MatView(Avv, PETSC_VIEWER_STDOUT_WORLD);
	    Mat         AdB;
	    Mat inverse = getApproximateInverseOpt(Avv,50);
//	    Mat inverse = getApproximateInverse(Avv);


//	    } else if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_DIAG) {
//	      ierr = MatGetDiagonal(A00,diag);CHKERRQ(ierr);
//	    } else SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Unknown MatSchurComplementAinvType: %D", ainvtype);

//	    VecView(diag, PETSC_VIEWER_STDOUT_WORLD);
	    ierr     = MatMatMult(inverse,Avp,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    /* Cannot really reuse Spmat in MatMatMult() because of MatAYPX() -->
	         MatAXPY() --> MatHeaderReplace() --> MatDestroy_XXX_MatMatMult()  */
	    ierr     = MatDestroy(S);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr     = MatMatMult(Apv,AdB,MAT_INITIAL_MATRIX,PETSC_DEFAULT,S);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	    /* TODO: when can we pass SAME_NONZERO_PATTERN? */
	    ierr     = MatAYPX(*S,-1,App,DIFFERENT_NONZERO_PATTERN);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	    ierr     = MatDestroy(&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	    return inverse;
}

//void create_s_complement_block(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat* S){
//	  PetscErrorCode ierr;
//
//	    Mat         AdB;
//	    Mat Avv_inverse_block;
//
//
//
////	    if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_LUMP) {
//
//
//	    const PetscScalar* values;
//	    const PetscInt* columns;
//	    PetscInt columns_number;
//	    int size;
//	    MatGetSize(amg_ns_supg_solver_data.Avv,&size,NULL);
//	    create_ns_supg_matrix(&(Avv_inverse_block), size, size);
//		ierr = MatMPIAIJSetPreallocation(Avv_inverse_block,3,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		ierr = MatSeqAIJSetPreallocation(Avv_inverse_block,3,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		ierr = MatSetUp(Avv_inverse_block); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    for(int i = 0; i<size; i+=3){
//	    	double block[9];
//	    	MatGetRow(Avv,i, &columns_number,&columns,&values);
//	    	for(int j = 0; j<columns_number; j++){
//	    		if(columns[j] % 3 == 0)
//	    			block[0] += fabs(values[j]);
//	    		else if(columns[j] % 3 == 1)
//	    			block[1] += fabs(values[j]);
//	    		else if(columns[j] % 3 == 2)
//	    			block[2] += fabs(values[j]);
//	    	}
//	    	MatRestoreRow(Avv, i, &columns_number, &columns, &values);
//	    	MatGetRow(Avv,i+1, &columns_number,&columns,&values);
//	    	for(int j = 0; j<columns_number; j++){
//	    		if(columns[j] % 3 == 0)
//	    			block[3] += fabs(values[j]);
//	    		else if(columns[j] % 3 == 1)
//	    			block[4] += fabs(values[j]);
//	    		else if(columns[j] % 3 == 2)
//	    			block[5] += fabs(values[j]);
//	    	}
//	    	MatRestoreRow(Avv, i+1, &columns_number, &columns, &values);
//	    	MatGetRow(Avv,i+2, &columns_number,&columns,&values);
//	    	for(int j = 0; j<columns_number; j++){
//	    		if(columns[j] % 3 == 0)
//	    			block[6] += fabs(values[j]);
//	    		else if(columns[j] % 3 == 1)
//	    			block[7] += fabs(values[j]);
//	    		else if(columns[j] % 3 == 2)
//	    			block[8] += fabs(values[j]);
//	    	}
//	    	MatRestoreRow(Avv, i+2, &columns_number, &columns, &values);
//
//	    	double inverse_factor = block[2]*(block[4]*block[6] - block[3]*block[7]) +
//	    			block[1]*(block[3]*block[8] - block[5]*block[6]) + block[0]*(block[5]*block[7] - block[4]*block[8]);
//	    	double inverse[9];
//	    	inverse[0] = (block[5]*block[7]-block[4]*block[8])/inverse_factor;
//	    	inverse[1] = (-block[2]*block[7]+block[1]*block[8])/inverse_factor;
//	    	inverse[2] = (block[2]*block[4]-block[1]*block[5])/inverse_factor;
//	    	inverse[3] = (-block[5]*block[6]+block[3]*block[8])/inverse_factor;
//	    	inverse[4] = (block[2]*block[6]-block[0]*block[8])/inverse_factor;
//	    	inverse[5] = (-block[2]*block[3]+block[0]*block[5])/inverse_factor;
//	    	inverse[6] = (block[4]*block[6]-block[3]*block[7])/inverse_factor;
//	    	inverse[7] = (-block[1]*block[6]+block[0]*block[7])/inverse_factor;
//	    	inverse[8] = (block[1]*block[3]-block[0]*block[4])/inverse_factor;
//	    	PetscInt idxm[3];
//	    	PetscInt idxn[3];
//	    	idxm[0] = i;
//	    	idxm[1] = i+1;
//	    	idxm[2] = i+2;
//	    	MatSetValues(Avv_inverse_block,3,idxm,3,idxm,inverse,INSERT_VALUES);
//	    }
//
//		ierr = MatAssemblyBegin(Avv_inverse_block,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		ierr = MatAssemblyEnd(Avv_inverse_block,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//
////	    } else if (ainvtype == MAT_SCHUR_COMPLEMENT_AINV_DIAG) {
////	      ierr = MatGetDiagonal(A00,diag);CHKERRQ(ierr);
////	    } else SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Unknown MatSchurComplementAinvType: %D", ainvtype);
//	    ierr = MatMatMult(Avv_inverse_block,Avp,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	    /* Cannot really reuse Spmat in MatMatMult() because of MatAYPX() -->
//	         MatAXPY() --> MatHeaderReplace() --> MatDestroy_XXX_MatMatMult()  */
//	    ierr     = MatDestroy(S);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    ierr     = MatMatMult(Apv,AdB,MAT_INITIAL_MATRIX,PETSC_DEFAULT,S);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	    /* TODO: when can we pass SAME_NONZERO_PATTERN? */
//	    ierr     = MatAYPX(*S,-1,App,DIFFERENT_NONZERO_PATTERN);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	    ierr     = MatDestroy(&Avv_inverse_block);
//	    ierr     = MatDestroy(&AdB);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//}
