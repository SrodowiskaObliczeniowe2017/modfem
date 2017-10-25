#include "VectorTransformationUtil.hpp"

void convert_vec_to_x(Vec solution_v, Vec solution_p, double* x, int Ndof)
{
	PetscScalar* values_v;
	PetscScalar* values_p;
	VecGetArray(solution_v,&values_v);
	VecGetArray(solution_p,&values_p);

	//velocity
	for(int i = 0; i< Ndof * 3 / 4; i+=3)
	{
		x[i + i/3] = (double) values_v[i];
		x[i+1 + i/3] = (double) values_v[i+1];
		x[i+2 + i/3] = (double) values_v[i+2];
	}
	//pressure
	for(int i = 0; i< Ndof / 4; i++)
	{
		x[4*i + 3] = (double) values_p[i];
	}

	VecRestoreArray(solution_v,&values_v);
	VecRestoreArray(solution_p,&values_p);
//	VecView(solution, PETSC_VIEWER_STDOUT_WORLD);
//	printf("===================\n");
//	for(int i = 0; i<Ndof; i++)
//		printf("%lf\n", x[i]);
}

void convert_vec_to_x(Vec solution, double* x, int Ndof)
{
	PetscScalar* values;
	VecGetArray(solution,&values);

	//velocity
	for(int i = 0; i< Ndof * 3 / 4; i+=3)
	{
		x[i + i/3] = (double) values[i];
		x[i+1 + i/3] = (double) values[i+1];
		x[i+2 + i/3] = (double) values[i+2];
	}
	//pressure
	for(int i = 0; i< Ndof / 4; i++)
	{
		x[4*i + 3] = (double) values[i + Ndof * 3 /4];
	}

	VecRestoreArray(solution,&values);
//	VecView(solution, PETSC_VIEWER_STDOUT_WORLD);
//	printf("===================\n");
//	for(int i = 0; i<Ndof; i++)
//		printf("%lf\n", x[i]);
}

void convert_vec_to_x_simple(Vec solution, double* x, int Ndof)
{
	PetscScalar* values;
	VecGetArray(solution,&values);

	for(int i = 0; i< Ndof; i++)
	{
		x[i] = (double) values[i];
	}
	VecRestoreArray(solution,&values);
}

void convert_x_to_vec(Vec solution_v, Vec solution_p, double* x, int Ndof)
{
	PetscScalar* values;
	int index;
	for(int i = 0; i< Ndof; i+=4)
	{
		index = i - i/4; VecSetValues(solution_v,1,&index,&x[i],INSERT_VALUES);
		index = i+1 - i/4; VecSetValues(solution_v,1,&index,&x[i+1],INSERT_VALUES);
		index = i+2 - i/4; VecSetValues(solution_v,1,&index,&x[i+2],INSERT_VALUES);
		index = i/4; VecSetValues(solution_p,1,&index,&x[i+3],INSERT_VALUES);
	}

//	VecView(solution, PETSC_VIEWER_STDOUT_WORLD);
//	printf("===================\n");
//	for(int i = 0; i<Ndof; i++)
//		printf("%lf\n", x[i]);
}

void convert_x_to_vec(Vec solution, double* x, int Ndof)
{
	PetscScalar* values;
	int index;
	for(int i = 0; i< Ndof; i+=4)
	{
		index = i - i/4; VecSetValues(solution,1,&index,&x[i],INSERT_VALUES);
		index = i+1 - i/4; VecSetValues(solution,1,&index,&x[i+1],INSERT_VALUES);
		index = i+2 - i/4; VecSetValues(solution,1,&index,&x[i+2],INSERT_VALUES);
		index = Ndof * 3 / 4 + i/4; VecSetValues(solution,1,&index,&x[i+3],INSERT_VALUES);
	}

//	VecView(solution, PETSC_VIEWER_STDOUT_WORLD);
//	printf("===================\n");
//	for(int i = 0; i<Ndof; i++)
//		printf("%lf\n", x[i]);
}

void convert_x_to_vec_simple(Vec solution, double* x, int Ndof){
	PetscScalar* values;
	int index;
	for(int i = 0; i< Ndof; i+=1)
	{
		index = i; VecSetValues(solution,1,&index,&x[i],INSERT_VALUES);
	}
}
