/*
 * mic_platrom.cpp
 *
 *  Created on: 16 lis 2015
 *      Author: pmaciol
 */
#include <stdio.h>
#include <unistd.h>
#include <omp.h>
#include <vector>
#include "fv_timer.h"
#include "mic.h"
#include "types.h"
#include "Log.h"
#include "Mesh.h"
#include "Field.h"
#include "BBox3D.h"

template<typename T>
inline void initarray(T* arr,const int size,const T value = 0) {
	for(int i=0;i<size;i++) arr[i] = value;
}

template<typename T>
inline void min_array(T* arr,T const * sol,const int size) {
	for(int i=0;i<size;i++) arr[i] = fv_min(arr[i],sol[i]);
}

template<typename T>
inline void max_array(T* arr,T const * sol,const int size) {
	for(int i=0;i<size;i++) arr[i] = fv_max(arr[i],sol[i]);
}

namespace FemViewer {
int init_platorm(int argc,char** argv)
{
	mfp_log_info("The number of available logic cores: %d\n",sysconf(_SC_NPROCESSORS_ONLN ));
#ifdef MIC
	int is_init = 1;
#else
	int is_init = 0;
#endif
	return is_init;
}

int init_data(int argc,char** argv)
{

	return 0;
}


BBox3D ExtractAxisAlignedMeshExtens(const Mesh* mesh_ptr,int do_parallel)
{
	using fvmath::Vec3d;

	BBox3D mesh_bbox;
	int num_inact = 0;
	int num_nodes = mesh_ptr->GetNumNodes();
	fv_timer t;
	t.start();

	#pragma omp parallel if(do_parallel) default(none) \
		firstprivate(mesh_ptr,num_nodes) shared(mesh_bbox,num_inact)
	{
		BBox3D lbbox;
		Vec3d lcoord;
		int lianno = 0;

		#pragma omp for schedule(dynamic,4)
		for (int i=1; i <= num_nodes; ++i)
		{
			if (mesh_ptr->GetNodeCoor(i,lcoord.v)) {
				lbbox += lcoord;
			}
			else {
				++lianno;
			}
		}

		#pragma omp critical(bbox)
		{
			mesh_bbox += lbbox;
			num_inact += lianno;
		}
	}
	t.stop();
	assert(mesh_bbox.isInitialized());
	mfp_log_info("Building%bounding box info in %lf sec",(do_parallel ? " in parallel " : " "), t.get_time_sec() );
	mfp_log_info("Found %d number of inactive nodes", num_inact);
	return mesh_bbox;
}

size_t ExtractAxisAlignedElementExtens(const Mesh* mesh_ptr,std::vector<BBox3D>& elem_bboxes,int do_parallel)
{
	using fvmath::Vec3d;

	int num_inact = 0;
	int num_elems = mesh_ptr->GetNumElems();

	elem_bboxes.reserve(num_elems);
	elem_bboxes.clear();

	fv_timer t;
	t.start();

	#pragma omp parallel if(do_parallel) default(none) \
		firstprivate(mesh_ptr,num_elems) shared(elem_bboxes)
	{
		BBox3D lbbox;
		CVec3d elCoords[MMC_MAXELVNO];

		size_t size = num_elems / omp_get_num_threads();
		if (size == 0) size = 1;
		std::vector<BBox3D> vbboxs(size);

		#pragma omp for schedule(dynamic)
		for (int iel=1; iel <= num_elems; ++iel)
		{
			if (mesh_ptr->GetElemStatus(iel) != MMC_ACTIVE) continue;
			int nno = mesh_ptr->GetElementCoordinates(iel,nullptr,elCoords[0].v);
			for (int ino=0;ino<nno;++ino) {
				lbbox += elCoords[ino];
			}
		}

		#pragma omp critical(elembbox)
		{
			elem_bboxes.insert(elem_bboxes.end(),vbboxs.begin(),vbboxs.end());
		}
	}
	t.stop();
	assert(elem_bboxes.size() > 0);
	mfp_log_info("Building%bounding element boxes in %lf sec",(do_parallel ? " in parallel " : " "), t.get_time_sec() );
	return elem_bboxes.size();
}

#define PDC_MAXEQ 5
// Calculate min and max values by sampling volume
int CalculateMinMaxValues(const Mesh* mesh_ptr,
		                  const Field* field_ptr,
		                  const int ctrl,
						  double*& min_max_array )
{
	if (!mesh_ptr || !field_ptr) return 0;

	double coords[3*MMC_MAXELVNO];
	double elmin[PDC_MAXEQ];
	double elmax[PDC_MAXEQ];
	std::vector<double> coeffs(APC_MAXELVD);
	std::vector<double> vSol(PDC_MAXEQ);

#if Math_Calculator == 1
	std::vector<fvmathParser::MathElement> mathElems;
	fvmathParser::MathCalculator::ONP(
				field_ptr->GetSolution().sFormula,
				field_ptr->Getsolution().nrEquetions,
				mathElems);
#endif

	int Pdeg[3],count = 0;
	int num_elems = mesh_ptr->GetNumElems();
	min_max_array = (double*)malloc(2*sizeof(double)*(num_elems+1));

	fv_timer t;
	t.start();
	for (int iel=1; iel <= num_elems; iel++)
	{
		// Check status
		if (mesh_ptr->GetElemStatus(iel) != MMC_ACTIVE) continue;
		// Get Coordinates
		int nno = mesh_ptr->GetElementCoordsd(iel, nullptr, coords);
		// Get base type and vector of degrees
		int base = field_ptr->GetDegreeVector(iel,Pdeg);
		// Get number of shape function
		int nrshp = field_ptr->GetNumberOfShapeFunc(iel,Pdeg);
		// Get element coefficients
		int nrdofs = field_ptr->GetElementDofs(iel,coeffs.data());
		// Loops over 10 point in each directions
		// Start at the first vertex in prism
		// Clear temporary arrays
		initarray(elmin,PDC_MAXEQ,FV_LARGE);
		initarray(elmax,PDC_MAXEQ,-FV_LARGE);
		CVec3d pt;
		if (nno == 6) { // prism
			//mfp_debug("Prism%d\n",iel);
			for (pt.z = -1.0; pt.z <= 1.001; pt.z += 0.2) {
				for (pt.y = 0.0; pt.y <= 1.001; pt.y+= 0.1) {
					for (pt.x = 0.0; pt.x <= 1.001 - pt.y; pt.x += 0.1) {
						field_ptr->CalculateSolution(ctrl,Pdeg[0],base,coords,pt.v,coeffs.data(),vSol.data());
#if Math_Calculator == 1
						fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol, sol[0]);
#endif
						// Element compare
						min_array(elmin,vSol.data(),PDC_MAXEQ);
						max_array(elmax,vSol.data(),PDC_MAXEQ);
						//printf("pt: {%lf %lf %lf}\n",pt.x,pt.y,pt.z);
						//printf("elmin: %lf elmax: %lf\n",elmin,elmax);
					}
				}
			}
		}
		else { // tetrahedron
			for (pt.z = 0.0; pt.z <= 1.001; pt.z += 0.1) {
				for (pt.y = 0.0; pt.y <= 1.001 - pt.z; pt.y += 0.1) {
					for (pt.x = 0.0; pt.x <= 1.001 - pt.y - pt.z; pt.x += 0.1) {
						field_ptr->CalculateSolution(ctrl,Pdeg[0],base,coords,pt.v,coeffs.data(),vSol.data());
#if Math_Calculator == 1
						fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol, sol[0]);
#endif
						// Element compare
						min_array(elmin,vSol.data(),PDC_MAXEQ);
						max_array(elmax,vSol.data(),PDC_MAXEQ);
					}
				}
			}
		}
		// Store global min/max values
		double* mn_mx_ptr = &min_max_array[0];
		min_array(mn_mx_ptr,elmin,PDC_MAXEQ);
		mn_mx_ptr = &min_max_array[PDC_MAXEQ];
		max_array(mn_mx_ptr,elmax,PDC_MAXEQ);

		// Store in global array
		mn_mx_ptr = &min_max_array[iel+1];
		for(int i=0;i<PDC_MAXEQ;i++) *mn_mx_ptr++ = elmin[i];
		for(int i=0;i<PDC_MAXEQ;i++) *mn_mx_ptr++ = elmax[i];

		// increment counter of added elems
		count++;
	}

	t.stop();
	if (count) mfp_log_debug("Found Min: %lf\tMax: %lf values in %d elems\n",min_max_array[0],min_max_array[PDC_MAXEQ],count);
	mfp_log_debug("Time for searching min/max value: %lf\n",t.get_time_sec() );

	return count;
}

int create_grid(double targetoccupancy,
		        const uint32_t numelems,
		        const BBox3D& gbox,
		        const BBox3D* boxes,
		        int** griddata,
		        int** eldata,
		        grid_t* gridinfo
		        )
{
	using namespace fvmath;
	size_t gridsize = gridinfo->resolution[0] * gridinfo->resolution[1] * gridinfo->resolution[2];
		//printf("cellCOunts= %d %d %d\tgridsize = %d\n",gridinfo->cellCount.x,gridinfo->cellCount.y,gridinfo->cellCount.z, gridsize);
	int* lgriddata = new int [gridsize+1];
	memset(lgriddata,0x0,sizeof(int)*(gridsize+1));
	int *leldata;
#pragma omp parallel default(none) shared(gridinfo,gbox,leldata) firstprivate(boxes,targetoccupancy,lgriddata,gridsize)
{
	double mnDimX = static_cast<double>(gbox.mn.x);
	double mnDimY = static_cast<double>(gbox.mn.y);
	double mnDimZ = static_cast<double>(gbox.mn.z);
	double mxDimX = static_cast<double>(gbox.mx.x);
	double mxDimY = static_cast<double>(gbox.mx.y);
	double mxDimZ = static_cast<double>(gbox.mx.z);
	double xlength = mxDimX - mnDimX;
	double ylength = mxDimY - mnDimY;
	double zlength = mxDimZ - mnDimZ;

#pragma omp single
	{
	gridinfo->min[0] = mnDimX;
	gridinfo->min[1] = mnDimY;
	gridinfo->min[2] = mnDimZ;
	gridinfo->max[0] = mxDimX;
	gridinfo->max[1] = mxDimY;
	gridinfo->max[2] = mxDimZ;
	double norm = std::cbrt((numelems/targetoccupancy)/(xlength*ylength*zlength));
	gridinfo->resolution[0]=static_cast<int>(fmax(1,std::roundf(xlength*norm)));
	gridinfo->resolution[1]=static_cast<int>(fmax(1,std::roundf(ylength*norm)));
	gridinfo->resolution[2]=static_cast<int>(fmax(1,std::roundf(zlength*norm)));
	gridinfo->cellSize[0]=xlength/gridinfo->resolution[0];
	gridinfo->cellSize[1]=ylength/gridinfo->resolution[1];
	gridinfo->cellSize[2]=zlength/gridinfo->resolution[2];
	}

	//cl_int3 minidx;
	//cl_int3 maxidx;
	CVec3i minidx, maxidx;
#pragma omp for
	for (uint32_t i=0;i<numelems;i++)
	{
		//count into griddata
		double mnX = boxes[i].mn.x;
		double mnY = boxes[i].mn.y;
		double mnZ = boxes[i].mn.z;
		double mxX = boxes[i].mx.x;
		double mxY = boxes[i].mx.y;
		double mxZ = boxes[i].mx.z;
		minidx.x=fmin(gridinfo->resolution[0]-1,(int)std::floor(((mnX-mnDimX)*gridinfo->resolution[0])/xlength));
		maxidx.x=fmin(gridinfo->resolution[0]-1,(int)std::floor(((mxX-mnDimX)*gridinfo->resolution[0])/xlength));
		minidx.y=fmin(gridinfo->resolution[1]-1,(int)std::floor(((mnY-mnDimY)*gridinfo->resolution[1])/ylength));
		maxidx.y=fmin(gridinfo->resolution[1]-1,(int)std::floor(((mxY-mnDimY)*gridinfo->resolution[1])/ylength));
		minidx.z=fmin(gridinfo->resolution[2]-1,(int)std::floor(((mnZ-mnDimZ)*gridinfo->resolution[2])/zlength));
		maxidx.z=fmin(gridinfo->resolution[2]-1,(int)std::floor(((mxZ-mnDimZ)*gridinfo->resolution[2])/zlength));
		for (int x=minidx.x;x<=maxidx.x;x++)
			for (int y=minidx.y;y<=maxidx.y;y++)
				for (int z=minidx.z;z<=maxidx.z;z++)
				{
					//if (cell_triangle_intersect(triangles[i],
					//	(gridinfo->minDimension.x+(x*gridinfo->cellSize.x)),(gridinfo->minDimension.x+((x+1)*gridinfo->cellSize.x)),
					//	(gridinfo->minDimension.y+(y*gridinfo->cellSize.y)),(gridinfo->minDimension.y+((y+1)*gridinfo->cellSize.y)),
					//	(gridinfo->minDimension.z+(z*gridinfo->cellSize.z)),(gridinfo->minDimension.z+((z+1)*gridinfo->cellSize.z)) ) == 1)
						(lgriddata[x + (gridinfo->resolution[0] * y) + (gridinfo->resolution[0] * gridinfo->resolution[1] * z)])++;
				}
 	}

#pragma omp single
	{
	for (uint32_t i=1; i<=gridsize; i++) lgriddata[i] += lgriddata[i-1];
    leldata = new int[ lgriddata[gridsize] ];
	}

#pragma omp for
    for (uint32_t i=0;i<numelems;i++)
	{
    	double mnX = boxes[i].mn.x;
    	double mnY = boxes[i].mn.y;
    	double mnZ = boxes[i].mn.z;
    	double mxX = boxes[i].mx.x;
    	double mxY = boxes[i].mx.y;
    	double mxZ = boxes[i].mx.z;
    	minidx.x=fmin(gridinfo->resolution[0]-1,(int)std::floor(((mnX-mnDimX)*gridinfo->resolution[0])/xlength));
    	maxidx.x=fmin(gridinfo->resolution[0]-1,(int)std::floor(((mxX-mnDimX)*gridinfo->resolution[0])/xlength));
    	minidx.y=fmin(gridinfo->resolution[1]-1,(int)std::floor(((mnY-mnDimY)*gridinfo->resolution[1])/ylength));
    	maxidx.y=fmin(gridinfo->resolution[1]-1,(int)std::floor(((mxY-mnDimY)*gridinfo->resolution[1])/ylength));
    	minidx.z=fmin(gridinfo->resolution[2]-1,(int)std::floor(((mnZ-mnDimZ)*gridinfo->resolution[2])/zlength));
    	maxidx.z=fmin(gridinfo->resolution[2]-1,(int)std::floor(((mxZ-mnDimZ)*gridinfo->resolution[2])/zlength));
		for (int x=minidx.x;x<=maxidx.x;x++)
			for (int y=minidx.y;y<=maxidx.y;y++)
				for (int z=minidx.z;z<=maxidx.z;z++)
				{
				//	if (cell_triangle_intersect(triangles[i],
				//			(gridinfo->minDimension.x+(x*gridinfo->cellSize.x)),(gridinfo->minDimension.x+((x+1)*gridinfo->cellSize.x)),
				//			(gridinfo->minDimension.y+(y*gridinfo->cellSize.y)),(gridinfo->minDimension.y+((y+1)*gridinfo->cellSize.y)),
				//			(gridinfo->minDimension.z+(z*gridinfo->cellSize.z)),(gridinfo->minDimension.z+((z+1)*gridinfo->cellSize.z)) ) == 1)
					lgriddata[x + (gridinfo->resolution[0]*y) + (gridinfo->resolution[0]*gridinfo->resolution[1] * z)]--;
					leldata[lgriddata[x + (gridinfo->resolution[0]*y) + (gridinfo->resolution[0]*gridinfo->resolution[1] * z)]]=i;
				}
	}
}// end pragma

	*griddata=lgriddata;
	*eldata=leldata;
	return 0;
}


}





