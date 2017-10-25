#include "../inc/FieldReaderDG.h"
#include "../inc/FileManager.h"

#include "../../include/IOManager.h"
#include "../../include/Enums.h"
#include "../../include/FieldElems.h"
#include "../../include/MeshElems.h"
#include "aph_intf.h"

#include "../../include/fv_config.h"

#include <string>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>





static int apr_get_pdeg_nrdofs(int base, int pdeg)
{
	int porder, porderz, numdofs=-1;
/* check */
#ifdef FV_DEBUG
	if(pdeg<0 || (base == APC_COMPLETE && pdeg > 100)){
		printf("Wrong pdeg %d in apr_get_pdeg_nrdofs\n",pdeg);
		return( -1);	
	}
#endif
	
	if(base==APC_TENSOR_DG){

		/* decipher horizontal and vertical orders of approximation */
		porderz = pdeg/100;
		porder  = pdeg%100;
		numdofs = (porderz+1)*(porder+1)*(porder+2)/2;

	}
	else if(base==APC_COMPLETE_DG){
		porder  = pdeg;
		numdofs = (porder+1)*(porder+2)*(porder+3)/6;
	}
	else {
		printf("Type of base in get_pdeg_nrdofs not valid for prisms!\n");
		exit(-1);
	}

	return( numdofs);
}

namespace IOmgr {

	inline bool mmr_el_status(const mmt_elems* pel, int nel)
	{
		return ((pel + nel)->type > 0);
	}

	FieldReaderDG inst_FiledReaderDG;

	FieldReaderDG& InstanceOfFieldReaderDG()
	{
		return  inst_FiledReaderDG;
	}

	FieldReaderDG::FieldReaderDG()
		: _name			("FieldReaderDG"),
		  _extension	("dat dmp grd"),
		  _initials		("field"),
		  _description	("dg")
	{
#ifdef FV_DEBUG
		std::cout << _name << " alive.\n" << std::endl;
#endif
		FileManager::GetInstance()->RegisterModules(this);
	}

	bool 
	FieldReaderDG::Read(const char* fname, BaseImporter& ibimpr)
	{
		std::ifstream ifs(fname,std::ios::in);
		if(!ifs.is_open() || !ifs.good()){
			std::cerr << "FieldReaderDG: can't open file: " << fname << std::endl;
			return false;
		}
#ifdef _DEBUG
		std::cout << "Reading file: " << fname << std::endl; 
#endif

		bool res = Read(ifs,ibimpr);
		ifs.close();
		return res;
	}

	bool FieldReaderDG::IsReadAble(const char *fname) const
	{
		if(BaseFileReader::IsReadAble(fname)){
			std::ifstream ifs(fname);
			if(ifs.is_open() && CheckInitials(fname)){
				ifs.close();
				return true;
			}
			ifs.close(); //?
		}
		return false;
	}

	bool 
	FieldReaderDG::Read(std::istream& is, BaseImporter& ibimpr)
	{

		if(!is.good()){
			std::cerr << "FieldReaderDG: can't read from stream!" << std::endl;
			return false;
		}

		// static field
		static apt_dg_field fieldp;
		int pdeg;

		// Read the number of solution vectors for each dof entity and
		// the type of shape functions
		is >> fieldp.nreq; 
		is >> fieldp.nr_sol;
		is >> pdeg; 


		// Reserve but onlu the mnel number
		fieldp.dof_ents = static_cast<apt_dg_dof_ent*>(ibimpr.Reserve(FemViewer::FIELD));
		if(!fieldp.dof_ents){
			std::cerr << "Dofs structures not allocated!\n";
			return( false);
		}

		// Get the number of elements in associated mesh
		int nmel = ibimpr.nElems();
		if(nmel<=0){
			std::cerr << "Error! no elements in mesh\n";
			return( false);
		}

		// Get handle to elements in associated mesh
		const mmt_elems* pelem = static_cast<const mmt_elems*>(ibimpr.GetItem());
		//pelem++; // increment to the first element
		//FV_ASSERT(pelem != 0);
		// Read degree of polynomial for each element (or uniform p if pdeg <0 )
/*		int pdeg;
		is >> pdeg;	*/	
		
		if(pdeg<0) {
			pdeg = -pdeg;
			fieldp.uniform = 1;
			if(pdeg < 100) fieldp.base = FemViewer::COMPLETE;
			if(fieldp.base == FemViewer::TENSOR && pdeg<100){
				std::cerr << "WARNING: piecewise constant approximation in xy-plane; p="
					<< pdeg << " < 100 in input file\n";
			}
			if( ((fieldp.base == FemViewer::COMPLETE) && (pdeg > APC_MAXELP_COMP)) ||
				((fieldp.base == FemViewer::TENSOR)   && (pdeg > APC_MAXELP_TENS)) ) {
				   std::cerr << "Too big degree of polynomial for complete basis: " << pdeg
					   << std::endl;
				   return( false);
			}

			int iaux = 0;
			for(int nel = 1; nel <= nmel; nel++) {
				if(mmr_el_status(pelem,nel)){ // if ACTIVE
					iaux++;
					fieldp.dof_ents[nel].pdeg = pdeg;
				}
				else {
					fieldp.dof_ents[nel].pdeg = -1;//APC_NO_DOFS;
				}
			}
		}
		else {
			int iaux = 0;
			fieldp.uniform = 0;
			if(pdeg<100) fieldp.base = FemViewer::COMPLETE;
			else fieldp.base = FemViewer::TENSOR;
			for(int nel=1; nel<=nmel; nel++) {
				if(mmr_el_status(pelem,nel)) {
					iaux++;
					if(iaux == 1)  
						fieldp.dof_ents[nel].pdeg = pdeg;
					else 
						is >> fieldp.dof_ents[nel].pdeg;
				} else {
					fieldp.dof_ents[nel].pdeg = -1;//APC_NO_DOFS;
				}
			} // end for
		} // else

		 // Create dofs data structure
		int base, num_dof;
		for(int nel=1; nel <= nmel; nel++){
			if(mmr_el_status(pelem,nel)){
				pdeg = fieldp.dof_ents[nel].pdeg;
				base = fieldp.base;
				num_dof = apr_get_pdeg_nrdofs(base,pdeg);

				fieldp.dof_ents[nel].vec_dof_1 = 
					(double *) malloc(num_dof*sizeof(double));
				if(!fieldp.dof_ents[nel].vec_dof_1){
					std::cerr << "Error! Can't allocale memory for dofs vec1!";
					return( false);
				}
				if(fieldp.nr_sol>1){
					fieldp.dof_ents[nel].vec_dof_2 = 
						(double *) malloc(num_dof*sizeof(double));
					if(!fieldp.dof_ents[nel].vec_dof_2){
						std::cerr << "Error! Can't allocale memory for dofs vec2!";
						return( false);
					}
				}	  
				if(fieldp.nr_sol>2){
					fieldp.dof_ents[nel].vec_dof_3 = 
						(double *) malloc(num_dof*sizeof(double));
					if(!fieldp.dof_ents[nel].vec_dof_3){
						std::cerr << "Error! Can't allocale memory for dofs vec3!";
						return( false);
					}
				}

				// Read values from a dump file
				for(int i = 0; i < num_dof; i++){
					is >> fieldp.dof_ents[nel].vec_dof_1[i];
				}
				if(fieldp.nr_sol>1){
					for(int i = 0; i <num_dof; i++){
						is >> fieldp.dof_ents[nel].vec_dof_2[i];
				}}
				if(fieldp.nr_sol>2){
					for(int i = 0; i <num_dof; i++){
						is >> fieldp.dof_ents[nel].vec_dof_3[i];
				}}

			} // end if
		} // end for

//#ifdef FV_DEBUG
//		std::cout << "Number of manged dof structures: " << i << std::endl;
//#endif

		ibimpr.ImportParams(&fieldp);

 	
		return true;


	}

	

	
}
