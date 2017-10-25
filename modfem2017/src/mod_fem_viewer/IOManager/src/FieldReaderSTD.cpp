#include "FieldReaderSTD.h"
#include "FileManager.h"
#include "io_config.h"
#include "../../include/Enums.h"
#include "../../include/FieldElems.h"
#include "../../include/MeshElems.h"
#include "aph_intf.h"

#include "../include/fv_compiler.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>




// some defines
//#ifndef FV_USE_EXTERNALS
//#define APC_MAXELP_COMP 9
//#define APC_MAXELP_TENS 909
//#else
#define APC_TRUE 1
#define APC_FALSE 0
//#endif


namespace IOmgr {

	inline bool mmr_el_status(const mmt_elems* pel, int nel)
	{
		bool res = ((pel + nel)->type > 0);
		return res;
	}

	

	FieldReaderSTD inst_FiledReaderSTD;

	FieldReaderSTD& InstanceOfFieldReaderSTD()
	{
		return  inst_FiledReaderSTD;
	}

	FieldReaderSTD::FieldReaderSTD()
		: _name("FieldReaderSTD"), _extension("dmp"),
		  _initials("field_std"), _description("Reader of filed data for STD problems")
	{
#ifdef FV_DEBUG
		std::cout << _name << " alive.\n";
#endif
		FileManager::GetInstance()->RegisterModules(this);
	}

	bool 
	FieldReaderSTD::Read(const char* fname, BaseImporter& ibimpr)
	{
		std::ifstream ifs(fname,std::ios::in);
		if(!ifs.is_open() || !ifs.good()){
			std::cerr << "Can't open file: " << fname << std::endl;
			return false;
		}
#ifdef _DEBUG
		std::cout << "Reading file: " << fname << std::endl; 
#endif

		bool res = Read(ifs,ibimpr);
		ifs.close();
		return res;
	}

	bool FieldReaderSTD::IsReadAble(const char *fname) const
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
	FieldReaderSTD::Read(std::istream& is, BaseImporter& ibimpr)
	{

		if(!is.good()){
			std::cerr << "FieldReaderSTD: can't read from stream!" << std::endl;
			return false;
		}

		static apt_std_field fieldp;
		fieldp.mesh_id = -1;
		fieldp.base = FemViewer::TENSOR;
		fieldp.pdeg = 101;
		fieldp.constr = APC_FALSE;
		fieldp.uniform = FemViewer::F_TRIA;
		fieldp.dof_ents = NULL;

		// Read the number of components inssolution vectors 
		// Degree of polynomials in each element
		// Number of nodes
		is >> fieldp.nreq; 
		is >> fieldp.nr_sol;
		is >> fieldp.nr_nodes;

		if(fieldp.nreq > F_MAXEQ) {
			std::cerr<< "number of equations: " << fieldp.nreq
				<< " is larger than maximal allowed: " << F_MAXEQ
				<< std::endl;
			return(false);
		}

		ibimpr.ImportParams(&fieldp);

		// Reserve
		//std::cout<<"Reader: before reserve\n";
		fieldp.dof_ents = static_cast<apt_std_dof_ent*>(ibimpr.Reserve( FemViewer::FIELDN ));
		//std::cout<<"Reader: after reserve\n";
		if(!fieldp.dof_ents){
			std::cerr << "Dofs structures not allocated!\n";
			return( false);
		}

		// Read in the coeficinet of solution
		int tmp; 
#ifdef FV_DEBUG
		int tmp1 = 0, counter = 0;
#endif
		for(int i = 1; i <= fieldp.nr_nodes; ++i) 
		{
			is >> tmp; // read index of node
			fieldp.dof_ents[tmp].vec_dof_1 = static_cast<double*>(malloc(fieldp.nreq*sizeof(double)));
			assert(fieldp.dof_ents[tmp].vec_dof_1 != NULL);
//#ifdef FV_DEBUG
//			if(tmp != (tmp1 +1)) {
//				std::cout<<"tmp= " << tmp << "tmp1= " << tmp1 << " ";
//				std::cout<< "wezel zwiazany: " << (tmp-1) << std::endl;
//				counter++;
//			}
//#endif
			if(fieldp.nr_sol>1) {
				fieldp.dof_ents[tmp].vec_dof_2 = static_cast<double*>(malloc(fieldp.nreq*sizeof(double)));
			}
			if(fieldp.nr_sol>2) {
				fieldp.dof_ents[tmp].vec_dof_3 = static_cast<double*>(malloc(fieldp.nreq*sizeof(double)));
			}

			// reading dofs's data
			for(int j = 0; j < fieldp.nreq; ++j) {
				is >> fieldp.dof_ents[tmp].vec_dof_1[j];
			}

			if(fieldp.nr_sol > 1) {
				for(int j = 0; j < fieldp.nreq; ++j) {
				 is >> fieldp.dof_ents[tmp].vec_dof_2[j];
				}
			}

			if(fieldp.nr_sol > 2) {
				for(int j = 0; j < fieldp.nreq; ++j) {
				 is >> fieldp.dof_ents[tmp].vec_dof_3[j];
				}
			}
#ifdef FV_DEBUG
			tmp1 = tmp;
#endif
		} // end while

#ifdef FV_DEBUG
		std::cout << "Number of constrain nodes: " << counter << std::endl;
		// chceking
		//FILE* fp = fopen("constrains.txt","w");
		//int cnt = 1;
		//for(int i=0;i<fieldp.nr_nodes;i++) {
		//	while(fieldp.dof_ents[cnt].vec_dof_1 == NULL) {
		//		cnt++;
		//	}
		//
		//	fprintf(fp,"pos= %d\t\tnode= %d\t\t value= %f\n",cnt,(i+1),fieldp.dof_ents[cnt].vec_dof_1[0]);
		//}
		//fclose(fp);
		//getchar();
#endif

		// After reading check if there are constrain nodes, and if
		// correct mesh data
		ibimpr.PostProcessing();		

 	
		return( true);


	}

	

	
}

#undef APC_TRUE
#undef APC_FALSE
