#include "../inc/io_config.h"
#include "../inc/MeshReader.h"
#include "../inc/FileManager.h"


#include "../../include/IOManager.h"
#include "../../include/MeshElems.h"
#include "../../include/Enums.h"
#include "../../utils/fv_assert.h"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>


static int *mmr_ivector(const int &nr, const std::string& text)
{
	int *v = (int*)malloc(nr*sizeof(int));
	if(!v){
		std::cerr << "Not enough space for allocating vector: " << text
			<< std::endl;
		return( NULL);
	}
	return( v);
}

namespace IOmgr {


	FMVReader implFMVReader;// = NULL;

	FMVReader& InstanceOfFMReader()
	{
		return  implFMVReader;
	}

	FMVReader::FMVReader() : _name("MeshReader"), _initials("mesh")
	{
#ifdef FV_DEBUG
		std::cout << _name << " alive.\n";
#endif
		FileManager::GetInstance()->RegisterModules(this);
	}

	bool 
	FMVReader::Read(const char* fname, BaseImporter& ibimpr)
	{
		std::ifstream ifs(fname,std::ios::in);
		if(!ifs.is_open() || !ifs.good()){
			std::cerr << "MeshReader: can't open file: " << fname << std::endl;
			return false;
		}
#ifdef FV_DEBUG
		std::cout << "Reading file: " << fname << std::endl; 
#endif

		bool res = Read(ifs,ibimpr);
		ifs.close();
		return res;
	}

	bool 
	FMVReader::IsReadAble(const char *fname) const
	{
		if(BaseFileReader::IsReadAble(fname)){
			std::ifstream ifs(fname);
			if(ifs.is_open() && CheckInitials(fname)){
				ifs.close();
				return true;
			}
			ifs.close();
		}
		return false;
	}

	bool 
	FMVReader::Read(std::istream& is, BaseImporter& ibimpr)
	{
#ifdef FV_DEBUG
		std::cout << "Reading mesh.\n"; 
#endif
		if(!is.good()){
			std::cerr << "MeshReader: can't read from stream!" << std::endl;
			return false;
		}

		// Local structure of mesh params
		static mmt_meshp mshp;

		// Read dimensions for mesh arrays
		is >> mshp.mxno; // maximal number of nodes
		is >> mshp.mxed; // maximal numder of edges
		is >> mshp.mxfa; // maximal number of faces
		is >> mshp.mxel; // maximal number of elements
		
		// Read the number of stored node structures and the pointer to 
		// the first free space 
		is >> mshp.nmno;
		is >> mshp.pfno;

		// Adjust the number of allocated structures
		if(mshp.mxno < mshp.nmno)
			mshp.mxno = mshp.nmno;

		// Allocate space for nodes' structures - only for maximal initialized
		mmt_nodes* pnode = static_cast<mmt_nodes*>(ibimpr.Reserve(FemViewer::NODE,mshp.nmno+1));
		pnode++;
       
		// Read nodes structures
		int iaux = 0;
		for(int i=0;i < mshp.nmno; ++i, ++pnode){

			// Set interprocessor ID to zero
			//pnode->ipid = 0;

			// read X-coor
			is >> pnode->x;
			if(pnode->x < -1e10){
				is >> pnode->y;
			} else {
				is >> pnode->y;
				is >> pnode->z;
				iaux++;
			}
		}

		// Set the number of active nodes
		mshp.nrno = iaux;

#ifdef FV_DEBUG
		// Report 
		std::cout << "Mesh:\n";
		std::cout << "Nodes   : allocated " << mshp.nmno << " structures, read "
			<< mshp.nmno << " structures for " << mshp.nrno << " active nodes\n";
#endif
		// Read the number of stored edge structures and the pointer to
		// first free space
		is >> mshp.nmed;
		is >> mshp.pfed;
   
		// Adjust the number of allocated structures
		if(mshp.mxed < mshp.nmed) 
			mshp.mxed = mshp.nmed;

		// Allocate space for nodes' structures - only for maximal initialized
		mmt_edges* pedge = static_cast<mmt_edges*>(ibimpr.Reserve(FemViewer::EDGE,mshp.nmed+1));
		pedge++;

		// Read edges structures
		iaux = 0;
		for(int i=0;i<mshp.nmed;++i, ++pedge){

			// Set interprocessor ID to zero
			//pedge->ipid = 0;

			is >> pedge->type;
			is >> pedge->node[0]; is >> pedge->node[1];
			if(pedge->type>0) iaux++;
		}

		// Set the number of acitve adges
		mshp.nred = iaux;
#ifdef FV_DEBUG
		// Report
		std::cout << "Edges   : allocated " << mshp.nmed << " structures, read " 
			<< mshp.nmed << " structures for " << mshp.nred << "edges\n";
#endif
		// Read the number of stored face structures and the pointer to 
		// first free space
		is >> mshp.nmfa;
		is >> mshp.pffa;

		// Adjust the number of allocated structures
		if(mshp.mxfa < mshp.nmfa) 
			mshp.mxfa = mshp.nmfa;

		// Allocate space for faces' structures - only for maximal initialized
		mmt_faces* pface = static_cast<mmt_faces*>(ibimpr.Reserve(FemViewer::FACE,mshp.nmfa+1));
		pface++;

		// Read face structures
		iaux = 0;
		for(int i = 0; i < mshp.nmfa; ++i, ++pface){

			// Set interprocessor ID to zero 
			//pface->ipid = 0;
			
			is >> pface->type;
			is >> pface->bc;
			is >> pface->neig[0]; is >> pface->neig[1];

			if(abs(pface->type) == FemViewer::F_TRIA){
				is >> pface->edge[0]; is >> pface->edge[1];
				is >> pface->edge[2];
			}
			else if(abs(pface->type) == FemViewer::F_QUAD){
				is >> pface->edge[0]; is >> pface->edge[1];
				is >> pface->edge[2]; is >> pface->edge[3];
			}

			// Initialize pointer to sons
			pface->sons = NULL;

			if(pface->type>0) {
				iaux++;
			}
			// Read sons for inactive faces 
			else if(pface->type<0) {
				const int ison = 4;
				pface->sons = mmr_ivector(ison,"face sons in read data");
				if(!pface->sons) {
					std::cout << "Can't allocate memory for sons elements!!\n";
					exit(-1);
				}
				//if(!pface->sons) throw( fv_exception("Can't allocate memory for sons in reading mesh procedure!"));
				for(int j=0;j<ison;j++){
					//int tmp; is >> tmp;
					is >> pface->sons[j];
				}
			}
		} // end for i

		// Set the number of active faces
		mshp.nrfa = iaux;
#ifdef FV_DEBUG
		// Report
		std::cout << "Faces   : allocated " << mshp.nmfa << " structures, read "
			<< mshp.nmfa << " structures for " << mshp.nrfa << " faces\n";
#endif
		// Read the number of stored element structures and the pointer to 
		// first free space
		is >> mshp.nmel; is >> mshp.pfel;
    
		// Adjust the number of allocated structures
		if(mshp.mxel < mshp.nmel) 
			mshp.mxel = mshp.nmel;

		// Allocate space for elements' structures - only for maximal initialized
		mmt_elems* pelem = static_cast<mmt_elems*>(ibimpr.Reserve(FemViewer::ELEM,mshp.nmel+1));
		pelem++; // increment to first element
		// Read element structures
		iaux = 0;
		for(int i = 0; i < mshp.nmel; ++i, ++pelem) {

			// Set interprocessor ID to zero
			//pelem->ipid = 0;

			is >> pelem->type;
			is >> pelem->mate;
			is >> pelem->fath;
			is >> pelem->refi;

#ifdef FV_DEBUG
			if(pelem->type == 0) {
				std::cout<< "Warning: element= " << (i+1) << " type=0, mate= "<< pelem->mate << " fath= " << pelem->fath
					<< " refi= " << pelem->refi << std::endl;
			}
#endif

			// Initialize pointer to sons
			//pelem->face = NULL;
			pelem->sons = NULL;

			// For active and inactive elements
			if(abs(pelem->type) == FemViewer::TETRA){
			}
			else if(abs(pelem->type) == FemViewer::PRIZM){
				const int ii = 5;
				pelem->face = mmr_ivector(ii,"element faces in read data");
				is >> pelem->face[0]; is >> pelem->face[1];
				is >> pelem->face[2]; is >> pelem->face[3];
				is >> pelem->face[4];
			}
			else if(abs(pelem->type) == FemViewer::BRICK){
			}

			// If element not a free space
			if(pelem->type != FemViewer::FREE){
				// For inactive elements
				if(pelem->type<0){
					// Refinement kind
					int ison;
					if(pelem->refi == FemViewer::REF_ISO) {
						ison = 8;
					} else {
						std::cerr << "Unknown refinement type for inactive element "
							<< i << std::endl;
						return( false);
					}
					pelem->sons = mmr_ivector(ison,"sons in read data");
					for(int j=0; j < ison; ++j){
						is >> pelem->sons[j];
					}
				}
				//fscanf(fp,"\n");
				else {
					// For active elements
					iaux++;
#ifdef FV_DEBUG
					if(pelem->refi != FemViewer::NOT_REF) {
						std::cerr << "Refinement type indicated for active element "
							<< i << std::endl;
						return( false);
					}
#endif
				} // else
			} // ednf if type == FREE
		} // endfor

		// Set the number of active elements
		mshp.nrel = iaux;

		// Copy params to mesh object
		ibimpr.ImportParams(&mshp);
#ifdef FV_DEBUG
		// Report 
		std::cout << "Elements: allocated " << mshp.nmel << " structures, read "
			<< mshp.nmel << " structures for " << mshp.nrel << " elements\n";
#endif
		return( true);
	
	}

} // end namespace
