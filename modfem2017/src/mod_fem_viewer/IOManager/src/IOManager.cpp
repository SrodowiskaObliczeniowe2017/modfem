#include "../../include/IOManager.h"
#include "../../include/BaseMesh.h"
#include "../../include/BaseField.h"

#include "../inc/FileManager.h"
#include "../inc/MeshImporter.h"
#include "../inc/FieldImporter.h"

#include<iostream>
#include<stdexcept>

namespace IOmgr {

	bool Init(const int argc,const char** argv)
	{
		FileManager::Init(argc,argv);
		return true;
	}

	bool Destroy()
	{
		FileManager::Destroy();
		return true;
	}

	template<class TImporter, class TItem>
	static bool Read(const char* fname, TItem* pItem)
	{
		TImporter importer(pItem);
		return FileManager::GetInstance()->ReadFile(fname, importer); 
	}
		
}



FV_Result FV_API
IOMgr_Initialize(const int& argc, const char** argv)
{
	try
	{
		if ( ! IOmgr::Init( argc,argv)) throw (std::runtime_error("Error while initializing FielManager!!!\n"));
	}
	catch(std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return FV_FAILED;
	}

	return FV_SUCCESS;
}

FV_Result FV_API
IOMgr_Destroy()
{
	if(!IOmgr::Destroy()) return FV_FAILED;
	return FV_SUCCESS;
}

FV_Result FV_API
IOMgr_FileExists(const char* fname, const StrDict* map, const unsigned int size, IO_Type* type)
{
	if (type!=NULL) *type = FV_NONE;
	if (fname == NULL) return FV_INVALID_FILE_PATH;

	// Recognize file by its name and then check if it can be opened
	if (type != NULL && size > 0){
		for(unsigned int i=0; i < size; ++i) {
			if(IOmgr::FileManager::GetInstance()->CheckInitials(fname, map[i].string)) {
				*type = static_cast<IO_Type>(map[i].key);

				break;
			}
		}
	}
	return IOmgr::FileManager::GetInstance()->IsFileExists(fname) ? FV_SUCCESS : FV_INVALID_FILE_PATH;
	
}



FV_Result FV_API
IOMgr_ReadFromFile(const char* fname, void *pItem, IO_Type type)
{
	using namespace FemViewer;
	if(!IOmgr::FileManager::GetInstance()->IsFileExists(fname))
		return FV_INVALID_FILE_PATH;

	bool res;
	switch(type)
	{
	case FV_MESH_PRISM:
		res = IOmgr::Read<IOmgr::MeshImporter,BaseMesh>(fname, reinterpret_cast<BaseMesh*>(pItem));
		break;
	case FV_FLD_DG:
	case FV_FLD_STD:
		res = IOmgr::Read<IOmgr::FieldImporter,BaseField>(fname, reinterpret_cast<BaseField*>(pItem));
		break;
	default:
		return FV_INVALID_TYPE;
		break;
	}

	if(!res) return FV_FAILED;
	else return FV_SUCCESS;
}
