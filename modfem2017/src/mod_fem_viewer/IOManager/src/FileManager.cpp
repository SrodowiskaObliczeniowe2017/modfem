#include "../inc/BaseReader.h"
#include "../inc/FileManager.h"
#include "../inc/MeshReader.h"
#include "../include/fv_compiler.h"
#include "../../utils/fv_exception.h"
#include <iostream>
#include <fstream>
#include <string>


namespace IOmgr {

	FileManager* FileManager::instance = 0;
	char** FileManager::rargv = 0;
	int FileManager::rargc = 0;

	void FileManager::Init(const int argc, const char** argv)
	{
		if(instance != NULL) {
			return;
		}

		instance = new FileManager();

		rargc = argc;
		rargv = const_cast<char**>(argv);
	}

	FileManager* FileManager::GetInstance()
	{
		if(instance == NULL) FileManager::Init(0,0);
		return instance;
	}

	FileManager::FileManager() /*: readers()*/{
#ifdef FV_DEBUG
			std::cout << "Creating FileManger\n";
#endif
		}


	void 
	FileManager::Destroy()
	{
#ifdef FV_DEBUG
		std::cout << "Destroing FileManager\n";
#endif
		delete instance;
		instance = NULL;
	}

	bool 
	FileManager::ReadFile(const char* fname, BaseImporter& bimpr)
	{
		std::vector<BaseFileReader*>::const_iterator it     =  readers.begin();
		std::vector<BaseFileReader*>::const_iterator it_end =  readers.end();

		for(; it != it_end; ++it){
#ifdef FV_DEBUG
			std::cout << "Query module: " << (*it)->GetName();
#endif
			if((*it)->IsReadAble(fname)){
				return((*it)->Read(fname, bimpr));

			}
#ifdef FV_DEBUG
			std::cout << " no\n"; 
#endif


		}
		return false;
	}

	void 
	FileManager::RegisterModules(BaseFileReader* bfr_)
	{
		assert(bfr_ != 0);

		readers.push_back(bfr_);
#ifdef FV_DEBUG
		std::cout << "FileManager: register module: " << bfr_->GetName() << std::endl;
#endif
	}

	bool 
	FileManager::IsFileExists(const char* fname) const
	{
		bool res;
		std::fstream file(fname, std::ios_base::in);

		if(file.is_open() && file.good()) {
			res = true;
		} else {
			res = false;
		}

		file.close();

    
		return res;
	}


	



}
