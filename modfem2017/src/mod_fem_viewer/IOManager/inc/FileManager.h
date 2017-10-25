#ifndef _FILE_MANGER_H_
#define _FILE_MANGER_H_

#include <vector>

#include "BaseImporter.h"
#include "BaseReader.h"

namespace IOmgr {


	class FileManager {
	protected:
		static FileManager* instance;
	public:
		static void Init(const int argc, const char** argv);
		static FileManager* GetInstance();
		static void Destroy();
		bool ReadFile(const char* fname, BaseImporter& bimpr);
		void RegisterModules(BaseFileReader* bfr);
		bool IsFileExists(const char* fname) const;
		bool CheckInitials(const char* fname,const char* initials) 
		{
			return BaseFileReader::CheckInitials(fname, initials);
		}
		
	private:
		static int rargc;
		static char** rargv;
		std::vector<BaseFileReader*> readers;
	private:
		//std::set<float> reader;
		FileManager(); 
		FileManager(const FileManager&);
		FileManager& operator=(FileManager&);
		
	};

}


#endif /* _FILE_MANGER_H_
*/