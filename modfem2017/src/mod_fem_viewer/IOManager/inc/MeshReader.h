#ifndef _MESH_READER_H_
#define _MESH_READER_H_

#include <string>

#include "BaseReader.h"
#include "BaseImporter.h"


namespace IOmgr {

	class FMVReader : public BaseFileReader {

	public:
		FMVReader();
		virtual ~FMVReader() {}
		const char* GetName() const { return _name; }
		const char* GetFileExtension() const { return "dmp dat"; }
		const char* GetInitilas() const { return "mesh"; }
		const char* GetDescription() const { return "Standard FEM mesh file format"; }
		
		bool Read(const char* fname, BaseImporter& ibimpr);
		bool IsReadAble(const char* fname) const;
	private:
		const char* _name;
		const char* _initials;
		bool Read(std::istream& is, BaseImporter& ibimpr);
		inline bool CheckInitials(const char* fname) const;
	};

	inline bool FMVReader::CheckInitials(const char *fname) const
	{
		return BaseFileReader::CheckInitials(fname,_initials);
	}

	extern FMVReader implFMVReader;

	FMVReader& InstanceOfFMReader();


}

#endif /* _MESH_READER_H_
*/