#ifndef _FIELD_READER_H_
#define _FIELD_READER_H_

#include "BaseReader.h"
#include<sstream>

namespace IOmgr {

	class FieldReaderDG : public BaseFileReader {
	public:
		// Constructor & destructor
		FieldReaderDG();
		virtual ~FieldReaderDG() {}
		
		// returns reader name
		const char* GetName() const 
		{ return( _name); }
		
		// returns file extension can read
		const char* GetFileExtension() const 
		{ return( _extension); }

		// returns reader initilas - from this the file is recognizable
		const char* GetInitials() const 
		{ return( _initials); }

		// returns reader desrciption
		const char* GetDescription() const 
		{ return( _description); }

		// read file
		bool Read(const char* fname, BaseImporter& ibimpr);

		// check if we can read a given file
		bool IsReadAble(const char *fname) const;
	
	private:
		// the name of the reader
		const char* _name;
		const char* _extension;
		const char* _initials;
		const char* _description;
	
	private:
		bool Read(std::istream& is, BaseImporter& ibimpr);

		inline bool CheckInitials(const char* fname) const;

	};

	inline bool FieldReaderDG::CheckInitials(const char *fname) const
	{
		return BaseFileReader::CheckInitials(fname,_initials);
	}

	extern FieldReaderDG inst_FieldReaderDG;
	FieldReaderDG& InstanceOfFieldReaderDG();
}

#endif /* _FIELD_READER_H_
*/
