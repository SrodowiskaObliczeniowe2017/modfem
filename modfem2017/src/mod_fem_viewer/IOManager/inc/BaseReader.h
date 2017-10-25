#ifndef _BASE_READER_H_
#define _BASE_READER_H_


namespace IOmgr {

	class BaseImporter;

	class BaseFileReader{
	public:
		/* Destructor
		**/
		virtual ~BaseFileReader() {}

		/* Returns the name of the reader
		**/
		virtual const char* GetName() const = 0;

		/* Returns file extension
		**/
		virtual const char* GetFileExtension() const = 0;

		/* Returns file description
		**/
		virtual const char* GetDescription() const = 0;

		/* Read the file
		**/
		virtual bool Read(const char* fname, BaseImporter& bimpr) = 0;

		/* Check if given file can be read
		**/
		virtual bool IsReadAble(const char* fname) const;

	protected:
		/* Find out if as given file has an extension that we can managed
		**/
		bool CompareExtensions(const char* fname, const char* ext) const;

		/* Extract the initials from filename and compare
		**/
	public:
		static bool CheckInitials(const char* path, const char *initials);// const;

	};


}

#endif /* _BASE_READER_H_
*/