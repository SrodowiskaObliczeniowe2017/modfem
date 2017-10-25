#include <algorithm>
#include <string>

#include "../inc/BaseReader.h"
#include "../utils/fv_str_utils.h"

#ifdef WIN32
# define DIR_SEPAR '\\'
#else
# define DIR_SEPAR '/'
#endif


namespace IOmgr {

	inline char tolower(char c)
	{
		return static_cast<char>(std::tolower(c));
	}

	bool BaseFileReader::IsReadAble(const char *fname) const
	{
		std::string _fname(fname), _ext;
		std::string::size_type pos(_fname.rfind("."));
		if(pos != std::string::npos){ 
			_ext = _fname.substr(pos+1, _fname.length()-pos-1);
			std::transform(_ext.begin(), _ext.end(), _ext.begin(), tolower);
		}
		return(std::string(GetFileExtension()).find(_ext) != std::string::npos);
	}

	bool BaseFileReader::CompareExtensions(const char* fname, const char* ext) const
	{
		std::string _fname(fname), _ext(ext);
		std::string::size_type pos(_fname.rfind("."));

		if(pos != std::string::npos && !_ext.empty()) {
			std::string ext_( _fname.substr(pos+1, _fname.length()-pos-1) );
			std::transform(ext_.begin(), ext_.end(), ext_.begin(), tolower);
			std::transform(_ext.begin(), _ext.end(), _ext.begin(), tolower);
    		return( (_ext.find(ext_) != std::string::npos) ? true : false);
		}
		return false; 
	}

	bool BaseFileReader::CheckInitials(const char* path, const char* initials)
	{
		std::string _path(path);
		std::string _initials(initials);
		std::string _fname;
		std::vector<std::string> vtokens;

		// Extract a filename from a path
		size_t pos = _path.rfind(DIR_SEPAR);
		if(pos != std::string::npos)
			_fname = _path.substr(pos+1, _path.length()-pos-1);

		// Split initials
		int size = trimString(_initials,vtokens), cnt = 0;
		for(int i(0);i<size;++i) {
			// check if initials are presents
			pos = _fname.find(vtokens[i]);
			if (pos!=std::string::npos) cnt++;
		}

		return( (size == cnt) && (size > 0));
	}
}
