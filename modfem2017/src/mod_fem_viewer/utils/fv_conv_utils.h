#ifndef _FV_CONV_UTL_H_
#define _FV_CONV_UTLS_H_

#include <ios>
#include <string>
#include <sstream>
#include "fv_exception.h"

template <typename Type>
inline
Type Str2T(std::string & str_)
{
	std::istringstream _isstr(str_);
	Type _value;
	_isstr >> _value;
	if(_isstr.fail())
		throw fv_exception("Error while parsing to type");
	return static_cast<Type>(_value);
}

template <typename Type>
inline
std::string T2Str(Type value_)
{
	std::stringstream _sstr;
	_sstr << value_;
	if(_sstr.fail())
		throw fv_exception("Error while parsing to string");
	return _sstr.str();
}

extern bool IsDig(std::string &str);
extern bool IsCyfr(std::string &str);

extern void int2bin(unsigned int nr,unsigned int size,char buff[]);
extern void  int2bin2(unsigned int n, unsigned int  buffer_size, char  *buffer);



#endif /* _FV_CONV_UTLS_H_
*/
