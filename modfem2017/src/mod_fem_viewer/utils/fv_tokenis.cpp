#include "fv_tokenis.h"
#include <cstdio>

using namespace std;

tokenis::tokenis(std::istream&isToAttach):is(isToAttach)
{
	init();
}

int 
tokenis::get()
{
	if(!is) return *this;
	int cr = is.get();
	if(cr == '\n') lineCounter++;
	return cr;
}

tokenis& 
tokenis::getline(signed char* buf, int size, char c)
{
	if(!is) return *this;
	is.getline((char *)buf, size, c);
	if(c=='\n' && !is.eof()) lineCounter++;
	return *this;
}

tokenis&
tokenis::getline(unsigned char* buf, int size, char cr)
{
	if(!is) return *this;
	is.getline((char *)buf, size, cr);
	if(cr == '\n' && !is.eof()) lineCounter++;
	return *this;
}

tokenis& 
tokenis::nextline()
{
	int cr;
	while(1){
		cr = get();
		if(cr == '\n' || cr == EOF) break;
	}
	return *this;
}

tokenis&
putback(tokenis& is,Token& aToken)
{
	if(aToken.empty()) return is;
	is.lastToken = aToken;
	is.tokenCounter--;
	if(((std::istream&)is).eof()) ((std::istream&)is).clear();
	return is;
}


/* EOF */
