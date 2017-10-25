#ifndef _FV_TOKENIS_H_
#define _FV_TOKENIS_H_

#include <iostream>
#include "fv_token.h"

class Token;

class tokenis
{
public:
	tokenis(std::istream& isToAttach);
	int get();                  
	tokenis& getline(signed char* buf, int size, char c= '\n');
	tokenis& getline(unsigned char* buf, int size, char c= '\n');
	tokenis& nextline();
	long getCounter()
	{ return lineCounter; }
	operator std::istream&()
	{ return is; }
	operator int()
	{ return is.good(); }
	void clear(int state= 0)
	{
#ifdef WIN32
		is.clear(state);
#else
		is.clear(static_cast<std::_Ios_Iostate>(state));
#endif
	}
	int eof()const{return is.eof();}
	int peek()
	{ return is.peek(); }
	int good() const
	{ return is.good(); }
	Token lastToken;
	friend tokenis& putback(tokenis& is,Token& aToken);
	long tokenCounter;
	long lineCounter;

protected:
	std::istream& is;

private:
	void init()
	{ lineCounter = 1; tokenCounter = 0; lastToken = ""; }
};


#endif

