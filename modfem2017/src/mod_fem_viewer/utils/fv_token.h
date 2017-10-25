#ifndef _FV_TOKEN_H_
#define _FV_TOKEN_H_

#include <istream>
#include <string>

class tokenis;

class Token : public std::string
{
public:
	Token() : std::string(), quoted(0) 
	{}
  
	Token(const char* a) : std::string(a), quoted(0) 
	{}
 
	Token& operator =(const Token& b){
		std::string::operator =(b);
		quoted = b.quoted;
		return *this;
	}

	friend std::istream& operator >>(std::istream& is,Token& aToken);
	friend tokenis& operator >>(tokenis&is,Token& aToken);
	int quoted;
	
	// static members and functions
	static const char* whiteSpaces;
	static const char* predefinedTokens;
	static int   commentChar;
	static int   stringChar;
	static void  defaultSyntax();
};

class lfToken : public Token
{
  public:
	  lfToken() : Token(){}
	  friend tokenis& operator >>(tokenis& is,lfToken& aToken);

	 void copy(std::string str){
		std::string::operator =(str);
	 }
};

class CppToken : public Token
{
public:
	CppToken& operator =(const Token& b){
		Token::operator =(b);
		return *this;
	}

	friend std::istream& operator >>(std::istream& is,CppToken& aToken);
	friend tokenis& operator >>(tokenis& is,CppToken& aToken);
};



#endif /* _FV_TOKEN_H_ 
*/
