#include "fv_token.h"
#include "fv_tokenis.h"
#include <string.h>
#include <iostream>
#include <cstdio>

using namespace std;
#ifndef EOF
# define EOF	(-1)
#endif

#define MAXBUF 1024
#ifdef WIN32
#define VERBATIMCHAR	'\\'
#else
#define VERBATIMCHAR	'/'
#endif
#define STRINGCHAR		'\"'

const char* Token::whiteSpaces	  =	" \t\r\n";
const char* Token::predefinedTokens =	"(),{}[]=*/#:;";
int   Token::commentChar	  = '#';
int   Token::stringChar		  = STRINGCHAR;

void 
Token::defaultSyntax()
{
	whiteSpaces			= " \t\r\n";
	predefinedTokens	= "(),{}[]=*/#:;";
	commentChar			= '#';
	stringChar			= STRINGCHAR;
}


tokenis&
operator >>(tokenis& is,Token& aToken)
{
	if(is.lastToken != ""){
		aToken = is.lastToken;
		is.lastToken = "";
		is.tokenCounter++;
		return is;
	}

	int cr, i = 0; 
	bool quoted = false, ready = false, isToken = false, strMode = false;
	char buf[MAXBUF];

	while(!ready){
		if(is.eof()) break;
		cr = is.peek();
		if(cr == EOF){
			is.clear(ios::eofbit);
			break;
		}
#if 0
	   if(quoted &&  c=='\\'){
			c=is.get();
			int c2=is.peek();
			if(c2==EOF){is.clear(ios::eofbit);break;}
			if(c2=='\\' || c2=='\"'){
				isToken=1;
				buf[i++]=(char)c2;
				c2=is.get();
			}else{
				isToken=1;
				buf[i++]=(char)c;
			}
			continue;
	   }
#endif

	if(isToken)
	{
		if(!strMode){
			if(cr == Token::commentChar){
				ready = !ready;   // skip to EOL or EOF
				while(1){
					cr = is.get();
					if(cr == '\n' || cr == EOF) break;
				}
			}
			if(strchr(Token::whiteSpaces,cr) || strchr(Token::predefinedTokens,cr)){
				ready = !ready;
				break;
			}
		} else {
			if(cr == Token::stringChar){
				is.get();
				ready = !ready;
				break;
			}
		}
		buf[i++] = (char)cr;
	}
	else{
		if(cr == Token::commentChar){
			// skip to EOL or EOF
			while(1){
				cr = is.get();
				if(cr == '\n' || cr == EOF) break;
			}
			continue;
		}
		if(strchr(Token::whiteSpaces,cr)){
			is.get();
			continue;
		}
		if(strchr(Token::predefinedTokens,cr)) ready = !ready;
		isToken = !isToken;
		if(cr == Token::stringChar){
			strMode = !strMode;
			quoted = !quoted;
			cr = is.get();
			continue;
		}
		buf[i++] =(char)cr;
	}
	cr = is.get();
	}
	buf[i] = 0;
	aToken = buf;
	aToken.quoted = quoted;
	is.tokenCounter++;
	return is;
}



tokenis& 
operator >>(tokenis& is, lfToken& aToken)
{
	char buf[MAXBUF];
	is.getline((signed char*)buf,MAXBUF);
	int cr;
	cr = is.peek();
	if(cr ==' \n') 
		is.get();
	aToken.copy(buf);
	return is;
}


std::istream&
operator >>(std::istream& is,Token& aToken)
{
 
	int cr, i = 0; 
	bool ready = false, isToken = false, strMode = false;
	char buf[MAXBUF];

	while(!ready){
		if(is.eof()) break;
		cr = is.peek();
		if(cr == EOF){
			is.clear(ios::eofbit);
			break;
		}

		if(isToken){
			if(!strMode){
				if(cr == Token::commentChar){
					ready = !ready;   // skip to EOL or EOF
					while(1){
						cr = is.get();
						if(cr == '\n' || cr == EOF) break;
					}
				}
				if(strchr(Token::whiteSpaces,cr) || strchr(Token::predefinedTokens,cr)){
					ready = !ready;
					break;
				}
			}else{
				if(cr == Token::stringChar){
					is.get();
					ready = !ready;
					break;
				}
			}
			buf[i++] = (char)cr;
		}
		else{
			if(cr == Token::commentChar){
				// skip to EOL or EOF
				while(1){
					cr = is.get();
					if(cr == '\n' || cr == EOF) break;
				}
				continue;
			}
			if(strchr(Token::whiteSpaces,cr)){
				is.get();
				continue;
			}
			if(strchr(Token::predefinedTokens,cr)) 
				ready = !ready;
			isToken = !isToken;
			if(cr==Token::stringChar){
				strMode = !strMode;
				cr = is.get();
				continue;
			}
			buf[i++] =(char)cr;
		}
		cr = is.get();
	}
	buf[i] = 0;
	aToken = buf;
 
	return is;
}

tokenis&
operator >>(tokenis& is,CppToken& aToken)
{
 if(is.lastToken!=""){
	aToken=is.lastToken;
	is.lastToken="";
	is.tokenCounter++;
	return is;
 }
 int quoted=0;
 //
 char buf[MAXBUF];
 int i=0;
 int c;
 int ready=0;
 int isToken=0;
 int stringMode=0;

 while(!ready)
 {
   if (is.eof()) break;
   c = is.peek();
   if(c==EOF){
	is.clear(ios::eofbit);
	break;
   }
   
   if(quoted &&  c=='\\'){
		c=is.get();
		int c2=is.get();
		if(c2==EOF){is.clear(ios::eofbit);break;}
		switch(c2){
			case 'n':buf[i++]='\n';break;
			case 'r':buf[i++]='\r';break;
			case 't':buf[i++]='\t';break;
			case '0':buf[i++]='\0';break;
			default:buf[i++]=(char)c2;break;
		}
		continue;
   }

   if (isToken)
   {
	 if(!stringMode){
		if(c==Token::commentChar){
			ready=1;   // skip to EOL or EOF
			while(1){c=is. get();if(c=='\n' || c==EOF)break;}
		}
		if (strchr(Token::whiteSpaces,c) || strchr(Token::predefinedTokens,c)) {ready=1;break;}
	 } else {
		if(c==Token::stringChar){is.get();ready=1;break;}
	 }
	 buf[i++]=(char)c;
   }
   else
   {
	 if(c==Token::commentChar){
			// skip to EOL or EOF
			while(1){c=is. get();if(c=='\n' || c==EOF)break;}
			continue;
	 }
	 if (strchr(Token::whiteSpaces,c)) {is.get();continue;}
	 if (strchr(Token::predefinedTokens,c)) ready=1;
	 isToken = 1;
	 if(c==Token::stringChar){
		stringMode=1;
		quoted=1;
		c = is.get();
		continue;
	 }
	 buf[i++] =(char) c;
   }
   c = is.get();
 }
 buf[i] = 0;
 aToken=buf;
 aToken.quoted=quoted;
 is.tokenCounter++;
 return is;
}



std::istream&
operator >>(std::istream& is,CppToken& aToken)
{
 //
 int quoted=0;
 //
 char buf[MAXBUF];
 int i=0;
 int c;
 int ready=0;
 int isToken=0;
 int stringMode=0;

 while(!ready)
 {
   if (is.eof()) break;
   c = is.peek();
   if(c==EOF){
	is.clear(ios::eofbit);
	break;
   }
   
   if(quoted &&  c=='\\'){
		c=is.get();
		int c2=is.get();
		if(c2==EOF){is.clear(ios::eofbit);break;}
		switch(c2){
			case 'n':buf[i++]='\n';break;
			case 'r':buf[i++]='\r';break;
			case 't':buf[i++]='\t';break;
			case '0':buf[i++]='\0';break;
			default:buf[i++]=(char)c2;break;
		}
		continue;
   }

   if (isToken)
   {
	 if(!stringMode){
		if(c==Token::commentChar){
			ready=1;   // skip to EOL or EOF
			while(1){c=is. get();if(c=='\n' || c==EOF)break;}
		}
		if (strchr(Token::whiteSpaces,c) || strchr(Token::predefinedTokens,c)) {ready=1;break;}
	 } else {
		if(c==Token::stringChar){is.get();ready=1;break;}
	 }
	 buf[i++]=(char)c;
   }
   else
   {
	 if(c==Token::commentChar){
			// skip to EOL or EOF
			while(1){c=is. get();if(c=='\n' || c==EOF)break;}
			continue;
	 }
	 if (strchr(Token::whiteSpaces,c)) {is.get();continue;}
	 if (strchr(Token::predefinedTokens,c)) ready=1;
	 isToken = 1;
	 if(c==Token::stringChar){
		stringMode=1;
		quoted=1;
		c = is.get();
		continue;
	 }
	 buf[i++] =(char) c;
   }
   c = is.get();
 }
 buf[i] = 0;
 aToken=buf;
 aToken.quoted=quoted;
 //
// is.tokenCounter++;
 return is;
}
