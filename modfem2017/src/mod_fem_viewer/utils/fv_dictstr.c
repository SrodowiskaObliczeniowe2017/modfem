#include "fv_dictstr.h"
#include <string.h>

#ifdef WIN32
# if defined(_MSC_VER) && (_MSC_VER >= 1400)
#	pragma warning(disable: 4996)
# endif
#endif

const char*
getString(const StrDict *table,unsigned int size,int key)
{
	unsigned int i;
	for(i = 0; i < size; i++)
		if(table[i].key == key)return table[i].string;
	return((const char*)0);
}

int 
getStrKey(const StrDict *table,unsigned int size,const char* str)
{
	unsigned int i;
	for(i = 0; i < size; i++)
	if(!strcmp(table[i].string,str)) return table[i].key;
	return(-1);
}

const char*
getDescription(const ErrDict *table,unsigned int size,int key)
{
	unsigned int i;
	for(i = 0; i < size; i++)
		if(table[i].key == key) return table[i].description;
	return((const char*)0);
}

int  
getKey(const ErrDict *table,unsigned int size,const char* str)
{
	unsigned int i;
	for(i=0;i<size;i++)
		if(!strcmp(table[i].description,str))return table[i].key;
	return(-1);
}

unsigned char 
getErrLevel(const ErrDict *table,unsigned int size,int key)
{
	unsigned int i;
	for(i=0;i<size;i++)
		if(table[i].key==key)return table[i].errlevel;
	return((char)0);
}

