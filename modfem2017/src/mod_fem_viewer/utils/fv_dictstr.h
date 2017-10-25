#if !defined(_DICTSTR_H_)
#define _DICTSTR_H_


#ifdef __cplusplus 
extern "C" {
#endif

#define FV_SIZEOF_ARRAY(arr) 	(sizeof(arr)/sizeof((arr)[0]))

	typedef struct tagChrDict{
		int  key;
		char chr;
	} ChrDict;
	
	typedef struct tagStrDict{
		int   key;
		const char *string;
	} StrDict;

	extern const char* getString(const StrDict* table,unsigned int size,int key);
	extern		 int   getStrKey(const StrDict *table,unsigned int size,const char *str);

	/* DICTIONARY FOR ERROR NAMES AND LEVELS */
	typedef struct tagErrDict{
		int key;
		unsigned char errlevel;
		char* description;
	} ErrDict;

	extern const char*   getDescription(const ErrDict *table,unsigned int size,int key);
	extern       int     getKey(const ErrDict *table,unsigned int size,const char *str);
	extern unsigned char getErrLevel(const ErrDict *table,unsigned int size,int key);

	typedef unsigned long int ulint_t;


#ifdef __cplusplus 
}
#endif

#endif /* _FV_DICTSTR_H_ 
*/
