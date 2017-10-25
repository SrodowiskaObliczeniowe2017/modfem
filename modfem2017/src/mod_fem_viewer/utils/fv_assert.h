#ifndef _FV_ASSERT_H_
#define _FV_ASSERT_H_

#include <stdio.h>
#include "fv_config.h"
#include "defs.h"


#ifdef FV_DEBUG
#include <signal.h>


#ifndef WIN32
# define fv_assert(Condition,file,line) { \
		bool Cond = Condition; 			\
		if(!Cond){						\
			fprintf(stderr,"ASSERTION: %s feild in %s, line %d\n",FV_STR(Condition),file,line); \
			::raise(SIGSEGV); 			\
		} 								\
	}

# else
# define WINDOWS_LEAN_AND_MEAN
#include <windows.h>
inline void fv_assert(bool Condition, const char* file, const int line)
{ 
	if(!Condition) { 
		wchar_t buf[1024];
		LPTSTR pbuf = (LPTSTR)buf;
		wsprintf(pbuf,TEXT("Assertion failed in file: %s line: %d"),file,line);
		MessageBox(0,pbuf,TEXT("Error!"),MB_ICONSTOP);
		exit(EXIT_FAILURE);
	}
}
#endif

#define FV_ASSERT(Condition)	fv_assert(Condition,__FILE__,__LINE__)

#else
#define FV_ASSERT(Condition)
#endif // Fv_DEBUF





#endif // _FD_ASSERT_H_
