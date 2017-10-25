#ifndef _FV_COMPILER_H_
#define _FV_COMPILER_H_
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include "fv_config.h"

#ifdef WIN32
# define WINDOWS_LEAN_AND_MEAN
# include<windows.h>
# include <crtdbg.h>
#else
#include<sys/types.h>
#include<dlfcn.h>
#endif

#ifdef WIN32
# ifndef fun_call
#   define fun_call _cdecl
# endif
# ifndef std_call
#   define std_call _stdcall
# endif
#else
# define fun_call
# define std_call
#endif

#if defined(WIN32) && defined(_MSC_VER)
# pragma warning(disable: 4996)
#endif

#include<cassert>
#include<stdexcept>

typedef unsigned int  	uint_t;
typedef unsigned char 	ubyte_t;
typedef unsigned char	uint8_t;
typedef long int	  	long_t;
typedef unsigned long 	ulong_t;

typedef unsigned int	id_t;
typedef int				size_type;

#ifndef __uint32_t_defined
typedef unsigned int		uint32_t;
# define __uint32_t_defined
#endif

#ifdef _USE_FV_LIB
# ifndef fv_const_t
#  	define fv_const_t const
# endif
#else
# define fv_const_t
#endif

// Typedefs of function pointers
typedef int  (fun_call *fv_const_t intfint)(int);
typedef int  (fun_call *fv_const_t intfcharp)(char*);
typedef int  (fun_call *fv_const_t intf2int)(int,int);
typedef int  (fun_call *fv_const_t intf2intintp)(int,int,int*);
typedef int  (fun_call *fv_const_t intf2int2intp)(int,int,int*,int*);
typedef int  (fun_call *fv_const_t intf2intintpdoublep)(int,int,int*,double*);
typedef void (fun_call *fv_const_t voidf2int4intp2doublep)(int,int,
	          int*,int*,int*,int*,double*,double*);
typedef int  (fun_call *fv_const_t intf2intdoublep)(int,int,double*);
typedef int  (fun_call *fv_const_t intfintcharp)(int,char*);
typedef int  (fun_call *fv_const_t intf3intdoublep)(int,int,int,double*);
typedef void (fun_call *fv_const_t voidf2int3intp)(int,int,int*,int*,int*);

typedef int  (fun_call *fv_const_t int_global_solutiom)(int,double*,int,int*,double*,double*,double*,double*,double*,int);
typedef double (fun_call *fv_const_t double_el_calc)(int,int,int*,int,
				double*,double*,double*,double*,double*,double*,double*,
				double*,double*,double*,double*,double*,double*);
typedef int (fun_call *fv_const_t intfchar5intcharpfnp)(char,int,int,int,int,int,char*,double(*)(int,double*,int));
typedef int (fun_call *fv_const_t intfintcharpfilep)(int,char*,FILE*);


#define FV_STATIC_ASSERT(Type,Num)	\
	static_assert( sizeof(Type) == sizeof(int) * Num \
			/*&&  std::is_standard_layout<Type>::value*/, \
			"Type does not satisfy contiguous storage requirements")

#endif /* _FV_COMPILER_H_
*/
