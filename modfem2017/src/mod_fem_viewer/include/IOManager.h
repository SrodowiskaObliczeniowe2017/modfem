#ifdef _WIN32
#pragma once
#pragma warning (disable: 4996)
#endif

#ifndef _MRCH_CUBE_H_
#define _MRCH_CUBE_H_

#include <stdio.h>
#include <stdlib.h>
#include "../utils/fv_dictstr.h"


#if defined(__cplusplus)
extern "C" {
#endif

#ifdef _WIN32
#  ifdef IOManager_EXPORTS
#    define DLL_MAPPING  __declspec(dllexport)
#  else
#    define DLL_MAPPING  __declspec(dllimport)
#  endif
#else 
#  define DLL_MAPPING 
#endif

#ifdef _WIN32
    #define FV_API __stdcall
#else
    #define FV_API
#endif

	typedef enum FV_Error_enum {
		FV_SUCCESS			    = 0,
		FV_FAILED				= -1,
		FV_INVALID_TYPE			= 1,
		FV_INVALID_FILE_PATH,
		FV_CANT_OPEN_FILE,
		FV_ERROR_FILE_NOT_FOUND,
	} FV_Result;

	typedef enum FV_IO_Items {
		FV_NONE				=  0x000,
		FV_MESH_PRISM		=  0x001,
		FV_MESH_TETRA		=  0x002,
		FV_MESH		    	=  FV_MESH_PRISM | FV_MESH_TETRA,
		FV_FLD_DG			=  0x010,
		FV_FLD_STD			=  0x020,
		FV_FIELD			=  FV_FLD_DG | FV_FLD_STD,
		FV_INIT_MESH		=  0x100,
		FV_INIT_FIELD_DG	=  0x200,
		FV_INIT_FIELD_STD	=  0x400,

	} IO_Type;

#define IS_MESH(w)			((w) & FV_MESH)
#define IS_FIELD(w)			((w) & FV_FIELD)
#define IS_DG_FIELD(w)		((w) & FV_FLD_DG)
	

	// this should init filemanager and etc
	DLL_MAPPING FV_Result FV_API
	IOMgr_Initialize(const int& argc, const char** argv);

	// this deactivate 
	DLL_MAPPING FV_Result FV_API
	IOMgr_Destroy();

	// this check file and return its type
	DLL_MAPPING FV_Result FV_API
	IOMgr_FileExists(const char* fname, const StrDict* map, const unsigned int size, IO_Type* type);

	// get the data from file
	DLL_MAPPING FV_Result FV_API
	IOMgr_ReadFromFile(const char* fname, void* pItem, IO_Type type);


#if defined(__cplusplus)
}
#endif

#endif /* _MRCH_CUBE_H_
*/
