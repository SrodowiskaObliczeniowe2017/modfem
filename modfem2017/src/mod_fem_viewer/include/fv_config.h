#ifndef _FV_CONFIG_H_
#define _FV_CONFIG_H_

#ifndef __OPENCL_VERSION__
#include <CL/cl.h>
#include <CL/cl_platform.h>
typedef cl_float3 float3;
typedef cl_float2 float2;
typedef cl_int3 int3;
typedef cl_int2 int2;

#endif
#include "types.h"

#ifndef NDEBUG
#elif defined(_DEBUG) || defined(DEBUG) || defined(__WXDEBUG__)
#  ifndef FV_DEBUG
#   define FV_DEBUG 1
#  endif
//# endif
#endif

// If the viewer must be statically linked
#ifdef MOD_FEM_VIEWER
# ifndef _USE_FV_LIB
#   define _USE_FV_LIB	1
# endif
#endif

// Allow to use IOManager
#ifndef _USE_FV_LIB
#  ifndef _USE_FV_IO_MGR
#    define _USE_FV_IO_MGR 1
#  endif
#endif

// Allow to use external modules
#ifndef _USE_FV_LIB
#  ifndef _USE_FV_EXT_MOD
#    define _USE_FV_EXT_MOD 1
#  endif
#endif

// some defines
/*
#ifndef APC_MAXELP_COMP 
#  define APC_MAXELP_COMP 7
#endif

#ifndef APC_MAXELP_TENS
#  define APC_MAXELP_TENS 707
#endif
*/

#ifndef _USE_FV_LIB
extern const char* dgDllName;
extern const char* stdDllName;
extern const char* stdDll2dName;
extern const char* stdhDllName;
#endif

/* Graphics module type */
//#define _USE_FV_WX 1
const int imageWidth  = 512;
#define DFLT_IMAGE_WIDTH 	512
const int imageHeight = 512;
#define DFLT_IMAGE_HEIGHT	512
// GUI interface
//#if defined(WIN32) && defined(VIEWER_WIN_GUI)
//# include "win/Common.h"
//# include "win/WinApp.h"
//#elif defined(VIEWER_FGLT_GUI)
//# include "glt/GlutWindow.h"
//#elif defined(VIEWER_WX_GUI)
//# include "wx/wxFemViewerApp.h"
//#else
//# error "Unknown GUI-type"
//#endif


/* Define application name */
#define MOD_FEM_VIEWER_NAME "FINITE ELEMENT MODEL VIEWER"
/* Define MAX_LEVEL_LOG */
#define MAX_LEVEL_LOG 	LogDEBUG
/* Define the threshold for the number of elements in mesh for parallel proceed
 */
#define NUM_ELEM_THRESHOLD	2
#define MAX_NUM_DOFS 128

#define IMG_WIDTH_CONSTRAINT	4
#define IMG_HEIGHT_CONSTRAINT	1

extern size_t tilesize[2];
//extern size_t lwsize[2];
extern frameconfig_s * frames_p;
extern frameconfig_s * current_frame_p;
extern frameconfig_s last_frame_s;
extern int framecount;
extern cl_int lightcount;
extern cl_int2 tiles;
extern struct light_s * lightlist;
extern float density;
extern char * outname;
extern char * scenefile;
extern int accelstruct;
extern int interactive;

extern const int nr_frames;
extern const bool initgl;
extern int img_width;
extern int img_height;
extern int gridcount;
extern int cellcount;

extern float cut_planes[3][4];
extern int current_plane_idx;
extern float dg_density;
extern int test_plane_cut;
extern float fixed_light_pos[3];

extern float bkg_color[3];
extern float node_id_color[3];
extern float elem_id_color[3];




#endif



