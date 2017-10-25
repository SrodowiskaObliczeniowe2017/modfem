#ifndef _ENUMS_TYPE_H_
#define _ENUMS_TYPE_H_

#include <cstring> // memset
#include "MathHelper.h"
#include "../../include/fv_config.h"


#ifndef LARGE_F
#define LARGE_F 10e10
#endif

namespace FemViewer
{

	typedef enum {
		REFIN = -1,
		FREE  =  0,
		
		PRIZM =  5,
		BRICK =  6,
		TETRA =  7,
	} ElemType;

	typedef enum {
		F_TRIA = 3,
		F_QUAD = 4,
	} FaceType;

	/*typedef enum {
		NITM = -1,
		NODE =  0, 
		EDGE =  1,
		FACE =  2,
		ELEM =  3,
		FIELD = 4,
		FIELDN = 5,
	} ElemItem;*/

	typedef enum {
		NOT_REF = 0,
		REF_ISO = 1,
	} RefinementType;

	typedef enum {
		TENSOR   = 0,
		COMPLETE = 1,
	} BaseType;

	typedef enum {
		INTERNAL = 0,
		EXTERNAL,
	} ModuleType;

	typedef enum {
		LINEAR = 0,
		HIGH_ORDER,
	} ApproximationType;

	typedef enum {
		NUM_VERTICES_LIST_WIREFRAME = 0,
		NUM_ELEMENTS_LIST_WIREFRAME,
		NUM_VERTICES_LIST_COLORMAP,
		NUM_ELEMENTS_LIST_COLORMAP,
		NUM_VERTICES_LIST_CUTTED_WIREFRAME,
		NUM_ELEMENTS_LIST_CUTTED_WIREFRAME,
		NUM_VERTICES_LIST_CUTTED_COLORMAP,
		NUM_ELEMENTS_LIST_CUTTED_COLORMAP,
		COLORBAR_LIST0,
		COLORBAR_LIST1,
		COLORBAR_LIST2,
		ID_COLORMAP_TEXT,
		ID_GRID,
		ID_AXES,
		NUM_TOTAL_LISTS
	 } GLListType;

	typedef enum {
		RASTERIZATION_GL = 0,
		RAYTRACE_GL_CL,
	} Render_t;

	typedef enum {
		Unknown 	= -1,
		MeshPrizm 	= 0,
		MeshHybrid 	= 1,
		MeshRemesh  = 2,
		FieldSTD 	= 0,
		FieldDG		= 1,
	} HandleType;

	enum eSelectionCategory
		{
			All = 0,
			Boundary,
			Cutted,
		};

	enum {
		vtxTriangle = 0,
		vtxQuad,
		vtxEdge,
		vtxAll
	};

	enum MouseMode {
						MOUSE_NONE,
						MOUSE_ROTATE,
						MOUSE_TRANSLATE,
						MOUSE_ZOOM
			};

	typedef enum TypeOfRenderer {
	   	Wireframe = 0,
		WireframeSlice,
		ColorMap,
		ColorMapSlice,
		ColorMapStd,
		ColorMapStdSlice,
		ColorMapBar,

		Totall
	} object_type;

	typedef enum RenderGLTypes {
		WIREFRAME_GL,
		WIREFRAME_CUTS_GL,
		COLORMAP_GL,
		COLORMAP_CUTS_GL,
		COLORMAP_BAR,

		NUM_DRAWS,
	} gl_object_type;

	/* Reference prismatic element */
	const double XlocPrizm[18] = {
				0.0, 0.0, -1.0,
				1.0, 0.0, -1.0,
				0.0, 1.0, -1.0,
				0.0, 0.0,  1.0,
				1.0, 0.0,  1.0,
				0.0, 1.0,  1.0
	};

	/* Reference tetrahedron element */
	const double XlocTetra[12] = {
				1.0, 0.0, 0.0,
				0.0, 1.0, 0.0,
		        0.0, 0.0, 1.0,
				0.0, 0.0, 0.0,
	};

	struct isect_info_t {
		float t;	  // t parameter for ray
		float u, v;   // u,v barycentric coordinates for ray/side intersection
		int   side;   // index of the side within the element
	};

    struct el_isect_info_t : isect_info_t {
		isect_info_t out;	// out info of intersection with element
	};

    enum VertexType {
    	INSIDE_SHAPE = 0,
    	EDGE_SHAPE = (1 << 0),
    	CORNER_SHAPE = (1 << 1),
    };





#pragma pack(push,1)

   typedef struct _Node_t {
      fvmath::CVec3d position;
      int  info; // Index in global mesh
      int  type; // State for outgoing edge in triangle
    } Node_t;

    struct BaseVertex {
    	typedef float info_t;
    	fvmath::CVec3f position;
    	info_t info;
    	//info_t type;

    };

    struct Vertex : public BaseVertex {
    	typedef float value_type;
    	float          value;
    	fvmath::CVec3f color;

    	Vertex()
    	{
    		position = fvmath::CVec3f(LARGE_F,LARGE_F,LARGE_F);
    		info = 0;
    		value = -LARGE_F;
    		color = fvmath::CVec3f(-LARGE_F,-LARGE_F,-LARGE_F);
    	}

    	Vertex(fvmath::Vec3f& pos_,int info_ = 0,float val = -LARGE_F)
    	{
    		position = pos_;
    		info = info_;
    		value = val;
    		color = fvmath::CVec3f(-LARGE_F,-LARGE_F,-LARGE_F);
    	}

    	Vertex(const Node_t& node)
    	{
        	BaseVertex::position = node.position;
        	BaseVertex::info = node.info;
        	BaseVertex::info = node.info;

    	}
    };
    struct MeshVertex {
    	//MeshVertex() : position(), info(0), value(-LARGE_F) {}
    	fvmath::CVec3d position;
    	double value;
    	int info;
    	int type;
    };

    typedef BaseVertex Origin;

    inline bool operator<(const Vertex& a, const Vertex& b) {
    	if (!fvmath::is_near( a.position.x, b.position.x, 0.001f)) return a.position.x < b.position.x;
    	if (!fvmath::is_near( a.position.y, b.position.y, 0.001f)) return a.position.y < b.position.y;
    	if (!fvmath::is_near( a.position.z, b.position.z, 0.001f)) return a.position.z < b.position.z;
    	return false;
    }

#pragma pack(pop)

struct BaseVertexInserterID {
	bool operator()(const BaseVertex& lhs,const BaseVertex& rhs) {
		return lhs.info < rhs.info;
	}
};

enum {
	HORIZONTAL=0,
	VERTICAL,
};

/* Rectangular color bar parameters */
struct colorbar_config_t {
  char type;		/* Type of color bar: g,G - gradient; f,F - flat */
  char orient;   	/* Orientation: h,H - horizontal; v,V - vertical */
  float x,y;            /* Lower left corner of rectangle bar */
  float w,h;            /* CoorBar's dimensions */
  double min,max;       /* Min/Max scalar values */
  double start,end;	/* Start/End values - should be Start >= Min AND End <= Max */
  float out_col[4];	/* Outside color on bar when Start > Min OR End < Max
						 * IF out_col[0] < 0 OR out_col[0] > 1, it means that color ashould be specified
						 * corespodly to first and last colors in given table of colors */
  float* colors;        /* Handle to array of colors */
  int size;		/* Number of colors in given table */
  int stride;		/* A stride between correspondig colors */
  float text_col[4];	/* Text color in labels */
  float theight;	/* TA height of the text field */
  void* font;		/* Font type */
  float edge_col[4];	/* Color of outlines */
  float line_width;	/* Line width */
  float out_width_coef;	/* Percent of color bar width for pre(su)fix rectangle fields */
  int num_ranges;	/* The number of labels >= 2 */
  //const char* tformat;  /*text format for displayed values */
};


struct Vertex2D {
  float x, y;
};

struct CVertex2D : public Vertex2D {
  float r,g,b;
  CVertex2D() { memset(this,0,sizeof(Vertex2D)); }
  CVertex2D(float xx,float yy,float rr,float gg,float bb)
  { x = xx; y = yy; r = rr; g = gg; b = bb; }
};

struct label_t {
  float xts,yts;
  char label[12];
  label_t() { memset(this,0,sizeof(label_t)); }
  label_t(float _xs,float _ys,const char* _label)
  : xts(_xs),yts(_ys)
  {
	//printf("Called ctr\n");
	strcpy(label,_label);
  }
  label_t(const label_t& rhs)
  {
	//printf("Called copy ctr\n");
	memcpy(this,&rhs,sizeof(label_t));
  }

 label_t(const label_t&& rhs)
  {
	//printf("Called move ctr\n");
	memcpy(this,&rhs,sizeof(label_t));
  }

};


//	typedef enum MenuEvent
//	{
//			IDM_SEPARATOR = -1,
//		#ifndef _USE_FV_LIB
//			IDM_OPEN_MESH = 0,
//			IDM_OPEN_FIELD,
//			IDM_REFRESH,
//		#else
//			IDM_REFRESH = 0,
//		#endif
//			IDM_RELOAD,
//			IDM_RESET,
//			IDM_QUIT,
//
//			IDM_VPERSP,
//			IDM_VORTHO,
//			IDM_VTOP,
//			IDM_VBOTTOM,
//			IDM_VFRONT,
//			IDM_VBACK,
//			IDM_VLEFT,
//			IDM_VRIGHT,
//			IDM_VDEFAULT,
//			IDM_VFULL,
//			IDM_VBBOX,
//			IDM_VFAST,
//			IDM_VNEW,
//			IDM_VNEXT,
//			IDM_VPREV,
//			IDM_VDUMP_CURR,
//			IDM_VDUMP_ALL,
//
//			IDM_CAXES,
//			IDM_CGRID,
//			IDM_CBKG_COLOR,
//			IDM_CLIGHT,
//			IDM_CLEGEND,
//			IDM_CRESET,
//			IDM_CSAVE,
//			IDM_CLEG_EDIT,
//			IDM_CMOD_APR,
//
//			IDM_RSOL_SET,
//			IDM_RDRAW_WIRE,
//			IDM_RDRAW_FILL,
//			IDM_RDRAW_CONT,
//			IDM_RDRAW_FLOODED,
//			IDM_RDRAW_CUT,
//			IDM_RCUT_SETS,
//			IDM_RSCREEN_SAVE,
//
//			IDM_HELP,
//		} MenuEvent_t;

} // end namespace

#endif /* _ENUMS_TYPE_H_
*/
