#ifndef _RENDERPARAMS_H_
#define _RENDERPARAMS_H_

#include "defs.h"
#include "Enums.h"
#include "Color.h"
#include "Light.h"
#include <cstring>
#include <ostream>
#include <type_traits>

namespace FemViewer {

struct BaseParams {
	/* Matrices */
	float proj[16]; // Projectiom matrix
	float view[16]; // Model-view matrix
	/* Colors */
	float wireframe_col[4]; // Color of edges of element when transparent mode is on
	float border_col[4]; // Color of edges of element when transparent mode is off
	float iso_color[4];		// Color of iso-lines
	float iso_values[MAX_ISO_VALUES][4]; // Table of iso values data
	int   num_breaks;	// Number of break-points (isovalues)
	int   num_tri_faces; // Total number of triangle faces in display model
	int   edges_on;		// Flag indicating whether element boundaries are drawn
	int   isolines_on;  // Flag indicating whether iso-lines are drawn
	float light_intensity[4];	// Position of a light
	float light_ambient[4]; // Color of a light

	//float light_pos[4]; // Light position
	BaseParams() {
		int i;
		memset(this,0x0,sizeof(BaseParams));
		// Wite wireframe color
		for (i=0;i<4;i++) wireframe_col[i] = 1.0f;
		// Black shading color
		border_col[3] = 1.f;
		// Gray isovaslue color
		for (i=0;i<3;i++) iso_color[i] = 0.1f;
		iso_color[i] = 1.0f;
		// Default position on a light
		for (i=0;i<4;i++) light_intensity[i] = 1.0f;
	    // White color of a light
		for (i=0;i<3;i++) light_ambient[i] = 1.0f;
	}

};


static_assert(
    (sizeof(BaseParams) == (sizeof(float) * 184) &&
     std::is_standard_layout<BaseParams>::value),

     "BaseParams does not satisfy contiguous storage requirements");

inline std::ostream& operator << (std::ostream& os, const BaseParams& rhs)
{
	os << "Basic shader parameters:\n";
	os << "\tWireframe color:\t\t{"
	   << rhs.wireframe_col[0] << ", "
	   << rhs.wireframe_col[1] << ", "
	   << rhs.wireframe_col[2] << ", "
	   << rhs.wireframe_col[3] << "}\n";
	os << "\tBorder color:\t\t{"
	   << rhs.border_col[0] << ", "
	   << rhs.border_col[1] << ", "
	   << rhs.border_col[2] << ", "
	   << rhs.border_col[3] << "}\n";
	os << "\tIsovalue color:\t\t{"
	   << rhs.iso_color[0] << ", "
	   << rhs.iso_color[1] << ", "
	   << rhs.iso_color[2] << ", "
	   << rhs.iso_color[3] << "}\n";
	os << "\tNumber of iso-breaks: " << rhs.num_breaks << "\n";
	for (int i=0;i<rhs.num_breaks;i++)
		os << "\t\tvalue: " << rhs.iso_values[i][3] << "\n";

}

struct RenderParams
{
	enum RenderType
	{
		eFull = 0,
		eBoundingBox,
		eFast
	};

	RenderType eRMode;

	Render_t   eRenderType;

	MouseMode  eMouseMode;
	bool       bSmoothNormals;    // Smooth normals, where applicable, and apply them to the primitive
	bool       bFacetFrame; // Render the boundaries of a facet
	bool       bDoubleSided;      // polygons are not culled but also drawn if seen from the back side
	ColorRGB   cBkgColor;
	ColorRGB   cEdgeColor;
	ColorRGB   cVertexIdColor;
	ColorRGB   cElemIdColor;
	int        iPrimitiveOptimizerValue; // Determine the maximum length of a glBegin/glEnd sequence
                                       //  default value of 100 is an acceptable performance/memory
                                       //  tradeoff.
  
	int        iRMode_Fast_Option;
	bool	   bShowNumOfVertices;
	bool       bShowNumOfElements;
	BaseParams sShaderParams;
	bool	   bDrawWireframe;
	bool 	   bDrawCutted;
	bool 	   bColorShading;
	Light 	   cLight;

  RenderParams()
	: eRMode                (eFull)
	, eRenderType			(RASTERIZATION_GL)
	, eMouseMode			(MOUSE_NONE)
	, bSmoothNormals		(false)
	, bFacetFrame			(false)
	, bDoubleSided			(false)
	, cBkgColor				()
	, cEdgeColor				()
	, cVertexIdColor(1.f,0.f,1.f)
	, cElemIdColor(0.f,1.f,1.f)
	, iPrimitiveOptimizerValue(100)
	, iRMode_Fast_Option	  (1)
	, bShowNumOfVertices	  (false)
	, bShowNumOfElements	  (false)
	, sShaderParams()
	, bDrawWireframe(false)
	, bDrawCutted(false)
	, bColorShading(true)
	, cLight()
	{
	  ;
    }

};

} // edn namespace FemViewer

#endif /* _RENDERPARAMS_H_
*/

