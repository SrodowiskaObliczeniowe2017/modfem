#ifndef _GRAPHICS_SETTINGS_H_
#define _GRAPHICS_SETTINGS_H_

#include "Enums.h"
#include "MathHelper.h"
#include "Color.h"
#include "Light.h"
#include "CubeGraphicSettings.h"

namespace FemViewer {

using namespace std;
using namespace fvmath;
class string;

class GraphicsSettings
{
public:
	  enum eMode 
	  { 
		  simpleMode_none,
		  simpleMode_bounding_box,
		  simpleMode_fast
	  };

	  GraphicsSettings();
	  ~GraphicsSettings();

	  void Defaults();
	  int  LoadFromFile(const std::string& sfileName, std::string& sError);
	  void SaveToFile(const std::string& pFilename, std::string& pError);

	  const char * path_to_meshfile;
	  const char * path_to_fieldfile;

	  eSelectionCategory eSelect;
	  Render_t	eRenderType;
	  ColorRGB BkgColor;
	  ColorRGB NodeIDColor;
	  ColorRGB ElemIDColor;
	  int	iBackgroundMode;
	  bool  bIsGridOn;
	  bool  bIsAxesOn;
	  bool  bIsLegendOn;
	  bool 	bEdgeOn;
	  bool	bShadingOn;
	  bool  bIsovalueLineOn;
	  bool  bLineCoolored;
	  bool  bIsSmothOn;
	  bool  bShowNumVertices;
	  bool  bShowNumElems;
	  bool  bDisplaySlices;
	  float fAmbientLight;
	  int   iLightModel;
	  int   iLightPositionLocal;
	  CVec3f vLightPos;
	  eMode eSmpMode;
	  int	iSimplicationMode_param;
	  bool	bIsOrthoOn;
	  bool  bIsBVHGridDraw;
	  bool  bIsBVHOctreeDraw;
	  bool  bIsRayTracing;
	  Light DirectionalLight;
	  bool  bSliceModel;
	  bool  bShowCutPlane[11];

};
} // end namespace FFemViewer


#endif /* _SETTINGS_H_
*/
