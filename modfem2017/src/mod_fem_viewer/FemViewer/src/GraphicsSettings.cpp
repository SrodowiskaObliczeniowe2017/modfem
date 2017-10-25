#include "fv_assert.h"
#include "fv_config.h"
#include "GraphicsSettings.h"

#include <string>
#include <algorithm>
#include <cstdio>

#ifdef WIN32
// Disable the "This function or variable may be unsafe. Consider using sprintf_s instead." etc.
#pragma warning( disable: 4996 )
#endif

namespace FemViewer {


// By default, distant light
GraphicsSettings::GraphicsSettings()
{
	//mfp_log_debug("graphicsSettings ctr");
	Defaults();
}

GraphicsSettings::~GraphicsSettings()
{
	//mfp_log_debug("GraphicsSettings dtr");
}

void GraphicsSettings::Defaults()
{
	eSelect = eSelectionCategory::All;
	eRenderType = RASTERIZATION_GL;
	BkgColor = ColorRGB(bkg_color);
	NodeIDColor = ColorRGB(node_id_color);
	ElemIDColor = ColorRGB(elem_id_color);
	iBackgroundMode = 0;
	bIsGridOn = false;
	bIsAxesOn = true;
	bIsLegendOn = true;
	bEdgeOn = true;
	bShadingOn = true;
	bIsovalueLineOn = true;
	bLineCoolored = true;
	bLineCoolored = false;
	bIsSmothOn = false;
	bShowNumVertices = false;
	bShowNumElems = false;
	bDisplaySlices = false;
	fAmbientLight = 0.3f;
	iLightModel = 1;
	iLightPositionLocal = 0;
	vLightPos = CVec3f(fixed_light_pos);
	eSmpMode = simpleMode_none;
	bIsOrthoOn = false;
	bIsBVHGridDraw = false;
	bIsBVHOctreeDraw = false;
	bIsRayTracing = false;
	DirectionalLight = Light();
	bSliceModel = true;
	for(int i=1;i< sizeof(bShowCutPlane)/sizeof(bShowCutPlane[0]); ++i) bShowCutPlane[i] = false;
}

// Return 0 if no error (see pError).  This does not override settings
// that arent listed in the file.
// return 1 if error, return 2 if no such file
int GraphicsSettings::LoadFromFile(const std::string& sfileName, std::string& sError)
{
 // FILE* file = fopen(sfileName.c_str(),"r");
 // if(!file) {
 //   sError = "Cannot open file : ";
	//sError += sfileName;
 //   return 2;
 // }

 // float lLightPositionX = vLightPos._x();
 // float lLightPositionY = vLightPos._y();
 // float lLightPositionZ = vLightPos._z();

 // char lLine[200];
 // while(fgets(lLine,200,file)) {

 //   if(lLine[0] == '\n' || lLine[0] == '#') {
 //     //Skip blank/comment lines
 //     continue;
 //   }
 //   std::string lLineStr(lLine);
 //   int lPosEqual = lLineStr.find("=");
 //   if(lPosEqual == std::string::npos) {
 //     sError = "Syntax error in configuration file: ";
 //     sError += sfileName;
 //     fclose(file);
 //     return 1;
 //   }

 //   std::string lVariable = lLineStr.substr(0,lPosEqual);
 //   std::string lParam    = lLineStr.substr(lPosEqual+1);

 //   if(lVariable == "background_mode") {
 //     sscanf(lParam.c_str(),"%d",&iBackgroundMode);
 //   }
 //   else if(lVariable == "background_r") {
 //     int lBackgroundR;
 //     sscanf(lParam.c_str(),"%d",&lBackgroundR);
 //     lBackgroundR = std::max(lBackgroundR, 0);
 //     lBackgroundR = std::min(lBackgroundR, 255);
 //     aBackgroundR = static_cast<float>(lBackgroundR)/255.0f;
 //   }
 //   else if(lVariable == "background_g") {
 //     int lBackgroundG;
 //     sscanf(lParam.c_str(),"%d",&lBackgroundG);
 //     lBackgroundG = std::max(lBackgroundG, 0);
 //     lBackgroundG = std::min(lBackgroundG, 255);
 //     aBackgroundG = static_cast<float>(lBackgroundG)/255.0f;
 //   }
 //   else if(lVariable == "background_b") {
 //     int lBackgroundB;
 //     sscanf(lParam.c_str(),"%d",&lBackgroundB);
 //     lBackgroundB = std::max(lBackgroundB, 0);
 //     lBackgroundB = std::min(lBackgroundB, 255);
 //     aBackgroundB = static_cast<float>(lBackgroundB)/255.0f;
 //   }
 //   else if(lVariable == "enable_grid") {
 //     int lFlagGrid;
 //     sscanf(lParam.c_str(),"%d",&lFlagGrid);
 //     aFlagGrid = (lFlagGrid == 1);
 //   }
 //   else if(lVariable == "enable_axes") {
 //     int lFlagAxes;
 //     sscanf(lParam.c_str(),"%d",&lFlagAxes);
 //     aFlagAxes = (lFlagAxes == 1);
 //   }
 //   else if(lVariable == "graphic_simplication") {
 //     int lSimplificationMode = -1;
 //     sscanf(lParam.c_str(),"%d",&lSimplificationMode);
 //     if (lSimplificationMode == 1) {
 //       aSimplicationMode = simplificationMode_bounding_box;
 //     }
 //     else if (lSimplificationMode == 2) {
 //       aSimplicationMode = simplificationMode_fast;
 //     }
 //     else if (lSimplificationMode == 3) {
 //       aSimplicationMode = simplificationMode_fast;
 //       aSimplicationMode_param = 1;
 //     }
 //     else {
 //       aSimplicationMode = simplificationMode_none;
 //     }
 //   }
 //   else if(lVariable == "light_model") {
 //     sscanf(lParam.c_str(),"%d",&aLightModel);
 //   }
 //   else if(lVariable == "light_position_x") {
 //     sscanf(lParam.c_str(),"%f",&lLightPositionX);
 //   }
 //   else if(lVariable == "light_position_y") {
 //     sscanf(lParam.c_str(),"%f",&lLightPositionY);
 //   }
 //   else if(lVariable == "light_position_z") {
 //     sscanf(lParam.c_str(),"%f",&lLightPositionZ);
 //   }
 //   else if(lVariable == "light_position_local") {
 //     sscanf(lParam.c_str(),"%d",&aLightPositionLocal);
 //   }
 //   else if(lVariable == "light_ambient") {
 //     sscanf(lParam.c_str(),"%f",&aLightAmbient);
 //   }
 // }
 // fclose(file);

 // aLightPosition = Vector3D(lLightPositionX, lLightPositionY, lLightPositionZ);

  return 0;
}

void GraphicsSettings::SaveToFile(const std::string& pFilename, std::string& pError)
{
 /* FILE* file = fopen(pFilename.c_str(),"w");
  if(!file) {
    pError = "Cannot write to file : ";
    pError += pFilename;
  }
  else {
    fprintf(file,"# Background color mode: integer value [0,4]\n");
    fprintf(file,"#   0 = user settings (see variable background_r,g,b), default=black\n");
    fprintf(file,"#   1 = white\n");
    fprintf(file,"#   2 = gray\n");
    fprintf(file,"#   3 = gray gradient\n");
    fprintf(file,"#   4 = black gradient\n");
    fprintf(file,"background_mode=%d\n\n",aBackgroundMode);
    fprintf(file,"# Background color (if background_mode is 0): R G B value [0,255]\n");
    fprintf(file,"background_r=%d\n",static_cast<int>(aBackgroundR*255.0f));
    fprintf(file,"background_g=%d\n",static_cast<int>(aBackgroundG*255.0f));
    fprintf(file,"background_b=%d\n",static_cast<int>(aBackgroundB*255.0f));
    fprintf(file,"# Show grid [0,1]\n");
    fprintf(file,"enable_grid=%d\n\n",aFlagGrid);
    fprintf(file,"# Show axes [0,1]\n");
    fprintf(file,"enable_axes=%d\n\n",aFlagAxes);
    fprintf(file,"# Graphic simplification [0,1,2]\n");
    fprintf(file,"#   0 = no simplification\n");
    fprintf(file,"#   1 = using bounding Boxes for simplification (-bbox)\n");
    fprintf(file,"#   2 = using a small subset of the primitives for simplification (-fast)\n");
    fprintf(file,"#  If 1 or 2 is selected, the 3D model is hidden while using the mouse to rotate, translate, zoom\n");
    fprintf(file,"#    thus, enabling a faster user interaction.  Use this for large models\n");
    switch (aSimplicationMode) {
    case simplificationMode_none:
      fprintf(file,"graphic_simplication=0\n\n");
      break;
    case simplificationMode_bounding_box:
      fprintf(file,"graphic_simplication=1\n\n");
      break;
    default:
      GLV_ASSERT(aSimplicationMode == simplificationMode_fast);
      if(aSimplicationMode_param == 1) {
        fprintf(file,"graphic_simplication=3\n\n");
      }
      else {
        fprintf(file,"graphic_simplication=2\n\n");
      }
      break;
    }
    fprintf(file,"# Light model [0,2]\n");
    fprintf(file,"#  0 = Flag shading, no lighting\n");
    fprintf(file,"#  1 = Light is camera  (default)\n");
    fprintf(file,"#  2 = Fixed , see light_x,y,z \n");
    fprintf(file,"light_model=%d\n",aLightModel);
    fprintf(file,"# Light position (if light_model is 2)\n");
    fprintf(file,"#  x,y,z = floating point position/vector  local : 0=infinite light 1=local light\n");
    fprintf(file,"light_position_x=%f\n",aLightPosition.x());
    fprintf(file,"light_position_y=%f\n",aLightPosition.y());
    fprintf(file,"light_position_z=%f\n",aLightPosition.z());
    fprintf(file,"light_position_local=%d\n\n",aLightPositionLocal);
    fprintf(file,"# Light ambient constant. floating point value [0,1]\n");
    fprintf(file,"#   if 0: unlit part of the model are totally black\n");
    fprintf(file,"light_ambient=%f\n\n",aLightAmbient);

    fclose(file);
  }*/
	assert(!"Not implemented");
}

} // end namespace FemViewer

