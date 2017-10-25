#include "MathHelper.h"
#include "ViewManager.h"
#include "BBox3D.h"
#include "CutPlane.h"
#include "ModelControler.h"
#include "RenderParams.h"
#include "VtxAccumulator.h"
#include "fv_txt_utls.h"
#include "fv_exception.h"
#include "fv_conv_utils.h"
#include "Tile.h"
//#include "RContext.h"
#include "defs.h"
#include "Enums.h"
#include "Log.h"

#ifdef _WIN32
#include"../win/resource.h"
#else
#include"../wx/resource.h"
#include<iostream>
#endif
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#ifndef WIN32
#include <unistd.h>
#endif

#include <GL/freeglut.h>

#define GL_GLEXT_PROTOTYPES
#include <GL/glext.h>

#define FV_GRID_DELTA 0.06250000f

namespace FemViewer {

// Handle to static ViewManager instance 
ViewManager& ViewManagerInst(void) {
	//mfp_debug("ViewManagerInst\n");
	static ViewManager viewMgrImpl;
	return viewMgrImpl;
}

// Constructor; initialize the viewport and load preferences
ViewManager::ViewManager(int destWidth, int destHeight)
: aBBox(fvmath::CVec3f(),fvmath::CVec3f(1.f,1.f,1.f))
, bGLinited(false)
, oAxes("Axes"),
currentView(),
aGraphicData(),
legend(),
oGrid("Grid"),
//oCube		   (0),
oLegend("Legend"),
iLastMouseX(0),
iLastMouseY(0),
eMoveMode(MOUSE_NONE),
bNewViewAdded(false),
eViewResetMode(ViewReset_FullReset),
vecStatusMsgs(),
sTitle(NULL),
UserSettings(),
iViewIndex(-1),
vecViews(),
iWndWidth(destWidth),
iWndHeight(destHeight),
fAspectRatio(1.f),
vboSupported(false),
vboUsed(true),
bGLEWInited(false),
bUseExternalModule(true),
uFrameCount(0),
bMeasurePerformance(true),
cTimer(),
dLastTime(0.0)
{
	//mfp_debug("ViewManager ctr\n");
	if (bMeasurePerformance) {
		cTimer.start();
		dLastTime = cTimer.get_time_sec();
	}
	fAspectRatio = (float)iWndWidth / (float)iWndHeight;
	assert(fAspectRatio > 0.0f);

	//pRenderer = new VtxAccumulator();
	//assert(pRenderer != NULL);
	grid3D[0] = 0;
	grid3D[1] = 0;
}

void ViewManager::Init(const BBox3D& pBBox)
{
	//mfp_debug("Init \n");
	assert(pBBox.isInitialized());
	if (!pBBox.isInitialized()) return;
	aBBox = pBBox;
	//std::cout << pBBox << std::endl;
	ColorRGB col(0.5f);
	Vec3D pt1;
	Vec3D pt2(1.0f, 0.0f, 0.0f);
	Vec3D pt3(0.0f, 1.0f, 0.0f);
	Vec3D pt4(0.0f, 0.0f, 1.0f);

	oAxes.AddArrowColored(pt1, col, pt2, col, 0.05f, 10);
	col.R = 0.0f;
	col.G = 0.5f; 

	oAxes.AddArrowColored(pt1, col, pt3, col, 0.05f, 10);
	col.G = 0.0f;
	col.B = 0.5f;
	oAxes.AddArrowColored(pt1, col, pt4, col, 0.05f, 10);

	col.R = 1.0f; col.G = col.B = 0.0f;
	oAxes.SetDrawColor(col);
	oAxes.AddText(1.02f, 0.02f ,0.02f,Object::fixed13,"\"X\"");
	col.G = 1.0f; col.R = col.B = 0.0f;
	oAxes.SetDrawColor(col);
	oAxes.AddText(0.02f, 1.02f ,0.02f,Object::fixed13,"\"Y\"");
	col.R = col.G = 0.0f; col.B = 1.0f;
	oAxes.SetDrawColor(col);
	oAxes.AddText(0.02f, 0.02f ,1.02f,Object::fixed13,"\"Z\"");

	// grid construction
	double maxDim = aBBox.maxDim();
	double gridSize = maxDim > 1.0 ? 0.2 : 0.2*maxDim;
	//mfp_log_debug("After init %f %f\n",maxDim, gridSize);
	fvmath::CVec3d orig(-4*maxDim,-4*maxDim,0.0);
	fvmath::CVec3d dims(8*maxDim,8*maxDim,0.0);
	fvmath::CVec3d density(gridSize,gridSize,0.0);

	gl_grid[0] = createQuadGrid(orig.v,dims.v,density,&gl_grid[1]);
	//const size_t
	//aBufferName[VBO_FACES] = createBuffer(NULL,offset + data_size,GL_ARRAY_BUFFER,GL_STATIC_DRAW);
	//if (glIsList(gl_grid)) glDeleteLists(gl_grid,1);
	//gl_grid = glGenLists(1);
	//glNewList(gl_grid, GL_COMPILE);
	// get maximal dimension of a model
/*
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glColor3f(0.3f,0.3f,0.5f);
	glBegin(GL_LINES);
	for (float y=-maxDim;y<maxDim;y=y+gridSize)
	{
		glBegin(GL_LINES);
	    glVertex3f( -maxDim, y,0);
	    glVertex3f( maxDim , y,0);
	    glEnd();
	}
	for (float x=-maxDim;x<maxDim;x=x+gridSize)
	{
	    glBegin(GL_LINES);
	    glVertex3f( x,  -maxDim,0);
	    glVertex3f( x,  maxDim,0);
	    glEnd();
	}
	glEnd();
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEndList();
*/

//
//	pt1[0] = 0.5f*(pBBox.Xmin() + pBBox.Xmax()) - maxDim;
//	pt1[1] = 0.5f*(pBBox.Ymin() + pBBox.Ymax()) - maxDim;
//	pt1[2] = 0.0f;
//
//
//	col.R = col.G = col.B = 0.3f;
//
//	pt1[1] = -1.0f;
//	pt2[1] =  1.0f;
//
//	pt3[0] = -1.0f;
//	pt3[1] =  0.0f;
//	pt4[0] =  1.0f;
//	pt4[2] =  0.0f;
//
//	for (int i=-16; i<=16; ++i) {
//
//		const float a = static_cast<float>(i)*FV_GRID_DELTA;
//		pt1[0] = pt2[0] = a;
//		pt3[1] = pt4[1] = a;
//		if ( i == 0) {
//			ColorRGB lYaxCol(0.5f,0.0f,0.0f), lXaxCol(0.0f,0.5f,0.0f);
//			oGrid.AddLineColored(pt1,lXaxCol,pt2,lXaxCol);
//			oGrid.AddLineColored(pt3,lYaxCol,pt4,lYaxCol);
//			continue;
//		}
//		oGrid.AddLineColored(pt1,col,pt2,col);
//		oGrid.AddLineColored(pt3,col,pt4,col);
//	}
//	mfp_log_debug("After init ViewManager\n");
	//FV_CHECK_ERROR_GL();
}

// Destructor
ViewManager::~ViewManager()
{
	//mfp_debug("ViewManager dtr\n");
	if (bGLinited) {
		mfp_debug("deleting buffer for grid");
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		if (glIsBuffer(gl_grid[0]) == GL_TRUE)
			glDeleteBuffers(1,gl_grid);
		bGLinited = false;
	}
}

// Reset 
void ViewManager::Reset()
{
	//pRenderer = new VtxAccumulator();
	//mfp_debug("VIEW_MGR: reset\n");
	//aBBox.mn = fvmath::CVec3f();
	//aBBox.mx = fvmath::CVec3f(1.f,1.f,1.f);
	oAxes.Reset();
	currentView = View();	// set default view
	aGraphicData.Reset();	// clear root object
	legend.Reset();		// set default legend settings
	bNewViewAdded = false;	// no new views
	vecViews.clear();		// clear all stored views
	//const static std::string v0("v0");
	UserSettings.Defaults();
	//this->SetTitle(&v0);
	this->ClearStatusMsg(); // clear all messages
	//glBindBuffer(GL_ARRAY_BUFFER, 0);
	if (bGLinited) {
		if (glIsBuffer(gl_grid[0]) == GL_TRUE)
			glDeleteBuffers(1,gl_grid);
		gl_grid[0] = 0;
		gl_grid[1] = 0;
		Init(aBBox);
	}
}

void ViewManager::AddView(const View& pView)
{
	currentView = pView;
	vecViews.push_back(pView);

	int nr_views = static_cast<int>(vecViews.size());

	iViewIndex    = nr_views - 1;
	bNewViewAdded = true;

	char msg[128];
	sprintf(msg,"Adding view %d", nr_views);
	AddStatusMsg(msg);
}

// Add a status message. (temporary text in the window)
void ViewManager::AddStatusMsg(const std::string& sMsg)
{
	if(vecStatusMsgs.size() == 5) {
	  vecStatusMsgs.erase(vecStatusMsgs.begin());
	}
	vecStatusMsgs.push_back(sMsg);
}

// Clear all status messages
void ViewManager::ClearStatusMsg()
{
	vecStatusMsgs.clear();
}

// prepare the display and renders the scene
void ViewManager::Display(const Tile& pTile)
{
	// Clear the window
	glClearColor(UserSettings.BkgColor.R,
				 UserSettings.BkgColor.G,
                 UserSettings.BkgColor.B,
                 0.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// Renders the background
	pTile.InitViewport();

	//pTile.InitOrtho2DMatrix(); ??

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);

	FV_ASSERT(UserSettings.iBackgroundMode >= 0);
	FV_ASSERT(UserSettings.iBackgroundMode <  5);

	if(UserSettings.iBackgroundMode < 3) {
		// Plain color
	} else if(UserSettings.iBackgroundMode == 3) {
		glBegin(GL_QUADS);
			glColor3f (0.7f, 0.7f, 0.7f);
			glVertex2f(0.0f, 0.0f);
			glColor3f (0.5f, 0.5f, 0.5f);
			glVertex2f(0.0f, 1.0f);
			glColor3f (0.4f, 0.4f, 0.4f);
			glVertex2f(1.0f, 1.0f);
			glColor3f (0.3f, 0.3f, 0.3f);
			glVertex2f(1.0f, 0.0f);
		glEnd();
	} else if(UserSettings.iBackgroundMode == 4) {
		glBegin(GL_QUADS);
			glColor3f (0.0f, 0.0f, 0.0f);
			glVertex2f(0.0f, 0.0f);
			glColor3f (0.0f, 0.0f, 0.0f);
			glVertex2f(0.0f, 1.0f);
			glColor3f (0.3f, 0.3f, 0.3f);
			glVertex2f(1.0f, 1.0f);
			glColor3f (0.2f, 0.2f, 0.2f);
			glVertex2f(1.0f, 0.0f);
		glEnd();
	}

	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);

	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();


	// Render the grid, axes, graphic data
	Render(pTile);

	// Reset openGL attribute state
	glPopAttrib();
}

// Callback: Called when the scene has to be rendered.
void ViewManager::DisplayCallback()
{
	Display(Tile(iWndWidth, iWndHeight, fAspectRatio));
}

// Return the graphic data structure; used to load input data
GraphicData& ViewManager::GetGraphicData()
{
	return aGraphicData;
}


// Callback: Called when the app is idle.
// We check if there is new data available.
// If so, we return true (that will tell the window to redraw)
bool ViewManager::IdleCallback()
{
	//const bool lNewDataParsed = aGraphicData.timerCallback();

  // Change currentView only if there is new data parsed
  // and that it wasn't a new view
  //if (lNewDataParsed && !bNewViewAdded) {
  //  if(eViewResetMode == ViewReset_BBoxOnly) {
  //    currentView = View(aGraphicData.GetGlobalBBox3D(), currentView);
  //  }
  //  else if(eViewResetMode == ViewReset_FullReset) {
  //    currentView = View(aGraphicData.GetGlobalBBox3D());
  //  }
  //}
	//mfp_debug("ViewManger: idle\n");
  bNewViewAdded = false;

  return false;
}

// Callback: Called when the keyboard is used
bool ViewManager::KeyboardCallback(unsigned char pKey)
{
  if(pKey == '1') {
    MenuCallback(IDM_VIEW_DEFAULT);
  }
  else if(pKey == '2') {
    MenuCallback(IDM_VIEW_TOP);
  }
  else if(pKey == '3') {
    MenuCallback(IDM_VIEW_BOTTOM);
  }
  else if(pKey == '4') {
    MenuCallback(IDM_VIEW_FRONT);
  }
  else if(pKey == '5') {
    MenuCallback(IDM_VIEW_BACK);
  }
  else if(pKey == '6') {
    MenuCallback(IDM_VIEW_LEFT);
  }
  else if(pKey == '7') {
    MenuCallback(IDM_VIEW_RIGHT);
  }
  else if(static_cast<int>(pKey) == 17) {
    // control-q
	  MenuCallback(IDM_FILE_EXIT);
  }
  else if(pKey == 'F') {
    MenuCallback(IDM_VIEW_FULL);
  }
  else if(pKey == 'D') {
	  MenuCallback(IDM_VIEW_BOUNDINGBOX);
  }
  else if(pKey == 'S') {
    MenuCallback(IDM_VIEW_FAST);
  }
  else if(pKey == 'V') {
	  MenuCallback(IDM_VIEW_DUMPCURRENTVIEW);
  }
  else if(pKey == '~') {
	  MenuCallback(IDM_VIEW_NEWVIEW);
  }
  else if(pKey == '>') {
	  MenuCallback(IDM_VIEW_NEWVIEW);
  }
  else if(pKey == '<') {
	  MenuCallback(IDM_VIEW_PREVIOUSVIEW);
  }
  else if(pKey == 'G') {
	  MenuCallback(IDM_CONFIGURE_GRID);
  }
  else if(pKey == 'A') {
	  MenuCallback(IDM_CONFIGURE_AXES);
  }
  else if(pKey == 'B') {
	  MenuCallback(IDM_CONFIGURE_BACKGROUNDCOLOR);
  }
  else if(pKey == 'L') {
	  MenuCallback(IDM_CONFIGURE_LIGHT);
  }
  else if(pKey == 'M'){
	  //MenuCallback("view_mesh");
  } else if(pKey == 'P'){
	  MenuCallback(IDM_VIEW_PERSPECTIVE);
  } else if(pKey == 'O'){
	  MenuCallback(IDM_VIEW_ORTOGRAPHIC);
  }
  return true;
}

// Callback : Called when special keys on the keyboard are used
bool ViewManager::KeyboardSpecialCallback(int pKey)
{
	return true;
}

// Callback: menus commands with "event" string
bool ViewManager::MenuCallback(const int& iEvent)
{
	//bool flag = false;
	if (iEvent == IDM_FILE_REFRESH) {
		ModelCtrlInst().Do(ModelCtrl::UPDATE | ModelCtrl::INIT_FIELD);
		ModelCtrlInst().Do();
	}
	else if (iEvent == IDM_FILE_RELOAD) {
		ModelCtrlInst().Do(ModelCtrl::CLEAR);
		ModelCtrlInst().Do();
	}
	else if (iEvent == IDM_FILE_CLEAR) {
		aGraphicData.Reset();
		ModelCtrlInst().Do(ModelCtrl::CLEAR);
		vecViews.clear();
		iViewIndex   = -1;
		currentView = View();
		return true;
	}
	else if (iEvent == IDM_VIEW_PERSPECTIVE) {
		currentView.SwitchProjection(View::ePerspective);
		AddStatusMsg("perspective view enabled");
	}
	else if (iEvent == IDM_VIEW_ORTOGRAPHIC) {
		currentView.SwitchProjection(View::eParallel);
		AddStatusMsg("Ortho view enabled");
	}
	else if (iEvent == IDM_VIEW_TOP) {
		currentView = View(currentView.GetType(),
				           aGraphicData.GetGlobalBBox3D(),
		                   CVec3f(0.0f, 1.0f,  0.0f),
		                   CVec3f(0.0f, 0.0f, -1.0f));
	}
	else if (iEvent == IDM_VIEW_BOTTOM) {
		currentView = View(currentView.GetType(),
		           	   	   aGraphicData.GetGlobalBBox3D(),
		                   CVec3f(0.0f, -1.0f,  0.0f),
		                   CVec3f(0.0f,  0.0f, -1.0f));
	}
	else if (iEvent == IDM_VIEW_FRONT) {
		currentView = View(currentView.GetType(),
				           aGraphicData.GetGlobalBBox3D(),
		                   CVec3f(0.0f, 0.0f, 1.0f),
		                   CVec3f(0.0f, 1.0f, 0.0f));
	}
	else if (iEvent == IDM_VIEW_BACK) {
		currentView = View(currentView.GetType(),
		           	   	   aGraphicData.GetGlobalBBox3D(),
		                   CVec3f(0.0f, 0.0f,-1.0f),
		                   CVec3f(0.0f, 1.0f, 0.0f));

	}
	else if (iEvent == IDM_VIEW_LEFT) {
		currentView = View(currentView.GetType(),
		           	   	   aGraphicData.GetGlobalBBox3D(),
		                   CVec3f(-1.0f, 0.0f, 0.0f),
		                   CVec3f( 0.0f, 1.0f, 0.0f));
	}
	else if (iEvent == IDM_VIEW_RIGHT) {
		currentView = View(currentView.GetType(),
				           aGraphicData.GetGlobalBBox3D(),
		                   CVec3d(1.0f, 0.0f, 0.0f),
		                   CVec3d(0.0f, 1.0f, 0.0f));
	}
	else if (iEvent == IDM_VIEW_DEFAULT) {
		currentView = View(aGraphicData.GetGlobalBBox3D(),currentView.GetType());
	}
	else if (iEvent == IDM_VIEW_FULL) {
		 UserSettings.eSmpMode = GraphicsSettings::simpleMode_none;
		 UserSettings.bLineCoolored = !UserSettings.bLineCoolored;
		 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	else if (iEvent == IDM_VIEW_BOUNDINGBOX) {
		UserSettings.eSmpMode = GraphicsSettings::simpleMode_bounding_box;
	}
	else if (iEvent == IDM_VIEW_FAST) {
		UserSettings.eSmpMode = GraphicsSettings::simpleMode_fast;
		UserSettings.bLineCoolored = !UserSettings.bLineCoolored;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else if (iEvent == IDM_VIEW_NEWVIEW) {
		AddView(currentView);
	}
	else if (iEvent == IDM_VIEW_NEXTVIEW) {
		int nViews = static_cast<int>(vecViews.size());

		if(nViews > 0) {

			iViewIndex = (iViewIndex + 1) % nViews;

			FV_ASSERT(iViewIndex >= 0);
			FV_ASSERT(iViewIndex <  nViews);
			currentView = vecViews[iViewIndex];

			char lMessage[128];
			sprintf(lMessage,"Using view %d/%d", iViewIndex+1, nViews);
			AddStatusMsg(lMessage);
		}
		else {
			AddStatusMsg("No user defined view");
		}
	}
	else if (iEvent == IDM_VIEW_PREVIOUSVIEW) {
		int nViews = static_cast<int>(vecViews.size());

		if (nViews > 0) {

			if (iViewIndex == -1) {
				iViewIndex = 1;
			}

			iViewIndex = (iViewIndex + nViews - 1) % nViews;

			FV_ASSERT(iViewIndex >= 0);
			FV_ASSERT(iViewIndex <  nViews);
			currentView = vecViews[iViewIndex];

			char lMessage[128];
			sprintf(lMessage,"Using view %d/%d", iViewIndex+1, nViews);
			AddStatusMsg(lMessage);
		}
		else {
			AddStatusMsg("No user defined view");
		}
	}
	else if (iEvent == IDM_VIEW_DUMPCURRENTVIEW) {
		currentView.dump(std::cerr);
		currentView.dump(*this);
	}
	else if (iEvent == IDM_VIEW_DUMPALLVIEWS) {
		int nViews = static_cast<int>(vecViews.size());

		if(nViews == 0) {
		  std::cerr << "No user defined views" << std::endl;
		  AddStatusMsg("No user defined views");
		}
		else {

		  std::vector<View>::const_iterator       it    = vecViews.begin();
		  const std::vector<View>::const_iterator it_e = vecViews.end  ();

		  while (it != it_e) {
			it->dump(std::cerr);
			++it;
		  }
		}
	}
	else if (iEvent == IDM_CONFIGURE_AXES) {
		UserSettings.bIsAxesOn = !UserSettings.bIsAxesOn;

		if(UserSettings.bIsAxesOn) {
			AddStatusMsg("Axes display enabled");
		} else {
			AddStatusMsg("Axes display disabled");
		}
	}
	else if (iEvent == IDM_CONFIGURE_GRID) {
		UserSettings.bIsGridOn = !UserSettings.bIsGridOn;
		if (UserSettings.bIsGridOn) {
			AddStatusMsg("Grid display enabled");
		} else {
			AddStatusMsg("Grid display disabled");
		}
	}
	else if (iEvent == IDM_CONFIGURE_CUT_PLANE) {
		UserSettings.bShowCutPlane[0] = !UserSettings.bShowCutPlane[0];
		if (UserSettings.bShowCutPlane[0]) {
			AddStatusMsg("Cut-plane display enabled");
		} else {
			AddStatusMsg("Cut-plane display disabled");
		}
	}
	else if (iEvent == IDM_CONFIGURE_BACKGROUNDCOLOR) {
		UserSettings.iBackgroundMode = (UserSettings.iBackgroundMode + 1) % 5;

		if(UserSettings.iBackgroundMode == 0) {
		  UserSettings.BkgColor.R = 0.0f;
		  UserSettings.BkgColor.G = 0.0f;
		  UserSettings.BkgColor.B = 0.0f;
		}
		else if(UserSettings.iBackgroundMode == 1) {
		  UserSettings.BkgColor.R = 1.0f;
		  UserSettings.BkgColor.G = 1.0f;
		  UserSettings.BkgColor.B = 1.0f;
		}
		else if(UserSettings.iBackgroundMode == 2) {
		  UserSettings.BkgColor.R = 70.0f/255.0f;
		  UserSettings.BkgColor.G = 70.0f/255.0f;
		  UserSettings.BkgColor.B = 70.0f/255.0f;
		}

	}
	else if (iEvent == IDM_CONFIGURE_LIGHT) {
		UserSettings.iLightModel = (UserSettings.iLightModel + 1) % 3;

		if (UserSettings.iLightModel == Light::Flat) {
			mfp_log_debug("Flat shading");
			UserSettings.DirectionalLight.Type() = Light::Flat;
			AddStatusMsg("Lighting model = Flat shading");
		}
		else if (UserSettings.iLightModel == Light::Camera) {
			mfp_log_debug("Set a light at a camera's position");
			// Light is camera
			UserSettings.DirectionalLight.Type() = Light::Camera;
			UserSettings.DirectionalLight.Position() = currentView.getPosition();
			AddStatusMsg("Lighting model = Light on camera");
		}
		else {
			FV_ASSERT(UserSettings.iLightModel == Light::Fixed);
			mfp_log_debug("Set a light at the fixed position");
			// Light is fixed
			UserSettings.DirectionalLight.Type() = Light::Fixed;
			UserSettings.DirectionalLight.Position() = UserSettings.vLightPos;
			AddStatusMsg("Lighting model = Fixed light");
		}
	}
	else if (iEvent == IDM_CONFIGURE_LEGEND) {

	}
	else if (iEvent == IDM_CONFIGURE_RESET) {

	}
	else if (iEvent == IDM_CONFIGURE_SAVECONFIGURATION) {

	}
	else if (iEvent == IDM_CONFIGURE_EDITLEGEND) {

	}
	else if (iEvent == IDM_CONFIGURE_CHANGEAPROXIMATIONMODULE) {
		if(bUseExternalModule) {
			AddStatusMsg("External modules enabled");
		} else {
			AddStatusMsg("External modules disabled");
		}
	}
	else if (iEvent == IDM_RENDER_SETSOLUTIONFORMULA) {

	}
	else if (iEvent == IDM_RENDER_WIREFRAME) {
		UserSettings.bEdgeOn = !UserSettings.bEdgeOn;
		if(UserSettings.bShadingOn) {
			AddStatusMsg("Edges display enabled");
		} else {
			AddStatusMsg("Edges display disabled");
		}

	}
	else if (iEvent == IDM_RENDER_SHADED) {
		UserSettings.bShadingOn = !UserSettings.bShadingOn;
		if(UserSettings.bShadingOn) {
			AddStatusMsg("Field display enabled");
		} else {
			AddStatusMsg("Field display disabled");
		}
	}
	else if (iEvent == IDM_RENDER_CONTOURLINES) {
		UserSettings.bIsovalueLineOn = !UserSettings.bIsovalueLineOn;
		if(UserSettings.bIsovalueLineOn) {
			AddStatusMsg("Isocontours display enabled");
		} else {
			AddStatusMsg("Isocontours display  disabled");
		}
	}
	else if (iEvent == IDM_RENDER_NUM_VERTEX) {
		UserSettings.bShowNumVertices = !UserSettings.bShowNumVertices;
		if(UserSettings.bShowNumVertices) {
			AddStatusMsg("Vertex id display enabled");
		} else {
			AddStatusMsg("Vertex id display disabled");
		}
	}
	else if (iEvent == IDM_RENDER_NUM_ELEM) {
		UserSettings.bShowNumElems = !UserSettings.bShowNumElems;
		if(UserSettings.bShowNumElems) {
			AddStatusMsg("Element id display enabled");
		} else {
			AddStatusMsg("Element id display disabled");
		}
	}
	else if (iEvent == IDM_RENDER_TRIM) {

	}
	else if (iEvent == IDM_RENDER_SETSECTION) {

	}
	else if (iEvent == IDM_RENDER_DUMPSCREEN) {

	}
//  else if(lEvent == "config_defaults") {
//	  UserSettings = GraphicsSettings();
//  }
//  else if(lEvent == "config_save") {
//
//    char* lHomeDir = getenv("HOME");
//    std::string lConfigFile(lHomeDir);
//    lConfigFile += "/.glvrc";
//
//    std::string lError;
//    UserSettings.saveToFile(lConfigFile, lError);
//
//    if(!lError.empty()) {
//      AddStatusMsg(lError);
//      std::cerr << lError << std::endl;
//    }
//  }
//  else if(lEvent == "view_zoom_100") {
//    currentView = View(aGraphicData.GetGlobalBBox3D(),
//                        Vec3D(0.0f, 0.0f, -1.0f),
//                        Vec3D(0.0f, 1.0f, 0.0f),
//						currentView.bOrthoOn);
//  }
//  else if(lEvent == "view_zoom_150") {
//    currentView = View(aGraphicData.GetGlobalBBox3D() * 3.0f,
//                        Vec3D(0.0f, 0.0f, -1.0f),
//                        Vec3D(0.0f, 1.0f, 0.0f),
//						currentView.bOrthoOn);
//  }
//  else if(lEvent == "view_zoom_200") {
//    currentView = View(aGraphicData.GetGlobalBBox3D() * 4.0f,
//                        Vec3D(0.0f, 0.0f, -1.0f),
//                        Vec3D(0.0f, 1.0f, 0.0f),
//						currentView.bOrthoOn);
//  }
//  else if(lEvent == "view_zoom_300") {
//    currentView = View(aGraphicData.GetGlobalBBox3D() * 5.0f,
//                        Vec3D(0.0f, 0.0f, -1.0f),
//                        Vec3D(0.0f, 1.0f, 0.0f),
//						currentView.bOrthoOn);
//  }
//  else if(lEvent == "view_back") {
//    currentView = View(aGraphicData.GetGlobalBBox3D(),
//                        Vec3D(0.0f, 0.0f, 1.0f),
//                        Vec3D(0.0f, 1.0f, 0.0f),
//						currentView.bOrthoOn);
//  }
//  else if(lEvent == "view_fullscreen") {
//    glutFullScreen();
//  }


  else {
    std::cerr << "Unhandled MENU " << iEvent << std::endl;
    FV_ASSERT(false);
  }

  return true;
}

// Callback: Mouse button.  Since the mouse's button 2 is linked
// to the contextual menus, we won't get callbacks for pButton=2
// pKeyState :   0x1 shift     0x2 control     0x4  alt
// pButton : 0 -lbtn; 1 - rbtn; 2 - mbtn
bool ViewManager::MouseButtonCallback(int pButton, int pState, unsigned pKeyState, int pX, int pY)
{
	ClearStatusMsg();

	pY = iWndHeight - pY;

	// Disable all modes if the mouse button is released
	eMoveMode = MOUSE_NONE;

	if(!pState) {

		// Mouse button is pressed:
		if(pButton == 0 && (pKeyState & 0x1)) {
			eMoveMode = MOUSE_ZOOM;

		}
#ifdef WIN32
		else if(pButton == 2) {
			eMoveMode = MOUSE_TRANSLATE;
		}
#endif
		else if(pButton == 1 && (pKeyState & 0x2)) {
			
			float value=0;
			
			// Read the Z-Buffer value
			glReadPixels(pX,pY,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&value);			

			// Init the camera matrices, for the current window size (no tiling here)
			GLint viewport[4];
			glGetIntegerv(GL_VIEWPORT,viewport);
			float aspect = viewport[0] / (float)viewport[1];
			currentView.initCamera(Tile(viewport[2],viewport[3],aspect));

			//float clickDistance = value*(currentView.getFarClipDistance()-currentView.getNearClipDistance())+currentView.getNearClipDistance();	
			
			// Obtain the computed camera matrices

			GLdouble modelview[16],projection[16];
			GLdouble lPickX,lPickY,lPickZ;
			
			glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
			glGetDoublev(GL_PROJECTION_MATRIX,projection);

			// Compute the reverse computation to obtain the _approximate_ XYZ from
			//  a X,Y and a Z-value
			gluUnProject(pX,pY,value,modelview,projection,viewport,&lPickX,&lPickY,&lPickZ);			

			// Only recenter the camera if the user clicked on something rendered
			//  value==1 is the "background"
			if(value < 1) {
				currentView.setCenter(CVec3d(lPickX,lPickY,lPickZ));
			}
			
		}
		else if(pButton == 1 || (pButton == 0 && (pKeyState & 0x2))) {
		  eMoveMode = MOUSE_TRANSLATE;
		}
		else if (pButton == 0){
		  eMoveMode = MOUSE_ROTATE;
		}
		else if(pButton == 4) {
		  currentView.scaleFOV(0.9f);
		}
		else if(pButton == 3) {
		  currentView.scaleFOV(1.0f/0.9f);
		}
	}

	// If the user used the mouse wheel; we want the state to
	// stay set at mouseMove_zoom
	// This way; in graphic simplification mode; the displays stays simplified
	// between successive mouse-wheel usages.
	if(pButton == 3 || pButton == 4) {
		eMoveMode = MOUSE_ZOOM;
	}

	iLastMouseX = pX;
	iLastMouseY = pY;

#ifndef WIN32
	usleep(1);
#endif

	return true;
}

bool ViewManager::MouseMotionCallback(int pX, int pY)
{
	pY = iWndHeight - pY;

	if(eMoveMode == MOUSE_ROTATE) {

    currentView.rotate(-0.01f* static_cast<float>(pX-iLastMouseX), 0.01f* static_cast<float>(pY-iLastMouseY));

  }
  else if(eMoveMode == MOUSE_TRANSLATE) {

    currentView.translate(static_cast<float>(pX-iLastMouseX)/static_cast<float>(iWndWidth),
                           static_cast<float>(pY-iLastMouseY)/static_cast<float>(iWndHeight));

  }
  else if(eMoveMode == MOUSE_ZOOM) {

    currentView.scaleFOV(1.0f+(0.01f* static_cast<float>(pY-iLastMouseY)));

  }
  // Mousewheel processing is done is the mouseButtonCallback

  iLastMouseX = pX;
  iLastMouseY = pY;

  return true;
}

// Callback: Called when the window was resized
void ViewManager::ReshapeCallback(int width, int height)
{
	//static int firstCall = 0;
	//mfp_debug("rshapeCallback\n");
	iWndWidth  = width;
	iWndHeight = height;
	fAspectRatio = (float)iWndWidth / (float)iWndHeight;
	//currentView.initCamera(Tile(width,height));
	//RenderManagerInst<GLCore>().Setup();
	//if (!firstCall) {
	legend.Create();
	//	++firstCall;
	//}
}

// Initialize OpenGL state machine
void ViewManager::InitGL(int width,int height)
{
//	// get OpenGL info
//	//mfp_debug("ViewMnager: >>>>>>>>>>> InitGL\n");
//	GLenum err = glewInit();
//		if (GLEW_OK != err)
//		{
//		  /* Problem: glewInit failed, something is seriously wrong. */
//		  fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
//
//		}
//    GLCore& glInfo = GLCore::instance();
//    glInfo.getInfo();
//    //glInfo.printSelf();
//
//#ifdef _WIN32
//    // check VBO is supported by your video card
//    if(glInfo.isExtensionSupported("GL_ARB_vertex_buffer_object"))
//    {
//        // get pointers to GL functions
//        glGenBuffersARB = (PFNGLGENBUFFERSARBPROC)wglGetProcAddress("glGenBuffersARB");
//        glBindBufferARB = (PFNGLBINDBUFFERARBPROC)wglGetProcAddress("glBindBufferARB");
//        glBufferDataARB = (PFNGLBUFFERDATAARBPROC)wglGetProcAddress("glBufferDataARB");
//        glBufferSubDataARB = (PFNGLBUFFERSUBDATAARBPROC)wglGetProcAddress("glBufferSubDataARB");
//        glDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC)wglGetProcAddress("glDeleteBuffersARB");
//        glGetBufferParameterivARB = (PFNGLGETBUFFERPARAMETERIVARBPROC)wglGetProcAddress("glGetBufferParameterivARB");
//        glMapBufferARB = (PFNGLMAPBUFFERARBPROC)wglGetProcAddress("glMapBufferARB");
//        glUnmapBufferARB = (PFNGLUNMAPBUFFERARBPROC)wglGetProcAddress("glUnmapBufferARB");
//
//        // check once again VBO extension
//        if(glGenBuffersARB && glBindBufferARB && glBufferDataARB && glBufferSubDataARB &&
//           glMapBufferARB && glUnmapBufferARB && glDeleteBuffersARB && glGetBufferParameterivARB)
//        {
//            vboSupported = vboUsed = true;
//            std::cout << "Video card supports GL_ARB_vertex_buffer_object." << std::endl;
//        }
//        else
//        {
//            vboSupported = vboUsed = false;
//            std::cout << "Video card does NOT support GL_ARB_vertex_buffer_object." << std::endl;
//        }
//    }
//
//#else // for linux, do not need to get function pointers, it is up-to-date
//    if(glInfo.isExtensionSupported("GL_ARB_vertex_buffer_object"))
//    {
//        vboSupported = vboUsed = true;
//        std::cout << "Video card supports GL_ARB_vertex_buffer_object." << std::endl;
//    }
//    else
//    {
//        vboSupported = vboUsed = false;
//        std::cout << "Video card does NOT support GL_ARB_vertex_buffer_object." << std::endl;
//    }
//#endif
//	this->ReshapeCallback(width,height);
//	glShadeModel(GL_SMOOTH);                        // shading mathod: GL_SMOOTH or GL_FLAT
//
//    // enable /disable features
//    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
//    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
//    glEnable(GL_DEPTH_TEST);
//    //glEnable(GL_BLEND);
//    //glDepthFunc (GL_LESS);
//    glDepthFunc (GL_LEQUAL);
//    //glDepthFunc (GL_EQUAL);
//    //glDepthFunc (GL_NOTEQUAL);
//	glClearColor(UserSettings.BkgColor.R,
//				 UserSettings.BkgColor.G,
//                 UserSettings.BkgColor.B,
//                 0);
//    glClearDepth(1.0);
//    glEnable(GL_CULL_FACE);
//    glFrontFace(GL_CCW);
//	glEnable(GL_LINE_SMOOTH);
//	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable(GL_BLEND);
//
//	glLineWidth(1.5125f);
//	//glDepthFunc(GL_ALWAYS);
//	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//	SetupLighting();
//	FV_CHECK_ERROR_GL();
//	//mfp_debug("After initGL\n");
}

void ViewManager::SetTitle(const std::string* strTile)
{
	sTitle = strTile;
}

void ViewManager::SetLegend(const Legend& pLegend)
{
	legend = pLegend;
}

/*RContext ViewManager::GetCurrentRCoontext()
{
	return RContext(*this);
}*/


//void ViewManager::SetViewResetMode(ViewResetMode pMode)
//{
//	eViewResetMode = pMode;
//}


//void ViewManager::DrawMesh(RenderParams& pParams)
//{
	//if(UserSettings.bIsMeshOn && oMesh != NULL) {
	//	glColor3f(1.0f , 1.0f, 1.0f);

	//	pParams.bSmoothNormals       = aGraphicData.bFlagSmoothing;
	//	pParams.iPrimitiveOptimizerValue = aGraphicData.iOptimizerValue;

	//	oMesh->Render(pParams);
	//}
//}


//void ViewManager::DrawField(RenderParams& pParams)
//{
	//if(UserSettings.bIsFieldOn && oField != NULL) {
	//	glColor3f(1.0f , 1.0f, 1.0f);

	//	pParams.bSmoothNormals       = aGraphicData.bFlagSmoothing;
	//	pParams.iPrimitiveOptimizerValue = aGraphicData.iOptimizerValue;
	//
	//	oField->Render(pParams);
	//}


//}

void ViewManager::InitGLEW(bool quiet)
{
//	if (!bGLEWInited) {
//		GLenum err = glewInit();
//		if (GLEW_OK != err) {
//			/* Problem: glewInit failed, something is seriously wrong. */
//			fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
//			exit(-1);
//		}
//		bGLEWInited = true;
//	}
//
//	if (!quiet) {
//		glInfo glInfo;
//		glInfo.getInfo();
//		glInfo.printSelf();
//		fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
//	}
}

void ViewManager::DrawAxes()
{
	//std::cout << "viewmgr:drawaxes\n";
	const BBox3D bbox = aGraphicData.GetGlobalBBox3D();
	const CVec3f& center = bbox.getCenter();
	float length = bbox.getCircumscribedSphereRadius();

	if(length == 0.0f) {
		length = 1.2f;
	} else {
		length *= 1.2f;
	}

	glPushMatrix();

	glTranslatef(center.x, center.y, center.z);
	glScalef(length, length, length);

	RenderParams rparams;
	oAxes.Render(rparams);  

	glPopMatrix();
}

void ViewManager::DrawGrid()
{
	// Disable depth information for grid rendering; they will always be "in the background
	//  and then never intrude over the object data
	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);

	const BBox3D bbox = aGraphicData.GetGlobalBBox3D();
	const CVec3f& center = bbox.getCenter();
	float length = bbox.getCircumscribedSphereRadius();

	glPushMatrix();

	glTranslatef(center.x, center.y - length, center.z);

	if(length == 0.0f) {
		length = 4.0f;
	} else {
		length *= 4.0f;
	}
	
	glScalef(length, length, length);

	RenderParams rparams;
	//rparams.bFacetFrame = true;
	oGrid.Render(rparams);		

	glPopMatrix();

	// Reactivate Z buffer
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
}

void ViewManager::DrawLegend(/*const Tile& pTile*/)
{
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	// to powinno byc zainicjowane juz wczseniej albo
	// w legendzie zeby bylo juz wiliczone!!!!!!!!!!
	if (this->legend.IsValid()) {
		//glPushAttrib(GL_ALL_ATTRIB_BITS);
		///glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		//glDisable(GL_CULL_FACE);
		// Notify Rendermanager

		this->legend.Draw();
		//glPopAttrib();

	}

	//glMatrixMode(GL_PROJECTION);
	//glPopMatrix ();
	//glMatrixMode(GL_MODELVIEW);
	//glPopMatrix ();

}

void ViewManager::DrawFPS()
{
	static std::string fps = "FPS: 0.0";
	// update fps every second
	++uFrameCount;
	double elapsedTime = cTimer.get_time_sec();
	if(elapsedTime >= 1.0)
	{
		std::stringstream ss;
	    ss << std::fixed << std::setprecision(1);
	    ss << (uFrameCount / elapsedTime) << " FPS" << std::ends; // update fps string
	    ss << std::resetiosflags(std::ios_base::fixed | std::ios_base::floatfield);
	    fps = ss.str();
	    uFrameCount = 0;                      // reset counter
	    cTimer.start();                  // restart timer
	 }
	// backup current model-view matrix
	//glMatrixMode(GL_MODELVIEW);
	//glPushMatrix();                     // save current modelview matrix
	//glLoadIdentity();                   // reset modelview matrix

	// set to 2D orthogonal projection
	//glMatrixMode(GL_PROJECTION);        // switch to projection matrix
	//glPushMatrix();                     // save current projection matrix
	//glLoadIdentity();                   // reset projection matrix
	//gluOrtho2D(0, iWndWidth, 0, iWndHeight); // set to orthogonal projection

	//glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
	GLfloat color[4] = {1, 1, 0, 1};
	int textWidth = (int)fps.size() * 8;
	//drawText(fps,iWndWidth - textWidth - 2,iWndHeight - 15,iWndWidth - 2,iWndHeight - 2,
	//		13, GLUT_BITMAP_8_BY_13, false,true);
	drawString(fps.c_str(), iWndWidth-textWidth, iWndHeight-13, color, GLUT_BITMAP_8_BY_13);

	// restore projection matrix
	//glPopMatrix();                      // restore to previous projection matrix

	// restore modelview matrix
	//glMatrixMode(GL_MODELVIEW);         // switch to modelview matrix
	//glPopMatrix();                      // restore to previous modelview matrix

}

void ViewManager::DrawCubeSides()
{
	// Disable depth information for grid rendering; they will always be "in the background
	//  and then never intrude over the object data
	//glDisable(GL_DEPTH_TEST);
	//glDepthMask(GL_FALSE);

	//glEnable(GL_CULL_FACE);
   // glFrontFace(GL_CW);
   // glCullFace(GL_FRONT);

	const BBox3D bbox = aGraphicData.GetGlobalBBox3D();
	const CVec3f& center = bbox.getCenter();
	float length = bbox.getCircumscribedSphereRadius();

	glPushMatrix();

	glTranslatef(center.x, center.y - length, center.z);

	if(length == 0.0f) {
		length = 1.2f;
	} else {
		length *= 1.2f;
	}
	
	glScalef(length, length, length);

	RenderParams rparams;
	//oCube->Render(rparams);

	glPopMatrix();

	// Reactivate Z buffer
	//glDepthMask(GL_TRUE);
	//glEnable(GL_DEPTH_TEST);
}

void ViewManager::DrawMouseMovementIndicators()
{
	// Viewpoint center is shown using a small orange axe; rendered only if
	// the mouse button is pressed

	currentView.drawCenter();
}

void ViewManager::DrawTextLayer(const Tile& pTile)
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	pTile.InitOrtho2DMatrix();



	// Render overlay text
	glColor3f(0.5f, 0.5f, 0.5f);
	std::string lStringText;

	std::vector<std::string>::const_iterator       lMessageIter    = vecStatusMsgs.begin();
	const std::vector<std::string>::const_iterator lMessageIterEnd = vecStatusMsgs.end  ();

	while (lMessageIter != lMessageIterEnd) {
		lStringText += *lMessageIter + "\n";
		++lMessageIter;
	}

	drawText(lStringText, 0.05f, 0.05f, 0.95f, 0.95f, 12, GLUT_BITMAP_8_BY_13, true, false);

	if(sTitle && !sTitle->empty()) {
		drawText(*sTitle, 0.0f, 0.0f, 1.0f, 1.0f, 24, GLUT_BITMAP_TIMES_ROMAN_24, true, true);
	}

	// Draw legend
	if (legend.IsValid() && UserSettings.bIsLegendOn)  {

		//mfp_debug("viewManger: drawLegend\n");
		DrawLegend(/*pTile*/);
		FV_CHECK_ERROR_GL();
	}

	if (bMeasurePerformance) {
		DrawFPS();
		FV_CHECK_ERROR_GL();
	}

	glMatrixMode(GL_PROJECTION);
	glPopMatrix ();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix ();
}

// Render the grid,axes an graphical data
void ViewManager::Render(const Tile& pTile)
{
	//FV_CHECK_ERROR_GL();
	currentView.initCamera(pTile);

	//currentView.setUp();
	/*{
		VtxAccumulator& pAccum(aGraphicData.GetRootObject().GetCurrentVtxAccumulator());
		pAccum.setUniformMatrix(currentView.getCameraMatrix().matrix.data());
	}*/
	if(UserSettings.bIsGridOn) {
		//DrawGrid();
		//glCallList(gl_grid);
		float col[] = {0.3f, 0.3f ,0.3f};
		drawGrid(gl_grid,col);
	}

	if(UserSettings.bIsBVHGridDraw) {
		float col[] = {1.f, 0.f ,1.f};
		float width = 1.3f;
		drawGrid(grid3D,col,width);
	}

	FV_CHECK_ERROR_GL();
	if(eMoveMode != MOUSE_NONE) {
		DrawMouseMovementIndicators();
	}
	FV_CHECK_ERROR_GL();
	// LIGHT
	if(UserSettings.iLightModel) {
		SetupLighting();
	}
	FV_CHECK_ERROR_GL();
	// Default settings (commands may override this)
	glLightModeli  (GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
	glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
	glEnable       (GL_COLOR_MATERIAL);
	glColor3f      (1.0f,1.0f,1.0f);
	glEnable	   (GL_CULL_FACE);

	 //Object rendering. only if the mouse button is not down OR graphic simplification
	 //is disabled.
	//if(UserSettings.bIsMeshOn || UserSettings.bIsFieldOn) {
	/*if(eMoveMode == MOUSE_NONE || UserSettings.eSmpMode == GraphicsSettings::simpleMode_none) {
		RenderParams lParams;
		lParams.eRenderType = RASTERIZATION_GL;
		lParams.eRMode = RenderParams::eFull;
		//DrawMesh(lParams);
		//DrawField(lParams);
		aGraphicData.Render(lParams);
	}
	else if(eMoveMode != MOUSE_NONE && UserSettings.eSmpMode == GraphicsSettings::simpleMode_fast) {
		RenderParams lParams;
		lParams.eRenderType = UserSettings.eRenderType;
		lParams.eRMode = RenderParams::eFast;
		lParams.iRMode_Fast_Option = UserSettings.iSimplicationMode_param;

		aGraphicData.Render(lParams);
	}
	else if(eMoveMode != MOUSE_NONE && UserSettings.eSmpMode == GraphicsSettings::simpleMode_bounding_box) {
		RenderParams lParams;
		lParams.eRenderType = RASTERIZATION_GL;
		lParams.eRMode = RenderParams::eBoundingBox;

		aGraphicData.Render(lParams);
	}*/
	RenderParams lParams;

	//lParams.bDrawWireframe = false;
	lParams.bDrawCutted = UserSettings.bDisplaySlices;
	lParams.bShowNumOfVertices = UserSettings.bShowNumVertices;
	lParams.bShowNumOfElements = UserSettings.bShowNumElems;
	lParams.eMouseMode = eMoveMode;
	memcpy(lParams.sShaderParams.proj,currentView.getProjectionMatrix().matrix.data(),16*sizeof(float));
	memcpy(lParams.sShaderParams.proj,currentView.getCameraMatrix().matrix.data(),16*sizeof(float));
	RenderManagerInst<GLCore>().Render(&lParams);
	FV_CHECK_ERROR_GL();
	if (UserSettings.bShowCutPlane[0]) {
		std::vector<CutPlane>& vcp = ModelCtrlInst().GetCutPlanes();
		for (size_t i = 0; i< vcp.size(); ++i) vcp[i].Draw();

	}
	FV_CHECK_ERROR_GL();


	glDisable(GL_COLOR_MATERIAL);

	if(UserSettings.iLightModel) {
		glDisable(GL_LIGHTING);
	}
	FV_CHECK_ERROR_GL();
	if(UserSettings.bIsAxesOn) {
		DrawAxes();
	}
	FV_CHECK_ERROR_GL();

	// Drawing 2D objects
	DrawTextLayer(pTile);
	FV_CHECK_ERROR_GL();


}


void ViewManager::SetupLighting()
{
	glEnable(GL_LIGHTING);

	// 70% of the object lighting is diffuse
	// 30% is ambient.
	GLfloat light_ambient[]  = {0.3f, 0.3f, 0.3f, 1.0f};
	GLfloat light_diffuse[]  = {0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat light_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
	GLfloat light_position[] = {1.0f, 1.0f, 1.0f, 0.0f};

	light_ambient[0] *= UserSettings.fAmbientLight;
	light_ambient[1] *= UserSettings.fAmbientLight;
	light_ambient[2] *= UserSettings.fAmbientLight;

	if (UserSettings.iLightModel != Light::Flat) {

		light_position[0] = UserSettings.DirectionalLight.Position().x;
		light_position[1] = UserSettings.DirectionalLight.Position().y;
		light_position[2] = UserSettings.DirectionalLight.Position().z;
		light_position[3] = 1.0f;
	}

	glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable(GL_LIGHT0);
}

void ViewManager::RenderByCastRays()
{
	// Init OCL device

}

void ViewManager::DrawAccelStruct(
		const float minb[],
		const float maxb[],
		const grid_t* grid)
{
	GLuint result= createGrid3D(minb,maxb,grid,grid3D);
	//printf("grid[0] = %u\ngrid[1] = %u\n",grid3D[0],grid3D[1]);
}


//void ViewManager::CreateCubeGrid()
//{
//	Reset(oCube);
//
//	Vec3D pt;
//	std::vector<Vec3D> vertices;
//
//	// vertex 0
//	pt[0] = -0.5f; pt[1] = -0.5f; pt[2] = -0.5f;
//	vertices.push_back(pt);
//	// vertex 1
//	pt[0] = 0.5f; pt[1] = -0.5f; pt[2] = -0.5f;
//	vertices.push_back(pt);
//	// vertex 2
//	pt[0] = 0.5f; pt[1] = 0.5f; pt[2] = -0.5f;
//	vertices.push_back(pt);
//	// vertex 3
//	pt[0] = -0.5f; pt[1] = 0.5f; pt[2] = -0.5f;
//	vertices.push_back(pt);
//	// vertex 4
//	pt[0] = -0.5f; pt[1] = -0.5f; pt[2] = 0.5f;
//	vertices.push_back(pt);
//	// vertex 5
//	pt[0] = 0.5f; pt[1] = -0.5f; pt[2] = 0.5f;
//	vertices.push_back(pt);
//	// vertex 6
//	pt[0] = 0.5f; pt[1] = 0.5f; pt[2] = 0.5f;
//	vertices.push_back(pt);
//	// vertex 7
//	pt[0] = -0.5f; pt[1] = 0.5f; pt[2] = 0.5f;
//	vertices.push_back(pt);
//
//	// add faces (0; 2; 5)  to object
//	// face 0
//	oCube->AddLine(vertices[0],vertices[1]);
//	oCube->AddLine(vertices[1],vertices[2]);
//	oCube->AddLine(vertices[2],vertices[3]);
//	oCube->AddLine(vertices[3],vertices[0]);
//
//	oCube->AddLine(vertices[4],vertices[5]);
//	oCube->AddLine(vertices[5],vertices[6]);
//	oCube->AddLine(vertices[6],vertices[7]);
//	oCube->AddLine(vertices[7],vertices[4]);
//
//	oCube->AddLine(vertices[0],vertices[4]);
//	oCube->AddLine(vertices[1],vertices[5]);
//	oCube->AddLine(vertices[2],vertices[6]);
//	oCube->AddLine(vertices[3],vertices[7]);
//	// face 2
//	//oCube->AddLine(vertices[1],vertices[0],vertices[4],vertices[5]);
//	//// face 5
//	//oCube->AddLine(vertices[4],vertices[0],vertices[3],vertices[7]);
//	
//shade
//}





} // end namespace FemViewer
