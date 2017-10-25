#ifndef _VIEW_MANAGER_H_
#define _VIEW_MANAGER_H_
#include"types.h"
#include"GraphicData.h"
#include"GraphicsSettings.h"
#include"Vec3D.h"
#include"BBox3D.h"
#include"View.h"
#include"Object.h"
#include"Legend.h"
#include"RContext.h"
#include"fv_dictstr.h"
#include"fv_timer.h"
#include "fv_config.h"
#include "Enums.h"
#include<vector>
#include<string>

namespace FemViewer
{
	class VtxAccumulator;
	class Tile;
	class GraphicsSettings;
	class BBox3D;

	class ViewManager
	{
	public:

		friend ViewManager& ViewManagerInst(void);
		friend class RContext;

		/* Ctr & Dtr
		**/
		ViewManager(int destWidth = DFLT_IMAGE_WIDTH, int destHeight = DFLT_IMAGE_HEIGHT);
	public:
		~ViewManager();
		/* Initialize member must be called after OpenGL init
		**/
		void Init(const BBox3D& pBBox);
		/* Reset ViewManger to its defaults
		**/
		void Reset();
		/* Clear OpenGL data
		**/
		//void Clear();
		/* Adds new view
		**/	
		void AddView(const View& pView);
		/* Generates status message on screen
		**/ 
		void AddStatusMsg(const std::string& msg);
		/* Clear all messages from window
		**/
		void ClearStatusMsg();
		/* Called every time if need to redraw content
		**/
		void Display(const Tile& pTile);

		void DisplayCallback();

		GraphicData& GetGraphicData();

		GraphicsSettings* GetSettings() { return &UserSettings; }
		const GraphicsSettings* GetSettings() const { return &UserSettings; }

		bool IdleCallback();

		bool KeyboardCallback(unsigned char key);

		bool KeyboardSpecialCallback(int key);

		bool MenuCallback(const int& iEvent);

		bool MouseButtonCallback(int button,
                             int state,
                             unsigned keystate,
                             int x,
                             int y);

		bool MouseMotionCallback(int pX, int pY);

		void ReshapeCallback(int pWidth, int pHeight);

		void InitGL(int width,int height);

		void SetTitle(const std::string* pTitle);

		void SetLegend(const Legend& pLegend);
		Legend& GetLegend() { return legend; }

		enum ViewResetMode { ViewReset_None, ViewReset_BBoxOnly, ViewReset_FullReset };
  
		void SetViewResetMode(ViewResetMode pMode) { eViewResetMode = pMode; }

		bool IsGLSupported() const { return bGLinited; }
		      bool& GLSupport() { return bGLinited; }
		const bool& GLSupport() const { return bGLinited; }

		bool IsVBOSuported() const { return vboSupported; }

		bool IsVBOUsed() const { return vboUsed; }

		//RContext GetCurrentRCoontext();
		const View& GetCurrentView() const { return currentView; }

		int GetWidth() const { return iWndWidth; }
		void SetWidth(int w) {
			iWndWidth = w;
			fAspectRatio = (float) iWndWidth / (float) iWndHeight;
		}
		int GetHeight() const { return iWndHeight; }
		void SetHeight(int h) {
			iWndHeight = h;
			fAspectRatio = (float) iWndWidth / (float) iWndHeight;
		}
		int GetMouseMode() const { return eMoveMode; }

		void DrawAccelStruct(const float [],const float [],const grid_t* grid);
	private:
		//VtxAccumulator* pRenderer;
		BBox3D aBBox;
		bool bGLinited;
		Object oAxes;
		View currentView;
		GraphicData aGraphicData;
		Legend legend;
		Object oGrid;
		Object oLegend;
		int iLastMouseX;
		int iLastMouseY;
		MouseMode eMoveMode;
		bool bNewViewAdded;
		ViewResetMode eViewResetMode;
		std::vector<std::string> vecStatusMsgs;
		const std::string* sTitle;
		GraphicsSettings UserSettings;
		int iViewIndex;
		std::vector<View> vecViews;
		int iWndWidth;
		int iWndHeight;
		float fAspectRatio;
		bool vboSupported;
		bool vboUsed;
		bool bGLEWInited;
		bool bUseExternalModule;
		unsigned int uFrameCount;
		bool bMeasurePerformance;
		fv_timer cTimer;
		double dLastTime;
		GLuint gl_grid[2], gl_axes;
		GLuint grid3D[2];

		void InitGLEW(bool quiet = true);
		void DrawAxes();
		void DrawGrid();
		void DrawLegend(/*const Tile& pTile*/);
		void DrawCubeSides();
		void DrawFPS();
		void DrawMouseMovementIndicators();
		void DrawTextLayer(const Tile& pTile);
		void Render(const Tile& pTile);
		void ResetCamera();
		void SetupLighting();
		void CreateCubeGrid();
		void RenderByCastRays();
	    /* Block use these
		**/
		ViewManager(const ViewManager& rhs);
		ViewManager& operator=(const ViewManager& rhs);
	};

	extern ViewManager& ViewManagerInst(void);

} // end namespace FemViewer

#endif /* _VIEW_MANAGER_H_ */
