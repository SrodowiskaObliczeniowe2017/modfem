/*
 * GlutWindow.h
 *
 *  Created on: 07-02-2012
 *      Author: dwg
 */

#ifndef GLUTWINDOW_H_
#define GLUTWINDOW_H_

#include <WindowFV.h>

namespace FemViewer {
namespace GlutGUI {

class FemvVieweApp;

class GlutWindow: public WindowFV
{
public:
	static void Init(int argc,char** argv);
	static GlutWindow& GetInstance(void);
	static void Run(void);
	/// Constructor
	GlutWindow(int posX,int posY,int width,int height,const char* title);
	virtual ~GlutWindow(void);

	virtual void Refresh(void);
	virtual bool Show(bool show = true);

	int MainLoop();
private:
	friend class FemViewerApp;
	/// Main window id
	int _wndId;

	/// Block use this
	GlutWindow(const GlutWindow&);
	GlutWindow& operator=(const GlutWindow&);

	static void IdleFunc(void);

	static void DisplayFunc(void);

	static void ResizeFunc(int width,int height);

	static void MouseFunc(int button,int state,int x,int y);

	static void MouseMotionFunc(int x,int y);

	static void MenuFunc(int option);

//	static void MenuViewFunc(int option);
//
//	static void MenuConfigFunc(int option);
//
//	static void MenuRenderFunc(int option);

	static void KeyboardFunc(unsigned char key,int x,int y);

};

}
}

#endif /* GLUTWINDOW_H_ */
