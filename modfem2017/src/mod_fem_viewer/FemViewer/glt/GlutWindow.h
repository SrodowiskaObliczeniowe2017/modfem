/*
 * GlutWindow.h
 *
 *  Created on: 07-02-2012
 *      Author: dwg
 */

#ifndef GLUTWINDOW_H_
#define GLUTWINDOW_H_

#include <WindowFV.h>
#include <GL/freeglut.h>

namespace FemViewer {
namespace GlutGUI {

static inline bool Init(int argc,char **argv) {
	if (! (glutGet(GLUT_INIT_STATE)==1)) {
		glutInit(&argc,argv);
	}
	return true;
}

class GlutWindow : public Window
{
public:
	static GlutWindow& GetInstance(void);
	static void Run(void);
	/// Constructor
	
	~GlutWindow(void);

	virtual void DoInit(void);
	virtual void RefreshScreen(void);
	virtual void Show(bool show = true);

	int MainLoop(void);
private:
	GlutWindow(int posX,int posY,int width,int height,const char* title);
	friend class FemViewerApp;
	/// Main window id
	int _wndId;

	void InitGL(void);
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
