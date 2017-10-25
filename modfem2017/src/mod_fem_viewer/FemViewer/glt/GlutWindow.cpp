/*
 * GlutWindow.cpp
 *
 *  Created on: 07-02-2012
 *      Author: dwg
 */

#include <GlutWindow.h>
#include <ViewManager.h>
#include "fv_inc.h"
#include "fv_config.h"
#include "Enums.h"
#include "../utils/fv_dictstr.h"
#include "../utils/fv_exception.h"
#include "../utils/fv_threads.h"
#include<iostream>

extern utd_critical_section crit_section;

namespace FemViewer {
namespace GlutGUI   {


static int close_it = 0;

// have to rewrite because of win/resource.h

static const StrDict MenuEntries[] =
{
	{Window::IDM_SEPARATOR,"--------------"},
	#ifndef _USE_FV_LIB
	{Window::IDM_OPEN_MESH,"Open mesh file"},
	{Window::IDM_OPEN_FIELD,"Open field file"},
	#endif
	{Window::IDM_REFRESH,"Refresh field data"},
	{Window::IDM_RELOAD,"Reload mesh and field data"},
	{Window::IDM_RESET,"Reset model"},
	{Window::IDM_QUIT,"Exit viewer"},
	{Window::IDM_VPERSP,"Perspective view"},
	{Window::IDM_VORTHO,"Parallel view"},
	{Window::IDM_VTOP,"Top"},
	{Window::IDM_VBOTTOM,"Bottom"},
	{Window::IDM_VFRONT,"Front"},
	{Window::IDM_VBACK,"back"},
	{Window::IDM_VLEFT,"Left"},
	{Window::IDM_VRIGHT,"Front"},
	{Window::IDM_VDEFAULT,"Default"},
	{Window::IDM_VFULL,"Normal"},
	{Window::IDM_VFAST,"Fast"},
	{Window::IDM_VBBOX,"Bounding box"},
	{Window::IDM_VNEW,"New"},
	{Window::IDM_VNEXT,"Next"},
	{Window::IDM_VPREV,"Preview"},
	{Window::IDM_VDUMP_CURR,"Dump current view"},
	{Window::IDM_VDUMP_ALL,"Dump all views"},
	{Window::IDM_CAXES,"Show axes"},
	{Window::IDM_CGRID,"Show grid"},
	{Window::IDM_CBKG_COLOR,"Switch bkg color"},
	{Window::IDM_CLIGHT,"Show legend"},
	{Window::IDM_CLEGEND,"Legend"},
	{Window::IDM_CRESET,"Defaults settings"},
	{Window::IDM_CSAVE,"Save settings"},
	{Window::IDM_CLEG_EDIT,"Edit legend"},
	{Window::IDM_CMOD_APR,"Chose module..."},
	{Window::IDM_RSOL_SET,"Solution settings..."},
	{Window::IDM_RDRAW_WIRE,"Wireframe"},
	{Window::IDM_RDRAW_FILL,"Solution"},
	{Window::IDM_RDRAW_CONT,"Contour lines"},
	{Window::IDM_RDRAW_FLOODED,"Colored edges"},
	{Window::IDM_RDRAW_CUT,"Slice"},
	{Window::IDM_RCUT_SETS,"Slice settings"},
	{Window::IDM_RSCREEN_SAVE,"Save screen"},
	{Window::IDM_HELP,"About"},

};

#define addMenu(event) glutAddMenuEntry(MenuEntries[event+1].string,MenuEntries[event+1].key); \
  std::cout << "dodałem: " << (event) << " z pozycja: " << MenuEntries[(event)+1].string << std::endl
#define addSepr		   addMenu(IDM_SEPARATOR)

//void GlutWindow::Init(int argc,char** argv)
//{
//	if (! (glutGet(GLUT_INIT_STATE)==1))
//		glutInit(&argc,argv);
//}

GlutWindow& GlutWindow::GetInstance()
{
	assert(_self != NULL);
	return static_cast<GlutWindow&>(*_self);
}

void GlutWindow::Run()
{
	if (_self != NULL) glutMainLoop();
}

GlutWindow::GlutWindow(int posX,int posY,int width,int height,const char* title)
: _wndId(-1)
{
	assert(this == _self);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(posX,posY);

	_wndId = glutCreateWindow(title);
	if(_wndId != 1) throw fv_exception("Error while creating GLUT Window!");

	glutDisplayFunc(DisplayFunc);
	glutReshapeFunc(ResizeFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MouseMotionFunc);
	glutKeyboardFunc(KeyboardFunc);


	const int menu_view = glutCreateMenu(MenuFunc);
	addMenu(IDM_VPERSP);
	addMenu(IDM_VORTHO);
	addSepr;
	addMenu(IDM_VTOP);
	addMenu(IDM_VBOTTOM);
	addMenu(IDM_VFRONT);
	addMenu(IDM_VBACK);
	addMenu(IDM_VLEFT);
	addMenu(IDM_VRIGHT);
	addSepr;
	addMenu(IDM_VFULL);
	addMenu(IDM_VFAST);
	addMenu(IDM_VBBOX);
	addSepr;
	addMenu(IDM_VNEW);
	addMenu(IDM_VNEXT);
	addMenu(IDM_VPREV);
	addSepr;
	addMenu(IDM_VDUMP_CURR);
	addMenu(IDM_VDUMP_ALL);
	const int menu_config = glutCreateMenu(MenuFunc);
	addMenu(IDM_CAXES);
	addMenu(IDM_CGRID);
	addMenu(IDM_CBKG_COLOR);
	addMenu(IDM_CLIGHT);
	addMenu(IDM_CLEGEND);
	addSepr;
	addMenu(IDM_CRESET);
	addMenu(IDM_CSAVE);
	addSepr;
	addMenu(IDM_CLEG_EDIT);
	addMenu(IDM_CMOD_APR);
	const int menu_render = glutCreateMenu(MenuFunc);
	addMenu(IDM_RSOL_SET);
	addSepr;
	addMenu(IDM_RDRAW_WIRE);
	addMenu(IDM_RDRAW_FILL);
	addMenu(IDM_RDRAW_CONT);
	addSepr;
	addMenu(IDM_RDRAW_FLOODED);
	addMenu(IDM_RDRAW_CUT);
	addMenu(IDM_RCUT_SETS);
	addSepr;
	//addMenu(IDM_RSCREEN_SAVE);
	const int menu_main = glutCreateMenu(MenuFunc);
	#ifndef _USE_FV_LIB
	addMenu(IDM_OPEN_MESH);
	addMenu(IDM_OPEN_FIELD);
	addSepr;
	#endif
	addMenu(IDM_RESET);
	addMenu(IDM_RELOAD);
	addSepr;
	glutAddSubMenu("View",menu_view);
	glutAddSubMenu("Configure",menu_config);
	glutAddSubMenu("Render",menu_render);
	addSepr;
	addMenu(IDM_QUIT);

	glutSetMenu(menu_main);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	glutIdleFunc(IdleFunc);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
				  GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	InitGL();
}

GlutWindow::~GlutWindow()
{
	// TODO Auto-generated destructor stub
}

void GlutWindow::DoInit(void)
{

}

void GlutWindow::RefreshScreen(void)
{
	glutPostRedisplay();
}

void GlutWindow::Show(bool show)
{
	if (show && _wndId)
	{
		glutSetWindow(_wndId);
	}
}

int GlutWindow::MainLoop(void)
{
	Run();
	return 0;
}

void GlutWindow::InitGL(void)
{
	glClearColor(0.0f,0.0f,0.0f,1.0f);
}


void GlutWindow::IdleFunc(void)
{
	ENTERCRITICALSECTION(crit_section);
	if (Window::GetInstance()->IsShutdown()) close_it++;
	else close_it = 0;
	LEAVECRITICALSECTION(crit_section);
	if (close_it) {
		close_it = 0;
		glutLeaveMainLoop();
	}
	else DisplayFunc();
}

void GlutWindow::DisplayFunc(void)
{
	ViewManagerInst().DisplayCallback();
	glutSwapBuffers();
}

void GlutWindow::ResizeFunc(int width,int height)
{
	glMatrixMode(GL_PROJECTION);

	// Resetujemy macierz projkecji
	glLoadIdentity();

	// Ustawiamy perspektywę
	gluPerspective(70.0, 1.0, 1.0, 20.0);

	// Korekta
	if(width <= height)
		glViewport(0, (height - width)/2, width, width);
	 else
		glViewport((width - height)/2, 0, height, height);

	// Ustawiamy macierz modelu
	glMatrixMode(GL_MODELVIEW);

	// Resetujemy macierz modelu
	glLoadIdentity();
}

void GlutWindow::MouseFunc(int button,int state,int x,int y)
{
	int modif = glutGetModifiers();
	ViewManagerInst().MouseButtonCallback(button,state,modif,x,y);
	glutPostRedisplay();
}

void GlutWindow::MouseMotionFunc(int x,int y)
{
	ViewManagerInst().MouseMotionCallback(x,y);
	glutPostRedisplay();
}

void GlutWindow::MenuFunc(int option)
{
	switch(option)
	{
	#ifndef _USE_FV_LIB
	case IDM_OPEN_MESH:
	case IDM_OPEN_FIELD:
		assert(!"Not implemented");
		break;
	#endif
	case IDM_REFRESH:
		ViewManagerInst().MenuCallback(IDM_RESET);
		break;
	case IDM_RELOAD:
	case IDM_RESET:
	case IDM_QUIT:
	{
		ENTERCRITICALSECTION(crit_section);
		Window::GetInstance()->SetShutdown();
		LEAVECRITICALSECTION(crit_section);
	} break;
	case IDM_VPERSP:
	case IDM_VORTHO:
	case IDM_VTOP:
	case IDM_VBOTTOM:
	case IDM_VFRONT:
	case IDM_VBACK:
	case IDM_VLEFT:
	case IDM_VRIGHT:
	case IDM_VDEFAULT:
	case IDM_VFULL:
	case IDM_VBBOX:
	case IDM_VFAST:
	case IDM_VNEW:
	case IDM_VNEXT:
	case IDM_VPREV:
	case IDM_VDUMP_CURR:
	case IDM_VDUMP_ALL:

	case IDM_CAXES:
	case IDM_CGRID:
	case IDM_CBKG_COLOR:
	case IDM_CLIGHT:
	case IDM_CLEGEND:
	case IDM_CRESET:
	case IDM_CSAVE:
	case IDM_CLEG_EDIT:
	case IDM_CMOD_APR:

	case IDM_RSOL_SET:
	case IDM_RDRAW_WIRE:
	case IDM_RDRAW_FILL:
	case IDM_RDRAW_CONT:
	case IDM_RDRAW_FLOODED:
	case IDM_RDRAW_CUT:
	case IDM_RCUT_SETS:
	case IDM_RSCREEN_SAVE:
		ViewManagerInst().MenuCallback(option);
		break;
	default:
		break;
	}

	glutPostRedisplay();
}

void GlutWindow::KeyboardFunc(unsigned char key,int x,int y)
{
	switch(key)
	{
	case 27:
		MenuFunc(IDM_QUIT);
	default:
		break;
	}

	glutPostRedisplay();
}


#undef addSepr
#undef addMenu
}
}
