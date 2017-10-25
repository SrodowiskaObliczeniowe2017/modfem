#ifndef _WINDOW_FV_H_
#define _WINDOW_FV_H_

#include "Common.h"


namespace FemViewer
{

// Forward declaretions
class ViewManager;
class ModelCtrl;

class mfvWindow //: public mfvBaseObject
{
//protected:
//	static Window * _self;
//public:
//	enum eMainOperation
//	{
//		IDM_SEPARATOR = -1,
//		#ifndef _USE_FV_LIB
//		IDM_OPEN_MESH = 0,
//		IDM_OPEN_FIELD,
//		IDM_REFRESH,
//		#else
//		IDM_REFRESH = 0,
//		#endif
//		IDM_RELOAD,
//		IDM_RESET,
//		IDM_QUIT,
//		IDM_VPERSP,
//		IDM_VORTHO,
//		IDM_VTOP,
//		IDM_VBOTTOM,
//		IDM_VFRONT,
//		IDM_VBACK,
//		IDM_VLEFT,
//		IDM_VRIGHT,
//		IDM_VDEFAULT,
//		IDM_VFULL,
//		IDM_VFAST,
//		IDM_VBBOX,
//		IDM_VNEW,
//		IDM_VNEXT,
//		IDM_VPREV,
//		IDM_VDUMP_CURR,
//		IDM_VDUMP_ALL,
//		IDM_CAXES,
//		IDM_CGRID,
//		IDM_CBKG_COLOR,
//		IDM_CLIGHT,
//		IDM_CLEGEND,
//		IDM_CRESET,
//		IDM_CSAVE,
//		IDM_CLEG_EDIT,
//		IDM_CMOD_APR,
//		IDM_RSOL_SET,
//		IDM_RDRAW_WIRE,
//		IDM_RDRAW_FILL,
//		IDM_RDRAW_CONT,
//		IDM_RDRAW_FLOODED,
//		IDM_RDRAW_CUT,
//		IDM_RCUT_SETS,
//		IDM_RSCREEN_SAVE,
//		IDM_HELP,
//		//IDM_SEPARATOR,
//	};
public:
	static bool init(int argc,char **argv);

	//int  IsShutdown() const { return this->_shutdown; }
	//void SetShutdown(){ this->_shutdown++; }
	// Destructor
	virtual ~mfvWindow(void);
protected:	//static ViewManager _vmgr;
	ViewManager & m_pview;
	ModelCtrl   & m_pmodel;
	// Constructor
	mfvWindow(void);





private:
	/// Not implemented, block to use it
	mfvWindow(const mfvWindow&);
	mfvWindow& operator =(const mfvWindow&);
};

} // end namespace FemViewer

#endif // WINDOW_H
