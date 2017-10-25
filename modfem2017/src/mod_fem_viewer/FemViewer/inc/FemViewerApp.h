/*
 * FemViewerApp.h
 *
 *  Created on: 26-01-2012
 *      Author: Pawel Maciol
 */

#ifndef FEMVIEWERAPP_H_
#define FEMVIEWERAPP_H_


namespace FemViewer
{

class WindowFV;
class ModelControler;
class ViewManger;

class FemViewerApp
{
public:
	/// Constructor
	/// @param argc number of application arguments (for glut initialization)
	/// @param argv application arguments (for glut initialization)
	/// @param posX render window position X coordinate
	/// @param posY render window position Y coordinate
	/// @param width render window height
	/// @param title render window name
	FemViewerApp(int argc, char** argv, int posX, int posY,
				 int width, int height, const char* title);

	/// Destructor
	~FemViewerApp(void);

	/// Start function
	int Run(void);

	/// Check initialization
	bool IsInitialized(void) const
	{
		return(_running);
	}

	/// Check if appliaction is running
	int IsRunning() const { return _running; }

	/// Close aplliaction
	bool Close();

protected:

	/// GUI window
	WindowFV* 		_pwnd;

	/// Model handle
	ModelControler* _pmodel;


private:
	/// Running state flag for threads
	volatile int 	_running;

	/// @note hide copy constructor from interface user
	FemViewerApp(const FemViewerApp&);

	/// @note hide assignment operator from interface user
	FemViewerApp& operator=(const FemViewerApp&);



};

}


#endif /* FEMVIEWERAPP_H_ */
