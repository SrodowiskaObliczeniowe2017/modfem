/*
 * RContext.h
 *
 *  Created on: 28 sty 2014
 *      Author: dwg
 */

#ifndef _RContext_h__
#define _RContext_h__


#include <vector>
#include "Enums.h"
#include "Camera.h"
#include "MathHelper.h"
//#include "Mesh.h"
#include "Color.h"


namespace FemViewer {

class Mesh;
class Field;
class Accelerator;
class Grid;
class mfvBaseObject;
class ViewManager;

template<typename T,typename U> class ArrayT;
template<typename T> class Array2;

class RContext
{
public:
	RContext(void);
	RContext(const ViewManager& VMgr,
			 Mesh* mesh_, // because we change elments
			 Field* field_);
	RContext(const RContext& rhs);
	~RContext(void);
	Camera camera;
	Mesh* mesh;
	Field* field;
	uint32_t imgw, imgh;
	float imageAspectRatio;
	ColorRGB backgroundColor;
	Array2<ColorRGB> *image;
	std::vector<mfvBaseObject *> objects;
	Grid *accelStruct;
};

}// end namespace FemViewer
#endif /* _RContext_h__ */
