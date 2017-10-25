/*
 * RContext.cpp
 *
 *  Created on: 28 sty 2014
 *      Author: dwg
 */
#include "fv_config.h"
#include "RContext.h"
#include "Geometry.h"
#include "Mesh.h"
#include "Field.h"
#include "Accelerators.h"
#include "ViewManager.h"


namespace FemViewer {

RContext::RContext()
: camera()
, mesh(NULL)
, field(NULL)
, imgw(imageWidth)
, imgh(imageHeight)
, backgroundColor(0.498f, 0.616f, 0.780f)
, image(NULL)
, objects()
, accelStruct(NULL)
{
	imageAspectRatio = imgw / (float)imgh;
	//accelStruct = new AccelerationStructure(this);
	//saccelStruct = new BVH(this);
}

RContext::RContext(const ViewManager& VMgr, Mesh* mesh_, Field* field_)
: camera(VMgr.currentView)
, mesh(mesh_)
, field(field_)
, imgw(VMgr.iWndWidth)
, imgh(VMgr.iWndHeight)
, image(NULL)
, backgroundColor(VMgr.UserSettings.BkgColor)
, objects()
, accelStruct(NULL)
{
	imageAspectRatio = imgw/(float)imgh;
}

RContext::RContext(const RContext& rc)
: camera(rc.camera)
, mesh(rc.mesh)
, field(rc.field)
, imgw(rc.imgw)
, imgh(rc.imgh)
, backgroundColor(rc.backgroundColor)
, image(NULL)
, accelStruct(NULL)
{
	imageAspectRatio = imgw / (float)imgh;
	//accelStruct = new AccelerationStructure(this);
	//saccelStruct = new BVH(this);
}

RContext::~RContext()
{
	while (!objects.empty()) {
		mfvBaseObject *obj = objects.back();
		objects.pop_back();
		delete obj;
	}
	if (accelStruct != NULL) {
		delete accelStruct; accelStruct = NULL;
	}
}

}// end namespace FemViewer



