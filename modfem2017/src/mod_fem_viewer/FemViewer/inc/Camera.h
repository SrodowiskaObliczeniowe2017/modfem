/*
 * Camera.h
 *
 *  Created on: 28 sty 2014
 *      Author: dwg
 */

#ifndef _Camera_h__
#define _Camera_h__

#include <cmath>
#include "MathHelper.h"
#include "Matrix.h"
#include "View.h"



namespace FemViewer {



template<typename T>
inline T degtorad(const T &angle) { return angle * M_PI / T(180); }

struct Camera {

  float near, far;
  float fov;
  float angle;
  fvmath::CVec3f pos, target, up;
  Matrix<float> c2w, w2c;

  Camera(void) { set(); }

  Camera(const View& view) : w2c(view.w2c) {
	  near = view.getNearClipDistance();
	  far  = view.getFarClipDistance();
	  fov  = view.FOVDegrees;
	  angle = atan(degtorad(fov * 0.5f));
	  pos = view.getPosition();
	  target = view.vCenter;
	  up = view.vUp;
	  c2w = w2c.Invers();
  }

  Camera(const Camera& cam) {
	  near = cam.near;
	  far = cam.far;
	  fov = cam.fov;
	  angle = cam.angle;
	  pos = cam.pos;
	  target = cam.target;
	  c2w = cam.c2w;
	  w2c = c2w.Invers();

  }

  void set(Matrix<float>& cam2world = Matrix<float>::Identity,
		  const float f = 90,
		  const float nea = 0.1f,
		  const float fr = 1000.0f) {
	  near = nea;
	  far = fr;
	  fov = f;
	  angle = atan(degtorad(fov * 0.5f));
	  pos = fvmath::CVec3f(1.f,1.f,1.f);
	  target = View::vDefaultCenter;
	  up = View::vDefaultUp;
	  c2w = cam2world;
	  w2c = c2w.Invers();
  }

  void set(fvmath::CVec3f eye_,fvmath::CVec3f at_,fvmath::CVec3f up_) {
	  pos = eye_;
	  target = at_;
	  up = up_;

  }
};
}// end namespace FemViewer
#endif /* _Camera_h__ */
