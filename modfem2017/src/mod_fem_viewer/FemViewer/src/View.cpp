#include "../../utils/fv_assert.h"

#ifdef FV_DEBUG
#define FV_CHECK_INVARIANTS invariants()
#else
#define FV_CHECK_INVARIANTS
#endif // #ifdef FD_BEBUG

#include "fv_inc.h"
#include "View.h"
#include "BBox3D.h"
#include "Matrix.h"
#include "Tile.h"
#include "ViewManager.h"
#include "Log.h"

#include <stdio.h>
#include <vector>
#define min_(a,b) (((a < b)) ? (a) : (b))
#define max_(a,b) (((a > b)) ? (a) : (b))

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

#define SQRT0 0.70710678f
namespace FemViewer {

using ::fabs;
class ViewManager;
const CVec3f View::vDefaultCenter              (0.0f, 0.0f, 0.0f);
const CVec3f View::vDefaultDirection           (-0.5f, -0.5f, -SQRT0);
const float View::fDeafultRadious			  (5.0f);
const float View::fDefaultDistance            (15.0f);//(5.0f);
const float View::fDefaultFOVDegrees          (45.0f);//(40.0f);
const CVec3f View::vDefaultUp                  (-0.5f, -0.5f, SQRT0);
const float View::fInitialObjectDistanceFactor(3.0f);//(3.0f);
const float View::fMaxFOVDegrees              (90.0f);//(60.0f);
const float View::fMinFOVDegrees              (0.001f);//(0.001f);
const bool  View::bOrthoState				  (false);

View::View()
: c2p(0.0f)
, w2c()
{
	//mfp_debug("View dflt ctr\n");
	init();
	FV_CHECK_INVARIANTS;
}

View::View(const View& pView)
: c2p(pView.c2p)
, w2c(pView.w2c)
{
	//mfp_log_debug("View ctr of copy\n");
	init(pView.type,
		pView.vCenter,
		pView.vDirection,
		pView.vUp,
		pView.Distance,
		pView.FOVDegrees
		);

	FV_CHECK_INVARIANTS;
}

View::View(const CVec3f& pCenter)
: w2c(0.0f)
, c2p(.0f)
{
	mfp_log_debug("View ctr with CVec3f\n");
	init(ePerspective,pCenter);

	FV_CHECK_INVARIANTS;
}

View::View(const CVec3f& pCenter,
           const CVec3f& pDirection)
: c2p(0.0f)
, w2c()
{
	mfp_debug("View(const CVec3f&,const CVec3f&)\n");
	init(ePerspective,pCenter, pDirection);
	FV_CHECK_INVARIANTS;
}

View::View(const int pType,
		   const CVec3f& pCenter,
           const CVec3f& pDirection,
           const CVec3f& pUp,
           const float  pDistance,
           const float  pFOVDegrees)
: c2p(0.0f)
, w2c()
{
	mfp_log_debug("View(const CVec3f& pCenter...)\n");
    init(pType,
       pCenter,
       pDirection,
       pUp,
       pDistance,
       pFOVDegrees);

  FV_CHECK_INVARIANTS;
}

View::View(const BBox3D& pBBox3D, const int pType)
: c2p(0.0f)
, w2c()
{
	mfp_log_debug("View(const BBox3D&)\n");
    init(pType);
	ajust(pBBox3D);

	FV_CHECK_INVARIANTS;
}

// Special constructor enabling to build a view, while keeping the position
// of the previous view.  But, the rotation distance is updated to be in
// sync with the given BBox
View::View(const BBox3D& pBBox3D,
           View& pView)
: c2p(0.0f)
, w2c()
{
  init(pView.type,
		  pView.vCenter,
		  pView.vDirection,
		  pView.vUp,
		  pView.Distance,
		  pView.FOVDegrees);

  ajust(pBBox3D);

  vCenter = pView.vCenter;

  FV_CHECK_INVARIANTS;
}

View::View(const int pType,
		   const BBox3D& pBBox3D,
           const CVec3f&    pDirection,
           const CVec3f&    pUp)
: c2p(0.0f)
, w2c()
{
  init (ePerspective,
		CVec3f(),
        pDirection,
        pUp);

  ajust(pBBox3D);

  FV_CHECK_INVARIANTS;
}

View::~View()
{
  FV_CHECK_INVARIANTS;
}

// Draw a small orange axis at the center of the View
void View::drawCenter() const
{
  FV_CHECK_INVARIANTS;

  glBegin(GL_LINES);

  glColor3f(0.7f,0.5f,0.1f);

  const float lSize = static_cast<float>(0.05*Distance);
  const float lX    = static_cast<float>(vCenter[0]);
  const float lY    = static_cast<float>(vCenter[1]);
  const float lZ    = static_cast<float>(vCenter[2]);

  glVertex3f(lX-lSize, lY      , lZ      );
  glVertex3f(lX+lSize, lY      , lZ      );
  glVertex3f(lX      , lY-lSize, lZ      );
  glVertex3f(lX      , lY+lSize, lZ      );
  glVertex3f(lX      , lY      , lZ-lSize);
  glVertex3f(lX      , lY      , lZ+lSize);

  glEnd();
  FV_CHECK_INVARIANTS;
}

// Dump the view on as a status message in the ViewManager
void View::dump(ViewManager& pViewManager) const
{
  FV_CHECK_INVARIANTS;

  char lViewString[1024];
  sprintf(lViewString," view center = %+f %+f %+f", vCenter[0], vCenter[1], vCenter[2]);
  pViewManager.AddStatusMsg(lViewString);
  sprintf(lViewString," view direction = %+f %+f %+f", vDirection[0], vDirection[1], vDirection[2]);
  pViewManager.AddStatusMsg(lViewString);
  sprintf(lViewString," view up direction =  %+f %+f %+f", vUp[0], vUp[1], vUp[2]);
  pViewManager.AddStatusMsg(lViewString);
  sprintf(lViewString," distance from center = %+f", Distance);
  pViewManager.AddStatusMsg(lViewString);
  sprintf(lViewString," field of view in degrees = %+f", FOVDegrees);
  pViewManager.AddStatusMsg(lViewString);

  FV_CHECK_INVARIANTS;
}

// Dump the view on an std::ostream
void View::dump(std::ostream& pOstream) const
{
  FV_CHECK_INVARIANTS;

  char lViewString[1024];
  sprintf(lViewString,"view %+f %+f %+f %+f %+f %+f %+f %+f %+f %+f %+f",
          vCenter[0], vCenter[1], vCenter[2],
          vDirection[0], vDirection[1], vDirection[2],
          vUp[0], vUp[1], vUp[2],
          Distance,
          FOVDegrees);

  pOstream << lViewString << std::endl;

  FV_CHECK_INVARIANTS;
}

// Returns the position of the camera
CVec3f View::getPosition() const
{
  FV_CHECK_INVARIANTS;

  return  vCenter - Distance*vDirection;
}

void View::initCamera(const Tile& pTile) // const
{
	FV_CHECK_INVARIANTS;

	if(type == eParallel) {
		pTile.InitiOrtho3DMatrix(FOVDegrees, Distance, 2.0*Distance,c2p.matrix.data());
	}
	else {
		pTile.InitFrustrumMatrix(FOVDegrees,getNearClipDistance(), getFarClipDistance(),c2p.matrix.data());
    }

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	const CVec3f lEye = vCenter - Distance*vDirection;
	w2c = Matrix<float>::LookAt(
			lEye,
			CVec3f(vCenter.x,vCenter.y,vCenter.z),
			CVec3f(vUp.x,vUp.y,vUp.z));
	glLoadMatrixf(w2c.matrix.data());
	//std::cout << w2c << std::endl;
	FV_CHECK_INVARIANTS;
}

void View::adjustToPlane(fvmath::CVec3f n,fvmath::CVec3f p,const float distane)
{
	// Just in case
	Normalize(n);
	CVec3f lEye = p + n*distane;
	w2c = Matrix<float>::LookAt(
			lEye, vCenter, vUp);
}

void View::setUp()
{
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();

	const CVec3f lEye = vCenter - Distance*vDirection;
	w2c = Matrix<float>::LookAt(
				lEye,
				CVec3f(vCenter.x,vCenter.y,vCenter.z),
				CVec3f(vUp.x,vUp.y,vUp.z));
	//glLoadMatrixf(w2c.matrix.data());

	//FV_CHECK_INVARIANTS;
}

View& View::operator=(const View& pView)
{
  FV_CHECK_INVARIANTS;

  if (&pView != this) {
	  c2p = pView.c2p;
	  w2c = pView.w2c;
	  init(pView.type,
    	 pView.vCenter,
         pView.vDirection,
         pView.vUp,
         pView.Distance,
         pView.FOVDegrees);
  }

  FV_CHECK_INVARIANTS;
  return *this;
}

void View::translate(const float pDeltaX, const float pDeltaY)
{
  FV_CHECK_INVARIANTS;

  const CVec3f lOrtho = vDirection * vUp;
  const float lCoef  = Distance*FOVDegrees/45.0f;

  vCenter -= lCoef*pDeltaX*lOrtho;
  vCenter -= lCoef*pDeltaY*vUp;

  FV_CHECK_INVARIANTS;
}

void View::rotate(const float pDeltaX, const float pDeltaY)
{
  FV_CHECK_INVARIANTS;

  const CVec3f lOrtho = vDirection * vUp;

  Matrix<float> lMatrixY;
  lMatrixY.rotateAbout(lOrtho.v, pDeltaY);

  vUp        = lMatrixY*vUp;
  vDirection = lMatrixY*vDirection;

  Matrix<float> lMatrixX;
  lMatrixX.rotateAbout(CVec3f(vUp), pDeltaX);

  vDirection = lMatrixX*vDirection;

  Normalize(vUp);
  Normalize(vDirection);


  FV_CHECK_INVARIANTS;
}

void View::scaleFOV(const float pFOVFactor)
{
  FV_CHECK_INVARIANTS;

  FOVDegrees *= pFOVFactor;
  FOVDegrees = max_(FOVDegrees, fMinFOVDegrees);
  FOVDegrees = min_(FOVDegrees, fMaxFOVDegrees);

  FV_CHECK_INVARIANTS;
}

void View::setCenter(const CVec3f& pCenter)
{
	vCenter = pCenter;
}

// Used internally by all the constructors using a BBox3D
// This is called after the View::init function in those constructors
void View::ajust(const BBox3D& pBBox3D)
{
  FV_CHECK_INVARIANTS;

  vCenter       = pBBox3D.getCenter();
  float fRadius = pBBox3D.getCircumscribedSphereRadius();

  if(fRadius == 0.0f) {
    fRadius = 1.0f;
  }

  Distance = fRadius*fInitialObjectDistanceFactor;

  FV_CHECK_INVARIANTS;
}

// Used internally by all the constructors
void View::init(const int    pProjType,
				const CVec3f& pCenter,
                const CVec3f& pDirection,
                const CVec3f& pUp,
                const float  pDistance,
                const float  pFOVDegrees
				)
{
	type = pProjType;
	vCenter = pCenter;
	vDirection = pDirection;
	vUp = pUp;
	//Radious = pDistance / fInitialObjectDistanceFactor;
	Distance = pDistance;
	FOVDegrees = pFOVDegrees;

	Normalize(vDirection);

	// Make sure that the vDirection and vUp are perpendicular
	CVec3f tmp = vDirection * vUp;
	if(Norm(tmp) < 1.0E-3f) {

		// If vUp is close to being parallel to vDirection, then we make an
    // arbitraty choice of vUp direction
    // We take the one of the CVec3f (1,0,0), (0,1,0) and (0,0,1)
    // that has the a dotProduct with vDirection closest to zero
    const float lX    = vDirection.x;
    const float lY    = vDirection.y;
    const float lZ    = vDirection.z;
    const float lAbsX = fabs(lX);
    const float lAbsY = fabs(lY);
    const float lAbsZ = fabs(lZ);

    // Use Gram-Schmidt to get vUp perpendicular to vDirection
    // using the fact that vDirection.dotProduct(vDirection) == 1.0f
    assert(fabs(Norm(vDirection) - 1.0f) < 1.0E-6f);

    if (lAbsX <= lAbsY && lAbsX <= lAbsZ) {
      // lV = CVec3f(1.0, 0.0, 0.0)
      // vUp = lV - ((lV.dotProduct(vDirection))/vDirection.dotProduct(vDirection))*vDirection;
      assert(lAbsX == fv_min(fv_max(lAbsX,lAbsY), lAbsZ));
      vUp = CVec3f(1.0f - lX*lX, -lX*lY, -lX*lZ);
    }
    else if (lAbsY <= lAbsZ) {
      // lV = CVec3f(0.0, 1.0, 0.0)
      // vUp = lV - ((lV.dotProduct(vDirection))/vDirection.dotProduct(vDirection))*vDirection;
      FV_ASSERT(lAbsY == min_(min_(lAbsX,lAbsY), lAbsZ));
      vUp = CVec3f(-lY*lX, 1.0f  - lY*lY, -lY*lZ);
    }
    else {
      // lV = CVec3f(0.0, 0.0, 1.0)
      // vUp = lV - ((lV.dotProduct(vDirection))/vDirection.dotProduct(vDirection))*vDirection;
      FV_ASSERT(lAbsZ == min_(min_(lAbsX,lAbsY), lAbsZ));
      vUp = CVec3f(-lZ*lX, -lZ*lY, 1.0f - lZ*lZ);
    }

    Normalize(vUp);

    assert(fabs(Dot(vUp, vDirection)) < 1.0e-6f);
    assert(fabs(Norm(vDirection * vUp) - 1.0f) < 1.0e-6f);
  }

  FV_CHECK_INVARIANTS;
}



// Accessor
const CVec3f& View::getDirection() const
{
  FV_CHECK_INVARIANTS;

  return vDirection;
}



// Accessor
float View::getDistance() const
{
  FV_CHECK_INVARIANTS;

  return Distance;
}

float View::getNearClipDistance() const
{
	return Distance/10.0f;
}

float View::getFarClipDistance() const
{
	return Distance*10.0f;
}

void View::SwitchProjection(int projection)
{
	switch(projection)
	{
	case ePerspective:
		type = projection;
		c2p = Matrixf(0.0f);
		break;
	case eParallel:
		type = projection;
		c2p = Matrixf::Identity;
		break;
	default:
		assert(!"Unsupported projection type");
	}
}

//void View::CalcViewVolume(float ViewVolume[7]) const
//{
//	TViewVolume& vv = *(TViewVolume*) ViewVolume;
//	vv.hw = fImageWidth;
//	vv.hh = fImageHeight;
//    vv.zn = fZNear;
//	vv.zf = fZFar;
//    vv.iez = tan(fFOV*0.00872665)/min(fImageWidth,fImageHeight);
//    vv.tsx = vv.tsy = 0.0f;
//	
//	// Scale the width or height up to match the aspect ratio
//    if (fImageWidth >= fImageHeight)
//		vv.hw *= fImageWidth / fImageHeight;
//	else
//		vv.hh *= fImageHeight / fImageWidth;
//
//}
//
//union Vec3f { float v3[3]; struct { float x, y, z; }; };
//
//void View::ConfigureView(float ViewToWorld[4][3], float ViewVolume[7]) const
//{
//	GLfloat g_mv[16], g_proj[16];
//	glGetFloatv(GL_MODELVIEW, g_mv);
//	glGetFloatv(GL_PROJECTION, g_proj);
//	// Create the ViewToWorld transformation by combining g_ViewOrientation
//    // and g_ViewCenter.
//    for (int j = 0, i=0; i<12;j++, i+=4)
//        *(Vec3f*) ViewToWorld[j] = *(Vec3f*) g_mv[i];
//    *(Vec3f*) ViewToWorld[3] = *(Vec3f*) &vCenter[0];
//    ViewToWorld[3][0] = ViewToWorld[3][1] =  ViewToWorld[3][2] = 0;
//
//    // Declare TViewVolume
//    //    HalfWidth, HalfHeight, ZNear, ZFar, InverseEyeZ, TanSkewX, TanSkewY
//
//	// Use a reference to ViewVolume[]
//    TViewVolume& vv = *(TViewVolume*) ViewVolume;
//
//    vv.tsx = vv.tsy = 0;    // not using skewed views
//    vv.hw = vv.hh = 0.5f * min(fImageWidth,fImageHeight);
//    
//	// Increase the wider dimension
//	if (fImageWidth >= fImageHeight)
//		vv.hw *= fImageWidth / fImageHeight;
//    else
//		vv.hh *= fImageHeight / fImageWidth;
//
//    // Map the bounding sphere from world space to view space
//	Vec3f SphereCenterView = {vCenter[0], vCenter[1], vCenter[2], };
//    Vec3MultInverseRotTrn( SphereCenterView.v3, g_BoundingSphere, ViewToWorld );
//
//    float& radius = g_BoundingSphere[3];
//    vv.zn = SphereCenterView.x + radius;
//    vv.zf = SphereCenterView.x - radius;
//    static const float kVerySmall = 1e-6f;
//    if (g_ViewAngle < kVerySmall)
//        vv.iez = 0;     // parallel view
//    else
//    {
//        // Perspective view
//        // Ensure the near clipping plane is not too close to the eyepoint
//        // and the ratio between the near and far clipping plane distances
//        // provides adequate resolution for the z-buffer.
//        float EyeZ = min(vv.hw,vv.hh) / tan(0.5f*g_ViewAngle);
//        static const float kMinNearDistance = 1e-4f;
//        if (vv.zn > EyeZ-kMinNearDistance)
//            vv.zn = EyeZ-kMinNearDistance;
//        static const float kMinNearFarFactor = 5e-4f;
//        if (vv.zn > (EyeZ-vv.zf)*kMinNearFarFactor)
//            vv.zn = (EyeZ-vv.zf)*kMinNearFarFactor;
//        vv.iez = 1 / EyeZ;
//    }
//}


#ifdef FV_CHECK_INVARIANTS
void View::invariants() const
{
//	FV_ASSERT(fRadious > 0.0f);
//	FV_ASSERT(fDistance >  0.0f);
//	FV_ASSERT(fFOVDegrees >= fMinFOVDegrees);
//	FV_ASSERT(fFOVDegrees <= fMaxFOVDegrees);
//	FV_ASSERT(fabs(vDirection.Norm() - 1.0f) < 1.0e-6f);
//	FV_ASSERT(fabs(vUp.Norm() - 1.0f) < 1.0e-6f);
//	//FV_ASSERT(fabs(dotProd(vUp,vDirection)) < 1.0e-6f);
//	FV_ASSERT(fabs(((vDirection % vUp).Norm() - 1.0f)) < 1.0e-6f);
}
#endif // #ifdef FD_BEBUG

} // end namespace FemViewer
