#ifndef _VIEW_H_
#define _VIEW_H_

#include "Vec3D.h"
#include "MathHelper.h"
#include "Matrix.h"
#include "BBox3D.h"

//#include <iostream>



namespace FemViewer {

class BBox3D;
class Tile;
class ViewManager;
class Camera;


class View
{
	friend class Camera;
	typedef fvmath::CVec3f vec3d;
public:
	/* Type of projection in view */
	enum ProjType {
		ePerspective = 0,
		eParallel
	};
	/* Default constructor */
	View ();
	/* Copy Constructor */
	View (const View&   pView);
	
	View (const CVec3f&     pCenter);

	View (const CVec3f&     pCenter,
		  const CVec3f&     pDirection);

	View (const int	    pProjType,
		  const CVec3f& pCenter,
          const CVec3f& pDirection,
          const CVec3f& pUp,
          const float   pDistance,
	      const float   pFOVDegrees);

	View (const BBox3D& pBBox3D, const int pType);

	View (const BBox3D& pBBox3D, View& pReferenceView);

	View (const int     pProjType,
		  const BBox3D& pBBox3D,
	      const CVec3f& pDirection,
	      const CVec3f& pUp
		  );
	/* Destructor */
	~View();

	void     drawCenter () const;

	void     dump       (ViewManager&  pViewManager) const;

	void     dump       (std::ostream& pOstream) const;

	CVec3f    getPosition() const;

	void  	initCamera(const Tile& pTile) /*const*/;
	void    adjustToPlane(fvmath::CVec3f n,fvmath::CVec3f p,const float distane);
	void    setUp();

	 	  Matrix<float>& getCameraMatrix()       { return w2c; }
	const Matrix<float>& getCameraMatrix() const { return w2c; }

		  Matrix<float>& getProjectionMatrix()       { return c2p; }
    const Matrix<float>& getProjectionMatrix() const { return c2p; }

	View&    operator=  (const View&   pView);

	void     rotate     (const float   pDeltaX,
                       const float   pDeltaY);

	void     scaleFOV   (const float   pFOVFactor);

	void     setCenter   (const CVec3f& pCenter);

	void     translate  (const float   pDeltaX,
		                   const float   pDeltaY);

	const    CVec3f& getDirection() const;
    int      GetType() const { return type; }
	//float	 GetRadious() const;
	float    getDistance() const;

	float    getNearClipDistance() const;
	float    getFarClipDistance () const;

	void SwitchProjection(int projection);
public:
	/* Static and const camera params */
	static const CVec3f	vDefaultCenter;
	static const CVec3f	vDefaultDirection;
	static const float  fDeafultRadious;
	static const float  fDefaultDistance;
	static const float  fDefaultFOVDegrees;
	static const CVec3f	vDefaultUp;
	static const float  fInitialObjectDistanceFactor;
	static const float  fMaxFOVDegrees;
	static const float  fMinFOVDegrees;
	static const float	Identity[4][3];
	static const bool   bOrthoState;

  void     ajust(const BBox3D&  pBBox3D);

  void     init (const int	pProjType    = ePerspective,
		         const CVec3f&     pCenter     = vDefaultCenter,
                 const CVec3f&     pDirection  = vDefaultDirection,
                 const CVec3f&     pUp         = vDefaultUp,
                 const float      pDistance   = fDefaultDistance,
                 const float      pFOVDegrees = fDefaultFOVDegrees);



#ifdef FV_CHECK_INVARIANTS
  void invariants() const;
#endif // #ifdef GLV_BEBUG
  Matrix<float> c2p;
  Matrix<float> w2c;
  int type;
  CVec3f vCenter;
  CVec3f vDirection;
  CVec3f vUp;
  //float Radious;
  float Distance;
  float FOVDegrees;
  //float BBoxRadious;

  //bool bOrthoOn;

};

} // end namespace FemViewer

#endif /* _VIEW_H_
*/
