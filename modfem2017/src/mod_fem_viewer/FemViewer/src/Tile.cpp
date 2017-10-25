#include "../../utils/fv_assert.h"
#include "fv_inc.h"
#include "Tile.h"
#include "Matrix.h"

#include <cmath>
#include <iostream>


namespace FemViewer {

#ifdef FV_DEBUG
#  define FV_CHECK_INVARIANTS invariants()
#else
#  define FV_CHECK_INVARIANTS
#endif
	
Tile::Tile(const int iImgWidth,
           const int iImgHeight,
           const float fAspectR)
  : fImageWidth(iImgWidth),
    fImageHeight(iImgHeight),
    fAspectRatio(fAspectR),
	fXmin(0.0f),
	fXmax(iImgWidth - 1.0f),
	fYmin(0.0f),
    fYmax(iImgHeight - 1.0f)
{
	//fMinHalfDim = fv_min(iImgWidth,iImgHeight) * 0.5f;
	  FV_ASSERT(iImgWidth  > 0);
	  FV_ASSERT(iImgHeight > 0);
	  FV_CHECK_INVARIANTS;
}

// A tile is not limited to the Image region, meaning
// that we can have:
//   pXMin < 0.0
//   pXMax > pImageWidth
//   pYMin < 0.0
//   pYMax > pImageHeight
Tile::Tile(const int iImgWidth,
           const int iImgHeight,
           const int iXMin,
           const int iXMax,
           const int iYMin,
           const int iYMax)
  : fImageWidth (iImgWidth),
	fImageHeight(iImgHeight),
    fXmin       (iXMin),
    fXmax       (iXMax),
	fYmin       (iYMin),
    fYmax       (iYMax)
{
	//fMinHalfDim = fv_min(iImgWidth,iImgHeight) / 2.0f;
	  
	FV_ASSERT(iImgHeight  > 0);
	FV_ASSERT(iImgWidth > 0);
	FV_ASSERT(iXMax > iXMin);
	FV_ASSERT(iYMax > iYMin);
	FV_CHECK_INVARIANTS;
}

Tile::~Tile()
{
  FV_CHECK_INVARIANTS;
}

void Tile::InitFrustrumMatrix(const float fFOV,
                              const float fZNear,
                              const float fZFar,
                              float mout[]) const
{
	FV_CHECK_INVARIANTS;
	FV_ASSERT(fFOV > 0.0);
	FV_ASSERT(fZNear > 0.0);
	FV_ASSERT(fZFar > fZNear);



	// Replace the call to gluPerspective to the
	// equivalent call to glFrustum in order to be
	// able to select a section of the original image
	//   gluPerspective(fFOVDegrees,
	//                  fImageWidth/fImageHeight,
	//                  fZNear,
	//                  fZFar);
	// The equivalent parameters to glFrustum are the following
	//   glFrustum(x_min, x_max, y_min, y_max, fZNear, fZfar);
	// with the following values:
	//const float aspectRatio = fImageWidth/fImageHeight;
	const float y_max        = fZNear * tan(fFOV*0.00872665);//M_PI/360);
	const float y_min        = -y_max;
	const float x_max        = y_max * fAspectRatio;
	const float x_min        = -x_max;

	// Compute new min/max to only display the tile
	// We have to keep in mind that OpenGL has the origin
	// at the lower left corner
	const float gltile_xmin = fXmin;
	const float gltile_xmax = fXmax + 1.0f;
	const float gltile_ymax = fImageHeight - fYmin;
	const float gltile_ymin = fImageHeight - fYmax - 1.0f;

	const float XMinTile = (gltile_xmin / fImageWidth )*(x_max-x_min) + x_min;
	const float XMaxTile = (gltile_xmax / fImageWidth )*(x_max-x_min) + x_min;
	const float YMinTile = (gltile_ymin / fImageHeight)*(y_max-y_min) + y_min;
	const float YMaxTile = (gltile_ymax / fImageHeight)*(y_max-y_min) + y_min;

	glMatrixMode(GL_PROJECTION);

	if (!mout) {
		glLoadIdentity();
		glFrustum(XMinTile, XMaxTile, YMinTile, YMaxTile, fZNear, fZFar);
	}
	else {
		mout[0] = 2.0f*fZNear / (XMaxTile - XMinTile);
		mout[5] = 2.0f*fZNear / (YMaxTile - YMinTile);
		mout[8] = (XMaxTile + XMinTile) / (XMaxTile - XMinTile);
		mout[9] = (YMaxTile + YMinTile) / (YMaxTile - YMinTile);
		mout[10] = -(fZFar + fZNear) / (fZFar - fZNear);
		mout[11] = -1.0f;
		mout[14] = -(2.0f*fZNear*fZFar) / (fZFar - fZNear);
		//glLoadIdentity();
		//MP.matrix.assign(mout,mout+16);
		glLoadMatrixf(mout);
		//glGetFloatv(GL_PROJECTION_MATRIX,MP.matrix.data());
		//std::cout << "MP: " << MP << std::endl;
	}
	FV_CHECK_INVARIANTS;
}

//void Tile::InitFrustrumMatrix(const float mt[]) const
//{
//	FV_ASSERT(mt! = NULL);
//	glMatrixMode(GL_PROJECTION);
//	glLoadMatrixf(mt);
//}

void Tile::InitOrtho2DMatrix() const
{
	FV_CHECK_INVARIANTS;

	// We have to keep in mind that OpenGL has the origin
	// at the lower left corner
	const float gltile_xmin = fXmin;
	const float gltile_xmax = fXmax + 1;
	const float gltile_ymax = fImageHeight - fYmin;
	const float gltile_ymin = fImageHeight - fYmax - 1;

	const float lXMinTile = (gltile_xmin / fImageWidth);
	const float lXMaxTile = (gltile_xmax / fImageWidth);
	const float lYMinTile = (gltile_ymin / fImageHeight);
	const float lYMaxTile = (gltile_ymax / fImageHeight);

	// To jest wczwesniej wywoływane
	// glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	gluOrtho2D(lXMinTile, lXMaxTile, lYMinTile, lYMaxTile);

	FV_CHECK_INVARIANTS;
}

void Tile::InitiOrtho3DMatrix(const float fFOV,
		                      const float fZNear,
		                      const float fZFar,
		                      float mout[]) const
{
	FV_CHECK_INVARIANTS;
	FV_ASSERT(fZNear > 0.0);
	FV_ASSERT(fZFar > fZNear);

	const float aspectRatio = fImageWidth/fImageHeight;
	const float y_max        = fZNear * tan(fFOV*0.00872665);//M_PI/360);
	const float y_min        = -y_max;
	const float x_max        = y_max * aspectRatio;
	const float x_min        = y_min * aspectRatio;

	// Compute new min/max to only display the tile
	// We have to keep in mind that OpenGL has the origin
	// at the lower left corner
	const float gltile_xmin = fXmin;
	const float gltile_xmax = fXmax + 1;
	const float gltile_ymax = fImageHeight - fYmin;
	const float gltile_ymin = fImageHeight - fYmax - 1;


	const float XMinTile = (gltile_xmin / fImageWidth )*(x_max-x_min) + x_min;
	const float XMaxTile = (gltile_xmax / fImageWidth )*(x_max-x_min) + x_min;
	const float YMinTile = (gltile_ymin / fImageHeight)*(y_max-y_min) + y_min;
	const float YMaxTile = (gltile_ymax / fImageHeight)*(y_max-y_min) + y_min;

	Matrixf MP;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (!mout) {

		glOrtho(XMinTile, XMaxTile, YMinTile, YMaxTile, -fZFar, fZFar);
	}
	else {
		mout[0] = 2.0f / (XMaxTile - XMinTile);
		mout[5] = 2.0f / (YMaxTile - YMinTile);
		mout[10] = 1.0f / (fZFar);
		mout[12] = -(XMaxTile + XMinTile) / (XMaxTile - XMinTile);
		mout[13] = -(YMaxTile + YMinTile) / (YMaxTile - YMinTile);
		glLoadMatrixf(mout);
	}

}

void Tile::InitViewport() const
{
  FV_CHECK_INVARIANTS;
  glViewport(0, 0, (GLsizei)(fXmax-fXmin+1), (GLsizei)(fYmax-fYmin+1));
  FV_CHECK_INVARIANTS;
}

#ifdef FV_DEBUG
void Tile::invariants() const
{
  FV_ASSERT(fImageWidth  >= 1.0);
  FV_ASSERT(fImageHeight >= 1.0);
  FV_ASSERT(fXmax > fXmin);
  FV_ASSERT(fYmax > fYmin);
}
#endif // #iFVef FV_BEBUG

} // end namespace FemViewer
