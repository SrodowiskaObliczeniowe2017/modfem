//#define NOMINMAX

#include "PrimitiveAccumulator.h"
#include "../../utils/fv_assert.h"
#include "../../include/fv_limits.h"
// system
#include <algorithm>
//#include <limits>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <stdexcept>


namespace FemViewer {
#ifndef M_PI
// PI is not defined in visual c headers
#define M_PI 3.14159265358979323846
#endif
#undef min
#undef max

#ifdef FV_DEBUG
inline bool Compare(const Vec3D& p1, const Vec3D& p2, const Vec3D& p3, const Vec3D& p4) {
	return (
		(p1 == p2) ||
		(p2 == p3) ||
		(p3 == p4) ||
		(p4 == p1) ||
		(p4 == p2) ||
		(p3 == p1)
		);
}
#endif

PrimitiveAccumulator::PrimitiveAccumulator(const bool pCreateSimplified)
  : aBBox3D                       (),
    aLines                             (),
    aLinesBBoxCounter                  (0),
    aLinesColored                      (),
    aLinesColoredBBoxCounter           (0),
    aPoints                            (),
    aPointsBBoxCounter                 (0),
    aPointsColored                     (),
    aPointsColoredBBoxCounter          (0),
    aQuads                             (),
    aQuadsBBoxCounter                  (0),
    aQuadsColored                      (),
    aQuadsColoredBBoxCounter           (0),
    aQuadsNormals                      (),
    aQuadsNormalsBBoxCounter           (0),
    aQuadsNormalsColored               (),
    aQuadsNormalsColoredBBoxCounter    (0),
    aSimplified                        (),
    aSimplifiedDirty                   (true),
    aTriangles                         (),
    aTrianglesBBoxCounter              (0),
    aTrianglesColored                  (),
    aTrianglesColoredBBoxCounter       (0),
    aTrianglesNormals                  (),
    aTrianglesNormalsBBoxCounter       (0),
    aTrianglesNormalsColored           (),
	aTrianglesNormalsColoredBBoxCounter(0),
    aTriElemsNormalColor					   (),
	aTriElemsNormalColorBBoxCounter		   (0),
	iPrimitiveOptimizerValue           (100)
{
  if (pCreateSimplified) {
    aSimplified = new PrimitiveAccumulator(false);
  }
  else {
    aSimplified = this;
  }
  FV_ASSERT(aSimplified != 0);
}

PrimitiveAccumulator::~PrimitiveAccumulator()
{
  FV_ASSERT(aSimplified != 0);

  if (aSimplified != this) {
    delete aSimplified;
  }
}

void PrimitiveAccumulator::AddArrow(const Vec3D& p1, const Vec3D& p2,
                                    const float tipprop, const int nrpolygons)
{
	FV_ASSERT(tipprop  > 0.0f);
	FV_ASSERT(tipprop <= 1.0f);
	FV_ASSERT(nrpolygons >= 1);

	AddLine(p1, p2);

	// Here it is no matter the 3rd params
	AddArrowTip(p1, p2, ColorRGB(), false, tipprop, nrpolygons);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddArrowColored(const Vec3D& p1, const ColorRGB& col1,
                                           const Vec3D& p2, const ColorRGB& col2,
                                           const float tipprop,
                                           const int nrpolygons)
{
	FV_ASSERT(tipprop  > 0.0f);
	FV_ASSERT(tipprop <= 1.0f);
	FV_ASSERT(nrpolygons >= 1);

	AddLineColored(p1, col1, p2, col2);
	AddArrowTip(p1, p2, col2, true, tipprop, nrpolygons);

	aSimplifiedDirty = true;
}

// We now do the equivalent of the glut commands in
// the PrimitiveAccumulator to take adventage of all
// it has to offer (BoundingBox updated, optimized diplay, etc)
//void PrimitiveAccumulator::AddSolidCube(const Vec3D& pCenter,
//                                        const float     pHalfSize)
//{
//  const float lMinX = pCenter._x() - pHalfSize;
//  const float lMinY = pCenter._y() - pHalfSize;
//  const float lMinZ = pCenter._z() - pHalfSize;
//  const float lMaxX = pCenter._x() + pHalfSize;
//  const float lMaxY = pCenter._y() + pHalfSize;
//  const float lMaxZ = pCenter._z() + pHalfSize;
//
//  Quad lQuad;
//
//  lQuad.pt[0] = Vec3D(lMinX, lMinY, lMinZ);
//  lQuad.pt[1] = Vec3D(lMaxX, lMinY, lMinZ);
//  lQuad.pt[2] = Vec3D(lMaxX, lMaxY, lMinZ);
//  lQuad.pt[3] = Vec3D(lMinX, lMaxY, lMinZ);
//  aQuads.push_back(lQuad);
//
//  lQuad.pt[0] = Vec3D(lMinX, lMinY, lMaxZ);
//  lQuad.pt[1] = Vec3D(lMaxX, lMinY, lMaxZ);
//  lQuad.pt[2] = Vec3D(lMaxX, lMaxY, lMaxZ);
//  lQuad.pt[3] = Vec3D(lMinX, lMaxY, lMaxZ);
//  aQuads.push_back(lQuad);
//
//  lQuad.pt[0] = Vec3D(lMinX, lMinY, lMinZ);
//  lQuad.pt[1] = Vec3D(lMaxX, lMinY, lMinZ);
//  lQuad.pt[2] = Vec3D(lMaxX, lMinY, lMaxZ);
//  lQuad.pt[3] = Vec3D(lMinX, lMinY, lMaxZ);
//  aQuads.push_back(lQuad       );
//
//  lQuad.pt[0] = Vec3D(lMinX, lMaxY, lMinZ);
//  lQuad.pt[1] = Vec3D(lMaxX, lMaxY, lMinZ);
//  lQuad.pt[2] = Vec3D(lMaxX, lMaxY, lMaxZ);
//  lQuad.pt[3] = Vec3D(lMinX, lMaxY, lMaxZ);
//  aQuads.push_back(lQuad);
//
//  lQuad.pt[0] = Vec3D(lMinX, lMinY, lMinZ);
//  lQuad.pt[1] = Vec3D(lMinX, lMaxY, lMinZ);
//  lQuad.pt[2] = Vec3D(lMinX, lMaxY, lMaxZ);
//  lQuad.pt[3] = Vec3D(lMinX, lMinY, lMaxZ);
//  aQuads.push_back(lQuad);
//
//  lQuad.pt[0] = Vec3D(lMaxX, lMinY, lMinZ);
//  lQuad.pt[1] = Vec3D(lMaxX, lMaxY, lMinZ);
//  lQuad.pt[2] = Vec3D(lMaxX, lMaxY, lMaxZ);
//  lQuad.pt[3] = Vec3D(lMaxX, lMinY, lMaxZ);
//  aQuads.push_back(lQuad);
//
//  aSimplifiedDirty = true;
//}

// We now do the equivalent of the glut commands in
// the PrimitiveAccumulator to take adventage of all
// it has to offer (BoundingBox updated, optimized diplay, etc)
//void PrimitiveAccumulator::addSolidSphere(const Vec3D& pCenter,
//                                          const float     pRadius,
//                                          const int       pSlices,
//                                          const int       pStacks)
//{
//
//  const double lDeltaTheta = 2.0f*M_PI/static_cast<double>(pSlices);
//  const double lDeltpthi   = M_PI/static_cast<double>(pStacks);
//
//  for (int i=0; i<pSlices; ++i) {
//    const double lTheta        = static_cast<double>(i)*lDeltaTheta;
//    const double lNextTheta    = lTheta + lDeltaTheta;
//    const double lCosTheta     = cos(lTheta);
//    const double lCosNextTheta = cos(lNextTheta);
//    const double lSinTheta     = sin(lTheta);
//    const double lSinNextTheta = sin(lNextTheta);
//
//    for (int j=0; j<pStacks; ++j) {
//      const double lPhi        = static_cast<double>(j)*lDeltpthi;
//      const double lNextPhi    = lPhi + lDeltpthi;
//      const double lCosPhi     = cos(lPhi);
//      const double lCosNextPhi = cos(lNextPhi);
//      const double lSinPhi     = sin(lPhi);
//      const double lSinNextPhi = sin(lNextPhi);
//
//      QuadNormals lQuad;
//      lQuad.aN1 = Vec3D(lCosNextTheta*lSinPhi    , lSinNextTheta*lSinPhi    , lCosPhi);
//      lQuad.aN2 = Vec3D(lCosTheta*lSinPhi        , lSinTheta*lSinPhi        , lCosPhi);
//      lQuad.aN3 = Vec3D(lCosTheta*lSinNextPhi    , lSinTheta*lSinNextPhi    , lCosNextPhi);
//      lQuad.aN4 = Vec3D(lCosNextTheta*lSinNextPhi, lSinNextTheta*lSinNextPhi, lCosNextPhi);
//
//      lQuad.pt[0] = pRadius*lQuad.aN1 + pCenter;
//      lQuad.pt[1] = pRadius*lQuad.aN2 + pCenter;
//      lQuad.pt[2] = pRadius*lQuad.aN3 + pCenter;
//      lQuad.pt[3] = pRadius*lQuad.aN4 + pCenter;
//
//      aQuadsNormals.push_back(lQuad);
//    }
//  }
//
//  aSimplifiedDirty = true;
//}

// We now do the equivalent of the glut commands in
// the PrimitiveAccumulator to take adventage of all
// it has to offer (BoundingBox updated, optimized diplay, etc)
//void PrimitiveAccumulator::addWireCube(const Vec3D& pCenter,
//                                       const float     pHalfSize)
//{
//  const float lMinX = pCenter._x() - pHalfSize;
//  const float lMinY = pCenter._y() - pHalfSize;
//  const float lMinZ = pCenter._z() - pHalfSize;
//  const float lMaxX = pCenter._x() + pHalfSize;
//  const float lMaxY = pCenter._y() + pHalfSize;
//  const float lMaxZ = pCenter._z() + pHalfSize;
//
//  Line lLine;
//  lLine.pt[0] = Vec3D(lMinX,lMinY,lMinZ); lLine.pt[1] = Vec3D(lMaxX,lMinY,lMinZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMinX,lMinY,lMinZ); lLine.pt[1] = Vec3D(lMinX,lMaxY,lMinZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMinX,lMinY,lMinZ); lLine.pt[1] = Vec3D(lMinX,lMinY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMaxX,lMinY,lMinZ); lLine.pt[1] = Vec3D(lMaxX,lMaxY,lMinZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMaxX,lMinY,lMinZ); lLine.pt[1] = Vec3D(lMaxX,lMinY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMaxX,lMaxY,lMinZ); lLine.pt[1] = Vec3D(lMinX,lMaxY,lMinZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMaxX,lMaxY,lMinZ); lLine.pt[1] = Vec3D(lMaxX,lMaxY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMinX,lMaxY,lMinZ); lLine.pt[1] = Vec3D(lMinX,lMaxY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMinX,lMinY,lMaxZ); lLine.pt[1] = Vec3D(lMaxX,lMinY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMinX,lMinY,lMaxZ); lLine.pt[1] = Vec3D(lMinX,lMaxY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMaxX,lMinY,lMaxZ); lLine.pt[1] = Vec3D(lMaxX,lMaxY,lMaxZ);
//  aLines.push_back(lLine);
//  lLine.pt[0] = Vec3D(lMinX,lMaxY,lMaxZ); lLine.pt[1] = Vec3D(lMaxX,lMaxY,lMaxZ);
//  aLines.push_back(lLine);
//
//  aSimplifiedDirty = true;
//}


void PrimitiveAccumulator::AddLine(const Vec3D& p1, const Vec3D& p2)
{
	Line l;
	l.pt[0] = p1; l.pt[1] = p2;

	aLines.push_back(l);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddLineColored(const Vec3D& p1, const ColorRGB& col1,
                                          const Vec3D& p2, const ColorRGB& col2)
{
	LineColored lc;
	lc.pt[0] = p1; lc.pt[1] = p2;
	lc.col1 = col1; lc.col2 = col2;

	aLinesColored.push_back(lc);
	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddPoint(const Vec3D& pt)
{
	Point p;
	p.pt = pt;

	aPoints.push_back(p);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddPointColored(const Vec3D& pt, const ColorRGB& col)
{
	PointColored p;
	p.pt = pt; p.col = col;

	aPointsColored.push_back(p);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddQuad(const Vec3D& p1, const Vec3D& p2,
                                   const Vec3D& p3, const Vec3D& p4)
{
	Quad q;
	q.pt[0] = p1; q.pt[1] = p2;
	q.pt[2] = p3; q.pt[3] = p4;

	aQuads.push_back(q);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddQuadColored(const Vec3D& p1, const ColorRGB& col1,
                                          const Vec3D& p2, const ColorRGB& col2,
                                          const Vec3D& p3, const ColorRGB& col3,
                                          const Vec3D& p4, const ColorRGB& col4)
{
	QuadColored qc;
	qc.pt[0] = p1; 
	qc.pt[1] = p2; 
	qc.pt[2] = p3; 
	qc.pt[3] = p4;

	//FV_ASSERT(!Compare(p1,p2,p3,p4));
	qc.col1 = col1; qc.col2 = col2; qc.col3 = col3; qc.col4 = col4;

	aQuadsColored.push_back(qc);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddTriangle(const Vec3D& p1, const Vec3D& p2,
                                       const Vec3D& p3)
{
	Triangle t;
	t.pt[0] = p1; t.pt[1] = p2; t.pt[2] = p3;

	aTriangles.push_back(t);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddTriangleColored(const Vec3D& p1, const ColorRGB& col1,
                                              const Vec3D& p2, const ColorRGB& col2,
                                              const Vec3D& p3, const ColorRGB& col3)
{
	TriangleColored tc;
	tc.pt[0] = p1;	tc.pt[1] = p2; tc.pt[2] = p3;
	tc.col1 = col1; tc.col2 = col2; tc.col3 = col3;

	aTrianglesColored.push_back(tc);

	aSimplifiedDirty = true;
}

void PrimitiveAccumulator::AddTriangleNormalColored(const float nv[], const float p1[], const float p2[], const float p3[],
									  const ColorRGB& col1, const ColorRGB& col2, const ColorRGB& col3)
{
	TriEelementNormalColor tnc;
	for(register int i(0);i<3;++i) {
		tnc.pt[0].v[i] = p1[i];
		tnc.pt[1].v[i] = p2[i];
		tnc.pt[2].v[i] = p3[i];
		tnc.nv.v[i]    = nv[i];
	}
	
	tnc.c1 = col1; tnc.c2 = col2; tnc.c3 = col3;
	
	aTriElemsNormalColor.push_back(tnc);

	aSimplifiedDirty = true;
}

#ifdef FV_DUMP_MEMORY_USAGE

template <class Container>
inline
std::string getStringSizeAndCptacity(const Container& pContainer)
{
  char lSizeAndCptacity[128];
  sprintf(lSizeAndCptacity,"%d/%d",
          sizeof(typename Container::value_type)*pContainer.size(),
          sizeof(typename Container::value_type)*pContainer.cptacity());
  return std::string(lSizeAndCptacity);
}

#endif // #ifdef  FV_DUMP_MEMORY_USAGE


// Dump in ASCII the caracteristics of the PrimitiveAccumulator
void PrimitiveAccumulator::DumpCharacteristics(std::ostream&       pOstream,
                                               const std::string&  pIndentation,
                                               const Matrix<float>&    pTransformation)
{
  pOstream << pIndentation << "PrimitiveAccumulator "  << std::endl;
  std::string lIndentation = pIndentation + "  ";

  // Update the bounding box. But we use the attribute
  // directly just after the call because we need it non-const
  GetBBox3D();
  aBBox3D.dumpCharacteristics(pOstream, lIndentation, pTransformation);

#ifdef GLV_DUMP_SIMPLICATION
  if (aSimplified != this) {
    pOstream << lIndentation << "Simplified "  << std::endl;
    aSimplified->dumpCharacteristics(pOstream, lIndentation + "  ", pTransformation);
  }
#endif // #ifdef GLV_DUMP_SIMPLICATION

#ifdef GLV_DUMP_MEMORY_USAGE
  {
    pOstream << lIndentation << "Memory used by the PrimitiveAccumulator = " << sizeof(*this) << std::endl;

    pOstream << lIndentation << "Memory used by ptoints                  = " << getStringSizeAndCptacity(ptoints                 ) << std::endl;
    pOstream << lIndentation << "Memory used by ptointsColored           = " << getStringSizeAndCptacity(ptointsColored          ) << std::endl;
    pOstream << lIndentation << "Memory used by aLines                   = " << getStringSizeAndCptacity(aLines                  ) << std::endl;
    pOstream << lIndentation << "Memory used by aLinesColored            = " << getStringSizeAndCptacity(aLinesColored           ) << std::endl;
    pOstream << lIndentation << "Memory used by aTriangles               = " << getStringSizeAndCptacity(aTriangles              ) << std::endl;
    pOstream << lIndentation << "Memory used by aTrianglesColored        = " << getStringSizeAndCptacity(aTrianglesColored       ) << std::endl;
    pOstream << lIndentation << "Memory used by aTrianglesNormals        = " << getStringSizeAndCptacity(aTrianglesNormals       ) << std::endl;
    pOstream << lIndentation << "Memory used by aTrianglesNormalsColored = " << getStringSizeAndCptacity(aTrianglesNormalsColored) << std::endl;
    pOstream << lIndentation << "Memory used by aQuads                   = " << getStringSizeAndCptacity(aQuads                  ) << std::endl;
    pOstream << lIndentation << "Memory used by aQuadsColored            = " << getStringSizeAndCptacity(aQuadsColored           ) << std::endl;
    pOstream << lIndentation << "Memory used by aQuadsNormals            = " << getStringSizeAndCptacity(aQuadsNormals           ) << std::endl;
    pOstream << lIndentation << "Memory used by aQuadsNormalsColored     = " << getStringSizeAndCptacity(aQuadsNormalsColored    ) << std::endl;

  }
#endif // #ifdef GLV_DUMP_MEMORY_USAGE

  SizeType lNbQuads            = aQuads           .size() + aQuadsNormals           .size();
  SizeType lNbQuadsColored     = aQuadsColored    .size() + aQuadsNormalsColored    .size();
  SizeType lNbTriangles        = aTriangles       .size() + aTrianglesNormals       .size();
  SizeType lNbTrianglesColored = aTrianglesColored.size() + aTrianglesNormalsColored.size();

  pOstream << lIndentation << "Number of point            = " << aPoints       .size() << std::endl;
  pOstream << lIndentation << "Number of point_colored    = " << aPointsColored.size() << std::endl;
  pOstream << lIndentation << "Number of line             = " << aLines        .size() << std::endl;
  pOstream << lIndentation << "Number of line_colored     = " << aLinesColored .size() << std::endl;
  pOstream << lIndentation << "Number of triangle         = " << lNbTriangles          << std::endl;
  pOstream << lIndentation << "Number of triangle_colored = " << lNbTrianglesColored   << std::endl;
  pOstream << lIndentation << "Number of quad             = " << lNbQuads              << std::endl;
  pOstream << lIndentation << "Number of quad_colored     = " << lNbQuadsColored       << std::endl;
}

const BBox3D& PrimitiveAccumulator::GetBBox3D() const
{

  // Scan only the new primitives for each type
  {
    const SizeType lSize = aLines.size();

    while (aLinesBBoxCounter < lSize) {
      const Line& lLine = aLines[aLinesBBoxCounter];
      aBBox3D += lLine.pt[0];
      aBBox3D += lLine.pt[1];
      ++aLinesBBoxCounter;
    }
  }

  {
    const SizeType lSize = aLinesColored.size();

    while (aLinesColoredBBoxCounter < lSize) {
      const Line& lLine = aLinesColored[aLinesColoredBBoxCounter];
	  if(fabs(lLine.pt[0][0]) > 1.0f){
		  printf("%d %f\n",(int)aLinesColoredBBoxCounter,lLine.pt[0][0]);
		  exit(0);
	  }
	  if(fabs(lLine.pt[0][1]) > 1.0f){
		  printf("%d %f\n",(int)aLinesColoredBBoxCounter,lLine.pt[0][1]);
		  exit(0);
	  }
      aBBox3D += lLine.pt[0];
      aBBox3D += lLine.pt[1];
      ++aLinesColoredBBoxCounter;
    }
  }

  {
    const SizeType lSize = aPoints.size();

    while (aPointsBBoxCounter < lSize) {
      const Point& lPoint = aPoints[aPointsBBoxCounter];
	  if(fabs(lPoint.pt[0]) > 1.0f){
		  printf("%d %f\n",aPointsBBoxCounter,lPoint.pt[0]);
		  exit(0);
	  }
      aBBox3D += lPoint.pt;
      ++aPointsBBoxCounter;
    }
  }

  {
    const SizeType lSize = aPointsColored.size();

    while (aPointsColoredBBoxCounter < lSize) {
      const Point& lPoint = aPointsColored[aPointsColoredBBoxCounter];
      aBBox3D += lPoint.pt;
      ++aPointsColoredBBoxCounter;
    }
  }

  {
    const SizeType lSize = aQuads.size();

    while (aQuadsBBoxCounter < lSize) {
      const Quad& lQuad = aQuads[aQuadsBBoxCounter];
      aBBox3D += lQuad.pt[0];
      aBBox3D += lQuad.pt[1];
      aBBox3D += lQuad.pt[2];
      aBBox3D += lQuad.pt[3];
      ++aQuadsBBoxCounter;
    }
  }

  {
    const SizeType lSize = aQuadsColored.size();

    while (aQuadsColoredBBoxCounter < lSize) {
      const Quad& lQuad = aQuadsColored[aQuadsColoredBBoxCounter];
      aBBox3D += lQuad.pt[0];
      aBBox3D += lQuad.pt[1];
      aBBox3D += lQuad.pt[2];
      aBBox3D += lQuad.pt[3];
      ++aQuadsColoredBBoxCounter;
    }
  }

  {
    const SizeType lSize = aQuadsNormals.size();

    while (aQuadsNormalsBBoxCounter < lSize) {
      const Quad& lQuad = aQuadsNormals[aQuadsNormalsBBoxCounter];
      aBBox3D += lQuad.pt[0];
      aBBox3D += lQuad.pt[1];
      aBBox3D += lQuad.pt[2];
      aBBox3D += lQuad.pt[3];
      ++aQuadsNormalsBBoxCounter;
    }
  }

  {
    const SizeType lSize = aQuadsNormalsColored.size();

    while (aQuadsNormalsColoredBBoxCounter < lSize) {
      const Quad& lQuad = aQuadsNormalsColored[aQuadsNormalsColoredBBoxCounter];
      aBBox3D += lQuad.pt[0];
      aBBox3D += lQuad.pt[1];
      aBBox3D += lQuad.pt[2];
      aBBox3D += lQuad.pt[3];
      ++aQuadsNormalsColoredBBoxCounter;
    }
  }

  {
    const SizeType lSize = aTriangles.size();

    while (aTrianglesBBoxCounter < lSize) {
      const Triangle& lTriangle = aTriangles[aTrianglesBBoxCounter];
      aBBox3D += lTriangle.pt[0];
      aBBox3D += lTriangle.pt[1];
      aBBox3D += lTriangle.pt[2];
      ++aTrianglesBBoxCounter;
    }
  }

  {
    const SizeType lSize = aTrianglesColored.size();

    while (aTrianglesColoredBBoxCounter < lSize) {
      const Triangle& lTriangle = aTrianglesColored[aTrianglesColoredBBoxCounter];
      aBBox3D += lTriangle.pt[0];
      aBBox3D += lTriangle.pt[1];
      aBBox3D += lTriangle.pt[2];
      ++aTrianglesColoredBBoxCounter;
    }
  }

  {
    const SizeType lSize = aTrianglesNormals.size();

    while (aTrianglesNormalsBBoxCounter < lSize) {
      const Triangle& lTriangle = aTrianglesNormals[aTrianglesNormalsBBoxCounter];
      aBBox3D += lTriangle.pt[0];
      aBBox3D += lTriangle.pt[1];
      aBBox3D += lTriangle.pt[2];
      ++aTrianglesNormalsBBoxCounter;
    }
  }

  {
	const SizeType lSize = aTrianglesNormalsColored.size();

    while (aTrianglesNormalsColoredBBoxCounter < lSize) {
      const Triangle& lTriangle = aTrianglesNormalsColored[aTrianglesNormalsColoredBBoxCounter];
      aBBox3D += lTriangle.pt[0];
      aBBox3D += lTriangle.pt[1];
      aBBox3D += lTriangle.pt[2];
      ++aTrianglesNormalsColoredBBoxCounter;
    }

	{
		const SizeType lSzie = aTriElemsNormalColor.size();

		while(aTriElemsNormalColorBBoxCounter < lSize) {
			const TriEelement& lTraingle = aTriElemsNormalColor[aTriElemsNormalColorBBoxCounter];
			aBBox3D += lTraingle.pt[0];
			aBBox3D += lTraingle.pt[1];
			aBBox3D += lTraingle.pt[2];
			++aTriElemsNormalColorBBoxCounter;
		}
	}

  }

  return aBBox3D;
}

void PrimitiveAccumulator::Render(const RenderParams& pParams)
{
  switch (pParams.eRMode)
  {
  case RenderParams::eFull:
    renderFull(pParams);
    break;
  case RenderParams::eBoundingBox:
    aBBox3D.render();
    break;
  case RenderParams::eFast:
    renderSimplified(pParams);
    break;
  default:
    FV_ASSERT(false);
  }
}



// Renders the content of aLines
void PrimitiveAccumulator::renderFacetsFrame(const RenderParams& pParams)
{
  const SizeType lNbFacets = aQuads.size()               +
                             aQuadsColored.size()        +
                             aQuadsNormals.size()        +
                             aQuadsNormalsColored.size() +
                             aTriangles.size()           +
                             aTrianglesColored.size()    +
                             aTrianglesNormals.size()    +
                             aTrianglesNormalsColored.size();

  // If facet drawing is enabled and we've got something to draw
  if(pParams.bFacetFrame && lNbFacets > 0) {

    // We push the attributes on the openGL attribute stack; so
    // we dont disturb the current values of color (current) or line width (line)
    glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_LIGHTING_BIT);
    // Warning: performance hit; when a lot of push(GL_LIGHTING_BIT) calls are made
    glDisable(GL_LIGHTING);
    glLineWidth(1.0f);
	glColor3f(pParams.cBkgColor.R, pParams.cBkgColor.G, pParams.cBkgColor.B);
	

    {
      Quads::const_iterator       lIterQuads    = aQuads.begin();
      const Quads::const_iterator lIterQuadsEnd = aQuads.end  ();

      while (lIterQuads != lIterQuadsEnd) {
        lIterQuads->renderFacetsFrame();
        ++lIterQuads;
      }
    }

    {
      QuadsColored::const_iterator       lIterQuads    = aQuadsColored.begin();
      const QuadsColored::const_iterator lIterQuadsEnd = aQuadsColored.end  ();

      while (lIterQuads != lIterQuadsEnd) {
        lIterQuads->renderFacetsFrame();
        ++lIterQuads;
      }
    }

    {
      QuadsNormals::const_iterator       lIterQuads    = aQuadsNormals.begin();
      const QuadsNormals::const_iterator lIterQuadsEnd = aQuadsNormals.end  ();

      while (lIterQuads != lIterQuadsEnd) {
        lIterQuads->renderFacetsFrame();
        ++lIterQuads;
      }
    }

    {
      QuadsNormalsColored::const_iterator       lIterQuads    = aQuadsNormalsColored.begin();
      const QuadsNormalsColored::const_iterator lIterQuadsEnd = aQuadsNormalsColored.end  ();

      while (lIterQuads != lIterQuadsEnd) {
        lIterQuads->renderFacetsFrame();
        ++lIterQuads;
      }
    }

    {
      Triangles::const_iterator       lIterTriangles    = aTriangles.begin();
      const Triangles::const_iterator lIterTrianglesEnd = aTriangles.end  ();

      while (lIterTriangles != lIterTrianglesEnd) {
        lIterTriangles->renderFacetsFrame();
        ++lIterTriangles;
      }
    }

    {
      TrianglesColored::const_iterator       lIterTriangles    = aTrianglesColored.begin();
      const TrianglesColored::const_iterator lIterTrianglesEnd = aTrianglesColored.end  ();

      while (lIterTriangles != lIterTrianglesEnd) {
        lIterTriangles->renderFacetsFrame();
        ++lIterTriangles;
      }
    }

    {
      TrianglesNormals::const_iterator       lIterTriangles    = aTrianglesNormals.begin();
      const TrianglesNormals::const_iterator lIterTrianglesEnd = aTrianglesNormals.end  ();

      while (lIterTriangles != lIterTrianglesEnd) {
        lIterTriangles->renderFacetsFrame();
        ++lIterTriangles;
      }
    }

    {
      TrianglesNormalsColored::const_iterator       lIterTriangles    = aTrianglesNormalsColored.begin();
      const TrianglesNormalsColored::const_iterator lIterTrianglesEnd = aTrianglesNormalsColored.end  ();

      while (lIterTriangles != lIterTrianglesEnd) {
        lIterTriangles->renderFacetsFrame();
        ++lIterTriangles;
      }
    }

    // Revert the lighting state and the line state
    glPopAttrib();
  }

  //glColor3f(0.7f, 0.7f, 0.7f);
}

void PrimitiveAccumulator::renderFull(const RenderParams& pParams)
{
	//std::cout<<"na poczatku render primitiw akumul\n";
  iPrimitiveOptimizerValue = pParams.iPrimitiveOptimizerValue;
  renderFacetsFrame            (pParams);
  renderLines                  ();
  renderPoints                 ();
  renderQuads                  ();
  renderQuadsNormals           ();
  renderTriangles              ();
  renderTrianglesNormals       ();
  renderLinesColored           ();
  renderPointsColored          ();
  renderQuadsColored           ();
  renderQuadsNormalsColored    ();
  renderTrianglesColored       ();
  renderTrianglesNormalsColored();
  renderTriElementNormalColor  ();
}

// Renders the content of aLines
void PrimitiveAccumulator::renderLines()
{
  int lPrimitiveCount = 0;
  if (aLines.size() > 0) {

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
	//glBlendFunc(GL_SRC0_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_BLEND);
    glBegin(GL_LINES);

    Lines::const_iterator       lIterLines    = aLines.begin();
    const Lines::const_iterator lIterLinesEnd = aLines.end  ();

    while (lIterLines != lIterLinesEnd) {

      const Line&     lLine = *lIterLines;
      const Vec3D& lP1   = lLine.pt[0];
      const Vec3D& lP2   = lLine.pt[1];

      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_LINES);
      }
      ++lIterLines;
    }

    glEnd();

    // Revert the lighting state
    glPopAttrib();

  }
}

// Renders the content of aLinesColored
void PrimitiveAccumulator::renderLinesColored()
{
	//return;
	//std::cout << "prim:linecolor\n";
  int lPrimitiveCount = 0;
  if (aLinesColored.size() > 0) {

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);

    LinesColored::const_iterator       lIterLinesColored    = aLinesColored.begin();
    const LinesColored::const_iterator lIterLinesColoredEnd = aLinesColored.end  ();
	//FILE* fp = fopen("siatka.txt","w");
    while (lIterLinesColored != lIterLinesColoredEnd) {

      const LineColored& lLineColored = *lIterLinesColored;
      const Vec3D&    lP1          = lLineColored.pt[0];
      const Vec3D&    lP2          = lLineColored.pt[1];
      const ColorRGB& lC1          = lLineColored.col1;
      const ColorRGB& lC2          = lLineColored.col2;
	  //try {
      glColor3f (lC1.R, lC1.G, lC1.B);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glColor3f (lC2.R, lC2.G, lC2.B);
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
	  //} catch(...)
	  //{
		  //std::cout<<"error in render linecolored\n";
		  //continue;
	  //}

	 // fprintf(fp,"%f %f %f %f %f %f\n",lP1._x(),lP1._y(),lP1._z(),
	//	  lP2._x(),lP2._y(),lP2._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_LINES);
      }
      ++lIterLinesColored;
    }

    glEnd();

	//fclose(fp);

    // Revert the lighting state
    glPopAttrib();

  }
}

// Renders the content of ptoints
void PrimitiveAccumulator::renderPoints()
{
  if (aPoints.size() > 0) {

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);

    Points::const_iterator       lIterPoints    = aPoints.begin();
    const Points::const_iterator lIterPointsEnd = aPoints.end  ();

    while (lIterPoints != lIterPointsEnd) {

      const Point&    lPoint = *lIterPoints;
      const Vec3D& lP     = lPoint.pt;

      glVertex3f(lP._x(), lP._y(), lP._z());

      ++lIterPoints;
    }

    glEnd();

    // Revert the lighting state
    glPopAttrib();

  }
}

// Renders the content of ptointsColored
void PrimitiveAccumulator::renderPointsColored()
{
  if (aPointsColored.size() > 0) {

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);

    PointsColored::const_iterator       lIterPointsColored    = aPointsColored.begin();
    const PointsColored::const_iterator lIterPointsColoredEnd = aPointsColored.end  ();

    while (lIterPointsColored != lIterPointsColoredEnd) {

      const PointColored&  lPointColored = *lIterPointsColored;
      const Vec3D&      lP            = lPointColored.pt;
      const ColorRGB&      lC            = lPointColored.col;

      glColor3f (lC.R, lC.G, lC.B);
      glVertex3f(lP._x(), lP._y(), lP._z());

      ++lIterPointsColored;
    }

    glEnd();

    // Revert the lighting state
    glPopAttrib();

  }
}

// Renders the content of aQuads
void PrimitiveAccumulator::renderQuads()
{
  int lPrimitiveCount = 0;

  if (aQuads.size() > 0) {

    glBegin(GL_QUADS);

    Quads::const_iterator       lIterQuads    = aQuads.begin();
    const Quads::const_iterator lIterQuadsEnd = aQuads.end  ();

    while (lIterQuads != lIterQuadsEnd) {

      const Quad&     lQuad = *lIterQuads;
      const Vec3D& lP1   = lQuad.pt[0];
      const Vec3D& lP2   = lQuad.pt[1];
      const Vec3D& lP3   = lQuad.pt[2];
      const Vec3D& lP4   = lQuad.pt[3];

      // Compute normals
      Vec3D lN1 = (lP2-lP1)%(lP4-lP1);
      Vec3D lN2 = (lP3-lP2)%(lP1-lP2);
      Vec3D lN3 = (lP4-lP3)%(lP2-lP3);
      Vec3D lN4 = (lP1-lP4)%(lP3-lP4);

      //try{
		  //if(lN1.Norm() == 0) {
			 // lN1.out("P1");
				//  exit(0);
		  //}
      lN1.normalize();
	  //if(lN2.Norm() == 0) {
			//  lN2.out("P2");
			//	  exit(0);
		 // }
      lN2.normalize();
	  //if(lN3.Norm() == 0) {
			//  lN3.out("P3");
			//	  exit(0);
		 // }
      lN3.normalize();
	  //if(lN4.Norm() == 0) {
			//  lN4.out("P4");
			//	  exit(0);
		 // }
      lN4.normalize();
	  //} catch(std::exception& ex) {
		 // std::cerr << ex.what();
		 // FV_ASSERT(false);
	  //}

      glNormal3f(lN1._x(), lN1._y(), lN1._z());
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f(lN2._x(), lN2._y(), lN2._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f(lN3._x(), lN3._y(), lN3._z());
      glVertex3f(lP3._x(), lP3._y(), lP3._z());
      glNormal3f(lN4._x(), lN4._y(), lN4._z());
      glVertex3f(lP4._x(), lP4._y(), lP4._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_QUADS);
      }
      ++lIterQuads;
    }

    glEnd();

  }
}

// Renders the content of aQuadsColored
void PrimitiveAccumulator::renderQuadsColored()
{
  int lPrimitiveCount = 0;
  if (aQuadsColored.size() > 0) {

    glBegin(GL_QUADS);

    QuadsColored::const_iterator       lIterQuadsColored    = aQuadsColored.begin();
    const QuadsColored::const_iterator lIterQuadsColoredEnd = aQuadsColored.end  ();

    while (lIterQuadsColored != lIterQuadsColoredEnd) {

      const QuadColored&  lQuadColored = *lIterQuadsColored;
      const Vec3D&     lP1          = lQuadColored.pt[0];
      const Vec3D&     lP2          = lQuadColored.pt[1];
      const Vec3D&     lP3          = lQuadColored.pt[2];
      const Vec3D&     lP4          = lQuadColored.pt[3];
      const ColorRGB&     lC1          = lQuadColored.col1;
      const ColorRGB&     lC2          = lQuadColored.col2;
      const ColorRGB&     lC3          = lQuadColored.col3;
      const ColorRGB&     lC4          = lQuadColored.col4;
	 
      // Compute normals
	 //(lP1-lP4).out("Roz1");
	// printf("\n");
	  //(lP3-lP4).out("Roz2");
	  //printf("\n");
	  FV_ASSERT(!Compare(lP1,lP2,lP3,lP4));
      Vec3D lN1 = (lP2-lP1) % (lP4-lP1);
      Vec3D lN2 = (lP3-lP2) % (lP1-lP2);
      Vec3D lN3 = (lP4-lP3) % (lP2-lP3);
      Vec3D lN4 = (lP1-lP4) % (lP3-lP4);



	 // FV_ASSERT(!Compare(lN1,lN2,lN3,lN4));
	  //lN1.out("N1");
	  //getchar();
	  
		 // if(lN1.Norm() == 0){
			//printf("stop\n");
			//lN1.out("lN1");
			//exit(0);
		 // }
      lN1.normalize();
	 // FV_ASSERT(lN2.Norm() != 0);
	  //if(lN2.Norm() == 0){
			//printf("stop\n");
			//lN2.out("lN2");
			//exit(0);
		 // }
      lN2.normalize();
	  //if(lN3.Norm() == 0){
			//printf("stop\n");
			//lN3.out("lN3");
			//exit(0);
		 // }
	 // FV_ASSERT(lN3.Norm() != 0);
      lN3.normalize();
	  //if(lN4.Norm() == 0){
			//printf("stop\n");
			//lN4.out("lN4");
			//exit(0);
		 // }
	 // FV_ASSERT(lN4.Norm() != 0);
      lN4.normalize();
	 // } catch(std::exception& ex) {
	//	  std::cerr << ex.what();
	//	  FV_ASSERT(false);
	 // }

      glNormal3f(lN1._x(), lN1._y(), lN1._z());
      glColor3f (lC1.R, lC1.G, lC1.B);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f(lN2._x(), lN2._y(), lN2._z());
      glColor3f (lC2.R, lC2.G, lC2.B);
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f(lN3._x(), lN3._y(), lN3._z());
      glColor3f (lC3.R, lC3.G, lC3.B);;
      glVertex3f(lP3._x(), lP3._y(), lP3._z());
      glNormal3f(lN4._x(), lN4._y(), lN4._z());
      glColor3f (lC4.R, lC4.G, lC4.B);
      glVertex3f(lP4._x(), lP4._y(), lP4._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_QUADS);
      }
      ++lIterQuadsColored;
    }

    glEnd();

  }
}

// Renders the content of aQuadsNormals
void PrimitiveAccumulator::renderQuadsNormals()
{
  int lPrimitiveCount = 0;
  if (aQuadsNormals.size() > 0) {

    glBegin(GL_QUADS);

    QuadsNormals::const_iterator       lIterQuadsNormals    = aQuadsNormals.begin();
    const QuadsNormals::const_iterator lIterQuadsNormalsEnd = aQuadsNormals.end  ();

    while (lIterQuadsNormals != lIterQuadsNormalsEnd) {

      const QuadNormals& lQuad = *lIterQuadsNormals;

      const Vec3D& lP1 = lQuad.pt[0];
      const Vec3D& lP2 = lQuad.pt[1];
      const Vec3D& lP3 = lQuad.pt[2];
      const Vec3D& lP4 = lQuad.pt[3];
      const Vec3D& lN1 = lQuad.n1;
      const Vec3D& lN2 = lQuad.n2;
      const Vec3D& lN3 = lQuad.n3;
      const Vec3D& lN4 = lQuad.n4;

      glNormal3f(lN1._x(), lN1._y(), lN1._z());
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f(lN2._x(), lN2._y(), lN2._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f(lN3._x(), lN3._y(), lN3._z());
      glVertex3f(lP3._x(), lP3._y(), lP3._z());
      glNormal3f(lN4._x(), lN4._y(), lN4._z());
      glVertex3f(lP4._x(), lP4._y(), lP4._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_QUADS);
      }
      ++lIterQuadsNormals;
    }

    glEnd();

  }
}

// Renders the content of aQuadsNormalsColored
void PrimitiveAccumulator::renderQuadsNormalsColored()
{
  int lPrimitiveCount = 0;
  if (aQuadsNormalsColored.size() > 0) {

    glBegin(GL_QUADS);

    QuadsNormalsColored::const_iterator       lIterQuadsNormalsColored    = aQuadsNormalsColored.begin();
    const QuadsNormalsColored::const_iterator lIterQuadsNormalsColoredEnd = aQuadsNormalsColored.end  ();

    while (lIterQuadsNormalsColored != lIterQuadsNormalsColoredEnd) {

      const QuadNormalsColored& lQuadColored = *lIterQuadsNormalsColored;

      const Vec3D& lP1 = lQuadColored.pt[0];
      const Vec3D& lP2 = lQuadColored.pt[1];
      const Vec3D& lP3 = lQuadColored.pt[2];
      const Vec3D& lP4 = lQuadColored.pt[3];
      const ColorRGB& lC1 = lQuadColored.col1;
      const ColorRGB& lC2 = lQuadColored.col2;
      const ColorRGB& lC3 = lQuadColored.col3;
      const ColorRGB& lC4 = lQuadColored.col4;
      const Vec3D& lN1 = lQuadColored.n1;
      const Vec3D& lN2 = lQuadColored.n2;
      const Vec3D& lN3 = lQuadColored.n3;
      const Vec3D& lN4 = lQuadColored.n4;

      glNormal3f(lN1._x(), lN1._y(), lN1._z());
      glColor3f (lC1.R, lC1.G, lC1.B);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f(lN2._x(), lN2._y(), lN2._z());
      glColor3f (lC2.R, lC2.G, lC2.B);
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f(lN3._x(), lN3._y(), lN3._z());
      glColor3f (lC3.R, lC3.G, lC3.B);
      glVertex3f(lP3._x(), lP3._y(), lP3._z());
      glNormal3f(lN4._x(), lN4._y(), lN4._z());
      glColor3f (lC4.R, lC4.G, lC4.B);
      glVertex3f(lP4._x(), lP4._y(), lP4._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_QUADS);
      }
      ++lIterQuadsNormalsColored;
    }

    glEnd();

  }
}

void PrimitiveAccumulator::renderSimplified(const RenderParams& pParams)
{
  // Check if we have to update the Simplified model
  if (aSimplifiedDirty) {
    if(pParams.iRMode_Fast_Option == 1) {
     //constructSimplified_heuristics();
    }
    else {
      //constructSimplified_points();
    }
  }

  RenderParams lParams = pParams;
  lParams.eRMode      = RenderParams::eFull;

  // We push the attributes on the openGL attribute stack; so
  // we dont disturb the current values of color (current) or line width (line)
  glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_LIGHTING_BIT);
  // Warning: performance hit; when a lot of push(GL_LIGHTING_BIT) calls are made
  glDisable(GL_LIGHTING);
  glPointSize(3);

  glColor3f(0.7, 0.5, 0.1);

  FV_ASSERT(aSimplified != 0);
  aSimplified->Render(lParams);

  // Revert the lighting state and the line state
  glPopAttrib();

  FV_ASSERT(aSimplified != 0);
  aSimplified->renderFull(lParams);
}

// Renders the content of aTriangles
void PrimitiveAccumulator::renderTriangles()
{
  int lPrimitiveCount = 0;

  if (aTriangles.size() > 0) {

    glBegin(GL_TRIANGLES);

    Triangles::const_iterator       lIterTriangles    = aTriangles.begin();
    const Triangles::const_iterator lIterTrianglesEnd = aTriangles.end  ();

    while (lIterTriangles != lIterTrianglesEnd) {

      const Triangle& lTriangle = *lIterTriangles;
      const Vec3D& lP1       = lTriangle.pt[0];
      const Vec3D& lP2       = lTriangle.pt[1];
      const Vec3D& lP3       = lTriangle.pt[2];

      // Compute normals
      Vec3D lN = (lP2-lP1) % (lP3-lP1);
      lN.normalize();

      glNormal3f( lN._x(),  lN._y(),  lN._z());
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f( lN._x(),  lN._y(),  lN._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f( lN._x(),  lN._y(),  lN._z());
      glVertex3f(lP3._x(), lP3._y(), lP3._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_TRIANGLES);
      }
      ++lIterTriangles;
    }

    glEnd();

  }
}

// Renders the content of aTrianglesColored
void PrimitiveAccumulator::renderTrianglesColored()
{
  int lPrimitiveCount = 0;

  if (aTrianglesColored.size() > 0) {
	  int i = 0;
	  
    glBegin(GL_TRIANGLES);

    TrianglesColored::const_iterator       lIterTrianglesColored    = aTrianglesColored.begin();
    const TrianglesColored::const_iterator lIterTrianglesColoredEnd = aTrianglesColored.end  ();

    while (lIterTrianglesColored != lIterTrianglesColoredEnd) {

      const TriangleColored& lTriangleColored = *lIterTrianglesColored;
      const Vec3D&        lP1          = lTriangleColored.pt[0];
      const Vec3D&        lP2          = lTriangleColored.pt[1];
      const Vec3D&        lP3          = lTriangleColored.pt[2];
      const ColorRGB&        lC1        = lTriangleColored.col1;
      const ColorRGB&        lC2          = lTriangleColored.col2;
      const ColorRGB&        lC3          = lTriangleColored.col3;

      // Compute normals
      Vec3D lN = (lP2-lP1) % (lP3-lP1);
	  try{
      lN.normalize();
	  } catch(std::runtime_error& er)
	  {
		  std::cout<< er.what() << "for i= "  << i << std::endl;
		  goto skok;
		 // exit(-1);
	  }

      glNormal3f( lN._x(),  lN._y(),  lN._z());
      glColor3f (lC1.R, lC1.G, lC1.B);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glColor3f (lC2.R, lC2.G, lC2.B);
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glColor3f (lC3.R, lC3.G, lC3.B);
      glVertex3f(lP3._x(), lP3._y(), lP3._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_TRIANGLES);
      }
skok:
      ++lIterTrianglesColored;
	  ++i;

	  
    }

    glEnd();
	  

  }
}

// Renders the content of aTrianglesNormals
void PrimitiveAccumulator::renderTrianglesNormals()
{
  int lPrimitiveCount = 0;

  if (aTrianglesNormals.size() > 0) {

    glBegin(GL_TRIANGLES);

    TrianglesNormals::const_iterator       lIterTrianglesNormals    = aTrianglesNormals.begin();
    const TrianglesNormals::const_iterator lIterTrianglesNormalsEnd = aTrianglesNormals.end  ();

    while (lIterTrianglesNormals != lIterTrianglesNormalsEnd) {

      const TriangleNormals& lTriangle = *lIterTrianglesNormals;

      const Vec3D& lP1 = lTriangle.pt[0];
      const Vec3D& lP2 = lTriangle.pt[1];
      const Vec3D& lP3 = lTriangle.pt[2];
      const Vec3D& lN1 = lTriangle.n1;
      const Vec3D& lN2 = lTriangle.n2;
      const Vec3D& lN3 = lTriangle.n3;

      glNormal3f(lN1._x(), lN1._y(), lN1._z());
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f(lN2._x(), lN2._y(), lN2._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f(lN3._x(), lN3._y(), lN3._z());
      glVertex3f(lP3._x(), lP3._y(), lP3._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_TRIANGLES);
      }
      ++lIterTrianglesNormals;
    }

    glEnd();

  }
}

// Renders the content of aTrianglesNormalsColored
void PrimitiveAccumulator::renderTrianglesNormalsColored()
{
  int lPrimitiveCount = 0;

  if (aTrianglesNormalsColored.size() > 0) {

    glBegin(GL_TRIANGLES);

    TrianglesNormalsColored::const_iterator       lIterTrianglesNormalsColored    = aTrianglesNormalsColored.begin();
    const TrianglesNormalsColored::const_iterator lIterTrianglesNormalsColoredEnd = aTrianglesNormalsColored.end  ();

    while (lIterTrianglesNormalsColored != lIterTrianglesNormalsColoredEnd) {

      const TriangleNormalsColored& lTriangleColored = *lIterTrianglesNormalsColored;

      const Vec3D& lP1 = lTriangleColored.pt[0];
      const Vec3D& lP2 = lTriangleColored.pt[1];
      const Vec3D& lP3 = lTriangleColored.pt[2];
	  const Vec3D& lN1 = lTriangleColored.n1;
      const Vec3D& lN2 = lTriangleColored.n2;
      const Vec3D& lN3 = lTriangleColored.n3;
      const ColorRGB& lC1 = lTriangleColored.col1;
      const ColorRGB& lC2 = lTriangleColored.col2;
      const ColorRGB& lC3 = lTriangleColored.col3;
      

      glNormal3f(lN1._x(), lN1._y(), lN1._z());
      glColor3f (lC1.R, lC1.G, lC1.B);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glNormal3f(lN2._x(), lN2._y(), lN2._z());
      glColor3f (lC2.R, lC2.G, lC2.B);
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glNormal3f(lN3._x(), lN3._y(), lN3._z());
      glColor3f (lC3.R, lC3.G, lC3.B);
      glVertex3f(lP3._x(), lP3._y(), lP3._z());

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_TRIANGLES);
      }
      ++lIterTrianglesNormalsColored;
    }

    glEnd();

  }
}

void PrimitiveAccumulator::renderTriElementNormalColor() const
{
	int lPrimitiveCount = 0;
	
	if (aTriElemsNormalColor.size() > 0) {

		glBegin(GL_TRIANGLES);

		TriElemsNormalColor::const_iterator       lIterTriElNC    = aTriElemsNormalColor.begin();
		const TriElemsNormalColor::const_iterator lIterTriElNCEnd = aTriElemsNormalColor.end  ();

		while (lIterTriElNC != lIterTriElNCEnd) {

			const TriEelementNormalColor& lTriangRef = *lIterTriElNC;

			const Vec3f& lP1 = lTriangRef.pt[0];
			const Vec3f& lP2 = lTriangRef.pt[1];
			const Vec3f& lP3 = lTriangRef.pt[2];
			const Vec3f& lN1 = lTriangRef.nv;
			//const Vec3d& lN2 = lTriangRef.n2;
			//const Vec3d& lN3 = lTriangRef.n3;
			const ColorRGB& lC1 = lTriangRef.c1;
			const ColorRGB& lC2 = lTriangRef.c2;
			const ColorRGB& lC3 = lTriangRef.c3;

			glNormal3fv(lN1.v);
			glColor3f  (lC1.R, lC1.G, lC1.B);
			glVertex3fv(lP1.v);
			//glNormal3dv(lN2.v3);
			glColor3f  (lC2.R, lC2.G, lC2.B);
			glVertex3fv(lP2.v);
			//glNormal3dv(lN3.v3);
			glColor3f  (lC3.R, lC3.G, lC3.B);
			glVertex3fv(lP3.v);
			
			lPrimitiveCount++;
			// This value sets the maximum length of the number of primitives
			// between a glBegin and a glEnd.  Some openGL driver optimize too
			// much and are slow to execute onto long primitive lists
			if(lPrimitiveCount % iPrimitiveOptimizerValue == 0) {
				glEnd();
				glBegin(GL_TRIANGLES);
			}
			++lIterTriElNC;
		}

		glEnd();

	}
}


struct SimplificationPoint
{
	SimplificationPoint(const float pX,
					  const float pY,
					  const float pZ)
	: fBarycenterDistance (std::numeric_limits<float>::max()/*FLT_MAX*/),
	  vPoint              (pX, pY, pZ),
	  vPrimiticeBarycenter()
	{}

	float  fBarycenterDistance;
	Vec3D  vPoint;
	Vec3D  vPrimiticeBarycenter;

};


template <class Container>
inline void scanForClosest(Container&           pContainer,
                           SimplificationPoint* pIterSPBegin,
                           SimplificationPoint* pIterSPEnd)
{
  const typename Container::size_type lSize = pContainer.size();

  for (typename Container::size_type i=0; i<lSize; ++i) {

    const typename Container::value_type&  lPrimitive = pContainer[i];
    SimplificationPoint*                   lIterSP    = pIterSPBegin;

    while (lIterSP != pIterSPEnd) {
      SimplificationPoint& lSP         = *lIterSP;
      const Vec3D       lBarycenter = lPrimitive .getBarycenter();
      const float       lDistance   = lBarycenter.getDistance(lSP.vPoint);
      if (lDistance < lSP.fBarycenterDistance) {
        lSP.vPrimiticeBarycenter = lBarycenter;
        lSP.fBarycenterDistance  = lDistance;
      }
      ++lIterSP;
    }
  }
}

void PrimitiveAccumulator::AddArrowTip(const Vec3D& p1,
                                       const Vec3D& p2,
                                       const ColorRGB& col,
                                       const bool colored,
                                       const float tipprop,
                                       const int nrpolygons)
{
	FV_ASSERT(tipprop  > 0.0f);
	FV_ASSERT(tipprop <= 1.0f);
	FV_ASSERT(nrpolygons >= 1);

	// Now construct the tip
	const float theta = 0.2f;
	const Vec3D vec = p2 - p1;
	const float length = vec.Norm();

	if (length > 0.0f) {
		const Vec3D direction = vec / length;

    // find a perpendicular vector to lDirection
    // using Gram-Schmidt on one of i,j or k. Pick
    // the one with the dot product with lDirection
    // closest to zero.
    Vec3D    lPerpendicularUnitary;
    const float absx = fabs(direction._x());
    const float absy = fabs(direction._y());
    const float absz = fabs(direction._z());
    if(absx < absy && absx < absz) {
		lPerpendicularUnitary = Vec3D(1.0f, 0.0f, 0.0f) - direction._x()*direction;
    } else if(absy < absz) {
		FV_ASSERT(absy <= absx);
		lPerpendicularUnitary = Vec3D(0.0f, 1.0f, 0.0f) - direction._y()*direction;
    } else {
		FV_ASSERT(absz <= absx);
		FV_ASSERT(absz <= absx);
		lPerpendicularUnitary = Vec3D(0.0f, 0.0f, 1.0f) - direction._z()*direction;
    }
    lPerpendicularUnitary.normalize();

    const float lTipRadius = theta * tipprop * length;
    Vec3D lPerpendicular = lTipRadius*lPerpendicularUnitary;

    // We rotate that lFirstTipPoint around lDirection
    // as many times as we want polygons
    Matrix<float> lTransformation;
    lTransformation.rotateAbout(direction.coord(), 2.0*M_PI/static_cast<float>(nrpolygons));

    const Vec3D lBaseTip    = (1.0f-tipprop)*vec + p1;
    const Vec3D lBaseNormal = -0.1f*direction;

    for(int i=0; i<=nrpolygons; ++i) {

		const Vec3D lNewPerpendicularUnitary = lTransformation*lPerpendicularUnitary;
		const Vec3D lNewPerpendicular        = lTipRadius*lPerpendicularUnitary;

		TriangleNormalsColored lTriangle;

		lTriangle.pt[0] = lBaseTip + lPerpendicular;
		lTriangle.pt[1] = lBaseTip + lNewPerpendicular;
		lTriangle.pt[2] = p2;
		lTriangle.col1 = col;
		lTriangle.col2 = col;
		lTriangle.col3 = col;
		lTriangle.n1 = lPerpendicularUnitary;
		lTriangle.n2 = lNewPerpendicularUnitary;
		lTriangle.n3 = direction;

      if(colored) {
		  aTrianglesNormalsColored.push_back(lTriangle);
      } else {
		  aTrianglesNormals.push_back(lTriangle);
      }

      lTriangle.pt[0] = lBaseTip + lPerpendicular;
      lTriangle.pt[1] = lBaseTip + lNewPerpendicular;
      lTriangle.pt[2] = lBaseTip;
      // Use normals to attenuate the intensity
      // of the base of the tip
      lTriangle.n1 = lBaseNormal;
      lTriangle.n2 = lBaseNormal;
      lTriangle.n3 = lBaseNormal;

      if(colored) {
		  aTrianglesNormalsColored.push_back(lTriangle);
      } else {
		  aTrianglesNormals.push_back(lTriangle);
      }

      lPerpendicularUnitary = lNewPerpendicularUnitary;
      lPerpendicular        = lNewPerpendicular;
    }
  }
}

//void PrimitiveAccumulator::constructSimplified_points()
//{
//
//  // The simplified model is a set of points located
//  // at the barycenter of a subset of the primitives.
//  // These primitives are chosen so that the
//  // object shpte can still be recognise by the user
//  // with a modest number of points.
//  // More precisely, they are the primitives "closest"
//  // to a set of points on a the circumscribed sphere to
//  // the BoundingBox.
//
//  if (!aSimplifiedDirty) {
//    // Just to be safe, should not pass here
//    // if nobody calls it when aSimplified is
//    // not used
//  }
//  else {
//
//    // Start from scratch
//    FV_ASSERT(aSimplified != 0);
//    if (aSimplified != this) {
//      delete aSimplified;
//    }
//    aSimplified = new PrimitiveAccumulator(false);
//
//    // We call getBoundingBox so that it gets updated
//    // if it needs to
//    const BBox3D& lBoundingBox = GetBBox3D();
//    const Vec3D     lCenter      = lBoundingBox.getCenter();
//    const float     lRadius      = lBoundingBox.getCircumscribedSphereRadius();
//
//    SimplificationPoint unitSphereSampling[] = {
//      SimplificationPoint(-1.0,0.0,0.0),
//      SimplificationPoint(-0.894427191,-0.4472135955,0.0),
//      SimplificationPoint(-0.894427191,0.0,-0.4472135955),
//      SimplificationPoint(-0.894427191,0.0,0.4472135955),
//      SimplificationPoint(-0.894427191,0.4472135955,0.0),
//      SimplificationPoint(-0.816496580928,-0.408248290464,-0.408248290464),
//      SimplificationPoint(-0.816496580928,-0.408248290464,0.408248290464),
//      SimplificationPoint(-0.816496580928,0.408248290464,-0.408248290464),
//      SimplificationPoint(-0.816496580928,0.408248290464,0.408248290464),
//      SimplificationPoint(-0.707106781187,-0.707106781187,0.0),
//      SimplificationPoint(-0.707106781187,0.0,-0.707106781187),
//      SimplificationPoint(-0.707106781187,0.0,0.707106781187),
//      SimplificationPoint(-0.707106781187,0.707106781187,0.0),
//      SimplificationPoint(-0.666666666667,-0.333333333333,-0.666666666667),
//      SimplificationPoint(-0.666666666667,-0.333333333333,0.666666666667),
//      SimplificationPoint(-0.666666666667,-0.666666666667,-0.333333333333),
//      SimplificationPoint(-0.666666666667,-0.666666666667,0.333333333333),
//      SimplificationPoint(-0.666666666667,0.333333333333,-0.666666666667),
//      SimplificationPoint(-0.666666666667,0.333333333333,0.666666666667),
//      SimplificationPoint(-0.666666666667,0.666666666667,-0.333333333333),
//      SimplificationPoint(-0.666666666667,0.666666666667,0.333333333333),
//      SimplificationPoint(-0.57735026919,-0.57735026919,-0.57735026919),
//      SimplificationPoint(-0.57735026919,-0.57735026919,0.57735026919),
//      SimplificationPoint(-0.57735026919,0.57735026919,-0.57735026919),
//      SimplificationPoint(-0.57735026919,0.57735026919,0.57735026919),
//      SimplificationPoint(-0.4472135955,-0.894427191,0.0),
//      SimplificationPoint(-0.4472135955,0.0,-0.894427191),
//      SimplificationPoint(-0.4472135955,0.0,0.894427191),
//      SimplificationPoint(-0.4472135955,0.894427191,0.0),
//      SimplificationPoint(-0.408248290464,-0.408248290464,-0.816496580928),
//      SimplificationPoint(-0.408248290464,-0.408248290464,0.816496580928),
//      SimplificationPoint(-0.408248290464,-0.816496580928,-0.408248290464),
//      SimplificationPoint(-0.408248290464,-0.816496580928,0.408248290464),
//      SimplificationPoint(-0.408248290464,0.408248290464,-0.816496580928),
//      SimplificationPoint(-0.408248290464,0.408248290464,0.816496580928),
//      SimplificationPoint(-0.408248290464,0.816496580928,-0.408248290464),
//      SimplificationPoint(-0.408248290464,0.816496580928,0.408248290464),
//      SimplificationPoint(-0.333333333333,-0.666666666667,-0.666666666667),
//      SimplificationPoint(-0.333333333333,-0.666666666667,0.666666666667),
//      SimplificationPoint(-0.333333333333,0.666666666667,-0.666666666667),
//      SimplificationPoint(-0.333333333333,0.666666666667,0.666666666667),
//      SimplificationPoint(0.0,-0.4472135955,-0.894427191),
//      SimplificationPoint(0.0,-0.4472135955,0.894427191),
//      SimplificationPoint(0.0,-0.707106781187,-0.707106781187),
//      SimplificationPoint(0.0,-0.707106781187,0.707106781187),
//      SimplificationPoint(0.0,-0.894427191,-0.4472135955),
//      SimplificationPoint(0.0,-0.894427191,0.4472135955),
//      SimplificationPoint(0.0,-1.0,0.0),
//      SimplificationPoint(0.0,0.0,-1.0),
//      SimplificationPoint(0.0,0.0,1.0),
//      SimplificationPoint(0.0,0.4472135955,-0.894427191),
//      SimplificationPoint(0.0,0.4472135955,0.894427191),
//      SimplificationPoint(0.0,0.707106781187,-0.707106781187),
//      SimplificationPoint(0.0,0.707106781187,0.707106781187),
//      SimplificationPoint(0.0,0.894427191,-0.4472135955),
//      SimplificationPoint(0.0,0.894427191,0.4472135955),
//      SimplificationPoint(0.0,1.0,0.0),
//      SimplificationPoint(0.333333333333,-0.666666666667,-0.666666666667),
//      SimplificationPoint(0.333333333333,-0.666666666667,0.666666666667),
//      SimplificationPoint(0.333333333333,0.666666666667,-0.666666666667),
//      SimplificationPoint(0.333333333333,0.666666666667,0.666666666667),
//      SimplificationPoint(0.408248290464,-0.408248290464,-0.816496580928),
//      SimplificationPoint(0.408248290464,-0.408248290464,0.816496580928),
//      SimplificationPoint(0.408248290464,-0.816496580928,-0.408248290464),
//      SimplificationPoint(0.408248290464,-0.816496580928,0.408248290464),
//      SimplificationPoint(0.408248290464,0.408248290464,-0.816496580928),
//      SimplificationPoint(0.408248290464,0.408248290464,0.816496580928),
//      SimplificationPoint(0.408248290464,0.816496580928,-0.408248290464),
//      SimplificationPoint(0.408248290464,0.816496580928,0.408248290464),
//      SimplificationPoint(0.4472135955,-0.894427191,0.0),
//      SimplificationPoint(0.4472135955,0.0,-0.894427191),
//      SimplificationPoint(0.4472135955,0.0,0.894427191),
//      SimplificationPoint(0.4472135955,0.894427191,0.0),
//      SimplificationPoint(0.57735026919,-0.57735026919,-0.57735026919),
//      SimplificationPoint(0.57735026919,-0.57735026919,0.57735026919),
//      SimplificationPoint(0.57735026919,0.57735026919,-0.57735026919),
//      SimplificationPoint(0.57735026919,0.57735026919,0.57735026919),
//      SimplificationPoint(0.666666666667,-0.333333333333,-0.666666666667),
//      SimplificationPoint(0.666666666667,-0.333333333333,0.666666666667),
//      SimplificationPoint(0.666666666667,-0.666666666667,-0.333333333333),
//      SimplificationPoint(0.666666666667,-0.666666666667,0.333333333333),
//      SimplificationPoint(0.666666666667,0.333333333333,-0.666666666667),
//      SimplificationPoint(0.666666666667,0.333333333333,0.666666666667),
//      SimplificationPoint(0.666666666667,0.666666666667,-0.333333333333),
//      SimplificationPoint(0.666666666667,0.666666666667,0.333333333333),
//      SimplificationPoint(0.707106781187,-0.707106781187,0.0),
//      SimplificationPoint(0.707106781187,0.0,-0.707106781187),
//      SimplificationPoint(0.707106781187,0.0,0.707106781187),
//      SimplificationPoint(0.707106781187,0.707106781187,0.0),
//      SimplificationPoint(0.816496580928,-0.408248290464,-0.408248290464),
//      SimplificationPoint(0.816496580928,-0.408248290464,0.408248290464),
//      SimplificationPoint(0.816496580928,0.408248290464,-0.408248290464),
//      SimplificationPoint(0.816496580928,0.408248290464,0.408248290464),
//      SimplificationPoint(0.894427191,-0.4472135955,0.0),
//      SimplificationPoint(0.894427191,0.0,-0.4472135955),
//      SimplificationPoint(0.894427191,0.0,0.4472135955),
//      SimplificationPoint(0.894427191,0.4472135955,0.0),
//      SimplificationPoint(1.0,0.0,0.0)
//    };
//
////     SimplificationPoint unitSphereSampling[] = {
////       SimplificationPoint(-1.0,0.0,0.0),
////       SimplificationPoint(-0.707106781187,-0.707106781187,0.0),
////       SimplificationPoint(-0.707106781187,0.0,-0.707106781187),
////       SimplificationPoint(-0.707106781187,0.0,0.707106781187),
////       SimplificationPoint(-0.707106781187,0.707106781187,0.0),
////       SimplificationPoint(-0.57735026919,-0.57735026919,-0.57735026919),
////       SimplificationPoint(-0.57735026919,-0.57735026919,0.57735026919),
////       SimplificationPoint(-0.57735026919,0.57735026919,-0.57735026919),
////       SimplificationPoint(-0.57735026919,0.57735026919,0.57735026919),
////       SimplificationPoint(0.0,-0.707106781187,-0.707106781187),
////       SimplificationPoint(0.0,-0.707106781187,0.707106781187),
////       SimplificationPoint(0.0,-1.0,0.0),
////       SimplificationPoint(0.0,0.0,-1.0),
////       SimplificationPoint(0.0,0.0,1.0),
////       SimplificationPoint(0.0,0.707106781187,-0.707106781187),
////       SimplificationPoint(0.0,0.707106781187,0.707106781187),
////       SimplificationPoint(0.0,1.0,0.0),
////       SimplificationPoint(0.57735026919,-0.57735026919,-0.57735026919),
////       SimplificationPoint(0.57735026919,-0.57735026919,0.57735026919),
////       SimplificationPoint(0.57735026919,0.57735026919,-0.57735026919),
////       SimplificationPoint(0.57735026919,0.57735026919,0.57735026919),
////       SimplificationPoint(0.707106781187,-0.707106781187,0.0),
////       SimplificationPoint(0.707106781187,0.0,-0.707106781187),
////       SimplificationPoint(0.707106781187,0.0,0.707106781187),
////       SimplificationPoint(0.707106781187,0.707106781187,0.0),
////       SimplificationPoint(1.0,0.0,0.0)
////     };
//
//    SimplificationPoint* lBegin = unitSphereSampling;
//    SimplificationPoint* lIter  = lBegin;
//    SimplificationPoint* lEnd   = unitSphereSampling + (sizeof(unitSphereSampling)/sizeof(SimplificationPoint));
//
//    // Scale the point from the unit sphere to the circumscribed sphere
//    FV_ASSERT(lBegin != 0);
//    FV_ASSERT(lEnd   != 0);
//
//    while (lIter != lEnd) {
//      FV_ASSERT(lIter != 0);
//      lIter->vPoint *= lRadius;
//      lIter->vPoint += lCenter;
//      ++lIter;
//    }
//
//    // Scan all the primitives to find the closest ones
//    // For now, we choose only one primitive in all the types
//    // for each point on the shere. Again, we'll have to see
//    // if if wouldn't be better to have a separate simplication
//    // for each type of primitive.  We'll see.
//    // Scan only the new primitives for each type
//    scanForClosest(aLines,                   lBegin, lEnd);
//    scanForClosest(aLinesColored,            lBegin, lEnd);
//    scanForClosest(ptoints,                  lBegin, lEnd);
//    scanForClosest(ptointsColored,           lBegin, lEnd);
//    scanForClosest(aQuads,                   lBegin, lEnd);
//    scanForClosest(aQuadsColored,            lBegin, lEnd);
//    scanForClosest(aQuadsNormals,            lBegin, lEnd);
//    scanForClosest(aQuadsNormalsColored,     lBegin, lEnd);
//    scanForClosest(aTriangles,               lBegin, lEnd);
//    scanForClosest(aTrianglesColored,        lBegin, lEnd);
//    scanForClosest(aTrianglesNormals,        lBegin, lEnd);
//    scanForClosest(aTrianglesNormalsColored, lBegin, lEnd);
//
//    // Now add Points to aSimplified for each primitive barycenter chosen
//    lIter = lBegin;
//
//    while (lIter != lEnd) {
//      FV_ASSERT(lIter != 0);
//      FV_ASSERT(aSimplified != 0);
//      aSimplified->AddPoint(lIter->vPrimiticeBarycenter);
//      ++lIter;
//    }
//  }
//
//  // aSimplified is not dirty anymore
//  aSimplifiedDirty = false;
//}
//
//#define SIMPLIFICATION_FACTOR 43
//#define SIMPLIFICATION_RANDOM_FACTOR 13
//void PrimitiveAccumulator::constructSimplified_heuristics()
//{
//  if (!aSimplifiedDirty) {
//    // Just to be safe, should not pass here
//    // if nobody calls it when aSimplified is
//    // not used
//  }
//  else {
//
//    // Start from scratch
//    FV_ASSERT(aSimplified != 0);
//    if (aSimplified != this) {
//      delete aSimplified;
//    }
//    aSimplified = new PrimitiveAccumulator(false);
//
//    // FOR OTHER PRIMITIVES (LINES, POINTS)
//    // WE SHOW 1/x
//    //  the factor is SIMPLIFICATION_FACTOR instead of 100, because that on structured data with a modulo
//    //  of 100, the chosent points does not seem "random"
//    {
//      size_t lNbPoints = ptoints.size();
//      for(size_t i=0;i<lNbPoints;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->AddPoint(ptoints[i].pt);
//      }
//    }
//    {
//      size_t lNbLines = aLines.size();
//      for(size_t i=0;i<lNbLines;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified-AddLine(aLines[i].pt[0],aLines[i].pt[1]);
//      }
//    }
//    {
//      size_t lNbTri = aTriangles.size();
//      for(size_t i=0;i<lNbTri;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->AddTriangle(aTriangles[i].pt[0],aTriangles[i].pt[1],aTriangles[i].pt[2]);
//      }
//    }
//    {
//      size_t lNbTri = aTrianglesColored.size();
//      for(int i=0;i<lNbTri;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->AddTriangle(aTrianglesColored[i].pt[0],aTrianglesColored[i].pt[1],aTrianglesColored[i].pt[2]);
//      }
//    }
//    {
//      size_t lNbTri = aTrianglesNormalsColored.size();
//      for(size_t i=0;i<lNbTri;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->AddTriangle(aTrianglesNormalsColored[i].pt[0],aTrianglesNormalsColored[i].pt[1],aTrianglesNormalsColored[i].pt[2]);
//      }
//    }
//    {
//      size_t lNbTri = aTrianglesNormals.size();
//      for(size_t i=0;i<lNbTri;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->AddTriangle(aTrianglesNormals[i].pt[0],aTrianglesNormals[i].pt[1],aTrianglesNormals[i].pt[2]);
//      }
//    }
//    {
//      size_t lNbQuads = aQuads.size();
//      for(size_t i=0;i<lNbQuads;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->AddQuad(aQuads[i].pt[0],aQuads[i].pt[1],aQuads[i].pt[2],aQuads[i].pt[3]);
//      }
//    }
//    {
//      size_t lNbQuads = aQuadsColored.size();
//      for(size_t i=0;i<lNbQuads;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->addQuad(aQuadsColored[i].pt[0],aQuadsColored[i].pt[1],aQuadsColored[i].pt[2],aQuadsColored[i].pt[3]);
//      }
//    }
//    {
//      size_t lNbQuads = aQuadsNormalsColored.size();
//      for(size_t i=0;i<lNbQuads;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->addQuad(aQuadsNormalsColored[i].pt[0],aQuadsNormalsColored[i].pt[1],aQuadsNormalsColored[i].pt[2],aQuadsNormalsColored[i].pt[3]);
//      }
//    }
//    {
//      size_t lNbQuads = aQuadsColored.size();
//      for(size_t i=0;i<lNbQuads;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
//        aSimplified->addQuad(aQuadsColored[i].pt[0],aQuadsColored[i].pt[1],aQuadsColored[i].pt[2],aQuadsColored[i].pt[3]);
//      }
//    }
//  }
//
//  // aSimplified is not dirty anymore
//  aSimplifiedDirty = false;
//}

} // end namespace FemViewer
