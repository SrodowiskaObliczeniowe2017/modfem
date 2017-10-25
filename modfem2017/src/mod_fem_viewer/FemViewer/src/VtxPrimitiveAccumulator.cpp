#include "VtxPrimitiveAccumulator.h"
#include "VtxAccumulator.h"
#include "PrimitiveAccumulator.h"


#include <cmath>
#include <iostream>
#include <stdio.h>

namespace FemViewer {

VtxPrimitiveAccumulator::VtxPrimitiveAccumulator(VtxAccumulator& pVertexes)
  : aBBox3D                (),
    aColors                     (pVertexes.getColors()),
    aLines                      (),
    aLinesBBoxCounter           (0),
    aPoints                     (),
    aPointsBBoxCounter          (0),
    aQuads                      (),
    aQuadsBBoxCounter           (0),
    aNormals                    (pVertexes.getNormals()),
//    aTexCoords                  (pVertexes.getTexCoords()),
    aSimplified                 (),
    aSimplifiedDirty            (true),
    aSimplifiedSelf             (false),
    aTriangles                  (),
    aTrianglesBBoxCounter       (0),
    aVertices                   (pVertexes.getVertices()),
    aPrimitiveOptimizerValue    (100),
    aVtxAccumulator          (pVertexes)
{
  aSimplified = new PrimitiveAccumulator(false);
  FV_ASSERT(aSimplified != 0);
}

VtxPrimitiveAccumulator::~VtxPrimitiveAccumulator()
{
  FV_ASSERT(aSimplified != 0);
  delete aSimplified;
}

// Returns false if one of the parameters
// is out of range
bool VtxPrimitiveAccumulator::addLine(const int pP1,
                                           const int pP2)
{
  const int  lSize      = static_cast<int>(aVertices.size());
  const bool lIndicesOk = (pP1 >= 0 && pP1 < lSize &&
                           pP2 >= 0 && pP2 < lSize);

  if (lIndicesOk) {
    FV_ASSERT(aVtxAccumulator.getVerticesFrozen());

    Line lLine;
    lLine.aP[0] = pP1;
    lLine.aP[1] = pP2;

    aLines.push_back(lLine);

    aSimplifiedDirty = true;
  }

  return lIndicesOk;
}

// Returns false if the parameters is out of range
bool VtxPrimitiveAccumulator::addPoint(const int pP)
{

  const int  lSize      = static_cast<int>(aVertices.size());
  const bool lIndicesOk = (pP >= 0 && pP < lSize);

  if (lIndicesOk) {
    FV_ASSERT(aVtxAccumulator.getVerticesFrozen());

    Point lPoint;
    lPoint.aP = pP;

    aPoints.push_back(lPoint);

    aSimplifiedDirty = true;
  }

  return lIndicesOk;
}

// Returns false if one of the parameters
// is out of range
bool VtxPrimitiveAccumulator::addQuad(const int pP1,
                                           const int pP2,
                                           const int pP3,
                                           const int pP4)
{
  const int  lSize      = static_cast<int>(aVertices.size());
  const bool lIndicesOk = (pP1 >= 0 && pP1 < lSize &&
                           pP2 >= 0 && pP2 < lSize &&
                           pP3 >= 0 && pP3 < lSize &&
                           pP4 >= 0 && pP4 < lSize);

  if (lIndicesOk) {
    FV_ASSERT(aVtxAccumulator.getVerticesFrozen());

    Quad lQuad;
    lQuad.aP[0] = pP1;
    lQuad.aP[1] = pP2;
    lQuad.aP[2] = pP3;
    lQuad.aP[3] = pP4;

    aQuads.push_back(lQuad);

    aSimplifiedDirty = true;
  }

  return lIndicesOk;
}

// Returns false if one of the parameters
// is out of range
bool VtxPrimitiveAccumulator::addTriangle(const int pP1,
                                               const int pP2,
                                               const int pP3)
{

  const int  lSize      = static_cast<int>(aVertices.size());
  const bool lIndicesOk = (pP1 >= 0 && pP1 < lSize &&
                           pP2 >= 0 && pP2 < lSize &&
                           pP3 >= 0 && pP3 < lSize);

  if (lIndicesOk) {
    FV_ASSERT(aVtxAccumulator.getVerticesFrozen());

    Triangle lTriangle;
    lTriangle.aP[0] = pP1;
    lTriangle.aP[1] = pP2;
    lTriangle.aP[2] = pP3;

    aTriangles.push_back(lTriangle);
    aSimplifiedDirty = true;
  }

  return lIndicesOk;
}

#ifdef FV_DUMP_MEMORY_USAGE

template <class Container>
inline
std::string getStringSizeAndCapacity(const Container& pContainer)
{
  char lSizeAndCapacity[128];
  sprintf(lSizeAndCapacity,"%d/%d",
          sizeof(typename Container::value_type)*pContainer.size(),
          sizeof(typename Container::value_type)*pContainer.capacity());
  return std::string(lSizeAndCapacity);
}

#endif // #ifdef  GLV_DUMP_MEMORY_USAGE

// Dump in ASCII the caracteristics of the VtxPrimitiveAccumulator
void VtxPrimitiveAccumulator::dumpCharacteristics(std::ostream&       pOstream,
                                                       const std::string&  pIndentation,
                                                       const Matrix<float>&    pTransformation)
{
  pOstream << pIndentation << "VtxPrimitiveAccumulator " << std::endl;
  std::string lIndentation = pIndentation + "  ";

  // Dump the data of the shared vertexes
  aVtxAccumulator.dumpCharacteristics(pOstream,lIndentation,pTransformation);

  // Update the bounding box. But we use the attribute
  // directly just after the call because we need it non-const
  getBBox3D();
  aBBox3D.dumpCharacteristics(pOstream, lIndentation, pTransformation);

#ifdef GLV_DUMP_SIMPLICATION
  if (!aSimplifiedSelf) {
    pOstream << lIndentation << "Simplified "  << std::endl;
    aSimplified->dumpCharacteristics(pOstream, lIndentation + "  ", pTransformation);
  }
#endif // #ifdef GLV_DUMP_SIMPLICATION

#ifdef GLV_DUMP_MEMORY_USAGE
  {

    pOstream << lIndentation << "Memory used by the VtxPrimitiveAccumulator = " << sizeof(*this) << std::endl;

    pOstream << lIndentation << "Memory used by aPoints    = " << getStringSizeAndCapacity(aPoints   ) << std::endl;
    pOstream << lIndentation << "Memory used by aLines     = " << getStringSizeAndCapacity(aLines    ) << std::endl;
    pOstream << lIndentation << "Memory used by aTriangles = " << getStringSizeAndCapacity(aTriangles) << std::endl;
    pOstream << lIndentation << "Memory used by aQuads     = " << getStringSizeAndCapacity(aQuads    ) << std::endl;

  }
#endif // #ifdef GLV_DUMP_MEMORY_USAGE

  pOstream << lIndentation << "Number of point_v    = " << aPoints    .size() << std::endl;
  pOstream << lIndentation << "Number of line_v     = " << aLines     .size() << std::endl;
  pOstream << lIndentation << "Number of triangle_v = " << aTriangles .size() << std::endl;
  pOstream << lIndentation << "Number of quad_v     = " << aQuads     .size() << std::endl;
}

const BBox3D& VtxPrimitiveAccumulator::getBBox3D() const
{

  // We could take for granted that all the vertices are
  // used by at least a primitive. But since there is a
  // chance that one vertex could not be used, we the
  // primitives instead.

  FV_ASSERT(aLinesBBoxCounter     <= aLines.size());
  FV_ASSERT(aPointsBBoxCounter    <= aPoints.size());
  FV_ASSERT(aQuadsBBoxCounter     <= aQuads.size());
  FV_ASSERT(aTrianglesBBoxCounter <= aTriangles.size());

  // Scan only the new primitives for each type
  {
    const SizeType lSize = aLines.size();

    while (aLinesBBoxCounter < lSize) {
      const Line& lLine = aLines[aLinesBBoxCounter];
      FV_ASSERT(lLine.aP[0] >= 0 && lLine.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lLine.aP[1] >= 0 && lLine.aP[1] < static_cast<int>(aVertices.size()));
      aBBox3D += aVertices[lLine.aP[0]].position;
      aBBox3D += aVertices[lLine.aP[1]].position;
      ++aLinesBBoxCounter;
    }
  }
	
  {
    const SizeType lSize = aPoints.size();

    while (aPointsBBoxCounter < lSize) {
      const Point& lPoint = aPoints[aPointsBBoxCounter];
      FV_ASSERT(lPoint.aP >= 0 && lPoint.aP < static_cast<int>(aVertices.size()));
      aBBox3D += aVertices[lPoint.aP].position;
      ++aPointsBBoxCounter;
    }
  }

  {
    const SizeType lSize = aQuads.size();

    while (aQuadsBBoxCounter < lSize) {
      const Quad& lQuad = aQuads[aQuadsBBoxCounter];
      FV_ASSERT(lQuad.aP[0] >= 0 && lQuad.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[1] >= 0 && lQuad.aP[1] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[2] >= 0 && lQuad.aP[2] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[3] >= 0 && lQuad.aP[3] < static_cast<int>(aVertices.size()));
      aBBox3D += aVertices[lQuad.aP[0]].position;;
      aBBox3D += aVertices[lQuad.aP[1]].position;
      aBBox3D += aVertices[lQuad.aP[2]].position;
      aBBox3D += aVertices[lQuad.aP[3]].position;
      ++aQuadsBBoxCounter;
    }
  }

  {
    const SizeType lSize = aTriangles.size();

    while (aTrianglesBBoxCounter < lSize) {
      const Triangle& lTriangle = aTriangles[aTrianglesBBoxCounter];
      FV_ASSERT(lTriangle.aP[0] >= 0 && lTriangle.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lTriangle.aP[1] >= 0 && lTriangle.aP[1] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lTriangle.aP[2] >= 0 && lTriangle.aP[2] < static_cast<int>(aVertices.size()));
      aBBox3D += aVertices[lTriangle.aP[0]].position;
      aBBox3D += aVertices[lTriangle.aP[1]].position;
      aBBox3D += aVertices[lTriangle.aP[2]].position;
      ++aTrianglesBBoxCounter;
    }
  }

  return aBBox3D;
}

void VtxPrimitiveAccumulator::computeNormals()
{
  FV_ASSERT(aNormals.size() != aVertices.size() || aVertices.size() == 0);

  const Vec3D lNullVector     (0.0f, 0.0f, 0.0f);
  const Vec3D lArbitraryNormal(1.0f, 0.0f, 0.0f);

  aNormals.resize(aVertices.size());
  std::fill(aNormals.begin(), aNormals.end(), lNullVector);

  Quads::const_iterator       lIterQuads    = aQuads.begin();
  const Quads::const_iterator lIterQuadsEnd = aQuads.end  ();

  while (lIterQuads != lIterQuadsEnd) {

    const Quad& lQuad = *lIterQuads;
    FV_ASSERT(lQuad.aP[0] >= 0 && lQuad.aP[0] < static_cast<int>(aVertices.size()));
    FV_ASSERT(lQuad.aP[1] >= 0 && lQuad.aP[1] < static_cast<int>(aVertices.size()));
    FV_ASSERT(lQuad.aP[2] >= 0 && lQuad.aP[2] < static_cast<int>(aVertices.size()));
    FV_ASSERT(lQuad.aP[3] >= 0 && lQuad.aP[3] < static_cast<int>(aVertices.size()));
    const Vec3D& lP1 = aVertices[lQuad.aP[0]].position.v;
    const Vec3D& lP2 = aVertices[lQuad.aP[1]].position.v;
    const Vec3D& lP3 = aVertices[lQuad.aP[2]].position.v;
    const Vec3D& lP4 = aVertices[lQuad.aP[3]].position.v;

    // Compute normals
    Vec3D lN1 = (lP2-lP1) % (lP4-lP1);
    Vec3D lN2 = (lP3-lP2) % (lP1-lP2);
    Vec3D lN3 = (lP4-lP3) % (lP2-lP3);
    Vec3D lN4 = (lP1-lP4) % (lP3-lP4);

    // We do not normalize the normals at this points
    // This plays the role of weights giving more
    // importance to bigger polygons
    aNormals[lQuad.aP[0]] += lN1;
    aNormals[lQuad.aP[1]] += lN2;
    aNormals[lQuad.aP[2]] += lN3;
    aNormals[lQuad.aP[3]] += lN4;

    ++lIterQuads;
  }

  Triangles::const_iterator       lIterTriangles    = aTriangles.begin();
  const Triangles::const_iterator lIterTrianglesEnd = aTriangles.end  ();

  while (lIterTriangles != lIterTrianglesEnd) {

    const Triangle& lTriangle = *lIterTriangles;
    FV_ASSERT(lTriangle.aP[0] >= 0 && lTriangle.aP[0] < static_cast<int>(aVertices.size()));
    FV_ASSERT(lTriangle.aP[1] >= 0 && lTriangle.aP[1] < static_cast<int>(aVertices.size()));
    FV_ASSERT(lTriangle.aP[2] >= 0 && lTriangle.aP[2] < static_cast<int>(aVertices.size()));
    const Vec3D& lP1 = aVertices[lTriangle.aP[0]].position.v;
    const Vec3D& lP2 = aVertices[lTriangle.aP[1]].position.v;
    const Vec3D& lP3 = aVertices[lTriangle.aP[2]].position.v;

    // Compute normals
    Vec3D lN = (lP2-lP1) % (lP3-lP1);

    // We do not normalize the normals at this points
    // This plays the role of weights giving more
    // importance to bigger polygons
    aNormals[lTriangle.aP[0]] += lN;
    aNormals[lTriangle.aP[1]] += lN;
    aNormals[lTriangle.aP[2]] += lN;

    ++lIterTriangles;
  }

  Normals::iterator       lIterNormals    = aNormals.begin();
  const Normals::iterator lIterNormalsEnd = aNormals.end  ();

  while (lIterNormals != lIterNormalsEnd) {
    // There is the possibility that a vertex doesn't
    // have adjacent Quads or Triangles. So we have to
    // check that
    Vec3D& lNormal = *lIterNormals;
    if (!(lNormal == lNullVector)) {
      lNormal.normalize();
    }
    else {
      // But there is also the possibility that two
      // adjacent Quad and/or Triangle are facing in
      // opposite directions. This would cause problems
      // later on, so we put an arbitrary unity vector
      // for the normal
      lNormal = lArbitraryNormal;
    }
    ++lIterNormals;
  }
}

void VtxPrimitiveAccumulator::render(const RenderParams& pParams)
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
    // Should not happen!
    FV_ASSERT(false);
  }
}

// Renders the content of aLines
void VtxPrimitiveAccumulator::renderFacetsFrame(const RenderParams& pParams)
{
  const SizeType lNbFacets = aQuads.size() + aTriangles.size();

  // If facet drawing is enabled and we've got something to draw
  if(pParams.bFacetFrame && lNbFacets > 0) {

    // We push the attributes on the openGL attribute stack; so
    // we dont disturb the current values of color (current) or line width (line)
    glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_LIGHTING_BIT);
    // Warning: performance hit; when a lot of push(GL_LIGHTING_BIT) calls are made
    glDisable(GL_LIGHTING);
    glLineWidth(1);
    glColor3f(pParams.cBkgColor.R, pParams.cBkgColor.G, pParams.cBkgColor.B);

    Quads::const_iterator       lIterQuads    = aQuads.begin();
    const Quads::const_iterator lIterQuadsEnd = aQuads.end  ();

    while (lIterQuads != lIterQuadsEnd) {

      const Quad& lQuad = *lIterQuads;
      FV_ASSERT(lQuad.aP[0] >= 0 && lQuad.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[1] >= 0 && lQuad.aP[1] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[2] >= 0 && lQuad.aP[2] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[3] >= 0 && lQuad.aP[3] < static_cast<int>(aVertices.size()));
      const Vec3D& lP1 = aVertices[lQuad.aP[0]].position.v;
      const Vec3D& lP2 = aVertices[lQuad.aP[1]].position.v;
      const Vec3D& lP3 = aVertices[lQuad.aP[2]].position.v;
      const Vec3D& lP4 = aVertices[lQuad.aP[3]].position.v;

      glBegin(GL_LINE_LOOP);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glVertex3f(lP3._x(), lP3._y(), lP3._z());
      glVertex3f(lP4._x(), lP4._y(), lP4._z());
      glEnd();

      ++lIterQuads;
    }

    Triangles::const_iterator       lIterTriangles    = aTriangles.begin();
    const Triangles::const_iterator lIterTrianglesEnd = aTriangles.end  ();

    while (lIterTriangles != lIterTrianglesEnd) {

      const Triangle& lTriangle = *lIterTriangles;
      FV_ASSERT(lTriangle.aP[0] >= 0 && lTriangle.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lTriangle.aP[1] >= 0 && lTriangle.aP[1] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lTriangle.aP[2] >= 0 && lTriangle.aP[2] < static_cast<int>(aVertices.size()));
      const Vec3D& lP1 = aVertices[lTriangle.aP[0]].position.v;
      const Vec3D& lP2 = aVertices[lTriangle.aP[1]].position.v;
      const Vec3D& lP3 = aVertices[lTriangle.aP[2]].position.v;

      glBegin(GL_LINE_LOOP);
      glVertex3f(lP1._x(), lP1._y(), lP1._z());
      glVertex3f(lP2._x(), lP2._y(), lP2._z());
      glVertex3f(lP3._x(), lP3._y(), lP3._z());
      glEnd();

      ++lIterTriangles;
    }

    // Revert the lighting state and the line state
    glPopAttrib();
  }
}

void VtxPrimitiveAccumulator::renderFull(const RenderParams& pParams)
{
  aPrimitiveOptimizerValue = pParams.iPrimitiveOptimizerValue;
  if (pParams.iRMode_Fast_Option == 1) {
    if (aNormals.size() != aVertices.size()) {
      computeNormals();
    }
  }
  renderFacetsFrame(pParams);
  renderLines      ();
  renderPoints     ();
  renderQuads      ();
  renderTriangles  ();
}

// Renders the content of aLines
void VtxPrimitiveAccumulator::renderLines()
{
  int lPrimitiveCount = 0;

  if (aLines.size() > 0) {

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);

    Lines::const_iterator       lIterLines    = aLines.begin();
    const Lines::const_iterator lIterLinesEnd = aLines.end  ();

    while (lIterLines != lIterLinesEnd) {

      const Line& lLine = *lIterLines;
      FV_ASSERT(lLine.aP[0] >= 0 && lLine.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lLine.aP[1] >= 0 && lLine.aP[1] < static_cast<int>(aVertices.size()));
      const Vec3D& lP1 = aVertices[lLine.aP[0]].position.v;
      const Vec3D& lP2 = aVertices[lLine.aP[1]].position.v;

/*      if(aTexCoords.size() == aVertices.size()) {
        FV_ASSERT(aTexCoords.size() != 0);
        const Tex2D& lC1 = aTexCoords[lLine.aP[0]];
        const Tex2D& lC2 = aTexCoords[lLine.aP[1]];

        glTexCoord2f(lC1.u(), lC1.v());
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
        glTexCoord2f(lC2.u(), lC2.v());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
      }
      else*/ if (aColors.size() == aVertices.size()) {
        FV_ASSERT(aColors.size() != 0);
        const Vec3D& lC1 = aColors[lLine.aP[0]];
        const Vec3D& lC2 = aColors[lLine.aP[1]];

        glColor3f (lC1._x(), lC1._y(), lC1._z());
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
        glColor3f (lC2._x(), lC2._y(), lC2._z());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
      }
      else {
        FV_ASSERT(aColors.size() == 0);
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
      }

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % aPrimitiveOptimizerValue == 0) {
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

// Renders the content of aPoints
void VtxPrimitiveAccumulator::renderPoints()
{
  if (aPoints.size() > 0) {

    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);

    Points::const_iterator       lIterPoints    = aPoints.begin();
    const Points::const_iterator lIterPointsEnd = aPoints.end  ();

    while (lIterPoints != lIterPointsEnd) {

      const Point& lPoint = *lIterPoints;
      FV_ASSERT(lPoint.aP >= 0 && lPoint.aP < static_cast<int>(aVertices.size()));
      const Vec3D& lP = aVertices[lPoint.aP].position.v;

      if (aColors.size() == aVertices.size()) {
        FV_ASSERT(aColors.size() != 0);
        const Vec3D& lC = aColors[lPoint.aP];

        glColor3f (lC._x(), lC._y(), lC._z());
        glVertex3f(lP._x(), lP._y(), lP._z());
      }
      else {
        FV_ASSERT(aColors.size() == 0);
        glVertex3f(lP._x(), lP._y(), lP._z());
      }

      ++lIterPoints;
    }

    glEnd();

    // Revert the lighting state
    glPopAttrib();

  }
}

// Renders the content of aQuads
void VtxPrimitiveAccumulator::renderQuads()
{
  int lPrimitiveCount = 0;

  if (aQuads.size() > 0) {

    glBegin(GL_QUADS);

    const Vec3D lNullVector     (0.0f, 0.0f, 0.0f);
    const Vec3D lArbitraryNormal(1.0f, 0.0f, 0.0f);

    Quads::const_iterator       lIterQuads    = aQuads.begin();
    const Quads::const_iterator lIterQuadsEnd = aQuads.end  ();

    while (lIterQuads != lIterQuadsEnd) {

      const Quad& lQuad = *lIterQuads;
      FV_ASSERT(lQuad.aP[0] >= 0 && lQuad.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[1] >= 0 && lQuad.aP[1] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[2] >= 0 && lQuad.aP[2] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lQuad.aP[3] >= 0 && lQuad.aP[3] < static_cast<int>(aVertices.size()));
      const Vec3D& lP1 = aVertices[lQuad.aP[0]].position.v;
      const Vec3D& lP2 = aVertices[lQuad.aP[1]].position.v;
      const Vec3D& lP3 = aVertices[lQuad.aP[2]].position.v;
      const Vec3D& lP4 = aVertices[lQuad.aP[3]].position.v;
      Vec3D        lN1;
      Vec3D        lN2;
      Vec3D        lN3;
      Vec3D        lN4;

      if (aNormals.size() == aVertices.size()) {
        FV_ASSERT(aNormals.size() != 0);
        lN1 = aNormals[lQuad.aP[0]];
        lN2 = aNormals[lQuad.aP[1]];
        lN3 = aNormals[lQuad.aP[2]];
        lN4 = aNormals[lQuad.aP[3]];
      }
      else {
        FV_ASSERT(aNormals.size() == 0);
        lN1 = (lP2-lP1) % (lP4-lP1);
        lN2 = (lP3-lP2) % (lP1-lP2);
        lN3 = (lP4-lP3) % (lP2-lP3);
        lN4 = (lP1-lP4) % (lP3-lP4);

        // For degenerated Quads (unfortunately, we can't
        // prevent this from happening since it depends
        // on the user's input)
        if (lN1 != lNullVector) {
          lN1.normalize();
        }
        else {
          lN1 = lArbitraryNormal;
        }

        if (lN2 != lNullVector) {
          lN2.normalize();
        }
        else {
          lN2 = lArbitraryNormal;
        }

        if (lN3 != lNullVector) {
          lN3.normalize();
        }
        else {
          lN3 = lArbitraryNormal;
        }

        if (lN4 != lNullVector) {
          lN4.normalize();
        }
        else {
          lN4 = lArbitraryNormal;
        }
      }

      FV_ASSERT(fabs(lN1.Norm() - 1.0f) < 1.0E-6);
      FV_ASSERT(fabs(lN2.Norm() - 1.0f) < 1.0E-6);
      FV_ASSERT(fabs(lN3.Norm() - 1.0f) < 1.0E-6);
      FV_ASSERT(fabs(lN4.Norm() - 1.0f) < 1.0E-6);

      //const Tex2D* lT1 = 0;
      //const Tex2D* lT2 = 0;
      //const Tex2D* lT3 = 0;
      //const Tex2D* lT4 = 0;
      /*if(aTexCoords.size() == aVertices.size()) {
        FV_ASSERT(aTexCoords.size() != 0);

        lT1 = &aTexCoords[lQuad.aP[0]];
        lT2 = &aTexCoords[lQuad.aP[1]];
        lT3 = &aTexCoords[lQuad.aP[2]];
        lT4 = &aTexCoords[lQuad.aP[3]];
      }*/

      if (aColors.size() == aVertices.size()) {
        FV_ASSERT(aColors.size() != 0);
        const Vec3D& lC1 = aColors[lQuad.aP[0]];
        const Vec3D& lC2 = aColors[lQuad.aP[1]];
        const Vec3D& lC3 = aColors[lQuad.aP[2]];
        const Vec3D& lC4 = aColors[lQuad.aP[3]];

        //if(lT1)
          //glTexCoord2f(lT1->u(), lT1->v());
        glNormal3f(lN1._x(), lN1._y(), lN1._z());
        glColor3f (lC1._x(), lC1._y(), lC1._z());
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
       // if(lT2)
        //  glTexCoord2f(lT2->u(), lT2->v());
        glNormal3f(lN2._x(), lN2._y(), lN2._z());
        glColor3f (lC2._x(), lC2._y(), lC2._z());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
      //  if(lT3)
       //   glTexCoord2f(lT3->u(), lT3->v());
        glNormal3f(lN3._x(), lN3._y(), lN3._z());
        glColor3f (lC3._x(), lC3._y(), lC3._z());
        glVertex3f(lP3._x(), lP3._y(), lP3._z());
       // if(lT4)
        //  glTexCoord2f(lT4->u(), lT4->v());
        glNormal3f(lN4._x(), lN4._y(), lN4._z());
        glColor3f (lC4._x(), lC4._y(), lC4._z());
        glVertex3f(lP4._x(), lP4._y(), lP4._z());
      }
      else {
        FV_ASSERT(aColors.size() == 0);
       //// if(lT1)
       //   glTexCoord2f(lT1->u(), lT1->v());
        glNormal3f(lN1._x(), lN1._y(), lN1._z());
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
        /*if(lT2)
          glTexCoord2f(lT2->u(), lT2->v());*/
        glNormal3f(lN2._x(), lN2._y(), lN2._z());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
        /*if(lT3)
          glTexCoord2f(lT3->u(), lT3->v());*/
        glNormal3f(lN3._x(), lN3._y(), lN3._z());
        glVertex3f(lP3._x(), lP3._y(), lP3._z());
       /* if(lT4)
          glTexCoord2f(lT4->u(), lT4->v());*/
        glNormal3f(lN4._x(), lN4._y(), lN4._z());
        glVertex3f(lP4._x(), lP4._y(), lP4._z());
      }

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % aPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_QUADS);
      }
      ++lIterQuads;
    }

    glEnd();

  }
}

void VtxPrimitiveAccumulator::renderSimplified(const RenderParams& pParams)
{
  // Check if we have to update the Simplified model
  if (aSimplifiedDirty) {
    if(pParams.eRenderType == 0) {
      constructSimplified_Points();
    }
    else {
      constructSimplified_Heuristics();
    }
  }

  RenderParams lParams = pParams;
  lParams.eRMode      = RenderParams::eFull;

  // If there is no need for simplication, we
  // just render *this
  if (aSimplifiedSelf) {
    renderFull(lParams);
  }
  else {

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
  }
}

// Renders the content of aTriangles
void VtxPrimitiveAccumulator::renderTriangles()
{
  int lPrimitiveCount = 0;

  if (aTriangles.size() > 0) {

    glBegin(GL_TRIANGLES);

    const Vec3D lNullVector     (0.0f, 0.0f, 0.0f);
    const Vec3D lArbitraryNormal(1.0f, 0.0f, 0.0f);

    Triangles::const_iterator       lIterTriangles    = aTriangles.begin();
    const Triangles::const_iterator lIterTrianglesEnd = aTriangles.end  ();

    while (lIterTriangles != lIterTrianglesEnd) {

      const Triangle& lTriangle = *lIterTriangles;
      FV_ASSERT(lTriangle.aP[0] >= 0 && lTriangle.aP[0] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lTriangle.aP[1] >= 0 && lTriangle.aP[1] < static_cast<int>(aVertices.size()));
      FV_ASSERT(lTriangle.aP[2] >= 0 && lTriangle.aP[2] < static_cast<int>(aVertices.size()));
      const Vec3D& lP1 = aVertices[lTriangle.aP[0]].position.v;
      const Vec3D& lP2 = aVertices[lTriangle.aP[1]].position.v;
      const Vec3D& lP3 = aVertices[lTriangle.aP[2]].position.v;
      Vec3D        lN1;
      Vec3D        lN2;
      Vec3D        lN3;

      if (aNormals.size() == aVertices.size()) {
        FV_ASSERT(aNormals.size() != 0);
        lN1 = aNormals[lTriangle.aP[0]];
        lN2 = aNormals[lTriangle.aP[1]];
        lN3 = aNormals[lTriangle.aP[2]];
      }
      else {
        FV_ASSERT(aNormals.size() == 0);
        lN1 = (lP2-lP1) % (lP3-lP1);

        // For degenerated Quads (unfortunately, we can't
        // prevent this from happening since it depends
        // on the user's input)
        if (lN1 != lNullVector) {
          lN1.normalize();
        }
        else {
          lN1 = lArbitraryNormal;
        }

        lN2 = lN1;
        lN3 = lN1;
      }

      FV_ASSERT(fabs(lN1.Norm() - 1.0f) < 1.0E-6);
      FV_ASSERT(fabs(lN2.Norm() - 1.0f) < 1.0E-6);
      FV_ASSERT(fabs(lN3.Norm() - 1.0f) < 1.0E-6);

      //const Tex2D* lT1 = 0;
      //const Tex2D* lT2 = 0;
      //const Tex2D* lT3 = 0;
      /*if(aTexCoords.size() == aVertices.size()) {
        FV_ASSERT(aTexCoords.size() != 0);

        lT1 = &aTexCoords[lTriangle.aP[0]];
        lT2 = &aTexCoords[lTriangle.aP[1]];
        lT3 = &aTexCoords[lTriangle.aP[2]];
      }*/

      if (aColors.size() == aVertices.size()) {
        FV_ASSERT(aColors.size() != 0);
        const Vec3D& lC1 = aColors[lTriangle.aP[0]];
        const Vec3D& lC2 = aColors[lTriangle.aP[1]];
        const Vec3D& lC3 = aColors[lTriangle.aP[2]];

        //if(lT1)
        //  glTexCoord2f(lT1->u(), lT1->v());
        glNormal3f(lN1._x(), lN1._y(), lN1._z());
        glColor3f (lC1._x(), lC1._y(), lC1._z());
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
/*        if(lT2)
          glTexCoord2f(lT2->u(), lT2->v())*/;
        glNormal3f(lN2._x(), lN2._y(), lN2._z());
        glColor3f (lC2._x(), lC2._y(), lC2._z());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
        //if(lT3)
        //  glTexCoord2f(lT3->u(), lT3->v());
        glNormal3f(lN3._x(), lN3._y(), lN3._z());
        glColor3f (lC3._x(), lC3._y(), lC3._z());
        glVertex3f(lP3._x(), lP3._y(), lP3._z());
      }
      else {
        FV_ASSERT(aColors.size() == 0);
        //if(lT1)
        //  glTexCoord2f(lT1->u(), lT1->v());
        glNormal3f(lN1._x(), lN1._y(), lN1._z());
        glVertex3f(lP1._x(), lP1._y(), lP1._z());
        //if(lT2)
        //  glTexCoord2f(lT2->u(), lT2->v());
        glNormal3f(lN2._x(), lN2._y(), lN2._z());
        glVertex3f(lP2._x(), lP2._y(), lP2._z());
        //if(lT3)
        //  glTexCoord2f(lT3->u(), lT3->v());
        glNormal3f(lN3._x(), lN3._y(), lN3._z());
        glVertex3f(lP3._x(), lP3._y(), lP3._z());
      }

      lPrimitiveCount++;
      // This value sets the maximum length of the number of primitives
      // between a glBegin and a glEnd.  Some openGL driver optimize too
      // much and are slow to execute onto long primitive lists
      if(lPrimitiveCount % aPrimitiveOptimizerValue == 0) {
        glEnd();
        glBegin(GL_TRIANGLES);
      }

      ++lIterTriangles;
    }

    glEnd();

  }
}

struct SimplificationPoint
{
	SimplificationPoint(
		const float pX, 
		const float pY,
        const float pZ)
	: aBarycenterDistance (FLT_MAX),
      aPoint              (pX, pY, pZ),
      aPrimiticeBarycenter()
    {}

  float     aBarycenterDistance;
  Vec3D  aPoint;
  Vec3D  aPrimiticeBarycenter;

};

template <class Container, class Vertices>
inline void scanForClosest(Container&           pContainer,
                           SimplificationPoint* pIterSPBegin,
                           SimplificationPoint* pIterSPEnd,
                           const Vertices&      aVertices)
{
  const typename Container::size_type lSize = pContainer.size();

  for (typename Container::size_type i=0; i<lSize; ++i) {

    const typename Container::value_type&  lPrimitive = pContainer[i];
    SimplificationPoint*                   lIterSP    = pIterSPBegin;

    while (lIterSP != pIterSPEnd) {
      SimplificationPoint& lSP         = *lIterSP;
      const Vec3D       lBarycenter = lPrimitive .getBarycenter(aVertices);
      const float          lDistance   = lBarycenter.getDistance(lSP.aPoint);
      if (lDistance < lSP.aBarycenterDistance) {
        lSP.aPrimiticeBarycenter = lBarycenter;
        lSP.aBarycenterDistance  = lDistance;
      }
      ++lIterSP;
    }
  }
}


void VtxPrimitiveAccumulator::constructSimplified_Points()
{

  // The simplified model is a set of points located
  // at the barycenter of a subset of the primitives.
  // These primitives are chosen so that the
  // object shape can still be recognise by the user
  // with a modest number of points.
  // More precisely, they are the primitives "closest"
  // to a set of points on a the circumscribed sphere to
  // the BoundingBox.

  if (!aSimplifiedDirty) {
    // Just to be safe, should not pass here
    // if nobody calls it when aSimplified is
    // not used
  }
  else {

    // Start from scratch
    FV_ASSERT(aSimplified != 0);
    delete aSimplified;
    aSimplified     = new PrimitiveAccumulator(false);
    aSimplifiedSelf = false;

    // We call getBoundingBox so that it gets updated
    // if it needs to
    const BBox3D& lBoundingBox = getBBox3D();
    const Vec3D     lCenter(lBoundingBox.getCenter().v);
    const float        lRadius      = lBoundingBox.getCircumscribedSphereRadius();

    SimplificationPoint unitSphereSampling[] = {
      SimplificationPoint(-1.0,0.0,0.0),
      SimplificationPoint(-0.894427191,-0.4472135955,0.0),
      SimplificationPoint(-0.894427191,0.0,-0.4472135955),
      SimplificationPoint(-0.894427191,0.0,0.4472135955),
      SimplificationPoint(-0.894427191,0.4472135955,0.0),
      SimplificationPoint(-0.816496580928,-0.408248290464,-0.408248290464),
      SimplificationPoint(-0.816496580928,-0.408248290464,0.408248290464),
      SimplificationPoint(-0.816496580928,0.408248290464,-0.408248290464),
      SimplificationPoint(-0.816496580928,0.408248290464,0.408248290464),
      SimplificationPoint(-0.707106781187,-0.707106781187,0.0),
      SimplificationPoint(-0.707106781187,0.0,-0.707106781187),
      SimplificationPoint(-0.707106781187,0.0,0.707106781187),
      SimplificationPoint(-0.707106781187,0.707106781187,0.0),
      SimplificationPoint(-0.666666666667,-0.333333333333,-0.666666666667),
      SimplificationPoint(-0.666666666667,-0.333333333333,0.666666666667),
      SimplificationPoint(-0.666666666667,-0.666666666667,-0.333333333333),
      SimplificationPoint(-0.666666666667,-0.666666666667,0.333333333333),
      SimplificationPoint(-0.666666666667,0.333333333333,-0.666666666667),
      SimplificationPoint(-0.666666666667,0.333333333333,0.666666666667),
      SimplificationPoint(-0.666666666667,0.666666666667,-0.333333333333),
      SimplificationPoint(-0.666666666667,0.666666666667,0.333333333333),
      SimplificationPoint(-0.57735026919,-0.57735026919,-0.57735026919),
      SimplificationPoint(-0.57735026919,-0.57735026919,0.57735026919),
      SimplificationPoint(-0.57735026919,0.57735026919,-0.57735026919),
      SimplificationPoint(-0.57735026919,0.57735026919,0.57735026919),
      SimplificationPoint(-0.4472135955,-0.894427191,0.0),
      SimplificationPoint(-0.4472135955,0.0,-0.894427191),
      SimplificationPoint(-0.4472135955,0.0,0.894427191),
      SimplificationPoint(-0.4472135955,0.894427191,0.0),
      SimplificationPoint(-0.408248290464,-0.408248290464,-0.816496580928),
      SimplificationPoint(-0.408248290464,-0.408248290464,0.816496580928),
      SimplificationPoint(-0.408248290464,-0.816496580928,-0.408248290464),
      SimplificationPoint(-0.408248290464,-0.816496580928,0.408248290464),
      SimplificationPoint(-0.408248290464,0.408248290464,-0.816496580928),
      SimplificationPoint(-0.408248290464,0.408248290464,0.816496580928),
      SimplificationPoint(-0.408248290464,0.816496580928,-0.408248290464),
      SimplificationPoint(-0.408248290464,0.816496580928,0.408248290464),
      SimplificationPoint(-0.333333333333,-0.666666666667,-0.666666666667),
      SimplificationPoint(-0.333333333333,-0.666666666667,0.666666666667),
      SimplificationPoint(-0.333333333333,0.666666666667,-0.666666666667),
      SimplificationPoint(-0.333333333333,0.666666666667,0.666666666667),
      SimplificationPoint(0.0,-0.4472135955,-0.894427191),
      SimplificationPoint(0.0,-0.4472135955,0.894427191),
      SimplificationPoint(0.0,-0.707106781187,-0.707106781187),
      SimplificationPoint(0.0,-0.707106781187,0.707106781187),
      SimplificationPoint(0.0,-0.894427191,-0.4472135955),
      SimplificationPoint(0.0,-0.894427191,0.4472135955),
      SimplificationPoint(0.0,-1.0,0.0),
      SimplificationPoint(0.0,0.0,-1.0),
      SimplificationPoint(0.0,0.0,1.0),
      SimplificationPoint(0.0,0.4472135955,-0.894427191),
      SimplificationPoint(0.0,0.4472135955,0.894427191),
      SimplificationPoint(0.0,0.707106781187,-0.707106781187),
      SimplificationPoint(0.0,0.707106781187,0.707106781187),
      SimplificationPoint(0.0,0.894427191,-0.4472135955),
      SimplificationPoint(0.0,0.894427191,0.4472135955),
      SimplificationPoint(0.0,1.0,0.0),
      SimplificationPoint(0.333333333333,-0.666666666667,-0.666666666667),
      SimplificationPoint(0.333333333333,-0.666666666667,0.666666666667),
      SimplificationPoint(0.333333333333,0.666666666667,-0.666666666667),
      SimplificationPoint(0.333333333333,0.666666666667,0.666666666667),
      SimplificationPoint(0.408248290464,-0.408248290464,-0.816496580928),
      SimplificationPoint(0.408248290464,-0.408248290464,0.816496580928),
      SimplificationPoint(0.408248290464,-0.816496580928,-0.408248290464),
      SimplificationPoint(0.408248290464,-0.816496580928,0.408248290464),
      SimplificationPoint(0.408248290464,0.408248290464,-0.816496580928),
      SimplificationPoint(0.408248290464,0.408248290464,0.816496580928),
      SimplificationPoint(0.408248290464,0.816496580928,-0.408248290464),
      SimplificationPoint(0.408248290464,0.816496580928,0.408248290464),
      SimplificationPoint(0.4472135955,-0.894427191,0.0),
      SimplificationPoint(0.4472135955,0.0,-0.894427191),
      SimplificationPoint(0.4472135955,0.0,0.894427191),
      SimplificationPoint(0.4472135955,0.894427191,0.0),
      SimplificationPoint(0.57735026919,-0.57735026919,-0.57735026919),
      SimplificationPoint(0.57735026919,-0.57735026919,0.57735026919),
      SimplificationPoint(0.57735026919,0.57735026919,-0.57735026919),
      SimplificationPoint(0.57735026919,0.57735026919,0.57735026919),
      SimplificationPoint(0.666666666667,-0.333333333333,-0.666666666667),
      SimplificationPoint(0.666666666667,-0.333333333333,0.666666666667),
      SimplificationPoint(0.666666666667,-0.666666666667,-0.333333333333),
      SimplificationPoint(0.666666666667,-0.666666666667,0.333333333333),
      SimplificationPoint(0.666666666667,0.333333333333,-0.666666666667),
      SimplificationPoint(0.666666666667,0.333333333333,0.666666666667),
      SimplificationPoint(0.666666666667,0.666666666667,-0.333333333333),
      SimplificationPoint(0.666666666667,0.666666666667,0.333333333333),
      SimplificationPoint(0.707106781187,-0.707106781187,0.0),
      SimplificationPoint(0.707106781187,0.0,-0.707106781187),
      SimplificationPoint(0.707106781187,0.0,0.707106781187),
      SimplificationPoint(0.707106781187,0.707106781187,0.0),
      SimplificationPoint(0.816496580928,-0.408248290464,-0.408248290464),
      SimplificationPoint(0.816496580928,-0.408248290464,0.408248290464),
      SimplificationPoint(0.816496580928,0.408248290464,-0.408248290464),
      SimplificationPoint(0.816496580928,0.408248290464,0.408248290464),
      SimplificationPoint(0.894427191,-0.4472135955,0.0),
      SimplificationPoint(0.894427191,0.0,-0.4472135955),
      SimplificationPoint(0.894427191,0.0,0.4472135955),
      SimplificationPoint(0.894427191,0.4472135955,0.0),
      SimplificationPoint(1.0,0.0,0.0)
    };

//     SimplificationPoint unitSphereSampling[] = {
//       SimplificationPoint(-1.0,0.0,0.0),
//       SimplificationPoint(-0.707106781187,-0.707106781187,0.0),
//       SimplificationPoint(-0.707106781187,0.0,-0.707106781187),
//       SimplificationPoint(-0.707106781187,0.0,0.707106781187),
//       SimplificationPoint(-0.707106781187,0.707106781187,0.0),
//       SimplificationPoint(-0.57735026919,-0.57735026919,-0.57735026919),
//       SimplificationPoint(-0.57735026919,-0.57735026919,0.57735026919),
//       SimplificationPoint(-0.57735026919,0.57735026919,-0.57735026919),
//       SimplificationPoint(-0.57735026919,0.57735026919,0.57735026919),
//       SimplificationPoint(0.0,-0.707106781187,-0.707106781187),
//       SimplificationPoint(0.0,-0.707106781187,0.707106781187),
//       SimplificationPoint(0.0,-1.0,0.0),
//       SimplificationPoint(0.0,0.0,-1.0),
//       SimplificationPoint(0.0,0.0,1.0),
//       SimplificationPoint(0.0,0.707106781187,-0.707106781187),
//       SimplificationPoint(0.0,0.707106781187,0.707106781187),
//       SimplificationPoint(0.0,1.0,0.0),
//       SimplificationPoint(0.57735026919,-0.57735026919,-0.57735026919),
//       SimplificationPoint(0.57735026919,-0.57735026919,0.57735026919),
//       SimplificationPoint(0.57735026919,0.57735026919,-0.57735026919),
//       SimplificationPoint(0.57735026919,0.57735026919,0.57735026919),
//       SimplificationPoint(0.707106781187,-0.707106781187,0.0),
//       SimplificationPoint(0.707106781187,0.0,-0.707106781187),
//       SimplificationPoint(0.707106781187,0.0,0.707106781187),
//       SimplificationPoint(0.707106781187,0.707106781187,0.0),
//       SimplificationPoint(1.0,0.0,0.0)
//     };

    SimplificationPoint* lBegin = unitSphereSampling;
    SimplificationPoint* lIter  = lBegin;
    SimplificationPoint* lEnd   = unitSphereSampling + (sizeof(unitSphereSampling)/sizeof(SimplificationPoint));

    // Scale the point from the unit sphere to the circumscribed sphere
    FV_ASSERT(lBegin != 0);
    FV_ASSERT(lEnd   != 0);

    while (lIter != lEnd) {
      FV_ASSERT(lIter != 0);
      lIter->aPoint *= lRadius;
      lIter->aPoint += lCenter;
      ++lIter;
    }

    // Scan all the primitives to find the closest ones
    // For now, we choose only one primitive in all the types
    // for each point on the shere. Again, we'll have to see
    // if if wouldn't be better to have a separate simplication
    // for each type of primitive.  We'll see.
    // Scan only the new primitives for each type

    scanForClosest(aLines,     lBegin, lEnd, aVertices);
    scanForClosest(aPoints,    lBegin, lEnd, aVertices);
    scanForClosest(aQuads,     lBegin, lEnd, aVertices);
    scanForClosest(aTriangles, lBegin, lEnd, aVertices);

    // Now add Points to aSimplified for each primitive barycenter chosen
    lIter = lBegin;

    while (lIter != lEnd) {
      FV_ASSERT(lIter != 0);
      FV_ASSERT(aSimplified != 0);
	  aSimplified->AddPoint(lIter->aPrimiticeBarycenter);
      ++lIter;
    }
  }

  // aSimplified is not dirty anymore
  aSimplifiedDirty = false;
}

#define SIMPLIFICATION_FACTOR 43
#define SIMPLIFICATION_RANDOM_FACTOR 13
void VtxPrimitiveAccumulator::constructSimplified_Heuristics()
{

  // The simplified model is a set of points located
  // at the barycenter of a subset of the primitives.
  // These primitives are chosen so that the
  // object shape can still be recognise by the user
  // with a modest number of points.
  // More precisely, they are the primitives "closest"
  // to a set of points on a the circumscribed sphere to
  // the BoundingBox.

  if (!aSimplifiedDirty) {
    // Just to be safe, should not pass here
    // if nobody calls it when aSimplified is
    // not used
  }
  else {

    // Start from scratch
    FV_ASSERT(aSimplified != 0);
    delete aSimplified;
    aSimplified     = new PrimitiveAccumulator(false);
    aSimplifiedSelf = false;

    // FOR TRIANGLES AND QUADS, WE BUILD A CONNECTIVITY AND DETECT
    // WHICH ARE "BOUNDARY" TRIANGLES, AND BETWEEN WHICH TRIANGLES
    // THERE IS A HARD BREAK (60 degr angle and more)
    //if(aTriangles.size() || aQuads.size()) {
    //  Connectivity lConnectivity;
    //  lConnectivity.reserveData(aVertices.size(),aTriangles.size()+aQuads.size());
    //  const int lNbTri = aTriangles.size();
    //  for(int i=0;i<lNbTri;i++) {
    //    lConnectivity.addPolygon(i,3,aTriangles[i].aP);
    //  }
    //  const int lNbQuad = aQuads.size();
    //  {for(int i=0;i<lNbQuad;i++) {
    //    lConnectivity.addPolygon(lNbTri+i,4,aQuads[i].aP);
    //  }}

    //  int lC[2] = { -1,-1 };
    //  const double lDetectionAngle = 0.52;

    //  {for(int i=0;i<lNbTri;i++) {
    //    Triangle& lTri = aTriangles[i];
    //    int match = lConnectivity.getEdge2Connectivity(lTri.aP[0],lTri.aP[1],lC);
    //    if(lC[0] != i && lC[1] != i) {
    //      continue;
    //    }
    //    double lAngle = 0;
    //    if(match == 2) {
    //      Triangle& lTri2 = lC[0]==i?aTriangles[lC[1]]:aTriangles[lC[0]];
    //      lAngle = lTri.getNormal(aVertices).getAngle(lTri2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2 || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lTri.aP[0]],aVertices[lTri.aP[1]]);
    //    }
    //    lAngle = 0;
    //    match = lConnectivity.getEdge2Connectivity(lTri.aP[1],lTri.aP[2],lC);
    //    if(match == 2) {
    //      Triangle& lTri2 = lC[0]==i?aTriangles[lC[1]]:aTriangles[lC[0]];
    //      lAngle = lTri.getNormal(aVertices).getAngle(lTri2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2  || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lTri.aP[1]],aVertices[lTri.aP[2]]);
    //    }
    //    lAngle = 0;
    //    match = lConnectivity.getEdge2Connectivity(lTri.aP[0],lTri.aP[2],lC);
    //    if(match == 2) {
    //      Triangle& lTri2 = lC[0]==i?aTriangles[lC[1]]:aTriangles[lC[0]];
    //      lAngle = lTri.getNormal(aVertices).getAngle(lTri2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2 || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lTri.aP[0]],aVertices[lTri.aP[2]]);
    //    }
    //  }}
    //  {for(int i=0;i<lNbQuad;i++) {
    //    Quad& lQuad = aQuads[i];
    //    int match = lConnectivity.getEdge2Connectivity(lQuad.aP[0],lQuad.aP[1],lC);
    //    if(lC[0] != -1) { lC[0] -= lNbTri; } if(lC[1] != -1) { lC[1] -= lNbTri; }
    //    if(lC[0] != i && lC[1] != i) {
    //      continue;
    //    }
    //    double lAngle = 0;
    //    if(match == 2) {
    //      Quad& lQuad2 = lC[0]==i?aQuads[lC[1]]:aQuads[lC[0]];
    //      lAngle = lQuad.getNormal(aVertices).getAngle(lQuad2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2 || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lQuad.aP[0]],aVertices[lQuad.aP[1]]);
    //    }
    //    lAngle = 0;
    //    match = lConnectivity.getEdge2Connectivity(lQuad.aP[1],lQuad.aP[2],lC);
    //    if(lC[0] != -1) { lC[0] -= lNbTri; } if(lC[1] != -1) { lC[1] -= lNbTri; }
    //    if(match == 2) {
    //      Quad& lQuad2 = lC[0]==i?aQuads[lC[1]]:aQuads[lC[0]];
    //      lAngle = lQuad.getNormal(aVertices).getAngle(lQuad2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2  || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lQuad.aP[1]],aVertices[lQuad.aP[2]]);
    //    }
    //    lAngle = 0;
    //    match = lConnectivity.getEdge2Connectivity(lQuad.aP[2],lQuad.aP[3],lC);
    //    if(lC[0] != -1) { lC[0] -= lNbTri; } if(lC[1] != -1) { lC[1] -= lNbTri; }
    //    if(match == 2) {
    //      Quad& lQuad2 = lC[0]==i?aQuads[lC[1]]:aQuads[lC[0]];
    //      lAngle = lQuad.getNormal(aVertices).getAngle(lQuad2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2 || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lQuad.aP[2]],aVertices[lQuad.aP[3]]);
    //    }
    //    lAngle = 0;
    //    match = lConnectivity.getEdge2Connectivity(lQuad.aP[0],lQuad.aP[3],lC);
    //    if(lC[0] != -1) { lC[0] -= lNbTri; } if(lC[1] != -1) { lC[1] -= lNbTri; }
    //    if(match == 2) {
    //      Quad& lQuad2 = lC[0]==i?aQuads[lC[1]]:aQuads[lC[0]];
    //      lAngle = lQuad.getNormal(aVertices).getAngle(lQuad2.getNormal(aVertices)) ;
    //    }
    //    if(match != 2 || lAngle > lDetectionAngle) {
    //      aSimplified->addLine(aVertices[lQuad.aP[0]],aVertices[lQuad.aP[3]]);
    //    }
    //  }}
    //}

    //// FOR OTHER PRIMITIVES (LINES, POINTS)
    //// WE SHOW 1/x
    ////  the factor is 87 instead of 100, because that on structured data with a modulo
    ////  of 100, the chosent points does not seem "random"
    //{
    //  int lNbPoints = aPoints.size();
    //  for(int i=0;i<lNbPoints;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
    //    aSimplified->addPoint(aVertices[aPoints[i].aP]);
    //  }
    //}
    //{
    //  int lNbLines = aLines.size();
    //  for(int i=0;i<lNbLines;i+=SIMPLIFICATION_FACTOR+rand()%SIMPLIFICATION_RANDOM_FACTOR) {
    //    aSimplified->addLine(aVertices[aLines[i].aP[0]],aVertices[aLines[i].aP[1]]);
    //  }
    //}
  }

  // aSimplified is not dirty anymore
  aSimplifiedDirty = false;
}

} // edn namespace FemViewer
