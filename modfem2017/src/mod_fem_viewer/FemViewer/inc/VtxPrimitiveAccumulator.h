#ifndef _VTX_PRIMITIVE_ACCUMULATOR_H_
#define _VTX_PRIMITIVE_ACCUMULATOR_H_

#include "fv_assert.h"
#include "BBox3D.h"
#include "fv_inc.h"
#include "RenderParams.h"
#include "VtxAccumulator.h"
#include <string>
#include <vector>

namespace FemViewer {
class PrimitiveAccumulator;

// Class used to accumulate OpenGL based on a mesh
// (vertices) and create an optimized order to display
// them
class VtxPrimitiveAccumulator
{
public:

  VtxPrimitiveAccumulator (VtxAccumulator& pVertexes);
  ~VtxPrimitiveAccumulator();

  bool  addLine            (const int           pVertexIndex1,
                            const int           pVertexIndex2);

  bool  addPoint           (const int           pVertexIndex);

  bool  addQuad            (const int           pVertexIndex1,
                            const int           pVertexIndex2,
                            const int           pVertexIndex3,
                            const int           pVertexIndex4);

  bool  addTriangle        (const int           pVertexIndex1,
                            const int           pVertexIndex2,
                            const int           pVertexIndex3);

  void  dumpCharacteristics(std::ostream&       pOstream,
                            const std::string&  pIndentation,
                            const Matrix<float>&    pTransformation);

  const  BBox3D&			getBBox3D() const;

  void						render  (const RenderParams& pParams);

private:

  typedef VtxAccumulator::Vertices Vertices;
  typedef VtxAccumulator::Normals Normals;
  // Block the use of those
  VtxPrimitiveAccumulator(const VtxPrimitiveAccumulator&);
  VtxPrimitiveAccumulator& operator=(const VtxPrimitiveAccumulator&);

  struct Line {
    int aP[2];

    Vec3D getBarycenter(const Vertices& pVertices) const {
      FV_ASSERT(aP[0] >= 0 && aP[0] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[1] >= 0 && aP[1] < static_cast<int>(pVertices.size()));
      return 0.5*Vec3D((pVertices[aP[0]].position + pVertices[aP[1]].position).v);
    }
  };

  struct Point {
    int aP;

    Vec3D getBarycenter(const Vertices& pVertices) const {
      FV_ASSERT(aP >= 0 && aP < static_cast<int>(pVertices.size()));
      return Vec3D(pVertices[aP].position.v);
    }
  };

  struct Quad {
    int aP[4];

    Vec3D getBarycenter(const Vertices& pVertices) const {
      FV_ASSERT(aP[0] >= 0 && aP[0] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[1] >= 0 && aP[1] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[2] >= 0 && aP[2] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[3] >= 0 && aP[3] < static_cast<int>(pVertices.size()));
      return 0.25*Vec3D((pVertices[aP[0]].position
    		  + pVertices[aP[1]].position
    		  + pVertices[aP[2]].position
    		  + pVertices[aP[3]].position).v);
    }
    Vec3D getNormal(const Vertices& pVertices) const {
      FV_ASSERT(aP[0] >= 0 && aP[0] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[1] >= 0 && aP[1] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[3] >= 0 && aP[3] < static_cast<int>(pVertices.size()));
      const Vec3D& lV0 = Vec3D(pVertices[aP[0]].position.v);
      const Vec3D& lV1 = Vec3D(pVertices[aP[1]].position.v);
      const Vec3D& lV2 = Vec3D(pVertices[aP[3]].position.v);
      return Vec3D(lV1-lV0) % (lV2-lV0);
    };
  };

  struct Triangle {
    int aP[3];

    Vec3D getBarycenter(const Vertices& pVertices) const {
      FV_ASSERT(aP[0] >= 0 && aP[0] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[1] >= 0 && aP[1] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[2] >= 0 && aP[2] < static_cast<int>(pVertices.size()));
      return (1.0/3.0)*Vec3D((pVertices[aP[0]].position + pVertices[aP[1]].position + pVertices[aP[2]].position).v);
    }
    Vec3D getNormal(const Vertices& pVertices) const {
      FV_ASSERT(aP[0] >= 0 && aP[0] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[1] >= 0 && aP[1] < static_cast<int>(pVertices.size()));
      FV_ASSERT(aP[2] >= 0 && aP[2] < static_cast<int>(pVertices.size()));
      const Vec3D& lV0 = Vec3D(pVertices[aP[0]].position.v);
      const Vec3D& lV1 = Vec3D(pVertices[aP[1]].position.v);
      const Vec3D& lV2 = Vec3D(pVertices[aP[2]].position.v);
      return Vec3D(lV1-lV0) % (lV2-lV0);
    };
  };

  typedef  std::vector<Line>      Lines;
  typedef  std::vector<Point>     Points;
  typedef  std::vector<Quad>      Quads;
  typedef  std::vector<Triangle>  Triangles;
  typedef  Lines::size_type       SizeType;


  void  computeNormals        ();
  void  constructSimplified_Points();
  void  constructSimplified_Heuristics();
  void  renderFacetsFrame     (const RenderParams& pParams);
  void  renderFull            (const RenderParams& pParams);
  void  renderLines           ();
  void  renderLinesColored    ();
  void  renderPoints          ();
  void  renderPointsColored   ();
  void  renderQuads           ();
  void  renderQuadsColored    ();
  void  renderSimplified      (const RenderParams& pParams);
  void  renderTriangles       ();
  void  renderTrianglesColored();


  mutable BBox3D		 aBBox3D;
  const VtxAccumulator::Colors& aColors;
  Lines                  aLines;
  mutable SizeType       aLinesBBoxCounter;
  Points                 aPoints;
  mutable SizeType       aPointsBBoxCounter;
  Quads                  aQuads;
  mutable SizeType       aQuadsBBoxCounter;
  VtxAccumulator::Normals& aNormals;
  //const VtxAccumulator::TexCoords& aTexCoords;
  PrimitiveAccumulator*  aSimplified;
  bool                   aSimplifiedDirty;
  bool                   aSimplifiedSelf;
  Triangles              aTriangles;
  mutable SizeType       aTrianglesBBoxCounter;
  const VtxAccumulator::Vertices& aVertices;
  int                    aPrimitiveOptimizerValue;

  const VtxAccumulator&  aVtxAccumulator;
};

} // end namespace FemViewer


#endif /* _VTX_PRIMITIVE_ACCUMULATOR_H_
*/
