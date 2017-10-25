#ifndef _PRIMITIVE_ACCUMULATOR_H_
#define _PRIMITIVE_ACCUMULATOR_H_

#include "fv_inc.h"
#include "BBox3D.h"
#include "RenderParams.h"
#include "MathHelper.h"

#include <string>
#include <vector>


namespace FemViewer {

using namespace fvmath;


// Class used to accumulate OpenGL
// primitives and create an optimized order to
// display them
class PrimitiveAccumulator
{
public:

	PrimitiveAccumulator (const bool pCreateSimplified);
	~PrimitiveAccumulator();

	void AddArrow(const Vec3D& p1,const Vec3D& p2,
                  const float tipprop, const int nrpolygons);

	void AddArrowColored(const Vec3D& p1, const ColorRGB& col1,
                         const Vec3D& p2, const ColorRGB& col2,
                         const float tipprop, const int nrpolygons);

  //void  AddSolidCube       (const Vec3D& pCenter,
  //                          const float     pHalfSize);

  //void  AddSolidSphere     (const Vec3D& pCenter,
  //                          const float     pRadius,
  //                          const int       pSlices,
  //                          const int       pStacks);

  //void  AddWireCube        (const Vec3D& pCenter,
  //                          const float     pHalfSize);

	void AddLine(const Vec3D& p1, const Vec3D& p2);

	void AddLineColored(const Vec3D& p1, const ColorRGB& col1,
                        const Vec3D& p2, const ColorRGB& col2);

	void AddPoint(const Vec3D& pt);

	void AddPointColored(const Vec3D& pt, const ColorRGB& col);

	void AddQuad(const Vec3D& p1, const Vec3D& p2,
                 const Vec3D& p3, const Vec3D& p4);

	void AddQuadColored(const Vec3D& p1, const ColorRGB& col1,
                        const Vec3D& p2, const ColorRGB& col2,
                        const Vec3D& p3, const ColorRGB& col3,
                        const Vec3D& p4, const ColorRGB& col4);

	void AddTriangle(const Vec3D& p1, const Vec3D& p2,
                     const Vec3D& p3);

	void AddTriangleColored(const Vec3D& p1, const ColorRGB& col1,
                            const Vec3D& p2, const ColorRGB& col2,
                            const Vec3D& p3, const ColorRGB& col3);

	void AddTriangleNormalColored(const float nv[], const float p1[], const float p2[], const float p3[],
								  const ColorRGB& col1, const ColorRGB& col2, const ColorRGB& col3);

	void DumpCharacteristics(std::ostream& pOstream,
                             const std::string& pIndentation,
                             const Matrix<float>& pTransformation);

	const BBox3D&  GetBBox3D() const;
	
	void Render(const RenderParams& pParams);

private:

	// Block the use of those
	PrimitiveAccumulator();
	PrimitiveAccumulator(const PrimitiveAccumulator&);
	PrimitiveAccumulator& operator=(const PrimitiveAccumulator&);

	struct Line {
		Vec3D pt[2];

		Vec3D getBarycenter() const {
			return 0.5*(pt[0] + pt[1]);
		}
	};

	struct LineColored : public Line {
		ColorRGB col1;
		ColorRGB col2;
	};

	struct Point {
		Vec3D pt;

		Vec3D getBarycenter() const {
			return pt;
		}
	};

	struct PointColored : public Point {
		ColorRGB col;
	};

	struct Quad {
		Vec3D pt[4];

		Vec3D getBarycenter() const {
			return 0.25*(pt[0] + pt[1] + pt[2] + pt[3]);
		}

		void renderFacetsFrame() const {
			glBegin(GL_LINE_LOOP);
				glVertex3f(pt[0]._x(), pt[0]._y(), pt[0]._z());
				glVertex3f(pt[1]._x(), pt[1]._y(), pt[1]._z());
				glVertex3f(pt[2]._x(), pt[2]._y(), pt[2]._z());
				glVertex3f(pt[3]._x(), pt[3]._y(), pt[3]._z());
			glEnd();
		}
	};

	struct QuadColored : public Quad {
		ColorRGB col1, col2, col3, col4;
	};


	struct QuadNormals : public Quad {
		Vec3D n1, n2, n3, n4;
	};

	struct QuadNormalsColored : public QuadNormals {
		ColorRGB col1, col2, col3, col4;;
	};

	struct TriEelement {
		Vec3f pt[3];

		Vec3f getBarycenter() const {
			return (1.0f/3.0f)*(pt[0] + pt[1] + pt[2]);
		}

		void renderFacetFrame() const {
			glBegin(GL_LINE_LOOP);
				glVertex3fv(pt[0].v);
				glVertex3fv(pt[1].v);
				glVertex3fv(pt[2].v);
			glEnd();
		}
	};

	//struct TriEelemntColor : public TriEelemnt {
	//	ColorRGB c1, c2, c3;
	//};

	struct TriEelementNormal : public TriEelement {
		Vec3f nv;
	};

	struct TriEelementNormalColor : public TriEelementNormal {
		ColorRGB c1, c2, c3;
	};

	struct Triangle {
		Vec3D pt[3];

		Vec3D getBarycenter() const {
			return (1.0f/3.0f)*(pt[0] + pt[1] + pt[2]);
		}

		void renderFacetsFrame() const {
			glBegin(GL_LINE_LOOP);
				glVertex3f(pt[0]._x(), pt[0]._y(), pt[0]._z());
				glVertex3f(pt[1]._x(), pt[1]._y(), pt[1]._z());
				glVertex3f(pt[2]._x(), pt[2]._y(), pt[2]._z());
			glEnd();
		}
	};

	struct TriangleColored : public Triangle {
		ColorRGB col1, col2, col3;
	};

	struct TriangleNormals : public Triangle {
		Vec3D n1, n2, n3;
	};

	struct TriangleNormalsColored : public TriangleNormals {
		ColorRGB col1, col2, col3;
	};


	void AddArrowTip(const Vec3D& p1, const Vec3D& p2, const ColorRGB& col,
                     const bool colored, const float tipprop, const int nrpolyons);

	//void  constructSimplified_points    ();
	//void  constructSimplified_heuristics();
	void  renderFacetsFrame            (const RenderParams& pParams);
	void  renderFull                   (const RenderParams& pParams);
	void  renderLines                  ();
	void  renderLinesColored           ();
	void  renderPoints                 ();
	void  renderPointsColored          ();
	void  renderQuads                  ();
	void  renderQuadsColored           ();
	void  renderQuadsNormals           ();
	void  renderQuadsNormalsColored    ();
	void  renderSimplified             (const RenderParams& pParams);
	void  renderTriangles              ();
	void  renderTrianglesColored       ();
	void  renderTrianglesNormals       ();
	void  renderTrianglesNormalsColored();
	//void  renderTriangNormalColorf    () const;
	void  renderTriElementNormalColor  () const;

	

  typedef  std::vector<Line>                    Lines;
  typedef  std::vector<LineColored>             LinesColored;
  typedef  std::vector<Point>                   Points;
  typedef  std::vector<PointColored>            PointsColored;
  typedef  std::vector<Quad>                    Quads;
  typedef  std::vector<QuadColored>             QuadsColored;
  typedef  std::vector<QuadNormals>             QuadsNormals;
  typedef  std::vector<QuadNormalsColored>      QuadsNormalsColored;
  typedef  std::vector<Triangle>                Triangles;
  typedef  std::vector<TriangleColored>         TrianglesColored;
  typedef  std::vector<TriangleNormals>         TrianglesNormals;
  typedef  std::vector<TriangleNormalsColored>  TrianglesNormalsColored;
  typedef  std::vector<TriEelementNormalColor>	TriElemsNormalColor;
  typedef  Lines::size_type                     SizeType;

  mutable BBox3D		   aBBox3D;
  Lines                    aLines;
  mutable SizeType         aLinesBBoxCounter;
  LinesColored             aLinesColored;
  mutable SizeType         aLinesColoredBBoxCounter;
  Points                   aPoints;
  mutable SizeType         aPointsBBoxCounter;
  PointsColored            aPointsColored;
  mutable SizeType         aPointsColoredBBoxCounter;
  Quads                    aQuads;
  mutable SizeType         aQuadsBBoxCounter;
  QuadsColored             aQuadsColored;
  mutable SizeType         aQuadsColoredBBoxCounter;
  QuadsNormals             aQuadsNormals;
  mutable SizeType         aQuadsNormalsBBoxCounter;
  QuadsNormalsColored      aQuadsNormalsColored;
  mutable SizeType         aQuadsNormalsColoredBBoxCounter;
  PrimitiveAccumulator*    aSimplified;
  bool                     aSimplifiedDirty;
  Triangles                aTriangles;
  mutable SizeType         aTrianglesBBoxCounter;
  TrianglesColored         aTrianglesColored;
  mutable SizeType         aTrianglesColoredBBoxCounter;
  TrianglesNormals         aTrianglesNormals;
  mutable SizeType         aTrianglesNormalsBBoxCounter;
  TrianglesNormalsColored  aTrianglesNormalsColored;
  mutable SizeType         aTrianglesNormalsColoredBBoxCounter;
  TriElemsNormalColor	   aTriElemsNormalColor;
  mutable SizeType         aTriElemsNormalColorBBoxCounter;
  
  int                      iPrimitiveOptimizerValue;

};

} // end namespace FemViewer
#endif /* _PRIMITIVE_ACCUMULATOR_H_
*/
