#ifndef _OBJECT_H_
#define _OBJECT_H_

#include "fv_inc.h"
#include "BBox3D.h"
#include "RenderParams.h"
#include "MathHelper.h"

/* Include system */
#include <map>
#include <string>
#include <vector>


namespace FemViewer {

class PrimitiveAccumulator;
class VtxAccumulator;
class VtxPrimitiveAccumulator;
class Vec3D;
class ColorRGB;

// Class used to manage the data associated with
// an Object displayed in OpenGL
class Object
{
	friend class ViewManager;
	public:
		explicit Object(const char* name = "Unknown",bool usevbo = false);
		~Object();
		void Reset();
		void AddMeshNode(const fvmath::Vec3f& vertex);
		void AddIndexOfNode(const unsigned int idx,const int dest);
		void AddColorOfNode(const fvmath::Vec3f& col);
		void AddNormalToNode(const fvmath::Vec3f& norm); 
		void AddEdge(const unsigned int Edge[]);

		void AddArrow(const Vec3D& p1,const Vec3D& p2,
                      const float tipprop, const int nrpolygons);

		void AddArrowColored(const Vec3D& p1, const ColorRGB& col1,
                             const Vec3D& p2, const ColorRGB& col2,
                             const float tipprop, const int nrpolygons);

		void AddPoint(const Vec3D& pt);

		void AddPointColored(const Vec3D& pt, const ColorRGB& col);

		void AddLine(const Vec3D& p1, const Vec3D& p2);

		void AddLineColored(const Vec3D& p1, const ColorRGB& col1,
							const Vec3D& p2, const ColorRGB& col2);

		void AddTriangle(const Vec3D& p1, const Vec3D& p2,
						 const Vec3D& p3);

		void AddTriangleColored(const Vec3D& p1, const ColorRGB& col1,
								const Vec3D& p2, const ColorRGB& col2,
								const Vec3D& p3, const ColorRGB& col3);

		void AddTriangleNormalColored(const float nv[], const float p1[], const float p2[], const float p3[],
									  const ColorRGB& col1, const ColorRGB& col2, const ColorRGB& col3);

		void AddQuad(const Vec3D& p1, const Vec3D& p2,
					 const Vec3D& p3, const Vec3D& p4);

		void AddQuadColored(const Vec3D& p1, const ColorRGB& col1,
							const Vec3D& p2, const ColorRGB& col2,
							const Vec3D& p3, const ColorRGB& col3,
							const Vec3D& p4, const ColorRGB& col4);

		/*void AddText(const Vec3D& pt, const std::string& text, const char* font);*/
		enum eFont {
			fixed13,
			fixed15,
			times10,
			times24,
			helvetica10,
			helvetica12,
			helvetica18,
		};
		void AddText(const float x, const float y, const float z,
			eFont font, const std::string& stext);

		Object*& AddNewObject(const std::string& sName, int* pIdx = NULL);
		Object*  ResetSubObject(const unsigned int idx_);
		void     ClearSubObject(const unsigned int idx_);

		void SwitcObject(const int idx, const bool hide);

		/*void AddToobject(const int idx);*/

		void SetDrawColor(const ColorRGB& col);

		void SetPointSize(const float size);

		void SetLineWidth(const float size);

		
				
		void DeleteDisplayLists();

		void DumpCharacteristics(std::ostream& pOstream,
                                 const std::string& pIndentation,
                                 const Matrix<float>& pTransformation);

		const BBox3D& GetBBox3D() const;

		void Render(RenderParams&  pParams, bool hide=false);

		const std::string& GetName() const { return aName; }


private:
	// Block the use of those
	Object(const Object&);
	Object& operator=(const Object&);

	typedef  std::vector<std::string>                 Commands;
	typedef  std::map<std::string, Object*>           IndexNamedObjects;
	typedef  std::vector<PrimitiveAccumulator*>       PrimitiveAccumulators;
	typedef  std::vector< std::pair<bool,Object*> >   SubObjects;
	typedef  std::vector<VtxAccumulator*>             VtxAccumulators;
	typedef  std::vector<VtxPrimitiveAccumulator*>    VtxPrimitiveAccumulators;


	void ConstructDisplayList(RenderParams& pParams);

	enum enum_exec_cmd {
		eExecPrimAccum,
		eExecVertexPrimAccum,
	};

	void executeCommand(const std::string& pCommand, RenderParams& pParams) const;

	PrimitiveAccumulator& GetCurrentPrimitiveAccumulator();
public:
	VtxAccumulator& GetCurrentVtxAccumulator();
private:
	VtxPrimitiveAccumulator& GetCurrentVtxPrimitiveAccumulator();

	GLuint& GetGLDisplayList(const RenderParams& pParams);


	mutable BBox3D	             aBoundingBox;
	Commands                       aCommands;
	bool                           aFrozen;
	GLuint                         aGLDisplayListBoundingBox;
	GLuint                         aGLDisplayListFast;
	GLuint                         aGLDisplayListFull;
	GLuint						 aGLvbo[2];
	std::string                    aName;
	bool                           aNewPrimitiveAccumulatorNeeded;
	bool                           aNewVtxAccumulatorNeeded;
	bool                           aNewVtxPrimitiveAccumulatorNeeded;
	PrimitiveAccumulators          aPrimitiveAccumulators;
	int                            aRawModeArrowTipNbPolygons;
	float                          aRawModeArrowTipProportion;
	SubObjects                     aSubObjects;
	VtxAccumulators                aVtxAccumulators;
	VtxPrimitiveAccumulators		 aVtxPrimitiveAccumulators;
	int							 aIndex;


};

} // end namespace FemViewer
#endif /* _OBJECT_H_
*/
