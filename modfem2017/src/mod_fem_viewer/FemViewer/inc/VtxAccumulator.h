#ifndef _VTXACCUMULATOR_H_
#define _VTXACCUMULATOR_H_

#include "fv_inc.h"
#include "fv_assert.h"
#include "Enums.h"
#include "BBox3D.h"
#include "Matrix.h"
#include "RenderParams.h"
#include "MathHelper.h"
#include "Shader.h"
#include<cstddef>	// for size_t
#include <string>
#include <vector>

namespace FemViewer {

extern GLuint g_GLList;
extern bool initGLLists();
extern GLuint createGLLists(const GLuint type,int size);

class Legend;
class ViewManager;






class PrimitiveAccumulator;
class Light;
class VtxAccumulator;
typedef void (VtxAccumulator::*drawMethod)(const RenderParams&);

// Class used to accumulate OpenGL based on a mesh
// (vertices) .  Contains the base data for VertexedPrimitiveAccumulator
//  : vertexes, colors and normals.
// This data may be shared between multiple VertexedPrimitiveAccumulator


class VtxAccumulator
{
  enum {
	UBO_PROJ = 0,
	UBO_PARAMS,
	VBO_VERTICES,
	VBO_FACES,
	VBO_PIXEL,
	IBO_EDGES,
	IBO_FACES,
	ALL_BO
  };

  static int counter;

public:
  void setUniformMatrix(const GLfloat *pData,const int pType = 0);
  void setUniformParams();
  typedef  GLuint index_type_t;
public:
 explicit VtxAccumulator();
  ~VtxAccumulator();
  bool  init(const size_t nVerts,const size_t nEdges, const int type = LINEAR);
  bool  reset(int type = -1);
  void  update();
  void  create(/*const int nType = 0*/);
  void  createVertexBuffer();
  void  createIndexBuffer();
  void  createEdges(GLuint *pVAO);
  void  createTriangles(GLuint *pVA0);
  void  createTriangleStrips(GLuint *pVA0);
  void  createText(const GLuint type,const float *color);
  void  createGLCL();
  void  cleanUp();
  void  addColor           (const Vec3D&     pColor);
  //void  addVertex          (const Vec3D&     pVertex);
  void  addVertex(const Vertex& pVertex,const int dest = vtxEdge);
  bool  isVertexInitialized(const size_t index) const {
	  return !(aVertices.at(index).color.x < 0.0f);
  }
  void  addNode			   (const fvmath::Vec3f& pNode) { 
	  aNodes.push_back(pNode);
	  aBBox += pNode;
  }
  void  addEdge			   (const int s,const int e) {
	  assert(s != e);
	  if (s != e) {
		  aEdIndices.push_back(s);
		  aEdIndices.push_back(e);
		 // aEdgeCount += 2;
	  }
  }

  void addNodeColor(const fvmath::Vec3f& color) {
	  aNodeColors.push_back(color);
  }

  void addNormal2Node(const fvmath::Vec3f& normal) {
	  aNodeNormals.push_back(normal);
  }

  void addIndexOfNode(const unsigned int index,const int type = vtxTriangle) {
	  if (type == vtxTriangle) aTrIndices.push_back(index);
	  else aQuIndices.push_back(index);
  }

  void addLocalIndexOfNode(const index_type_t index,const int type = vtxTriangle) {
	  if (type == vtxTriangle) aTrLocIndices.push_back(index);
	  else aQuLocIndices.push_back(index);
  }

  void addPrimitiveCount(const unsigned int count,const int dest = vtxTriangle) {
	  if (dest) aQFaceCounts.push_back(count);
	  else aTFaceCounts.push_back(count);
  }

//  void addBaseIndex(const ubyte_t index, const int dest = vtxTriangle) {
//	  if (dest) aQuBaseIndices.push_back(index);
//	  else aTrBaseIndices.push_back(index);
//  }

  void dumpCharacteristics(std::ostream& pOstream,
                           const std::string&  pIndentation,
                           const Matrix<float>&  pTransformation) const;


  bool  getVerticesFrozen  () const;
  bool  getColorsFrozen    () const;

  bool  freezeColors       ();
  void  freezeVertices     ();

  typedef  std::size_t	size_type;
  typedef  std::vector<Vec3D>  Colors;
  typedef  std::vector<Vec3D>  Normals;
  //typedef  std::vector<Vec3D>  Vertices;
  typedef  std::vector<Vertex> Vertices;
  typedef  fvmath::Vec3f	   NodeCoords;
  typedef  fvmath::Vec3f	   NodeColor;
  typedef  std::vector<NodeCoords> Nodes;
  typedef  std::vector<NodeColor> NodeColors;
  typedef  std::vector<NodeCoords> NodeNormals;
  typedef  Nodes::size_type     NodesSizeType;
  typedef  std::vector<unsigned int> Indices;
  typedef  std::vector<GLsizei> Counts;
  typedef  std::vector<GLint>	BaseIndices;
  typedef  std::vector<Shader *> Shaders;
  typedef  std::vector<index_type_t> LocalIndices;


  Vertices&   getVertices(const int type = vtxEdge) {
	  if (type == vtxEdge) return aVertices;
	  return (type == vtxTriangle ? aTrVertices : aQuVertices);
  }

  Indices& getIndices(const int type = vtxEdge) {
	  if (type == vtxEdge) return aEdIndices;
	  return (type == vtxTriangle) ? aTrIndices : aQuIndices;
  }

  LocalIndices& getLocalIndices(const int type = vtxTriangle) {
	  if (type) return aQuLocIndices;
	  else return aTrLocIndices;
  }
  NodeColors& getNodeColors() { return aNodeColors; }

  size_t& getVertexCounter() { return aVertexCounter; }
  // Record the number of main nodes
  //void Record() { aNumVerts = aVertices.size(); }
  GLuint* getPixelBufferHandle() { return &aBufferName[VBO_PIXEL]; }

  BBox3D getBBox3D() { return aBBox; }


protected:
  const ViewManager* viewmgr_ptr;
  drawMethod		 renderer;
  int 				 aType;
  Vertices           aVertices;
  Indices			 aEdIndices;
  Vertices			 aTrVertices;
  Indices			 aTrIndices;
  Vertices			 aQuVertices;
  Indices			 aQuIndices;
  bool			     aGLinited;
  GLuint			 aVAOIds[vtxAll];
  GLuint			 aBufferName[ALL_BO];
  Counts			 aTFaceCounts;
  Counts			 aQFaceCounts;
  std::vector<int64_t> aCounts;
  unsigned int	     aNumTriFaces;
  const Matrixf&	 aProjection;
  const Matrixf&	 aModelView;
  const Light&		 aLight;
  const Legend&		 aLegend;
  Colors             aColors;
  bool               aColorsFrozen;
  Normals          aNormals;
  bool               aSimplifiedDirty;


  bool                   aVerticesFrozen;
  Nodes					 aNodes;
  NodeColors			 aNodeColors;
  NodeNormals			 aNodeNormals;


  //bool					 aVBOUsed;

  //BBox3D				 aBBox;
  CVec3f				 aEdgeColor;
  float					 aEdgeThickness;
  //unsigned int 			 aEdgeCount;
  //GLuint 				 aProgramId[2];
  GLuint				 aUniforms[5];
  EdgeShader* edges_renderer;
  TriangleShader* faces_renderer;

  LocalIndices           aTrLocIndices;
  LocalIndices			 aQuLocIndices;
  //BaseIndices			 aTrBaseIndices, aQuBaseIndices;
  //GLuint				 aTxtID[ID_ALL];
  size_type				 aVertexCounter;
  BBox3D				 aBBox;

  public:
    mutable bool  aUpdateRenderParams;
	bool  aRenderEdges;
	bool  aRenderIsoLines;
	bool  aRenderTraingleStrips;
	bool  aRenderTraingles;
	const Colors& getColors() const;
	const Vertices& getVertices() const;
	const Nodes& getNodes() const { return aNodes; }

	//const Indices& getIndices() const { return aIndices; }
	//bool& UseVBO() { return aVBOUsed; }
	Normals& getNormals();
	void render(const RenderParams& pPram);
	//const BBox3D& getBBox3D() const { return aBBox; }
	void setEdgeColor(const Vec3f& color) { aEdgeColor = color; }
	void setEdgeWidth(const float width) { aEdgeThickness = width; }
  private:
	void init_gl(bool quiet = true);
	void UpdateRenderParams() const;

	void renderTriangleStrips(const RenderParams& pPr);
	void renderTriangles(const RenderParams& pPr);
	void renderEdges(const RenderParams& pPr);
	void renderIdle(const RenderParams& pPr);
	void renderGLCL(const RenderParams& pPr);


  private:

    // Block the use of those
    VtxAccumulator(const VtxAccumulator&);
    VtxAccumulator& operator=(const VtxAccumulator&);
};

} // end namespace FemViewer

#endif /* _VTXACCUMULATOR_H_
*/
