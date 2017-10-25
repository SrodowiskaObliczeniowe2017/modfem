#include "defs.h"
#include "../../utils/fv_assert.h"
#include "fv_txt_utls.h"
#include "Log.h"
#include "ocl.h"
#include "VtxAccumulator.h"
#include "ViewManager.h"
#include "Legend.h"
#include "Light.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <cstdio>
#include <cstring>
#include <memory>
#include <GL/freeglut.h>

#define MAX_LEGEND_ITEMS 5

namespace FemViewer {

GLuint g_GLList;

bool initGLLists()
{
	if (glIsList(g_GLList) == GL_FALSE) {
		//mfp_log_debug("Generating %u GL lists\n",ID_ALL);
		g_GLList = glGenLists(NUM_TOTAL_LISTS);
	}

	return (glIsList(g_GLList) == GL_TRUE);
}

void deleteGLList(const GLuint type, const int size)
{
	if (glIsList(g_GLList+type)) glDeleteLists(g_GLList+type,size);
}

GLuint createGLLists(const GLuint type,const int size)
{
	deleteGLList(type,size);
	return glGenLists(size);
}

template<typename T>
bool OnList(const std::vector<T>& List, const T& val)
{
	for (size_t i(0); i < List.size(); ++i) {
		if (val == List[i]) return true;
	}

	return false;
}

int VtxAccumulator::counter;

void VtxAccumulator::setUniformMatrix(const GLfloat* pData,
		                              const int pType)
{
	static const size_t lSize = 16 * sizeof(GLfloat);
	if (aBufferName[UBO_PROJ] == 0) return;

	glBindBuffer(GL_UNIFORM_BUFFER, aBufferName[UBO_PROJ]);
	glBufferSubData(GL_UNIFORM_BUFFER, pType ? lSize : 0, lSize, (const GLvoid*)pData);
	//glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void VtxAccumulator::setUniformParams()
{
	BaseParams lgData;
	//mfp_debug("In uniform params\n");
	int num = aLegend.PackValuesIntoArray<float>(lgData.iso_values, MAX_ISO_VALUES);
	if (num > 2) {
		if (glIsBuffer(aBufferName[UBO_PARAMS]) == GL_TRUE) {
			glBindBuffer(GL_UNIFORM_BUFFER, 0);
			glDeleteBuffers(1, &aBufferName[UBO_PARAMS]);
		}

		//mfp_debug("Setting data for color map sizeof baseParams %d\n",sizeof(BaseParams));
		lgData.wireframe_col[0] = 1;//aEdgeColor.x;
		lgData.wireframe_col[1] = 1;//aEdgeColor.y;
		lgData.wireframe_col[2] = 1;//aEdgeColor.z;
		lgData.wireframe_col[3] = 1.0f;
		lgData.border_col[3]	= 1.0f;
		lgData.iso_color[0]   	= 0.1f;
		lgData.iso_color[1]		= 0.1f;
		lgData.iso_color[2]		= 0.1f;
		lgData.border_col[3]    = 1.0f;

		lgData.num_breaks = num;

		aBufferName[UBO_PARAMS] = createBuffer(&lgData,  sizeof(BaseParams), GL_UNIFORM_BUFFER, GL_STATIC_DRAW);
		if (aBufferName[UBO_PARAMS] == 0) {
			fprintf(stderr,"Error in creating uniform buffer for legend values\n");
			exit(-1);
		}

		glBindBufferRange(GL_UNIFORM_BUFFER, Shader::IsoValuesBlkIdx, aBufferName[UBO_PARAMS],
						0, sizeof(BaseParams));
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
		//delete[] pValues;
		aUpdateRenderParams = !aUpdateRenderParams;

	}
}

VtxAccumulator::VtxAccumulator()
: viewmgr_ptr(nullptr)
, renderer(nullptr)
, aType(LINEAR)
, aGLinited(false)
, edges_renderer(nullptr)
, faces_renderer(nullptr)
, aProjection(ViewManagerInst().GetCurrentView().getProjectionMatrix())
, aModelView(ViewManagerInst().GetCurrentView().getCameraMatrix())
, aLight(ViewManagerInst().GetSettings()->DirectionalLight)
, aLegend(ViewManagerInst().GetLegend())
, aColors()
, aColorsFrozen(false)
//, aNormals                    ()
, aSimplifiedDirty            (true)
, aNumTriFaces				  (0)
, aVerticesFrozen             (false),
//	aNodes						(),
	aNodeColors					(),
	aNodeNormals				(),
	//aEdIndices(),

//	aUBOId						(0),
//	aVBOUsed					(false),
	aUpdateRenderParams			(false),
	aRenderEdges				(false),
	aRenderIsoLines				(false),
	aRenderTraingleStrips		(false),
	aRenderTraingles			(false),
	aEdgeColor					(0.f, 0.f, 0.f),
	aEdgeThickness				(1.5f),
	//aEdgeCount					(0),
	//aProgramId					(),

	aVertexCounter(0)
{
	viewmgr_ptr = & ViewManagerInst();
	renderer = &VtxAccumulator::renderIdle;
	//mfp_log_debug("ctrt\n");
	memset(aVAOIds,0x0,sizeof(aVAOIds));
	memset(aBufferName,0x0,sizeof(aBufferName));
	//initGLEW();
	++counter;
}

VtxAccumulator::~VtxAccumulator()
{
	//mfp_log_debug("Dtr %d\n",counter);
	reset();
	--counter;

}

bool VtxAccumulator::init(const size_t nVerts, const size_t nEdges,const int nType)
{
	//mfp_log_debug("init\n");
	// just in case
	if (nVerts <= 0 || nEdges <= 0) return false;
	// Reset
	reset();
	// Set type
	if (nType == LINEAR || nType == HIGH_ORDER) aType = nType;
	else aType = LINEAR;
	// clear date
	aVertices.reserve(nVerts);
	// reserve buffer for indices
	aEdIndices.reserve(nEdges * 2);
	// Try init GL stuff
	init_gl();

	return true;
}

bool VtxAccumulator::reset(int type)
{
	//mfp_debug("Reseting for type %d\n",type);
	switch(type)
	{
	default:
		aNodes.clear();
		aNodeColors.clear();
		aNodeNormals.clear();
		aBBox.Reset();
		aVertices.clear();
		aEdIndices.clear();
		//aTrBaseIndices.clear();
		//aQuBaseIndices.clear();
		FV_FREE_PTR(edges_renderer);
		aVertexCounter = 0;
		glDeleteVertexArrays(FV_SIZEOF_ARR(aVAOIds), aVAOIds);
		glDeleteBuffers(FV_SIZEOF_ARR(aBufferName), aBufferName);
	case 1: // data for high order interpolation
		aTFaceCounts.clear();
		aQFaceCounts.clear();
		aCounts.clear();
		aTrVertices.clear();
		aQuVertices.clear();
		glDeleteVertexArrays(1,&aVAOIds[vtxQuad]);
		glDeleteBuffers(1,&aBufferName[VBO_FACES]);
	case 0: // data for linear interpolation
		aTrIndices.clear();
		aQuIndices.clear();
		aNumTriFaces = 0;
		FV_FREE_PTR(faces_renderer);
		glDeleteVertexArrays(1,&aVAOIds[vtxTriangle]);
		glDeleteBuffers(1,&aBufferName[IBO_FACES]);
		cleanUp();
	}

	return true;
}

void VtxAccumulator::update()
{
	if (viewmgr_ptr->GetSettings()->bShadingOn) {
		renderer = (aType == LINEAR) ? &VtxAccumulator::renderTriangles : &VtxAccumulator::renderTriangleStrips;
	} else if (viewmgr_ptr->GetSettings()->bEdgeOn) {
		renderer = &VtxAccumulator::renderEdges;
	} else {
		renderer = &VtxAccumulator::renderIdle;
	}

	if (viewmgr_ptr->GetSettings()->eRenderType == RAYTRACE_GL_CL) {
		renderer = &VtxAccumulator::renderGLCL;
	}
}

void VtxAccumulator::create(/*const int nType*/)
{
	// Init OpenGL stuff
	init_gl();
	// Create UBO for projection matrice
	size_t sizeBytes = 2 * 16 * sizeof(GLfloat);
	aBufferName[UBO_PROJ] = createBuffer(NULL, sizeBytes, GL_UNIFORM_BUFFER, GL_DYNAMIC_DRAW);
	if (aBufferName[UBO_PROJ] == 0) {
		fprintf(stderr,"Error in creating uniform buffer for matrices\n");
		exit(-1);
	}
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	glBindBufferRange(GL_UNIFORM_BUFFER, Shader::ProjectionBlkIdx, aBufferName[UBO_PROJ],
			0, sizeBytes);
	// Create UBO for rendering parameters
	setUniformParams();
	//FV_CHECK_ERROR_GL();
	// Create edges
	createEdges(&aVAOIds[vtxEdge]);
	// Create triangles
	if (aType) createTriangleStrips(&aVAOIds[vtxQuad]);
	else createTriangles(&aVAOIds[vtxTriangle]);
	//FV_CHECK_ERROR_GL();
	float Vcol[] = {1.f, 0.f, 1.f };
	float Ecol[] = {0.f, 1.f, 1.f };
	createText(NUM_VERTICES_LIST_COLORMAP,Vcol);
	createText(NUM_ELEMENTS_LIST_COLORMAP,Ecol);

	createGLCL();
}

void VtxAccumulator::createEdges(GLuint *pVAO)
{
	// Load vertices into GPU memory
	assert(aVertices.size() > 0);
	size_t data_size = aVertices.size() * sizeof(aVertices[0]);
	aBufferName[VBO_VERTICES] = createBuffer(aVertices.data(),data_size,GL_ARRAY_BUFFER,GL_STATIC_DRAW);
	//mfp_debug("aBufferName[VBO_VERTICES] = %u\n",aBufferName[VBO_VERTICES]);
	if (aBufferName[VBO_VERTICES] == 0) {
		fprintf(stderr,"Error in creating buffer for vertices\n");
		throw "Can't allocate OpenGL Buffer for vertices!";
	}

	// Create and load edge index data
	assert(aEdIndices.size() > 0);
	//mfp_debug("liczba indeksow: %u i verteksów: %u\n",aEdIndices.size(),aVertices.size());
	data_size = aEdIndices.size() * sizeof(aEdIndices[0]);
	aBufferName[IBO_EDGES] = createBuffer(aEdIndices.data(), data_size,GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
	if (aBufferName[IBO_EDGES] == 0) throw "Can't create IBO for edges!";
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	// Generate and bind VAO
	glGenVertexArrays(1, pVAO);
	glBindVertexArray(*pVAO);
	glBindBuffer(GL_ARRAY_BUFFER, aBufferName[VBO_VERTICES]);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(aVertices[0]), (const GLvoid*)0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, aBufferName[IBO_EDGES]);

	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(aIndices[0])*aIndices.size(),&aIndices[0],GL_STATIC_DRAW);
	glBindVertexArray(0);

	FV_CHECK_ERROR_GL();

}
// We just only use VBO for edges
void VtxAccumulator::createTriangles(GLuint *pVAO)
{
	//mfp_debug("Creating tri: %u qua: %u\n",aTrIndices.size()/3,aQuIndices.size()/3);
	// We assume that createEdges just have created VBO buffer
	// so we use it, but we have to create IBO for triangles an load data into it.
	assert(aTrIndices.size() > 0);
	aNumTriFaces = aTrIndices.size() / 3;
	size_t offset = aTrIndices.size() * sizeof(GLuint);
	size_t data_size = aQuIndices.size() * sizeof(GLuint);
	aBufferName[IBO_FACES] = createBuffer(NULL,offset+data_size,GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
	if (aBufferName[IBO_FACES] == 0) throw "Can't create IBO for faces!";
	// Load data
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, offset, (const GLvoid *)&aTrIndices[0]);
	if (!aQuIndices.empty()) {//mfp_debug("qIndices size: %u\n",aQuIndices.size());
		glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, offset, data_size, (const GLvoid *)&aQuIndices[0]);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	FV_CHECK_ERROR_GL();
	// Generate VAO
	glGenVertexArrays(1, pVAO);
	glBindVertexArray(*pVAO);
	//FV_CHECK_ERROR_GL();
	// Bind to position + color
	glBindBuffer(GL_ARRAY_BUFFER, aBufferName[VBO_VERTICES]);
	//FV_CHECK_ERROR_GL();
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	//FV_CHECK_ERROR_GL();
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);
	//FV_CHECK_ERROR_GL();
	//mfp_debug("aBufferName[VBO_VERTICES] = %u\n",aBufferName[VBO_VERTICES]);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)(4*sizeof(GLfloat)));
	//FV_CHECK_ERROR_GL();
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, aBufferName[IBO_FACES]);
	//FV_CHECK_ERROR_GL();
	glBindVertexArray(0);
	//mfp_debug("After creating triangles\n");
	//FV_CHECK_ERROR_GL();
}


#define BUFFER_OFFSET(offset)	((GLchar *)NULL + (offset))
void VtxAccumulator::createTriangleStrips(GLuint *pVAO)
{
	//mfp_debug("Creating triangle-strips\n");
	// This function is for triangle strips.
	// We use different VBO then for edges, and must load vertices into it
	//assert(aTrVertices.size() > 0);
	size_t offset = aTrVertices.size() * sizeof(aTrVertices[0]);
	size_t data_size = aQuVertices.size() * sizeof(aTrVertices[0]);
	aBufferName[VBO_FACES] = createBuffer(NULL,offset + data_size,GL_ARRAY_BUFFER,GL_STATIC_DRAW);
	if (aBufferName[VBO_FACES] == 0) {
		fprintf(stderr,"Error in creating buffer for vertices\n");
		exit(-1);
	}
	glBufferSubData(GL_ARRAY_BUFFER, 0, offset, (const GLvoid *)&aTrVertices[0]);
	if (!aQuVertices.empty()) {
		glBufferSubData(GL_ARRAY_BUFFER, offset, data_size, (const GLvoid *)&aQuVertices[0]);
	}
	FV_CHECK_ERROR_GL();
	//size_t ntri_in_faces = 0;
	//for (size_t i(0); i < aTFaceCounts.size(); ++i) ntri_in_faces += aTFaceCounts[i] + 2;
	aNumTriFaces = aTFaceCounts.size();
	//for (int i(0);i<aTrLocIndices.size();++i) std::cout << "TrLocInd" << i << " " << aTrLocIndices[i] << std::endl;
	//for (int i(0);i<aQuLocIndices.size();++i) std::cout << "QuLocInd" << i << " " << aQuLocIndices[i] << std::endl;
	//mfp_debug("TFCounts: %u\n",aTFaceCounts.size());
	// Merge two vectors into one
	if (!aQFaceCounts.empty()) {
		aTFaceCounts.insert(aTFaceCounts.end(),aQFaceCounts.begin(),aQFaceCounts.end());
	}
	//mfp_debug("Before creating %u %u\n",aTFaceCounts.size(), aQFaceCounts.size());
	//for(int i(0); i < aQFaceCounts.size(); ++i) std::cout << "aTRF" << i << " = " << aTFaceCounts[i] << std::endl;
	// Now, we must prepare a storage for offsets of strips

	aCounts.reserve(aTFaceCounts.size());
	size_t ptr(0);
	//aCounts.push_back(ptr);
	Counts::iterator it   = aTFaceCounts.begin();//aTFaceCounts.begin();
	Counts::iterator it_e = aTFaceCounts.end();//aTFaceCounts.end();

	while(it != it_e) {
		//mfp_debug("Adding offset = %u\n",ptr);
		aCounts.push_back(ptr);
		ptr += sizeof(GLuint)*(*it);
		++it;
	}
	//mfp_debug("After creating %u\n",aCounts.size());
	// Merge base indices of vertex
	/*if (!aQuBaseIndices.empty()) {
		const int offset = static_cast<int>(aTrVertices.size());
		//mfp_debug("baseQuindicesoffset  = %u\n",aQuBaseIndices.size());
		for (size_t i(0); i < aQuBaseIndices.size(); ++i) {
			//mfp_debug("offet = %d BaseQu%u = %d\n",offset,i,aQuBaseIndices[i]);
			int base_index = offset + aQuBaseIndices[i];
			aTrBaseIndices.push_back(base_index);
		}
	}*/
	if (!aQuIndices.empty()) {
		const size_t lSize = aTrVertices.size();
		for(size_t i(0); i < aQuIndices.size();++i) aQuIndices[i] += lSize;
	}
	//mfp_debug("After merging\n");
	//for (int i(0);i<aCounts.size();++i) std::cout << "aCounts" << i << " " << aCounts[i] << std::endl;
	// Create IBO for faces
	//assert(aTrIndices.size() > 0);
	//2
	//offset = aTrLocIndices.size() * sizeof(aTrLocIndices[0]);
	//data_size = aQuLocIndices.size() * sizeof(aTrLocIndices[0]);
	offset = aTrIndices.size() * sizeof(aTrIndices[0]);
	data_size = aQuIndices.size() * sizeof(aTrIndices[0]);
	//for(int i(0); i < aTrIndices.size(); ++i)
	//	std::cout<< " " << aTrIndices[i];
	//data_size = aQuIndices.size() * sizeof(aTrIndices[0]);
	aBufferName[IBO_FACES] = createBuffer(NULL,offset+data_size,GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
	if (aBufferName[IBO_FACES] == 0) throw "Can't create IBO for faces!";
	// Load data
	//3
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, offset, &aTrIndices[0]);
	if (!aQuIndices.empty()) {
		glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, offset, data_size, &aQuIndices[0]);
	}
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	FV_CHECK_ERROR_GL();
	glGenVertexArrays(1, pVAO);
	glBindVertexArray(*pVAO);
	glBindBuffer(GL_ARRAY_BUFFER, aBufferName[VBO_FACES]);
	//offset = aVertices.size()*sizeof(Vertex);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);
	glEnableVertexAttribArray(1);
	//size_t offset = 4*sizeof(GLfloat);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)(4*sizeof(GLfloat)));
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, aBufferName[IBO_FACES]);
	glBindVertexArray(0);
	//mfp_debug("After creating strips\n");
}

// create a text as a displaylist
void VtxAccumulator::createText(const GLuint type,const float* color)
{
	//mfp_log_debug("createText: %s\n",type == ID_VERTEX ? "VERTEX" : "ELEMENT");
	if (aVertices.empty()) return;
	assert(aVertexCounter <= aVertices.size());

	char buff[32];
	size_t start_id     = (type == NUM_VERTICES_LIST_COLORMAP) ? 0 : aVertexCounter;
	const size_t end_id = (type == NUM_VERTICES_LIST_COLORMAP) ? aVertexCounter : aVertices.size();
	Vertex * p = aVertices.data();
	FV_CHECK_ERROR_GL();
	glNewList(g_GLList+type, GL_COMPILE);
	FV_CHECK_ERROR_GL();
	glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT | GL_ENABLE_BIT); // lighting and color mask
	glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
	glDisable(GL_TEXTURE_2D);

	glColor3fv(color);
	for (size_t iv(start_id); iv < end_id; ++iv) {
		glRasterPos3f(p[iv].position.x,p[iv].position.y,p[iv].position.z);

		const char * pbuff = &buff[0];
		sprintf(buff,"%d",(int)p[iv].info);

		while (*pbuff) {
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13,*pbuff++);
		}

	}
	//mfp_debug("creating GLList: %d\n",type);
	// Restore
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_LIGHTING);
	glPopAttrib();

	glEndList();

	FV_CHECK_ERROR_GL();
}

void VtxAccumulator::createGLCL()
{
	// Specify memory size
	size_t w_screen(viewmgr_ptr->GetWidth());
	size_t h_screen(viewmgr_ptr->GetHeight());
	size_t size_bytes = w_screen * h_screen * sizeof(GLubyte) * 4;
	// Create PBO
	aBufferName[VBO_PIXEL] = createBuffer(NULL, size_bytes, GL_PIXEL_UNPACK_BUFFER, GL_STREAM_DRAW);
	if (aBufferName[VBO_PIXEL] == 0) throw "Can't create PBO for GL-CL operations\n";
	glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
	// Invoke OpenCL
	bindOpenCL(aBufferName[VBO_PIXEL],w_screen,h_screen);


}

void VtxAccumulator::cleanUp()
{
	//mfp_debug("CleanUp\n");
	if (ViewManagerInst().IsGLSupported()) {
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		if (glIsVertexArray(aVAOIds[0]) == GL_TRUE) {
			//mfp_debug("Deleting %d objects in vertex arrays",FV_SIZEOF_ARR(aVAOIds));
			glDeleteVertexArrays(FV_SIZEOF_ARR(aVAOIds), &aVAOIds[0]);
		}
		//mfp_debug("here1\n");
		if (glIsBuffer(aBufferName[0]) == GL_TRUE) {
			//mfp_debug("Deleting %d objects in vertex buffers",FV_SIZEOF_ARR(aBufferName));
			glDeleteBuffers(FV_SIZEOF_ARR(aBufferName), &aBufferName[0]);
		}
		//mfp_debug("here\n");
		memset(aVAOIds,0x0,sizeof(aVAOIds));
		memset(aBufferName,0x0,sizeof(aBufferName));

		g_GLList = createGLLists(NUM_VERTICES_LIST_COLORMAP,2);
	}


	//aVBOUsed = false;
	aBBox.Reset();
	//FV_FREE_PTR(edges_renderer);
	//FV_FREE_PTR(faces_renderer);
}



void VtxAccumulator::addColor(const Vec3D& pColor)
{
  aColors.push_back(pColor);

  aSimplifiedDirty = true;
}

//void VtxAccumulator::addTex(const Tex2D& pTex)
//{
//  aTexCoords.push_back(pTex);
//}

void VtxAccumulator::addVertex(const Vertex& pVertex,const int dest)
{
	if (dest == vtxTriangle) aTrVertices.push_back(pVertex);
	else if (dest == vtxQuad) aQuVertices.push_back(pVertex);
	else aVertices.push_back(pVertex);
	aBBox += pVertex.position;
	aSimplifiedDirty = true;
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

#endif // #ifdef  FV_DUMP_MEMORY_USAGE

// Dump in ASCII the caracteristics of the VtxAccumulator
void VtxAccumulator::dumpCharacteristics(std::ostream&       pOstream,
                                            const std::string&  pIndentation,
                                            const Matrix<float>&    pTransformation) const
{
  pOstream << pIndentation << "VtxAccumulator " << std::endl;
  std::string lIndentation = pIndentation + "  ";

  if (aVerticesFrozen) {
    pOstream << lIndentation << "Vertices frozen = true" << std::endl;
  }
  else {
    pOstream << lIndentation << "Vertices frozen = false" << std::endl;
  }

  if (aColorsFrozen) {
    pOstream << lIndentation << "Colors frozen   = true" << std::endl;
  }
  else {
    pOstream << lIndentation << "Colors frozen   = false" << std::endl;
  }

 /* if (aTexCoordsFrozen) {
    pOstream << lIndentation << "TexCoords frozen   = true" << std::endl;
  }
  else {
    pOstream << lIndentation << "TexCoords frozen   = false" << std::endl;
  }*/

#ifdef FV_DUMP_MEMORY_USAGE
  {
    pOstream << lIndentation << "Memory used by the VtxAccumulator = " << sizeof(*this) << std::endl;

    pOstream << lIndentation << "Memory used by aColors    = " << getStringSizeAndCapacity(aColors  ) << std::endl;
    pOstream << lIndentation << "Memory used by aNormals   = " << getStringSizeAndCapacity(aNormals ) << std::endl;
    pOstream << lIndentation << "Memory used by aVertices  = " << getStringSizeAndCapacity(aVertices) << std::endl;
    pOstream << lIndentation << "Memory used by aTexCoords = " << getStringSizeAndCapacity(aTexCoords) << std::endl;

  }
#endif // #ifdef GLV_DUMP_MEMORY_USAGE

  pOstream << lIndentation << "Number of vertex  = " << aVertices  .size() << std::endl;
  pOstream << lIndentation << "Number of color_v = " << aColors    .size() << std::endl;
}

// Tell the VtxAccumulator that no more colors
// will be added to *this.
// Returns false if the number of colors is different
// from the number of vertices (this i bad and should
// be used to generate an error message
bool VtxAccumulator::freezeColors()
{
  FV_ASSERT(aVerticesFrozen);
  aColorsFrozen = true;

  const bool lColorsOk = (aColors.size() == aVertices.size());

  if (!lColorsOk) {
    // Remove the colors since they're not usable
    std::vector<Vec3D> lTmp;
    aColors.swap(lTmp);
  }

  return lColorsOk;
}

// Tell the VtxAccumulator that no more colors
// will be added to *this.
// Returns false if the number of colors is different
// from the number of vertices (this i bad and should
// be used to generate an error message
//bool VtxAccumulator::freezeTexCoords()
//{
//  FV_ASSERT(aVerticesFrozen);
//  aTexCoordsFrozen = true;
//
//  const bool lCoordsOk = (aTexCoords.size() == aVertices.size());
//
//  if (!lCoordsOk) {
//    // Remove the coords since they're not usable
//    std::vector<Tex2D> lTmp;
//    aTexCoords.swap(lTmp);
//  }
//
//  return lCoordsOk;
//}

// Tell the VtxAccumulator that no more vertices
// will be added to *this.
void VtxAccumulator::freezeVertices()
{
  FV_ASSERT(!aColorsFrozen);
  aVerticesFrozen = true;
}

// Return a reference to the internal data
const VtxAccumulator::Colors& VtxAccumulator::getColors() const
{
  return aColors;
}

// Return a reference to the internal data
//const VtxAccumulator::TexCoords& VtxAccumulator::getTexCoords() const
//{
//  return aTexCoords;
//}

// Return a reference to the internal data
const VtxAccumulator::Vertices& VtxAccumulator::getVertices() const
{
  return aVertices;
}

// Return a reference to the internal data
//  non-const because of it's computed into the VertexedPrimitiveAccumulators
VtxAccumulator::Normals& VtxAccumulator::getNormals()
{
  return aNormals;
}

// Return the state of the frozen flag; used for GLV_ASSERTs
bool VtxAccumulator::getVerticesFrozen() const
{
  return aVerticesFrozen;
}

// Return the state of the frozen flag; used for GLV_ASSERTs
bool VtxAccumulator::getColorsFrozen() const
{
  return aColorsFrozen;
}

void VtxAccumulator::init_gl(bool quiet)
{
	//mfp_debug("init_gl\n");
	if (!viewmgr_ptr->IsGLSupported()) return;
	if (!aGLinited) {
		// Init shader program for edges
		//mfp_debug("before create edge shader\n");
		edges_renderer = new EdgeShader();
		bool res = edges_renderer->Init();
		// Init shader for faces
		// 1 - DG; 0 - STD
		if (aType) {
			renderer = &VtxAccumulator::renderTriangleStrips;
			faces_renderer = new TriStripsVGFShader();
		}
		else {
			renderer = &VtxAccumulator::renderTriangles;
			faces_renderer = new TriangleShader();
		}
		bool res2 = faces_renderer->Init();
		if (viewmgr_ptr->GetSettings()->eRenderType == RAYTRACE_GL_CL) {
			renderer = &VtxAccumulator::renderGLCL;
		}
		// Generate list for text
		//mfp_debug("Before checing GLLists\n");
		bool res3 = initGLLists();
	}

}

void VtxAccumulator::UpdateRenderParams() const
{
	size_t num = aLegend.GetColors().size();
	Light llight(aLight);
	llight.Position() = aModelView * aLight.Position();
	RenderParams lgData;
	double start = aLegend.GetColors()[0].value;
	double delta = aLegend.GetColors()[1].value - start;

	//mfp_log_debug("Loading legend data: %f %f %d size of l",start,delta,num - 2);
	// Send yo GPU memory
	glBindBuffer(GL_UNIFORM_BUFFER, aBufferName[UBO_PARAMS]);                    // activate vbo id to use
	glBufferData(GL_UNIFORM_BUFFER, sizeof(BaseParams), &lgData.sShaderParams, GL_STATIC_DRAW); // upload data to video card
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	aUpdateRenderParams = !aUpdateRenderParams;
}

void VtxAccumulator::render(const RenderParams& pParams)
{
	setUniformMatrix(aModelView.matrix.data(),0);
	setUniformMatrix(aProjection.matrix.data(),1);
	(this->*renderer)(pParams);
	// If no mouse movement
	if (!viewmgr_ptr->GetMouseMode()) {
		// Try render vertex id
		if (viewmgr_ptr->GetSettings()->bShowNumVertices) {
			//mfp_debug("calling vertex: %u",sizeof(Vertex));
			if (glIsList(g_GLList+NUM_VERTICES_LIST_COLORMAP)) glCallList(g_GLList+NUM_VERTICES_LIST_COLORMAP);
		}
		// Try render element id
		if (viewmgr_ptr->GetSettings()->bShowNumElems) {
			if (glIsList(g_GLList+NUM_ELEMENTS_LIST_COLORMAP)) glCallList(g_GLList+NUM_ELEMENTS_LIST_COLORMAP);
		}
		FV_CHECK_ERROR_GL();
	}
}

void VtxAccumulator::renderTriangleStrips(const RenderParams& pParams)
{
	// Specify the number of triangle faces
	const GLsizei count = aCounts.size();
	const int index_begin = 0;

	//GLint baseIndex = 56;
	//mfp_debug("Render strips %u\n",count);
	Light llight(aLight);
	llight.Position() = aModelView * aLight.Position();
	//FV_CHECK_ERROR_GL();
	//glDisable(GL_CULL_FACE);
	faces_renderer->Enable();
	glBindVertexArray(aVAOIds[vtxQuad]);

	//mfp_debug("aVAOIds = %u\n",aVAOIds[vtxQuad]);
	//mfp_debug("nTriangles = %u\n",aNumTriFaces);
	faces_renderer->SetNumOfTriangles(aNumTriFaces /*1000*/);
	faces_renderer->SetLight(aLight);
	faces_renderer->EdgesOn(viewmgr_ptr->GetSettings()->bEdgeOn);
	faces_renderer->IsoLinesOn(viewmgr_ptr->GetSettings()->bIsovalueLineOn);
	//glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT,(const GLvoid*)0);
	//glDrawElements(GL_TRIANGLE_STRIP, 63, GL_UNSIGNED_INT,(const GLvoid*)0);
	//mfp_debug("counts = %d offset = %d baseIndex = %d\n",aTFaceCounts[index_begin],aCounts[index_begin],aTrBaseIndices[index_begin]);
//	glMultiDrawElementsBaseVertex(GL_TRIANGLE_STRIP,&aTFaceCounts[index_begin],
//			GL_UNSIGNED_INT,
//			(GLvoid **)&aCounts[index_begin], count,
//			//(GLvoid **)&my, count,
//			(GLint *)&aTrBaseIndices[index_begin]
//			//&baseIndex
//			                         );
	glMultiDrawElements(GL_TRIANGLE_STRIP,&aTFaceCounts[index_begin],
				GL_UNSIGNED_INT,
				(const GLvoid **)&aCounts[index_begin], count
				//(GLvoid **)&my, count,
				//(GLint *)&aTrBaseIndices[index_begin]
				//&baseIndex
				);

	FV_CHECK_ERROR_GL();
	glBindVertexArray(0);
	glUseProgram(0);
}

void VtxAccumulator::renderTriangles(const RenderParams& pPr)
{
	// Specify the number of triangle faces
	const GLsizei count = aTrIndices.size() + aQuIndices.size();
	faces_renderer->Enable();
	glBindVertexArray(aVAOIds[vtxTriangle]);
	Light llight(aLight);
	llight.Position() = aModelView * aLight.Position();
	//FV_CHECK_ERROR_GL();
	//mfp_debug("nTriangles = %u %d\n",aTrIndices.size(),count);
	faces_renderer->SetNumOfTriangles(aNumTriFaces);
	faces_renderer->SetLight(aLight);
	faces_renderer->EdgesOn(viewmgr_ptr->GetSettings()->bEdgeOn);
	faces_renderer->IsoLinesOn(viewmgr_ptr->GetSettings()->bIsovalueLineOn);
	//FV_CHECK_ERROR_GL();

	glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT,(const GLvoid *)0);
	FV_CHECK_ERROR_GL();

	glBindVertexArray(0);
	glUseProgram(0);
	//FV_CHECK_ERROR_GL();
}

void VtxAccumulator::renderEdges(const RenderParams& pParams)
{
	//mfp_debug("Render edges %d %d\n",aEdIndices.size(),aEdgeCount);
	//glColor3f(1,1,0);
	edges_renderer->Enable();
	//glUseProgram(aProgamId[0]);

	glBindVertexArray(aVAOIds[vtxEdge]);
	//std::cout << aProjection << std::endl;
	//std::cout << aModelView << std::endl;
	//aEdgeShader.SetMatrixMP(aProjection);
	//aEdgeShader.SetMatrixMV(aModelView);
	//aEdgeShader.SetColor(aEdgeColor.v);

	glDrawElements(GL_LINES,aEdIndices.size(),GL_UNSIGNED_INT,(const GLvoid*)0);
	//glDrawArrays(GL_TRIANGLES, 0 ,3);
	FV_CHECK_ERROR_GL();
	//glDisableVertexAttribArray(0);
	glBindVertexArray(0);
	glUseProgram(0);

}

void VtxAccumulator::renderIdle(const RenderParams& pParams)
{
	;
}

void VtxAccumulator::renderGLCL(const RenderParams& pPr)
{
	;
}

} // end namespace FemViewer
