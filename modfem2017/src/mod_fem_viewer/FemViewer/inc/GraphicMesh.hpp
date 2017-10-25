#ifndef GRAPHIC_MESH_
#define GRAPHIC_MESH_

#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <array>
#include <string.h>
#include "defs.h"
#include "Enums.h"
#include "RenderParams.h"
#include "Matrix.h"
#include "MathHelper.h"
#include "fv_inc.h"
#include "fv_txt_utls.h"
#include "Mesh.h"
#include "Field.h"
#include "Light.h"
#include "Shader.h"
#include"ViewManager.h"
#include "GLRenderBase.h"

namespace FemViewer {

enum {
	EDGE_RENDERER,
	TRIANGLE_RENDERER,
	TRAINGLE_STRIP_RENDERER,
};



class GraphicsSettings;

template<>
class RenderObject<TypeOfRenderer::Wireframe> : public GLRenderBase {

public:
	RenderObject();
	~RenderObject();

	bool Create(void* arg_ptr);
	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings);
	virtual void DoRender(RenderParams* rparams);

private:
	RenderObject(const RenderObject&);
	RenderObject& operator=(const RenderObject&);

	GLShaderProgram m_program;
	GLuint m_dispLists[2];

	bool m_dirty;
};

template<>
class RenderObject<TypeOfRenderer::WireframeSlice> : public GLRenderBase {

public:
	RenderObject(const std::vector<CutPlane>& cuts);
	~RenderObject();

	bool Create(void* arg_ptr);
	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings);
	virtual void DoRender(RenderParams* rparams);

private:
	GLShaderProgram m_program;
	const std::vector<CutPlane>& m_cut_planes;
	std::vector<GLsizei> m_cutted_counts;
	std::vector<int64_t> m_cutted_addresses;
	GLuint m_dispLists[2];

	std::vector<Vertex> m_vertices;
	std::vector<Origin> m_origins;

	bool m_dirty;

};


template<>
class RenderObject<TypeOfRenderer::ColorMap> : public GLRenderBase {

public:
	RenderObject();
	~RenderObject();

	bool Create(void* arg_ptr);
	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings);
	virtual void DoRender(RenderParams* rparams);

private:
	GLShaderProgram m_program;
	std::vector<GLsizei> m_counts;
	std::vector<int64_t> m_addresses;
	GLuint m_dispLists[2];

	bool m_dirty;
};



template<>
class RenderObject<TypeOfRenderer::ColorMapSlice> : public GLRenderBase {

public:
	RenderObject(const std::vector<CutPlane>& cuts);
	~RenderObject();


	bool Create(void* arg_ptr);
	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings);
	virtual void DoRender(RenderParams* rparams);

private:
	GLShaderProgram m_program;
	const std::vector<CutPlane>& m_cut_planes;
	GLuint m_dispLists[2];
    bool m_dirty;

};

template<>
class RenderObject<TypeOfRenderer::ColorMapStd> : public GLRenderBase {

public:
	RenderObject();
	~RenderObject();

	bool Create(void* arg_ptr);
	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings);
	virtual void DoRender(RenderParams* rparams);
private:
	GLShaderProgram m_program;
	int m_num_tri_faces;
	GLuint m_dispLists[2];
	bool m_dirty;

};


template<>
class RenderObject<TypeOfRenderer::ColorMapStdSlice> : public GLRenderBase {

public:
	RenderObject(const std::vector<CutPlane>& cuts);
	~RenderObject();

	bool Create(void* arg_ptr);
	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings);
	virtual void DoRender(RenderParams* rparams);
private:
	GLShaderProgram m_program;
	const std::vector<CutPlane>& m_cut_planes;
	GLuint m_dispLists[2];
	bool m_dirty;
};

enum RenderCore {
	CORE_GL = 0,
	CORE_OPENCL
};

// Forward declaration
template< typename TCore >
class RenderManager;

template< typename TCore>
RenderManager<TCore>& RenderManagerInst(void)
{
	static RenderManager<TCore> render_mgr;
	return render_mgr;
}

template< typename TCore >
using RenderCaller = void (RenderManager<TCore>::*)(RenderParams*);

// Class for graphic mesh
template< typename TCore >
class RenderManager
{
	friend RenderManager<TCore>& RenderManagerInst<>(void);

	typedef std::unique_ptr<GLRenderBase> Hnd2Render;
	typedef std::array<Hnd2Render,NUM_DRAWS> arRenders;
	typedef arRenders::iterator itRenders;

	typedef enum {
		mesh_items = 0,
		field_items,
		idlle,
		num_items,
	} render_type_index;

public:
	~RenderManager(void);

private:
	// Ctr
	RenderManager();
	// Not allowed
	RenderManager(const RenderManager&);
	RenderManager& operator=(const RenderManager&);

	//TCore& m_core;
	arRenders m_renders;
	GLuint m_uboID;
	int m_currMethodId;
	Light& m_light;
	RenderCaller<TCore> m_renderCall;
	bool m_needSetup;

public:
	void Setup();
	void Update(BaseParams* pParams);
	void SetUniformParams(BaseParams* params);
	void SetUniformMatrix(GLfloat* matrix,const int type);
	void SetMatrix(GLfloat matrix[32]);
	void AddObject2Render(int objType,const std::vector<CutPlane>* pCutPlanes=nullptr);
	bool InitFromMeshData(Mesh* pMesh,std::vector<CutPlane>* pCutPlanes=nullptr);
	bool InitFromFieldData(Field* pField,std::vector<CutPlane>* pCutPlanes=nullptr);

	bool IsObjectInitialized(std::size_t objType) const {
		return (m_renders.at(objType) ? m_renders.at(objType)->IsEnabled() : false);
	}
	void SetObjectVisibility(std::size_t objType, bool flag) { m_renders.at(objType)->SetVisible(true); }

	void Render(RenderParams* pParams);
private:
	void RenderColorMap(RenderParams* pParams);
	void RenderWireframe(RenderParams* pParams);
	void RenderIdle(RenderParams* pParams);
public:
	void Clear();

};


template< typename TCore >
RenderManager<TCore>::RenderManager()
: m_renders()
, m_uboID(0)
, m_currMethodId(field_items)
, m_light(ViewManagerInst().GetSettings()->DirectionalLight)
, m_renderCall(nullptr)
, m_needSetup(true)
{
	// Create uniform buffer for parameters
	//mfp_log_debug("RenderManager ctr");
	size_t size_bytes = sizeof(BaseParams);
	m_uboID = createBuffer(NULL, size_bytes, GL_UNIFORM_BUFFER, GL_STATIC_DRAW);
	if (m_uboID == 0) throw fv_exception("Can't allocate OpenGL memory for uniform parameters!");
	glBindBufferRange(GL_UNIFORM_BUFFER, GLShaderProgram::ParametersBlockIndex, m_uboID, 0, size_bytes);
}

template< typename TCore >
RenderManager<TCore>::~RenderManager()
{
	//mfp_log_debug("RenderManager dtr");
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	if (m_uboID != 0) {
		glDeleteBuffers(1,&m_uboID);
	}

}

template< typename TCore >
void RenderManager<TCore>::Setup()
{
	m_needSetup = true;
	GraphicsSettings* settings_ptr = ViewManagerInst().GetSettings();
	for(itRenders it=m_renders.begin();it != m_renders.end();++it)
	{
		if (*it) (*it)->Setup(settings_ptr);
	}
}

template< typename TCore >
void RenderManager<TCore>::Update(BaseParams* pParams)
{
	//mfp_log_debug("in Update\n");
	// Update parameters
	if (pParams == nullptr) return;

	ViewManager& vm = ViewManagerInst();
	SetMatrix(pParams->proj);

	// Update parameters
	pParams->edges_on = vm.GetSettings()->bEdgeOn;
	pParams->isolines_on = vm.GetSettings()->bIsovalueLineOn;
	//std::cout << "ISOLINES STATUS: " << pParams->sShaderParams.isolines_on << std::endl;
	//pParams->num_tri_faces = num_triangles;
	pParams->num_breaks = vm.GetLegend().PackValuesIntoArray<float>(
			pParams->iso_values,MAX_ISO_VALUES);
	//std::cout << "Num breaks: " << pParams->sShaderParams.num_breaks << std::endl;

	// Setup light parameters
	// Diffuse
	int i;
	for (i=0;i<3;i++) pParams->light_intensity[i] = m_light.DiffuseIntensity()[i];
	pParams->light_intensity[i] = 1.0f;
	// Ambient
	for (i=0;i<3;i++) pParams->light_ambient[i] = m_light.AmbientIntensity()[i];
	pParams->light_ambient[i] = 1.0f;

	// Upload parameters to GPU
	glBindBuffer(GL_UNIFORM_BUFFER, m_uboID);
	GLintptr   offset = static_cast<GLintptr>(32*sizeof(GLfloat));
	GLsizeiptr size   = static_cast<GLsizeiptr>(sizeof(BaseParams)-offset);
	// Load to second part of basic parameters omit the matrices
	glBufferSubData(GL_UNIFORM_BUFFER, offset, size,reinterpret_cast<const GLvoid*>(pParams->wireframe_col));
	glBindBuffer(GL_UNIFORM_BUFFER, 0);

	if (vm.GetSettings()->bShadingOn) {
		m_renderCall = &RenderManager::RenderColorMap;
	}
	else if (vm.GetSettings()->bEdgeOn) {
		m_renderCall = &RenderManager::RenderWireframe;
	}
	else {
		m_renderCall = &RenderManager::RenderIdle;
	}

	FV_CHECK_ERROR_GL();
}

template< typename TCore >
void RenderManager<TCore>::SetUniformParams(BaseParams* params)
{
//	//mfp_debug("In uniform params\n");
	//if (!params) return;
//
	//assert(m_uboID!=0);
	//glBindBuffer(1, m_uboID);
//	glBindBufferRange(GL_UNIFORM_BUFFER, Shader::IsoValuesBlkIdx, m_VBO[UBO_PARAMS],
//					0, sizeof(BaseParams));
//	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	//return false;
}

template< typename TCore >
void RenderManager<TCore>::SetUniformMatrix(GLfloat* pData,const int type)
{
//	if (m_mID == 0) return;
//
//    glGetFloatv(type,pData);
//
//	static const int offset = 16*sizeof(GLfloat);
//	glBindBuffer(GL_UNIFORM_BUFFER, m_mID);
//	glBufferSubData(GL_UNIFORM_BUFFER, (type == GL_PROJECTION_MATRIX ? 0 : offset), offset, reinterpret_cast<const GLvoid*>(pData));
//	glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

template< typename TCore >
void RenderManager<TCore>::SetMatrix(GLfloat matrices[32])
{
	if (m_uboID == 0) return;

	glGetFloatv(GL_PROJECTION_MATRIX, &matrices[0]);
	glGetFloatv(GL_MODELVIEW_MATRIX, &matrices[16]);

	glBindBuffer(GL_UNIFORM_BUFFER, m_uboID);
	glBufferSubData(GL_UNIFORM_BUFFER, 0, 32*sizeof(GLfloat), reinterpret_cast<const GLvoid*>(matrices));
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
}



template< typename TCore >
bool RenderManager<TCore>::InitFromMeshData(Mesh* pMesh,std::vector<CutPlane>* pCutPlanes)
{
	//mfp_debug("Initialization from mesh data\n");
	if (!pMesh) return false;

	// Specify a type
	int type = (pCutPlanes != nullptr ? WireframeSlice : Wireframe);

	// Add object
	AddObject2Render(type,pCutPlanes);

	// Initialize a data
	//mfp_log_debug("Before create wireframe mesh\n");
	//m_graphicsObjects.back()->Create((void*)pMesh);
	if (m_renders.at(type)->IsEnabled()) m_renders.at(type)->Clear();
	m_renders.at(type)->Create(pMesh);
	m_renders.at(type)->SetVisible(true);
	return true;
}


template< typename TCore >
bool RenderManager<TCore>::InitFromFieldData(Field* pField,std::vector<CutPlane>* pCutPlanes)
{
	//mfp_debug("Initialization from field data\n");
	if (!pField) return false;

	// Base type
	int type = (pField->GetApproximationType()  == FieldDG ? ColorMap : ColorMapStd);
	int index = COLORMAP_GL;
	// Upgrade to slice type if necessary
	if (pCutPlanes) { type++; index = COLORMAP_CUTS_GL; }

	// Add object if there is not any
	AddObject2Render(type,pCutPlanes);

	// Clear if it was already initiated
	//mfp_log_debug("Before create in InitFromField\n");
	if (m_renders[index]->IsEnabled()) m_renders[index]->Clear();
	m_renders[index]->Create(pField);
	m_renders[index]->SetVisible(true);

	return true;
}


template< typename TCore >
void RenderManager<TCore>::Render(RenderParams* pParams)
{
	//return;
	static Matrixf mat;


	//SetUniformMatrix(mat.matrix.data(),GL_PROJECTION_MATRIX);
	//SetUniformMatrix(mat.matrix.data(),GL_MODELVIEW_MATRIX);
	//Update(&pParams->sShaderParams);
	SetMatrix(pParams->sShaderParams.proj);

	(this->*m_renderCall)(pParams);

	glBindVertexArray(0);
	glUseProgram(0);

	//mfp_debug("draw num of vertices\n");
	//if (pParams->bShowNumOfVertices)
	//{
		//m_core.renderGLList(NUM_VERTICES_LIST);
	//}




}

template< typename TCore >
void RenderManager<TCore>::RenderColorMap(RenderParams* pParams)
{
	static Matrixf mat;

	if (m_light.Type() == Light::Camera) {
		pParams->cLight.Position() = mat * m_light.Position();
	}

	if (pParams->bDrawCutted) {
		m_renders[COLORMAP_CUTS_GL]->Render(pParams);
	}
	else {
		//mfp_log_debug("Rendering colormap");
		m_renders[COLORMAP_GL]->Render(pParams);
	}
}

template< typename TCore >
void RenderManager<TCore>::RenderWireframe(RenderParams* pParams)
{
	if (pParams->bDrawCutted) {
		m_renders[WIREFRAME_CUTS_GL]->Render(pParams);
	}
	else {
		m_renders[WIREFRAME_GL]->Render(pParams);
	}
}

template< typename TCore >
void RenderManager<TCore>::RenderIdle(RenderParams* pParams)
{
	;
}


template< typename TCore >
void RenderManager<TCore>::Clear()
{
	//mfp_log_debug("In clear method\n");
	//glDeleteBuffers(1,&m_uboID);

	for (auto& it : m_renders) it.reset();

	m_renderCall = &RenderManager::RenderIdle;
	m_needSetup = false;
}



template< typename TCore >
void RenderManager<TCore>::AddObject2Render(int objType,const std::vector<CutPlane>* cuts_ptr)
{
	int index = -1;
	switch (objType) {
	case Wireframe: index = WIREFRAME_GL; break;
	case WireframeSlice: index = WIREFRAME_CUTS_GL; break;
	case ColorMap:
	case ColorMapStd: index = COLORMAP_GL; break;
	case ColorMapSlice:
	case ColorMapStdSlice: index = COLORMAP_CUTS_GL; break;
	}

	if (index < 0) throw fv_exception("Failed to add unknown type of renderer!");

	if (m_renders.at(index)) {
		//mfp_log_debug("Renderer (%d) already is present\n",objType);
	}
	else {
		switch (objType) {
		case Wireframe:
			//mfp_log_debug("In AddObject2Render: wireframe\n");
			m_renders[WIREFRAME_GL].reset(new RenderObject<Wireframe>()); break;
		case WireframeSlice:
			//mfp_log_debug("In AddObject2Render: wireframeSlice\n");
			m_renders[WIREFRAME_CUTS_GL].reset(new RenderObject<WireframeSlice>(*cuts_ptr)); break;
		case ColorMap:
			//mfp_log_debug("In AddObject2Render: ColorMap\n");
			m_renders[COLORMAP_GL].reset(new RenderObject<ColorMap>()); break;
		case ColorMapStd:
			//mfp_log_debug("In AddObject2Render: ColorMapStd\n");
			m_renders[COLORMAP_GL].reset(new RenderObject<ColorMapStd>()); break;
		case ColorMapSlice:
			//mfp_log_debug("In AddObject2Render: ColorMapSlice\n");
			m_renders[COLORMAP_CUTS_GL].reset(new RenderObject<ColorMapSlice>(*cuts_ptr)); break;
		case ColorMapStdSlice:
			//mfp_log_debug("In AddObject2Render: ColorMapStdSlice\n");
			m_renders[COLORMAP_CUTS_GL].reset(new RenderObject<ColorMapStdSlice>(*cuts_ptr)); break;
		default:
			assert(!"uknown render type");
		}
	}

	//mfp_debug("Here is something\n");
	//if (isOk) m_graphicsObjects.push_back(RenderHandle(tmpPtr));
	//mfp_debug("No, everithibng isd oK\n");
}













} // end namespace FemViewer

#endif /* _GRAPH_ELEMENT_H_ */
