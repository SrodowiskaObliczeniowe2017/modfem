/*
 * GraphicMesh.cpp
 *
 *  Created on: 29 lip 2015
 *      Author: pmaciol
 */
#include"GraphicMesh.hpp"
#include"ModelControler.h"

namespace FemViewer {


RenderObject<TypeOfRenderer::Wireframe>::RenderObject()
: GLRenderBase(TypeOfRenderer::Wireframe)
, m_program()
, m_dispLists({0,0})
, m_dirty(false)
{
	//mfp_log_debug("RenderObject ctr type: %d",this->m_type);
	bool res = m_program.AddShader(GLShader::GetPath(EDGE_VERT),GL_VERTEX_SHADER);
	res = res && m_program.AddShader(GLShader::GetPath(EDGE_FRAG), GL_FRAGMENT_SHADER);
	if (!res) throw fv_exception("Can't load shaders for wireframe rendering!");
	else m_program.Complete();
	//mfp_log_debug("After init shader\n");
}


RenderObject<TypeOfRenderer::Wireframe>::~RenderObject()
{
	//mfp_log_debug("RenderObject dtr type: %d",this->m_type);
	//m_coreGL.deleteGLLists(m_dispLists[0],1);
	//m_coreGL.deleteGLLists(m_dispLists[1],1);
}

//void RenderObject<EdgeShader,TypeOfRenderer::Wireframe>::Clear()
//{
//	Renderable::Clear();
//}


bool RenderObject<TypeOfRenderer::Wireframe>::Create(void* arg_ptr)
{
	///mfp_debug("Create Wireframe objecrt: %d\n", (arg_ptr ? 1 : 0));
	if (!arg_ptr) return false;

	std::vector<BaseVertex> vertices;
	std::vector<unsigned> indices;

	Mesh* pMesh = reinterpret_cast<Mesh*>(arg_ptr);
	unsigned int size = pMesh->RenderWireframe(vertices,indices);
	assert(size != 0);

	/*std::cout <<"Mesh vertices:\n";
	for (const auto& v: m_vertices) {
		std::cout << v.position << " id: " << v.info << std::endl;
	}
	std:cout << "Edde indices:\n";
	for (size_t i=0; i < indices.size(); i+=2) {
		std::cout << "Edge" << i << " = {" << indices[i] << ", " << indices[i+1] << std::endl;
	}*/
	//mfp_debug("Number of center of faces: %u\n number of center of elems: %u",size,origins.size()- size);
	//mfp_debug("sizeof Vertex: %u size of origin: %u sizeof Node: %u\n",sizeof(Vertex),sizeof(Origin),sizeof(Node_t));

	m_dirty = (size != 0);
	if (m_dirty) {

		m_vboId = createBuffer(vertices.data(),sizeof(vertices[0])*vertices.size(),GL_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_vboId == 0) return false;

		m_count = indices.size();
		m_iboId = createBuffer(indices.data(),sizeof(unsigned)*indices.size(),GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_iboId == 0) return false;

		glGenVertexArrays(1, &m_vaoId);
		glBindVertexArray(m_vaoId);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboId);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboId);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (const GLvoid*)0);
		glBindVertexArray(0);

		// Create index tags for vertices;
		//mfp_debug("Before creating vertex ids for wireframe");
		m_dispLists[0] = m_coreGL.createGLText(&vertices[0],sizeof(BaseVertex),size,
				ViewManagerInst().GetSettings()->NodeIDColor, NUM_VERTICES_LIST_WIREFRAME,GL_FALSE);
		// Create index tags for faces
		size = pMesh->GetElementOrigins(vertices);
		assert(size != 0);
		//mfp_debug("Before creating element ids for wireframe %d",num_fa_origs);
		// Create index tags for elements
		m_dispLists[1] = m_coreGL.createGLText(&vertices[0],sizeof(BaseVertex),size,
				ViewManagerInst().GetSettings()->ElemIDColor, NUM_ELEMENTS_LIST_WIREFRAME,GL_FALSE);

		assert(m_dispLists[0]!=0&&m_dispLists[1]!=0);
		//mfp_log_debug("disp0=%u disp1=%u",m_dispLists[0],m_dispLists[1]);
	}

	return m_dirty;
}

void RenderObject<TypeOfRenderer::Wireframe>::DoSetup(const GraphicsSettings* settings)
{
	bool result = true;
	if (!m_dirty) {
		result = Create(ModelCtrlInst().GetCurrentMesh());
	}

	/*mfp_log_debug("In DoSetup function of wireframe");
	if (m_dispLists[0] == 0) {
		if (m_vertices.empty()) this->Create(ModelCtrlInst().GetCurrentMesh());
		m_dispLists[0] = m_coreGL.createGLText((BaseVertex*)&m_vertices[0],sizeof(Vertex),m_vertices.size(),NUM_VERTICES_LIST_WIREFRAME,GL_FALSE);
	}
	assert(m_dispLists[0]!=0);

	if (m_dispLists[1] == 0) {
		Mesh* pmesh = ModelCtrlInst().GetCurrentMesh();
	    size_t size = pmesh->GetElementOrigins(m_vertices);
	    assert(size != 0);
		m_dispLists[1] = m_dispLists[1] = m_coreGL.createGLText((BaseVertex*)&m_vertices[0],sizeof(Origin),size,NUM_ELEMENTS_LIST_WIREFRAME,GL_FALSE);
	}
	assert(m_dispLists[1]!=0);

	mfp_log_debug("After DoSetup %d %d\n",m_dispLists[0],m_dispLists[1]);*/
}

void RenderObject<TypeOfRenderer::Wireframe>::DoRender(RenderParams* rparams)
{
	m_program.Enable();
	glBindVertexArray(m_vaoId);
	glDrawElements(GL_LINES, m_count, GL_UNSIGNED_INT, (const GLvoid*)0);

	glBindVertexArray(0);
	glUseProgram(0);

	if (!rparams->eMouseMode) {
		if (rparams->bShowNumOfVertices) {
			//mfp_log_debug("drawing num vertices with %d list",m_dispLists[0]);
			glCallList(m_dispLists[0]);
		}
		if (rparams->bShowNumOfElements) {
			//mfp_log_debug("drawing num elements with %d list",m_dispLists[1]);
			glCallList(m_dispLists[1]);
		}
	}
}



RenderObject<TypeOfRenderer::WireframeSlice>::RenderObject(const std::vector<CutPlane>& cuts)
: GLRenderBase(TypeOfRenderer::WireframeSlice)
, m_program()
, m_cut_planes(cuts)
, m_cutted_counts()
, m_cutted_addresses()
, m_vertices()
, m_origins()
, m_dirty(false)
{
	//mfp_log_debug("RenderObject ctr type: %d",this->m_type);
	bool res = m_program.AddShader(GLShader::GetPath(EDGE_VERT),GL_VERTEX_SHADER);
	res = res && m_program.AddShader(GLShader::GetPath(EDGE_FRAG), GL_FRAGMENT_SHADER);
	if (!res) throw fv_exception("Can't load shaders for wireframe cuts rendering!");
	else m_program.Complete();
}

RenderObject<TypeOfRenderer::WireframeSlice>::~RenderObject()
{
	//mfp_log_debug("RenderObject dtr type: %d",this->m_type);
	//m_coreGL.deleteGLLists(m_dispLists[0],1);
	//m_coreGL.deleteGLLists(m_dispLists[1],1);
}

bool RenderObject<TypeOfRenderer::WireframeSlice>::Create(void* arg_ptr)
{
	//mfp_debug("Create slice\n");
	if (!arg_ptr) return false;

    Mesh* pMesh = reinterpret_cast<Mesh*>(arg_ptr);
	//std::vector<Origin> origins;
	//std::vector<Vertex> vertices;
	std::vector<unsigned> indices;

	m_cutted_counts.clear();
	unsigned int num_counts = pMesh->RenderWireframeCuts(m_cut_planes,m_vertices,m_origins,indices,m_cutted_counts);
	//return false;
	//printf("inideices %u counts: %u origins: %u, vertices: %u \n",
	//			indices.size(),m_cutted_counts.size(),m_origins.size(),m_vertices.size());

	size_t p(0);
	std::vector<int>::const_iterator id_p = m_cutted_counts.begin();
	//std::cout << "Adresses:\n";
	m_cutted_addresses.reserve(indices.size());
	m_cutted_addresses.clear();

	while (id_p != m_cutted_counts.end()) {
		m_cutted_addresses.push_back(p);
		//std::cout << "\t" << p << std::endl;
		p += sizeof(unsigned int)*(*id_p);
		++id_p;
	}

	m_count = static_cast<GLsizei>(m_cutted_counts.size());

	if (m_vertices.empty()) return false;

	m_vboId = createBuffer(m_vertices.data(),sizeof(Vertex)*m_vertices.size(),GL_ARRAY_BUFFER,GL_STATIC_DRAW);
	if (m_vboId == 0) return false;

	m_iboId= createBuffer(indices.data(),sizeof(unsigned)*indices.size(),GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
	if (m_iboId == 0) return false;


	glGenVertexArrays(1, &m_vaoId);
	glBindVertexArray(m_vaoId);
	glBindBuffer(GL_ARRAY_BUFFER, m_vboId);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboId);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(m_vertices[0]), (const GLvoid*)0);
	glBindVertexArray(0);

	// Create index tags for cutted vertices;
	m_dispLists[0] = m_coreGL.createGLText((BaseVertex*)&m_vertices[0],sizeof(m_vertices[0]),m_vertices.size(),
			ViewManagerInst().GetSettings()->NodeIDColor,NUM_VERTICES_LIST_CUTTED_WIREFRAME,GL_FALSE);
	// Create index tags for cutted faces
	m_dispLists[1] = m_coreGL.createGLText((BaseVertex*)&m_origins[0],sizeof(Origin),m_origins.size(),
			ViewManagerInst().GetSettings()->ElemIDColor, NUM_ELEMENTS_LIST_CUTTED_WIREFRAME,GL_TRUE);

	assert(m_dispLists[0] != 0 && m_dispLists[1] != 0);

	m_dirty = true;
	return m_dirty;
}

void RenderObject<TypeOfRenderer::WireframeSlice>::DoSetup(const GraphicsSettings* settings)
{
	/*
	if (m_dispLists[0] == 0)
		m_dispLists[0] = m_coreGL.createGLText((BaseVertex*)&m_vertices[0],sizeof(Vertex),m_vertices.size(),NUM_VERTICES_LIST_WIREFRAME,GL_FALSE);
	assert(m_dispLists[0]!=0);

	if (m_dispLists[1] == 0)
		m_dispLists[1] = m_dispLists[1] = m_coreGL.createGLText((BaseVertex*)&m_origins[0],sizeof(Origin),m_origins.size(),NUM_ELEMENTS_LIST_WIREFRAME,GL_FALSE);
	assert(m_dispLists[1]!=0);
	*/
}

void RenderObject<TypeOfRenderer::WireframeSlice>::DoRender(RenderParams* rparams)
{
	m_program.Enable();
	glBindVertexArray(m_vaoId);

	glMultiDrawElements(GL_LINE_LOOP, &m_cutted_counts[0], GL_UNSIGNED_INT,
			reinterpret_cast<const GLvoid**>(&m_cutted_addresses[0]), m_cutted_counts.size());

	glBindVertexArray(0);
	glUseProgram(0);

	if (!rparams->eMouseMode) {
		if (rparams->bShowNumOfVertices && m_dispLists[0]!=0) glCallList(m_dispLists[0]);
		if (rparams->bShowNumOfElements && m_dispLists[1]!=0) glCallList(m_dispLists[1]);
	}
}

/********************** TypeOfRenderer::ColorMap DG*********************/

RenderObject<TypeOfRenderer::ColorMap>::RenderObject()
: GLRenderBase(TypeOfRenderer::ColorMap)
, m_program()
, m_counts()
, m_addresses()
, m_dispLists({0,0})
{
	//mfp_log_debug("RenderObject ctr type: %d",this->m_type);
	const char* path = GLShader::GetPath(TRI_VERT);
	//mfp_log_debug("VERTEX PATH: %s",path);
	bool res = m_program.AddShader(GLShader::GetPath(TRI_VERT),GL_VERTEX_SHADER);
	path = GLShader::GetPath(TRI_FRAG);
	//mfp_log_debug("FRAGMENT PATH: %s",path);
	res = res && m_program.AddShader(GLShader::GetPath(TRI_FRAG), GL_FRAGMENT_SHADER);
	path = GLShader::GetPath(TRISTRIP_GEOM);
	///mfp_log_debug("GEOMETRY PATH: %s",path);
	res = res && m_program.AddShader(GLShader::GetPath(TRISTRIP_GEOM), GL_GEOMETRY_SHADER);
	if (!res) throw fv_exception("Can't load shaders for color map rendering!");
	else m_program.Complete();
}

RenderObject<TypeOfRenderer::ColorMap>::~RenderObject()
{
	//mfp_log_debug("RenderObject dtr type: %d",this->m_type);
	//m_coreGL.deleteGLLists(m_dispLists[0],1);
	//m_coreGL.deleteGLLists(m_dispLists[1],1);
}


bool RenderObject<TypeOfRenderer::ColorMap>::Create(void* arg_ptr)
{
	//mfp_debug("Init from field\n");
	if (!arg_ptr) return false;
	else {
		Field* pField = reinterpret_cast<Field*>(arg_ptr);

		// Temporary containers
		std::vector<Vertex> vertices;
		std::vector<unsigned> indices;

		// Clear member containers
		m_counts.clear();
		m_addresses.clear();

		m_count = pField->RenderDG(vertices,indices,m_counts,m_addresses);
		assert(m_count != 0);
		if (m_count == 0) return false;

		m_vboId = createBuffer(vertices.data(),sizeof(Vertex)*vertices.size(),GL_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_vboId == 0) return false;

		m_iboId = createBuffer(indices.data(),sizeof(unsigned)*indices.size(),GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_iboId == 0) return false;

		glGenVertexArrays(1, &m_vaoId);
		glBindVertexArray(m_vaoId);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboId);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboId);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (const GLvoid*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)(4*sizeof(GLfloat)));
		glBindVertexArray(0);

		m_dirty = true;
	}


	if (m_dirty) {
		std::vector<BaseVertex> origins, nodes;
		unsigned size = ModelCtrlInst().GetCurrentMesh()->GetOriginsOfBoundaryFaces(origins,nodes);
		//std::cout << "Boundary vertices size: " << size << "\n";
		//for (const auto& v: nodes) std::cout << v.position << " id: " << v.info << std::endl;
		//std::cout << "Face(Elements) centers\n";
		//for (const auto& v: origins) std::cout << v.position << " id: " << v.info << std::endl;
		assert(size != 0);
		m_dispLists[0] = m_coreGL.createGLText(&nodes[0],sizeof(BaseVertex),size,
				ViewManagerInst().GetSettings()->NodeIDColor,NUM_VERTICES_LIST_COLORMAP,GL_FALSE);
		assert(m_dispLists[0]!=0);
		m_dispLists[1] = m_coreGL.createGLText(&origins[0],sizeof(BaseVertex),origins.size(),
				ViewManagerInst().GetSettings()->ElemIDColor,NUM_ELEMENTS_LIST_COLORMAP,GL_FALSE);
		assert(m_dispLists[1]!=0);

		//mfp_log_debug("disp0=%u disp1=%u",m_dispLists[0],m_dispLists[1]);
	}

	return m_dirty;
}


void RenderObject<TypeOfRenderer::ColorMap>::DoSetup(const GraphicsSettings* settings)
{
	/*
	std::vector<BaseVertex> origins, vertices;
	Mesh* pmesh = ModelCtrlInst().GetCurrentMesh();
	int num_vrts = pmesh->GetOriginsOfBoundaryFaces(origins,vertices);
    assert(num_vrts > 0);
	if (m_dispLists[0] == 0) {
		m_dispLists[0] = m_coreGL.createGLText((BaseVertex*)&vertices[0],sizeof(BaseVertex),vertices.size(),NUM_VERTICES_LIST_COLORMAP,GL_FALSE);
	}
	assert(m_dispLists[0]!=0);

	if (m_dispLists[1] == 0){
		m_dispLists[1] = m_coreGL.createGLText((BaseVertex*)&origins[0],sizeof(Origin),origins.size(),NUM_ELEMENTS_LIST_COLORMAP,GL_TRUE);
	}
	assert(m_dispLists[1]!=0);
	*/
}

void RenderObject<TypeOfRenderer::ColorMap>::DoRender(RenderParams* rparams)
{

	//mfp_debug("renderer colorMap %u",m_counts.size());
	//cout << "Light position: " << rparams.cLight.Position();
	m_program.Enable();
	FV_CHECK_ERROR_GL();
	m_program.SetUniform("posLight",rparams->cLight.Position().v);
	FV_CHECK_ERROR_GL();
	glBindVertexArray(m_vaoId);
	glMultiDrawElements(GL_TRIANGLE_STRIP, &m_counts[0], GL_UNSIGNED_INT,
			reinterpret_cast<const GLvoid**>(&m_addresses[0]), m_counts.size());

	glBindVertexArray(0);
	glUseProgram(0);

	if (!rparams->eMouseMode) {
		//mfp_log_debug("in static mouse mode %d",rparams->bShowNumOfVertices);
		if (rparams->bShowNumOfVertices) {
			//mfp_log_debug("Rendering vertices list\n");
			//m_coreGL.renderGLList(NUM_VERTICES_LIST_COLORMAP);
			glCallList(m_dispLists[0]);
			FV_CHECK_ERROR_GL();
		}
		if (rparams->bShowNumOfElements) {
			//mfp_log_debug("Rendering element list %d\n");
			//m_coreGL.renderGLList(NUM_ELEMENTS_LIST_COLORMAP);
			glCallList(m_dispLists[1]);
			FV_CHECK_ERROR_GL();
		}
	}

}

/********************** TypeOfRenderer::ColorMapSlice *********************/

RenderObject<TypeOfRenderer::ColorMapSlice>::RenderObject(const std::vector<CutPlane>& cuts)
: GLRenderBase(TypeOfRenderer::ColorMapSlice)
, m_program()
, m_cut_planes(cuts)
, m_dispLists({0,0})
{
	//mfp_log_debug("RenderObject ctr type: %d",this->m_type);
	bool res = m_program.AddShader(GLShader::GetPath(TRI_VERT),GL_VERTEX_SHADER);
	res = res && m_program.AddShader(GLShader::GetPath(TRI_FRAG), GL_FRAGMENT_SHADER);
	res = res && m_program.AddShader(GLShader::GetPath(TRISTRIP_GEOM), GL_GEOMETRY_SHADER);
	if (!res) throw fv_exception("Can't load shaders for color map cuts rendering!");
	else m_program.Complete();
}

RenderObject<TypeOfRenderer::ColorMapSlice>::~RenderObject()
{
	//mfp_log_debug("RenderObject dtr type: %d",this->m_type);
	//m_coreGL.deleteGLLists(m_dispLists[0],1);
	//m_coreGL.deleteGLLists(m_dispLists[1],1);
}

bool RenderObject<TypeOfRenderer::ColorMapSlice>::Create(void* arg_ptr)
{
	//mfp_debug("Init from field\n");
	if (!arg_ptr) return false;
	else {
		Field* pField = reinterpret_cast<Field*>(arg_ptr);

		// Temporary containers
		std::vector<Vertex> vertices;
		std::vector<unsigned int> indices;

		// Clear member containers
		m_count = pField->RenderDGCutted(m_cut_planes,vertices,indices/*,m_counts_fd,m_addresses_fd*/);
		assert(m_count!=0);
		//else pField->RenderSTD();
		//mfp_debug("Vertices: %u\tindices: %u\n",vertices.size(),indices.size());
		//int i=0;
		//for (const auto& vt : vertices)
		//	std::cout << "vertex " << i++ << " = " << vt.position << std::endl;
		//i = 0;
		//for (const auto& vt : indices)
		//		std::cout << "index " << i++ << " = " << vt << std::endl;
		if (vertices.empty()) return false;

		m_vboId = createBuffer(vertices.data(),sizeof(Vertex)*vertices.size(),GL_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_vboId == 0) return false;

		m_iboId = createBuffer(indices.data(),sizeof(unsigned int)*indices.size(),GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_iboId == 0) return false;

		glGenVertexArrays(1, &m_vaoId);
		glBindVertexArray(m_vaoId);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboId);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboId);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (const GLvoid*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)(4*sizeof(GLfloat)));
		glBindVertexArray(0);

		m_dirty = true;
	}

	if (m_dirty) {

		std::vector<BaseVertex> vertices;
		size_t num_origs = ModelCtrlInst().GetCurrentMesh()->GetVerticesOfCuts(m_cut_planes,vertices);

		if (vertices.size() - num_origs > 0) {
			m_dispLists[0] = m_coreGL.createGLText(&vertices[num_origs],sizeof(BaseVertex),vertices.size() - num_origs,
						ViewManagerInst().GetSettings()->NodeIDColor,NUM_VERTICES_LIST_CUTTED_COLORMAP,GL_FALSE);
			assert(m_dispLists[0]!=0);
		}

		if (num_origs) {
			m_dispLists[1] = m_coreGL.createGLText(&vertices[0],sizeof(BaseVertex),num_origs,
					ViewManagerInst().GetSettings()->ElemIDColor,NUM_ELEMENTS_LIST_CUTTED_COLORMAP,GL_FALSE);
			assert(m_dispLists[1]!=0);
		}
	}
	//mfp_log_debug("COLORMAP_SLICE: obtained: %u %u",m_dispLists[0],m_dispLists[1]);
	return m_dirty;
}

void RenderObject<TypeOfRenderer::ColorMapSlice>::DoSetup(const GraphicsSettings* settings)
{

}

void RenderObject<TypeOfRenderer::ColorMapSlice>::DoRender(RenderParams* rparams)
{
	m_program.Enable();
	m_program.SetUniform("posLight",rparams->cLight.Position().v);

	glDisable(GL_CULL_FACE);
	glBindVertexArray(m_vaoId);
	//mfp_log_debug("Number of lines to draw: %u\n",m_indices.size());
	glDrawElements(GL_TRIANGLES, m_count, GL_UNSIGNED_INT, (const GLvoid*)0);

	glBindVertexArray(0);
	glUseProgram(0);

	if (!rparams->eMouseMode) {
		if (rparams->bShowNumOfVertices && m_dispLists[0]!=0) { /*mfp_debug("Draw gllist=%d",m_dispLists[0]);*/glCallList(m_dispLists[0]); }
		if (rparams->bShowNumOfElements && m_dispLists[1]!=0) { /*mfp_debug("Draw gllist=%d",m_dispLists[1]);*/glCallList(m_dispLists[1]); }
	}

}


/********************** TypeOfRenderer::ColorMap STD*********************/

RenderObject<TypeOfRenderer::ColorMapStd>::RenderObject()
: GLRenderBase(TypeOfRenderer::ColorMapStd)
, m_program()
, m_num_tri_faces(0)
, m_dispLists({0,0})
{
	//mfp_log_debug("RenderObject ctr type: %d",this->m_type);
	const char* path = GLShader::GetPath(TRI_VERT);
	//mfp_log_debug("VERTEX PATH: %s",path);
	bool res = m_program.AddShader(GLShader::GetPath(TRI_VERT),GL_VERTEX_SHADER);
	path = GLShader::GetPath(TRI_FRAG);
	//mfp_log_debug("FRAGMENT PATH: %s",path);
	res = res && m_program.AddShader(GLShader::GetPath(TRI_FRAG), GL_FRAGMENT_SHADER);
	path = GLShader::GetPath(TRI_GEOM);
	//mfp_log_debug("GEOMETRY PATH: %s",path);
	res = res && m_program.AddShader(GLShader::GetPath(TRI_GEOM), GL_GEOMETRY_SHADER);
	if (!res) throw fv_exception("Can't load shaders for color map rendering!");
	else m_program.Complete();;

}

RenderObject<TypeOfRenderer::ColorMapStd>::~RenderObject()
{
	//mfp_log_debug("RenderObject dtr type: %d",this->m_type);
	//m_coreGL.deleteGLLists(m_dispLists[0],1);
	//m_coreGL.deleteGLLists(m_dispLists[1],1);

}


bool RenderObject<TypeOfRenderer::ColorMapStd>::Create(void* arg_ptr)
{
	//mfp_debug("Init from field TriangleShader\n");
	if (!arg_ptr) return false;
	else {
		Field* pField = reinterpret_cast<Field*>(arg_ptr);
		// Temporary containers
		std::vector<Vertex> vertices;
		std::vector<unsigned> indices;

		m_num_tri_faces = pField->RenderSTD(vertices,indices);

		m_count = indices.size();
		if (m_count == 0) return false;

		m_vboId = createBuffer(vertices.data(),sizeof(Vertex)*vertices.size(),GL_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_vboId == 0) return false;

		m_iboId = createBuffer(indices.data(),sizeof(unsigned)*indices.size(),GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_iboId == 0) return false;

		//std::cout <<"VERTICES:\n";
		//for (const auto& v : vertices) std::cout << v.position;
		//std::cout <<"INDICES:\n";
		//for (const auto& v : indices) std::cout << v << " ";

		glGenVertexArrays(1, &m_vaoId);
		glBindVertexArray(m_vaoId);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboId);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboId);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (const GLvoid*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)(4*sizeof(GLfloat)));
		glBindVertexArray(0);

		m_dirty = true;
	}

	if (m_dirty) {
		std::vector<BaseVertex> origins, nodes;
		unsigned size = ModelCtrlInst().GetCurrentMesh()->GetOriginsOfBoundaryFaces(origins,nodes);
		//std::cout << "Boundary vertices:\n";
		//for (const auto& v: nodes) std::cout << v.position << " id: " << v.info << std::endl;
		//std::cout << "Face(Elements) centers\n";
		//for (const auto& v: origins) std::cout << v.position << " id: " << v.info << std::endl;
		assert(size != 0);
		m_dispLists[0] = m_coreGL.createGLText(&nodes[0],sizeof(BaseVertex),size,
				ViewManagerInst().GetSettings()->NodeIDColor,NUM_VERTICES_LIST_COLORMAP,GL_FALSE);
		assert(m_dispLists[0]!=0);
		m_dispLists[1] = m_coreGL.createGLText(&origins[0],sizeof(BaseVertex),origins.size(),
				ViewManagerInst().GetSettings()->ElemIDColor,NUM_ELEMENTS_LIST_COLORMAP,GL_FALSE);
		assert(m_dispLists[1]!=0);

		//mfp_log_debug("disp0=%u disp1=%u",m_dispLists[0],m_dispLists[1]);
	}

	return m_dirty;
}

void RenderObject<TypeOfRenderer::ColorMapStd>::DoSetup(const GraphicsSettings* settings)
{
}

void RenderObject<TypeOfRenderer::ColorMapStd>::DoRender(RenderParams* rparams)
{
	m_program.Enable();
	m_program.SetUniform("posLight",rparams->cLight.Position().v);
	m_program.SetUniform("nTriangles",m_num_tri_faces);

	//mfp_log_debug("Just before render %d",m_count);
	glBindVertexArray(m_vaoId);
	//mfp_log_debug("Number of lines to draw: %u\n",m_indices.size());
	glDrawElements(GL_TRIANGLES, m_count, GL_UNSIGNED_INT, (const GLvoid*)0);

	glBindVertexArray(0);
	glUseProgram(0);

	if (!rparams->eMouseMode) {
		if (rparams->bShowNumOfVertices) glCallList(m_dispLists[0]);
		if (rparams->bShowNumOfElements) glCallList(m_dispLists[1]);
	}

}
/********************** TypeOfRenderer::ColorMapSlice *********************/

RenderObject<TypeOfRenderer::ColorMapStdSlice>::RenderObject(const std::vector<CutPlane>& cuts)
: GLRenderBase(TypeOfRenderer::ColorMapStdSlice)
, m_program()
, m_cut_planes(cuts)
, m_dispLists({0,0})
{
	//mfp_log_debug("RenderObject ctr type: %d",this->m_type);
	bool res = m_program.AddShader(GLShader::GetPath(TRI_VERT),GL_VERTEX_SHADER);
	res = res && m_program.AddShader(GLShader::GetPath(TRI_FRAG), GL_FRAGMENT_SHADER);
	res = res && m_program.AddShader(GLShader::GetPath(TRISTRIP_GEOM), GL_GEOMETRY_SHADER);
	if (!res) throw fv_exception("Can't load shaders for color map cuts rendering!");
	else m_program.Complete();
}

RenderObject<TypeOfRenderer::ColorMapStdSlice>::~RenderObject()
{
	//mfp_log_debug("RenderObject dtr type: %d",this->m_type);
	//m_coreGL.deleteGLLists(m_dispLists[0],1);
	//m_coreGL.deleteGLLists(m_dispLists[1],1);
}

bool RenderObject<TypeOfRenderer::ColorMapStdSlice>::Create(void* arg_ptr)
{
	//mfp_debug("Init from field\n");
	if (!arg_ptr) return false;
	else {
		Field* pField = reinterpret_cast<Field*>(arg_ptr);

		// Temporary containers
		std::vector<Vertex> vertices;
		std::vector<unsigned int> indices;

		// Clear member containers
		m_count = pField->RenderDGCutted(m_cut_planes,vertices,indices/*,m_counts_fd,m_addresses_fd*/);
		assert(m_count!=0);
		//else pField->RenderSTD();
		//mfp_debug("Vertices: %u\tindices: %u\n",vertices.size(),indices.size());
		//int i=0;
		//for (const auto& vt : vertices)
		//	std::cout << "vertex " << i++ << " = " << vt.position << std::endl;
		//i = 0;
		//for (const auto& vt : indices)
		//		std::cout << "index " << i++ << " = " << vt << std::endl;
		if (vertices.empty()) return false;

		m_vboId = createBuffer(vertices.data(),sizeof(Vertex)*vertices.size(),GL_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_vboId == 0) return false;

		m_iboId = createBuffer(indices.data(),sizeof(unsigned int)*indices.size(),GL_ELEMENT_ARRAY_BUFFER,GL_STATIC_DRAW);
		if (m_iboId == 0) return false;

		glGenVertexArrays(1, &m_vaoId);
		glBindVertexArray(m_vaoId);
		glBindBuffer(GL_ARRAY_BUFFER, m_vboId);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iboId);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(vertices[0]), (const GLvoid*)0);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)(4*sizeof(GLfloat)));
		glBindVertexArray(0);

		m_dirty = true;
	}

	if (m_dirty) {

		std::vector<BaseVertex> vertices;
		size_t num_origs = ModelCtrlInst().GetCurrentMesh()->GetVerticesOfCuts(m_cut_planes,vertices);

		if (vertices.size() - num_origs > 0) {
			m_dispLists[0] = m_coreGL.createGLText(&vertices[num_origs],sizeof(BaseVertex),vertices.size() - num_origs,
						ViewManagerInst().GetSettings()->NodeIDColor,NUM_VERTICES_LIST_CUTTED_COLORMAP,GL_FALSE);
			assert(m_dispLists[0]!=0);
		}

		if (num_origs) {
			m_dispLists[1] = m_coreGL.createGLText(&vertices[0],sizeof(BaseVertex),num_origs,
					ViewManagerInst().GetSettings()->ElemIDColor,NUM_ELEMENTS_LIST_CUTTED_COLORMAP,GL_FALSE);
			assert(m_dispLists[1]!=0);
		}
	}
	return m_dirty;
}

void RenderObject<TypeOfRenderer::ColorMapStdSlice>::DoSetup(const GraphicsSettings* settings)
{
}

void RenderObject<TypeOfRenderer::ColorMapStdSlice>::DoRender(RenderParams* rparams)
{
	m_program.Enable();
	// Set up a light position
	m_program.SetUniform("posLight",rparams->cLight.Position().v);
	//m_program.SetUniform("nTriangles",m_num_tri_faces);

	glDisable(GL_CULL_FACE);
	glBindVertexArray(m_vaoId);
	//mfp_log_debug("Number of lines to draw: %u\n",m_count);
	glDrawElements(GL_TRIANGLES, m_count, GL_UNSIGNED_INT, (const GLvoid*)0);

	glBindVertexArray(0);
	glUseProgram(0);

	if (!rparams->eMouseMode) {
		if (rparams->bShowNumOfVertices && m_dispLists[0]!=0) glCallList(m_dispLists[0]);
		if (rparams->bShowNumOfElements && m_dispLists[1]!=0) glCallList(m_dispLists[1]);
	}
}

}


