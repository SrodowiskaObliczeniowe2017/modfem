/*
 * GLRenderBase.cpp
 *
 *  Created on: 14 paź 2015
 *      Author: pmaciol
 */
#include "GraphicsSettings.h"
#include "ViewManager.h"
#include "GLRenderBase.h"

namespace FemViewer {



GLRenderBase::GLRenderBase(int type, bool visible)
: m_coreGL(GLCore::instance())
//, m_viewMgr(ViewManagerInst())
, m_type(type)
//, m_dirty(false)
, m_visible(visible)
, m_vboId(0)
, m_iboId(0)
, m_vaoId(0)
, m_count(0)
{
	//mfp_log_debug("Renderable dtr type: %d",m_type);
}

GLRenderBase::~GLRenderBase()
{
	//mfp_log_debug("Renderable dtr type: %d",m_type);
	Clear();
}

void GLRenderBase::Clear()
{
	//m_dirty = false;
	m_visible = false;

	if (m_vboId != 0) {
		//mfp_debug("Deleting VBO : %d\n",m_vboId);
		glDeleteBuffers(1, &m_vboId);
		m_vboId = 0;
	}
	if (m_iboId != 0) {
		//mfp_debug("Deleting IBO : %d\n",m_iboId);
		glDeleteBuffers(1, &m_iboId);
		m_iboId = 0;
	}
	if (m_vaoId != 0) {
		//mfp_debug("Deleting VAO : %d\n",m_vaoId);
		glDeleteVertexArrays(1, &m_vaoId);
		m_vaoId = 0;
	}

	m_count = 0;

}

void GLRenderBase::Render(RenderParams* rparams)
{
	/*if (m_visible)*/ DoRender(rparams);
}

void GLRenderBase::Setup(const GraphicsSettings* gsettings) {
	DoSetup(gsettings);
}

}// end namespace FemViewer



