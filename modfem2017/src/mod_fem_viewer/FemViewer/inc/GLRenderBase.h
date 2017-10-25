/*
 * GLRenderBase.h
 *
 *  Created on: 14 paź 2015
 *      Author: pmaciol
 */

#ifndef MOD_FEM_VIEWER_FEMVIEWER_INC_GLRENDERBASE_H_
#define MOD_FEM_VIEWER_FEMVIEWER_INC_GLRENDERBASE_H_

#include"fv_txt_utls.h"


class GLcore;
namespace FemViewer {

// Forward decl.
class ViewManager;
class GraphicsSettings;

class GLRenderBase {

	template<typename> friend class RenderManager;

public:
	GLRenderBase(int type,bool visible=true);
	virtual ~GLRenderBase();

	virtual void Clear();
	virtual bool Create(void* argPtr) = 0;

	void Render(RenderParams* rparams);
    void Setup(const GraphicsSettings* gsettings);

	int  Type() { return m_type; }
	virtual bool IsEnabled() const = 0;
	//bool SetEnabled(bool state) { m_dirty = state; }
	bool IsVisible() { return m_visible; }
	void SetVisible(bool flag) { m_visible = flag; }

protected:

	virtual void DoSetup(const GraphicsSettings* gsettings) = 0;
	virtual void DoRender(RenderParams* rparams) = 0;

	GLCore& m_coreGL;
	//ViewManager& m_viewMgr;
	const int m_type;
	bool m_visible;
	GLuint m_vboId, m_iboId;
	GLuint m_vaoId;
	GLsizei m_count;

private:
	GLRenderBase(const GLRenderBase&);
	GLRenderBase& operator=(const GLRenderBase&);
};


template< int RenderType>
class RenderObject : public GLRenderBase {

public:
	RenderObject();
	~RenderObject();

	bool Create(void* arg_ptr) {};

	void Render(const RenderParams* rparams) {mfp_log_debug("Empty render");};


};


}// end namepsace FemViewer

#endif /* MOD_FEM_VIEWER_FEMVIEWER_INC_GLRENDERBASE_H_ */
