/*
 * ColorBar.h
 *
 *  Created on: 17-09-2013
 *      Author: Paweł Macioł
 */

#ifndef COLORBAR_H_
#define COLORBAR_H_

#include "fv_inc.h"
#include "fv_txt_utls.h"
#include "Enums.h"
#include "Log.h"
#include "LegendValue.h"
#include "GLRenderBase.h"
#include <vector>
#include <cassert>
#include <iostream>
#include <GL/gl.h>

namespace FemViewer {

class Legend;

enum { left, up, };
enum { flat, gradient };



class ColorBar : public GLRenderBase
{
public:
	// Ctr & Dsr
	explicit ColorBar (int type = HORIZONTAL);

	virtual ~ColorBar (void) {
		mfp_log_debug("ColorBar dtr");
		Destroy();
	}
	virtual void Draw (void) {
		//glCallList(m_dispList);
		m_coreGL.renderGLList(m_dispList);
	}
	virtual void Reset (void) { Destroy(); Create(m_colorMap.get());  }
	//virtual void Update (void) {}
	virtual bool Create (void* argPtr);

	bool IsEnabled() const { return m_dirty; }

protected:
	virtual void DoSetup(const GraphicsSettings* settings) {}
	virtual void DoRender(RenderParams* rparams) {};

	unsigned int m_dispList;
	colorbar_config_t m_cbCfg;

	std::shared_ptr<Legend> m_colorMap;

	void Destroy (void);


private:
	// Not implemented
	ColorBar(const ColorBar&);
	ColorBar& operator=(const ColorBar&);

	static void Defaults(ColorBar& cb);

	bool m_dirty;

};

}


#endif /* COLORBAR_H_ */
