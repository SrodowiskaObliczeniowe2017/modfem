/*
 * ColorBar.cpp
 *
 *  Created on: 17-09-2013
 *      Author: Paweł Macioł
 */

#include "Enums.h"
#include "ColorBar.h"
#include "Legend.h"
#include "fv_txt_utls.h"
#include "Log.h"
#include <GL/glut.h>

namespace FemViewer {

ColorBar::ColorBar(int type)
: GLRenderBase(type)
, m_dispList(0)
, m_dirty(false)
{
	Defaults(*this);
}


// TODO: How to inform about orient ??
bool ColorBar::Create (void* argPtr)
{
	//mfp_log_debug("In create colorbar\n");
	if (!argPtr) return false;
	if (m_dirty) Destroy();

	Legend& leg = *reinterpret_cast<Legend*>(argPtr);
	m_cbCfg.type = leg.GetLegendType() ? 'g': 'f';
    m_cbCfg.orient = 'h';
    m_cbCfg.x = 0.05f;
    m_cbCfg.y = 0.05f;
	m_cbCfg.w = 0.9f;
    m_cbCfg.h = 0.03f;
	m_cbCfg.min = leg.GetMin();
    m_cbCfg.max = leg.GetMax();
	m_cbCfg.start = leg.GetStart();
	m_cbCfg.end = leg.GetEnd();
	if (leg.GetOutsideBlack()) {
		int i=0;
		for(i=0;i<3;++i) m_cbCfg.out_col[i] = 0.0f;
		m_cbCfg.out_col[i] = 1.0f;
	}
	else m_cbCfg.out_col[0] = -1.0f;
	std::vector<float> colors(3*(leg.GetRangesSize()+1));
	colors.clear();
	for (std::vector<colormap_value>::iterator it = leg.GetColors().begin(); it != leg.GetColors().end(); ++it)
	{
		if (it->type & OUTSIDE) {
			//mfp_log_debug("Warning there is OUTSIDE color");
			continue;
		}
		if (it->type & LVT_COLOR) {
			colors.push_back(it->rgb.r);
			colors.push_back(it->rgb.g);
			colors.push_back(it->rgb.b);
			//mfp_log_debug("Add color: {%f %f %f}",it->rgb.r,it->rgb.g,it->rgb.b);
		}
	}
    m_cbCfg.colors = colors.data();
	m_cbCfg.size = static_cast<int>(colors.size()/3);
	//mfp_log_debug("Ther is %lu colors",colors.size()/3);
	m_cbCfg.stride = 0;
	{
		int i=0;
		for (;i<3;++i) m_cbCfg.text_col[i] = 0.8f;
		m_cbCfg.text_col[i] = 1.0f;
	}

	m_cbCfg.font = GLUT_BITMAP_HELVETICA_10;
	m_cbCfg.theight = 0.02f;
	{
		int i=0;
		for (;i<3;++i) m_cbCfg.edge_col[i] = 0.8f;
		m_cbCfg.edge_col[i] = 1.0f;
	}
	m_cbCfg.line_width = 1.2f;
	m_cbCfg.out_width_coef = 0.1f;
	m_cbCfg.num_ranges = 5;

	m_dispList = createColorBar(&m_cbCfg, &m_vboId);
	m_dirty = (m_dispList != 0);
	return m_dirty;
}

void ColorBar::Destroy (void)
{
	//mfp_debug("ColorBar: Destroy\n");
	m_coreGL.deleteBuffers(&m_vboId, 1);
	m_coreGL.deleteGLLists(m_dispList,1);

}

void ColorBar::Defaults(ColorBar& cb)
{
#define ARRAY_SIZE(arr) sizeof(arr)/sizeof(arr[0]);
  //printf("In horizonatal kolora bar");
  static float colors[5][3] = {
	{0.0f,0.0f,1.0f}, {0.0f,1.0f,1.0f}, {1.0f,1.0f,0.0f}, {1.0f,1.0f,0.0f}, {1.0f,0.0f,0.0f}
  };

  switch(cb.m_type)
  {

  case HORIZONTAL:
	  cb.m_cbCfg.type = 'g';
	  cb.m_cbCfg.orient = 'h';
	  cb.m_cbCfg.x = 0.05f;
	  cb.m_cbCfg.y = 0.05f;
	  cb.m_cbCfg.w = 0.9f;
	  cb.m_cbCfg.h = 0.03f;
	  cb.m_cbCfg.min = 0.0;
	  cb.m_cbCfg.max = 1.0;// maxValue;
	  cb.m_cbCfg.start = 0.0;
	  cb.m_cbCfg.end = 1.0;
	  cb.m_cbCfg.out_col[0]=0.0f; cb.m_cbCfg.out_col[1]=0.0f; cb.m_cbCfg.out_col[2]=0.0f; cb.m_cbCfg.out_col[3]=1.f;
	  cb.m_cbCfg.colors = &colors[0][0];
	  cb.m_cbCfg.size = ARRAY_SIZE(colors);
	  cb.m_cbCfg.stride = 0;
	  cb.m_cbCfg.text_col[0]=0.8f;cb.m_cbCfg.text_col[1]=0.8f;cb.m_cbCfg.text_col[2]=0.8f;cb.m_cbCfg.text_col[3]=1.f;
	  cb.m_cbCfg.font = GLUT_BITMAP_HELVETICA_10;
	  cb.m_cbCfg.theight = 0.02f;
	  cb.m_cbCfg.edge_col[0]=0.8f;cb.m_cbCfg.edge_col[1]=0.8f;cb.m_cbCfg.edge_col[2]=0.8f;cb.m_cbCfg.edge_col[3]=1.0f;
	  cb.m_cbCfg.line_width = 1.2f;
	  cb.m_cbCfg.out_width_coef = 0.1f;
	  cb.m_cbCfg.num_ranges = 5;
	  break;
  case VERTICAL:
	  cb.m_cbCfg.type = 'f';
	  cb.m_cbCfg.orient = 'v';
	  cb.m_cbCfg.x = 0.9f;
	  cb.m_cbCfg.y = 0.05f;
	  cb.m_cbCfg.w = 0.03f;
	  cb.m_cbCfg.h = 0.9f;
	  cb.m_cbCfg.min = 0.0;
	  cb.m_cbCfg.max = 1.0;// maxValue;
	  cb.m_cbCfg.start = 0.0;
	  cb.m_cbCfg.end = 1.0;
	  cb.m_cbCfg.out_col[0]=0.0f; cb.m_cbCfg.out_col[1]=0.0f; cb.m_cbCfg.out_col[2]=0.0f; cb.m_cbCfg.out_col[3]=1.f;
	  cb.m_cbCfg.colors = &colors[0][0];
	  cb.m_cbCfg.size = ARRAY_SIZE(colors);
	  cb.m_cbCfg.stride = 0;
	  cb.m_cbCfg.text_col[0]=0.8f;cb.m_cbCfg.text_col[1]=0.8f;cb.m_cbCfg.text_col[2]=0.8f;cb.m_cbCfg.text_col[3]=1.f;
	  cb.m_cbCfg.font = GLUT_BITMAP_HELVETICA_10;
	  cb.m_cbCfg.theight = 0.02f;
	  cb.m_cbCfg.edge_col[0]=0.8f;cb.m_cbCfg.edge_col[1]=0.8f;cb.m_cbCfg.edge_col[2]=0.8f;cb.m_cbCfg.edge_col[3]=1.0f;
	  cb.m_cbCfg.line_width = 1.2f;
	  cb.m_cbCfg.out_width_coef = 0.1f;
	  cb.m_cbCfg.num_ranges = 5;
	  break;
  }
#undef ARRAY_SIZE
}





}// end namespace FemViewer


