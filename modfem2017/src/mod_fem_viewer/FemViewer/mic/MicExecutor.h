/*
 * MicExecutor.h
 *
 *  Created on: 17 lis 2015
 *      Author: pmaciol
 */

#ifndef MOD_FEM_VIEWER_FEMVIEWER_MIC_MICEXECUTOR_H_
#define MOD_FEM_VIEWER_FEMVIEWER_MIC_MICEXECUTOR_H_

#include<memory>
#include "mic.h"
#include "fv_float.h"
#include "GraphicElem.hpp"

namespace FemViewer {
namespace MIC {

class MicRenderer;

class MicExecutor {
public:
	MicExecutor(int type = mic_CPU);
	virtual ~MicExecutor();

	int GetCore() const { return m_core; }

	bool Execute(void* (*pfn)(void*));
	void SetRenderingInegrationType(vol_rend_integr_t type);
protected:
	template<int Target>
	int SelectMeshElements(void* data);


private:
	MicRenderer* m_render;
	int m_core;
	vol_rend_integr_t m_integrtype;
	ScalarValueType m_epsilon;
	ScalarValueType m_hstep;

	std::shared_ptr<HostTransferFunction> m_tf;

	bool m_initialized;
	bool m_skipempty;

};

} /* namespace MIC */
} /* namespace FemViewer */

#endif /* MOD_FEM_VIEWER_FEMVIEWER_MIC_MICEXECUTOR_H_ */
