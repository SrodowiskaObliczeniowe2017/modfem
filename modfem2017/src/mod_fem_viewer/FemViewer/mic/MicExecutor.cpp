/*
 * MicExecutor.cpp
 *
 *  Created on: 17 lis 2015
 *      Author: pmaciol
 */

#include <mic/MicExecutor.h>
#include <mic/MicRenderer.h>
#include <omp.h>
#include "mmh_intf.h"

namespace FemViewer {
namespace MIC {

MicExecutor::MicExecutor(int type)
: m_render(nullptr)
, m_core(type)
, m_integrtype(mic_Riemann)
, m_epsilon(0.001)
, m_hstep(0.01)
, m_tf(new HostTransferFunction())
, m_initialized(true)
, m_skipempty(true)
{
	// TODO Auto-generated constructor stub
	if (m_core == mic_XEONPHI)
		m_render = new MicRenderer();
}

MicExecutor::~MicExecutor() {
	// TODO Auto-generated destructor stub
	if (m_render) delete m_render;
}

bool MicExecutor::Execute(void* (*pfn)(void*))
{
	return false;
}

void MicExecutor::SetRenderingInegrationType(vol_rend_integr_t type)
{
	if (m_integrtype == type) return;
	m_integrtype = type;
	// Notify chnage of integration type
}

template<>
int MicExecutor::SelectMeshElements<mic_CPU>(void *data)
{

}


} /* namespace MIC */
} /* namespace FemViewer */
