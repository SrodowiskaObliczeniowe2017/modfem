/*
 * mic_operator.h
 *
 *  Created on: 16 lis 2015
 *      Author: pmaciol
 */

#ifndef MOD_FEM_VIEWER_FEMVIEWER_MIC_MIC_OPERATOR_H_
#define MOD_FEM_VIEWER_FEMVIEWER_MIC_MIC_OPERATOR_H_

namespace FemViewer {
namespace MIC {

class MicOperator {
public:
	static int initMIC(bool quiet=false);
	static void shutdownMIC();

protected:
	int m_numCores;
	int m_targetIdCore;

private:
	MicOperator();
	MicOperator(const MicOperator&);
	MicOperator& operator=(const MicOperator&);
};

}
}



#endif /* MOD_FEM_VIEWER_FEMVIEWER_MIC_MIC_OPERATOR_H_ */
