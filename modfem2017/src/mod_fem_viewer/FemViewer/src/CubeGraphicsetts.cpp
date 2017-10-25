#include "CubeGraphicSettings.h"

namespace FemViewer
{
	const int CubeGprahicSetts::CubeDims::max_cells = 50;

	CubeGprahicSetts::CubeGprahicSetts()
	: dims(),
	  originPt(),
	  bkgColor(1.f, 1.f, 1.f, 1.f),
	  lineColor(.3f),
	  lineWidth(2.f),
	  pointColor(.6f),
	  pointSize(2.f),
	  cubeMask(CUBE_NONE)
	{ }

	void CubeGprahicSetts::defaults()
	{
		bkgColor.R = 1.f;
		bkgColor.G = 1.f;
		bkgColor.B = 1.f;
		bkgColor.A = 1.f;

		lineColor = .3f;
		lineWidth = 2.f;
		
		pointColor = .6f;
		pointSize  = 2.f;

		cubeMask = CUBE_NONE;
	}

}// end namespace FemViewer
