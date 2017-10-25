#include "GraphicData.h"
#include "../utils/fv_assert.h"
#include "../utils/fv_str_utils.h"
#include "Matrix.h"
#include "ModelControler.h"
#include "Object.h"
#include "RenderParams.h"
//#include <fcntl.h>
//#include <cstdio>

namespace FemViewer {



// Constructor; initialize with default values
GraphicData::GraphicData()
  : bFlagSmoothing(false),
    iOptimizerValue(100),
	_oRootObject("RootObject")
{
	//mfp_log_debug("Graphicdata ctr\n");
}

// Destructor.  Delete allocated openGL display lists
GraphicData::~GraphicData()
{
	//mfp_debug("Graphicdata dtr\n");
}

// Delete the display lists of the root Model and all its SubModels.
// Used primarily to force the reconstruction of the display lists when
// different OpenGL contexts can't share them
void GraphicData::DeleteDisplayLists()
{
	_oRootObject.DeleteDisplayLists();
}

// Dump in ASCII the caracteristics of the GraphicData
void GraphicData::DumpCharacteristics(std::ostream& os,
                                      const std::string&  pIndentation)
{
	os << std::endl;
	os << pIndentation << "*>  Characteristics" << std::endl;
	os << std::endl;

	//FV_ASSERT(pRootModel != 0);
	//pRootModel->DumpCharacteristics(os, pIndentation, Matrix());
	_oRootObject.DumpCharacteristics(os, pIndentation, Matrix<float>());

	os << std::endl;
	os << pIndentation << "<*  Characteristics" << std::endl;
	os << std::endl;
}



// Enable smoothing (only applicable to some primitives)
void GraphicData::EnableSmoothingMode()
{
	bFlagSmoothing = true;
}



// Returns the global BoundingBox based on the BoundingBox of all the objects
const BBox3D& GraphicData::GetGlobalBBox3D() const
{
	return _oRootObject.GetBBox3D();
}



// Render the accumulated graphic data.
void GraphicData::Render(RenderParams& pParams)
{

	// Default values
	glColor3f(1.0f , 1.0f, 1.0f);

	pParams.bSmoothNormals       = bFlagSmoothing;
	pParams.iPrimitiveOptimizerValue = iOptimizerValue;

	//pRootModel->Render(pParams);
	_oRootObject.Render(pParams);
}

//  Clear all accumulated graphic data

void GraphicData::Reset()
{
	_oRootObject.Reset();
}

void GraphicData::SetOptimizerValue(const int OptimalValue)
{
	iOptimizerValue = OptimalValue;
}


} // end namespace FemViewer

