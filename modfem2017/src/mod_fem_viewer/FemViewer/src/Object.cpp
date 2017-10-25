#include "fv_inc.h"
#include "fv_txt_utls.h"
#include "defs.h"
#include "Log.h"
#include "Object.h"
#include "Matrix.h"
#include "PrimitiveAccumulator.h"
#include "../../utils/fv_str_utils.h"
#include "../../utils/fv_assert.h"
#include "Vec3D.h"
#include "Color.h"
#include "VtxAccumulator.h"
#include "VtxPrimitiveAccumulator.h"
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>

#include "../../include/ElemType.h"

namespace FemViewer {


// Should be put in a header file
static std::string extractCommandWord(const std::string&      pLine,
                               std::string::size_type& pEndWord)
{
	size_t start_word = pLine.find_first_not_of(" \t");
	
	if(start_word == std::string::npos){
		pEndWord = std::string::npos;
		return "";
	}

	pEndWord = pLine.find_first_of(" \n\t",start_word+1);
	return pLine.substr(start_word, pEndWord-start_word);

}


Object::Object(const char* name,bool usevbo)
  : aBoundingBox                          (),
    aCommands                             (),
    aFrozen                               (false),
    aGLDisplayListBoundingBox             (0),
    aGLDisplayListFast                    (0),
    aGLDisplayListFull                    (0),
    aName                                 (name),
    aNewPrimitiveAccumulatorNeeded        (true),
    aNewVtxAccumulatorNeeded           (true),
    aNewVtxPrimitiveAccumulatorNeeded(true),
    aPrimitiveAccumulators                (),
    //aRawMode                              (rawMode_not_in_raw_section),
    aRawModeArrowTipNbPolygons            (-1),
    aRawModeArrowTipProportion            (-1.0f),
    aSubObjects                           (),
    aVtxAccumulators                   (),
    aVtxPrimitiveAccumulators        ()
{
	//mfp_log_debug("ctr for %s\n",aName.c_str());
	//logd2(("Constructor Object: %s\n",aName.c_str()));
}

Object::~Object()
{
	//mfp_debug("Object dtr for %s\n",aName.c_str());

	if (aGLDisplayListFull != 0) {
		glDeleteLists(aGLDisplayListFull, 1);
	}
	if (aGLDisplayListBoundingBox != 0) {
		glDeleteLists(aGLDisplayListBoundingBox, 1);
	}
	if (aGLDisplayListFast != 0) {
		glDeleteLists(aGLDisplayListFast, 1);
	}

	  for (size_t i(0); i < aSubObjects.size(); ++i)
		  FV_FREE_PTR(aSubObjects[i].second);

	  for (size_t i(0); i < aPrimitiveAccumulators.size(); ++i)
		  FV_FREE_PTR(aPrimitiveAccumulators[i]);
	  
  
  

  
    for (size_t i(0); i < aVtxPrimitiveAccumulators.size(); ++i)
		  FV_FREE_PTR(aVtxPrimitiveAccumulators[i]);
	  
  

  
    for (size_t i(0); i < aVtxAccumulators.size(); ++i)
		  FV_FREE_PTR(aVtxAccumulators[i]);
	aSubObjects.clear();
    aPrimitiveAccumulators.clear();
	aVtxPrimitiveAccumulators.clear();
	aVtxAccumulators.clear();
    aCommands.clear();

}

// This function resets the object
// Writen becauseof feild destructor when exiting
void Object::Reset()
{
	//mfp_log_debug("Reset function: %s\n",aName.c_str());
	if (aGLDisplayListFull != 0) {
		glDeleteLists(aGLDisplayListFull, 1);
	}
	if (aGLDisplayListBoundingBox != 0) {
		glDeleteLists(aGLDisplayListBoundingBox, 1);
	}
	if (aGLDisplayListFast != 0) {
		glDeleteLists(aGLDisplayListFast, 1);
	}

	for (size_t i(0); i < aSubObjects.size(); ++i)
		  FV_FREE_PTR(aSubObjects[i].second);
	  
    for (size_t i(0); i < aPrimitiveAccumulators.size(); ++i)
		  FV_FREE_PTR(aPrimitiveAccumulators[i]);
	  
    for (size_t i(0); i < aVtxPrimitiveAccumulators.size(); ++i)
		  FV_FREE_PTR(aVtxPrimitiveAccumulators[i]);
	  
    for (size_t i(0); i < aVtxAccumulators.size(); ++i)
		  FV_FREE_PTR(aVtxAccumulators[i]);
	// Clear containers
	aSubObjects.clear();
    aPrimitiveAccumulators.clear();
	aVtxPrimitiveAccumulators.clear();
	aVtxAccumulators.clear();
    aCommands.clear();

	aNewPrimitiveAccumulatorNeeded = true;
    aNewVtxAccumulatorNeeded = true;
    aNewVtxPrimitiveAccumulatorNeeded = true;
}



// Delete the display lists of the current Object and all its SubObjects.
// Used primarily to force the reconstruction of the display lists when
// different OpenGL contexts can't share them
void Object::DeleteDisplayLists()
{

  if (aGLDisplayListFull != 0) {
    glDeleteLists(aGLDisplayListFull, 1);
    aGLDisplayListFull = 0;
  }
  if (aGLDisplayListBoundingBox != 0) {
    glDeleteLists(aGLDisplayListBoundingBox, 1);
    aGLDisplayListBoundingBox = 0;
  }
  if (aGLDisplayListFast != 0) {
    glDeleteLists(aGLDisplayListFast, 1);
    aGLDisplayListFast = 0;
  }

  SubObjects::iterator       lIter    = aSubObjects.begin();
  const SubObjects::iterator lIterEnd = aSubObjects.end  ();

  while (lIter != lIterEnd) {

    // The sub-Object might have been deleted
    Object* lSubObject = (*lIter).second;

    if (lSubObject != 0) {
      lSubObject->DeleteDisplayLists();
    }
    ++lIter;
  }
}

// Dump in ASCII the caracteristics of the Object and its parts
void Object::DumpCharacteristics(std::ostream&       pOstream,
                                 const std::string&  pIndentation,
                                 const Matrix<float>&    pTransformation)
{
  pOstream << pIndentation << "Object " << std::endl;
  std::string lIndentation = pIndentation + "  ";

  // Update the bounding box if necessary by calling getBoundingBox
  // instead of accessing directly aBoundingBox
  GetBBox3D();
  aBoundingBox.dumpCharacteristics(pOstream, lIndentation, pTransformation);

  pOstream << lIndentation << "Name               = " << aName << std::endl;

  if (aFrozen) {
    pOstream << lIndentation << "Frozen             = true" << std::endl;
  }
  else {
    pOstream << lIndentation << "Frozen             = false" << std::endl;
  }

#ifdef FV_DUMP_MEMORY_USAGE
  {
    typedef  std::vector<std::string>::size_type  SizeType;

    pOstream << lIndentation << "Memory used by the Object = " << sizeof(*this) << std::endl;
    pOstream << lIndentation << "Memory used by aName      = "
             << aName.size() << "/"
             << aName.capacity() << std::endl;

    Commands::const_iterator       lIterCommands     = aCommands.begin();
    const Commands::const_iterator lIterCommandsEnd  = aCommands.end  ();
    SizeType                       lCommandsSize     = sizeof(std::string)*aCommands.size    ();
    SizeType                       lCommandsCapacity = sizeof(std::string)*aCommands.capacity();

    while (lIterCommands != lIterCommandsEnd) {
      lCommandsSize     += lIterCommands->size();
      lCommandsCapacity += lIterCommands->capacity();
      ++lIterCommands;
    }

    pOstream << lIndentation << "Memory used by aCommands  = "
             << lCommandsSize << "/"
             << lCommandsCapacity << std::endl;
  }
#endif // #ifdef GLV_DUMP_MEMORY_USAGE

  pOstream << lIndentation << "Number of commands = " << aCommands.size() << std::endl;

  if (aCommands.size() > 0) {
    pOstream << lIndentation << "Commands (**): " << std::endl;
  }

  lIndentation += "  ";

  // Output every command and the dumpCharacteristics() on
  // commands that support it
  Commands::const_iterator       lIterCommands    = aCommands.begin();
  const Commands::const_iterator lIterCommandsEnd = aCommands.end  ();

  // Take a local copy of pTransformation in order to update
  // it when we encounter a gltranslate or glscale command
  Matrix<float> lTransformation = pTransformation;

  while (lIterCommands != lIterCommandsEnd) {

    std::string::size_type lEndWord = std::string::npos;
    std::string            lCommand;// = extractCommandWord(*lIterCommands, lEndWord);

    std::string lParameters;
    if(lEndWord != std::string::npos) {
      lParameters = trimString(lIterCommands->substr(lEndWord+1), " \t\n");
    }

    FV_ASSERT(!lCommand.empty());

    pOstream << lIndentation << "** " << *lIterCommands << std::endl;

    if (lCommand == "execute_primitive_accumulator_id") {
      FV_ASSERT(countWords(lParameters) == 1);
      int lPrimitiveAccumulatorId = atoi(lParameters.c_str());
      FV_ASSERT(lPrimitiveAccumulatorId >= 0);
      FV_ASSERT(lPrimitiveAccumulatorId <  static_cast<int>(aPrimitiveAccumulators.size()));
      FV_ASSERT(aPrimitiveAccumulators[lPrimitiveAccumulatorId] != 0);

      aPrimitiveAccumulators[lPrimitiveAccumulatorId]->DumpCharacteristics(pOstream,
                                                                           lIndentation + "  ",
                                                                           lTransformation);

    }
    else if (lCommand == "execute_vertex_primitive_accumulator_id") {
      FV_ASSERT(countWords(lParameters) == 1);
      int lVtxPrimitiveAccumulatorId = atoi(lParameters.c_str());
      FV_ASSERT(lVtxPrimitiveAccumulatorId >= 0);
      FV_ASSERT(lVtxPrimitiveAccumulatorId <  static_cast<int>(aVtxPrimitiveAccumulators.size()));
      FV_ASSERT(aVtxPrimitiveAccumulators[lVtxPrimitiveAccumulatorId] != 0);

      aVtxPrimitiveAccumulators[lVtxPrimitiveAccumulatorId]->dumpCharacteristics(pOstream,
                                                                                           lIndentation + "  ",
                                                                                           lTransformation);

    }
	else if (lCommand == "execute_vertex_accumulator_id") {
	}
    else if (lCommand == "execute_subobjects_id") {
      FV_ASSERT(countWords(lParameters) == 1);
      int lSubObjectId = atoi(lParameters.c_str());
      FV_ASSERT(lSubObjectId >= 0);
      FV_ASSERT(lSubObjectId <  static_cast<int>(aSubObjects.size()));

      // The sub-Object might have been deleted
      if (aSubObjects[lSubObjectId].second != 0) {
        aSubObjects[lSubObjectId].second->DumpCharacteristics(pOstream,
                                                       lIndentation + "  ",
                                                       lTransformation);
      }
      else {
        pOstream << lIndentation << "  Object deleted" << std::endl;
      }
    }
    // DIRECT OPENGL CALLS
    else if(lCommand == "gltranslate") {
      float tx, ty, tz;
      sscanf(lParameters.c_str(), "%f %f %f", &tx, &ty, &tz);

      // Add the translation to the transformation matrix lTM
      lTransformation.translate(Vec3D(tx, ty, tz));
    }
    else if(lCommand == "glscale") {
      float sx, sy, sz;
      sscanf(lParameters.c_str(), "%f %f %f", &sx, &sy, &sz);

      // Add the scaling to the transformation matrix lTM
      lTransformation.scale(Vec3D(sx, sy, sz));

    }

    ++lIterCommands;
  }
}

// Add mesh node
void Object::AddMeshNode(const fvmath::Vec3f& vertex)
{
	VtxAccumulator& vtxa = GetCurrentVtxAccumulator();
	vtxa.addNode(vertex);
}

void Object::AddIndexOfNode(const unsigned int index,
							const int dest)
{
	VtxAccumulator& vtxa = GetCurrentVtxAccumulator();
	vtxa.addIndexOfNode(index,dest);
}

void Object::AddColorOfNode(const fvmath::Vec3f& col)
{
	VtxAccumulator& vtxa = GetCurrentVtxAccumulator();
	vtxa.addNodeColor(col);
}

void Object::AddNormalToNode(const fvmath::Vec3f& norm)
{
	VtxAccumulator& vtxa = GetCurrentVtxAccumulator();
	vtxa.addNormal2Node(norm);
}

void Object::AddEdge(const unsigned int Edge[])
{
	VtxAccumulator& vtxa = GetCurrentVtxAccumulator();
	vtxa.addEdge(Edge[0],Edge[1]);
}

void Object::AddArrow(const Vec3D& p1,const Vec3D& p2,
                  const float tipprop, const int nrpolygons)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddArrow(p1,p2,tipprop,nrpolygons);
}

void Object::AddArrowColored(const Vec3D& p1, const ColorRGB& col1,
                             const Vec3D& p2, const ColorRGB& col2,
                             const float tipprop, const int nrpolygons)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddArrowColored(p1,col1,p2,col2,tipprop,nrpolygons);
}

void Object::AddPoint(const Vec3D& pt)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
    pa.AddPoint(pt);
}

void Object::AddPointColored(const Vec3D& pt, const ColorRGB& col)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
    pa.AddPointColored(pt,col);
}

void Object::AddLine(const Vec3D& p1, const Vec3D& p2)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddLine(p1,p2);
}

void Object::AddLineColored(const Vec3D& p1, const ColorRGB& col1,
                            const Vec3D& p2, const ColorRGB& col2)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddLineColored(p1,col1,p2,col2);
}

void Object::AddTriangle(const Vec3D& p1, const Vec3D& p2,
						 const Vec3D& p3)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddTriangle(p1,p2,p3);
}

void Object::AddTriangleColored(const Vec3D& p1, const ColorRGB& col1,
								const Vec3D& p2, const ColorRGB& col2,
								const Vec3D& p3, const ColorRGB& col3)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddTriangleColored(p1,col1,p2,col2,p3,col3);
}

void Object::AddTriangleNormalColored(const float nv[], const float p1[], const float p2[], const float p3[],
									  const ColorRGB& col1, const ColorRGB& col2, const ColorRGB& col3)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddTriangleNormalColored(nv,p1,p2,p3,col1,col2,col3);
}

void Object::AddQuad(const Vec3D& p1, const Vec3D& p2,
					 const Vec3D& p3, const Vec3D& p4)
{
	//FV_ASSERT( !Compare(p1,p2,p3,p4) );

	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddQuad(p1,p2,p3,p4);
}

void Object::AddQuadColored(const Vec3D& p1, const ColorRGB& col1,
							const Vec3D& p2, const ColorRGB& col2,
							const Vec3D& p3, const ColorRGB& col3,
							const Vec3D& p4, const ColorRGB& col4)
{
	PrimitiveAccumulator& pa = GetCurrentPrimitiveAccumulator();
	pa.AddQuadColored(p1,col1,p2,col2,p3,col3,p4,col4);
}
//
//void Object::AddText(const Vec3D& pt, const std::string& text, const char* font)
void Object::AddText(const float x, const float y, const float z, eFont font, const std::string& stext)
{
	FV_ASSERT(stext != "");
	char buff[256];
	sprintf(buff,"%f %f %f ",x,y,z);
	std::string label("text ");
	label += buff;

	switch(font)
	{
	case fixed13:
		label += "fixed13 ";
		break;
	case fixed15:
		label += "fixed15 ";
		break;
	case times10:
		label += "times10 ";
		break;
	case times24:
		label += "times24 ";
		break;
	case helvetica10:
		label += "helvetica10 ";
		break;
	case helvetica12:
		label += "helvetica12 ";
		break;
	case helvetica18:
		label += "helvetica24 ";
		break;
	}

	label += stext;
	aCommands.push_back(label);
 
}

Object*& Object::AddNewObject(const std::string& sName, int* pIdx)
{
	aSubObjects.push_back(std::make_pair(false,new Object(sName.c_str())));
    char cmd[64];
	size_t count = aSubObjects.size()-1;
    sprintf(cmd,"execute_subobjects_id %d", (int)count);
    aCommands.push_back(cmd);

	if (pIdx != NULL) *pIdx = static_cast<int>(count);

	return aSubObjects.back().second;
}

Object* Object::ResetSubObject(const unsigned int idx_)
{
	Object* pObj = NULL;
	try {
		std::string sName_ = aSubObjects.at(idx_).second->GetName();
		delete aSubObjects.at(idx_).second;
		pObj = new Object(sName_.c_str());
		aSubObjects[idx_].first = false;
		aSubObjects[idx_].second = pObj;
	} catch (std::out_of_range& /*ex*/)
	{
		std::cerr << "Error in reset subobject! index: " << idx_ << " is out of range\n";
	}
	

	return pObj;
}

void Object::ClearSubObject(const unsigned int idx_)
{
	try {
		Object*& pObj = aSubObjects.at(idx_).second;
		if( pObj != NULL) {
			delete pObj;
			pObj = NULL;
		}
	} catch (std::out_of_range& /*ex*/)
	{
		std::cerr << "Error in clear subobject! index: " << idx_ << " is out of range\n";
	}

}
// idx  < 0 - switch freez on/off for all subobjects
// idx >= 0 - switch freez in/off for object of idx and its all subobjects
// if idx is out of range - nothing is done
void Object::SwitcObject(const int idx, const bool hide)
{
	aSubObjects[idx].first = hide;
}


void Object::SetDrawColor(const ColorRGB& col)
{
	char cmd[64];
	sprintf(cmd,"color %f %f %f",col.R,col.G,col.B);
	aCommands.push_back(std::string(cmd));
	//aNewPrimitiveAccumulatorNeeded = true;
   // aNewVtxPrimitiveAccumulatorNeeded = true;

}

void Object::SetPointSize(const float size)
{
	char cmd[24];
	sprintf(cmd,"pointsize %f",size);
	aCommands.push_back(std::string(cmd));
	//aNewPrimitiveAccumulatorNeeded = true;
   // aNewVtxPrimitiveAccumulatorNeeded = true;

}

void Object::SetLineWidth(const float size)
{
	char cmd[24];
	sprintf(cmd,"linewidth %f",size);
	aCommands.push_back(std::string(cmd));
	//aNewPrimitiveAccumulatorNeeded = true;
    //aNewVtxPrimitiveAccumulatorNeeded = true;

}


const BBox3D& Object::GetBBox3D() const
{
	// If the object is frozen, then the BoundingBox is frozen too
	if (!aFrozen) {

    // We have to use the commands so that we can provide
    // an accurate BoundingBox by taking into account
    // the commands gltranslate and glscale. Glut
    // primitives are also taken into account as well as
    // the text command.
    aBoundingBox.Reset();

    // At first, we set the transformation to the identity
    Matrix<float> lTM;

    Commands::const_iterator       lIterCommands    = aCommands.begin();
    const Commands::const_iterator lIterCommandsEnd = aCommands.end  ();

    while (lIterCommands != lIterCommandsEnd) {

      std::string::size_type lEndWord = std::string::npos;
      std::string            lCommand = extractCommandWord(*lIterCommands, lEndWord);

      std::string lParameters;
      if(lEndWord != std::string::npos) {
        lParameters = trimString(lIterCommands->substr(lEndWord+1), " \t\n");
      }

      FV_ASSERT(!lCommand.empty());

      if (lCommand == "execute_primitive_accumulator_id") {
		  //std::cout<<"w exec primitw akumul\n";
        FV_ASSERT(countWords(lParameters) == 1);
        int lPrimitiveAccumulatorId = atoi(lParameters.c_str());
        FV_ASSERT(lPrimitiveAccumulatorId >= 0);
        FV_ASSERT(lPrimitiveAccumulatorId <  static_cast<int>(aPrimitiveAccumulators.size()));
        FV_ASSERT(aPrimitiveAccumulators[lPrimitiveAccumulatorId] != 0);
        const BBox3D lBoundingBox = aPrimitiveAccumulators[lPrimitiveAccumulatorId]->GetBBox3D();
		//std::cout<<"w exec primitw akumul po global\n";
        // Apply the transformation to the Bounding box
        // and add it to the object bounding box
        aBoundingBox += lTM * lBoundingBox;

      }
      else if (lCommand == "execute_vertex_primitive_accumulator_id") {
        FV_ASSERT(countWords(lParameters) == 1);
        int lVtxPrimitiveAccumulatorId = atoi(lParameters.c_str());
        FV_ASSERT(lVtxPrimitiveAccumulatorId >= 0);
        FV_ASSERT(lVtxPrimitiveAccumulatorId <  static_cast<int>(aVtxPrimitiveAccumulators.size()));
        FV_ASSERT(aVtxPrimitiveAccumulators[lVtxPrimitiveAccumulatorId] != 0);

        const BBox3D lBoundingBox = aVtxPrimitiveAccumulators[lVtxPrimitiveAccumulatorId]->getBBox3D();

        // Apply the transformation to the Bounding box
        // and add it to the object bounding box
        aBoundingBox += lTM * lBoundingBox;

      }
	  else if (lCommand == "execute_vertex_accumulator_id") {
        FV_ASSERT(countWords(lParameters) == 1);
        int lVtxAccumulatorId = atoi(lParameters.c_str());
        FV_ASSERT(lVtxAccumulatorId >= 0);
        FV_ASSERT(lVtxAccumulatorId <  static_cast<int>(aVtxAccumulators.size()));
        FV_ASSERT(aVtxAccumulators[lVtxAccumulatorId] != 0);

        const BBox3D lBoundingBox = aVtxAccumulators[lVtxAccumulatorId]->getBBox3D();

        // Apply the transformation to the Bounding box
        // and add it to the object bounding box
        aBoundingBox += lTM * lBoundingBox;

      }
      else if (lCommand == "execute_subobjects_id") {
        FV_ASSERT(countWords(lParameters) == 1);
        int lSubObjectId = atoi(lParameters.c_str());
        FV_ASSERT(lSubObjectId >= 0);
        FV_ASSERT(lSubObjectId <  static_cast<int>(aSubObjects.size()));

        // The sub-Object might have been deleted
        if (aSubObjects[lSubObjectId].second != 0) {
          BBox3D lBoundingBox = aSubObjects[lSubObjectId].second->GetBBox3D();

          // Apply the transformation to the Bounding box
          // and add it to the object bounding box
          aBoundingBox += lTM * lBoundingBox;
        }

      }
      // DIRECT OPENGL CALLS
      else if(lCommand == "gltranslate") {
        float tx, ty, tz;
        sscanf(lParameters.c_str(), "%f %f %f", &tx, &ty, &tz);

        // Add the translation to the transformation matrix lTM
        lTM.translate(Vec3D(tx, ty, tz));

      }
      else if(lCommand == "glscale") {
        float sx, sy, sz;
        sscanf(lParameters.c_str(), "%f %f %f", &sx, &sy, &sz);

        // Add the scaling to the transformation matrix lTM
        lTM.scale(Vec3D(sx, sy, sz));

      }
      else if(lCommand == "text") {

        // Add just the position of the text in the BoundingBox
        float x, y, z;

        sscanf(lParameters.c_str(), "%f %f %f", &x, &y, &z);

        aBoundingBox += lTM * Vec3D(x, y, z);

      }

      ++lIterCommands;
    }

  }

  return aBoundingBox;
}

void Object::Render(RenderParams& pParams,bool hide)
{
  if (hide) return;
  GLuint& lGLDisplayList = GetGLDisplayList(pParams);

  if (lGLDisplayList == 0) {
    if (aFrozen) {
      ConstructDisplayList(pParams);
    }
  }

  //if(aTextures.isDefined()) {
  //  aTextures.startTexturing();
  //}

  if (lGLDisplayList != 0) {
    glCallList(lGLDisplayList);
  }
  else {
    // The Object is not frozen *or* GL doesn't want
    // to return a valid aGLDisplayList*. So we just
    // execute the commands
    Commands::const_iterator       lIterCommands    = aCommands.begin();
    const Commands::const_iterator lIterCommandsEnd = aCommands.end  ();

	/*std::cout << (aCommands.size() > 0) ? "comm\n" : "nokom\n";*/

    while (lIterCommands != lIterCommandsEnd) {
	  //std::cout << *lIterCommands << std::endl; 
      executeCommand(*lIterCommands, pParams);
      ++lIterCommands;
    }
  }

  //if(aTextures.isDefined()) {
  //  aTextures.endTexturing();
  //}

}

void Object::ConstructDisplayList(RenderParams& pParams)
{
  // First, we make sure that all the display lists for that
  // RenderParameters::RenderMode for all the children are constructed
  // We have to loop on the aCommands and not directly on aSubObjects
  // because some commands change the attributes of pParams and
  // influence the rendering

  Commands::const_iterator       lIterCommands    = aCommands.begin();
  const Commands::const_iterator lIterCommandsEnd = aCommands.end  ();
  bool                           lAllChildrenOk   = true;
  RenderParams               lParams          = pParams;

  while (lIterCommands != lIterCommandsEnd) {

    std::string::size_type lEndWord = std::string::npos;
    std::string            lCommand = extractCommandWord(*lIterCommands, lEndWord);

    std::string lParameters;
    if(lEndWord != std::string::npos) {
      lParameters = trimString(lIterCommands->substr(lEndWord+1), " \t\n");
    }

    if (lCommand == "execute_subobjects_id") {

      FV_ASSERT(countWords(lParameters) == 1);
      int lSubObjectId = atoi(lParameters.c_str());
      FV_ASSERT(lSubObjectId >= 0);
      FV_ASSERT(lSubObjectId <  static_cast<int>(aSubObjects.size()));

      Object* lSubObject = aSubObjects[lSubObjectId].second;

      // The sub-Object might have been deleted
      if (lSubObject != 0) {

        // this->aFrozen is true, then all the SubObjects and their
        // children should have aFrozen set to true. In other words,
        // an Object can't be frozen if one of its parts are still
        // subject to change. This is automatic in the recursive
        // parsing order, but we make sure here that we don't have a bug.
        FV_ASSERT(lSubObject->aFrozen);

        GLuint& lGLDisplayListSubObject = lSubObject->GetGLDisplayList(lParams);

        if (lGLDisplayListSubObject == 0 || glIsList(lGLDisplayListSubObject) == GL_FALSE) {
          lSubObject->ConstructDisplayList(lParams);
        }

        // If the construction of one display list fails,
        // then we are not in a position to construct the
        // current display list
        lAllChildrenOk = (lAllChildrenOk                         &&
                          lGLDisplayListSubObject           != 0 &&
                          glIsList(lGLDisplayListSubObject) == GL_TRUE);

      }
    }
    else if(lCommand == "draw_facetboundary_enable") {
      lParams.bFacetFrame = true;
    }
    else if(lCommand == "draw_facetboundary_disable") {
      lParams.bFacetFrame = false;
    }

    ++lIterCommands;
  }


  if (!lAllChildrenOk) {
    // We can't create a display list  for the current Object
    // because one of the children wasn't able to create its own
    // display list (OpenGL ressources exhausted?)
    FV_ASSERT(GetGLDisplayList(pParams) == 0);
  }
  else {

    GLuint& lGLDisplayList = GetGLDisplayList(pParams);

    // Starts recording lGLDisplayList
    // Reuse the already allocated aGLDisplayList if possibility

    if (lGLDisplayList == 0) {
      lGLDisplayList = glGenLists(1);
    }

    // There is a possibility that GL refuses to
    // generate a new display list
    if (glIsList(lGLDisplayList) == GL_FALSE) {
      // For the rest of the code, if for some obscure reason
      // glIsList(aGLDisplayList) return GL_FALSE even if
      // glGenLists returned a non-zero value, then we
      // consider the display list id as not initialised
      lGLDisplayList = 0;
    }

    if (lGLDisplayList != 0) {

      glNewList(lGLDisplayList, GL_COMPILE);

      lIterCommands = aCommands.begin();

      while (lIterCommands != lIterCommandsEnd) {
        executeCommand(*lIterCommands, pParams);
        ++lIterCommands;
      }

      glEndList();
    }

    FV_ASSERT(lGLDisplayList == 0 || glIsList(lGLDisplayList) == GL_TRUE);
  }
}

void Object::executeCommand(const std::string& pCommand, RenderParams& pParams) const
{

  std::string::size_type lEndWord = std::string::npos;
  std::string            lCommand = extractCommandWord(pCommand, lEndWord);

  std::string lParameters;
  if(lEndWord != std::string::npos) {
    lParameters = trimString(pCommand.substr(lEndWord+1), " \t\n");
  }

  FV_ASSERT(!lCommand.empty());

  if (lCommand == "execute_primitive_accumulator_id") {
    FV_ASSERT(countWords(lParameters) == 1);
    int lPrimitiveAccumulatorId = atoi(lParameters.c_str());
    FV_ASSERT(lPrimitiveAccumulatorId >= 0);
    FV_ASSERT(lPrimitiveAccumulatorId <  static_cast<int>(aPrimitiveAccumulators.size()));
    FV_ASSERT(aPrimitiveAccumulators[lPrimitiveAccumulatorId] != 0);

    aPrimitiveAccumulators[lPrimitiveAccumulatorId]->Render(pParams);
	//FV_CHECK_ERROR_GL();
  }

  else if (lCommand == "execute_vertex_primitive_accumulator_id") {
    FV_ASSERT(countWords(lParameters) == 1);
    int lVtxPrimitiveAccumulatorId = atoi(lParameters.c_str());
    FV_ASSERT(lVtxPrimitiveAccumulatorId >= 0);
    FV_ASSERT(lVtxPrimitiveAccumulatorId <  static_cast<int>(aVtxPrimitiveAccumulators.size()));
    FV_ASSERT(aVtxPrimitiveAccumulators[lVtxPrimitiveAccumulatorId] != 0);

    aVtxPrimitiveAccumulators[lVtxPrimitiveAccumulatorId]->render(pParams);
	//FV_CHECK_ERROR_GL();
  }

  else if (lCommand == "execute_vertex_accumulator_id") {
	  FV_ASSERT(countWords(lParameters) == 1);
      int lVtxAccumulatorId = atoi(lParameters.c_str());
      FV_ASSERT(lVtxAccumulatorId >= 0);
      FV_ASSERT(lVtxAccumulatorId <  static_cast<int>(aVtxAccumulators.size()));
      FV_ASSERT(aVtxAccumulators[lVtxAccumulatorId] != 0);

      aVtxAccumulators[lVtxAccumulatorId]->render(pParams);
  }

  else if (lCommand == "execute_subobjects_id") {
    FV_ASSERT(countWords(lParameters) == 1);
    int lSubObjectId = atoi(lParameters.c_str());
    FV_ASSERT(lSubObjectId >= 0);
    FV_ASSERT(lSubObjectId <  static_cast<int>(aSubObjects.size()));

    // The sub-Object might have been deleted
    if (aSubObjects[lSubObjectId].second != 0) {

      // We don't want the transformations of that
      // sub-Object to influence the transformation
      // of other sub-Objects and the parents ;-)
      glPushMatrix();
      glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT | GL_POLYGON_BIT | GL_ENABLE_BIT);

      aSubObjects[lSubObjectId].second->Render(pParams,aSubObjects[lSubObjectId].first);
	  //FV_CHECK_ERROR_GL();
      glPopAttrib();
      glPopMatrix();
    }

  }
  // DIRECT OPENGL CALLS
  else if(lCommand == "color") {
    float r,g,b;
    sscanf(lParameters.c_str(), "%f %f %f", &r, &g, &b);
    glColorMaterial(GL_FRONT_AND_BACK,GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glColor3f(r, g, b);

  }
  //else if(lCommand == "glpushmatrix") {
  //  FV_ASSERT(lParameters == "");
  //  glPushMatrix();

  //}
  //else if(lCommand == "glpopmatrix") {
  //  FV_ASSERT(lParameters == "");
  //  glPopMatrix();

  //}
  //else if(lCommand == "glbegin_triangles") {
  //  FV_ASSERT(lParameters == "");
  //  glBegin(GL_TRIANGLES);

  //}
  //else if(lCommand == "glbegin_lines") {
  //  FV_ASSERT(lParameters == "");
  //  glBegin(GL_LINES);

  //}
  //else if(lCommand == "glbegin_points") {
  //  FV_ASSERT(lParameters == "");
  //  glBegin(GL_POINTS);

  //}
  //else if(lCommand == "glend") {
  //  FV_ASSERT(lParameters == "");
  //  glEnd();

  //}
  //else if(lCommand == "glvertex") {
  //  float x,y,z;
  //  sscanf(lParameters.c_str(), "%f %f %f", &x, &y, &z);
  //  glVertex3f(x, y, z);

  //}
  //else if(lCommand == "gltranslate") {
  //  float x,y,z;
  //  sscanf(lParameters.c_str(), "%f %f %f", &x, &y, &z);
  //  glTranslatef(x, y, z);

  //}
  //else if(lCommand == "glscale") {
  //  float x,y,z;
  //  sscanf(lParameters.c_str(), "%f %f %f", &x, &y, &z);
  //  glScalef(x, y, z);

  //}
  else if(lCommand == "pointsize") {
    float lSize;
    sscanf(lParameters.c_str(), "%f", &lSize);
    glPointSize(lSize);

  }
  else if(lCommand == "linewidth") {
    float lSize;
    sscanf(lParameters.c_str(), "%f", &lSize);
    glLineWidth(lSize);

  //}
  //else if(lCommand == "texturing_linear") {
  //  // Not implemented at this level
  //}
  //// DRAWING PARAMETERS
  //else if(lCommand == "draw_single_sided") {
  //  FV_ASSERT(lParameters == "");
  //  glEnable(GL_CULL_FACE);
  //  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);

  //}
  //else if(lCommand == "draw_double_sided") {
  //  FV_ASSERT(lParameters == "");
  //  glDisable(GL_CULL_FACE);
  //  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);

  //}
  //else if(lCommand == "glenable_polygonoffset_fill") {
  //  FV_ASSERT(lParameters == "");
  //  glEnable(GL_POLYGON_OFFSET_FILL);
  //  glPolygonOffset(1,1);

  //}
  //else if(lCommand == "gldisable_polygonoffset_fill") {
  //  FV_ASSERT(lParameters == "");
  //  glDisable(GL_POLYGON_OFFSET_FILL);

  //}
  //else if(lCommand == "draw_facetboundary_enable") {
  //  pParams.bFacetFrame = true;
  //  sscanf(lParameters.c_str(), "%f %f %f", &pParams.BkgColor.R, &pParams.BkgColor.G, &pParams.BkgColor.B);
  //  glEnable(GL_POLYGON_OFFSET_FILL);
  //  glPolygonOffset(1,1);

  //}
  //else if(lCommand == "draw_facetboundary_disable") {
  //  pParams.bFacetFrame = false;

  }
  else if(lCommand == "text") {
    char lParamFont[64];
    float x,y,z;

    sscanf(lParameters.c_str(), "%f %f %f %s", &x, &y, &z, lParamFont);

    FV_ASSERT(lParameters.find("\"") != std::string::npos);
    std::string::size_type lStart = lParameters.find("\"") + 1;

    FV_ASSERT(lParameters.find("\"",lStart) != std::string::npos);
    std::string::size_type lEnd   = lParameters.find("\"",lStart);

    std::string lText = lParameters.substr(lStart,lEnd-lStart);

    // Make sure we dont draw with lighting enabled
    // but.. we dont want to mess around with the lighting settings
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);

    void* lFont = GLUT_BITMAP_8_BY_13;
    if(!strcmp(lParamFont,"fixed13")) {
      lFont = GLUT_BITMAP_8_BY_13;
    }
    else if(!strcmp(lParamFont,"fixed15")) {
      lFont = GLUT_BITMAP_9_BY_15;
    }
    else if(!strcmp(lParamFont,"times10")) {
      lFont = GLUT_BITMAP_TIMES_ROMAN_10;
    }
    else if(!strcmp(lParamFont,"times24")) {
      lFont = GLUT_BITMAP_TIMES_ROMAN_24;
    }
    else if(!strcmp(lParamFont,"helvetica10")) {
      lFont = GLUT_BITMAP_HELVETICA_10;
    }
    else if(!strcmp(lParamFont,"helvetica12")) {
      lFont = GLUT_BITMAP_HELVETICA_12;
    }
    else if(!strcmp(lParamFont,"helvetica18")) {
      lFont = GLUT_BITMAP_HELVETICA_18;
    }

    drawText3D(lText, x, y, z, lFont);

    // Revert to previous lighting settings
    glPopAttrib();

  }
  //else if(lCommand == "texture") {
  //  std::string lFilename = lParameters;
  //}
  else {
    FV_ASSERT(false);
  }
}

PrimitiveAccumulator& Object::GetCurrentPrimitiveAccumulator()
{
  if (aNewPrimitiveAccumulatorNeeded) {

    aPrimitiveAccumulators.push_back(new PrimitiveAccumulator(true));

    char lCommand[64];
    sprintf(lCommand, "execute_primitive_accumulator_id %d", (int)aPrimitiveAccumulators.size()-1);
    aCommands.push_back(lCommand);

    aNewPrimitiveAccumulatorNeeded = false;
  }

  FV_ASSERT(!aPrimitiveAccumulators.empty());
  FV_ASSERT(aPrimitiveAccumulators.back() != 0);
  return *(aPrimitiveAccumulators.back());
}

VtxAccumulator& Object::GetCurrentVtxAccumulator()
{
  if (aNewVtxAccumulatorNeeded) {

    aVtxAccumulators.push_back(new VtxAccumulator);

	char lCommand[64];
    sprintf(lCommand, "execute_vertex_accumulator_id %d", aVtxAccumulators.size()-1);
	aCommands.push_back(lCommand);

    aNewVtxAccumulatorNeeded = false;
  }

  FV_ASSERT(!aVtxAccumulators.empty());
  FV_ASSERT(aVtxAccumulators.back() != 0);
  return *(aVtxAccumulators.back());
}

VtxPrimitiveAccumulator& Object::GetCurrentVtxPrimitiveAccumulator()
{
  if (aNewVtxPrimitiveAccumulatorNeeded) {

    aVtxPrimitiveAccumulators.push_back(new VtxPrimitiveAccumulator(GetCurrentVtxAccumulator()));

    char lCommand[64];
    sprintf(lCommand, "execute_vertex_primitive_accumulator_id %d", static_cast<int>(aVtxPrimitiveAccumulators.size()-1));
    aCommands.push_back(lCommand);

    aNewVtxPrimitiveAccumulatorNeeded = false;
  }

  FV_ASSERT(!aVtxPrimitiveAccumulators.empty());
  FV_ASSERT(aVtxPrimitiveAccumulators.back() != 0);
  return *(aVtxPrimitiveAccumulators.back());
}

GLuint& Object::GetGLDisplayList(const RenderParams& pParams)
{
	switch (pParams.eRMode)
	{
	case RenderParams::eFull:
		return aGLDisplayListFull;
		break;
	case RenderParams::eBoundingBox:
		return aGLDisplayListBoundingBox;
		break;
	case RenderParams::eFast:
		return aGLDisplayListFast;
		break;
	default:
		FV_ASSERT(false);
    break;
  }
  return aGLDisplayListFast;
}

} //end namespace FemViewer
