#ifndef _GRAPHIC_DATA_H_
#define _GRAPHIC_DATA_H_

#include <string>
#include <iostream>
#include "Object.h"

namespace FemViewer {

class BBox3D;
class RenderParams;
class ModelCtrl;

class GraphicData
{
public:

	GraphicData();
	~GraphicData();

	void DeleteDisplayLists();

	void DumpCharacteristics(std::ostream& os,
                             const std::string& pIndentation);

	void EnableSmoothingMode();

	const BBox3D& GetGlobalBBox3D() const;


	void Render(RenderParams&  pParams);

	void Reset();

	void SetOptimizerValue(const int OptimalValue);


	Object& GetRootObject() { return _oRootObject; }
	const Object& GetRootObject() const { return _oRootObject; }

	public:

		bool    bFlagSmoothing;
		int     iOptimizerValue;
	private:
		Object _oRootObject;

	private:
		// Block the use of those
		GraphicData(const GraphicData&);
		GraphicData& operator=(const GraphicData&);
};

} // end namespace FemViewer


#endif /* _GRAPHIC_DATA_H_
*/ 
