#ifndef _BASE_IMPORTER_H_
#define _BASE_IMPORTER_H_

#include "../../include/Enums.h"

namespace IOmgr {

	class BaseImporter{
	public:
		/* Destructor
		**/
		virtual ~BaseImporter() {}

		/* Reserves need memory
		**/
		virtual void* Reserve(FemViewer::ElemItem /*elem_type*/, const unsigned int = 0 /*size*/) = 0;
		//{ return NULL; }

		/* launch post-import processing
		**/
		virtual void PostProcessing () = 0;

		/* Get the handle to the item
		**/
		virtual const void* GetItem () const = 0;
		
		/* Returns number of elems
		**/
		virtual int nElems() const = 0;

		/* Import parameters
		**/
		virtual void ImportParams(const void*) = 0;

		/* Import vertex
		**/
		virtual void AddVertex(const void* /*v_*/) {} 

		/* Import edge
		**/
		virtual void AddEdge(const void* /*e_*/) {} 

		/* Import face
		**/
		virtual void AddFace(const void* /*f_*/) {}

		/* Import element
		**/
		virtual void AddElement(const void* /*el_*/) {} 

	};	
}
#endif /* _BASE_IMPORTER_H_
*/
