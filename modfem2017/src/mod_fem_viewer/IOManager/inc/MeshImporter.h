#ifndef _MESH_IMPORTER_H_
#define _MESH_IMPORTER_H_

#include "fv_compiler.h"

#include "BaseImporter.h"
#include "BaseMesh.h"

namespace IOmgr {

	/*template <class TMesh>*/
	class MeshImporter : public BaseImporter {

	public:

		MeshImporter(FemViewer::BaseMesh* pmesh_ ) : _pmesh(pmesh_) {}

		~MeshImporter() {}

		void* Reserve(FemViewer::ElemItem elem_type, const unsigned int size = 0)
		{
			assert(elem_type >= 0);
			assert(size    > 0);
			return _pmesh->Reserve(elem_type,size);
		}

		void PostProcessing()
		{
			_pmesh->Prepare();
		}

		const void* GetItem () const
		{
			return( NULL);// _pmesh->GetItem();
		}

		void ImportParams(const void *p_)
		{
			_pmesh->ImportParams(p_);
		}

		void AddNode(const void* p_)
		{
			_pmesh->AddNode(p_);
		}

		void AddEdge(const void* p_)
		{
			_pmesh->AddEdge(p_);
		}

		void AddFace(const void* p_)
		{
			_pmesh->AddFace(p_);
		}

		void AddElem(const void* p_)
		{
			_pmesh->AddElem(p_);
		}

		int nElems () const
		{ return _pmesh->nElems(); }

		
	private:
		FemViewer::BaseMesh* _pmesh;
	};
}

#endif /* _MESH_IMPORTER_H_
*/
