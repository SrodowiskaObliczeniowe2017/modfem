#ifndef _FIELD_IMPORTER_H_
#define _FIELD_IMPORTER_H_



#include "BaseImporter.h"
#include "../../include/Enums.h"
#include "../../include/BaseField.h"
#include "../include/fv_compiler.h"


namespace IOmgr {

	//class Mesh;

	class FieldImporter : public BaseImporter {

	public:

		FieldImporter(FemViewer::BaseField* pfield_ ) : _pfield(pfield_) {}

		~FieldImporter() {}

		void* Reserve(FemViewer::ElemItem eltype, const unsigned int /*size = 0*/)
		{
			assert(eltype == FemViewer::FIELD || eltype == FemViewer::FIELDN);
			return _pfield->Reserve();
			/*else if(eltype == FIELDN) return _pfield*/
		}

		int nElems() const
		{
			return _pfield->nElems();
		}

		void ImportParams(const void* p_)
		{
			_pfield->ImportParams(p_);
		}

		void PostProcessing()
		{
			_pfield->Prepare();
		}

		//void AddField(const void* p_)
		//{
		//	FV_ASSERT(p_ != 0);
		//	_pfield->ImportField(p_);
		//}

		const void* GetItem () const
		{
			return _pfield->GetItem();
		}

		
	private:
		FemViewer::BaseField* _pfield;
	};
}

#endif /* _FIELD_IMPORTER_H_
*/
