#ifndef _FV_MODEL_H_
#define _FV_MODEL_H_

#include <set>
#include <vector>
#include <string>
#include <cassert>

#include "Legend.h"
#include "CutPlane.h"
#include "Object.h"
#include "ViewManager.h"
#include "GrphElement.hpp"
#include "Mesh.h"


#include "../../utils/fv_dictstr.h"
#include "../../include/fv_compiler.h"

namespace FemViewer {

	class BBox3D;
	class BaseField; //?
	class Matrix;
	class ColorRGB;
	class Object;
	class GraphElement;
	class ViewManager;


	const int MeshKey  = 0;
	const int FieldKey = 1;
	
	class BaseInfo {
	
	public:
	
		BaseInfo(u_int status_=0) 
			: _status(status_)  {}

		bool operator==(const BaseInfo& rhs_) {
			return _status == rhs_.status();
		}
		
		// get and set flags of status
		u_int status()  const { return _status; }

		void set_flags(u_int flags_) { _status = flags_; }

		// is flag set
		bool is_flag_set(u_int flag_) const { return ( _status & flag_ ) > 0; }
		// set flag
		void set_flag(u_int flag_) { _status |= flag_; }
		// unset flag
		void unset_flag(u_int flag_) { _status &= ~flag_; }
		// change flag
		void change(u_int flag_, bool f_) 
		{ _status = (f_) ? _status | flag_ : _status & ~flag_; }  

	  protected:
		u_int _status;
	};

	class MeshInfo : public BaseInfo {
	public:
		enum { Wireframe = 1, };

	public:
		MeshInfo(u_int status_=0) : BaseInfo(status_) {}

		MeshInfo(const MeshInfo& rhs_) {
			_status = rhs_.status();
		}

		MeshInfo& operator=(const MeshInfo& rhs_) {
			this->_status = rhs_.status();
			return *this;
		}

		bool is_wireframe() const { return is_flag_set(Wireframe); }
		void set_wireframe(bool f_) { change(Wireframe,f_); }
	};

	class FieldInfo : public BaseInfo {
	public:
		enum {
			Flooded			 = 1,
			Contoured		 = 2,
			ColoredContoured = 4,
			IsoSurfaced		 = 8,
			Vectored		 = 16,
		};
	public:
		FieldInfo(u_int status_=0) : BaseInfo(status_) {}

		FieldInfo(const FieldInfo& rhs_) {
			_status = rhs_.status();
		}

		FieldInfo& operator=(const FieldInfo& rhs_) {
			this->_status = rhs_.status();
			return *this;
		}

		bool is_flooded() const { return is_flag_set(Flooded); }
		void set_flooded(bool f_) { change(Flooded,f_); }

		bool is_contoured() const { return is_flag_set(Contoured); }
		void set_contoured(bool f_) { change(Contoured,f_); }

		bool is_coloredcontoured() const { return is_flag_set(ColoredContoured); }
		void set_coloredcontoured(bool f_) { change(ColoredContoured,f_); }

		bool is_isosurfaced() const { return is_flag_set(IsoSurfaced); }
		void set_isosurfaced(bool f_) { change(IsoSurfaced,f_); }

		bool is_vectored() const { return is_flag_set(Vectored); }
		void set_vectored(bool f_) { change(Vectored,f_); }
	};

	struct VisualInfo : public BaseInfo {
	public: 
		enum {
			Hidden	= 1,
			Cutted	= 2,
			Changed	= 4,
			Deleted = 8,
		};
	public:
		VisualInfo(u_int status_=0) : BaseInfo(status_) {}
		VisualInfo(const VisualInfo& rhs_) {
			_status = rhs_.status();
		}

		VisualInfo& operator=(const VisualInfo& rhs_) {
			this->_status = rhs_.status();
			return *this;
		}

		bool is_hidden() const { return is_flag_set(Hidden); }
		void set_hidden(bool f_) { change(Hidden,f_); }

		bool is_cutted() const { return is_flag_set(Cutted); }
		void set_cutted(bool f_) { change(Cutted,f_); }

		bool is_changed() const { return is_flag_set(Changed); }
		void set_changed(bool f_) { change(Changed,f_); }

		bool is_deleted() const { return is_flag_set(Deleted); }
		void set_deleted(bool f_) { change(Deleted,f_); }
	};

	class BaseHandle {
	public:
		explicit BaseHandle(int idx_=-1, bool act_=false) 
			: _idx(idx_), _active(act_) {}

		inline BaseHandle& operator=(const BaseHandle& rhs_) {
			_idx = rhs_.idx(); return *this;
		}

		void activate(const bool flg_) { _active = flg_; }

		// get _idx
		int idx() const { return _idx; }
		// get active
		bool active() const { return _active; }

		// check if _idx != -1
		bool is_valid() const { return (_idx != -1); }

		// check if active
		//bool is_active() const { return _active; }

		// reset
		void reset() { _idx=-1; _active=false; }
		// invalidate

		// operators
		bool operator==(const BaseHandle& rhs_) const
		{ return _idx == rhs_.idx(); }

		bool operator!=(const BaseHandle& rhs_) const
		{ return _idx != rhs_._idx; }

		bool operator <(const BaseHandle& rhs_) const
		{ return _idx < rhs_._idx; }

		// friends
		friend std::ostream& operator<<(std::ostream& os_,const BaseHandle& rhs_);

	private:
		int  _idx;
		bool _active;
	};

	inline std::ostream& operator<<(std::ostream& os_,const BaseHandle& rhs_)
	{
		return os_ << rhs_.idx();
	}

	struct MeshHandle : public BaseHandle {
		explicit MeshHandle(int idx_=-1) : BaseHandle(idx_) {}
	};

	struct FieldHandle : public BaseHandle {
		explicit FieldHandle(int idx_=-1) : BaseHandle(idx_) {}
	};

	struct grh_data 
	{
		Object*    objroot;
		const char *name;
		std::pair<int,Object*> object;
		bool        change;

		Mesh::arBndElems vecdt;

		explicit grh_data(Object* obj_=NULL,const char* nm_="unknown",int id_=-1,Object *ob_=NULL,bool flg_=true)
			: objroot(obj_), name(nm_), object(id_,ob_), change(flg_), vecdt() 
		{}

		/*explicit*/ grh_data(const grh_data& rhs_) {
			objroot = rhs_.objroot;
			name = rhs_.name;
			object = rhs_.object;
			change = rhs_.change;
			vecdt  = rhs_.vecdt;
		}

		grh_data& operator=(const grh_data& rhs_) {
			objroot = rhs_.objroot;
			name = rhs_.name;
			object = rhs_.object;
			change = rhs_.change;
			vecdt  = rhs_.vecdt;
			return *this;
		}

		~grh_data() { clear(); }

		void clear()
		{
			//objroot->ClearSubObject(object.first);
			//object.second = NULL;
			//change = true;
			vecdt.clear();
		}

		void reset() 
		{
			object.second = objroot->ResetSubObject(object.first);
			change = true;
			vecdt.clear();
		}
	};

	template<class TStatus, class TData, size_t nTSize> 
	class ModelItem : public BaseHandle {
	private:
		Object*    root;
		TData*     elem;
		VisualInfo visInfo;
		TStatus	   attributes;
		std::vector<grh_data> vecData;
		std::vector<int> idxs;

	public:	

		explicit ModelItem( Object* rt_=NULL,
			        TData* elm_=NULL,
					size_t size=nTSize,
					u_int vis_=0,
					u_int attr_=0,
					int idx_=-1,
					bool act_=false) 
			: BaseHandle(idx_,act_),
			  root(rt_),
		      elem(elm_),
			  visInfo(vis_),
			  attributes(attr_),
			  vecData(nTSize),
			  idxs()			  
		{}

		ModelItem& operator=(const ModelItem& rhs_)
		{
			BaseHandle::operator=(rhs_);
			root = rhs_.root;
			elem = rhs_.elem;
			visInfo = rhs_.visInfo;
			attributes = rhs_.attributes;
			vecData = rhs_.vecData;
			idxs = rhs_.idxs;
			return *this;
		}

		~ModelItem()
		{
			vecData.clear();
			idxs.clear();
		}

//		void ResetItem(const int idx_)
//		{
//			assert(idx_ != -1);
//
//			if( idx_ >= vecData.size()) return;
//			vecData[idx_].reset(root);
//		}

		void RemoveItem(const int idx_) 
		{
			assert( idx_ != -1);

			if( idx_ >= vecData.size()) return;
			vecData.erase(vecData.begin()+idx_);
		}

		void Reset()
		{
			this->reset();

			std::vector<grh_data>::iterator             itrb = vecData.begin();
			const std::vector<grh_data>::const_iterator itre = vecData.end();

			for(;itrb != itre; ++itrb)
			{
				itrb->reset();
			}
		}

		inline const TData* GetElem() const { return elem; }
		inline       TData* GetElem() { return elem; }

		inline const VisualInfo& GetVisualInfo() const { return visInfo; }
		inline       VisualInfo& GetVisualInfo() { return visInfo; }

		inline const TStatus& GetAttributes() const { return attributes; }
		inline		 TStatus& GetAttributes() { return attributes; }

		inline const std::vector<grh_data>& GetGraphicData() const { return vecData; }
		inline		 std::vector<grh_data>& GetGraphicData() { return vecData; }

		inline const std::vector<int>& GetIndexes() const { return idxs; }
		inline		 std::vector<int>& GetIndexes() { return idxs; }

	};
 

		


	class Model {
	public:
		typedef struct _struct_sol {
			static const int max_sols = 100;
			bool is_init;
			int nr_sols;
			int nr_equs;
			int nr_curr_sol;
			std::string formula;
			_struct_sol() 
				: is_init(false), 
				  nr_sols(1), 
				  nr_equs(1), 
				  nr_curr_sol(0), 
				  formula("v0")
			{}

			_struct_sol(const _struct_sol& rhs_)
				: is_init(rhs_.is_init), 
				  nr_sols(rhs_.nr_sols), 
				  nr_equs(rhs_.nr_equs), 
				  nr_curr_sol(rhs_.nr_curr_sol), 
				  formula(rhs_.formula)
			{
				assert(nr_sols >= 1 && nr_sols <= max_sols);
				assert(nr_equs >= 1 && nr_equs <= 10);
				assert(nr_curr_sol >= 0 && nr_curr_sol <= max_sols-1);
			}

		} Solution, *SolutionPtr;

		typedef struct _contener {
			Legend*     legend_ptr;
			SolutionPtr sol_ptr;
			CutPlane*	plane_ptr;
			Mesh::arBndElems* data_ptr;
			Object*		object_ptr;
		} contener, *contener_ptr;

		enum GeometryType {
			NoneType   = 0,
			MeshTypes   = 1,
			FieldTypes = 2,
			AllTypes   = 3,
		};

	private:
		static Model*		_self;
		static ViewManager* _vmgr;  // pointer to ViewManager
		static Object*		_root;
		static Legend		_legend;
		static CutPlane		_cutplane;
		static int			_mshlast;
		static int			_fldlast;
		static contener		_package;
		static Solution		_solution;

		Model();
	public:
		static void Init(ViewManager& vmgr_);
		static void Init() { if( ! _self) _self = new Model(); }
		static void Destroy();
		static void Reset() { Destroy(); Init(); }
		static void Refresh(GeometryType Type);
		static Legend& GetLegend() { return _legend; }
		static CutPlane& GetPalne() { return _cutplane; }
		static Solution& GetSolution() { return _solution; }
		static Model* GetInstance() { if(!_self) Init(); return _self; }
		static void SetViewManager(ViewManager* pVMgr) { _vmgr = pVMgr; Init(*_vmgr); }
		enum RenderAttributes {
			wireframe = 0x01,
			flooded   = 0x02,
			contoured = 0x04,
			cutted	  = 0x08,
		};

		u_int RenderOptions;
		
		bool isMeshes;
		bool isFields;

		bool isWireframe;
		bool isFlooded;
		bool isContoured;
		bool isColorContoured;
		bool isIsoSurfaced;
		bool isCutted;


		//CutPlane plane;
		//int		 nrVar;

		bool isPlaneChange;
		bool isContourChange;
		bool isSolutionChange;



		~Model();

		int  AddMesh(Mesh* mesh_);
		int  AddField(BaseField* pFiled);

		void ChangeCuttingPlane();
		void ChangeSolution();
		void RemoveMesh(const u_int idx);
		void RemoveField(const u_int idx);

		bool ActivateMeshObj(const u_int idx_, const bool flg_);
		bool ActivateFieldObj(const u_int idx_, const bool flg_);

		//static bool IsModelReady() const;

		void DeleteDisplayLists();

		void DumpCharacteristics(std::ostream& os, 
			const std::string& pIndentation, Matrix pTransformation);
	
		/* Returns the model BBox */
		const BBox3D& GetModelBBox3D() const {
			FV_ASSERT(_vmgr != NULL);
			return _vmgr->GetGraphicData().GetGlobalBBox3D();
		}

		void Update();
		void Redraw(bool redrawAll = true);

		const BBox3D& GetGlobalBBox3D() { return _root->GetBBox3D(); }
		

		void reset();
		void clear();


		int GetMeshNumber() const { return static_cast<int>(vMeshes.size()); }
		int GetFieldNumber() const;

		void SetGoemetry(Object* pRootObject);
		//void SetViewManager(ViewManager* pvmgr);

		BaseField* GetCurrentField();
	private:
		
		//void RemoveItem(Items& vItem, const int idx);

		typedef ModelItem<MeshInfo,Mesh,0> MeshObj;
		typedef ModelItem<FieldInfo,BaseField,0> FieldObj;
		typedef std::vector<MeshObj>        vecMeshes;
		typedef std::vector<FieldObj>		vecFields;
		typedef vecMeshes::iterator	      vecMeshIter;
		typedef vecFields::iterator		  vecFieldIter;
		typedef vecMeshes::const_iterator vecMeshConstIter;
		typedef vecFields::const_iterator vecFieldConstIter;
		template< class TMObj > 
			static inline bool InitModelObj(TMObj& ref_);

		void RenderMeshes(const std::vector<MeshObj*>& vMeshes2render_);
		void RenderFields(const std::vector<FieldObj*>& vFields2render_);

		bool CheckForChanges() const;

		//template< typename TMObj = vecMeshes >
		//	inline TMobj& GetActives(const TMObj& vModelObjs)
		//{
		//	TMObj tmpMdls;
		//	for (size_t i(0);i<vModelObjs.size();++i)
		//	{
		//		if (vModelObjs[i].is_active()) tmpMdls.push_back(
		//	}
		//}

		//template< typename TMObj >
		//	bool CheckVisualState(std::vector<TMObj&>& vHdls_,const u_int state_);

		//template< typename TMObj >
		//	bool CheckAttribState(std::vector<TMObj&>& vHdls_,const u_int state_);
	
		vecMeshes  vMeshes;
		vecFields  vFields;
		
	};

	template<class TMObj>
	inline bool Model::InitModelObj(TMObj& ref_)
	{
		return InitModelObj< typename TMObj::value_type >(ref_);
	}

	template<>
	inline	bool Model::InitModelObj<Model::MeshObj>(MeshObj& ref_)
	{
		int id;
		Object* obj = _root->AddNewObject("Mesh_Wireframe",&id);
		std::cout<<"first id= " << id << "obj= " << (void*)obj <<"\n";
		ref_.GetGraphicData().push_back( grh_data(_root,"Mesh_Wireframe",id,obj,true));
		
		obj = _root->AddNewObject("Mesh_Cutted_Wireframe",&id);
		std::cout<<"second id= " << id << "obj= " << (void*)obj <<"\n";
		ref_.GetGraphicData().push_back( grh_data(_root,"Mesh_Cutted_Wireframe",id,obj,true));
		std::cout <<"After MeghObj Tempate\n";
		return true;
	}

	template<>
	inline	bool Model::InitModelObj<Model::FieldObj>(FieldObj& ref_)
	{
		int id;
		Object* obj = _root->AddNewObject("Filed_Flooded",&id);
		ref_.GetGraphicData().push_back( grh_data(_root,"Field_Flooded",id,obj,true));

		obj = _root->AddNewObject("Filed_Contoured",&id);
		ref_.GetGraphicData().push_back( grh_data(_root,"Field_Contoured",id,obj,true));

		obj = _root->AddNewObject("Filed_Cutted_Flooded",&id);
		ref_.GetGraphicData().push_back( grh_data(_root,"Field_Cutted_Flooded",id,obj,true));

		obj = _root->AddNewObject("Filed_Cutted_Contoured",&id);
		ref_.GetGraphicData().push_back( grh_data(_root,"Field_Cutted_Contoure",id,obj,true));

		return true;

	}



	//template<typename TMObj,typename TAssign>
	//	inline void Model::AssignIndex(const int idx_)
	//{
	//	assert(idx_ > -1);
	//	typedef typename T::value_type lObj;
	//	lObj.push_back(idx_);
	//}

	//template<typename TMOb, class j>
	//	inline bool Model::Activate(const int idx_, const bool flg_)
	//{
	//	assert(idx_ > -1);

	//	typedef TMObj::value_type lObj;
	//	try 
	//	{
	//		lObj.at(idx_).active(flg_);
	//		return true;
	//	} 
	//	catch ( std::out_of_range& ex)
	//	{
	//		std::cerr << "Out of ragne for "<< idx_ " in activation\n!";
	//	}
	//	
	//	return false;
	//}

	//template< typename TMObj, typename TVObj >
	//	bool Model::CheckActive(std::vector<TMObj&>& vHdls_) /*const*/
	//{
	//	for(int i(0);i<static_cast<int>(TVObj.size());++i)
	//	{
	//		TMObj& lHdl = TVObj[i];
	//		if (lHdl.active()) vHdls_->push_back(lHdl);
	//	}

	//	return (vHdls_->size() == 0);
	//}


	//template< typename TMObj >
	//	bool CheckVisualState(std::vector<TMObj&>& vHdls_, const u_int state_)
	//{
	//	std::vectro<TMObj&> vhdlTmp;
	//	for(u_int i(0);i<static_cast<u_int>(vHdls_.size());++i)
	//	{
	//		TMObj& obj = vHdls_[i];
	//		if (obj.visInfo.is_flag_set(state_)) vhdlTmp = vHdls_[i];
	//	}

	//	vHdls_ = vhdlTmp;

	//	return (vhdlTmp.szie() == 0);
	//}

	//template< typename TMObj >
	//	bool CheckAttribState(std::vector<TMObj&>& vHdls_, const u_int state_)
	//{
	//	std::vectro<TMObj&> vhdlTmp;
	//	for(u_int i(0);i<static_cast<u_int>(vHdls_.size());++i)
	//	{
	//		TMObj& obj = vHdls_[i];
	//		if (obj.attributes.is_flag_set(state_)) vhdlTmp = vHdls_[i];
	//	}

	//	*vHdls_ = vhdlTmp;

	//	return (vhdlTmp.szie() == 0);
	//}

	


} // end namespace FemViewer

#endif /* _FV_MODEL_H_ 
*/
