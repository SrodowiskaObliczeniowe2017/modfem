/*
 * ElemId.hpp
 *
 *  Created on: 24-09-2011
 *      Author: Paweł Macioł
 */

#ifndef ELEMID_HPP_
#define ELEMID_HPP_

#include"defs.h"
#include"fv_compiler.h"
#include<stdexcept>
#include<cstring>//memset
#include<type_traits>


namespace FemViewer {

#define EL_ID_ATTRIB(Tid) struct{ Tid eid : 24, etype : 1, faces : 5, active : 1, bound : 1; }
#define FA_ID_ATTRIB(Tid) struct{ Tid fid : 31, ftype : 1; }

template<typename T>
struct BaseId {
	union {
		T id;
		union {
			EL_ID_ATTRIB(T);
			FA_ID_ATTRIB(T);
		};
	};
};


template<typename T>
class Id : public BaseId<T>
{
public:
	Id(const T id_ = T(0)) { this->id = id_; }
	Id(const BaseId<T>& ref) { this->id = ref.id; }
	Id& operator=(const BaseId<T>& rhs) { this->id = rhs.id; return *this; }
	//virtual ~Id(){} // +8 bytes for vtable
	inline bool is_bit(T bit) const { return IS_BIT_SET(this->id,bit); }
	inline bool is_bit_not(T bit) const {return !is_bit(bit); }
	//inline virtual int type() const = 0;
};

template<typename T>
class ElemId : public Id<T>
{
public:
	typedef T value_type;

	ElemId(const T id_ = T(0)) : Id<T>(id_) {}
	ElemId(const BaseId<T>& ref) : Id<T>(ref) {}
	ElemId(const T eid_,
		   const int type_,
		   const bool f0 = false,
		   const bool f1 = false,
		   const bool f2 = false,
		   const bool f3 = false,
		   const bool f4 = false,
		   const bool act = false,
		   const bool bnd = false);
	ElemId& operator=(const BaseId<T>& rhs) { this->id = rhs.id; return *this; }

	// 0 - tetrahedron
    // 1 - prism
	bool is_prism() const { return this->is_bit(ELTYPE_BIT_POS(this->id)); }
	bool is_tetra() const { return !this->is_prism(); }

	bool operator < (const Id<T>& rhs) const {
		return (this->eid < rhs.eid);
	}

	//auto Index() const { return this->eid; }
};

template<typename T>
ElemId<T>::ElemId(const T eid_, const int type_,
		   const bool f0, const bool f1, const bool f2,
		   const bool f3, const bool f4, const bool act, const bool bnd)
: Id<T>()
{
	this->id = eid_;
	if (type_) SET_PRISM(this->id);
	if (f0) SET_FACE0(this->id);
	if (f1) SET_FACE1(this->id);
	if (f2) SET_FACE2(this->id);
	if (f3) SET_FACE3(this->id);
	if (f4) SET_FACE4(this->id);
	if (act) SET_ACTIVE(this->id);
	if (bnd) SET_BOUND(this->id);
}


template<typename T>
struct FaceId : public Id<T> {
public:
	typedef T value_type;
	FaceId(const T id_ = T(0)) : Id<T>(id_) {}
	FaceId(const BaseId<T>& ref) : Id<T>(ref) {}
	FaceId& operator=(const BaseId<T>& rhs) { this->id = rhs.id; return *this; }

	// 0 - traingle face
	// 1 - quad face
	bool is_triangle() const { return this->is_bit(FACE_TYPE_POS(this->id)); }
	bool is_quad() 		 const { return !this->is_triangle(); }

	bool operator < (const Id<T>& rhs) const {
		return (this->fid < rhs.fid);
	}
};

#undef EL_ID_ATTRIB
#undef FA_ID_ATTRIB

template <typename T>
struct CompareBndAct {
	bool operator() (T* el1, T* el2)
	{
		return (ELEM_ID(el1->id) < ELEM_ID(el2->id));
	}
};

template<typename T>
bool compare_func(T* it1, T* it2) { return (ELEM_ID(it1->el_id) < ELEM_ID(it2->el_id)); }

struct ElemInfo : public ElemId<id_t>
{
	int nodes[7];

	ElemInfo() : ElemId<id_t>() {
		memset(nodes,0x0,sizeof(nodes));
	}

	bool operator==(const ElemInfo& rh) {
		bool result = (nodes[0] == rh.nodes[0]) ? this->eid == rh.eid : false;
		return result;
	}

};

FV_STATIC_ASSERT(ElemInfo,8);


}// end namespace

#endif /* ELEMID_HPP_ */
