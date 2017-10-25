/*
 * ArrayT.hpp
 *
 *  Created on: 2011-04-23
 *      Author: pawel
 */

#ifndef ARRAYT_HPP_
#define ARRAYT_HPP_
#include <cstddef>
#include <algorithm>
//#include <stack>
#include <vector>
#include <stdexcept>
#include"fv_compiler.h"
#include"fv_exception.h"
#include"ElemId.hpp"
#include "uth_log.h"
#ifdef FV_DEBUG
#include<iostream>
#endif

#define BOUND_ACTIVE 	0
#define ACTIVE_BOUND	1

namespace FemViewer {
template< typename IdT, typename T = ElemId<IdT> >
class ArrayT
{
public:
	typedef T				value_type;
	typedef T*				iterator;
	typedef const T*		const_iterator;
	typedef T&				reference;
	typedef const T&		const_reference;
	typedef std::size_t 	size_type;

private:
	static const size_type default_size = 8;
	T*			_data;
	size_type 	_max_size; 			// real size in indexes
	size_type 	_next_idx;      	// number of used indexes
	size_type	_count;				// number of specific elements
	uint8_t	_order_segregate; 	// order for segregating
	iterator& 	_sepIter1;			// internal iterator over boundary elements
	iterator&   _sepIter2;			// internal iterator over active   elements

	inline T*	_alloc(size_type n) { return reinterpret_cast<T*>( new char [ n * sizeof(T) ]); }
	inline void _check(size_type n) { if (n > _max_size) {
		char buf[24];
		sprintf(buf,"index: %u size %u\n",n,_max_size);
		throw std::out_of_range(buf); }
	}

public:

	explicit ArrayT(size_type size_ = default_size)
	: _data( _alloc(size_) )
	, _max_size(size_)
	, _next_idx(0)
	, _count(0)
	, _sepIter1(_data)
	, _sepIter2(_data)
	{
		//mfp_debug("ArrayT ctr\n");
		size_type i(0);
		try {
			for(;i<_max_size;++i) new (_data+i) T();
		} catch (...)
		{
			for(;i>0;--i) (_data+i)->~T();
			FV_FREE_ARR((char*)_data);
			throw;
		}
	}

	template<class U>
	ArrayT(const ArrayT<U>& ar)
	: _data(ar._data),
	  _max_size(ar._max_size),
	  _next_idx(ar._next_idx),
	  _count(ar._count),
	  _order_segregate(ar._order_segregate),
	  _sepIter1(ar._sepIter1),
	  _sepIter2(ar._sepIter1)
	{;}


	~ArrayT() {
		//mfp_debug("ArrayT dtr\n");
		clear();
	}

	// iterators
	iterator begin() 			  { return _data; }
	const_iterator begin()  const { return _data; }
	iterator end()				  { return _data + _next_idx; }
	const_iterator end()  	const { return _data + _next_idx; }
	/*iterator end()				  { return _data + _max_size; }
	const_iterator end()  	const { return _data + _max_size; }*/


	size_type capacity() const { return _max_size; }
	size_type size() 	 const { return _next_idx; }
	size_type nbound() 	const { assert(_next_idx >= _count); return (_next_idx - _count); }
	size_type nactive()	const {return _count; }

	uint_t order() const { return _order_segregate; }
	void set_order(uint_t order_);
	bool is_empty() const { return (_next_idx == 0); }


	void clear()
	{
		for( size_type i(0); i < _max_size; ++i)
			(_data + i)->~T();
		//std::cout<<"array clear before\n";
		FV_FREE_ARR((char*)_data);
		_max_size 	= 0;
		_next_idx 	= 0;
		_count		= 0;
		_order_segregate = 0x00;
		_sepIter1 = _sepIter2 	= 0;
	}

	void reserve(size_type capacity_);
	void resize(size_type  nsize_);
	void insert(const T& it_,size_type idx_);
	void push_back(const T& it) { insert(it,_next_idx); }
	void segregate(uint_t order_=BOUND_ACTIVE);

	reference fisrt() 				{ return _data[0]; }
	const_reference fisrt() const 	{ return _data[0]; }

	reference last() 				{ assert(_next_idx > 0); return _data[_next_idx-1]; }
	const_reference last() const 	{ assert(_next_idx > 0); return _data[_next_idx-1]; }

	reference at_pos(const size_type id_) {
		this->_check(id_);
		return *(this->begin() + id_);
	}

	const_reference at_pos(const size_type id_) const {
		assert((id_ >= 0) && (id_ < _max_size));
		return *(this->_begin() + id_);
	}

	reference operator[](const size_type id_) 	 	{ /*std::cout<<"setpr= " << _sepIter1 << std::endl*/;return *(_sepIter1 + id_); }
	const_reference operator[](const size_type id_) const  	{ /*std::cout<<"setpr= " << _sepIter1 << std::endl;*/return *(_sepIter1 + id_); }

	reference operator()(const size_type id_) { return *(_sepIter2 + id_); }
	const_reference operator()(const size_type id_) const { return *(_sepIter2 + id_); }

	reference next_boundary(const size_type idx_);
	const_reference next_boundary(const size_type idx_) const;

	reference next_active(const size_type idx_);
	const_reference next_active(const size_type idx_) const;
#ifdef FV_DEBUG
	void dumpArray() const
	{
		std::cout << "\tArray max size:\t\t" << _max_size
				  << "\n\tArray next index\t\t"  << _next_idx 
				  << "\n\tArray specific els:\t\t" << _count 
				  << std::endl;
	}
#endif
private:
	// anlarge capacity
	void _extend_array();
	// Block use of these
	//template<class U>
	//ArrayT(const ArrayT<U>&);

	template<class U>
	ArrayT& operator=(const ArrayT<U>&);

};

template<typename IdT,typename T>
void ArrayT<IdT,T>::set_order(uint_t order_)
{
	_order_segregate = order_;
	segregate();
}


template<typename IdT,typename T>
void ArrayT<IdT,T>::reserve(size_type size_)
{
	assert(size_ > 0);
	if (size_ <= _max_size) return;

	try {
		T* ptr(NULL);
		if ((ptr = this->_alloc(size_)) != NULL)
		{
			for(size_type i(0); i<_next_idx; ++i) {
				new(ptr + i) T(_data[i]);
				(_data+i)->~T();
			}
			FV_FREE_ARR((char*)_data);
			_data = ptr;
			_max_size  = size_;
		}

	}
	catch (std::bad_alloc&)
	{
		throw "Can't reserve memory for array!";
	}
}

template<typename IdT,typename T>
void ArrayT<IdT,T>::resize(size_type size_)
{
	assert(size_ >= 0);
	this->reserve(size_);
	for (size_type i = _next_idx; i< size_;++i)
		new (_data+i) T();
	for(size_type i = size_; i< _next_idx; ++i)
		(_data+i)->~T();
}

template<typename IdT,typename T>
void ArrayT<IdT,T>::insert(const T& it_,size_type idx_)
{
	//mfp_debug("insert index: %u\n",idx_);
	this->_check(idx_);

	if (_max_size == _next_idx) _extend_array();

	if (idx_ == _next_idx) {
		new(_data+idx_) T(it_);
	}
	else {
		new(_data + _next_idx) T(_data[_next_idx-1]);
		for(size_type i = _next_idx -1; i > idx_; --i)
			_data[i] = _data[i-1];
		_data[idx_] = it_;
	}
	++_next_idx;
	if (IS_ACTIVE(it_.id)) ++_count;
}

template<typename IdT,typename T>
void ArrayT<IdT,T>::segregate(uint_t order_)
{
	// Empty container
	if (_count == 0) return;

	// Invalid order of segregation
	if (order_!=BOUND_ACTIVE || order_!= ACTIVE_BOUND)
		throw fv_exception("Invalid order for segregation");

	_order_segregate 	  = order_;
	_sepIter1 = _sepIter2 = begin();

	iterator it_ (begin()), _it (&last());
	assert(it_ < _it);

	IdT first_mask(0), second_mask(0);

	if (_order_segregate == BOUND_ACTIVE)
	{
		SET_BOUND(first_mask);
		SET_ACTIVE(second_mask);
	}
	else {
		SET_ACTIVE(first_mask);
		SET_BOUND(second_mask);
	}
	std::cout<<"przed pierwszym: it_= " << it_ << " "<< _it << "=_it\n";
	// First segregate over first_mask
	while (it_ < _it) {

		while (it_->is_bit(first_mask)     && it_ < _it) ++it_;
		while (_it->is_bit_not(first_mask) && it_ < _it) --_it;
		if (it_ < _it) {
			std::swap(*it_,*_it);
			++it_;
			--_it;
		}
	}
	//std::cout<<"po pierwszym: it_= " << it_ << " "<< _it << "=_it\n";

	//std::cout<<"ArraySegregate: after boundary\n";

	// Separate boundary&active from boundary only
	_it = (it_->is_bit(first_mask)) ? it_ : --it_;

	if (_order_segregate == BOUND_ACTIVE) {
		_sepIter1 = begin();
		_sepIter2 = _it + 1;
	} else {
		_sepIter1 = _it + 1;
		_sepIter2 = begin();
	}

	iterator it(_it);
	it_ = begin();

	while (it_ < it) {

		while (it_->is_bit_not(second_mask) && it_ < it) ++it_;
		while (it->is_bit(second_mask)      && it_ < it) --it;
		if (it_ < it) {
			std::swap(*it_,*it);
			++it_;
			--it;
		}
	}

	//std::cout<<"ArraySegregate: after boundary2\n";

	it = it_->is_bit_not(second_mask) ? it_ : --it_;
	//std::cout<<"begin= " << begin() << std::endl;
	//std::cout<<"it_= " << it_ << std::endl;
	// Sort over boundary
	//std::sort(begin(),it_);

	//std::cout<<"it= " << it << std::endl;
	//std::cout<<"_it= " << _it << std::endl;
	// Sort over boundary&active
	std::sort(++it,_it++);

	// Separate active first
	it = it_ = _it;
	_it = &last();
	//std::cout<<"it= " << it << std::endl;
	//std::cout<<"_it= " << _it << std::endl;
	//std::cout<<"ArraySegregate: after boundary2\n";
	while (it_ < _it) {

		while (it_->is_bit(second_mask) && it_ < _it) ++it_;

		while (_it->is_bit_not(second_mask) && it_ < _it) --_it;
		if (it_ < _it) {
			std::swap(*it_,*_it);
			++it_;
			--_it;
		}
	}

	//std::cout<<"ArraySegregate: after boundary3\n";

	// Sort over active
	if (it_->is_bit_not(second_mask)) --it_;
	_it = it_ + 1;
	std::sort(it,it_); // why not qsort?
	std::sort(_it,&last());
}

template<typename IdT,typename T>
T& ArrayT<IdT,T>::next_boundary(const size_type idx_)
{
	assert(idx_ < _next_idx);
	assert(_next_idx > 0);
	return *(this->_sepIter1 + idx_ + 1);;
}

template<typename IdT,typename T>
const T& ArrayT<IdT,T>::next_boundary(const size_type idx_) const {
	assert(idx_ < _next_idx);
	assert(_next_idx > 0);
	return *(this->_sepIter1 + idx_ + 1);
}

template<typename IdT,typename T>
T& ArrayT<IdT,T>::next_active(const size_type idx_) {
	assert(idx_ < _next_idx);
	assert(_next_idx > 0);
	return *(this->_sepIter2 + idx_ + 1);
}

template<typename IdT,typename T>
const T& ArrayT<IdT,T>::next_active(const size_type idx_) const {
	assert(idx_ < _next_idx);
	assert(_next_idx > 0);
	return *(this->_sepIter2 + idx_ + 1);
}

template<typename IdT,typename T>
void ArrayT<IdT,T>::_extend_array()
{

	size_type n(1);
	if (_max_size > 0)
		n = _max_size << 1;
	this->reserve(n);
}

#undef BONUD_ACTIVE
#undef ACTIVE_BOUND

}
#endif /* ARRAYT_HPP_ */
