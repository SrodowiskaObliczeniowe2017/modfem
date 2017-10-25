/*
 * BaseHandle.h
 *
 *  Created on: 07-03-2012
 *      Author: dwg
 */

#ifndef BASEHANDLE_H_
#define BASEHANDLE_H_

// System
#include<sstream>
#include<iostream>
#include<stdexcept>
#include<cassert>
// Internals
#include "Enums.h"

namespace FemViewer {



class BaseHandle
{
public:

	explicit BaseHandle(HandleType htype_=Unknown,int idx_=-1,bool act_=false)
		: _htype(htype_), _idx(idx_), _act(act_) {;}

	virtual ~BaseHandle(){;}

	inline BaseHandle& operator=(const BaseHandle& rhs_) {
		_htype = rhs_.type();
		_idx   = rhs_.idx();
		_act   = rhs_.is_active();
		return *this;
	}

	// get type
	HandleType& type() { return _htype; }
	const HandleType& type() const { return _htype; }
	// get _idx
	int& idx()		{ return _idx; }
	const int& idx() const { return _idx; }

	// get active
	bool is_active() const { return _act; }

	// change activation
	void activate(const bool flg_) { _act = flg_; }

	// check if _idx != -1
	bool is_valid_idx() const { return (_idx != -1); }

	// reset
	virtual void Reset() { _htype = Unknown; _idx=-1; _act=false; }

	// is valid type
	bool is_valid_type(const BaseHandle& rhs_) const {
		return _htype == rhs_._htype;
	}
	// operators
	bool operator==(const BaseHandle& rhs_) const {
		return (_idx == rhs_.idx()) && is_valid_type(rhs_); }

	bool operator!=(const BaseHandle& rhs_) const
	{ return (_idx != rhs_.idx()) && is_valid_type(rhs_); }

	bool operator <(const BaseHandle& rhs_) const
	{ return (_idx < rhs_.idx()) && is_valid_type(rhs_); }


protected:
	friend std::ostream& operator<<(std::ostream& os_,const BaseHandle& rhs_);
	HandleType _htype;
	int  	_idx;
	bool 	_act;

};

inline std::ostream& operator<<(std::ostream& os_,const BaseHandle& rhs_)
{
	return (os_ << rhs_._idx);
}

template < typename T,unsigned N >
class TArrayPtrs
{
public: typedef std::size_t		size_type;


private: // members________________________

	T* _ptrs;
	unsigned int _size;
	size_type	 _count;
	static const unsigned int _dim = N;

private: // methods________________________
	inline void _check(unsigned n) const {
		if (n >= _size) {
			std::stringstream ss("Index ");
			ss << n << " is over range [" << _size << "]";
			throw std::out_of_range(ss.str());
		}
	}
public: // methods__________________________

	explicit TArrayPtrs(unsigned int size_= _dim)
	: _ptrs(new T[size_]), _size(size_), _count(0) {
		assert(_ptrs!=(T*)0);
		this->Clear();
	}

	~TArrayPtrs() { delete [] _ptrs; }

	typedef T*				iterator;
	typedef const T*		const_iterator;


	iterator 	   begin()		{ return _ptrs; }
	const_iterator begin()const { return _ptrs; }
	iterator	   end()		{ return _ptrs + _size; }
	const_iterator end()  const { return _ptrs + _size; }
	unsigned Size() const { return _size; }
	unsigned Count() const { return _count; }
	void 	 Resize(unsigned new_size_) {
		delete _ptrs;
		_ptrs = new T[new_size_];
		_size = new_size_;
	}

		  T& operator [](unsigned idx_) 	  { return _ptrs[idx_]; }
	const T& operator [](unsigned idx_) const { return _ptrs[idx_]; }


	T& at(unsigned idx_) {
		this->_check(idx_);
		return (*this)[idx_];
	}

	bool Register(T const ptr_) {
		if (Exists(ptr_)) return false;
		register unsigned i(0);
		do {
			if (!(*this)[i]) {
				(*this)[i] = ptr_;
				_count++;
				//std::cout<< "addning new object " << *ptr_ << std::endl;
			    return true;
			}
		} while(++i<_size);
		return false;
	}

	void Unregister(unsigned idx_) {
		this->_check(idx_);
		(*this)[idx_] = (T)0;
	}

	void Clear() { for(unsigned i(0);i<_size;++i) (*this)[i] = (T)0; }

private:

	inline bool Exists(T const ptr_)
	{
		if (_count>0) {
			register unsigned int i=0;
			while (i<_size) {
				if (!(*this)[i])   {i++; continue;}
				if (*(_ptrs[i++]) == *ptr_) return true;
			}
		}
		return false;
	}
	/// Block use of those
	template < typename U,unsigned M >
	TArrayPtrs(const TArrayPtrs<U,M>&);
	template < typename U,unsigned M >
	TArrayPtrs<T,N>& operator=(const TArrayPtrs<U,M>&);

};




}// end namespace FemViewer





#endif /* BASEHANDLE_H_ */

