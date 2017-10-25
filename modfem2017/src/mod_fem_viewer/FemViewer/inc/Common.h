#ifndef __COMMON_H
#define __COMMON_H

#include<cstdlib>
#include"fv_compiler.h"

namespace FemViewer {

template< class T>
class mfvSingleton : private T {
  private:
	mfvSingleton(void) {/**/}
	~mfvSingleton(void) {/**/}

  public:
	static T& instance(void) {
		static mfvSingleton<T> obj;
		return(obj);
	}
};

template< class T>
class mfvSingleton<T*> : private T {
  private:
	static mfvSingleton<T*> * _objp;
	mfvSingleton(void) {
		atexit(destroy);
	}
	~mfvSingleton(void) {/**/}

  public:
	static void destroy(void) {
		if (_objp) { delete _objp; _objp= NULL; }
	}
	static T& instance(void) {
		if (!_objp) {
			_objp = new mfvSingleton<T>();
		}
		return(*_objp);
	}
};

template<class T>
mfvSingleton<T*>* mfvSingleton<T*>::_objp = NULL;
/*
// Base Object Handle
class mfvBaseObject {
  protected:
	// Ctr & Dtr
	mfvBaseObject(void){}
	virtual ~mfvBaseObject(void){}
};*/

template<typename T,uint32_t N>
class Vec {
  public:
	T data[N];
	Vec()    		{ memset(data,0x0, sizeof(T) * N); }
	Vec(T x) 		{ for (uint32_t i(0);i<N;++i) data[i] = x; }
	Vec(const T *x) { for (uint32_t i(0);i<N;++i) data[i] = x[i]; };
};

// dimension
template<uint32_t D, class C>
class Array {
public:
	Array(const Vec<size_t, D> &dim) : dimension(dim), arraySize(1)
	{
		for (uint32_t i(0); i < D; ++i) arraySize *= dim[i];
		data = new C[arraySize];
	}
	virtual ~Array() { delete [] data; }
	void clear() { memset(this->data, 0x0, sizeof(C) * this->arraySize); }
	void set(const C &c) { for (uint32_t i = 0; i < arraySize; ++i) data[i] = c; }
	Vec<size_t, D> dimension;
	size_t arraySize;
	C *data;
};

template<class C>
C bilinearInterpolate1(
	const double &s, const double &t,
	const C &a, const C &b, const C &c, const C &d)
{
	C e = a * (1 - s) + b * s;
	C f = c * (1 - s) + d * s;
	return e * (1 - t) + f * t;
}

// dimension
template<class C>
class Array2 : public Array<2, C>
{
public:
	Array2(const Vec<size_t, 2> &dim) : Array<2, C>(dim) {}
	template<typename T> C interpolate(const T &s, const T &t) const
	{
		assert(s <= 1 && t <= 1);
		double x = s * (this->dimension[0] - 1), y = t * (this->dimension[1] - 1);
		uint32_t xi = static_cast<uint32_t>(x), yi = static_cast<uint32_t>(y);
		double dx = x - xi, dy = y - yi;
		const C &a = (*this)[yi][xi];
		const C &b = (*this)[yi][std::min(xi + 1, this->dimension[0] - 1)];
		const C &c = (*this)[std::min(yi + 1, this->dimension[1] - 1)][xi];
		const C &d = (*this)[std::min(yi + 1, this->dimension[1] - 1)][std::min(xi + 1, this->dimension[0] - 1)];
		return bilinearInterpolate1<C>(dx, dy, a, b, c, d);
	}
	C* operator [] (size_t j) { return this->data + j * this->dimension[0]; }
	C* operator [] (size_t j) const { return this->data + j * this->dimension[0]; }
};

}// end namespace FemViewer



#endif
