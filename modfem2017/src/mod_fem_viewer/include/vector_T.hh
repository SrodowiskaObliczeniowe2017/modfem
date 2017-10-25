/*
 * vector_T.hh
 *
 *  Created on: 11-01-2012
 *      Author: Paweł Macioł
 */

#ifndef VECTOR_T_hh__
#define VECTOR_T_hh__

#define AVEC2(T)    union { struct { T x,y      }; T data_[2]; }
#define AVEC3(T)    union { struct { T x,y,z;   }; T data_[3]; }
#define AVEC4(T)	union { struct { T x,y,z,w; }; T data_[4]; }

template,typename T,int N> struct VecT { T data_[N] };
template<> struct VecT<T,2> { AVEC2(T); };
template<> struct VecT<T,3> { AVEC3(T); };
template<> struct VecT<T,4> { AVEC4(T); };


template<typename T,int N> class CVecT: protected VecT<T,N> {

public:
	typedef typename VecT<T,N>	base;
	typedef typename CVecT<T,N>	vec_type;

	explicit inline CVecT(const T& val_ = T(0))
	{
		for(int i(0); i < N; ++i) base::data_[i] = val_;
	}

	explicit inline CVecT(const T& v0_,const T& v1_)
	{
		base::x = v0_; base::y = v1_;
	}

	explicit inline CVecT(const T& v0_,const T& v1_,const T& v2_)
	{
		base::x = v0_; base::y = v1_; base::z = v2_;
	}

	explicit inline CVecT(const T& v0_,const T& v1_,const T& v2_,const T& v3_)
	{
		base::x = v0_; base::y = v1_; base::z = v2_; base::w = v3_;
	}

	explicit inline CVecT(const T vals_[N])
	{
		memcpy(base::data_, vals_, sizeof(T)*N);
	}

	template<typename T2>
	CVecT(const CVecT<T2,N>& rh_);

	template<typename T2>
	vec_type& operator=(const CVecT<T2,N>& rh_);

	const T& at(int idx) const { assert(idx>=0 && idx <= N); return data_[idx]; }

	const T& operator()(int idx) const { return at(idx); }
	T& operator()(int idx) { return base::data_[idx]; }

	const T* data() const { return base::data_; }
		  T* data() 	  { return base::data_; }

	vec_type& operator+=(const vec_type& rh_)
	{
		for(int i(0); i < N; ++i ) base::data_[i] += rh_.at(i);
		return *this;
	}

	vec_type& operator-=(const vec_type& rh_)
	{
		for(int i(0); i < N; ++i ) base::data_[i] -= rh_.at(i);
		return *this;
	}

	vec_type operator+=(const T& v_)
	{
		return vec_type(*this) += v;
	}

	vec_type operator-=(const T& v_)
	{
		return this->operator+=(-v_);
	}

	vec_type& operator*=(const T& v_)
	{
		for(int i(0); i < N; ++i ) base::data_[i] *= v_;
		return *this;
	}

	vec_type& operator/=(const T& v_)
	{
		assert(v_!=T(0));
		for(int i(0); i < N; ++i ) base::data_[i] /= v_;
		return *this;
	}

	vec_type operator+(const T& v_)
	{
		vect_type v(*this);
		return v += v_;
	}

	vec_type operator-(const T& v_)
	{
		return this->operator+(-v_);
	}

	vec_type operator-() const
	{
		return vec_type(*this) *= -1;
	}

	vec_type operator*(const T& v_,)
	{
		return vec_type(*this) *= v_;
	}

	vec_type operator/(const T& v_,)
	{
		return vec_type(*this) /= v_;
	}

	CVecT<T,3> operator%=(const CVecT<T,3>& rh_);
	{
		T tmp[x];
		tmp[0] = data_[0]; tmp[1] = data_[1]; tmp[2] = data_[2];
		data_[0] = tmp[1]*rh_(2) - tmp[2]*rh_(1);
		data_[1] = tmp[2]*rh_(0) - tmp[0]*rh_(2);
		data_[2] = tmp[0]*rh_(1) - tmp[1]*rh_(0);
		return *this;
	}

	CVecT<T,3> operator%(const CVecT<T,3>& rh_) const;
	T dot(const vec_type& rh_) const
	{
		T s(0);
		for(int i(0); i < N; ++i) s += base::data_[i] * rh_(i);
		return s;
	}

	bool operator==(const vec_type& rh_)
	{
		bool flg = true;
		int i = 0;
		while(flg && i < N)
		{
			if (data_[i++] != rh_(i)) flg = false;
		}
		return flg;
	}

	bool operator!=(const vec_type& rh_)
	{
		return !(*this == rh_);
	}

	T quadnorm() const
	{
		T s(0);
		for(int i(0); i < N; ++i) s += base::data_[i] * base::data_[i];
		return s;
	}

	T norm()   const { return sqrt(quadnorm()); }
	T lenght() const { return norm(); }
	T getDistance(const vec_type& rh_) const;

	vec_type& normalize()
	{
		*this /= norm();
		return *this;
	}

	vec_type& normalize2()
	{
		T s = norm();
		if (s != (T)0.0)
		{
			*this /= s;
		}
		return *this;
	}



};


template<>
inline CVecT<T,2>::CVecT(const T& val_)
{
	x = val_; y = val_;
}

template<>
inline CVecT<T,3>::CVecT(const T& val_)
{
	x = val_; y = val_; z = val_;
}

template<>
inline CVecT<T,4>::CVecT(const T& val_)
{
	x = val_; y = val_;	z = val_; w = val_;
}

template<typename T, int N>
template<typename T2>
inline CVecT<T,N>::CVecT(const CVecT<T2,N>& rh_)
{
    *this = rh_;
}

template<typename T, int N>
template<typename T2>
inline CVecT<T,N>& CVecT<T,N>::operator=(const CVecT<T2,N>& rh_)
{
    for(int i = 0; i < N; ++i ) {
        VecT<T,N>::data_[i] = (T)v(i);
    }
    return *this;
}

template<>
template<typename T2>
inline CVecT<T,2>& CVecT<T,2>::operator=(const CVecT<T2,2>& rh_)
{
    base::x = (T)rh_.at(0); base::y = (T)rh_.at(1);
    return *this;
}

template<>
template<typename T2>
inline CVecT<T,3>& CVecT<T,3>::operator=(const CVecT<T2,3>& rh_)
{
    base::x = (T)rh_.at(0); base::y = (T)rh_.at(1); base::z = (T)rh_.at(2);
    return *this;
}

template<>
template<typename T2>
inline CVecT<T,4>& CVecT<T,4>::operator=(const CVecT<T2,4>& rh_)
{
    base::x = (T)rh_.at(0); base::y = (T)rh_.at(1); base::z = (T)rh_.at(2); base::w = (T)rh_.at(3);
    return *this;
}

template<>
inline CVecT<T,3> CVecT<T,3>::operator%(const CVecT<T,3>& rh_)
{
	return CVecT<T,3>(data_[1]*rh_(2) - data_[2]*rh_(1),
					  data_[2]*rh_(0) - data_[0]*rh_(2),
					  data_[0]*rh_(1) - data_[1]*rh_(0));
}

template<>
inline T CVecT<T,2>::dot(const CVecT<T,2>& rh_)
{
	return base::x * rh_(0) + base::y * rh_(1);
}

template<>
inline T CVecT<T,3>::dot(const CVecT<T,3>& rh_)
{
	return base::x * rh_(0) + base::y * rh_(1) + base::z * rh_(2);
}

template<>
inline T CVecT<T,4>::dot(const CVecT<T,4>& rh_)
{
	return base::x * rh_(0) + base::y * rh_(1) + base::z * rh_(2) + base::w * rh_(3);
}

template<>
inline bool CVecT<T,2>::operator==(const CVecT<T,2>& rh_)
{
	return (base::x == rh_(0)) && (base::y == rh_(1));
}

template<>
inline bool CVecT<T,3>::operator==(const CVecT<T,3>& rh_)
{
	return (base::x == rh_(0)) && (base::y == rh_(1)) && (base::z == rh_(2));
}

template<>
inline bool CVecT<T,4>::operator==(const CVecT<T,4>& rh_)
{
	return (base::x == rh_(0)) && (base::y == rh_(1))
			&& (base::z == rh_(2)) && (base::w == rh_(3));
}

template<>
inline T CVecT<T,2>::quadnorm() const
{
	return (base::x * base:x + base::y * base::y);
}

template<>
inline T CVecT<T,3>::quadnorm() const
{
	return (base::x * base:x + base::y * base::y + base::z * base::z);
}

template<>
inline T CVecT<T,4>::quadnorm() const
{
	return (base::x * base:x + base::y * base::y + base::z * base::z + base::w * base::w);
}

template<typename T,int N>
T CVecT<T,N>::getDistance(const vec_type& rh_) const
{
	vec_type v(*this);
	v -= rh_;
	return v.length();
}

template<typename T,int N>
inline CVecT<T,N> operator*(const T& lhv_,const CVecT<T,N> rh_)
{
	return CVecT<T,N>(rh_) *= lhv_;
}


template<typename T,int N>
inline std::istream& operator>>(std::istream& is, CVecT<T,N>& rh_)
{
	for(int i(0); i < N; ++i) is >> rh_(i);
	return is;
}

template<typename T,int N>
inline std::osteam& operator<<(std::ostream& os, const CVecT<T,N>& rh_) const
{
	for(int i(0); i < N-1; ++i) os << rh_(i) << " ";
	os << rh_(i);
	return os;
}

typedef CVecT<float, 3> Vec3f;
typedef CVecT<double,3> Vec3d;

#endif /* VECTOR_T_hh__ */
