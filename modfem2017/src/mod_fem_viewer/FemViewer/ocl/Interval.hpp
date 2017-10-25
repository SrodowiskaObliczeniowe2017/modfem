#ifndef _Interval_hpp__
#define _Interval_hpp__

#include <ElVis/Core/Float.h>
#include <ElVis/Core/Cuda.h>

#if !defined(__CUDACC__)
#include <iostream>
#include <limits>
#endif

namespace FemViewer
{
    template<typename T>
    class Interval
    {
        public:
            Interval() :
                m_low(1.0),
                m_high(-1.0)
            {
            }

            Interval(const T& v) :
                _low(v),
                _high(v)
            {
            }

            ELVIS_DEVICE Interval(const T& low, const T& high) :
                m_low(low),
                m_high(high)
            {
            }

            ELVIS_DEVICE Interval(const Interval& rhs) :
                m_low(rhs.m_low),
                m_high(rhs.m_high)
            {
            }

            ELVIS_DEVICE ~Interval() {}

            ELVIS_DEVICE Interval& operator=(const Interval& rhs)
            {
                m_low = rhs.m_low;
                m_high = rhs.m_high;
                return *this;
            }

            ELVIS_DEVICE void Combine(const Interval<T>& rhs)
            {
                if( rhs.GetLow() < GetLow() ) SetLow(rhs.GetLow());
                if( rhs.GetHigh() > GetHigh() ) SetHigh(rhs.GetHigh());
            }

            ELVIS_DEVICE const T& GetLow() const { return m_low; }
            ELVIS_DEVICE const T& GetHigh() const { return m_high; }
            ELVIS_DEVICE void SetLow(const T& low) { m_low = low; }
            ELVIS_DEVICE void SetHigh(const T& high) { m_high = high; }
            ELVIS_DEVICE void Set(const T& low, const T& high) { SetLow(low); SetHigh(high); }
            ELVIS_DEVICE bool Contains(const T& value) const { return value >= m_low && value <= m_high; }

            ELVIS_DEVICE bool IsEmpty() const { return m_low > m_high; }
            ELVIS_DEVICE T GetWidth() const { return m_high - m_low; }
            ELVIS_DEVICE T GetMidpoint() const { return (m_low + m_high)*MAKE_FLOAT(.5); }
            ELVIS_DEVICE T GetMax() const { return fmaxf(fabsf(m_low), fabsf(m_high)); }

            ELVIS_DEVICE void Union(const Interval<T>& rhs)
            {
                if( !IsEmpty() && !rhs.IsEmpty() )
                {
                    SetLow(fminf(GetLow(), rhs.GetLow()));
                    SetHigh(fmaxf(GetHigh(), rhs.GetHigh()));
                }
                else if( IsEmpty() && !rhs.IsEmpty() )
                {
                    (*this) = rhs;
                }
            }

            ELVIS_DEVICE void operator+=(const Interval<T>& rhs)
            {
                m_low += rhs.m_low;
                m_high += rhs.m_high;
            }

        private:
            T _low;
            T _high;
    };

    template<typename T>
    ELVIS_DEVICE bool Subset(const Interval<T>& a, const Interval<T>& b)
    {
        return b.Contains(a.GetLow()) && b.Contains(a.GetHigh());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> fabs(const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() ) return rhs;
        return Interval<T>(fabs(rhs.GetLow()), fabs(rhs.GetHigh()));
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator+(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return Interval<T>();
        }

        return Interval<T>(lhs.GetLow() + rhs.GetLow(), lhs.GetHigh() + rhs.GetHigh());
    }

    // Assumes 0 is not in rhs.
    template<typename T>
    ELVIS_DEVICE Interval<T> operator/(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return Interval<T>();
        }

        T a = lhs.GetLow();
        T b = lhs.GetHigh();
        T c = rhs.GetLow();
        T d = rhs.GetHigh();
        T v0 = a/c;
        T v1 = a/d;
        T v2 = b/c;
        T v3 = b/d;

        Interval<T> result(fminf(fminf(fminf(v0, v1), v2), v3), fmaxf(fmaxf(fmaxf(v0, v1), v2), v3));
        return result;
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator/(const Interval<T>& lhs, const double& rhs)
    {
        if( lhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        Interval<T> newRhs(rhs, rhs);
        return lhs/newRhs;
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator/(const Interval<T>& lhs, const float& rhs)
    {
        if( lhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        Interval<T> newRhs(rhs, rhs);
        return lhs/newRhs;
    }


    template<typename T>
    ELVIS_DEVICE Interval<T> operator/(const double& lhs, const Interval<T>& rhs)
    {
        if( rhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        Interval<T> newLhs(lhs, lhs);
        return newLhs/rhs;
        //return Interval<T>(lhs/rhs.GetHigh(), lhs/rhs.GetLow());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator/(const float& lhs, const Interval<T>& rhs)
    {
        if( rhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        Interval<T> newLhs(lhs, lhs);
        return newLhs/rhs;
        //return Interval<T>(lhs/rhs.GetHigh(), lhs/rhs.GetLow());
    }

    template<typename T>
    ELVIS_DEVICE void Divide(const Interval<T>& lhs, const Interval<T>& rhs, Interval<T>& out1, Interval<T>& out2)
    {
        out1 = Interval<T>();
        out2 = Interval<T>();

        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return;
        }

        if( rhs.Contains(0.0) )
        {
            if( rhs.GetLow() < 0.0 &&
                rhs.GetHigh() > 0.0 )
            {
                Interval<double> rhs1(-ELVIS_FLOAT_MAX,1.0/rhs.GetLow());
                out1 = lhs*rhs1;

                Interval<double> rhs2(1.0/rhs.GetHigh(),ELVIS_FLOAT_MAX);
                out2 = lhs*rhs2;
            }
            else if( rhs.GetLow() == 0.0 )
            {
                Interval<double> rhs2(1.0/rhs.GetHigh(),ELVIS_FLOAT_MAX);
                out1 = lhs*rhs2;
            }
            else
            {
                Interval<double> rhs1(-ELVIS_FLOAT_MAX,1.0/rhs.GetLow());
                out1 = lhs*rhs1;
            }
        }
        else
        {
            out1 = lhs/rhs;
        }
    }

    template<typename T>
    ELVIS_DEVICE void Divide(const T& lhs, const Interval<T>& rhs, Interval<T>& out1, Interval<T>& out2)
    {
        Interval<T> newLhs(lhs, lhs);
        Divide(newLhs, rhs, out1, out2);
    }

    template<typename T>
    ELVIS_DEVICE void Divide(const Interval<T>& lhs,const T& rhs,  Interval<T>& out1, Interval<T>& out2)
    {
        Interval<T> newRhs(rhs, rhs);
        Divide(lhs, newRhs, out1, out2);
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator+(const Interval<T>& lhs, const double& rhs)
    {
        if( lhs.IsEmpty() )
        {
            return Interval<T>();
        }

        return Interval<T>(lhs.GetLow() + rhs, lhs.GetHigh() + rhs);
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator+(const double& lhs, const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        return Interval<T>(static_cast<T>(lhs) + rhs.GetLow(), static_cast<T>(lhs) + rhs.GetHigh());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator+(const Interval<T>& lhs, const float& rhs)
    {
        if( lhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        return Interval<T>(lhs.GetLow() + static_cast<T>(rhs), lhs.GetHigh() + static_cast<T>(rhs));
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator+(const float& lhs, const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        return Interval<T>(static_cast<T>(lhs) + rhs.GetLow(), static_cast<T>(lhs) + rhs.GetHigh());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator-(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        return Interval<T>(lhs.GetLow() - rhs.GetHigh(), lhs.GetHigh() - rhs.GetLow());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator-(const Interval<T>& lhs, const float& rhs)
    {
        if( lhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        return Interval<T>(lhs.GetLow() - rhs, lhs.GetHigh()-rhs);
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator-(const Interval<T>& lhs, const double& rhs)
    {
        if( lhs.IsEmpty()  )
        {
            return Interval<T>();
        }
        return Interval<T>(lhs.GetLow() - rhs, lhs.GetHigh()-rhs);
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator-(const float& lhs, const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        return Interval<T>(lhs - rhs.GetHigh(), lhs-rhs.GetLow());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator-(const double& lhs, const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        return Interval<T>(lhs - rhs.GetHigh(), lhs-rhs.GetLow());
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator-(const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        Interval<T> result(-rhs.GetHigh(), -rhs.GetLow());
        return result;
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator*(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return Interval<T>();
        }

        T t1 = lhs.GetLow() * rhs.GetLow();
        T t2 = lhs.GetLow() * rhs.GetHigh();
        T t3 = lhs.GetHigh() * rhs.GetLow();
        T t4 = lhs.GetHigh() * rhs.GetHigh();

        T low = fminf(t1, fminf(t2, fminf(t3, t4)));
        T high = fmaxf(t1, fmaxf(t2, fmaxf(t3, t4)));

        return Interval<T>(low, high);
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator*(const Interval<T>& lhs, const double& rhs)
    {
        Interval<T> temp(rhs, rhs);
        return lhs*temp;
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator*(const double& lhs, const Interval<T>& rhs)
    {
        Interval<T> temp(static_cast<T>(lhs), static_cast<T>(lhs));
        return temp*rhs;
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator*(const Interval<T>& lhs, const float& rhs)
    {
        Interval<T> temp(rhs, rhs);
        return lhs*temp;
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> operator*(const float& lhs, const Interval<T>& rhs)
    {
        Interval<T> temp(lhs, lhs);
        return temp*rhs;
    }

    template<typename T>
    ELVIS_DEVICE bool Overlaps(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return false;
        }

        return rhs.GetHigh() >= lhs.GetLow() &&
            rhs.GetLow() <= lhs.GetHigh();
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> Intersection(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return Interval<T>();
        }

        return Interval<T>(fmaxf(lhs.GetLow(), rhs.GetLow()),
            fminf(lhs.GetHigh(), rhs.GetHigh()));
    }

    template<typename T>
    ELVIS_DEVICE Interval<T> Union(const Interval<T>& lhs, const Interval<T>& rhs)
    {
        if( lhs.IsEmpty() || rhs.IsEmpty() )
        {
            return Interval<T>();
        }

        Interval<T> result = lhs;
        result.Union(rhs);
        return result;
    }

    // Bring in the other definitions of exp to prevent this method from hiding them.
    using ::exp;

    template<typename T>
    ELVIS_DEVICE Interval<T> exp(const Interval<T>& rhs)
    {
        if( rhs.IsEmpty() )
        {
            return Interval<T>();
        }
        Interval<T> result(exp(rhs.GetLow()), exp(rhs.GetHigh()));
        return result;
    }

    #if !defined(__CUDACC__)
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const Interval<T>& interval)
    {
        os << "[" << interval.GetLow() << ", " << interval.GetHigh() << "]";
        return os;
    }
    #endif
}


#endif /* _Interval_hpp__
*/

