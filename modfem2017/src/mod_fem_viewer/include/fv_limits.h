#ifndef _FV_LIMITS_H_
#define _FV_LIMITS_H_

#if _MSC_VER < 1600 

// Auto-enable limits for visual c++
#if (defined(WIN32) && defined(_MSC_VER))
#define FV_MISSING_LIMITS_HEADER_FILE
#define FV_MISSING_STD_MIN_MAX
#undef min
#undef max
#endif

#ifdef FV_MISSING_STD_MIN_MAX
#undef min
#undef max

namespace std {
	
	template <class T> T min(T a,T b) 
	{ return( a < b ? a : b); };

	template <class T> T max(T a,T b) 
	{ return( a > b ? a : b); };
};

#endif //FV_MISSING_STD_MIN_MAX

#ifdef FV_MISSING_LIMITS_HEADER_FILE

#ifdef FV_MISSING_CFLOAT_HEADER_FILE
#include <float.h>
#else // #ifdef FV_MISSING_CFLOAT_HEADER_FILE
#include <cfloat>
#endif // #ifdef FV_MISSING_CFLOAT_HEADER_FILE

namespace std {
	template <class T>
	struct numeric_limits {};

	template <>
	struct numeric_limits<float>
	{
		static double min() throw() 
		{ return FLT_MIN; }
    
		static double max() throw() 
		{ return FLT_MAX;}
	};

	template <>
	struct numeric_limits<double>
	{
		static double min() throw() 
		{ return FLT_MIN; }
    
		static double max() throw() 
		{ return FLT_MAX;}
  };
}

#else // #ifdef FV_MISSING_LIMITS_HEADER_FILE

#include <limits>

#endif // #ifdef FV_MISSING_LIMITS_HEADER_FILE

#endif

#endif /* _FV_LIMITS_H_
*/
