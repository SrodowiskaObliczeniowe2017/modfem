#ifndef _GRAPHIC_ELEM_H_
#define _GRAPHIC_ELEM_H_

#include "fv_float.h"
#include "defs.h"
#include "MathHelper.h"
#include "Interval.hpp"
#include "Color.h"
//#include "Fire.h"

//#include <sstream>
#include <vector>
#include <map>
#include <string.h>

namespace FemViewer {
using fvmath::Interval;
typedef enum {
	eDensity,
	eRed,
	eGreen,
	eBlue
} tf_chanel_t;

class TransferFunction
{
public:
	TransferFunction() {}

	void Init() {
		memset(this,0x0,sizeof(*this));
	}
#define SET(start,length) \
	begin = start; \
	end = start + (length-1); \
	count = length; break

	void GetBreakpointsForChannel(tf_chanel_t channel, ScalarValueType*& begin, ScalarValueType*& end, int& count) const
	{
		switch(channel)
		{
		case eDensity:
			SET(m_dbrks_ptr,m_numd_bkrs);
		case eRed:
			SET(m_rbrks_ptr,m_numr_brks);
		case eGreen:
			SET(m_gbrks_ptr,m_numb_brks);
		case eBlue:
			SET(m_bbrks_ptr,m_numb_brks);
		default:
			SET(m_dbrks_ptr,m_numd_bkrs);
		}
	}

	void GetValuesForChannel(tf_chanel_t channel, ScalarValueType*& begin, ScalarValueType*& end, int& count) const
	{
		switch(channel)
		{
		case eDensity:
			SET(m_dvals_ptr,m_numd_bkrs);
		case eRed:
			SET(m_rvals_ptr,m_numr_brks);
		case eGreen:
			SET(m_gbrks_ptr,m_numb_brks);
		case eBlue:
			SET(m_bbrks_ptr,m_numb_brks);
		default:
			SET(m_dbrks_ptr,m_numd_bkrs);
		}
	}

#undef SET

	bool ColorContainsAtLeastOneBreakpoint(const Interval<ScalarValueType>& range) const
	{
		return RangeContainsAtLeastOneBreakpoint(eGreen, range) ||
				RangeContainsAtLeastOneBreakpoint(eRed, range) ||
				RangeContainsAtLeastOneBreakpoint(eBlue, range);
	}

	bool RangeContainsAtLeastOneBreakpoint(tf_chanel_t channel, const Interval<ScalarValueType>& range) const
    {
		ScalarValueType* breakpointStart = 0;
		ScalarValueType* breakpointEnd = 0;
		int numBreakpoints = 0;
		GetBreakpointsForChannel(channel, breakpointStart, breakpointEnd, numBreakpoints);

		for(int i = 0; i < numBreakpoints; ++i)
		{
			if( range.Contains(breakpointStart[i]) )
			{
				return true;
			}
		}

		return false;
    }

	bool IntersectsRange(tf_chanel_t channel, const Interval<ScalarValueType>& range) const
	{
		ScalarValueType* brs_ptr = nullptr;
		ScalarValueType* bre_ptr = nullptr;
		int count = 0;
		GetBreakpointsForChannel(channel, brs_ptr, bre_ptr, count);

		return (*bre_ptr) >= range.GetLow() && (*brs_ptr) <= range.GetHigh();
	}

//            ELVIS_DEVICE Interval<ScalarValueType> Sample(TransferFunctionChannel channel, const Interval<ScalarValueType>& range) const
//            {
//                Interval<ScalarValueType> result;

//                ScalarValueType* breakpointStart = 0;
//                ScalarValueType* breakpointEnd = 0;
//                int numBreakpoints = 0;
//                GetBreakpointsForChannel(channel, breakpointStart, breakpointEnd, numBreakpoints);

//                // Find the location of the breakpoint.
//                int index0 = 0;
//                int index1 = 0;
//                for(int i = 0; i < numBreakpoints-1; ++i)
//                {
//                    bool test0 = range.GetLow() >= breakpointStart[i];
//                    test0 &= (range.GetLow() < breakpointStart[i+1]);

//                    if( test0 ) index0 = i+1;

//                    bool test1 = range.GetHigh() >= breakpointStart[i];
//                    test1 &= (range.GetHigh() < breakpointStart[i+1]);

//                    if( test1 ) index1 = i+1;
//                }
//                if( range.GetLow() > breakpointStart[numBreakpoints-1] )
//                {
//                    index0 = numBreakpoints;
//                }
//                if( range.GetHigh() > breakpointStart[numBreakpoints-1] )
//                {
//                    index1 = numBreakpoints;
//                }

//                ScalarValueType* valueBegin = 0;
//                ScalarValueType* valueEnd = 0;
//                int count = 0;
//                GetValuesForChannel(channel, valueBegin, valueEnd, count);

//                if( index0 == 0 )
//                {
//                    result.SetLow(valueBegin[0]);
//                }
//                else if( index0 == numBreakpoints )
//                {
//                    result.SetLow(valueBegin[index0-1]);
//                }
//                else
//                {
//                    ScalarValueType v0 = valueBegin[index0-1];
//                    ScalarValueType v1 = valueBegin[index0];

//                    ScalarValueType s0 = breakpointStart[index0-1];
//                    ScalarValueType s1 = breakpointStart[index0];

//                    if( v0 == v1 )
//                    {
//                        result.SetLow(v0);
//                    }
//                    else
//                    {
//                        ScalarValueType scale = MAKE_FLOAT(1.0)/(s1-s0);
//                        result.SetLow(scale*v1*(range.GetLow()-s0) + scale*v0*(s1-range.GetLow()));
//                    }
//                }


//                if( index1 == 0 )
//                {
//                    result.SetHigh(valueBegin[0]);
//                }
//                else if( index1 == numBreakpoints )
//                {
//                    result.SetHigh(valueBegin[index1-1]);
//                }
//                else
//                {
//                    ScalarValueType v0 = valueBegin[index1-1];
//                    ScalarValueType v1 = valueBegin[index1];

//                    ScalarValueType s0 = breakpointStart[index1-1];
//                    ScalarValueType s1 = breakpointStart[index1];

//                    if( v0 == v1 )
//                    {
//                        result.SetHigh(v0);
//                    }
//                    else
//                    {
//                        ScalarValueType scale = MAKE_FLOAT(1.0)/(s1-s0);
//                        result.SetHigh(scale*v1*(range.GetHigh()-s0) + scale*v0*(s1-range.GetHigh()));
//                    }
//                }


//                for(int i = index0; i < index1; ++i)
//                {
//                    if( valueBegin[i] < result.GetLow() )
//                    {
//                        result.SetLow(valueBegin[i]);
//                    }
//                    if( valueBegin[i] > result.GetHigh() )
//                    {
//                        result.SetHigh(valueBegin[i]);
//                    }
//                }
//                return result;
//            }

	// Not as optimal as the above version is supposed to be, but it appears it has a bug, and I
	// need something correct right now.
	Interval<ScalarValueType> Sample(tf_chanel_t channel, const Interval<ScalarValueType>& range) const
    {
		Interval<ScalarValueType> result;
		ScalarValueType left = Sample(channel, range.GetLow());
		ScalarValueType right = Sample(channel, range.GetHigh());
		result.Set(fv_min(left, right), fv_max(left, right));

		ScalarValueType* brs_ptr = nullptr;
		ScalarValueType* bre_ptr = nullptr;
		int numBreakpoints = 0;
		GetBreakpointsForChannel(channel, brs_ptr, bre_ptr, numBreakpoints);
		ScalarValueType* vals_ptr = nullptr;
		ScalarValueType* vale_ptr = nullptr;
		int count = 0;
		GetValuesForChannel(channel, vals_ptr, vals_ptr, count);

		for(int i = 0; i < numBreakpoints; ++i)
		{
			if( range.Contains(brs_ptr[i]) )
			{
				result.SetLow(fv_min(result.GetLow(), vals_ptr[i]));
				result.SetHigh(fv_max(result.GetHigh(), vals_ptr[i]));
			}
		}
		return result;
    }

	ScalarValueType Sample(tf_chanel_t channel, const ScalarValueType& s) const
	{
		ScalarValueType* brs_ptr = nullptr;
		ScalarValueType* bre_ptr = nullptr;
		int numBreakpoints = 0;
		GetBreakpointsForChannel(channel, brs_ptr, bre_ptr, numBreakpoints);

		// Find the location of the breakpoint.
		int index = 0;
		for(int i = 0; i < numBreakpoints-1; ++i)
		{
			bool test = s >= brs_ptr[i];
			test &= (s < brs_ptr[i+1]);

			if (test) index = i+1;
		}

		if( s >= brs_ptr[numBreakpoints-1] )
		{
			index = numBreakpoints;
		}

		ScalarValueType* vals_ptr = nullptr;
		ScalarValueType* vale_ptr = nullptr;
		int count = 0;
		GetValuesForChannel(channel, vals_ptr, vale_ptr, count);

		ScalarValueType result = FLOAT_CONVERT(0.0);
		if( index == 0 )
		{
			result = vals_ptr[0];
		}
		else if( index == numBreakpoints )
		{
			result = vals_ptr[index-1];
		}
		else
		{
			ScalarValueType v0 = vals_ptr[index-1];
			ScalarValueType v1 = vals_ptr[index];

			ScalarValueType s0 = brs_ptr[index-1];
			ScalarValueType s1 = brs_ptr[index];

			if( v0 == v1 )
			{
				result = v0;
			}
			else
			{
				ScalarValueType scale = FLOAT_CONVERT(1.0)/(s1-s0);
				result = scale*v1*(s-s0) + scale*v0*(s1-s);
			}
        }

		return result;
	}

    ScalarValueType3 SampleColor(const ScalarValueType& s) const
    {
    	ScalarValueType3 result;
    	result.x = Sample(eRed, s);
    	result.y = Sample(eGreen, s);
    	result.z = Sample(eBlue, s);

    	return result;
    }

    ScalarValueType GetMaxValue(tf_chanel_t channel) const
    {
    	ScalarValueType* begin = 0;
    	ScalarValueType* end = 0;
    	int count = 0;
    	GetValuesForChannel(channel, begin, end, count);
    	ScalarValueType result = FLOAT_CONVERT(0.0);

    	for(int i = 0; i < count; ++i)
    	{
    		result = fv_max(begin[i], result);
    	}
    	return result;
    }

    ScalarValueType GetMaxValue(tf_chanel_t channel, const Interval<ScalarValueType>& range) const
    {
    	ScalarValueType* valueBegin = 0;
    	ScalarValueType* valueEnd = 0;
    	int count = 0;
    	GetValuesForChannel(channel, valueBegin, valueEnd, count);

    	ScalarValueType* breakpointBegin = 0;
    	ScalarValueType* breakpointEnd = 0;
    	GetBreakpointsForChannel(channel, breakpointBegin, breakpointEnd, count);

    	ScalarValueType result = FLOAT_CONVERT(0.0);

    	for(int i = 0; i < count; ++i)
    	{
    		if( range.Contains(breakpointBegin[i]) )
    		{
    			result = fmaxf(valueBegin[i], result);
    		}
    	}
    	result = fv_max(Sample(channel, range.GetLow()), result);
    	result = fv_max(Sample(channel, range.GetHigh()), result);
    	return result;
    }

    ScalarValueType*& DensityBreakpoints() { return m_dbrks_ptr; }
    ScalarValueType*& RedBreakpoints() { return m_rbrks_ptr; }
    ScalarValueType*& GreenBreakpoints() { return m_gbrks_ptr; }
    ScalarValueType*& BlueBreakpoints() { return m_bbrks_ptr; }

    ScalarValueType*& DensityValues() { return m_dvals_ptr; }
    ScalarValueType*& RedValues() { return m_rvals_ptr; }
    ScalarValueType*& GreenValues() { return m_gvals_ptr; }
    ScalarValueType*& BlueValues() { return m_bvals_ptr; }

    int& NumDensityBreakpoints() { return m_numd_bkrs; }
    int& NumRedBreakpoints() { return m_numr_brks; }
    int& NumGreenBreakpoints() { return m_numg_brks; }
    int& NumBlueBreakpoints() { return m_numb_brks; }

private:
    TransferFunction(const TransferFunction&);
    TransferFunction& operator=(const TransferFunction&);

    ScalarValueType* m_dbrks_ptr;
    ScalarValueType* m_rbrks_ptr;
    ScalarValueType* m_gbrks_ptr;
    ScalarValueType* m_bbrks_ptr;

    ScalarValueType* m_dvals_ptr;
    ScalarValueType* m_rvals_ptr;
    ScalarValueType* m_gvals_ptr;
    ScalarValueType* m_bvals_ptr;

    int m_numd_bkrs;
    int m_numr_brks;
    int m_numg_brks;
    int m_numb_brks;

};


struct Breakpoint
{
	ColorRGB Col;
	mfvFloat_t Density;
	mfvFloat_t Scalar;
};

class HostTransferFunction
{
public:
	HostTransferFunction();
	TransferFunction GetOptixObject();

    void SetBreakpoint(double s, const ColorRGBA& c);
    void SetBreakpoint(double s, const ColorRGBA& c, const mfvFloat_t& density);

    bool IsValid() const
    {
    	return m_breakpoints.size() >= 2;
    }

    const std::map<double, Breakpoint>& GetBreakpoints() const { return m_breakpoints; }

    void Clear();

    //void CopyToOptix(optixu::Context context, OptiXBuffer<ElVisFloat>& buffer, OptiXBuffer<ElVisFloat>& values, TransferFunctionChannel channel);

    bool& Dirty() { return m_dirty; }

private:
    HostTransferFunction(const HostTransferFunction&);
    HostTransferFunction& operator=(const HostTransferFunction&);

    void UpdateBreakpoints(std::map<double, double>& container);

    void SynchronizeDeviceIfNeeded();
    void SynchronizeOptiXIfNeeded();
    void FreeDeviceMemory();
    void AllocateDeviceMemory();
    void CopyToDeviceMemory();

    std::map<double, Breakpoint> m_breakpoints;

    TransferFunction m_localDeviceTransferFunction;
    bool m_dirty;
};

} // end namespace FemViewer

#endif /* _GAPHIC_ELEM_H_
  */
