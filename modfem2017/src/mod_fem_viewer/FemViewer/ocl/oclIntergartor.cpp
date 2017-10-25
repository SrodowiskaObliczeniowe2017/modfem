#include <CL/cl.h>
#include "MathHelper.h"

typedef float3 mfvPoint;

// Estimate a range of scalar field
int EstimateRangeCuda(uint elemId,
		              uint elemType,
		              int  fieldId,
                      const mfvPoint& p0,
                      const mfvPoint& p1,
                      ElVis::Interval<ElVisFloat>& result)
{
    ReferencePoint r0, r1;
    ElVisError e0 = ConvertWorldToReferenceSpaceCuda(elementId, elementType, p0, ElVis::eReferencePointIsInvalid, r0);
    ElVisError e1 = ConvertWorldToReferenceSpaceCuda(elementId, elementType, p1, ElVis::eReferencePointIsInvalid, r1);

    IntervalPoint referenceInterval(r0, r1);
    IntervalPoint worldInterval(p0, p1);
    ElVisError e2 = SampleScalarFieldAtReferencePointCuda(elementId, elementType, fieldId,
                                                          worldInterval, referenceInterval, result);

    return eNoError;
}

