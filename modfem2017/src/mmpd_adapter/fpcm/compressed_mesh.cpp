#include "compressed_mesh.hpp"

#include <cmath>

namespace fpcm
{
namespace CompressedMesh
{
int EncodingProcessor(const double* const real_points, const int n_points, __restrict__ PTID * encoded_points, CoordMesh& out_mesh)
{
    mf_check_mem(real_points);
    mf_check_mem(encoded_points);
    // encoding
    const double* p_in=real_points;
    if(encoded_points!=NULL) {
        for(int i=0; i < n_points; ++i,p_in+=3) {
            encoded_points[i] = out_mesh.encode(p_in);
        }
    }
    else {
     for(int i=0; i < n_points; ++i,p_in+=3) {
            out_mesh.encode(p_in);
        }
    }
    // checking
    double tmp[3];
    p_in = real_points;
    for(int i=0; i < n_points; ++i,p_in+=3) {
        out_mesh.getPointCoords(encoded_points[i],tmp);
        mf_test(fabs(tmp[0]-*p_in) < point_accuracy(), "Precision X check failure for %d(PTID=%d) (%lf != %lf)!", i,encoded_points[i], tmp[0],*p_in);
        mf_test(fabs(tmp[1]-*(p_in+1)) < point_accuracy(), "Precision Y check failure for %d(PTID=%d) (%lf != %lf)!",i,encoded_points[i], tmp[1],*(p_in+1));
        mf_test(fabs(tmp[2]-*(p_in+2)) < point_accuracy(), "Precision Z check failure for %d(PTID=%d)(%lf != %lf)!",i,encoded_points[i], tmp[2],*(p_in+2));
    }

    mf_check_mem(encoded_points);

    return n_points;
}

int TransformProcessor(const PTID *points_in, const int n_size, const CoordMesh& from,  const CoordMesh& to, PTID * points_out )
{

    return 0;
}


} // namespace
} // namespace

