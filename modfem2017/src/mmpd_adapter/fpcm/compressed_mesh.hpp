#ifndef COMPRESSED_MESH_HPP
#define COMPRESSED_MESH_HPP

/// \file compressed_mesh.hpp
/// \brief compressed_mesh - Implementation of algorithm for compressing and storing mesh in compressed form for fast and memory efficient parallel handling.
/// \author kazimierz.michalik@agh.edu.pl
///

#include <inttypes.h>
#include <limits>
#include <algorithm>

#include <omp.h>

#include "coucal.h"
#include "uth_log.h"

namespace fpcm
{
namespace CompressedMesh
{

typedef int32_t PTID; //point ID
typedef PTID GID;    // globalid

inline double point_accuracy() { return 10e-3; }

struct CoordMesh{
    double  origin[3],d[3];
    PTID binary_oper[3]; //pos_in_line,line_in_layer
    int8_t in_line_digits,line_digits,layer_digits,shift;


    static const PTID digits = std::numeric_limits<PTID>::digits;

    int createMeshBase(const double* const real_points, const int n_points)
    {

        mf_check_mem(real_points);

        {
            double min[3]={std::numeric_limits<double>::max()},
                    max[3]={std::numeric_limits<double>::min()};
            for(int i=0; i < n_points; ++i) {
                for(int j=0;j<3;++j) {
                    if(min[j] > real_points[3*i+j]) {
                        min[j] = real_points[3*i+j];
                    }
                    else if(max[j] < real_points[3*i+j]) {
                        max[j] = real_points[3*i+j];
                    }
                }
            }


            const double span[3] = { max[0] - min[0], max[1] - min[1], max[2] - min[2] } ;
            const int smallest_dim = std::min_element(span,span+3) - span;
    //        const int biggest_dim = std::max_element(span, span+3) - span;
            int ratio[3] = { int(span[0]/span[smallest_dim]), int(span[1]/span[smallest_dim]), int(span[2]/span[smallest_dim]) };
            int ratio_sum = ratio[0]+ratio[1]+ratio[2];

            mf_log_info("Compressed mesh ratios %d : %d : %d",ratio[0],ratio[1],ratio[2]);

            // NOTE: if failing with accuracy in smallest dim increase min_digits
            const int min_digits=6;
            int n_bits = std::numeric_limits<PTID>::digits-(3*min_digits);

    //        // NOTE: we assume, that in each direction at least 5-level adaptaption is posiible (2^5)
    //        const int min_n_adapts_in_dim = 5;
    //        n_bits-=min_n_adapts_in_dim*3;
    //        ratio[smallest_dim]=min_n_adapts_in_dim;
    //        n_bits-=min_n_adapts_in_dim;

            double avg = double(n_bits)/double(ratio_sum);

            mf_log_info("Avg= %lf",avg);

            for(int i=0; i<3 ;++i) {
                ratio[i] = min_digits + double(ratio[i])*avg;
    //                    + (min_n_adapts_in_dim-1);
            }

            ratio_sum = ratio[0]+ratio[1]+ratio[2];
            assert( ratio_sum <= std::numeric_limits<PTID>::digits );

            while(ratio_sum < std::numeric_limits<PTID>::digits) {
                ++ratio[smallest_dim];
                ++ratio_sum;
            }

            for(int i=0; i < 3; ++i) {
                origin[i] = min[i];
                d[i] = span[i] / ((1<<ratio[i])-1);
            }

            setLengths(ratio[0],ratio[1]);


            mf_log_info("Min= %lf, %lf, %lf",min[0],min[1],min[2]);
            mf_log_info("Max= %lf, %lf, %lf",max[0],max[1],max[2]);
            mf_log_info("Span= %lf, %lf, %lf",span[0],span[1],span[2]);
            mf_log_info("Smallest dim=%d",smallest_dim);
            mf_log_info("Compressed mesh ratios %d : %d : %d",ratio[0],ratio[1],ratio[2]);
            mf_log_info("pt per line=%d, lines=%d, layers=%d",in_line_digits,line_digits,layer_digits);
            mf_log_info("binary_oper[3] = %d, %d, %d",binary_oper[0],binary_oper[1],binary_oper[2]);
            mf_log_info("d[3]=%lf, %lf, %lf",d[0], d[1], d[2]);

            mf_check(min[0] == origin[0], "Error in encoding processor!");
            mf_check(min[1] == origin[1], "Error in encoding processor!");
            mf_check(min[2] == origin[2], "Error in encoding processor!");
            mf_check(max[0] == origin[0]+d[0]*((1<<ratio[0])-1), "Error in encoding processor!");
            mf_check(max[1] == origin[1]+d[1]*((1<<ratio[1])-1), "Error in encoding processor!");
            mf_check(max[2] == origin[2]+d[2]*((1<<ratio[2])-1), "Error in encoding processor!");
        }
    }

    void setLengths(int lineLengthAs2Pow, int layerLengthAs2Pow)
    {
        in_line_digits = lineLengthAs2Pow;
        line_digits = layerLengthAs2Pow;
        layer_digits = digits - (lineLengthAs2Pow+layerLengthAs2Pow);
        shift = in_line_digits + line_digits;

        binary_oper[0]=0;
        binary_oper[1]=0;
        binary_oper[2]=0;

        for(int i=0;i<in_line_digits; ++i) {
            binary_oper[0]|=(1<<i);
        }

        for(int i=0;i<line_digits; ++i) {
            binary_oper[1]|=(1<<i);
        }
        binary_oper[1] = binary_oper[1] << lineLengthAs2Pow;

        for(int i=0;i<layer_digits; ++i) {
            binary_oper[2]|=(1<<i);
        }
        binary_oper[2] = (binary_oper[2] << (lineLengthAs2Pow+layerLengthAs2Pow));
    }

    inline void getPointCoords(const PTID pt,   double* coords) const {
        coords[0]=origin[0] + d[0]*(pt & (binary_oper[0])) ;
        coords[1]=origin[1] + d[1]*((pt & (binary_oper[1]))>>in_line_digits) ;
        coords[2]=origin[2] + d[2]*(pt >> (in_line_digits+line_digits));

    }

    inline PTID encode(const double* const coords) const  {

        assert(PTID((coords[0]-origin[0]) / d[0]) < (1<<in_line_digits));
        assert(PTID((coords[1]-origin[1]) / d[2])>>in_line_digits < (1<<line_digits));
        assert(PTID((coords[2]-origin[2]) / d[2])>>(in_line_digits+line_digits) < (1<<layer_digits));

        return    PTID((coords[0]-origin[0]) / d[0]) // pos in line
                | (PTID((coords[1]-origin[1]) / d[1]) << in_line_digits) // which line
                | (PTID((coords[2]-origin[2]) / d[2]) << (line_digits+in_line_digits)); // which layer
    }

    inline PTID middle(const PTID p1, const PTID p2) const {
        return PTID((p1 & (binary_oper[0]))/2 + (p2 & (binary_oper[0]))/2)
                + PTID((p1 & (binary_oper[1]))/2 + (p2 & (binary_oper[1]))/2)
                + PTID((p1 & (binary_oper[2]))/2 + (p2 & (binary_oper[2]))/2);
    }

    inline double X(const PTID pt) const {
        return origin[0] + d[0]*double(pt & binary_oper[0]);
    }

    inline double Y(const PTID pt) const {
        return origin[1] + d[1]*((pt>>in_line_digits) & (binary_oper[0]));
    }

    inline double Z(const PTID pt) const {
        return  origin[2] + d[2]*(pt >> shift);
    }

};

int EncodingProcessor(const double* const real_points, const int n_points, __restrict__ PTID * encoded_points, const CoordMesh& out_mesh);

int TransformProcessor(const PTID *points_in, const int n_size, const CoordMesh& from,  const CoordMesh& to, PTID * points_out );

}
}
#endif
