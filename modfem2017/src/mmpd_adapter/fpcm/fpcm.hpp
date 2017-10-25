#ifndef FCPM_HPP
#define FCPM_HPP

/// \file fcpm.hpp
/// \brief Fast Compressed Parallel Mesh
/// \author kazimierz.michalik@agh.edu.pl
///

#include "es.hpp"
#include "compressed_mesh.hpp"
#include "dbg.h"

using fpcm::CompressedMesh::PTID;
using namespace fpcm::ES;

namespace fpcm
{

template<int8_t Tmax_el_verts,int8_t Tmax_fa_verts>
class FastParallelCompressedMesh
{
public:
    typedef CompressedMesh::PTID PTID;
    typedef ES::GUID GUID;

    static const int8_t max_el_verts = Tmax_el_verts;
    static const int8_t max_fa_verts = Tmax_fa_verts;
    typedef Verticable<Tmax_el_verts> ElemVerticable;
    typedef Verticable<Tmax_fa_verts> FaceVerticable;
    typedef Verticable<2>             EdgeVerticable;
    // Storages
    // Mesh vertices (in compressed form)
    CompressedMesh::CoordMesh compressed_mesh;
    std::vector<PTID> used_vertices;

    // Mesh edges.
    ES::GUID n_edges;
    std::vector<EdgeVerticable > init_edges;

    // Mesh faces.
    ES::GUID n_faces;
    std::vector<FaceVerticable> init_faces;
    std::vector<Neigborable> init_faces_neigs;

    // Mesh elements.
    ES::GUID n_elements;
    std::vector<ElemVerticable> init_elems;
    std::vector<Materialable>   init_elems_materials;

    int8_t max_gen,max_gen_diff;


    int encodePoints(const double* const real_points, const int n_points)
    {
        used_vertices.resize(n_points);
        return CompressedMesh::EncodingProcessor(real_points,n_points,
                                                 used_vertices.data(),
                                                 compressed_mesh);
    }

    PTID*    points_begin() { return used_vertices.data(); }
    /// NOTE: this is 'one-past-the-end' iterator.
    /// Dereferencing this pointer has undefined behaviour.
    const PTID*    points_end() const {return used_vertices.data()
                + used_vertices.size(); }
    inline void     coords(const PTID id,double* coords_ar) const
    {
        return compressed_mesh.getPointCoords(id,coords_ar);
    }
    PTID        points_size() const { return used_vertices.size(); }

private:
    void defaults() {
        n_edges=0;
        n_faces=0;
        n_elements=0;
        max_gen=1;
        max_gen_diff=1;
    }
public:
    void clear() {
        defaults();
        used_vertices.clear();
        init_edges.clear();
        init_faces.clear();
        init_faces_neigs.clear();
        init_elems.clear();
        init_elems_materials.clear();
    }

    ~FastParallelCompressedMesh() {
        clear();
    }


    int in_memory_size() const {
        return sizeof(compressed_mesh) +
                container_in_memory_size(used_vertices) +
                container_in_memory_size(init_edges) +
                container_in_memory_size(init_faces) +
                container_in_memory_size(init_faces_neigs) +
                container_in_memory_size(init_elems) +
                container_in_memory_size(init_elems_materials);
    }

    // Processors


};
}

#endif
