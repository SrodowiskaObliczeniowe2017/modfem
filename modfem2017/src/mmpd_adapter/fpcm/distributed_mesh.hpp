#ifndef _DISTRIBUTED_MESH_HPP_
#define _DISTRIBUTED_MESH_HPP_

#include <mpi.h>


#include "compressed_mesh.hpp"

namespace fpcm
{
namespace DistributedMesh
{
using namespace fpcm::CompressedMesh;

/// Local (in-process) representation of distributed mesh.
/// In keeps infomration about topology of local mesh AND its' direct neighbours
/// in other processes.

const int MAX_EL_NEIG = 5;
const int NO_NEIG = 0;

struct Element
{
    GID neigs[MAX_EL_NEIG];
};

enum MeshEntityType { EVertex=0, EEdge, EFace, EElement, EN_MeshEntityTypes };


/// \brief DistributedMesh is responsible for logic behind parallel mesh operations.
/// It decides which mesh entities are transfered and between which processes.
/// Thus this class is responsible for mesh entities recognition, migration orders and others.
struct DistributedMesh
{
    typedef coucal HashMap;

    // Unique global identyfication of all mesh entities
    CoordMesh global_mesh_base;

    // Interprocess face neigbouring info.
    HashMap extern_faces_neighs; //keeps GID of foregin elements connecting with our external faces

    // Local ownership info (inc. overlap mesh entities)
    HashMap foreign_owners[EN_MeshEntityTypes]; //GID -> owning proc id for all local and neighbouring mesh entities

    /// [GID -> neighs GIDs], elements that has no BC or foreign neigh (not stored locally), not affected by ownership.
    HashMap elements_core,
    /// [GID -> neighs GIDs], elements that are stored locally and has at least one foreign neigh, or BC neigh, not affected by ownership.
            elements_boundary;

    /// [GID -> owner proc id], should affect only 'elements_boundary'.
    HashMap overlap_elements_owners;


    HashMap boundary_vts; // only boundary vertices may have foreign owner GID -> owner

    void Init(int* elem_GIDs,  //out: local element ids in reference to given numeration
              const double* const vertices_coordinates = NULL,
              const int n_vertices = 0,
              const double* const elems_center_pts = NULL,
              int* elem_neigs = NULL, // in  [MAX_EL_NEIG*n_elems] in: local numeration, out: GIDs numbers
              const int* elem_domain = NULL, // which element to which subdomain
              const int n_elems = 0
              ) {

        /// Initializing hashes.
        extern_faces_neighs = coucal_new(0);
        elements_core = coucal_new(0);
        elements_boundary = coucal_new(0);

        for(int i=0; i < EN_MeshEntityTypes; ++i) {
            foreign_owners[i] = coucal_new(0);
        }

        Clear();

        /// Computing and distributing global mesh base.
        if(MPI_Comm_rank(MPI_COMM_WORLD,NULL) == 0) {
            mf_check_mem(vertices_coordinates);

             global_mesh_base.createMeshBase(vertices_coordinates, n_vertices);
        }

        MPI_Bcast(global_mesh_base.origin, 3, MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(global_mesh_base.d, 3, MPI_DOUBLE,0,MPI_COMM_WORLD);

        // assuming PTID == MPI_INT
        MPI_Bcast(global_mesh_base.binary_oper, 3, MPI_INT,0,MPI_COMM_WORLD);

        MPI_Bcast(& global_mesh_base.in_line_digits, 1, MPI_BYTE,0,MPI_COMM_WORLD);
        MPI_Bcast(& global_mesh_base.line_digits, 1, MPI_BYTE,0,MPI_COMM_WORLD);
        MPI_Bcast(& global_mesh_base.layer_digits, 1, MPI_BYTE,0,MPI_COMM_WORLD);
        MPI_Bcast(& global_mesh_base.shift, 1, MPI_BYTE,0,MPI_COMM_WORLD);


        /// Converting neigs numbers into GIDs.
        for(int e=0; e < n_elems; ++e) {
            elem_GIDs[e] = global_mesh_base.encode( & elems_center_pts[3*e] );
            mf_check(elem_GIDs[e] != NO_NEIG,"Incorrect GID for %d!",e);
        }

        for(int n=0; n < MAX_EL_NEIG*n_elems; ++n) {
            if(elem_neigs[n] != NO_NEIG) {
                elem_neigs[n]=elem_GIDs[ elem_neigs[n] ];
            }
        }


        /// Distributing elements.
        for(int e=0; e < n_elems; ++e) {
            Element *p_el = new Element();

            //for(int n=0; n < MAX_EL_NEIG; ++n) {
            //
            //}

            //coucal_add(elements,global_mesh_base.encode(& elems_center_pts[3*e]),);
        }

        // barrier of distributed initialization
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void Clear() {
        for(int i=0; i < EN_MeshEntityTypes; ++ i) {
            if(coucal_created(foreign_owners[i]) ) {
                coucal_delete(& foreign_owners[i]);
                foreign_owners[i] = coucal_new(0);
            }
        }

        if(coucal_created(extern_faces_neighs) ) {
            coucal_delete(& extern_faces_neighs);
            extern_faces_neighs = coucal_new(0);
        }

    }




    DistributedMesh() {
    }





//    // are those needed? pergaps it will be easier to compute them?
//    // it should be moved to mmpd_adapter from fpcm!
//    HashMap ve_GIDs; // local number -> GID
//    HashMap ed_GIDs; // local number -> GID
//    HashMap fa_GIDs; // local number -> GID
//    HashMap el_GIDs; // local number -> GID

    template<MeshEntityType TType>
    int getOwner(const GID id) const {
         int *p_owner=NULL;
         coucal_read_pvoid(foreign_owners[TType],&id,(void**) &p_owner);
         return p_owner==NULL ? MPI_Comm_rank(MPI_COMM_WORLD,NULL) : *p_owner;
    }

    template<MeshEntityType TType>
    void setOwner(const GID id, const int owner) {
         coucal_add(foreign_owners[TType],id,owner);
    }

    void DomainDecompositionProcessor();
    void OverlapProcessor();
    void AdaptationProcessor();

};

}
} //!namespace


#endif // _DISTRIBUTED_MESH_HPP_
