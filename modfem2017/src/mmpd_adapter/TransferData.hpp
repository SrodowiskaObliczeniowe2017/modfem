#ifndef TRANSFERER_DATA_HPP_
#define TRANSFERER_DATA_HPP_

#include <vector>
#include <map>

#include "EntitiesOwnership.hpp"


namespace mmpt{

enum TRANSFER_POLICY {
    TRANSFER_MOVE,
    TRANSFER_COPY,
    TRANSFER_MOVE_RR // Response Required
};


/// Transfer vertices
/// [ [coord1,coord2, coord3] ,....]

/// Transfer element:
/// [vertices]

struct TransferElement {
    int type;
    int vertices[MMC_MAXELVNO]; // ids of transfer vertices
    int BCs[MMC_MAXELFAC];
    int id_at_prev_owner;
    int groupID;

    TransferElement() ;
    TransferElement(const int t,const int *v,const int* bc,const int old_id, const int gID)
        : type(t), id_at_prev_owner(old_id), groupID(gID)
    {
        if(type == MMC_PRISM) {
            memcpy(vertices,v,sizeof(int) * MMC_MAXELVNO);
            memcpy(BCs,bc,sizeof(int)* MMC_MAXELFAC);
        }
        else {
            memcpy(vertices,v,sizeof(int) * 4);
            memcpy(BCs,bc,sizeof(int)* 4);
        }
    }
};

//struct MeshPart
//{
//    std::vector<double>     vertices_coords;
//    std::vector<TransferElement> elems;
//};

///
/// \brief The TransferOrder class is responsible encapsulating a
/// single transfer information about source, destination, affected mesh entities.
///
class TransferOrder
{
    // source+dest+policy + Elements.size() = 4
    static const int n_ints_to_serialize_empty=4;
public:
    TransferOrder() : Source_id(-1),Destination_id(-1),t_policy(TRANSFER_COPY)
    {}
    TransferOrder(const TransferOrder& o)
        : Source_id(o.Source_id),Destination_id(o.Destination_id),
          t_policy(o.t_policy),Elements(o.Elements)
    {}
    TransferOrder(const int Sourceproc_id, const int Destproc_id,
                  const TRANSFER_POLICY Policy,
                  const int* elements = NULL, const int n_elements=0);


    int Source_id,Destination_id;
    TRANSFER_POLICY t_policy;
    std::vector<int> Elements;

    /// \return number of int's need for serialization of TransferOrder
    int n_ints_to_serialize() const {
        return n_ints_to_serialize_empty+Elements.size();
    }

    /// \return pointer to the 'next after readed' int
    int* serialize(int * buffer) const {
        buffer[0] = Source_id;
        buffer[1] = Destination_id;
        buffer[2] = t_policy;
        buffer[3] = Elements.size();
        memcpy(buffer+n_ints_to_serialize_empty, Elements.data(), Elements.size() * sizeof(int) );
        return buffer+n_ints_to_serialize();
    }

    /// \return pointer to the 'next after readed' int
    int* deserialize(int * buffer) {
        Source_id      = buffer[0];
        Destination_id = buffer[1];
        t_policy       = static_cast<TRANSFER_POLICY>(buffer[2]);
        Elements.resize(buffer[3]);
        memcpy(Elements.data(), buffer+n_ints_to_serialize_empty, Elements.size() * sizeof(int) );
        return buffer+n_ints_to_serialize();
    }

    std::vector<TransferOrder>  relative_transfer_orders;
};


bool    operator==(const TransferOrder& first, const TransferOrder& second);

enum FaceSubdomainSide {
    FaceSide0=0,
    FaceSide1=1,
    FaceSideBoth=2,
    FaceSideNone=-1
};

/// and also a result of the transfer.
struct TransferResult
{
    TransferResult() : wasApplied(false) {}
    // Update to the parallel mesh.
    // mmpt_owners_map  ovl_nodes;
    // std::vector<char> bnd_vtx_flags;


    // Update to the sequential mesh.
    typedef std::map<int,FaceSubdomainSide> mmpt_faceSideMap;
    mmpt_faceSideMap    new_sub_bdn_faces;
    std::vector<double> new_vertices_coords;
    std::vector<int> new_vertices_transfered_ids;
    std::vector<TransferElement> new_elems;
    mutable bool wasApplied;

    std::vector<EntOwn>      transfer_vertices_ownerships;
    std::vector<int>         ownerships2vtx_trans_id;
};

}

#endif // TRANSFERER_DATA_HPP_
