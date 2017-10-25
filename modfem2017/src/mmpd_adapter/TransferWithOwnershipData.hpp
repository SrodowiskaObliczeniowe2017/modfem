#ifndef TRANSFERERWITHOWNERSHIP_DATA_HPP_
#define TRANSFERERWITHOWNERSHIP_DATA_HPP_

#include "TransferData.hpp"

namespace mmpt
{

    /// and also a result of the transfer.
class TransferWithOwnershipResult
    : public TransferResult
{
    public:
    // Update to the parallel mesh.
    // mmpt_owners_map  ovl_nodes;
    // std::vector<char> bnd_vtx_flags;

};

}

#endif // TRANSFERERWITHOWNERSHIP_DATA_HPP_
