#ifndef COMPOSITE_TRANSFER_ORDER_HPP_
#define COMPOSITE_TRANSFER_ORDER_HPP_

#include "TransferWithOwnershipData.hpp"

namespace mmpt{

class CompositeTransferWithOwnershipResult
: public TransferWithOwnershipResult
{
public:
    std::vector<TransferOrder>  transfer_orders;
};

}

#endif
