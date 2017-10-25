#ifndef COMPOSITE_TRANSFERER_WITH_OWNERSHIP_HPP_
#define COMPOSITE_TRANSFERER_WITH_OWNERSHIP_HPP_

#include "TransfererWithOwnership.hpp"
#include "CompositeTransferData.hpp"

namespace mmpt
{

class CompositeTransfererWithOwnership
: public TransfererWithOwnership
{
public:
    CompositeTransfererWithOwnership(int my_id, const mmpt_mesh & parallel_mesh);
    using Transferer::doTransfer;
    const TransferResult& doTransfer(const TransferOrder& order);

protected:

    void doMassTransferInit(const std::vector<TransferOrder>& orders);
    void doMassTransferClear(const std::vector<TransferOrder>& orders,
                                     std::vector<const TransferResult*> &results);

    void doSendRelativeOrders(const TransferOrder& order,
                              CompositeTransferWithOwnershipResult & r);
    void doRecvRelativeOrders(const TransferOrder & order,
                              CompositeTransferWithOwnershipResult & r);

};

}

#endif // COMPOSITE_TRANSFERER_WITH_OWNERSHIP_HPP_
