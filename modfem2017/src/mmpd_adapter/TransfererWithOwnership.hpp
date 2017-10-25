#ifndef TRANSFERER_WITH_OWNERSHIP_HPP_
#define TRANSFERER_WITH_OWNERSHIP_HPP_


#include "mmph_mesh.h"
#include "Transferer.hpp"
#include "TransferWithOwnershipData.hpp"


namespace mmpt
{

///
/// \brief The Transferer class is reponsible for realization, providing a results and storing
/// of a transfers between subdomains according to the given Transfer orders.
///
class TransfererWithOwnership : public mmpt::Transferer
{
public:
    TransfererWithOwnership(int my_id, const mmpt_mesh & parallel_mesh);
    using Transferer::doTransfer;
    const TransferResult& doTransfer(const TransferOrder& order);


protected:

    void doSendOwnership(const TransferOrder& order,TransferWithOwnershipResult& r) const ;
    void doRecvOwnership(const TransferOrder& order,TransferWithOwnershipResult& r)  ;

    void doMassTransferInit(const std::vector<TransferOrder>& orders);
    void doMassTransferClear(const std::vector<TransferOrder>& orders,
                                     std::vector<const TransferResult*> &results);


    mutable std::map<int,EntOwn> contested_vrts_ownerships;
};

}//!namespace

#endif //TRANSFERER_WITH_OWNERSHIP_HPP_
