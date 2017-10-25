#ifndef TRANSFERER_HPP_
#define TRANSFERER_HPP_

#include <vector>
#include <map>

#include "mmh_intf.h"
#include "mmph_mesh.h"
#include "TransferData.hpp"

namespace mmpt
{

enum MSG_IDS {
    MMPC_MSG_ID_COPYING,
    MMPC_TRANSFER_MSG_ID,
    MMPC_DISTIRB_MSG_ID,
    MMPC_ADAPT_REF_MSG_ID,
    MMPC_ADAPT_PROJ_READY_MSG_ID,
    MMPC_ADAPT_HS_MSG_ID,
    MMPC_MSG_ID_OWNERSHIP,
    MMPC_MSG_ID_COMPOSITE_ORDER,
    MMPC_MSG_ID_MESH_BASE
};



///
/// \brief The Transferer class is reponsible for realization, providing a results and storing
/// of a transfers between subdomains according to the given Transfer orders.
///
class Transferer
{
public:
    Transferer(int my_id, const mmpt_mesh & parallel_mesh);
    virtual ~Transferer();

    virtual const TransferResult& doTransfer(const TransferOrder& order);
    void doTransfer(const std::vector<TransferOrder>& orders, std::vector<const TransferResult*> &results);


protected:

    void doSend(const TransferOrder& order,TransferResult& r) const ;
    void doRecv(const TransferOrder& order,TransferResult& r) const ;

    virtual void doMassTransferInit(const std::vector<TransferOrder>& orders);
    virtual void doMassTransferClear(const std::vector<TransferOrder>& orders,
                                     std::vector<const TransferResult*> &results);

    const mmpt_mesh & pmesh;
    int my_subdomain_id;
    int last_applaied_transfer;
    std::vector<TransferOrder>  transfers_history;
    std::vector<TransferResult*> transfers_results;
    const  std::vector<TransferOrder>* current_orders;
};

}//!namespace

#endif //TRANSFERER_HPP_
