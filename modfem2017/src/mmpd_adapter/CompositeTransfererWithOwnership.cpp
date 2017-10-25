#include "CompositeTransfererWithOwnership.hpp"
#include "CompositeTransferData.hpp"

namespace mmpt
{

CompositeTransfererWithOwnership::CompositeTransfererWithOwnership(int my_id, const mmpt_mesh & parallel_mesh)
: TransfererWithOwnership(my_id, parallel_mesh)
{}

void CompositeTransfererWithOwnership::doMassTransferInit(const std::vector<TransferOrder>& orders)
{
    TransfererWithOwnership::doMassTransferInit(orders);
}

void CompositeTransfererWithOwnership::doMassTransferClear(const std::vector<TransferOrder>& orders,
                                 std::vector<const TransferResult*> &results)
{
    TransfererWithOwnership::doMassTransferClear(orders,results);
}

const TransferResult& CompositeTransfererWithOwnership::doTransfer(const TransferOrder& order)
{

    assert(order.Source_id > 0 || order.Source_id==PCC_ANY_PROC);
    assert(order.Source_id <= pcr_nr_proc() || order.Source_id==PCC_ANY_PROC);
    assert(order.Destination_id > 0);
    assert(order.Destination_id <= pcr_nr_proc());
    assert(order.Destination_id != order.Source_id);

    transfers_history.push_back(order);
    CompositeTransferWithOwnershipResult& r = * new CompositeTransferWithOwnershipResult();
    transfers_results.push_back(&r);

    if(pcr_my_proc_id() == order.Source_id) {

        doSend(order,r);
        doSendOwnership(order,r);
        doSendRelativeOrders(order,r);

    }
    else if(pcr_my_proc_id() == order.Destination_id) {

        doRecv(order,r);
        doRecvOwnership(order,r);
        doRecvRelativeOrders(order,r);
    }

    return r;
}

void CompositeTransfererWithOwnership::doSendRelativeOrders(const TransferOrder& order,
                          CompositeTransferWithOwnershipResult & r)
{
    const int bufid = pcr_send_buffer_open(MMPC_MSG_ID_COMPOSITE_ORDER,PCC_DEFAULT_BUFFER_SIZE);
    const int n_transfer_orders = order.relative_transfer_orders.size();
    pcr_buffer_pack_int(MMPC_MSG_ID_COMPOSITE_ORDER,bufid,1, &n_transfer_orders);
    if(n_transfer_orders > 0) {
        int n_total_ints=0;
        for(int i=0; i < order.relative_transfer_orders.size(); ++i) {
            n_total_ints+=order.relative_transfer_orders[i].n_ints_to_serialize();
        }
        pcr_buffer_pack_int(MMPC_MSG_ID_COMPOSITE_ORDER,bufid,1, &n_total_ints );

        std::vector<int>    internal_buf(n_total_ints);
        int * buf_ptr = internal_buf.data();
        for(int o=0; o < order.relative_transfer_orders.size(); ++o) {
            buf_ptr = order.relative_transfer_orders[o].serialize(buf_ptr);
        }

        pcr_buffer_pack_int(MMPC_MSG_ID_COMPOSITE_ORDER,bufid,
                        internal_buf.size(),
                        internal_buf.data() );
    }
    pcr_buffer_send(MMPC_MSG_ID_COMPOSITE_ORDER,bufid,order.Destination_id);
}


void CompositeTransfererWithOwnership::doRecvRelativeOrders(const TransferOrder & order,
                          CompositeTransferWithOwnershipResult & r)
{
    const int bufid = pcr_buffer_receive(MMPC_MSG_ID_COMPOSITE_ORDER, order.Source_id,PCC_DEFAULT_BUFFER_SIZE);

    int n_rel_orders=0;
    pcr_buffer_unpack_int(MMPC_MSG_ID_COMPOSITE_ORDER,bufid,1, &n_rel_orders );

    if(n_rel_orders > 0) {
        r.transfer_orders.resize(n_rel_orders);

        int n_total_ints = 0;
        pcr_buffer_unpack_int(MMPC_MSG_ID_COMPOSITE_ORDER,bufid,1, &n_total_ints );
        std::vector<int>    internal_buf(n_total_ints);

        int * buf_ptr =  internal_buf.data();
        pcr_buffer_unpack_int(MMPC_MSG_ID_COMPOSITE_ORDER,bufid, n_total_ints, buf_ptr );
        pcr_recv_buffer_close(MMPC_MSG_ID_COMPOSITE_ORDER,bufid);

        for(int o=0; o < n_rel_orders; ++o) {
            buf_ptr = r.transfer_orders[o].deserialize(buf_ptr);
        }
    }
}

} // namespace
