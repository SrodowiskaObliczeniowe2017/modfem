

#include "TransfererWithOwnership.hpp"

namespace mmpt
{

TransfererWithOwnership::TransfererWithOwnership(int my_id, const mmpt_mesh & parallel_mesh)
: Transferer(my_id, parallel_mesh)
{

}

void TransfererWithOwnership::doMassTransferInit(const std::vector<TransferOrder>& orders)
{
    Transferer::doMassTransferInit(orders);
    // create correct ownership for all vertices in order pack
    contested_vrts_ownerships.clear();

    pmesh.findContestedVerticesAssigmentToProcs(orders,
                                                contested_vrts_ownerships);

    mf_log_info("Mass Transfer init");

}

void TransfererWithOwnership::doMassTransferClear(const std::vector<TransferOrder>& orders,
                                 std::vector<const TransferResult*> &results)
{
    Transferer::doMassTransferClear(orders,results);
    contested_vrts_ownerships.clear();
    mf_log_info("Mass Transfer clear");
}

const TransferResult& TransfererWithOwnership::doTransfer(const TransferOrder& order)
{
    assert(order.Source_id > 0 || order.Source_id==PCC_ANY_PROC);
    assert(order.Source_id <= pcr_nr_proc() || order.Source_id==PCC_ANY_PROC);
    assert(order.Destination_id > 0);
    assert(order.Destination_id <= pcr_nr_proc());
    assert(order.Destination_id != order.Source_id);

    transfers_history.push_back(order);
    TransferWithOwnershipResult& r = * new TransferWithOwnershipResult();
    transfers_results.push_back(&r);

    if(pcr_my_proc_id() == order.Source_id) {

        doSend(order,r);
        doSendOwnership(order,r);
    }
    else if(pcr_my_proc_id() == order.Destination_id) {

        doRecv(order,r);
        doRecvOwnership(order,r);
    }

    return r;
}

void TransfererWithOwnership::doSendOwnership(const TransferOrder& order,TransferWithOwnershipResult& r) const
{
    // if transfer order is COPY then
    // (by assumption)
    // we only copy from process with lower id
    // and lower-id process will keep ownership.
    //


    //  those vertices will be kept on both processes,
    // so the one with lower id will have ownership over them.
    // Send global id of vertices to other process to get it updated.

    // Ownership of elements.
    if(order.t_policy == TRANSFER_MOVE) {

        //if(order.Destination_id > order.Source_id) {

            std::vector<int> contested_vertices;
            const int bufid = pcr_send_buffer_open(MMPC_MSG_ID_OWNERSHIP,PCC_DEFAULT_BUFFER_SIZE);

            if(contested_vrts_ownerships.empty()) {
                // gather ownership for only this order
                pmesh.findContestedExternalVerticesOfElements(order.Elements,
                                                              contested_vertices);
                // keep ownership for this proc
                // send local id of vertices with kept ownership to other proc
                int n_vertd_ids = contested_vertices.size();
                pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid,1,& n_vertd_ids);
                pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid,n_vertd_ids, contested_vertices.data() );
            }
            else {
                std::vector<EntOwn> ownerships_to_transfer;
                std::vector<int>       ownerships2vtx_trans_id;
                for(int i=0; i < r.new_vertices_transfered_ids.size(); ++i) {
                    std::map<int, EntOwn>::iterator it
                            = contested_vrts_ownerships.find(r.new_vertices_transfered_ids[i]);
                    if(it!=contested_vrts_ownerships.end()) {
                        if(it->second.id_at_owner == EntOwnerships::UNDEFINED) {
                           it->second.id_at_owner =  i + 1; // only for initial distribution
                        }
                        ownerships_to_transfer.push_back(it->second);
                        ownerships2vtx_trans_id.push_back(r.new_vertices_transfered_ids[i]);
//                        mf_log_info("doSendOwnership(M): vtx %d sending as (%d,%d)",
//                                    glob_id2vtx_trans_id.back(),
//                                    glob_ids_to_transfer.back().owner,
//                                    glob_ids_to_transfer.back().id_at_owner);
                    }

                }

                int n_vertd_ids = - ownerships_to_transfer.size(); // note the '-' sign
                pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid,1,& n_vertd_ids);
                pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid, EntOwn::size_in_ints*ownerships_to_transfer.size(), reinterpret_cast<const int*>(ownerships_to_transfer.data()) );
                pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid, ownerships2vtx_trans_id.size(), ownerships2vtx_trans_id.data() );
            }

            pcr_buffer_send(MMPC_MSG_ID_OWNERSHIP,bufid,order.Destination_id);

        //}
        }
    else if(order.t_policy == TRANSFER_COPY) {
        std::vector<EntOwn> glob_ids_to_transfer;
        std::vector<int>       glob_id2vtx_trans_id;
        EntOwn tmp_gid;

        for(int v=0; v < r.new_vertices_transfered_ids.size(); ++v) {
            pmesh.ownership.GetOwnership<MMC_NODE>(r.new_vertices_transfered_ids[v],tmp_gid);

            glob_ids_to_transfer.push_back(tmp_gid);
            glob_id2vtx_trans_id.push_back(r.new_vertices_transfered_ids[v]);

//            mf_log_info("doSendOwnership(C): vtx %d sending as (%d,%d)",
//                        glob_id2vtx_trans_id.back(),
//                        glob_ids_to_transfer.back().owner,
//                        glob_ids_to_transfer.back().id_at_owner);

    }

        const int bufid = pcr_send_buffer_open(MMPC_MSG_ID_OWNERSHIP,PCC_DEFAULT_BUFFER_SIZE);
        int n_vertd_ids = glob_ids_to_transfer.size();
        pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid,1,& n_vertd_ids);
        pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid, EntOwn::size_in_ints*glob_ids_to_transfer.size(), reinterpret_cast<const int*>(glob_ids_to_transfer.data()));
        pcr_buffer_pack_int(MMPC_MSG_ID_OWNERSHIP,bufid, glob_id2vtx_trans_id.size(), glob_id2vtx_trans_id.data() );
        pcr_buffer_send(MMPC_MSG_ID_OWNERSHIP,bufid,order.Destination_id);
    } //! else if
}

void TransfererWithOwnership::doRecvOwnership(const TransferOrder& order,TransferWithOwnershipResult& r)
{
    // Ownership of vertices.
    if(order.t_policy == TRANSFER_MOVE) {
         //if(order.Destination_id > order.Source_id) {
            // keep ownership for this proc
            // send local id of vertices with kept ownership to other proc
            const int bufid = pcr_buffer_receive(MMPC_MSG_ID_OWNERSHIP, order.Source_id,PCC_DEFAULT_BUFFER_SIZE);
            int n_vertd_ids = 0;
            pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid,1,& n_vertd_ids);

            if(n_vertd_ids > 0) {
                r.ownerships2vtx_trans_id.resize(n_vertd_ids);
                pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid,n_vertd_ids, r.ownerships2vtx_trans_id.data() );

                EntOwn foreign_vtx_id(order.Source_id,EntOwnerships::UNDEFINED);

                for(int v=0; v < r.ownerships2vtx_trans_id.size(); ++v) {
                    //mf_log_info("Contested MOVE vertex: %d", r.glob_id2vtx_trans_id[v]);
                    foreign_vtx_id.id_at_owner = r.ownerships2vtx_trans_id[v];
                    mf_check_debug(foreign_vtx_id.id_at_owner>0,"Wrong id(%d)",foreign_vtx_id.id_at_owner);

                    r.transfer_vertices_ownerships.push_back(foreign_vtx_id);
                }
            }
            else { // n_vert_ids < 0
                n_vertd_ids*=-1; // back to positive values

                assert(r.transfer_vertices_ownerships.empty());
                r.transfer_vertices_ownerships.resize(n_vertd_ids);
                pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid, EntOwn::size_in_ints*n_vertd_ids, reinterpret_cast<int*>(r.transfer_vertices_ownerships.data()));
                r.ownerships2vtx_trans_id.resize(n_vertd_ids);
                pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid, n_vertd_ids,  r.ownerships2vtx_trans_id.data() );
            }
            pcr_recv_buffer_close(MMPC_MSG_ID_OWNERSHIP,bufid);
        //}
        // else {
        //    mf_fatal_err("Not implemented!");
        //}

    }
    else { /// TRANSFER_COPY
        // All vertices are supposed to be cloned.
        // Ownership stays at source proc.
         const int bufid = pcr_buffer_receive(MMPC_MSG_ID_OWNERSHIP, order.Source_id,PCC_DEFAULT_BUFFER_SIZE);
         int n_vertd_ids = 0;
         pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid,1,& n_vertd_ids);
         assert(r.transfer_vertices_ownerships.empty());
         r.transfer_vertices_ownerships.resize(n_vertd_ids);
         pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid, EntOwn::size_in_ints*n_vertd_ids, reinterpret_cast<int*>(r.transfer_vertices_ownerships.data()) );
         r.ownerships2vtx_trans_id.resize(n_vertd_ids);
         pcr_buffer_unpack_int(MMPC_MSG_ID_OWNERSHIP,bufid, n_vertd_ids,  r.ownerships2vtx_trans_id.data() );
         pcr_recv_buffer_close(MMPC_MSG_ID_OWNERSHIP,bufid);

    }

    assert(r.ownerships2vtx_trans_id.size() == r.transfer_vertices_ownerships.size());
}

}//!namespace

