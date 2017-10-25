#define BOOST_TEST_MODULE ModFEM_mmpl_adapter_test
#include <boost/test/included/unit_test.hpp>

#include <boost/filesystem.hpp>

#include "mmh_intf.h"
#include "uth_io_files.h"

#include "mmph_intf.h"
#include "pch_intf.h"
#include "mmph_adapter.h"
#include "CompositeTransfererWithOwnership.hpp"


BOOST_AUTO_TEST_SUITE( ModFEM_mmpl_adapter_test )

std::string ctest_dir = "ctest/modfem_test/";
char  io_name[1024] = "modfem_test_mpi_out.txt";
FILE* f_out=NULL;
int n_proc=0, my_proc_id=0;

int load_2_prism()
{
    int mesh_id=-1;
    if(pcr_is_this_master()) {
        /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
        #ifdef MODFEM_NEW_MPI
          // Using MODFEM_NEW_MPI (mmpd_adapter)
          mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
          if(mesh_id > 0) {
           int n_read = utr_io_read_mesh(mesh_id,
                                         (ctest_dir+"../../meshes/cube/no_BC/dat/").c_str(),
                                         "2_prism_mesh.dat",
                                         MMC_MOD_FEM_PRISM_DATA);
           BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
          }
        #else
          // Using default mpi (mmpd_prism)
          // Create full mesh file path.
          sprintf(arg, "%s/%s", Work_dir, Mesh_file);
          // Initialize sequential/shared memory  mesh.
          mesh_id = mmr_init_mesh(mesh_desc, arg, stdout);
        #endif
        // !MODFEM_NEW_MPI

        BOOST_CHECK(mmr_get_nr_elem(mesh_id) == 2);
    }
    else if(pcv_my_proc_id == 2) {
        mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
    }
    return mesh_id;
}

int load_6_tetra()
{
    std::string dir = ctest_dir+"../../meshes/cube/no_BC/nas/";


    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     dir.c_str(),
                                     "mesh.nas",
                                     MMC_NASTRAN_DATA);
       BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
      }
    #else
      // Using default mpi (mmpd_prism)
      // Create full mesh file path.
      sprintf(arg, "%s/%s", Work_dir, Mesh_file);
      // Initialize sequential/shared memory  mesh.
      mesh_id = mmr_init_mesh(mesh_desc, arg, stdout);
    #endif
    /// !MODFEM_NEW_MPI


    BOOST_CHECK(mesh_id > 0);
    return mesh_id;
}

void print_mesh(const int Mesh_id)
{
    mf_log_info(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
    mf_log_info("Mesh %d at this point:", Mesh_id);
    int owner, id_at_owner;
    int nel=0;
    while(0 != (nel=mmr_get_next_elem_all(Mesh_id,nel))) {
        owner = mmpr_el_owner(Mesh_id, nel);
        id_at_owner = mmpr_el_id_at_owner(Mesh_id,nel);
        mf_log_info("Element %d (%d,%d)",
                    nel,owner, id_at_owner);
    }

    int nno=0;
    double coords[3]={0.0};
    while(0 != (nno=mmr_get_next_node_all(Mesh_id,nno))) {
        owner = mmpr_ve_owner(Mesh_id, nno);
        id_at_owner = mmpr_ve_id_at_owner(Mesh_id,nno);
        mmr_node_coor(Mesh_id,nno,coords);
        mf_log_info("Node %d (%d,%d) [%lf,%lf,%lf]",
                    nno,owner, id_at_owner,coords[0],coords[1],coords[2]);
    }
    mf_log_info("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
}

BOOST_AUTO_TEST_CASE( setup )
{
    utv_log_out = f_out = stdout;
    int n_loops=0;
    while( !boost::filesystem::exists(ctest_dir+"test_beacon.txt") && (n_loops<5) ) {
        ctest_dir = "../"+ctest_dir;
        ++n_loops;
    }

    BOOST_CHECK( boost::filesystem::exists(ctest_dir+"test_beacon.txt") );

    BOOST_CHECK( 0 <= pcr_init_parallel(& boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv,
                      "",io_name,
                      &f_out, & n_proc, &my_proc_id) );

    BOOST_REQUIRE_MESSAGE(n_proc >= 2, "ModFEM_mmpl_adapter_test have to be lunched as [at least] 2 process mpi (e.g. mpirun -np 2 MOD_FEM_test_mmpl_adapter");

}

BOOST_AUTO_TEST_CASE( simple_transfer_to_empty_2nd_proc_move )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        int el_id = 2;

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        Transferer t(my_proc_id,pmesh);

        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const TransferResult &tr = t.doTransfer(to);

        if(pcr_is_this_master()) {
            BOOST_CHECK( tr.new_elems.size() == 0 );
            BOOST_CHECK( tr.new_vertices_coords.size() == 0 );
        }
        else if(pcv_my_proc_id == 2) {
            BOOST_CHECK( tr.new_elems.size() == 1 );
            BOOST_CHECK( tr.new_vertices_coords.size() == 18);
        }
    }
}


BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_move )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        int el_id = 2;
        const int elgroup = 12;

        if(my_proc_id == 1) {
            mmr_el_set_groupID(mesh_id,el_id,elgroup);
        }

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        Transferer t(my_proc_id,pmesh);

        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const TransferResult &tr = t.doTransfer(to);

        pmesh.apply(to,tr);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);
            BOOST_CHECK( mmr_el_groupID(mesh_id,1) == elgroup);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);
                ++i;
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_move_2 )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        int el_id = 1;

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        Transferer t(my_proc_id,pmesh);

        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const TransferResult &tr = t.doTransfer(to);

        pmesh.apply(to,tr);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);
                ++i;
            }
        }
    }
}

/// END OF UNIT TESTS FOR TRANSFERER-move functionality


BOOST_AUTO_TEST_CASE( simple_transfer_to_empty_2nd_proc_copy)
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        int el_id = 2;

        BOOST_REQUIRE(mmpr_select_mesh(mesh_id) != NULL);

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        Transferer t(my_proc_id,pmesh);

        TransferOrder to(1,2,TRANSFER_COPY, &el_id, 1);

        const TransferResult &tr = t.doTransfer(to);

        if(pcr_is_this_master()) {
            BOOST_CHECK( tr.new_elems.size() == 0 );
            BOOST_CHECK( tr.new_vertices_coords.size() == 0 );
        }
        else if(pcv_my_proc_id == 2) {
            BOOST_CHECK( tr.new_elems.size() == 1 );
            BOOST_CHECK( tr.new_vertices_coords.size() == 18);
        }
    }

}
//
//

BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_copy)
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        int el_id = 2;

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        Transferer t(my_proc_id,pmesh);

        TransferOrder to(1,2,TRANSFER_COPY, &el_id, 1);

        const TransferResult &tr = t.doTransfer(to);

        pmesh.apply(to,tr);



        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 2);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 8);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 9);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);
                ++i;
            }
        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);
                ++i;
            }
        }
    }

}
//
//
BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_copy_w_ownership)
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        int el_id = 2;

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        TransfererWithOwnership t(my_proc_id,pmesh);

        TransferOrder to(1,2,TRANSFER_COPY, &el_id, 1);

        const TransferWithOwnershipResult &tr = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to));

        pmesh.apply(to,tr);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);



        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 2);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 8);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 9);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);
                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);


        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);
                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),1);

        }

        int ve=0;
        while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),1);
        }

    }

}


BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_move_2_w_ownership_v1 )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        TransfererWithOwnership t(my_proc_id,pmesh);

        int el_id = 2;
        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const TransferWithOwnershipResult &tr = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to));

        pmesh.apply(to,tr);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),my_proc_id);// 4
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),my_proc_id);//8

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),7);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),6); // loc

        }
    }
}
//
//
//

BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_move_2_w_ownership_v2 )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();



        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        TransfererWithOwnership t(my_proc_id,pmesh);

        int el_id = 1;
        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const TransferWithOwnershipResult &tr = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to));

        pmesh.apply(to,tr);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),my_proc_id); // 2
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),my_proc_id); // 6
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),1);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),2); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),5); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),7);

        }
    }
}


BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_to_empty_2nd_proc_move_2_w_ownership_overlap_v1 )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
       const int mesh_id = load_2_prism();

        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        TransfererWithOwnership t(my_proc_id,pmesh);

        int el_id = 2;
        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const TransferWithOwnershipResult &tr = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to));

        pmesh.apply(to,tr);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 1);
            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 6);
            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK(i == iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),my_proc_id);// 4
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),my_proc_id);//8

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),7);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),6); // loc

        }

        // sending overlap
        int ovl_el_id = 1;
        TransferOrder to_ovl(1,2,TRANSFER_COPY, &ovl_el_id, 1);

        const TransferWithOwnershipResult &tr_ovl = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to_ovl));

        pmesh.apply(to_ovl,tr_ovl);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap in this proc!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            // overlap in this proc!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_ELEMENT>(2),ovl_el_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),my_proc_id);// 4
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),my_proc_id);//8
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(7),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(8),1);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1); // owner proc 1
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),3); // owner proc 1
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5); // owner proc 1
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),7); // owner proc 1
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),6); // loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(7),2); // owner proc 1
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(8),6);  // owner proc 1
        }
    }
}
//
//
//

BOOST_AUTO_TEST_CASE( no_masstransfer )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();


        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        CompositeTransfererWithOwnership t(my_proc_id,pmesh);

        int el_id = 1;
        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);

        const CompositeTransferWithOwnershipResult &tr = static_cast<const CompositeTransferWithOwnershipResult &>(t.doTransfer(to));

        pmesh.apply(to,tr);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),my_proc_id); // 2
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),my_proc_id); // 6
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),1);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),2); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),5); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),7);

        }

        int ovl_el_id = 2;
        TransferOrder to_ovl(1,2,TRANSFER_COPY, &ovl_el_id, 1);

        const TransferWithOwnershipResult &tr_ovl = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to_ovl));

        pmesh.apply(to_ovl,tr_ovl);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap in this proc
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            // overlap in this proc
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_ELEMENT>(2),ovl_el_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),my_proc_id); // 2
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),my_proc_id); // 6
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(7),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(8),1);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),2); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),5); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),7);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(7),4);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(8),8);
        }
    }
}

BOOST_AUTO_TEST_CASE( masstransfer )
{
    if(my_proc_id < 3 && (n_proc>=2)) {
        const int mesh_id = load_2_prism();



        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

        using namespace mmpt;

        CompositeTransfererWithOwnership t(my_proc_id,pmesh);

        int el_id = 1;

        std::vector<TransferOrder> t_orders;
        TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);
        t_orders.push_back(to);

        std::vector<const TransferResult*>    t_results;

        t.doTransfer(t_orders,t_results);

        pmesh.apply(t_orders[0],
                static_cast<const TransferWithOwnershipResult&>(*t_results[0]) );

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        if(pcr_is_this_master()) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }
            pmesh.neighbours.insert(2);
        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),my_proc_id); // 2
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),my_proc_id); // 6
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),1);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),2); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),5); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),7);

        }

        int ovl_el_id = 2;
        TransferOrder to_ovl(1,2,TRANSFER_COPY, &ovl_el_id, 1);

        const TransferWithOwnershipResult &tr_ovl = static_cast<const TransferWithOwnershipResult &>(t.doTransfer(to_ovl));

        pmesh.apply(to_ovl,tr_ovl);

        BOOST_CHECK_EQUAL(pmesh.ownership.getLocalOwnerID(),my_proc_id);

        BOOST_CHECK_EQUAL(pmesh.neighbours.size(),1);

        if(pcr_is_this_master()) {
            BOOST_CHECK_EQUAL(* pmesh.neighbours.begin(), 2);

            // no overlap in this proc
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),my_proc_id);

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(pmesh.mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(ve),my_proc_id);
            }

        }
        else if(pcv_my_proc_id == 2) {
            BOOST_CHECK_EQUAL(* pmesh.neighbours.begin(), 1);
            // overlap in this proc
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(1),my_proc_id);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_ELEMENT>(2),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_ELEMENT>(2),ovl_el_id);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(2),my_proc_id); // 2
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(3),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(4),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(5),my_proc_id); // 6
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(6),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(7),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetOwner<MMC_NODE>(8),1);

            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(1),1);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(2),2); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(3),3);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(4),5);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(5),5); //loc
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(6),7);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(7),4);
            BOOST_CHECK_EQUAL(pmesh.ownership.GetId_at_owner<MMC_NODE>(8),8);
        }
    }
}

//BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_6_tetra_copy)
//{
//    if(my_proc_id < 3 && (n_proc>=2)) {
//        int mesh_id=-1;
//        if(pcr_is_this_master()) {
//            /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
//            #ifdef MODFEM_NEW_MPI
//              // Using MODFEM_NEW_MPI (mmpd_adapter)
//              mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//              if(mesh_id > 0) {
//                   int n_read = utr_io_read_mesh(mesh_id,
//                                             (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
//                                             "mesh.nas",
//                                             MMC_NASTRAN_DATA);
//               BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
//              }
//            #else
//              // Using default mpi (mmpd_prism)
//              // Create full mesh file path.
//              sprintf(arg, "%s/%s", Work_dir, Mesh_file);
//              // Initialize sequential/shared memory  mesh.
//              mesh_id = mmr_init_mesh(mesh_desc, arg, stdout);
//            #endif
//            // !MODFEM_NEW_MPI
//
//            BOOST_CHECK(mmr_get_nr_elem(mesh_id) == 6);
//        }
//        else if(pcv_my_proc_id == 2) {
//            mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//        }
//
//        for(int i=1; i <=6; ++i) {
//
//            int el_id = i;
//
//            mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);
//
//            using namespace mmpt;
//
//            Transferer t(my_proc_id,n_proc+1,pmesh);
//
//            TransferOrder to(1,2,TRANSFER_COPY, &el_id, 1);
//
//            const TransferResult &tr = t.doTransfer(to);
//
//            pmesh.apply(to,tr);
//
//
//
//            if(pcr_is_this_master()) {
//                // no overlap!
//                BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 6);
//                BOOST_CHECK( mmr_get_nr_node(mesh_id) == 8);
////                BOOST_CHECK( mmr_get_nr_face(mesh_id) == 9);
//
//                int iNode=0, j=1;
//                while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
//                    BOOST_CHECK(j == iNode);
//                    ++j;
//                }
//            }
//            else if(pcv_my_proc_id == 2) {
//                // no overlap!
//                BOOST_CHECK( mmr_get_nr_elem(mesh_id) == i);
//                BOOST_CHECK_MESSAGE( mmr_get_nr_node(mesh_id) >= 2+i,
//                                     "Incorrect number of vertices for " << i << ": " << mmr_get_nr_node(mesh_id));
////                BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);
//
//                int iNode=0, j=1;
//                while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
//                    BOOST_CHECK(j == iNode);
//                    ++j;
//                }
//
//                switch(i) {
//                case 1: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),4);
//                    break;
//                case 2:
//                case 3: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),6);
//                    break;
//                case 4:
//                case 5:
//                case 6: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),8);
//                    break;
//
//                }
//            }
//         }
//
//        if(pcr_is_this_master()) {
//            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 6);
//            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
////            BOOST_CHECK_EQUAL( mmr_get_nr_edge(mesh_id) , 0);
////            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 0);
//
//        }
//        else if(pcv_my_proc_id == 2) {
//            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 6);
//            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
////            BOOST_CHECK( mmr_get_nr_edge(mesh_id) == 0);
////            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 0);
//
//        }
//    }
//
//}
//
//BOOST_AUTO_TEST_CASE( simple_transfer_and_apply_6_tetra_move)
//{
//    if(my_proc_id < 3 && (n_proc>=2)) {
//        int mesh_id=-1;
//        if(pcr_is_this_master()) {
//            /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
//            #ifdef MODFEM_NEW_MPI
//              // Using MODFEM_NEW_MPI (mmpd_adapter)
//              mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//              if(mesh_id > 0) {
//                   int n_read = utr_io_read_mesh(mesh_id,
//                                             (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
//                                             "mesh.nas",
//                                             MMC_NASTRAN_DATA);
//               BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
//              }
//            #else
//              // Using default mpi (mmpd_prism)
//              // Create full mesh file path.
//              sprintf(arg, "%s/%s", Work_dir, Mesh_file);
//              // Initialize sequential/shared memory  mesh.
//              mesh_id = mmr_init_mesh(mesh_desc, arg, stdout);
//            #endif
//            // !MODFEM_NEW_MPI
//
//            BOOST_CHECK(mmr_get_nr_elem(mesh_id) == 6);
//        }
//        else if(pcv_my_proc_id == 2) {
//            mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//        }
//
//        mmpt_mesh& pmesh = *mmpr_select_mesh(mesh_id);

//        using namespace mmpt;
//
//        const int n_faces = mmr_get_nr_face(mesh_id);
//
//        // Transferring all elements from 1 to 2 in a loop.
//        for(int i=1; i <=6; ++i) {
//
//            int el_id = i;
//
//            Transferer t(my_proc_id,n_proc+1,pmesh);
//
//            TransferOrder to(1,2,TRANSFER_MOVE, &el_id, 1);
//
//            const TransferResult &tr = t.doTransfer(to);
//
//            BOOST_CHECK(!tr.new_elems.empty() || !to.Elements.empty());
//
//            pmesh.apply(to,tr);
//
//
//            if(pcr_is_this_master()) {
//                // no overlap!
//                BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 6-i);
//                BOOST_CHECK( mmr_get_nr_face(mesh_id) <= n_faces-i);
//
////                switch(i) {
////                case 1:
////                case 2:BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),8);
////                    break;
////                case 3: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),7);
////                    break;
////                case 4:BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),6);
////                    break;
////                case 5:BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),5);
////                    break;
////                case 6: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),4);
////                    break;
////                }
//
//            }
//            else if(pcv_my_proc_id == 2) {
//                // no overlap!
//                BOOST_CHECK( mmr_get_nr_elem(mesh_id) == i);
//                BOOST_CHECK_MESSAGE( mmr_get_nr_node(mesh_id) >= 2+i,
//                                     "Incorrect number of vertices for " << i << ": " << mmr_get_nr_node(mesh_id));
////                BOOST_CHECK( mmr_get_nr_face(mesh_id) == 5);
//
//                int iNode=0, j=1;
//                while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
//                    BOOST_CHECK(j == iNode);
//                    ++j;
//                }
//
//                switch(i) {
//                case 1: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),4);
//                    break;
//                case 2:
//                case 3: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),6);
//                    break;
//                case 4:
//                case 5:
//                case 6: BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),8);
//                    break;
//
//                }
//            }
//         }
//
//        if(pcr_is_this_master()) {
//            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 0);
//            BOOST_CHECK( mmr_get_nr_node(mesh_id) == 0);
//            BOOST_CHECK( mmr_get_nr_edge(mesh_id) == 0);
//            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 0);
//
//        }
//        else if(pcv_my_proc_id == 2) {
//            BOOST_CHECK( mmr_get_nr_elem(mesh_id) == 6);
//            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
////            BOOST_CHECK( mmr_get_nr_edge(mesh_id) == 0);
////            BOOST_CHECK( mmr_get_nr_face(mesh_id) == 0);
//
//        }
//    }
//
//}



BOOST_AUTO_TEST_CASE( simple_domain_decomposition_2_dom_2_prism )
{
    if(pcv_nr_proc == 2) {
   const int mesh_id = load_2_prism();

    BOOST_CHECK( 0 < mmpr_init_mesh(MMC_INIT_GEN_LEVEL,mesh_id,pcr_nr_proc(),pcr_my_proc_id()) );

    BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 1);

    print_mesh(mesh_id);

    if(pcv_my_proc_id == 1) {
            int neig_proc=0;
            mmpr_neighbouring_procs(mesh_id, &neig_proc);
            BOOST_CHECK_EQUAL( neig_proc, 2 );

            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

            BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,1), 1);
            BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,2), 2);

            BOOST_CHECK_EQUAL(mmpr_el_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_el_id_at_owner(mesh_id,2),1);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,7),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,8),2);

            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3); //loc
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),6);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),5);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),6); //loc
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,7),7);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,8),3);
    }
    else if(pcv_my_proc_id == 2) {
        int neig_proc=0;
        mmpr_neighbouring_procs(mesh_id, &neig_proc);
        BOOST_CHECK_EQUAL( neig_proc, 1 );

        // no overlap!
        BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
        BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
        BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

        BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,1), 2);
        BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,2), 1);

        BOOST_CHECK_EQUAL(mmpr_el_id_at_owner(mesh_id,1),1);

        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),2);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),1);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),2);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,7),1);
        BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,8),1);

        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),3);
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3); //loc
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),5);
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),7);
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),6); //loc
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,7),2);
        BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,8),6);
    }
    }
}

BOOST_AUTO_TEST_CASE( simple_dd_t4_2_procs )
{
    if(pcv_nr_proc == 2) {
        const int mesh_id = load_6_tetra();

        BOOST_CHECK( 0 < mmpr_init_mesh(MMC_INIT_GEN_LEVEL,mesh_id,pcr_nr_proc(),pcr_my_proc_id()) );

        BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 1);

        if( pcr_is_this_master() ) {
            int neig_proc=0;
            mmpr_neighbouring_procs(mesh_id, &neig_proc);
            BOOST_CHECK_EQUAL( neig_proc, 2 );

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),5);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),7);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),4);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,8),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,8),8);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,7),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,7),5);
        }
        else if(pcv_my_proc_id == 2) {
            int neig_proc=0;
            mmpr_neighbouring_procs(mesh_id, &neig_proc);
            BOOST_CHECK_EQUAL( neig_proc, 1 );

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),5);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),7);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),PCC_MASTER_PROC_ID);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),PCC_MASTER_PROC_ID);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),PCC_MASTER_PROC_ID);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),pcv_my_proc_id);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),4);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),pcv_my_proc_id);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),PCC_MASTER_PROC_ID);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),8);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,7),PCC_MASTER_PROC_ID);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,7),1);
        }
    }
}

BOOST_AUTO_TEST_CASE( simple_dd_t4_6_procs )
{
    if(pcv_nr_proc == 6) {
       const int mesh_id = load_6_tetra();


        BOOST_CHECK(mesh_id > 0);

        BOOST_CHECK( 0 < mmpr_init_mesh(MMC_INIT_GEN_LEVEL,mesh_id,pcr_nr_proc(),pcr_my_proc_id()) );

        int correct_neigs[6]={0};
        switch(pcv_my_proc_id) {
            case 1:
            correct_neigs[0]=2;
            correct_neigs[1]=3;
            correct_neigs[2]=4;
            correct_neigs[3]=5;
            BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 4);

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),3);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),6);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),5);
            break;
            case 2:
            correct_neigs[0]=1;
            correct_neigs[1]=3;
            BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 2);

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),2);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),5);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),3);    //
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),2);          //
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),4);    //
            break;
            case 3:
            correct_neigs[0]=1;
            correct_neigs[1]=2;
            correct_neigs[2]=4;
            correct_neigs[3]=6;
            BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 4);

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),4);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),7);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),4);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),3);

            break;
            case 4:
            correct_neigs[0]=1;
            correct_neigs[1]=5;
            BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 2);

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),2);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),5);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1), 1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1), 2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2), 1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2), 5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3), 4);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3), 3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4), 4);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4), 4);
            break;
            case 5:
            correct_neigs[0]=1;
            correct_neigs[1]=2;
            correct_neigs[2]=4;
            correct_neigs[3]=6;
            BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 4);

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),4);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),7);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),4);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),4);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),4);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),3);

            break;
            case 6:
            correct_neigs[0]=1;
            correct_neigs[1]=2;
            correct_neigs[2]=3;
            correct_neigs[3]=4;
            correct_neigs[4]=5;
            BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 5);

            BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),3);
            BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),6);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),3);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),5);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),4);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),4);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),4);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),2);
            break;
        }

        int neig_proc[6]={0};
        mmpr_neighbouring_procs(mesh_id, neig_proc);
        BOOST_CHECK_EQUAL_COLLECTIONS(neig_proc,
                                      neig_proc+mmpr_n_neighbouring_procs(mesh_id),
                                      correct_neigs,
                                      correct_neigs+mmpr_n_neighbouring_procs(mesh_id));
    }
}

BOOST_AUTO_TEST_CASE( adapt_2_prism )
{
    if(pcv_nr_proc == 2) {
        const int mesh_id = load_2_prism();

        BOOST_CHECK( 0 < mmpr_init_mesh(MMC_INIT_GEN_LEVEL,mesh_id,pcr_nr_proc(),pcr_my_proc_id()) );

        BOOST_CHECK_EQUAL( mmpr_n_neighbouring_procs(mesh_id), 1);

        if(pcv_my_proc_id == 1) {
                int neig_proc=0;
                mmpr_neighbouring_procs(mesh_id, &neig_proc);
                BOOST_CHECK_EQUAL( neig_proc, 2 );

                // no overlap!
                BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
                BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
                BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

                BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,1), 1);
                BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,2), 2);

                BOOST_CHECK_EQUAL(mmpr_el_id_at_owner(mesh_id,1),1);
                BOOST_CHECK_EQUAL(mmpr_el_id_at_owner(mesh_id,2),1);

                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),2);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),1);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,7),1);
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,8),2);

                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),2);
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3); //loc
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),6);
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),5);
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),6); //loc
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,7),7);
                BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,8),3);
        }
        else if(pcv_my_proc_id == 2) {
            int neig_proc=0;
            mmpr_neighbouring_procs(mesh_id, &neig_proc);
            BOOST_CHECK_EQUAL( neig_proc, 1 );

            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

            BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,1), 2);
            BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,2), 1);

            BOOST_CHECK_EQUAL(mmpr_el_id_at_owner(mesh_id,1),1);

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),2);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,7),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,8),1);

            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),3);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3); //loc
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),5);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),7);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),6); //loc
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,7),2);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,8),6);
        }

        int* ref_list = (int*)calloc(1,sizeof(int) );
        ref_list[0]=1;
        int n_ref = 1;

        print_mesh(mesh_id);

        mmr_init_ref(mesh_id);

        mmpr_init_ref(mesh_id);

        mmpr_check_mesh(mesh_id);

        mmpr_update_ref_list(mesh_id,&n_ref,&ref_list);


        for(int i=0; i < n_ref; ++i) {
            mmr_refine_el(mesh_id,ref_list[i]);
        }


        mmpr_is_ready_for_proj_dof_ref(mesh_id);

        mmr_is_ready_for_proj_dof_ref(mesh_id);
        mmr_final_ref(mesh_id);

        mmpr_final_ref(mesh_id);

        BOOST_CHECK_EQUAL(mmr_get_nr_elem(mesh_id),16);
        BOOST_CHECK_EQUAL(mmr_get_max_elem_id(mesh_id),18);

        BOOST_CHECK_EQUAL(mmr_get_nr_node(mesh_id),2*(6+12)-9);

        print_mesh(mesh_id);
    }
}

BOOST_AUTO_TEST_CASE( setup_end )
{
    pcr_exit_parallel();
}

BOOST_AUTO_TEST_SUITE_END()
