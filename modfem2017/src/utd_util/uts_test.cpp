#define BOOST_TEST_MODULE ModFEM_test_util
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <string>

#include "mmh_intf.h"
#include "uth_log.h"
#include "uth_mat.h"
#include "uth_bc.h"
#include "uth_io_compression.h"
#include "uth_io_files.h"


#ifdef PARALLEL
#include "mmph_intf.h"
#include "pch_intf.h"
#endif // PARALLEL

std::string ctest_dir = "ctest/modfem_test/";


BOOST_AUTO_TEST_SUITE( Utilities_materials )

BOOST_AUTO_TEST_CASE( Utilities_test_setup )
{
    int n_loops=0;
    while( !boost::filesystem::exists(ctest_dir+"test_beacon.txt") && (n_loops<5) ) {
        ctest_dir = "../"+ctest_dir;
        ++n_loops;
    }

    BOOST_CHECK_MESSAGE( boost::filesystem::exists(ctest_dir+"test_beacon.txt"),
                         "Invalid path: " << ctest_dir+"test_beacon.txt");
}

BOOST_AUTO_TEST_CASE( materials_normal )
{
    std::string dir = ctest_dir+"materials";
    BOOST_CHECK( UTE_FAIL != utr_mat_read(dir.c_str(),"materials.dat",stdout) );
    BOOST_CHECK( 3 == utr_mat_get_n_materials() );
}
BOOST_AUTO_TEST_CASE( materials_ill_formed_files )
{
    std::string dir = ctest_dir+"materials";
    BOOST_CHECK( UTE_FAIL == utr_mat_read(dir.c_str(),"no_such_file.dat",stdout) );
    BOOST_CHECK( UTE_SUCCESS == utr_mat_read(dir.c_str(),"empty_file.dat",stdout) );
}
BOOST_AUTO_TEST_CASE( materials_additonal_files )
{
    std::string dir = ctest_dir+"materials";
    BOOST_CHECK( UTE_FAIL != utr_mat_read(dir.c_str(),"1_material.dat",stdout) );
    BOOST_CHECK( 4 == utr_mat_get_n_materials() );
}
BOOST_AUTO_TEST_CASE( materials_normal_zero_id_bug_test )
{
    utr_mat_clear_all();
    std::string dir2 = ctest_dir+"materials/zero_material";
    BOOST_CHECK( UTE_FAIL != utr_mat_read(dir2.c_str(),"materials.dat",stdout) );
    BOOST_CHECK_EQUAL( 1 , utr_mat_get_n_materials() );
    BOOST_CHECK_EQUAL( 1 , utr_mat_get_matID(1) );
    BOOST_CHECK_EQUAL( 1 , utr_bc_get_blockID(1) );
}
BOOST_AUTO_TEST_CASE( materials_no_groups )
{
    utr_mat_clear_all();
    std::string dir2 = ctest_dir+"materials/zero_material";
    BOOST_CHECK_EQUAL( UTE_FAIL , utr_mat_read(dir2.c_str(),"no_such_file.dat",stdout) );
    BOOST_CHECK_EQUAL( 1 , utr_mat_get_n_materials() );
    BOOST_CHECK_EQUAL( 1 , utr_mat_get_matID(1) );
    BOOST_CHECK_EQUAL( UTE_GROUP_ID_DEFAULT_BLOCK , utr_bc_get_blockID(1) );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( Utilities_io_compression )

BOOST_AUTO_TEST_CASE( backward_compatiliblity_no_compression )
{
    std::string dir = ctest_dir+"compression/";
    char decompressed_file[255]={0};
    BOOST_CHECK( 0 == utr_io_decompress_file(dir.c_str(),"LDC_z_005_20x20_nocomp.nas",decompressed_file) );
    BOOST_CHECK( boost::filesystem::exists(dir+decompressed_file) );
}

BOOST_AUTO_TEST_CASE( decompression )
{
    std::string dir = ctest_dir+"compression";
    char decompressed_file[255]={0};
    BOOST_CHECK( utr_io_decompress_file(dir.c_str(),"LDC_z_005_20x20.nas.zip",decompressed_file) );
    BOOST_CHECK( boost::filesystem::exists(dir+"/LDC_z_005_20x20.nas") );
    BOOST_CHECK( boost::filesystem::exists(dir+"/"+decompressed_file) );

}

BOOST_AUTO_TEST_CASE( compression_dir_with_filename )
{
    std::string dir = ctest_dir+"compression";

    boost::filesystem::copy(dir+"/field90.dmp.ref", dir+"/field90.dmp");

    BOOST_CHECK( boost::filesystem::exists(dir+"/field90.dmp") );
    BOOST_CHECK( utr_io_compress_file(dir.c_str(),"field90.dmp") );
    BOOST_CHECK( boost::filesystem::exists(dir+"/field90.dmp.zip") );

    boost::filesystem::remove(dir+"/field90.dmp.zip");
    boost::filesystem::remove(dir+"/field90.dmp");

    //boost::filesystem::rename(dir+"/field90.dmp", dir+"/field90.dmp.ref");

}

BOOST_AUTO_TEST_CASE( decompression_filename_only )
{
    std::string file = ctest_dir+"compression/LDC_z_005_20x20.nas";
    char decompressed_file[255]={0};
    BOOST_CHECK( utr_io_decompress_file( NULL, (file+".zip").c_str() ,decompressed_file) );
    BOOST_CHECK( boost::filesystem::exists(file) );
    BOOST_CHECK( boost::filesystem::exists(decompressed_file) );
}


BOOST_AUTO_TEST_CASE( compression_filename_only )
{
    std::string dir = ctest_dir+"compression";

    boost::filesystem::copy((dir+"/field90.dmp.ref"), (dir+"/field90.dmp") );

    BOOST_CHECK( boost::filesystem::exists(dir+"/field90.dmp") );
    BOOST_CHECK( utr_io_compress_file(NULL,(dir+"/field90.dmp").c_str()) );
    BOOST_CHECK( boost::filesystem::exists(dir+"/field90.dmp.zip") );

    boost::filesystem::remove(dir+"/field90.dmp.zip");
    boost::filesystem::remove(dir+"/field90.dmp");

    //boost::filesystem::rename(dir+"/field90.dmp", dir+"/field90.dmp.ref");

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( Utilities_io_mesh_files )

BOOST_AUTO_TEST_CASE( mesh_nas_simple )
{
    std::string dir = ctest_dir+"../../meshes/cube/no_BC/nas/";
    std::string filename = "mesh.nas";

    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     dir.c_str(),
                                     filename.c_str(),
                                     MMC_NASTRAN_DATA);
       BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
      }
    #else
      // Using default mpi (mmpd_prism)
      // Create full mesh file path.
      char arg[1000]={0};
      sprintf(arg, "%s/%s", dir.c_str(), filename.c_str());
      // Initialize sequential/shared memory  mesh.
      mesh_id = mmr_init_mesh('n', arg, stdout);
    #endif
    /// !MODFEM_NEW_MPI
}

BOOST_AUTO_TEST_CASE( mesh_nas )
{
    std::string dir = ctest_dir+"../../meshes/cube/LDC/nas/";
    std::string filename = "LDC_z_001_100x100.nas";

    char decompressed_file[255]={0};
    utr_io_decompress_file(dir.c_str(),filename.c_str(),decompressed_file);
    filename = decompressed_file;

    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     dir.c_str(),
                                     filename.c_str(),
                                     MMC_NASTRAN_DATA);
       BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
      }
    #else
      // Using default mpi (mmpd_prism)
      // Create full mesh file path.
      char arg[1000]={0};
      sprintf(arg, "%s/%s", dir.c_str(), filename.c_str());
      // Initialize sequential/shared memory  mesh.
      mesh_id = mmr_init_mesh('n', arg, stdout);
    #endif
    /// !MODFEM_NEW_MPI
}

//BOOST_AUTO_TEST_CASE( mesh_jk_simple )
//{
//    std::string dir = ctest_dir+"../../meshes/cube/LDC/jk/struct/scale1/";
//    std::string filename = "mesh_jk_11x11.dat";

//    int mesh_id=-1;
//    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
//    #ifdef MODFEM_NEW_MPI
//      // Using MODFEM_NEW_MPI (mmpd_adapter)
//      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//      if(mesh_id > 0) {
//       int n_read = utr_io_read_mesh(mesh_id,
//                                     dir.c_str(),
//                                     filename.c_str(),
//                                     MMC_NASTRAN_DATA);
//       BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
//      }
//    #else
//      // Using default mpi (mmpd_prism)
//      // Create full mesh file path.
//      sprintf(arg, "%s/%s", Work_dir, Mesh_file);
//      // Initialize sequential/shared memory  mesh.
//      mesh_id = mmr_init_mesh(mesh_desc, arg, stdout);
//    #endif
//    /// !MODFEM_NEW_MPI
//}

//BOOST_AUTO_TEST_CASE( mesh_jk )
//{
//    std::string dir = ctest_dir+"../../meshes/cube/LDC/nas/";
//    std::string filename = "LDC_z_001_100x100.nas";

//    utr_io_decompress_file(dir.c_str(),filename.c_str());

//    int mesh_id=-1;
//    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
//    #ifdef MODFEM_NEW_MPI
//      // Using MODFEM_NEW_MPI (mmpd_adapter)
//      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//      if(mesh_id > 0) {
//       int n_read = utr_io_read_mesh(mesh_id,
//                                     dir.c_str(),
//                                     filename.c_str(),
//                                     MMC_NASTRAN_DATA);
//       BOOST_CHECK_MESSAGE(n_read > 0,"Mesh read not ok returned " << n_read);
//      }
//    #else
//      // Using default mpi (mmpd_prism)
//      // Create full mesh file path.
//      sprintf(arg, "%s/%s", Work_dir, Mesh_file);
//      // Initialize sequential/shared memory  mesh.
//      mesh_id = mmr_init_mesh(mesh_desc, arg, stdout);
//    #endif
//    /// !MODFEM_NEW_MPI
//}

/// TEST FOR PARALLEL MESH READING



BOOST_AUTO_TEST_CASE( parllel_2prism_mesh )
{
    std::string dir = ctest_dir+"../../meshes/cube/no_BC/dat/";
    std::string filename = "2_prism_mesh.dat";
    char mtype = 'p';

#ifdef PARALLEL
    char  io_name[1024] = "modfem_test_mpi_out.txt";
    FILE* f_out=NULL;
    int n_proc=0, my_proc_id=0;
    pcr_init_parallel(& boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv,
                      "",io_name,
                      &f_out, & n_proc, &my_proc_id);
#endif // PARALLEL

    const int mesh_id = utr_io_initialize_mesh(utv_log_out, dir.c_str(),
                           mtype, filename.c_str());

#ifdef PARALLEL
    if(pcr_is_this_master() && (pcv_nr_proc > 1))  {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 1);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 6);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 5);

            BOOST_CHECK_EQUAL( mmpr_el_owner(mesh_id,1), pcr_my_proc_id());

            int ve=0;
            while(0 != (ve=mmr_get_next_node_all(mesh_id,ve)) ) {
                BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,ve),pcr_my_proc_id());
            }

        }
        else if(pcv_my_proc_id == 2) {
            // no overlap!
            BOOST_CHECK_EQUAL( mmr_get_nr_elem(mesh_id) , 2);
            BOOST_CHECK_EQUAL( mmr_get_nr_node(mesh_id) , 8);
            BOOST_CHECK_EQUAL( mmr_get_nr_face(mesh_id) , 9);

            int iNode=0, i=1;
            while( (iNode = mmr_get_next_node_all(mesh_id,iNode)) != 0 ) {
                BOOST_CHECK_EQUAL(i , iNode);

                ++i;
            }

            BOOST_CHECK_EQUAL(mmpr_el_owner(mesh_id,1),pcr_my_proc_id());

            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,2),pcr_my_proc_id()); // 2
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,3),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,4),1);
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,5),pcr_my_proc_id()); // 6
            BOOST_CHECK_EQUAL(mmpr_ve_owner(mesh_id,6),1);

            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,1),1);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,2),2); //loc
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,3),3);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,4),5);
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,5),5); //loc
            BOOST_CHECK_EQUAL(mmpr_ve_id_at_owner(mesh_id,6),7);

        }
#endif // PARALLEL
}

BOOST_AUTO_TEST_SUITE_END()
