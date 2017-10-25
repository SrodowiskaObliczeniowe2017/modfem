#define BOOST_TEST_MODULE ModFEM_ddl_parmetis_test
#include <boost/test/included/unit_test.hpp>

#include <boost/filesystem.hpp>
#include <vector>

#include "ddh_intf.h"
#include "mmh_intf.h"
#include "uth_io_files.h"
#include "uth_io_compression.h"
#include "uth_log.h"

std::string ctest_dir = "ctest/modfem_test/";

BOOST_AUTO_TEST_SUITE( ModFEM_ddl_parmetis_test )

BOOST_AUTO_TEST_CASE( ddl_parmetis_test_setup )
{
    int n_loops=0;
    while( !boost::filesystem::exists(ctest_dir+"test_beacon.txt") && (n_loops<5) ) {
        ctest_dir = "../"+ctest_dir;
        ++n_loops;
    }

    BOOST_CHECK_MESSAGE( boost::filesystem::exists(ctest_dir+"test_beacon.txt"),
                         "Invalid path: " << ctest_dir+"test_beacon.txt");
}


BOOST_AUTO_TEST_CASE( simple_domain_decomposition_2_dom )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
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

    BOOST_CHECK(mesh_id == 1);

    const int n_sdm = 2;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  6);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);


    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),NULL,NULL,NULL);

    BOOST_WARN(n_sdm_sizes[0] == 3);
    BOOST_WARN(n_sdm_sizes[1] == 3);

}

BOOST_AUTO_TEST_CASE( simple_domain_decomposition_3_dom )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
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

    BOOST_CHECK(mesh_id == 2);

    const int n_sdm = 3;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  6);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);


    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),NULL, NULL,NULL);

    BOOST_WARN_MESSAGE(n_sdm_sizes[0] == 2, "Why " << n_sdm_sizes[0] << "? Should be 2.");
    BOOST_WARN_MESSAGE(n_sdm_sizes[1] == 2, "Why " << n_sdm_sizes[1] << "? Should be 2.");
    BOOST_WARN_MESSAGE(n_sdm_sizes[2] == 2, "Why " << n_sdm_sizes[2] << "? Should be 2.");

}

BOOST_AUTO_TEST_CASE( simple_domain_decomposition_6_dom )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
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

    BOOST_CHECK(mesh_id == 3);

    const int n_sdm = 6;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  6);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);


    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),NULL, NULL,NULL);

    BOOST_WARN(n_sdm_sizes[0] == 1);
    BOOST_WARN(n_sdm_sizes[1] == 1);
    BOOST_WARN(n_sdm_sizes[2] == 1);
    BOOST_WARN(n_sdm_sizes[3] == 1);
    BOOST_WARN(n_sdm_sizes[4] == 1);
    BOOST_WARN(n_sdm_sizes[5] == 1);

}

BOOST_AUTO_TEST_CASE( simple_domain_decomposition_2_dom_2_prism )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
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
    /// !MODFEM_NEW_MPI

    BOOST_CHECK(mesh_id == 4);

    const int n_sdm = 2;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  2);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);


    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),NULL, NULL,NULL);

    BOOST_WARN(n_sdm_sizes[0] == 1);
    BOOST_WARN(n_sdm_sizes[1] == 1);

}

BOOST_AUTO_TEST_CASE( dd_2_dom )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.

    std::string filename = "LDC_z_001_100x100.nas",
                workdir = ctest_dir+"../../meshes/cube/LDC/nas/";

    char decompressed_file[255]={0};
    utr_io_decompress_file(workdir.c_str(),filename.c_str(),decompressed_file);
    filename = decompressed_file;


    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     workdir.c_str(),
                                     filename.c_str(),
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

    BOOST_CHECK(mesh_id == 5);

    const int n_sdm = 2;
    const int n_e = mmr_get_nr_elem(mesh_id);

//    BOOST_CHECK(n_e ==  6);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);


    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),NULL, NULL,NULL);

    BOOST_WARN(n_sdm_sizes[0] == n_e/2);
    BOOST_WARN(n_sdm_sizes[1] == n_e/2);

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( ModFEM_ddl_parmetis_test_overlap )

BOOST_AUTO_TEST_CASE( simple_domain_decomposition_2_dom_2_prism_overlap )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
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
    /// !MODFEM_NEW_MPI

    BOOST_CHECK(mesh_id == 6);

    const int n_sdm = 2;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  2);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);


    std::vector<int> overlap_sizes, overlap_elems;

    overlap_sizes.resize(n_sdm);
    overlap_elems.resize(n_e);

    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),
                          overlap_sizes.data(), overlap_elems.data(),NULL);

    BOOST_CHECK(n_sdm_sizes[0] == 1);
    BOOST_CHECK(n_sdm_sizes[1] == 1);

    BOOST_CHECK(elems_ids[0] == 1);
    BOOST_CHECK(elems_ids[1] == 2);

    BOOST_CHECK(overlap_sizes[0] == 0);
    BOOST_CHECK(overlap_sizes[1] == 1);

    // first proc do not have overlap, but second do
    BOOST_CHECK_MESSAGE(overlap_elems[0] == 1, "Overlap elem should be 1, but is " << overlap_elems[0]);

    int* core = elems_ids.data();
    int core_size=0;
    int *overlap = overlap_elems.data();
    int overlap_size=0;
    std::vector<int>    result(overlap_elems.size());
    int* presult = result.data();

    for(int p=0; p < n_sdm; ++p) {
        core_size = n_sdm_sizes[p];
        overlap_size = overlap_sizes[p];

        BOOST_CHECK(core_size > 0);
        if(p==0) {
            BOOST_CHECK(overlap_size == 0);
        }
        else {
            BOOST_CHECK(overlap_size > 0);
        }
        // intersection should be empty
        BOOST_CHECK(presult == std::set_intersection(core, core+core_size,
                                        overlap, overlap+overlap_size,
                                        presult ) );
        core+=core_size;
        overlap+=overlap_size;
    }

}

BOOST_AUTO_TEST_CASE( simple_domain_decomposition_2_dom_overlap )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
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

    BOOST_CHECK(mesh_id == 7);

    const int n_sdm = 2;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  6);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);

    std::vector<int> overlap_sizes, overlap_elems;

    overlap_sizes.resize(n_sdm);
    overlap_elems.resize(n_e);

    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),
                          overlap_sizes.data(), overlap_elems.data(),NULL);

    BOOST_WARN(n_sdm_sizes[0] == 3);
    BOOST_WARN(n_sdm_sizes[1] == 3);

    BOOST_WARN(overlap_sizes[1] == 3);

    BOOST_CHECK(overlap_sizes[0] == 0);

    int* core = elems_ids.data();
    int core_size=0;
    int *overlap = overlap_elems.data();
    int overlap_size=0;
    std::vector<int>    result(overlap_elems.size());
    int* presult = result.data();

    for(int p=0; p < n_sdm; ++p) {
        core_size = n_sdm_sizes[p];
        overlap_size = overlap_sizes[p];

        BOOST_CHECK(core_size > 0);
        if(p==0) {
            BOOST_CHECK(overlap_size == 0);
        }
        else {
            BOOST_CHECK(overlap_size > 0);
        }
        // intersection should be empty
        BOOST_CHECK(presult == std::set_intersection(core, core+core_size,
                                        overlap, overlap+overlap_size,
                                        presult ) );
        core+=core_size;
        overlap+=overlap_size;
    }


}

BOOST_AUTO_TEST_CASE( simple_domain_decomposition_6_dom_overlap )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.
    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     (ctest_dir+"../../meshes/cube/no_BC/nas/").c_str(),
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

    BOOST_CHECK(mesh_id == 8);

    const int n_sdm = 6;
    const int n_e = mmr_get_nr_elem(mesh_id);

    BOOST_CHECK(n_e ==  6);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);

    std::vector<int> overlap_sizes, overlap_elems;

    overlap_sizes.resize(n_sdm);
    overlap_elems.resize(n_e);

    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),
                          overlap_sizes.data(), overlap_elems.data(),NULL);


    for(int i=0; i < n_sdm; ++i) {
        BOOST_CHECK(n_sdm_sizes[i] == 1);
        BOOST_CHECK(elems_ids[i] == i+1);
    }


    BOOST_CHECK(overlap_sizes[0] == 0);

    BOOST_CHECK(overlap_sizes[1] == 0);

    BOOST_CHECK(overlap_sizes[2] == 2);
    BOOST_CHECK(overlap_elems[0] == 1);
    BOOST_CHECK(overlap_elems[1] == 2);

    BOOST_CHECK(overlap_sizes[3] == 0);

    BOOST_CHECK(overlap_sizes[4] == 2);
    BOOST_CHECK(overlap_elems[2] == 1);
    BOOST_CHECK(overlap_elems[3] == 4);

    BOOST_CHECK(overlap_sizes[5] == 2);
    BOOST_CHECK(overlap_elems[4] == 3);
    BOOST_CHECK(overlap_elems[5] == 5);



    int* core = elems_ids.data();
    int core_size=0;
    int *overlap = overlap_elems.data();
    int overlap_size=0;
    std::vector<int>    result(overlap_elems.size());
    int* presult = result.data();

    for(int p=0; p < n_sdm; ++p) {
        core_size = n_sdm_sizes[p];
        overlap_size = overlap_sizes[p];

        BOOST_CHECK(core_size > 0);

        // intersection should be empty
        BOOST_CHECK(presult == std::set_intersection(core, core+core_size,
                                        overlap, overlap+overlap_size,
                                        presult ) );
        core+=core_size;
        overlap+=overlap_size;
    }
}

BOOST_AUTO_TEST_CASE( dd_tetra_with_overlap )
{
    //std::string dir = ctest_dir+"meshes";
    int mesh_id=-1;
    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.

    std::string filename = "LDC_z_001_100x100.nas",
                workdir = ctest_dir+"../../meshes/cube/LDC/nas/";

    char decompressed_file[255]={0};
    utr_io_decompress_file(workdir.c_str(),filename.c_str(),decompressed_file);
    filename = decompressed_file;


    #ifdef MODFEM_NEW_MPI
      // Using MODFEM_NEW_MPI (mmpd_adapter)
      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
      if(mesh_id > 0) {
       int n_read = utr_io_read_mesh(mesh_id,
                                     workdir.c_str(),
                                     filename.c_str(),
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

    BOOST_CHECK(mesh_id == 9);

    const int n_sdm = 6;
    const int n_e = mmr_get_nr_elem(mesh_id);

    std::vector<int> n_sdm_sizes,elems_ids;

    n_sdm_sizes.resize(n_sdm+1);
    elems_ids.resize(n_e+1);

    std::vector<int> overlap_sizes, overlap_elems;

    overlap_sizes.resize(n_sdm);
    overlap_elems.resize(n_e);

    ddr_create_subdomains_scheme(mesh_id,DDC_DEFAULT,
                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),
                          overlap_sizes.data(), overlap_elems.data(),NULL);

    BOOST_CHECK(overlap_sizes[0] == 0);

    int* core = elems_ids.data();
    int core_size=0;
    int *overlap = overlap_elems.data();
    int overlap_size=0;
    std::vector<int>    result(overlap_elems.size());
    int* presult = result.data();

    for(int p=0; p < n_sdm; ++p) {
        core_size = n_sdm_sizes[p];
        overlap_size = overlap_sizes[p];

        BOOST_CHECK(core_size > 0);
        if(p==0) {
            BOOST_CHECK(overlap_size == 0);
        }
        else {
            BOOST_WARN_MESSAGE(overlap_size > 0, "Empty overlap for subdomain " << p);
        }
        // intersection should be empty
        BOOST_CHECK(presult == std::set_intersection(core, core+core_size,
                                        overlap, overlap+overlap_size,
                                        presult ) );
        core+=core_size;
        overlap+=overlap_size;
    }
}

// TODO: uncommend when jk mesh reading will be fixed
//BOOST_AUTO_TEST_CASE( dd_prism_with_overlap )
//{
//    //std::string dir = ctest_dir+"meshes";
//    int mesh_id=-1;
//    /// Here MODFEM_NEW_MPI define that switches old and new version of I/O handling.

//    std::string filename = "mesh_jk_61x61.dat",
//                workdir = ctest_dir+"../../meshes/cube/LDC/jk/struct/scale1/";


//    #ifdef MODFEM_NEW_MPI
//      // Using MODFEM_NEW_MPI (mmpd_adapter)
//      mesh_id = mmr_init_mesh2(stdout,0,0,0,0);
//      if(mesh_id > 0) {
//       int n_read = utr_io_read_mesh(mesh_id,
//                                     workdir.c_str(),
//                                     filename.c_str(),
//                                     MMC_GRADMESH_DATA);
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

//    BOOST_CHECK(mesh_id == 10);

//    const int n_sdm = 2;
//    const int n_e = mmr_get_nr_elem(mesh_id);

//    BOOST_CHECK(n_e ==  6);

//    std::vector<int> n_sdm_sizes,elems_ids;

//    n_sdm_sizes.resize(n_sdm+1);
//    elems_ids.resize(n_e+1);

//    std::vector<int> overlap_sizes, overlap_elems;

//    overlap_sizes.resize(n_sdm);
//    overlap_elems.resize(n_e);

//    ddr_create_subdomains(mesh_id,DDC_DEFAULT,
//                          n_sdm,n_sdm_sizes.data(),elems_ids.data(),
//                          overlap_sizes.data(), overlap_elems.data(),NULL);

//BOOST_CHECK(overlap_sizes[0] == 0);

//int* core = elems_ids.data();
//int core_size=0;
//int *overlap = overlap_elems.data();
//int overlap_size=0;
//std::vector<int>    result(overlap_elems.size());
//int* presult = result.data();

//for(int p=1; p < n_sdm; ++p) {
//    core_size = n_sdm_sizes[p];
//    overlap_size = overlap_sizes[p];

//    BOOST_CHECK(core_size > 0);
//    BOOST_CHECK(overlap_size > 0);

//    // intersection should be empty
//    BOOST_CHECK(presult == std::set_intersection(core, core+core_size,
//                                    overlap, overlap+overlap_size,
//                                    presult ) );
//    core+=core_size;
//    overlap+=overlap_size;
//}

//}

BOOST_AUTO_TEST_SUITE_END()
