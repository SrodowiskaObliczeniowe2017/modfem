#define BOOST_TEST_MODULE ModFEM_pcs_mpi_safe_test
#include <boost/test/included/unit_test.hpp>

#include "../include/pch_intf.h"
#include "../include/mmph_intf.h"
#include "../include/uth_intf.h"
#include "../include/mmh_intf.h"
#include "../include/uth_log.h"

const int PCC_MSG_ID_TEST = 1351;
const int PCC_MSG_ID_TEST_BUF =1333;


int prev_id= 0;
int next_id= 0;

char  io_name[1024] = "modfem_test_mpi_out.txt";
FILE* f_out=NULL;
int n_proc=0, my_proc_id=0;

BOOST_AUTO_TEST_SUITE( ModFEM_pcs_mpi_safe_test )

BOOST_AUTO_TEST_CASE( init_communication_testing )
{

    pcr_init_parallel(& boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv,
                      "",io_name,
                      &f_out, & n_proc, &my_proc_id);

    prev_id = pcv_my_proc_id == 1 ? pcv_nr_proc : (pcv_my_proc_id-1);
    next_id = pcv_my_proc_id == pcv_nr_proc ? 1 : (pcv_my_proc_id+1);

//    pcr_exit_parallel();

}

BOOST_AUTO_TEST_CASE( communication_testing_2_point )
{
    /// Automatic testing procedure.
    /// pcr_test is intented to mf_test some selected cases, and assure, that
    /// 1) pcl_mpi_safe works in this environment
    /// 2) error handling is correct

        if(pcv_nr_proc < 2)
            return;

        // 1. Basic mf_test:
        // 2-point communication.

            int number_in = 658687;
            int num_array_in[10] = {646,64684,64,686,545,7475646,44242,7656,43434,64685};
            if(pcv_my_proc_id == 1 ) {
                const int other_id = (pcv_my_proc_id+1)%(pcv_nr_proc + 1);
                mf_check(other_id > 0 && other_id <= pcv_nr_proc,"Wrong other id!");
                // Plain.
                mf_test(PCC_OK == pcr_send_int(other_id,PCC_MSG_ID_TEST,1,& number_in), "Plain communication mf_test failed!");
                // Buffered.
                const int buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF,other_id,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                int numbers_out[10] = {0};
                mf_test(PCC_OK == pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF,buf,10,numbers_out), "Unpacking failed!");
                mf_test(PCC_OK == pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF,buf), "Buffer closeing failed!");

                for(int i=0; i < 10; ++i) {
                    mf_test(numbers_out[i] == num_array_in[i], "Wrong numbers unpacked!");
                }

            }
            else if(pcv_my_proc_id == 2){
                // Plain.
                int number_out = 0;
                const int other_id = (pcv_my_proc_id-1)%(pcv_nr_proc + 1);
                pcr_receive_int(other_id,PCC_MSG_ID_TEST,1,& number_out);
                mf_test(number_in == number_out, "Failed basic 2-point communication mf_test!");
                // Buffered.
                const int buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF,buf,10,num_array_in),"Packig failed!");
                mf_test(PCC_OK == pcr_buffer_send(PCC_MSG_ID_TEST_BUF,buf,other_id), "Sending failed!");
            }

        pcr_barrier();
        mf_log_info("pcl_mpi_safe: Basic test passed.");
        //pcr_exit_parallel();
}

BOOST_AUTO_TEST_CASE( communication_testing_reduction )
{
        // 2. Advanced.
        // Reduction


            int number_in=pcv_my_proc_id;
            int number_out = 0;

            pcr_allreduce_max_int(1,&number_in, &number_out);
            mf_test(number_out== pcv_nr_proc, "Wrong result of pcr_allreduce_max_int!  (result=%d, expected %d)", number_out,  pcv_nr_proc );

            pcr_allreduce_sum_int(1, &number_in, &number_out);
            int sum = 0;
            for(int i=1; i <= pcv_nr_proc; ++i) {
                sum += i;
            }
            mf_test(number_out == sum, "Wrong result of pcr_allreduce_sum_int (result=%d, expected %d) !", number_out,  sum);

    //        delete [] number_in;
    //        delete [] number_out;
//            pcr_exit_parallel();

}

BOOST_AUTO_TEST_CASE( communication_testing_reduction_token_ring)
{
        // Token ring.

            // Send a token with value = 0.
            // Each node increases val by 1.
            const int n_total_loops = 5;
            int n_loops = n_total_loops; // (this is random value)
            // Token have to make n_loops rounds.
            int token = 0;


            if(pcv_my_proc_id == PCC_MASTER_PROC_ID) {
                // Plain.
                // sending token!
                pcr_send_int(next_id,PCC_MSG_ID_TEST,1, &token);

                do {
                    pcr_receive_int(prev_id,PCC_MSG_ID_TEST,1, &token);
                    ++token;
                    pcr_send_int(next_id,PCC_MSG_ID_TEST,1, &token);
                } while(--n_loops);

                pcr_receive_int(prev_id,PCC_MSG_ID_TEST,1, &token);
                ++token;

                mf_test(token == ((n_total_loops+1) * pcv_nr_proc ), "Token ring test value have error!");
            }
            else { // Other proc then master.
                // Plain.
                token=0;
                 do{
                    pcr_receive_int(prev_id,PCC_MSG_ID_TEST,1, &token);
                    ++token;
                    pcr_send_int(next_id,PCC_MSG_ID_TEST,1, &token);
                }while(n_loops--);
            }

            // Buffered.
            // Send a token with value = 0.
            // Each node increases val by 1 and adds current val into buffer (so buffer became bigger and bigger).
            n_loops = n_total_loops;

            // loop
            int expected_size=1;
            std::vector<int> numbers(PCC_DEFAULT_BUFFER_SIZE,0);
            numbers[0]=0;
            int buf=-1;
            int l=0;
            for(;l<n_total_loops;++l) {
                // omit at start
                if(pcv_my_proc_id != PCC_MASTER_PROC_ID
                        || l != 0) {
                    numbers.clear();
                    buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF,prev_id,PCC_DEFAULT_BUFFER_SIZE);
                    expected_size = l*pcv_nr_proc + (pcv_my_proc_id-1);
                    mf_test(expected_size > 0, "Invalid expected size(%d)!",expected_size);
                    numbers.resize(expected_size+1);
                    pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF,buf,expected_size, numbers.data());
                    pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF,buf);

                    token = numbers[expected_size-1];
                    mf_test(token == (expected_size-1), "Buffered token ring mf_test intermidate value error!");
                    ++token;
                    ++expected_size;
                    numbers[expected_size-1]=token;
                }
                // omit at ending
                //if(pcv_my_proc_id != pcv_nr_proc || l < (n_total_loops-1)) {
                    buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF,PCC_DEFAULT_BUFFER_SIZE);
                    pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF,buf,expected_size,numbers.data());
                    pcr_buffer_send(PCC_MSG_ID_TEST_BUF,buf,next_id);
                //}
            };

            // handle last request
            if(pcv_my_proc_id == PCC_MASTER_PROC_ID) {
                numbers.clear();
                buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF,prev_id,PCC_DEFAULT_BUFFER_SIZE);
                expected_size = l*pcv_nr_proc + (pcv_my_proc_id-1);
                mf_test(expected_size > 0, "Invalid expected size(%d)!",expected_size);
                numbers.resize(expected_size+1);
                pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF,buf,expected_size, numbers.data());
                pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF,buf);

                token = numbers[expected_size-1];
                mf_test(token == (expected_size-1), "Buffered token ring mf_test intermidate value error!");
                ++token;
                ++expected_size;
                numbers[expected_size-1]=token;

                // checking results
                for(int ll=0; ll < n_total_loops; ++ll) {
                    for(int p=0; p < pcv_nr_proc; ++p) {
                        mf_test(numbers[ll*pcv_nr_proc+p] == ll*pcv_nr_proc+p,"Token ring buffer mf_test value have error!");
                    }
                }

            }




        pcr_barrier();
        mf_log_info("pcl_mpi_safe: Advanced test passed.");

        // pcr_exit_parallel();
}

BOOST_AUTO_TEST_CASE( communication_testing_B_cast )
{
        // 4. Bcast
        // Plain.

            int numbers[5]={0};
            int numbers2[5]={84648,684684,684684648,64882,22222};
            if(pcv_my_proc_id == PCC_MASTER_PROC_ID) {
                std::copy(numbers2,numbers2+5,numbers);
            }
            pcr_bcast_int(PCC_MASTER_PROC_ID,5,numbers);


             BOOST_CHECK_EQUAL_COLLECTIONS(numbers,numbers+5,numbers2,numbers2+5);

            mf_log_info("pcl_mpi_safe: plain Bcast test passed.");

            //pcr_exit_parallel();
}
BOOST_AUTO_TEST_CASE( communication_testing_B_cast_buffered )
{

        // Buffered
    int numbers[5]={0};
    int numbers2[5]={84648,684684,684684648,64882,22222};
            int buf = -1;
            if(pcv_my_proc_id == PCC_MASTER_PROC_ID) {
                buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF,PCC_DEFAULT_BUFFER_SIZE);
                pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF,buf,5,numbers);
            }
            int numbers3[5]={0};
            int buf_out = pcr_buffer_bcast(PCC_MSG_ID_TEST_BUF,buf,PCC_MASTER_PROC_ID);
            if(pcv_my_proc_id  != PCC_MASTER_PROC_ID) {
                pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF,buf_out,5,numbers3);
                BOOST_CHECK_EQUAL_COLLECTIONS(numbers,numbers+5,numbers3,numbers3+5);
            }


        pcr_barrier();
        mf_log_info("pcl_mpi_safe: buffered Bcast test passed.");


}

BOOST_AUTO_TEST_CASE( communication_testing_2_point_out_of_order )
{
    /// Automatic testing procedure.
    /// pcr_test is intented to mf_test some selected cases, and assure, that
    /// 1) pcl_mpi_safe works in this environment
    /// 2) error handling is correct

        if(pcv_nr_proc < 2)
            return;

        // 1. Basic mf_test:
        // 2-point communication.

            int number_in = 658687;
            int num_array_in[10] = {646,64684,64,686,545,7475646,44242,7656,43434,64685};
            if(pcv_my_proc_id == 1 ) {
                const int other_id = (pcv_my_proc_id+1)%(pcv_nr_proc + 1);
                mf_check(other_id > 0 && other_id <= pcv_nr_proc,"Wrong other id!");
                // Plain.
                double tmp = 19283713.123123;
                // mixing different types and ids
                mf_test(PCC_OK == pcr_send_double(other_id,PCC_MSG_ID_TEST*2,1,& tmp), "Plain communication mf_test failed!");
                mf_test(PCC_OK == pcr_send_int(other_id,PCC_MSG_ID_TEST*2,1,& number_in), "Plain communication mf_test failed!");
                mf_test(PCC_OK == pcr_send_int(other_id,PCC_MSG_ID_TEST,1,& number_in), "Plain communication mf_test failed!");
                mf_test(PCC_OK == pcr_send_double(other_id,PCC_MSG_ID_TEST,1,& tmp), "Plain communication mf_test failed!");

                // Buffered.
                 int buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF*2,other_id,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                int numbers_out[10] = {0};
                mf_test(PCC_OK == pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF*2,buf,10,numbers_out), "Unpacking failed!");
                mf_test(PCC_OK == pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF*2,buf), "Buffer closeing failed!");

                // Buffered 2.
                  buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF*4,PCC_ANY_PROC,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF*4,buf,10,numbers_out), "Unpacking failed!");
                mf_test(PCC_OK == pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF*4,buf), "Buffer closeing failed!");

                // Buffered 3.
                  buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF*3,other_id,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF*3,buf,10,numbers_out), "Unpacking failed!");
                mf_test(PCC_OK == pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF*3,buf), "Buffer closeing failed!");

                // Buffered 4.
                  buf = pcr_buffer_receive(PCC_MSG_ID_TEST_BUF,PCC_ANY_PROC,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_unpack_int(PCC_MSG_ID_TEST_BUF,buf,10,numbers_out), "Unpacking failed!");
                mf_test(PCC_OK == pcr_recv_buffer_close(PCC_MSG_ID_TEST_BUF,buf), "Buffer closeing failed!");

                for(int i=0; i < 10; ++i) {
                    mf_test(numbers_out[i] == num_array_in[i], "Wrong numbers unpacked!");
                }

            }
            else if(pcv_my_proc_id == 2){
                // Plain.
                int number_out = 0;
                const int other_id = (pcv_my_proc_id-1)%(pcv_nr_proc + 1);
                pcr_receive_int(other_id,PCC_MSG_ID_TEST,1,& number_out);
                pcr_receive_int(PCC_ANY_PROC,PCC_MSG_ID_TEST*2,1,& number_out);
                mf_test(number_in == number_out, "Failed basic 2-point communication mf_test!");

                double recv_tmp=0.0;
                pcr_receive_double(other_id,PCC_MSG_ID_TEST,1,& recv_tmp);
                pcr_receive_double(PCC_ANY_PROC,PCC_MSG_ID_TEST*2,1,& recv_tmp);

                // Buffered 1.
                 int buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF,buf,10,num_array_in),"Packig failed!");
                mf_test(PCC_OK == pcr_buffer_send(PCC_MSG_ID_TEST_BUF,buf,other_id), "Sending failed!");

                // Buffered 2.
                  buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF*2,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF*2,buf,10,num_array_in),"Packig failed!");
                mf_test(PCC_OK == pcr_buffer_send(PCC_MSG_ID_TEST_BUF*2,buf,other_id), "Sending failed!");

                // Buffered 3.
                  buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF*3,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF*3,buf,10,num_array_in),"Packig failed!");
                mf_test(PCC_OK == pcr_buffer_send(PCC_MSG_ID_TEST_BUF*3,buf,other_id), "Sending failed!");

                // Buffered 4.
                  buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF*4,PCC_DEFAULT_BUFFER_SIZE);
                mf_test(buf > -1, "Buffer creation faied!");
                mf_test(PCC_OK == pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF*4,buf,10,num_array_in),"Packig failed!");
                mf_test(PCC_OK == pcr_buffer_send(PCC_MSG_ID_TEST_BUF*4,buf,other_id), "Sending failed!");
            }

        pcr_barrier();
        mf_log_info("pcl_mpi_safe: Basic test passed.");
        //pcr_exit_parallel();
}

    // These tests below are disabled by default, since, they checks what happens when error is introduced by user.
    // If those errors will be enabled, the pcr_test() show correct information about error, but program will terminate.
    // It can be changed, if needed.
    //    // 5. Forced overflow.
    //    {
    //        std::unique_ptr<int>  big_data(new int[PCC_DEFAULT_BUFFER_SIZE]);
    //        int result = PCC_OK;

    //        // Plain
    //        for(int i=0; i < 1000; ++i) {
    //            result = pcr_send_int(next_id,PCC_MSG_ID_TEST,PCC_DEFAULT_BUFFER_SIZE,big_data.get());
    //        }
    //        mf_test(result != PCC_OK, "Unexpected overflow mf_test success; it has to fail!");

    //        // Buffered.
    //        int buf = pcr_send_buffer_open(PCC_MSG_ID_TEST_BUF,PCC_DEFAULT_BUFFER_SIZE);
    //        for(int i=0; i < 10000; ++i) {
    //            result = pcr_buffer_pack_int(PCC_MSG_ID_TEST_BUF,buf,PCC_DEFAULT_BUFFER_SIZE,big_data.get());
    //        }
    //        mf_test(result != PCC_OK, "Unexpected buffer overflow mf_test success; it has to fail!");
    //        mf_test(PCC_OK != pcr_buffer_send(PCC_MSG_ID_TEST_BUF,buf,next_id),"Unexpected buffer overflow mf_test success; it has to fail!");

    //    }
    //    mf_log_info("pcl_mpi_safe: Forced overflow test passed.");

    //#if (__cplusplus > 199711L) // for <random> support
    //    {
    //        // 3. Invalid usage.
    //        // Plain.
    //        std::default_random_engine generator;
    //        std::uniform_int_distribution<int> distribution(-10000,10000000);
    //        int dest=distribution(generator),
    //                msg=distribution(generator),
    //                msg2 = distribution(generator),
    //                n_num=distribution(generator),
    //                nums=distribution(generator);
    //        mf_test(PCC_OK != pcr_send_int(dest,msg,n_num,&nums),"Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_receive_int(dest,msg2,n_num,&nums), "Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_send_int(dest,msg,n_num,NULL),"Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_receive_int(dest,msg2,n_num,NULL), "Passed mf_test, which is intended to fail!");

    //        // Buffered.
    //        int buf = pcr_send_buffer_open(msg,n_num);
    //        mf_test(PCC_OK != pcr_buffer_pack_int(msg2,buf,n_num,&nums),"Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_buffer_send(msg2,buf,dest),"Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_buffer_pack_int(msg2,buf,n_num,&nums),"Passed mf_test, which is intended to fail!");

    //        buf = pcr_send_buffer_open(msg,n_num);
    //        mf_test(PCC_OK != pcr_buffer_pack_int(msg2,buf,n_num,NULL),"Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_buffer_send(msg2,buf,dest),"Passed mf_test, which is intended to fail!");
    //        mf_test(PCC_OK != pcr_buffer_pack_int(msg2,buf,n_num,NULL),"Passed mf_test, which is intended to fail!");
    //    }
    //    mf_log_info("pcl_mpi_safe: Invalid usage test passed.");
    //#endif


BOOST_AUTO_TEST_CASE( communication_end )
{
        pcr_exit_parallel();
}

BOOST_AUTO_TEST_SUITE_END()
