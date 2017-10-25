/*
 * fv_threads.c
 *
 *  Created on: 31-03-2012
 *      Author: Paweł Macioł
 * Description:
 * 		Implementation.
 */

#include "fv_threads.h"
#ifdef _WIN32
# include <process.h>
#endif

int fv_internal_error;

#if _WIN32
    //Create thread
    fv_thread fv_start_thread(fv_thread_id thread_id, fv_thread_routine func, void *data){
        return _beginthreadex(0, 0, (unsigned (_stdcall *)(void*))func, data, 0, thread_id);
    }

    //Wait for thread to finish
    void fv_end_thread(fv_thread thread){
        WaitForSingleObject(thread, INFINITE);
        CloseHandle(thread);
    }

    //Destroy thread
    void fv_destroy_thread(fv_thread thread){
        TerminateThread(thread, 0);
        CloseHandle(thread);
    }

    //Wait for multiple threads
    void fv_wait_for_threads(const fv_thread * threads, int num) {
		int i;
        WaitForMultipleObjects(num, threads, TRUE, INFINITE);

        for(i = 0; i < num; i++)
            CloseHandle(threads[i]);
    }

#else
    //Create thread
    fv_thread fv_start_thread(fv_thread_id *thread_id,fv_thread_routine func, void * data){
        fv_internal_error = pthread_create(thread_id, NULL, func, data);
        return *thread_id;
    }

    //Wait for thread to finish
    void fv_end_thread(fv_thread thread){
        pthread_join(thread, NULL);
    }

    //Destroy thread
    void fv_destroy_thread(fv_thread thread){
        pthread_cancel(thread);
    }

    //Wait for multiple threads
    void fv_wait_for_threads(const fv_thread* threads, int num){
        int i=0; for(; i < num; i++)
        	fv_end_thread(threads[i]);
    }

    int fv_wait_for_wakeup(fv_signal_t *signal,fv_critical_section *critical) {
    	return pthread_cond_wait(signal,critical);
    }

    int fv_wakeup(fv_signal_t *signal) {
    	return pthread_cond_signal(signal);
    }

#endif

