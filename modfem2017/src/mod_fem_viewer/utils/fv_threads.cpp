/*
 * fv_threads.c
 *
 *  Created on: 31-03-2012
 *      Author: Paweł Macioł
 * Description:
 * 		Implementation.
 */

#include "uth_log.h"
#include "fv_threads.h"
#include "fv_exception.h"
#include <iostream>

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

    Thread::Thread() : _handle(0), _id(0)
    {
    	INITCRITICALSECTION(_crtsec);
    	INITSIGNAL(_signal);
    }

    Thread::~Thread()
    {
    	DELETECRITICALSECTION(_crtsec);
    	DELETESIGNAL(_signal);
    }

    void Thread::Start()
    {
    	//mfp_debug("Just before running new thread");
    	_handle = fv_start_thread(&this->_id,this->ThreadFunc,this);
    	if (_handle == 0) {
    		throw "Can't start thread\n";
    	}
    }

    int Thread::Run()
    {
    	return 0;
    }

    void Thread::Stop()
    {
    	fv_wait_for_wakeup(&_signal,&_crtsec);
    }

    void* Thread::ThreadFunc(void* data)
    {
    	try {
    		reinterpret_cast<Thread*>(data)->Run();
    	}
    	catch (const fv_exception &e) {
    		std::cerr << e.what();
    		return (void*)-1;
    	}

    	return (void*)0;
    }

    void Thread::Terminate()
    {
    	fv_destroy_thread(_handle);
    	//fv_end_thread(_handle);
    }

    void Thread::WaitForNotify()
    {
    	fv_wait_for_wakeup(&_signal,&_crtsec);
    }

    int Thread::Notify()
    {
    	return fv_wakeup(&_signal);
    }



