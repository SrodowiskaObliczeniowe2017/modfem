/*
 * fv_threads.h
 *
 *  Created on: 31-03-2012
 *      Author: Paweł Macioł
 *  Description:
 *  	Simple library for creation and managment of threads.
 */

#ifndef fv_threads_H_
#define fv_threads_H_

/* Defines */
#ifdef _WIN32
  // Windows threads use different data structures
  #include <windows.h>
  typedef HANDLE   fv_thread;
  typedef unsigned fv_thread_id;
  typedef unsigned (__stdcall * fv_thread_routine)(void *);
  typedef CRITICAL_SECTION fv_critical_section;
  
  #define fv_thread_return unsigned __stdcall
  #define INITCRITICALSECTION(crit_section)	 InitializeCriticalSection(&(crit_section));
  #define DELETECRITICALSECTION(crit_section) DeleteCriticalSection(&(crit_section));
  #define ENTERCRITICALSECTION(crit_section) EnterCriticalSection(&(crit_section));
  #define LEAVECRITICALSECTION(crit_section) LeaveCriticalSection(&(crit_section));

#else

  // Includes POSIX thread headers for Linux thread support
  #include <pthread.h>
  typedef pthread_t fv_thread;
  typedef pthread_t fv_thread_id;
  typedef void *(*fv_thread_routine)(void *);
  typedef pthread_mutex_t fv_critical_section;
  typedef pthread_cond_t  fv_signal_t;
  
  #define fv_thread_return void*
  #define INITCRITICALSECTION(crit_section)	  pthread_mutex_init(&(crit_section),NULL);
  #define DELETECRITICALSECTION(crit_section) pthread_mutex_destroy(&(crit_section));
  #define ENTERCRITICALSECTION(crit_section)  pthread_mutex_lock  (&(crit_section));
  #define LEAVECRITICALSECTION(crit_section)  pthread_mutex_unlock(&(crit_section));
  #define INITSIGNAL(signal_var) 			  pthread_cond_init(&(signal_var), NULL);
  #define DELETESIGNAL(signal_var)			  pthread_cond_destroy(&(signal_var));

#endif

#ifdef __cplusplus

  /*struct Event
  {
    protected:
	  void* _handle;

    public:
	  explicit Event();
	  ~Event();

	  bool Wait(const unsigned int time = -1);
	  bool TryWait();
	  void Signal();
  };
  */

/* Very simple lock */
  class Lock
  {
    protected:
	  fv_critical_section _sec;
	  bool _is;

    public:
	  Lock() : _is(false) {
		  INITCRITICALSECTION(this->_sec);
	  }
	  ~Lock() {
		  DELETECRITICALSECTION(this->_sec);
	  }

	  void lock() {
		  if (!_is) {
			  ENTERCRITICALSECTION(this->_sec);
			  _is = true;
		  }
	  }
	  void unlock() {
		  _is = false;
		  LEAVECRITICALSECTION(this->_sec);
	  }
  };


  class Thread
  {
    public:
	  Thread();
	  virtual ~Thread();

	  virtual void Start();
	  virtual void Stop();
	  virtual int  Run();
	  virtual int  Check() const { return 1; }

	  static void* ThreadFunc(void *data);

	  void Terminate();
	  void WaitForNotify();
	  int Notify();

    protected:
	  fv_thread		_handle;
	  fv_thread_id	_id;
	  fv_critical_section _crtsec;
	  fv_signal_t	_signal;

  };



    extern "C" {
#endif

    extern int utd_internal_error;

    //Create thread.
    extern fv_thread fv_start_thread(fv_thread_id* thread_id, fv_thread_routine func, void *data);

    //Wait for thread to finish.
    extern void fv_end_thread(fv_thread thread);

    //Destroy thread.
    extern void fv_destroy_thread(fv_thread thread);

    //Wait for multiple threads.
    extern void fv_wait_for_threads(const fv_thread *threads, int num);

    //Waiting for wake signal
    extern int fv_wait_for_wakeup(fv_signal_t *signal,fv_critical_section *critical);

    //Sending wake signal
    extern int fv_wakeup(fv_signal_t *signal);

#ifdef __cplusplus
} //extern "C"
#endif


#endif /* fv_threads_H_
  */
