#ifndef AppImpl_H_
#define AppImpl_H_

#include "Log.h"
#include "fv_threads.h"
#include "fv_config.h"
#include "pch_intf.h"

#include <unistd.h>
#include <memory>

namespace FemViewer {

/* The instatnt type */
enum InstantType {
	MODULE,
	APPLICATION
};

/* Types of module */
enum CoreType {
	ALONE = 0,
	MASTER,
	SLAVE,
	COUNT
};


/* Forward declaration of start function */
template< typename T,CoreType Type > struct Start;

int run_application(int argc, char** argv, bool initgl, Thread* parent = nullptr);

class ModuleThread : public Thread
{
public:
	int 	_shutdown;
	int  	_argc;
	char** 	_argv;
	bool	_running;

	explicit ModuleThread(int argc,char** argv);
	virtual ~ModuleThread();

	int Check() const { return _shutdown; }

	void Start();
	int  Run();
	int  Update();
	void Stop();
};

/* AppImpl class for application */
template<class T> class AppImpl
{
	typedef T Core;

  public:

	static int run(int argc,char **argv)
	{
		//mfp_log_debug("W inicjajcji2\n");
		return Core::Init(argc,argv);
	}
};

/* Partial specialization of AppImpl for module instances */
template<class T>
class AppImpl<T*>
{
	typedef T Core;

  public:

	static int run(int argc,char** argv);
	static void close(void);

	ModuleThread _thread;
	~AppImpl() {}

  protected:
	static std::unique_ptr<AppImpl<T*> > _self;

	AppImpl(int argc,char** argv) : _thread(argc,argv) {
		mfp_log_debug("ctr of pointer types\n");
		atexit(close);
	}

  private:
	AppImpl(const AppImpl&);
	void operator=(const AppImpl&);
};


template<class T>
std::unique_ptr<AppImpl<T*> > AppImpl<T*>::_self;

/* Start function for STAND_ALONE module */
template< typename T >
struct Start<T,ALONE> {
	static fv_thread_return run(void *pParams)
	{
		AppImpl<T*> * lptr = reinterpret_cast<AppImpl<T*> *>(pParams);
		// Try init logger
		T::Init(lptr->_argc,lptr->_argv, lptr);
		ENTERCRITICALSECTION(lptr->_crtsec);
		//Log::GetInstance().Close();
		//mfp_debug("Wake up main thread\n");
		fv_wakeup(&(lptr->_signal));
		LEAVECRITICALSECTION(lptr->_crtsec);
		lptr->destroy_self();
		return (fv_thread_return)0;
	}
};

/* Start function for PARALLEL_MASTER module */
template< typename T >
struct Start<T,MASTER> {
	static fv_thread_return run(void *pParams) {
		AppImpl<T> * lptr = reinterpret_cast<AppImpl<T> *>(pParams);
		T::Init(lptr->_argc,lptr->_argv);
		//Log::GetInstance().Close();
		lptr->close();
		return (fv_thread_return)0;
	}
};

/* Start function for PARALLEL_SLAVES modules */
template< typename T >
struct Start<T,SLAVE> {
	static fv_thread_return run(void *pParams) {
		AppImpl<T*> * lptr = static_cast<AppImpl<T*> *>(pParams);
		T::Init2(lptr->_argc,lptr->_argv);
		//mfvWindow::init(lptr->_argc,lptr->_argv);
		//Log::GetInstance().Close();
		lptr->destroy_self();
		return (fv_thread_return)0;
	}
};

#define NO_GUI 1

template<class T>
int AppImpl<T*>::run(int argc, char **argv)
{
	if (!_self.get())
	{
		try {
			mfp_log_debug("init graphics\n");
			_self.reset(new AppImpl<T*>(argc,argv));
			_self->_thread.Start();
		}
		catch (const std::bad_alloc & e) {
			std::string msg("Error during graphic application init: ");
			msg.append( e.what());
			throw msg.c_str();
		}
	}
	else {
		close();
	}


	return 0;
}

template<class T>
void AppImpl<T*>::close(void)
{
	mfp_debug("close\n");
	if (_self.get() != nullptr && !(_self->_thread._shutdown))
	{
		// Signal end of the work
		_self->_thread.Stop();

		// Sllep for one msecond and destroy kernel
		::sleep(100);
		mfp_debug("before releasing\n");
		_self.release();
	}
}


}// end namespace FemViewer

#endif /* AppImpl_H_ */
