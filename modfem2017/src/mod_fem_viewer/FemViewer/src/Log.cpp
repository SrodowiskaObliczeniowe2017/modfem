#include <cstdarg>
#include <cassert>
#include <string>
#include<iostream>
/* parallel comunication */
//#ifdef PARALELL
#ifdef PARALLEL
#include <mpi.h>
#include "pch_intf.h"
#endif
#include "fv_config.h"
#include "Log.h"
#include "WindowFV.h"
// namespaces
namespace FemViewer {

// initializations
std::unique_ptr<Log> Log::m_selfp(nullptr);
const eLevel Log::m_dfErrLevel = MAX_LEVEL_LOG;
const int 	 Log::m_dfLogDest = /*LOGCONSOLE |*/ LOGFILE;
const int MAX_LOG_STRING = 1024;

FILE * Log::m_dfLogFile = fopen("execution.txt", "w");
FILE * Log::m_dfLogConsole = stderr;


bool Log::Init(FILE* fp,mfvWindow* dlgp)
{
	if (m_selfp.get() == nullptr)
	{
		FILE *lfp, *lcp;
		int mode  = m_dfLogDest;
		// Configure logger's mode
		if (fp == stdout || fp == stderr) {
			mode |= LOGCONSOLE;
			lcp = fp;
			lfp = m_dfLogFile;
			assert(lfp != nullptr);
		}
		else if (fp != nullptr) {
			mode |= LOGFILE | LOGEXTERNAL;
			lcp = m_dfLogConsole;
			lfp = fp;
		}
		m_selfp = std::unique_ptr<Log>(new Log(mode,lcp,lfp,dlgp));
		atexit(Log::Destroy);
	}
	return (m_selfp.get() != nullptr);
}

void Log::Destroy()
{
	if (m_selfp.get() != nullptr) {
		mfp_log_debug("Deleting Logger\n");
		m_selfp.reset();
	}
}

Log& Log::GetInstance()
{
	assert(m_selfp.get() != nullptr);
	return *m_selfp;
}

Log::Log(int mode,FILE* console,FILE* file,mfvWindow * dlgp)
: m_log()
, m_mode(mode)
, m_consolep(console)
, m_filep(file)
, m_dlgp(dlgp)
{
	//mfp_debug("Ctr\n");
	if (m_filep != nullptr) StartLogToFile();
}

Log::~Log()
{
	//mfp_debug("Dtr");
	Close();
	//mfp_debug("After Dtr\n");
}

void Log::Msg(const char* str)
{
	if (m_mode & LOGCONSOLE) {
		assert(m_consolep != nullptr);
		fprintf(m_consolep,str); fflush(m_consolep);
	}
	if (m_mode & LOGFILE) {
		m_log.flush();
		m_log << str;
	}
}

void Log::SetMode(int mode)
{
	m_mode = mode;
}

void Log::SetLevel(eLevel level)
{
	//_level = level;
}

void Log::Close()
{
	if (m_mode & LOGFILE && m_filep != nullptr) {
		//mfp_debug("In close log\n");
		EndLogToFile();
		fprintf(m_filep,"%s\n",m_log.str().c_str());
		fflush(m_filep);
		if (~(m_mode & LOGEXTERNAL)) fclose(m_filep);
	}
}

void Log::StartLogToFile()
{
	//_slogstrm.open(LogFileName,std::ios::in | std::ios::out);
	//if(_slogstrm.fail()) return;
	//mfp_debug("starting loggining\n");
	m_log.clear();
	m_log << "______Strat Log_______\n\n";
}

void Log::EndLogToFile()
{
	m_log << "______End Log________\n";
	//mfp_debug("ending loggining\n");
}

#ifdef PARALLEL
int Log::sendlog(const char *str) {
	return 1;
}

int Log::updatelog(const char *str) {
	return 1;
}
#endif

void log(const int level,const char* msg, ...)
{
	if (level > Log::Level()) return;
    char buffer[MAX_LOG_STRING];

    // do the formating
    va_list valist;
    va_start(valist, msg);
#ifdef _WIN32
	_vsnprintf(buffer,MAX_LOG_STRING,msg,valist);
#else
    vsprintf(buffer, msg, valist);
#endif
    va_end(valist);

    Log::GetInstance().Msg(buffer);
}

#ifdef PARALLEL
	void printParallel(const char* msg, ...)
	{
		// Sanity check
		//assert("pre: MPI controller is NULL!" && (comm != NULL) );
		//assert("pre: format argument is NULL!" && (format != NULL) )
		if( ::pcr_my_proc_id() == ::pcr_print_master())
		{
		    log(msg);
		    //printf("\nHello\n");
		    //fflush(stdout);
		}

		 ::pcr_barrier();
	}

	void printSynchronizedParallel(const char* msg, ...)
	{
		MPI_Request *request;
		// Sanity checks
		int procId  = ::pcr_my_proc_id();
		int numProc = ::pcr_nr_proc();
		char buffer[MAX_LOG_STRING];
		sprintf(buffer,"[%d]: %s",procId,msg);
		int* nullmsg = NULL;

		if (procId == ::pcr_print_master()) {
			// print message
		    log(buffer);

		    // signal next process (if any) to print
		    if (numProc > 1)
		    {
		    	MPI_Isend(nullmsg, 0, MPI_INT, 1, 0, MPI_COMM_WORLD,request);
		    }
		}
		else if (procId == numProc) {
		    // Block until previous process completes
		    ::pcr_receive_int(procId-1,0,0,nullmsg);

		    // print message
		    log(buffer);
		} // END last rank
		else {
		    // Block until previous process completes
		    ::pcr_receive_int(procId-1,0,0,nullmsg);

		    // print message
		    log(buffer);

		    // signal next process to print
		    MPI_Isend(nullmsg, 0, MPI_INT, procId, 0, MPI_COMM_WORLD,request);
		}

		::pcr_barrier();
	}
#endif

}// end namespace FemViewer
