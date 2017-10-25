#ifndef _FV_LOG_H_
#define _FV_LOG_H_

#ifdef PARALLEL
#include"pch_intf.h"
#endif
#include"uth_log.h"
#include<memory>
#include<cstdio>
#include<fstream>
#include<sstream>

namespace FemViewer {

	// Destination of messages
	enum {
		NOWERE	   = 0x0000,  // nowhere
		LOGCONSOLE = 0x0001,  // direct messages to stdout/stderr
		LOGFILE    = 0x0010,  // direct messages to internal log file (at destructur write to file)
		LOGDIALOG  = 0x0100,  // direct messages to dialog window
		LOGEXTERNAL = 0x1000, // whether log messages are written to external or internal logfile
	};

	// Levels of messages
	enum eLevel {
		LogERROR = 0,		// log errors
		LogWARNING,			// log errors + warnings
		LogINFO,			// log error + warnings + infos
		LogDEBUG,			// + debug info 0
		LogDEBUG1,			// + debug info 1, etc
	};

	// Forward declaration
	class mfvWindow;

	class Log
	{
		static std::unique_ptr<Log> m_selfp;
		static const eLevel m_dfErrLevel;	// Default error level logginig
		static const int    m_dfLogDest;		// Default destination for log messages
		static FILE * m_dfLogFile;			// Handle to default logfile
		static FILE * m_dfLogConsole;	// Default log destinaltion console (stdout or stderr)

	public:
		static bool Init(FILE *fp = stdout,mfvWindow *dlgp = nullptr);
		static bool IsInit(void) { return m_selfp != nullptr; }
		static const eLevel& Level() { return m_dfErrLevel; }
		static void Destroy();
		static Log& GetInstance();
		~Log();
		void Msg(const char* str);
		void SetMode(int mode);
		void SetLevel(eLevel level);
		void Close();

	protected:
		std::ostringstream m_log;
		int		  	  m_mode;
		FILE*		  m_consolep;
		FILE*		  m_filep;
		mfvWindow*    m_dlgp;

		void StartLogToFile();
		void EndLogToFile();

	private:
		// Ctr
		Log(int mode,FILE* console,FILE* file,mfvWindow * dlgp);

		// These are not allowed
		Log(const Log& rhs);
		Log& operator=(const Log& rhs);

#ifdef PARALLEL
		int 	  _myprocid;
		int 	  _mymasterid;
		int (*logto)(const char* str);
		static int sendlog(const char* str);
		static int updatelog(const char* str);
#endif
	};

void log(const int level, const char *msg, ...);
// this macro is for showing only the filename, not full file path
#define MFP_FILE_NAME 	MF_FILE_NAME
#ifdef NDEBUG
//#define mfp_log_err(M, ...)
//#define mfp_log_warn(M, ...)
//#define mfp_log_info(M, ...)
#define mfp_log_debug(M, ...)
#else
//#define mfp_log_err(M, ...)		log(FemViewer::LogERROR, "ERROR: " M "\n", ##__VA_ARGS__)
//#define mfp_log_warn(M, ...)	log(FemViewer::LogWARNING,"WARNING: " M "\n", ##__VA_ARGS__)
//#define mfp_log_info(M, ...)	log(FemViewer::LogINFO,"INFO: " M "\n", ##__VA_ARGS__)
#define mfp_log_debug(M, ...)	log(FemViewer::LogDEBUG, "DEBUG %s:%d: " M "\n",MFP_FILE_NAME, __LINE__, ##__VA_ARGS__)
#endif

#define mfp_log_err(M, ...)		log(FemViewer::LogERROR, "ERROR: " M "\n", ##__VA_ARGS__)
#define mfp_log_warn(M, ...)	log(FemViewer::LogWARNING,"WARNING: " M "\n", ##__VA_ARGS__)
#define mfp_log_info(M, ...)	log(FemViewer::LogINFO,"INFO: " M "\n", ##__VA_ARGS__)

#ifdef PARALLEL
	void printParallel(const char*, ...);
	void printSynchronizedParallel(const char*, ...);
#endif


}// end namespace FeemViewer


#endif
