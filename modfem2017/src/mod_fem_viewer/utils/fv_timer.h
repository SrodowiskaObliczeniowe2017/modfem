/*
 * fv_timer.h
 *
 *  Created on: 23-10-2013
 *      Author: Paweł Macioł
 */

#ifndef FV_TIMER_H_
#define FV_TIMER_H_

#ifdef WIN32
#include <windows.h>
#include <stdlib.h>
#else
#include <sys/time.h>
#include <stdlib.h>
#endif

class fv_timer {
  public:
	fv_timer();
	~fv_timer();

	inline void start();
	inline void stemp();
	inline void stop();

	double get_time_sec()	{ return get_time_us() * 1e-6; }
	double get_time_ms()	{ return get_time_us() * 1e-3; }
	double get_time_us();

	double get_duration_sec() { return get_duration_us() * 1e-6; }
	double get_duration_ms()  { return get_duration_us() * 1e-3; }
	double get_duration_us();

  private:
	double stime_us;
	double etime_us;
	double ltime_us;
	int    running;

#ifdef WIN32
    LARGE_INTEGER freq;                    // ticks per second
    LARGE_INTEGER scount;                  //
    LARGE_INTEGER ecount;                  //
    LARGE_INTEGER lcount;
    double milion_over_freq;
#else
    timeval scount;                         //
    timeval ecount;                           //
    timeval lcount;
#endif
};

void fv_timer::start()
{
#ifdef WIN32
	QueryPerformanceCounter(&_start);
#else
    gettimeofday(&scount, NULL);
#endif
    lcount = scount;
}

void fv_timer::stemp()
{
#ifdef WIN32
    QueryPerformanceCounter(&lcount);
#else
    gettimeofday(&lcount, NULL);
#endif
}

void fv_timer::stop()
{
#ifdef WIN32
    QueryPerformanceCounter(&ecount);
#else
    gettimeofday(&ecount, NULL);
#endif
    lcount = ecount;
    running = 0;
}


#endif /* FV_TIMER_H_ */
