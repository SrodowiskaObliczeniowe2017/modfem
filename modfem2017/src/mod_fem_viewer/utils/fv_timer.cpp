/*
 * fv_timer.cpp
 *
 *  Created on: 23-10-2013
 *      Author: Paweł Macioł
 */

#include "fv_timer.h"


fv_timer::fv_timer() : stime_us(0), etime_us(0), running(1)
{
#ifdef WIN32
    QueryPerformanceFrequency(&freq);
    scount.QuadPart = 0;
    ecount.QuadPart = 0;
    lcount.QuadPart = 0;
    milion_over_freq = 1000000.0 / (double)freq.QuadPart;
#else
    scount.tv_sec = scount.tv_usec = 0;
    ecount.tv_sec = ecount.tv_usec = 0;
    lcount.tv_sec = ecount.tv_usec = 0;
#endif
}

fv_timer::~fv_timer()
{
}


double fv_timer::get_time_us()
{
#ifdef WIN32
    if(running)
        QueryPerformanceCounter(&ecount);

    stime_us = scount.QuadPart * milion_over_freq;
    etime_us = ecount.QuadPart * milion_over_freq;
#else
    if(running)
        gettimeofday(&ecount, NULL);

    stime_us = (scount.tv_sec * 1000000.0) + scount.tv_usec;
    etime_us = (ecount.tv_sec * 1000000.0) + ecount.tv_usec;
#endif

    return etime_us - stime_us;
}

double fv_timer::get_duration_us()
{
#ifdef WIN32
    if(running)
        QueryPerformanceCounter(&ecount);

    ltime_us = lcount.QuadPart * milion_over_freq;
    etime_us = ecount.QuadPart * milion_over_freq;
#else
    if(running)
        gettimeofday(&ecount, NULL);

    ltime_us = (lcount.tv_sec * 1000000.0) + lcount.tv_usec;
    etime_us = (ecount.tv_sec * 1000000.0) + ecount.tv_usec;
#endif
    lcount = ecount;
    return etime_us - ltime_us;
}
