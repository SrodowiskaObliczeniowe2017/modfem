#include "../inc/timer.h"
#include <cstdio>

#ifdef WIN32

myTimer::myTimer(void)
{
	QueryPerformanceFrequency(&m_iFrequecy);
	QueryPerformanceCounter(&m_iLastQuery);
}

myTimer::~myTimer(void)
{
	m_iFrequecy.QuadPart  = 0;
	m_iLastQuery.QuadPart = 0;
}

void myTimer::Upade(void)
{
	LARGE_INTEGER	tmpTimer;
	time_t			iTempTimeDate;

	QueryPerformanceCounter(&tmpTimer);
	m_iDelta.QuadPart = tmpTimer.QuadPart - m_iLastQuery.QuadPart;

	/* Keep current clock state */
	m_iLastQuery.QuadPart = tmpTimer.QuadPart;

	/* Get current date and time */
	time(&iTempTimeDate);
	m_pTime = localtime(&iTempTimeDate);

}

float myTimer::getDelta(void)
{
	return((static_cast<float>(m_iDelta.QuadPart))/(static_cast<float>(m_iFrequecy.QuadPart)));
}


void myTimer::print()
{
	printf("Measured time in ms: %lf",getDelta());
}

#endif
