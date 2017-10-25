#ifndef _FV_TIMER_H_ 
#define _FV_TIMER_H_ 

#ifndef WIN32
#else
#define WINDOWS_LEAN_AND_MEAN
#include<windows.h>
#include<time.h>

class myTimer {
protected:
	LARGE_INTEGER m_iFrequecy;
	LARGE_INTEGER m_iLastQuery;
	LARGE_INTEGER m_iDelta;

	tm* m_pTime;
public:
	myTimer(void);
	~myTimer(void);

	void Upade(void);
	float getDelta(void);
	void print();

};

#endif

#endif /* _FV_TIMER_H_ 
*/
