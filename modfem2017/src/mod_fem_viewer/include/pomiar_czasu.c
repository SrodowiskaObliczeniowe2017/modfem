#include"pomiar_czasu.h"

// system
#include<stdio.h>
#include<time.h>
#include<sys/timeb.h>

// pomocnicze zmienne 
#ifdef WIN32
static clock_t ct;
struct _timeb tb;
static LARGE_INTEGER iFrequency;
static LARGE_INTEGER iLastQuery;

static int iIsFrequencyInit;


/*---------------------------------------------------------
  czas_C - zwraca czas CPU w sekundach od rozpoczecia dzialania lub -1 jesli
           funkcja clock() nie jest prawidlowo zaimplementowana w systemie 
---------------------------------------------------------*/
double czas_C()
{        
	clock_t timer;

	timer = clock();

	if(timer==-1) return( (double) timer);
	else return ( (double) timer / (double) CLOCKS_PER_SEC);
}

/*---------------------------------------------------------
inicjuj_czas - inicjuje pomiar czasu nadajac wartosci 
               poczatkowe odpowiednim zmiennym
---------------------------------------------------------*/
void inicjuj_czas()
{
	if(!iIsFrequencyInit)
	{
		QueryPerformanceFrequency( &iFrequency);
		iIsFrequencyInit++;
	}

	QueryPerformanceCounter( &iLastQuery);
		
	ct = clock();
	_ftime_s( &tb);
}

/*---------------------------------------------------------
  czas_zegara - procedura zwraca czas zegarowy w sekundach od momentu
                inicjacji (jako liczbe podwojnej precyzji) 
---------------------------------------------------------*/
double czas_zegara()
{
	struct _timeb ltb;
	time_t lt;

	_ftime_s( &ltb);
	lt  = (ltb.time - tb.time)*1000;
	lt += ltb.millitm - tb.millitm;

	return( (double)lt / 1.e3); // rozdielczo�� do milisekund
}

/*---------------------------------------------------------
  czas_CPU - procedura zwraca czas CPU w sekundach od momentu
             inicjacji (jako liczbe podwojnej precyzji) 
 ---------------------------------------------------------*/
double czas_CPU()
{
	LARGE_INTEGER iCurrentQuery;
	double		  cputime;

	QueryPerformanceCounter( &iCurrentQuery);
	cputime = (double)(iCurrentQuery.QuadPart - iLastQuery.QuadPart) / (double)iFrequency.QuadPart;
  
	return( cputime);
}

/*---------------------------------------------------------
 drukuj_czas - pomiar i wydruk czasu CPU i zegarowego 
               od momentu inicjacji pomiaru czasu
 ---------------------------------------------------------*/
void drukuj_czas()
{ 
	clock_t	timer;
	LARGE_INTEGER iCurrentQuery;
	struct _timeb ltb;
	time_t lt;
	double stdtime, cputime, daytime;
  
	timer = clock();
	QueryPerformanceCounter( &iCurrentQuery);
	_ftime_s( &ltb);

	lt  = (ltb.time - tb.time)*1000;
	lt += ltb.millitm - tb.millitm;

	daytime = (double)lt / 1.e3;
	cputime = (double)(iCurrentQuery.QuadPart - iLastQuery.QuadPart) / (double)(iFrequency.QuadPart);
	stdtime = (double)(timer - ct) / (double) CLOCKS_PER_SEC;
	  
	printf("czas standardowy = %.12lf [s]\n",stdtime);
	printf("czas CPU         = %.12lf [s]\n",cputime);
	printf("czas zegarowy    = %.12lf [s]\n",daytime);
	printf("rozdzielczosc zegara: %d\n",iFrequency);
  
}
#else
#endif
