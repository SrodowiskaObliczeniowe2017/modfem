/***
Pawe� Macio�
Time mesasurment - windows version
***/

#ifndef _pomiar_czasu_H_
#define _pomiar_czasu_H_

#include "fv_compiler.h"

#ifdef __cplusplus
extern "C" 
{
#endif

/*---------------------------------------------------------
  czas_C - zwraca czas CPU w sekundach od rozpoczecia dzialania lub -1 jesli
           funkcja clock() nie jest prawidlowo zaimplementowana w systemie 
---------------------------------------------------------*/
double czas_C();

/*---------------------------------------------------------
  inicjuj_czas - inicjuje pomiar czasu 
---------------------------------------------------------*/
void inicjuj_czas();

/*---------------------------------------------------------
  czas_zegara - procedura zwraca czas zegarowy w sekundach od momentu
                inicjacji (jako liczbe podwojnej precyzji) 
---------------------------------------------------------*/
double czas_zegara();

/*---------------------------------------------------------
  czas_CPU - procedura zwraca czas CPU w sekundach od momentu
             inicjacji (jako liczbe podwojnej precyzji) 
 ---------------------------------------------------------*/
double czas_CPU();

/*---------------------------------------------------------
  drukuj_czas - pomiar i wydruk czasu standardowego C, CPU i zegarowego 
               w sekundach od momentu inicjacji pomiaru czasu
 ---------------------------------------------------------*/
void drukuj_czas();

#ifdef __cplusplus 
}
#endif

#endif /* _pomiar_czasu_H_
*/
