#include"fv_sysloc.h"
#include<stdio.h>
//#include<crtdefs.h>
#include<locale.h>

void check_locale()
{
	_locale_t currloc = _get_current_locale();
	#ifdef _DEBUG
	puts("Checking decimal point...\n");
	#endif

	lconv * sysloc = localeconv();
	if(currloc->locinfo->lconv->decimal_point == (char*)(',')) {
	//if(sysloc->decimal_point == (char*)(',')) {
		#ifdef _DEBUG
		puts(">>>>>>>>>>>>>>>>>>>Warrning: Decimal point changing to '.'!\n");
		#endif
		setlocale(LC_NUMERIC, "C");
	}
}
