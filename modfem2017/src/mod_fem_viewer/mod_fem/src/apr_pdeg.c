#include<stdio.h>
#include"mod_fem.h"
#include"../include/aph_intf.h"


int apr_get_pdeg_nrdofs(int base, int pdeg)
{
	int porder, porderz, numdofs=-1;
/* check */
#ifdef DEBUG
	if(pdeg<0 || (base == APC_COMPLETE && pdeg > 100)){
		printf("Wrong pdeg %d in apr_get_pdeg_nrdofs\n",pdeg);
		return( -1);	
	}
#endif
	
	if(base==APC_TENSOR_DG){

		/* decipher horizontal and vertical orders of approximation */
		porderz = pdeg/100;
		porder  = pdeg%100;
		numdofs = (porderz+1)*(porder+1)*(porder+2)/2;

	}
	else if(base==APC_COMPLETE_DG){
		porder  = pdeg;
		numdofs = (porder+1)*(porder+2)*(porder+3)/6;
	}
	else {
		printf("Type of base in get_pdeg_nrdofs not valid for prisms!\n");
	}

	return( numdofs);
}
