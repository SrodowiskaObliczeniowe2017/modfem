/************************************************************************
File pds_conv_diff_std_util.c - approximation dependent utility procedures
                               for convection-diffusion equations

Contains routines:
  pdr_limit_slope - perform slope limiting
 ------------------------------  			
History:    
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* header files for the problem dependent module for Conv_Diff's equation */
#include "../../include/pdh_conv_diff.h"	

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

#include "pdh_intf.h"

/* interface of the mesh manipulation module */
#include "mmh_intf.h"	

/* interface for all approximation modules */
#include "aph_intf.h"	

/* utilities - including simple time measurement library */
#include "uth_intf.h"


/*---------------------------------------------------------
pdr_slope_limit - to limit the slope of linear solution
---------------------------------------------------------*/
int pdr_slope_limit( /* returns: >0 - success, <=0 - failure */
	int Problem_id	/* in: data structure to be used  */
	)
{

  return(1);

}

