/*
 * fv_helper_routines.c
 *
 *  Created on: 18 cze 2014
 *      Author: Paweł Macioł
 */
#include "fv_helper_routines.h"
#include <stdlib.h>

/* Fill in an array with values */
void fvFillArray(float *outArray, const int size, const float *value)
{
	int i;
	const float invf = 1.0f / RAND_MAX;

	for (i = 0; i < size; ++i) {
		outArray[i] = (! value) ? invf * rand() : *value;
	}
}




