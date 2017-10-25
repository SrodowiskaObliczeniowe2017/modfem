/*
 * Consolidator.cpp
 *
 *  Created on: 2011-04-19
 *      Author: Paweł Macioł
 */

#include "Consolidator.h"

Consolidator* pConsolidator = NULL;

Consolidator& ConsolInstance()
{
	if(!pConsolidator) pConsolidator = new Consolidator();
	return *pConsolidator;
}

