/***************************************************************************
 *
 * Copyright 2015 Ian Liu Rodrigues and Edson Borin
 *
 ***************************************************************************/

#include "interpol.h"

float interpol_linear(float x0, float x1, float y0, float y1, float x)
{
	return (y1 - y0) * (x - x0) / (x1 - x0) + y0;
}


//x0, x1 , y0         , y1           , x
//k , k+1, tr->data[k], tr->data[k+1], t*idt+j-tau

//(tr->data[k+1] - tr->data[k]) * (t*idt+j-tau - k) + tr->data[k]
