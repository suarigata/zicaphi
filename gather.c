/***************************************************************************
 *
 * Copyright 2015 Ian Liu Rodrigues and Edson Borin
 *
 ***************************************************************************/

#include "gather.h"

#include <utils.h>

void gather_get_traces_for_aperture(int *m0s,
		su_trace_t *traces, int ntr,
		aperture_t *ap)
{
	vector_init(ap->traces);
	float x0, y0;
	float x, y;
	su_get_receiver(&traces[*m0s], &x0, &y0);
	su_get_receiver(traces, &x, &y);
	float dist = SQR(x0 - x) + SQR(y0 - y);
	while (ntr-- && dist > SQR(ap->ap_m)) {
		traces++;
		su_get_receiver(traces, &x, &y);
		dist = SQR(x0 - x) + SQR(y0 - y);
		(*m0s)--;
	}
	ntr++;
	if (ntr == 0)
		return;
	while (ntr-- && dist <= SQR(ap->ap_m)) {
		vector_push(ap->traces, traces);
		traces++;
		su_get_receiver(traces, &x, &y);
		dist = SQR(x0 - x) + SQR(y0 - y);
	}
}

