/***************************************************************************
 *
 * Copyright 2015 Ian Liu Rodrigues and Edson Borin
 *
 ***************************************************************************/

#ifndef SEMBLANCE_H__
#define SEMBLANCE_H__

#include <vector.h>
#include <su.h>

typedef struct aperture aperture_t;

struct aperture {
	float ap_m, ap_h, ap_t;
	vector_t(su_trace_t*) traces;
};

float semblance_2d(aperture_t *ap,
		float A, float B, float C,
		int t0s, float m0x, float m0y,
		float *stack);

#endif /* SEMBLANCE_H__ */
