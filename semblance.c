/***************************************************************************
 *
 * Copyright 2015 Ian Liu Rodrigues and Edson Borin
 *
 ***************************************************************************/

#include "semblance.h"

#include <math.h>
#include <string.h>
#include <vector.h>
#include <utils.h>
#include <interpol.h>
#include <x86intrin.h>

static float time2d (float A, float B, float C, float t0,
                float m0x, float m0y,
                float mx, float my,
                float hx, float hy) __attribute__ ((const));

//XXX pure function

static float time2d(float A, float B, float C, float t0,
		float m0x, float m0y,
		float mx, float my,
		float hx, float hy)
{
	float _mx = mx - m0x;
	float _my = my - m0y;
	float _m2 = _mx*_mx + _my*_my;
	float t2 = t0;// + A * sqrt(_m2);
	t2 *= t2;
	//t2 += B * _m2;
	t2 += C * (hx*hx + hy*hy);
	//printf("%f\n",t2);
	if (t2 < 0)
		return -1;
	else
		return sqrt(t2);
}

static float get_halfoffset(su_trace_t *tr)
{
	float hx = (float)(tr->gx - tr->sx) / 2;
	float hy = (float)(tr->gy - tr->sy) / 2;
	return sqrt(hx * hx + hy * hy);
}

static float get_midpoint(su_trace_t *tr)
{
	float mx = (float)(tr->gx + tr->sx) / 2;
	float my = (float)(tr->gy + tr->sy) / 2;
	return sqrt(mx * mx + my * my);
}

static float get_scalco(su_trace_t *tr)
{
        if (tr->scalco == 0)
                return 1;
        if (tr->scalco > 0)
                return tr->scalco;
        return 1.0f / tr->scalco;
}


/* POINTER ALIASING, RESTRICT
float semblance_2d(aperture_t* __restrict ap,
		float A, float B, float C,
		int t0s, float m0x, float m0y,
		float* __restrict stack)*/
float semblance_2d(aperture_t *ap,
		float A, float B, float C,
		int t0s, float m0x, float m0y,
		float *stack)
{
	su_trace_t *tr = vector_get(ap->traces, 0);
	float dt = (float) tr->dt / 1000000;
	float idt = 1 / dt;
	float t0 = t0s * dt;
	int tau = MAX((int)(ap->ap_t * idt), 0);
	int w = 2 * tau + 1;
	float num[w], den[w];
	memset(&num[0], 0, sizeof(num));
	memset(&den[0], 0, sizeof(den));
	int M = 0, skip = 0;
	float _stack = 0;
	
	float t1 = t0*t0; // XXX loop invariant code motion
	float _mx,_my,_m2,t2[8],t3[8],t;
	
	su_trace_t *trp[8];
	int ii,i;
	for (i = 0; i < ap->traces.len-8; i+=8) {
		float mx, my, hx, hy;
		float a;
		
		for(ii=0;ii<8;ii++){
			trp[ii] = vector_get(ap->traces, i+ii);
			a = get_scalco(trp[ii])*0.5;
			
			mx=(trp[ii]->gx+trp[ii]->sx)*a;
			my=(trp[ii]->gy+trp[ii]->sy)*a;
			hx=(trp[ii]->gx-trp[ii]->sx)*a;
			hy=(trp[ii]->gy-trp[ii]->sy)*a;
			
			_mx = mx - m0x; // XXX inline
			_my = my - m0y;
			_m2 = _mx*_mx + _my*_my;
			t2[ii] = t1 + C * (hx*hx + hy*hy);
			t3[ii]=0;
		}
		//*
		__m256 raizes;
		raizes = _mm256_load_ps(t2);
		raizes = _mm256_sqrt_ps(raizes);
		_mm256_store_ps(t3,raizes);
		//*/
		/*
		float t_idt_tau=t*idt; // common subexpression elimination
		int it = (int)(t_idt_tau);
		t_idt_tau-=tau;
		//*/

		for(ii=0;ii<8;ii++){
			t3[ii]=(t2[ii] < 0)?-1:t3[ii];
			int it = (int)(t3[ii] * idt);
			if (it - tau >= 0 && it + tau < trp[ii]->ns) {
				for (int j = 0; j < w; j++) {
					int k = it + j - tau;
					//float v=(tr->data[k+1] - tr->data[k]) * (t_idt_tau+j - k) + tr->data[k];
					float v=(trp[ii]->data[k+1] - trp[ii]->data[k]) * (t3[ii]*idt+j-tau - k) + trp[ii]->data[k];
					num[j] += v;
					den[j] += v*v;
					_stack += v;
				}
				M++;
			} else if (++skip == 2) {
				goto error;
			}
		}
	}
	
	for (; i < ap->traces.len-4; i+=4) {
		float mx, my, hx, hy;
		float a;
		
		for(ii=0;ii<4;ii++){
			trp[ii] = vector_get(ap->traces, i+ii);
			a = get_scalco(trp[ii])*0.5;
			
			mx=(trp[ii]->gx+trp[ii]->sx)*a;
			my=(trp[ii]->gy+trp[ii]->sy)*a;
			hx=(trp[ii]->gx-trp[ii]->sx)*a;
			hy=(trp[ii]->gy-trp[ii]->sy)*a;
			
			_mx = mx - m0x; // XXX inline
			_my = my - m0y;
			_m2 = _mx*_mx + _my*_my;
			t2[ii] = t1 + C * (hx*hx + hy*hy);
			t3[ii]=0;
		}
		//*
		__m128 raizes;
		raizes = _mm_load_ps(t2);
		raizes = _mm_sqrt_ps(raizes);
		_mm_store_ps(t3,raizes);
		//*/
		/*
		float t_idt_tau=t*idt; // common subexpression elimination
		int it = (int)(t_idt_tau);
		t_idt_tau-=tau;
		//*/

		for(ii=0;ii<4;ii++){
			t3[ii]=(t2[ii] < 0)?-1:t3[ii];
			int it = (int)(t3[ii] * idt);
			if (it - tau >= 0 && it + tau < trp[ii]->ns) {
				for (int j = 0; j < w; j++) {
					int k = it + j - tau;
					//float v=(tr->data[k+1] - tr->data[k]) * (t_idt_tau+j - k) + tr->data[k];
					float v=(trp[ii]->data[k+1] - trp[ii]->data[k]) * (t3[ii]*idt+j-tau - k) + trp[ii]->data[k];
					num[j] += v;
					den[j] += v*v;
					_stack += v;
				}
				M++;
			} else if (++skip == 2) {
				goto error;
			}
		}
	}
	
	// peeling
	for (i; i < ap->traces.len; i++) {
		tr = vector_get(ap->traces, i);
		float mx, my, hx, hy;
		float a= get_scalco(tr)*0.5;
		
		mx=(tr->gx+tr->sx)*a;
		my=(tr->gy+tr->sy)*a;
		hx=(tr->gx-tr->sx)*a;
		hy=(tr->gy-tr->sy)*a;
		
		float t = time2d(A, B, C, t0, m0x, m0y, mx, my, hx, hy);

		int it = (int)(t * idt);
		if (it - tau >= 0 && it + tau < tr->ns) {
			for (int j = 0; j < w; j++) {
				int k = it + j - tau;
				float v=(tr->data[k+1] - tr->data[k]) * (t*idt+j-tau - k) + tr->data[k];
				num[j] += v;
				den[j] += v*v;
				_stack += v;
			}
			M++;
		} else if (++skip == 2) {
			goto error;
		}
	}

	float sem = 0;
	float aux = 0;
	for (int j = 0; j < w; j++) {
		sem += num[j] * num[j];
		aux += den[j];
	}

	if (stack) {
		_stack /= M*w;
		*stack = _stack;
	}

	return sem / aux / M;

error:
	return 0;
}
