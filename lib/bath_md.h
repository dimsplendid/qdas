/***************************************************
 * bath_md.h
 *
 * Header file for bath_md.c.
 *
 * By Tseng Wei-Hsiang <dimsplendid@gmail.com>
 *
 *
 ***************************************************/

#pragma once
#include "qdas.h"
#include <complex.h>

/* Shared global variables */
/* exported functions */

int bath_md_init_params(const size_t nsize, const double beta,
                        const size_t bath_nparams, const double *bath_params);
int bath_md_free_params();
/* site lineshape function */
complex double md_C(double tau, complex double c[4096]);
complex double md_H(double tau, complex double h[4096]);
complex double md_G(double tau, complex double g[4096]);

double md_kernel_F(double tau, void *p);
double md_lambda0_f(complex double h[4096]);
void md_read(char *file_name, complex double array[4096]);

struct _md_par {
	double lambda0, Ea, Eb;
	double C_baab, C_aaba, C_bbba, C_aaab, C_bbab, C_bbaa, C_aaaa, C_bbbb;
	complex double c[4096];
	complex double h[4096];
	complex double g[4096];
};
typedef struct _md_par md_par;

/*
 * $Log$
 * Revision 0.1  2015/11/05 08:33:15 dimsplendid
 *   - The md fitting line shape function
 *
 */
