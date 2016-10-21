/***************************************************
 * bath_odbo.h
 *
 * Header file for bath_mdfit.c.
 *
 * By Tseng Wei-Hsiang <dimsplendid@gmail.com>
 *
 *
 ***************************************************/

#ifndef _BATH_FITMD_H
#define _BATH_FITMD_H 1

#include "qdas.h"

/* Shared global variables */

/* exported functions */

int bath_mdfit_init_params(const size_t nsize, const double beta,
			   const size_t bath_nparams, const double *bath_params);
int bath_mdfit_free_params();
/* site lineshape function */
gsl_complex mdfit_G(double tau,double beta);
gsl_complex mdfit_H(double tau,double beta);
gsl_complex mdfit_C(double tau,double beta);
double mdfit_kernel_F(double tau, void * p);
double lambda0_f(double beta);
#endif /* bath_mdfit.h */

/*
 * $Log$
 * Revision 0.1  2015/11/05 08:33:15  dimsplendid
 *   - An average form from MD simulation correlation fuction fitting by 王佑仁.
 *
 */
