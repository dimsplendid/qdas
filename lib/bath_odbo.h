/***************************************************
 * bath_odbo.h
 *
 * Header file for bath_odbo.c.
 *
 * By Tseng Wei-Hsiang <dimsplendid@gmail.com>
 *
 *
 ***************************************************/

#ifndef _BATH_ODBO_H
#define _BATH_ODBO_H 1

#include "qdas.h"

/* Shared global variables */
extern double BathODBOLambda;
extern double BathODBOGAMMA;

/* exported functions */

int bath_odbo_init_params(const size_t nsize, const double beta,
			   const size_t bath_nparams, const double *bath_params);
int bath_odbo_free_params();
/* site lineshape function */
double g_r(double tau,void * params);
double g_i(double tau,void * params);
double gg_r(double tau,void * params);
double gg_i(double tau,void * params);
double ggg_r(double tau,void * params);
double ggg_i(double tau,void * params);
/* exciton basis lineshape function                            */
/* Because these lineshape functions and reorginization energy */
/* do not depend on the site index, so the transfer sum can    */
/* just iosolate and be a product form.                        */
/* In this way, I just use a transform sum function.           */
double transCoeff(gsl_matrix * eigenvectors ,int ,int, int, int);
double kernel_F(double tau, void * p);

#endif /* bath_odbo.h */

/*
 * $Log$
 * Revision 0.1  2015/11/05 08:33:15 dimsplendid
 *   - A simple over-damped Brownian Oscillators Bath
 *
 */
