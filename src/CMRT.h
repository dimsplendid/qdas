/*****************************************************
 * [FILE]      CMRT.h
 * [AUTHOR]    Tseng Wei-Hsiang
 * [DATE]      20151015
 * ***************************************************/

#ifndef _CMRT_H
#define _CMRT_H

#include "gsl/gsl_matrix.h"
#include "qdas.h"
#include "params.h"
#include "aux.h"
#include "bath_odbo.h"
double char2double(const char *cline);
void plot_kernel(gsl_matrix * evec,gsl_vector * eval,double params[], double cut,int a,int b);
void cal_rate_matrix(gsl_matrix * evec,gsl_vector * eval,double params[], double cut,gsl_matrix * result);
void plot_C_r(void * params);
void plot_C_i(void * params);
#endif
