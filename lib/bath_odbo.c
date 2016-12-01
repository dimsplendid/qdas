/***************************************************
 * bath_odbo.c
 *
 * Part of the qdas package;
 * Over-damped Brownian Oscillators bath module.
 *
 * By Tseng Wei-Hsiang <dimsplendid@gmail.com>
 *
 ***************************************************/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex_math.h>

#include "qdas.h"
#include "aux.h"
#include "bath_odbo.h"

/* Global constant */
/* Integration constants */
#define NWSPACE (100000)
#define EPSABS (1e-5)
#define EPSREL (1e-5)
#define PI 3.14159265359
#define MAZCUT 1 // the sum of Mazbara term from 1 to n

// convert time from cm unit to fs
#define CM2FS (5309.1)

/* Define gsl complex macro */
#define ADD gsl_complex_add
#define MUL gsl_complex_mul
#define SUB gsl_complex_sub
#define INV gsl_complex_inverse
#define ADDR gsl_complex_add_real
#define SUBR gsl_complex_sub_real
#define MULR gsl_complex_mul_real
#define DIVR gsl_complex_div_real
#define EXP gsl_complex_exp
#define TAN gsl_complex_tan

double BathODBOLambda;
double BathODBOGAMMA;

/* lineshape functions */
double J(double omega, void * params){
  /* Over-damped Brownian oscillators */
  /* in this condition, I don't use this */
  double result;
  double lambda0, Gamma0;
  double *p = (double *) params;
  lambda0 = p[0];
  Gamma0 = p[1];
  result = (2.0*lambda0*omega*Gamma0)/(PI*omega*omega+Gamma0*Gamma0);
  return result;
}
/* I would use lambda_0 directly */
// double lambda_site(double tau,void *params){}

double transCoeff(gsl_matrix * evec ,int a,int b,int c,int d){
  int i;
  int size = evec->size1;
  double coeff,result ;
  double a_i,b_i,c_j,d_j;
  // printf("C_ai C_bi C_cj C_dj\n");
  // sum as i = j
  result = 0.0;
  for(i = 0;i<size;i++){
    a_i = gsl_matrix_get(evec,i,a);
    b_i = gsl_matrix_get(evec,i,b);
    c_j = gsl_matrix_get(evec,i,c);
    d_j = gsl_matrix_get(evec,i,d);
    // printf("%.2f %.2f %.2f %.2f\n",a_i,b_i,c_j,d_j);
    coeff = a_i * b_i * c_j * d_j;
    result += coeff;
  }
  // printf("C_%d%d%d%d = %.6f\n",a,b,c,d,result);
  return result;
}

double g_r(double tau,void * params){
  double result = 0.0;
  double beta,lambda0,Gamma0,cotx;
  double *p = (double *) params;
  beta = p[0];
  lambda0 = p[1];
  Gamma0 = p[2];
  cotx = p[3];

  // Mazbara term sum 1..MAZCUT
  // printf("test exp term: %.18f \n",beta);

  for(int i = 1; i < MAZCUT; i++){
    result += (1.0/i/PI)*\
      ((exp(-2.0*i*PI/beta)*tau-2.0*i*PI/beta*tau+1.0))/\
      (4*i*i*PI*PI/beta/beta+Gamma0*Gamma0);
  }
  // printf("The Mazbara term: %.18f\n",result );
  result = result*(-2.0*lambda0*Gamma0) + lambda0/Gamma0*cotx*\
    (exp(-Gamma0*tau)+Gamma0*tau-1.0); // unitless
  return result;
}
double g_i(double tau,void * params){
  double result = 0.0;
  // double beta;
  double lambda0,Gamma0;
  double *p = (double *) params;
  // beta = p[0]; // not used
  lambda0 = p[1];
  Gamma0 = p[2];

  result = -lambda0/Gamma0*(exp(-Gamma0*tau)+Gamma0*tau - 1.0); // unitless
  return result;
}
double gg_r(double tau,void * params){
  double result = 0.0;
  double beta,lambda0,Gamma0,cotx;
  double *p = (double *) params;
  beta = p[0];
  lambda0 = p[1];
  Gamma0 = p[2];
  cotx = p[3];
  // printf("params: %.3f %.3f %.3f %.3f\n",beta,lambda0,Gamma0,cotx);
  // printf("exp term: %.18f\n",-2*i*PI/beta*tau+1);
  // Mazbara term sum 1..MAZCUT
  for(int i = 1; i < MAZCUT; i++){
    result += (exp(-2.0*i*PI/beta*tau)+1.0)/\
      (4.0*i*i*PI*PI/beta/beta+Gamma0*Gamma0);
  }
  // printf("result: %.3f",result);
  result = result*(4.0*lambda0*Gamma0/beta) + lambda0*cotx*\
    (-exp(-Gamma0*tau)+1.0); // unit: cm^-1
  return result;
}
double gg_i(double tau,void * params){
  double result = 0.0;
  // double beta; // not used
  double lambda0,Gamma0;
  double *p = (double *) params;
  // beta = p[0]; // not used
  lambda0 = p[1];
  Gamma0 = p[2];

  result = lambda0*(exp(-Gamma0*tau)-1); // unit: cm^-1
  return result;
}
double ggg_r(double tau,void * params){
  double result = 0.0;
  double beta,lambda0,Gamma0,cotx;
  double *p = (double *) params;
  beta = p[0];
  lambda0 = p[1];
  Gamma0 = p[2];
  cotx = p[3];

  // Mazbara term sum 1..MAZCUT
  for(int i = 1; i < MAZCUT; i++){
    result += (4*i*PI/beta/beta)*\
      exp(-2*i*PI/beta*tau)/\
      (4*i*i*PI*PI/beta/beta+Gamma0*Gamma0);
  }
  result = result*(-2.0*lambda0*Gamma0) + lambda0*Gamma0*cotx*\
    exp(-Gamma0*tau); // unit cm^-2
  return result;
}
double ggg_i(double tau, void * params){
  double result = 0.0;
  // double beta; // not used
  double lambda0,Gamma0;
  double *p = (double *) params;
  // beta = p[0]; // not used
  lambda0 = p[1];
  Gamma0 = p[2];

  result = -lambda0*Gamma0*exp(-Gamma0*tau);// unit: cm^-2
  return result;
}
/*
* Use the bi-exponetial fitting for the lineshape function
* Tian, B. L.; Ding, J. J.; Xu, R. X.; Yan, Y.
* Biexponential Theory of Drude Dissipation via Hierarchical Quantum Master Equation.
* J. Chem. Phys. 2010, 133 (11), 1–6.
* C(t>0) ~ c1 exp(-r1 t) + c2 exp(-r2 t) + 2Δδ(t)
*/

gsl_complex C_biexp(double tau,void *p){
  gsl_complex result;
  double *params = (double *)p;

  // input params
  double beta   = params[0];
  double lambda = params[1];
  double r1     = params[2];

  // biexp fitting params
  double r2 = sqrt(42.0)/beta;
  double Delta = beta*lambda*r1/20.0;
  gsl_complex c1,c2;
  GSL_SET_COMPLEX(&c1,(2.0*lambda/beta)*(1.0-10*(Delta/lambda)*(Delta/lambda)\
                                      -2.45/((r2/r1)*(r2/r1)-1.0))\
                                      ,-lambda*r1);
  GSL_SET_COMPLEX(&c2,98*r2*Delta/(beta*beta)/(r2*r2-r1*r1),0.0);
  // C = c1*exp(-r1*t)+c2*exp(-r2*t);
  result = ADD(MULR(c1,exp(-r1*tau)),MULR(c2,exp(-r2*tau)));
  return result;
}

gsl_complex H_biexp(double tau,void *p){
  gsl_complex result;
  double *params = (double *)p;

  // input params
  double beta   = params[0];
  double lambda = params[1];
  double r1     = params[2];

  // biexp fitting params
  double r2 = sqrt(42.0)/beta;
  double Delta = beta*lambda*r1/20.0;
  gsl_complex c1,c2;
  GSL_SET_COMPLEX(&c1,(2.0*lambda/beta)*(1.0-10*(Delta/lambda)*(Delta/lambda)\
                                      -2.45/((r2/r1)*(r2/r1)-1.0))\
                                      ,-lambda*r1);
  GSL_SET_COMPLEX(&c2,98*r2*Delta/(beta*beta)/(r2*r2-r1*r1),0.0);
  // H = (c1-c1*exp(-r1*t))/r1+(c2-c2*exp(-r2*t))/r2+Delta;
  result = ADDR(ADD(DIVR(SUB(c1,MULR(c1,exp(-r1*tau))),r1),\
                    DIVR(SUB(c2,MULR(c2,exp(-r2*tau))),r2)),Delta);
  return result;
}

gsl_complex G_biexp(double tau,void *p){
  gsl_complex result;
  double *params = (double *)p;

  // input params
  double beta   = params[0];
  double lambda = params[1];
  double r1     = params[2];

  // biexp fitting params
  double r2 = sqrt(42.0)/beta;
  double Delta = beta*lambda*r1/20.0;
  gsl_complex c1,c2;
  GSL_SET_COMPLEX(&c1,(2.0*lambda/beta)*(1.0-10*(Delta/lambda)*(Delta/lambda)\
                                      -2.45/((r2/r1)*(r2/r1)-1.0))\
                                      ,-lambda*r1);
  GSL_SET_COMPLEX(&c2,98*r2*Delta/(beta*beta)/(r2*r2-r1*r1),0.0);
  // G = c1*(-1+exp(-r1*t))/r1^2+c2*(-1+exp(-r2*t))/r2^2+c1*t/r1+c2*t/r2+Delta*t;
  result = ADDR(ADD(\
                ADD(DIVR(MULR(c1,exp(-r1*tau)-1.0),r1*r1),MULR(c1,tau/r1)),\
                ADD(DIVR(MULR(c2,exp(-r2*tau)-1.0),r2*r2),MULR(c2,tau/r2))),\
                Delta*tau);
  return result;
}

double kernel_F(double tau, void * p){
  // exp()[Re{}cos()-Im{}sin()]
  // double imag_part;
  double real_part;
  double * params = (double*)p;
  double lambda0 = params[1];
  double Ea = params[4];
  double Eb = params[5];
  double C_baab = params[6];
  double C_aaba = params[7];
  double C_bbba = params[8];
  double C_aaab = params[9];
  double C_bbab = params[10];
  double C_bbaa = params[11];
  double C_aaaa = params[12];
  double C_bbbb = params[13];
  /*
  result = exp((2*C_bbaa-C_aaaa-C_bbbb)*g_r(tau,params))\
    *((((C_baab*ggg_r(tau,params)-(C_aaba-C_bbba)*(C_aaab-C_bbab)*gg_r(tau,params)*gg_r(tau,params))\
    +((C_aaba-C_bbba)*gg_i(tau,params)-2.0*C_bbba*lambda0)*((C_aaab-C_bbab)*gg_i(tau,params)-2.0*C_bbab*lambda0))\
    *cos((2.0*C_bbaa-C_aaaa-C_bbbb)*g_i(tau,params)+(2.0*(C_bbaa-C_bbbb)*lambda0+(Eb-Ea))*tau))\
    +((-C_baab*ggg_i(tau,params)+(C_aaba-C_bbba)*gg_r(tau,params)*((C_aaab-C_bbab)*gg_i(tau,params)-2.0*C_bbab*lambda0)\
    +(C_aaab-C_bbab)*gg_r(tau,params)*((C_aaba-C_bbba)*gg_i(tau,params)-2.0*C_bbab*lambda0))\
    *sin((2.0*C_bbaa-C_aaaa-C_bbbb)*g_i(tau,params)+(2.0*(C_bbaa-C_bbbb)*lambda0+(Eb-Ea))*tau)));
  */
  double G_r,G_i,H_r,H_i,C_r,C_i,cosx,sinx,Re_p,Im_p,expx;
  gsl_complex G,H,C;
  /*
  * ODBO analytic need to modified
  G_r = g_r(tau,params);
  G_i = g_i(tau,params);
  H_r = gg_r(tau,params);
  H_i = gg_i(tau,params);
  C_r = ggg_r(tau,params);
  C_i = ggg_i(tau,params);
  */
  // ODBO biexp fitting
  G = G_biexp(tau,params);
  H = H_biexp(tau,params);
  C = C_biexp(tau,params);
  G_r = GSL_REAL(G);
  G_i = GSL_IMAG(G);
  H_r = GSL_REAL(H);
  H_i = GSL_IMAG(H);
  C_r = GSL_REAL(C);
  C_i = GSL_IMAG(C);
  expx = exp((2.0*C_bbaa - C_aaaa - C_bbbb)*G_r);
  cosx = cos((2.0*C_bbaa - C_aaaa - C_bbbb)*G_i + (2.0*(C_bbaa - C_bbbb)*lambda0+Eb-Ea)*tau);
  sinx = sin((2.0*C_bbaa - C_aaaa - C_bbbb)*G_i + (2.0*(C_bbaa - C_bbbb)*lambda0+Eb-Ea)*tau);

  Re_p = C_baab*C_r - ((C_aaba - C_bbba)*H_r*(C_aaab - C_bbab)*H_r \
        -(((C_aaba-C_bbba)*H_i - 2.0*C_bbba*lambda0)*((C_aaab-C_bbab)*H_i - 2.0*C_bbab*lambda0)));
  Im_p = C_baab*C_i - ((C_aaba - C_bbba)*H_r*((C_aaab-C_bbab)*H_i - 2.0*C_bbab*lambda0)\
        + ((C_aaba-C_bbba)*H_i - 2.0*C_bbba*lambda0)*(C_aaab - C_bbab)*H_r);
  real_part = expx * (Re_p * cosx - Im_p * sinx);


  // imag_part = expx * (Re_p * sinx + Im_p * cosx);
  return real_part; // unitless
}

/* interface functions */
int bath_odbo_init_params(const size_t nsize, const double beta,
			   const size_t bath_nparams, const double *bath_params)
{
  int idx;
  idx = 0;
  printf("Use Over-damped Brownian Oscillators bath:\n");
  printf("\n");
  printf("          /inf \n");
  printf("g_nm(t) = | dw (J_nm/w^2)[(coth(beta*w/2)(1-cos(wt))+i(sin(wt)-wt)]\n");
  printf("          /0 \n");
  printf("\n");
  printf("J_nm = 2*lambda0/pi*w*Gamma0/(w^2+Gamma0^2)\n");
  printf("\n");
  printf("For each site:\n");
  printf("%12s %12s\n","lambda (cm^-1)","tau (cm^-1)");
  BathODBOLambda = bath_params[idx];
  BathODBOGAMMA = bath_params[idx+1];
  printf("%12.4f %12.4f\n",BathODBOLambda,BathODBOGAMMA);
  printf("\n");

  return 0;
}

int bath_odbo_free_params(){
  return 0;
}
#ifdef MAIN
void plot_C(void * params){
	double t;
  // double tmp_Ci;
  gsl_complex c;
	for(uint32_t i = 0; i < 3767; i++){
		t = ((double)i)/10000.0;
    /*
    tmp_Ci = ggg_i(t,params);
    if (tmp_Ci < 0.0)
		  printf("%.18f,%.18f%.18fi\n",t*CM2FS,ggg_r(t,params),tmp_Ci);
    else
      printf("%.18f,%.18f+%.18fi\n",t*CM2FS,ggg_r(t,params),tmp_Ci);
    */
    c = C_biexp(t,params);
    printf("%.18f, %.18f %+.18fi\n",t*CM2FS,GSL_REAL(c),GSL_IMAG(c));
	}
}
void plot_H(void * params){
  double t;
  gsl_complex h;
  // double tmp_Hi;
	for(uint32_t i = 0; i < 3767; i++){
		t = ((double)i)/10000.0;
    /*
    tmp_Hi = gg_i(t,params);
    if (tmp_Hi < 0.0)
      printf("%.18f,%.18f%.18fi\n",t*CM2FS,gg_r(t,params),tmp_Hi);
    else
		  printf("%.18f,%.18f+%.18fi\n",t*CM2FS,gg_r(t,params),tmp_Hi);
    */
    h = H_biexp(t,params);
    printf("%.18f, %.18f %+.18fi\n",t*CM2FS,GSL_REAL(h),GSL_IMAG(h));
	}
}
void plot_G(void * params){
  double t;
  gsl_complex g;
  // double tmp_Gi;
	for(uint32_t i = 0; i < 3767; i++){
		t = ((double)i)/10000.0;
    /*
    tmp_Gi = g_i(t,params);
    if (tmp_Gi < 0.0)
      printf("%.18f\t%.18f%.18fi\n",t*CM2FS,g_r(t,params),tmp_Gi);
    else
		  printf("%.18f\t%.18f+%.18fi\n",t*CM2FS,g_r(t,params),tmp_Gi);
    */
    g = G_biexp(t,params);
    printf("%.18f, %.18f %+.18fi\n",t*CM2FS,GSL_REAL(g),GSL_IMAG(g));
	}
}

int main(void){
  // double beta = 0.00478683;// unit: cm // 300K
  double beta = 0.01865; // 77 K
  double lambda0 = 70.0; // reorginization energy cm^-1
  double Gamma0 = 40.0; // cm^-1 lineshape brodening?
  double cotx = 1/tan(beta*Gamma0/2);
  double p[] = {beta,lambda0,Gamma0,cotx};
  printf("Print the lineshape function for check");
  printf("C: t(fs)\tC(cm^-2)\n");
  plot_C(p);

  printf("H: t(fs)\tH(cm^-1)\n");
  plot_H(p);
  printf("G: t(fs)\tG(unitless)\n");
  plot_G(p);

  return 0;
}
#endif
/*
 * $Log$
 * Revision 1.2  2007/01/30 08:33:15  platin
 *   - use more reasonable coupling constant factor in the bath_htgau module,
 *     note that this is still a high-T approximation.
 *   - more informative messages in 2c3ppes-htgau.c
 *
 * Revision 1.1.1.1  2006/05/24 00:42:18  platin
 *
 *   - initial import of the qdas package.
 *   - qdas stands for quantum dynamics and spectroscopy.
 *   - basic code inherited from the "lineshape" package.
 *
 *
 *
 */
