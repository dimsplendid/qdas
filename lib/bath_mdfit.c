/***************************************************
 * bath_mdfit.c
 *
 * Part of the qdas package;
 * Time correation function from MD fitting by 王佑仁
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
#include "bath_mdfit.h"
#include "mdfit_par.h"

/* Define gsl macro */
#define ADD gsl_complex_add
#define MUL gsl_complex_mul
#define SUB gsl_complex_sub
#define INV gsl_complex_inverse
#define ADDR gsl_complex_add_real
#define SUBR gsl_complex_sub_real
#define MULR gsl_complex_mul_real
#define EXP gsl_complex_exp
#define TAN gsl_complex_tan
/* Global constant */
/* Integration constants */
#define NWSPACE (100000)
#define EPSABS (1e-5)
#define EPSREL (1e-5)
#define PI 3.14159265359

// convert time from cm unit to fs
#define CM2FS (5309.1)
// change the unit to cm^-2, = <C(0)> = varience of energy, [cm-2]
// #define UNITTRANS 20000
#define UNITTRANS 6.660363243376534e+02
double mdfit_kernel_F(double tau, void * p){
  // exp()[Re{}cos()-Im{}sin()]
  double real_part;
  double t_ps;
  double * params = (double*)p;
  double beta = params[0];
  double lambda0 = params[1];
  double tan_C = params[2]; // beta*h_bar/2 unit: ps
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
  double G_r,G_i,H_r,H_i,C_r,C_i,cosx,sinx,Re_p,Im_p,expx;

  gsl_complex G,H,C;
  t_ps = tau*CM2FS/1000; // transfer time from cm to ps
  G = MULR(mdfit_G(t_ps,beta),UNITTRANS*1e6/CM2FS/CM2FS);// transfer [cm^-2 ps^2] to [unitless]
  H = MULR(mdfit_H(t_ps,beta),UNITTRANS*1000/CM2FS);// transfer [cm^-2 ps] to [cm^-1]
  C = MULR(mdfit_C(t_ps,beta),UNITTRANS);// transfer [cm^-2] to [cm^-2]

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
        - (((C_aaba-C_bbba)*H_i - 2.0*C_bbba*lambda0)*((C_aaab-C_bbab)*H_i - 2.0*C_bbab*lambda0)));
  Im_p = C_baab*C_i - ((C_aaba - C_bbba)*H_r*((C_aaab-C_bbab)*H_i - 2.0*C_bbab*lambda0)\
        + ((C_aaba-C_bbba)*H_i - 2.0*C_bbba*lambda0)*(C_aaab - C_bbab)*H_r);
  real_part = expx * (Re_p * cosx - Im_p * sinx);
  // double imag_part;
  //imag_part = expx * (Re_p * sinx + Im_p * cosx);
  return real_part;
}
gsl_complex mdfit_G(double t,double tan_C){
  gsl_complex A,B,G,Gr,Gi,tmp,tanx;
  GSL_SET_COMPLEX(&G,0.0,0.0);
  GSL_SET_COMPLEX(&Gr,0.0,0.0);
  GSL_SET_COMPLEX(&Gi,0.0,0.0);
  for(uint32_t i = 0; i < FIT_SIZE; i++){
    GSL_SET_COMPLEX(&A,A_r[i],A_i[i]);
    GSL_SET_COMPLEX(&B,B_r[i],B_i[i]);
    tmp=MUL(A,MUL(MUL(INV(B),INV(B)),SUBR(SUB(EXP(MULR(B,t)),MULR(B,t)),1)));
    tanx=TAN(MULR(B,tan_C));
    Gr=ADD(Gr,tmp);
    Gi=ADD(Gi,MUL(tmp,tanx));
  }
  GSL_SET_COMPLEX(&G,GSL_REAL(Gr),GSL_REAL(Gi));
  return G;
}
gsl_complex mdfit_H(double t, double tan_C){
  gsl_complex A,B,Hr,Hi,H,tmp,tanx;
  GSL_SET_COMPLEX(&H,0.0,0.0);
  GSL_SET_COMPLEX(&Hi,0.0,0.0);
  GSL_SET_COMPLEX(&Hr,0.0,0.0);
  for(uint32_t i = 0; i < FIT_SIZE; i++){
    GSL_SET_COMPLEX(&A,A_r[i],A_i[i]);
    GSL_SET_COMPLEX(&B,B_r[i],B_i[i]);
    tmp=MUL(A,MUL(INV(B),SUBR(EXP(MULR(B,t)),1)));
    tanx=TAN(MULR(B,tan_C));
    Hr=ADD(Hr,tmp);
    Hi=ADD(Hi,MUL(tmp,tanx));
  }
  GSL_SET_COMPLEX(&H,GSL_REAL(Hr),GSL_REAL(Hi));
  return H;
}
gsl_complex mdfit_C(double t, double tan_C){
  gsl_complex A,B,Cr,Ci,C,tmp,tanx;
  /* initialize */
  GSL_SET_COMPLEX(&C,0.0,0.0);
  GSL_SET_COMPLEX(&Cr,0.0,0.0);
  GSL_SET_COMPLEX(&Ci,0.0,0.0);
  for(uint32_t i = 0;i < FIT_SIZE; i++){
    GSL_SET_COMPLEX(&A,A_r[i],A_i[i]);
    GSL_SET_COMPLEX(&B,B_r[i],B_i[i]);
    //printf("A[%d]: ",i);
    //gsl_complex_print(A);
    //printf("B[%d]: ",i);
    //gsl_complex_print(B);
    tmp=MUL(A,EXP(MULR(B,t)));
    // printf("Cr part: ");
    // gsl_complex_print(tmp);
    /* What's wrong with imagenary part???? */

    tanx=TAN(MULR(B,tan_C));
    // printf("tanx: ");
    // gsl_complex_print(tanx);
    Cr=ADD(Cr,tmp);
    Ci=ADD(Ci,MUL(tmp,tanx));
    // printf("Ci part:");
    // gsl_complex_print(MUL(tmp,tanx));
  }
  GSL_SET_COMPLEX(&C,GSL_REAL(Cr),GSL_REAL(Ci));
  return C;
}

double lambda0_f(double tan_C){
  gsl_complex A,B,Hi,tmp,tanx;
  GSL_SET_COMPLEX(&Hi,0.0,0.0);
  for(uint32_t i = 0; i < FIT_SIZE; i++){
    GSL_SET_COMPLEX(&A,A_r[i],A_i[i]);
    GSL_SET_COMPLEX(&B,B_r[i],B_i[i]);
    tmp=MULR(MUL(A,INV(B)),-1.0);
    tanx=TAN(MULR(B,tan_C));
    Hi=ADD(Hi,MUL(tmp,tanx));
  }
  return fabs(UNITTRANS*1000/CM2FS*GSL_REAL(Hi));
}

int bath_mdfit_init_params(const size_t nsize, const double beta,
			   const size_t bath_nparams, const double *bath_params){
  printf("Use MD fitting bath:\n");
  printf("\n");

  return 0;
}

int bath_mdfit_free_params(){
  return 0;
}

#ifdef MAIN
void plot_C_r(void * params){
	double t, tan_C;
	double * p = (double *)params;
	tan_C = p[0];
	for(uint32_t i = 0; i < 2000; i++){
		t = ((double)i)/1000.0;
		printf("%.18f\t%.18f\n",t,GSL_REAL(mdfit_G(t,tan_C)));
	}
}
void plot_C_i(void * params){
	double t, tan_C;
	double * p = (double *)params;
	tan_C = p[0];
	for(uint32_t i = 0; i < 2000; i++){
		t = ((double)i)/1000.0;
		printf("%.18f\t%.18f\n",t,GSL_IMAG(mdfit_G(t,tan_C)));
	}
}

int main(void){
  double beta = 0.00478683;// unit: cm // 300K
  double tan_C = 2.65459449960518*beta; // beta*h_bar/2 unit: ps
  double * p = &tan_C;

  printf("test mdfit bath:\n");
  printf("Gr: t(ps)\tGr(t)\n");
  plot_C_r(p);
  printf("Gi: t(ps)\tGi(t)\n");
  plot_C_i(p);
  return 0;
}
#endif
