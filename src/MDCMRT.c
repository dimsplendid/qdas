/******************************************************
 * [FILE]     MDCMRT.c
 * [AUTHOR]   Tseng Wei-Hsiang
 * [DATE]     20161019
 * Use the CMRT to calculate the exciton transfer rate
 ******************************************************/
#include <stdio.h>
#include <stdlib.h>
// #include <errno.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_trig.h>
// #include <gsl/gsl_sf_bessel.h>

#include "MDCMRT.h"
/* Global constants */

/* Global variables */
const char *program_name;

/* Global constant */
/* Integration constants */
#define NWSPACE (100000000)
#define EPSABS (1e-5)
#define EPSREL (1e-5)
#define PI 3.14159265359


// convert time from cm unit to fs
#define CM2FS (5309.1)
/* Macros */

/* Aux function */
void print_input_MDCMRT(qdas_keys * Keys)
{
	/* the simple form from params.c */
	// int i,j;
	int n;
	n=Keys->nsize;
	printf("\n");
	printf("NSIZE = %d\n",n);
	printf("\n");
	printf("\n");
	printf("-----------------------\n");
	printf("| Hamiltonian (cm^-1) |\n");
	printf("-----------------------\n");
	printf("\n");
	gsl_matrix_print(Keys->He);
	printf("\n");
	printf("------------------------------\n");
	printf("|Beta|\n");
	printf("------------------------------\n");
	printf("%.6f",Keys->beta);
	printf("\n");
}

/* Display usage information and exit. */
static void usage()
{
		printf( "\
    Usage: %s KEYFILE \n\
    Use the parameters described in the KEYFILE to\n\
    calculate MD-CMRT exciton transfer rate\n\
    \n",
		program_name);
}

/* Main program */
int main(int argc, char *argv[]){
	char *key_file_name=NULL;
	//char *cache_file_name=NULL;
	/* variables for the options */

	/* main keyword structure */
	qdas_keys keys;
	/* set program name for messages. */
	program_name = argv[0];
	/* HEADER MESSAGE */
	printf("\n");
	printf("%s\n",program_name);
	printf("\n");
	printf("        MD-CMRT for rate\n");
	printf("        Copyright(C) 2016.\n");
	printf("        Tseng Wei-Hsiang <dimsplendid@gmail.com>.\n");
	/* parse the arguments... */

	/* check if we have a key file,
	 * if yes, use the key file to initialize all parameters */
	if(key_file_name == NULL) {
		key_file_name=(char*)strdup(argv[1]);
	}
	else {
		usage();
		exit(EXIT_FAILURE);
	}
	if(key_file_name != NULL) {
		params_init(key_file_name,&keys);
	}
	else {
		usage();
		exit(EXIT_FAILURE);
	}
	/* parameters read and initialized */
	print_input_MDCMRT(&keys);

	int nsize = keys.nsize;
	gsl_matrix * rate = gsl_matrix_alloc(nsize,nsize);
	gsl_matrix * He = gsl_matrix_alloc(nsize,nsize);
	gsl_matrix_memcpy(He, keys.He);
	double beta = keys.beta;
	/* eigenvector eigen value */
	gsl_vector *eval = gsl_vector_alloc (nsize);
	gsl_matrix *evec = gsl_matrix_alloc (nsize, nsize);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (nsize);

	gsl_eigen_symmv (He, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	{
	  int i;
	  for (i = 0; i < nsize; i++)
	    {
	      double eval_i = gsl_vector_get (eval, i);
	      gsl_vector_view evec_i = gsl_matrix_column (evec, i);
	      printf ("eigenvalue = %g\n", eval_i);
	      printf ("eigenvector = \n");
	      gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
	    }
	}
	printf("eigenvector: \n");
	gsl_matrix_print(evec);
	printf("------------------------------------\n");
	double tan_C = 2.65459449960518*beta; // beta*h_bar/2 unit: ps
	double lambda0 = lambda0_f(tan_C)*8000; // change unit to cm-1, need check

	double p[] = {beta,lambda0,tan_C}; //define parameters
	printf("reorginization energy: %.18f\n",lambda0);

	printf("rate matrix: \n");
	gsl_matrix * rate_matrix = gsl_matrix_alloc(nsize,nsize);
	cal_rate_matrix(evec,eval,p,rate_matrix);
	gsl_matrix_free (rate_matrix);

	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	gsl_matrix_free (He);
	gsl_matrix_free (rate);
	return 0;
}

// aux function
void plot_kernel(gsl_matrix * evec,gsl_vector * eval,double params[],int a,int b){
  int i = 0;
  double p[14];
  double t;
  printf("a: %d, b: %d\n",a,b);
  double C_baab = transCoeff(evec,b,a,a,b);
  double C_aaba = transCoeff(evec,a,a,b,a);
  double C_bbba = transCoeff(evec,b,b,b,a);
  double C_aaab = transCoeff(evec,a,a,a,b);
  double C_bbab = transCoeff(evec,b,b,a,b);
  double C_bbaa = transCoeff(evec,b,b,a,a);
  double C_aaaa = transCoeff(evec,a,a,a,a);
  double C_bbbb = transCoeff(evec,b,b,b,b);
  p[0] = params[0];
  p[1] = params[1];
  p[2] = params[2];
  p[3] = params[3];
  p[4] = gsl_vector_get(eval,a);
  p[5] = gsl_vector_get(eval,b);
  p[6] = C_baab;
  p[7] = C_aaba;
  p[8] = C_bbba;
  p[9] = C_aaab;
  p[10] = C_bbab;
  p[11] = C_bbaa;
  p[12] = C_aaaa;
  p[13] = C_bbbb;
  printf("plot ...\n");
  for(i = 0;i < 1000; i+=1){
    t = ((double)i) / 10000.0;
    printf("%.18f %.18f\n", t*CM2FS,mdfit_kernel_F(t,p));
    // printf("%.18f %.18f\n", t,kernel_F(t,p));
  }

  printf("\n");
}
void cal_rate_matrix(gsl_matrix * evec,gsl_vector * eval,double params[] ,gsl_matrix * result){
  double p[14];
  int size = evec->size1;
  int a,b;
  gsl_function F;
  double ra=0.0,rb=0.2;
  double element, abserr;
  for(a = 0;a < size; a++ ){
    for(b = 0; b < size; b++){
      p[0] = params[0];
			p[1] = params[1];
			p[2] = params[2];
      p[4] = gsl_vector_get(eval,a);
      p[5] = gsl_vector_get(eval,b);
      p[6] = transCoeff(evec,b,a,a,b);
      p[7] = transCoeff(evec,a,a,b,a);
      p[8] = transCoeff(evec,b,b,b,a);
      p[9] = transCoeff(evec,a,a,a,b);
      p[10] = transCoeff(evec,b,b,a,b);
      p[11] = transCoeff(evec,b,b,a,a);
      p[12] = transCoeff(evec,a,a,a,a);
      p[13] = transCoeff(evec,b,b,b,b);

      F.function = &mdfit_kernel_F;
      F.params = p;

      gsl_integration_workspace *w = gsl_integration_workspace_alloc(NWSPACE);

      gsl_integration_qags(&F,ra,rb,EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &element, & abserr);
      gsl_matrix_set(result,a,b,2.0*element/CM2FS);
    }
  }

  gsl_matrix_print(result);
  printf("\n");
}
