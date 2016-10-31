/******************************************************
 * [FILE]     CMRT.c
 * [AUTHOR]   Tseng Wei-Hsiang
 * [DATE]     20151015
 * Use the CMRT to calculate the exciton transfer rate
 ******************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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

#include "CMRT.h"
/* Global constants */

/* Global variables */
const char *program_name;
extern double BathODBOLambda;
extern double BathODBOGAMMA;
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
void print_input_CMRT(qdas_keys * Keys)
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
	printf("| lambda_0\tGamma_0\tBeta|\n");
	printf("------------------------------\n");
	printf("   %.3f\t%.3f\t%.6f",BathODBOLambda,BathODBOGAMMA,Keys->beta);
	printf("\n");
}

/* Display usage information and exit. */
static void usage()
{
		printf( "\
    Usage: %s [OPTIONS] KEYFILE \n\
    Use the parameters described in the KEYFILE to\n\
    calculate CMRT exciton transfer rate\n\
    \n\
    OPTIONS:\n\
    \n\
    --nocache  Do not cache polarizability data.\n\
    --resume   Resume calculation if a .qcache file is found.\n",
		program_name);
}

/* the following functions handle cache files */
/* read the cache file; if not successful, return 0; *iter=0 too, if not succeful */
// int read_ptavg_cache(char *fname, size_t *iter, gsl_vector_complex *Pt)
// {
// 	FILE *input_file;
// 	input_file=fopen(fname, "r");
// 	if(input_file == NULL) {
// 			*iter=0;
// 			return 0;
// 	}
// 	int err = fread(iter, sizeof(size_t), 1, input_file);
// 	gsl_vector_complex_fread(input_file,Pt);
// 	fclose(input_file);
// 	return *iter;
// }
// void save_ptavg_cache(char *fname, size_t iter, gsl_vector_complex *Pt)
// {
// 	FILE *output_file;
// 	output_file=fopen(fname, "w");
// 	if(output_file == NULL) {
// 			fprintf(stderr, "error while opening \"%s\" for writing: %s \n",
// 							fname,strerror(errno));
// 			exit(errno);
// 	}
// 	fwrite(&iter, sizeof(size_t), 1, output_file);
// 	gsl_vector_complex_fwrite(output_file,Pt);
// 	fclose(output_file);
// }

/* Main program */
int main(int argc, char *argv[]){
	int i;
	char *key_file_name=NULL;
	char *cache_file_name=NULL;
	/* variables for the options */
	// int resume=0;
	// int cache=1;
	/* main keyword structure */
	qdas_keys keys;
	/* set program name for messages. */
	program_name = argv[0];
	/* HEADER MESSAGE */
	printf("\n");
	printf("%s\n",program_name);
	printf("\n");
	printf("        CMRT for rate\n");
	printf("        Copyright(C) 2015.\n");
	printf("        Tseng Wei-Hsiang <dimsplendid@gmail.com>.\n");
	/* parse the arguments... */
	double kernel_cut = 500; // default cut
	for(i=1;i<argc;i++) {
				if(argv[i][0] != '-') {
						if(key_file_name == NULL) {
								key_file_name=(char*)strdup(argv[i]);
								cache_file_name=(char*)malloc((strlen(key_file_name)+8)*sizeof(char));
								strcpy(cache_file_name,key_file_name);
								strcpy(cache_file_name+strlen(key_file_name),".qcache");
						}
						else {
								usage();
								exit(EXIT_FAILURE);
						}
				}
				else {
						// if(!strncmp(argv[i],"--resume",8)) {
						// 		/* resume previously interrupted calculation
						// 		 * when a cache file is found */
						// 		resume = 1;
						// }
						// else if(!strncmp(argv[i],"--nocache",9)) {
						// 		/* do not save cache file. */
						// 		cache = 0;
						// }
						if(!strncmp(argv[i],"--cut",5)){
								i++;
								kernel_cut = char2double(argv[i]);
						}
						else {
								/* new options go here */
								usage();
								exit(EXIT_FAILURE);
						}
				}
	}
	/* check if we have a key file,
	 * if yes, use the key file to initialize all parameters */
	if(key_file_name != NULL) {
			params_init(key_file_name,&keys);
	}
	else {
			usage();
			exit(EXIT_FAILURE);
	}
	/* parameters read and initialized */
	print_input_CMRT(&keys);

	int nsize = keys.nsize;
	gsl_matrix * rate = gsl_matrix_alloc(nsize,nsize);
	gsl_matrix * He = gsl_matrix_alloc(nsize,nsize);
	gsl_matrix_memcpy(He, keys.He);
	double beta = keys.beta;
	double lambda0 = BathODBOLambda;
	double Gamma0 = BathODBOGAMMA;
	double cotx = cos(beta*Gamma0/2)/sin(beta*Gamma0/2);
	/* eigenvector eigen value */
	gsl_vector *eval = gsl_vector_alloc (nsize);
	gsl_matrix *evec = gsl_matrix_alloc (nsize, nsize);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (nsize);

	gsl_eigen_symmv (He, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);
	printf("eigenvector: \n");
	gsl_matrix_print(evec);
	printf("------------------------------------\n");
	double p[] = {beta,lambda0,Gamma0,cotx}; //define parameters

	printf("rate matrix: \n");
	gsl_matrix * rate_matrix = gsl_matrix_alloc(nsize,nsize);
	cal_rate_matrix(evec,eval,p, kernel_cut,rate_matrix);
	gsl_matrix_free (rate_matrix);

	gsl_vector_free (eval);
	gsl_matrix_free (evec);
	gsl_matrix_free (He);
	gsl_matrix_free (rate);
	return 0;
}

// aux function
double char2double(const char *cline)
{

// converts character table to the double


	double value,decimal,divided;
	int number, negflag;

	char c;

	c = *cline;
	value = 0;

	negflag = 0;
	if (c=='-') {
		negflag = 1;
		cline++;
		c = *cline;
		}

        while ((c != 0) && (c != '.')) {


	   number = c - '0';			// i:s number
           value += number;

	   cline++;
	   c = *cline;

	   if ((c != 0) && (c != '.'))	{
		value *= 10;	// more values to come
		}
 	}

// decimals:


	if (c != '.') { 			// no decimals
		if (negflag != 0) value = 0 - value;		//negativity
		return value;
		}

	cline++;		// jump the decimalpoint


	divided = 10;
	c = *cline;
	while (c != 0) {
	   number = c - '0';

  	   if (number != 0) {
	   decimal = number;
	   decimal /= divided;
	   value += decimal;
	}

	   divided *= 10;
	   cline++;
	   c = *cline;
	}


	if (negflag != 0) value = 0 - value;		// negativity
        return value;

}

// aux function
void plot_kernel(gsl_matrix * evec,gsl_vector * eval,double params[], double cut,int a,int b){
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
    printf("%.18f %.18f\n", t*CM2FS,kernel_F(t,p));
    // printf("%.18f %.18f\n", t,kernel_F(t,p));
  }

  printf("\n");
}
void cal_rate_matrix(gsl_matrix * evec,gsl_vector * eval,double params[], double cut,gsl_matrix * result){
  double p[14];
  int size = evec->size1;
  gsl_function F;
  double ra=0.0,rb=0.2;
  double element, abserr;
  for(uint_fast32_t a = 0;a < size; a++ ){
    for(uint_fast32_t b = 0; b < size; b++){
      p[0] = params[0];
      p[1] = params[1];
      p[2] = params[2];
      p[3] = params[3];
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

      F.function = &kernel_F;
      F.params = p;

      gsl_integration_workspace *w = gsl_integration_workspace_alloc(NWSPACE);

      gsl_integration_qags(&F,ra,rb,EPSABS/100.0, EPSREL/100.0, NWSPACE, w, &element, & abserr); // int from ra to rb
			// gsl_integration_qagiu (&F,ra, EPSABS, EPSREL, NWSPACE, w, &element, &abserr); // int from ra to inf
      gsl_matrix_set(result,a,b,2.0*element/CM2FS);
			// check
			// printf("H(%lu,%lu) = %.6f\n",a,b,element);
    }
  }

  gsl_matrix_print(result);
  printf("\n");
}
