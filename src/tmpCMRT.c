/******************************************************
 * [FILE]     MDCMRT.c
 * [AUTHOR]   Tseng Wei-Hsiang
 * [DATE]     20161019
 * Use the CMRT to calculate the exciton transfer rate
 ******************************************************/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
// #include <errno.h>
#include <ctype.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <string.h>
// #include <gsl/gsl_sf_bessel.h>

#include "tmpCMRT.h"
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
void print_input_tmpCMRT(qdas_keys *Keys) {
    /* the simple form from params.c */
    // int i,j;
    int n;
    n = Keys->nsize;
    printf("\n");
    printf("NSIZE = %d\n", n);
    printf("\n");
    printf("\n");
    printf("-----------------------\n");
    printf("| Hamiltonian (cm^-1) |\n");
    printf("-----------------------\n");
    printf("\n");
    gsl_matrix_print(Keys->He);
    printf("\n");
}

/* Display usage information and exit. */
static void usage() {
    printf("\
    Usage: %s KEYFILE \n\
    Use the parameters described in the KEYFILE and data in data/*.csv to\n\
    calculate MD-CMRT exciton transfer rate\n\
    \n",
           program_name);
}

/* Main program */
int main(int argc, char *argv[]) {
    char *key_file_name = NULL;
    // char *cache_file_name=NULL;
    /* variables for the options */

    /* main keyword structure */
    qdas_keys keys;
    /* set program name for messages. */
    program_name = argv[0];
    /* HEADER MESSAGE */
    printf("\n");
    printf("%s\n", program_name);
    printf("\n");
    printf("        MD-CMRT for rate\n");
    printf("        Copyright(C) 2016.\n");
    printf("        Tseng Wei-Hsiang <dimsplendid@gmail.com>.\n");
    /* parse the arguments... */

    /* check if we have a key file,
     * if yes, use the key file to initialize all parameters */
    if (key_file_name == NULL) {
        key_file_name = (char *)strdup(argv[1]);
    } else {
        usage();
        exit(EXIT_FAILURE);
    }
    if (key_file_name != NULL) {
        params_init(key_file_name, &keys);
    } else {
        usage();
        exit(EXIT_FAILURE);
    }
    /* parameters read and initialized */
    print_input_tmpCMRT(&keys);

    int nsize = keys.nsize;
    gsl_matrix *rate = gsl_matrix_alloc(nsize, nsize);
    gsl_matrix *He = gsl_matrix_alloc(nsize, nsize);
    gsl_matrix_memcpy(He, keys.He);
    /* eigenvector eigen value */
    gsl_vector *eval = gsl_vector_alloc(nsize);
    gsl_matrix *evec = gsl_matrix_alloc(nsize, nsize);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(nsize);

    gsl_eigen_symmv(He, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);
    
    printf("eigenvalue");
    for (uint32_t i = 0; i < nsize; i++) {
        printf("eigenvalue = %g\n", gsl_vector_get(eval, i));
    }
    printf("eigenvector: \n");
    gsl_matrix_print(evec);
    printf("------------------------------------\n");
    

    md_par *p = calloc(sizeof(md_par), 1); // define parameters

    printf("rate matrix: \n");
    gsl_matrix *rate_matrix = gsl_matrix_alloc(nsize, nsize);
    cal_rate_matrix(evec, eval, p, rate_matrix);
    
	free(p);
	gsl_matrix_free(rate_matrix);

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(He);
    gsl_matrix_free(rate);
    return 0;
}

// aux function
void cal_rate_matrix(gsl_matrix *evec, gsl_vector *eval, md_par *p,
                     gsl_matrix *result) {
    int size = evec->size1;
    gsl_function F;
    double ra = 0.0, rb = 0.2;
    double element, abserr;
	md_read("../data/g.csv", p->g);
	md_read("../data/h.csv", p->h);
	md_read("../data/c.csv", p->c);
	p->lambda0 = md_lambda0_f(p->h);
    for (uint32_t a = 0; a < size; a++) {
        for (uint32_t b = 0; b < size; b++) {
            p->Ea = gsl_vector_get(eval, a);
            p->Eb = gsl_vector_get(eval, b);
            p->C_baab = transCoeff(evec, b, a, a, b);
            p->C_aaba = transCoeff(evec, a, a, b, a);
            p->C_bbba = transCoeff(evec, b, b, b, a);
            p->C_aaab = transCoeff(evec, a, a, a, b);
            p->C_bbab = transCoeff(evec, b, b, a, b);
            p->C_bbaa = transCoeff(evec, b, b, a, a);
            p->C_aaaa = transCoeff(evec, a, a, a, a);
            p->C_bbbb = transCoeff(evec, b, b, b, b);

            F.function = &md_kernel_F;
            F.params = p;

            gsl_integration_workspace *w =
                gsl_integration_workspace_alloc(NWSPACE);

            gsl_integration_qags(&F, ra, rb, EPSABS / 100.0, EPSREL / 100.0,
                                 NWSPACE, w, &element, &abserr);
            gsl_matrix_set(result, a, b, 2.0 * element / CM2FS);
            gsl_integration_workspace_free(w);
        }
    }
    gsl_matrix_print(result);
    printf("\n");
}
