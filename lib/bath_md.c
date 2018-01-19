/***************************************************
 * bath_md.c
 *
 * Part of the qdas package;
 * md bath module.
 *
 * By Tseng Wei-Hsiang <dimsplendid@gmail.com>
 *
 ***************************************************/

#define _GNU_SOURCE

#include <ctype.h>
#include <errno.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "aux.h"
#include "bath_md.h"
#include "bath_odbo.h"
#include "qdas.h"

/* Global constant */
/* Integration constants */
#define NWSPACE (100000)
#define EPSABS (1e-5)
#define EPSREL (1e-5)
#define PI 3.14159265359
#define MAZCUT 1 // the sum of Mazbara term from 1 to n

// convert time from cm unit to fs
// CM2FS (fs/cm)
#define CM2FS (5309.1)

double md_kernel_F(double tau, void *p) {
    // exp()[Re{}cos()-Im{}sin()]
    // double imag_part;
	// printf("time: %.4f cm\n", tau);
    double real_part;
    md_par *params = (md_par *)p;
    double lambda0 = params->lambda0;
    double Ea = params->Ea;
    double Eb = params->Eb;
    double C_baab = params->C_baab;
    double C_aaba = params->C_aaba;
    double C_bbba = params->C_bbba;
    double C_aaab = params->C_aaab;
    double C_bbab = params->C_bbab;
    double C_bbaa = params->C_bbaa;
    double C_aaaa = params->C_aaaa;
    double C_bbbb = params->C_bbbb;
    double G_r, G_i, H_r, H_i, C_r, C_i, cosx, sinx, Re_p, Im_p, expx;
    complex double G, H, C;
    G = md_G(tau, params->g);
    H = md_H(tau, params->h);
    C = md_C(tau, params->c);
    G_r = creal(G);
    G_i = cimag(G);
    H_r = creal(H);
    H_i = cimag(H);
    C_r = creal(C);
    C_i = cimag(C);
    expx = exp((2.0 * C_bbaa - C_aaaa - C_bbbb) * G_r);
    cosx = cos((2.0 * C_bbaa - C_aaaa - C_bbbb) * G_i +
               (2.0 * (C_bbaa - C_bbbb) * lambda0 + Eb - Ea) * tau);
    sinx = sin((2.0 * C_bbaa - C_aaaa - C_bbbb) * G_i +
               (2.0 * (C_bbaa - C_bbbb) * lambda0 + Eb - Ea) * tau);

    Re_p =
        C_baab * C_r - ((C_aaba - C_bbba) * H_r * (C_aaab - C_bbab) * H_r -
                        (((C_aaba - C_bbba) * H_i - 2.0 * C_bbba * lambda0) *
                         ((C_aaab - C_bbab) * H_i - 2.0 * C_bbab * lambda0)));
    Im_p = C_baab * C_i - ((C_aaba - C_bbba) * H_r * ((C_aaab - C_bbab) * H_i -
                                                      2.0 * C_bbab * lambda0) +
                           ((C_aaba - C_bbba) * H_i - 2.0 * C_bbba * lambda0) *
                               (C_aaab - C_bbab) * H_r);
    real_part = expx * (Re_p * cosx - Im_p * sinx);

    // imag_part = expx * (Re_p * sinx + Im_p * cosx);
    return real_part; // unitless
}

static complex double discrete2continuous(double tau, complex double array[]) {
    tau = tau * CM2FS;    // unit: fs
	if (tau > 4092.0) {
		tau = 4092.0;
	}
    double time_step = 1; // unit: fs
    uint32_t index = (uint32_t)(tau / time_step);
    double residual = (tau / time_step) - (double)index;
    complex double result =
        array[index] * (1 - residual) + array[index + 1] * residual;
    return result;
}

complex double md_C(double tau, complex double c[4096]) {
    // time unit: cm
    return discrete2continuous(tau, c);
}
complex double md_G(double tau, complex double g[4096]) {
    return discrete2continuous(tau, g);
}
complex double md_H(double tau, complex double h[4096]) {
    return discrete2continuous(tau, h);
}
double md_lambda0_f(complex double h[4096]) {
    return creal(h[4092]);
}
void md_read(char *file_name, complex double array[4096]) {
    double real, imag;
    FILE *file = NULL;
    if ((file = fopen(file_name, "r")) == NULL) {
        printf("open file error !\n");
        exit(EXIT_FAILURE);
    }
    for (uint32_t i = 0;
         i < 4096 && (fscanf(file, "%lf%lfi", &real, &imag) != EOF); i++) {
        array[i] = real + imag * I;
    }
}

/* interface functions */
int bath_md_init_params(const size_t nsize, const double beta,
                        const size_t bath_nparams, const double *bath_params) {
    int idx;
    idx = 0;
    printf("Use mdfitting bath:\n");
    printf("\n");
    printf(" data/c.csv, data/g.csv, data/h.csv ");
    printf("\n");
    printf("\n");
    printf("\n");

    return 0;
}

int bath_md_free_params() {
    return 0;
}
#ifdef MAIN
typedef complex double (*lineshape_f)(double tau, complex double data[4096]);
void md_plot(complex double data[4096], lineshape_f F) {
    double t = 0.0;
    complex double c;
    for (uint32_t i = 0; i < 3767; i++) {
        t = ((double)i) / 10000.0;
        c = F(t, data);
        printf("%.18f, %.18f%+.18fi\n", t * CM2FS, creal(c), cimag(c));
    }
}

int main(void) {
    complex double c[4096], h[4096], g[4096];
    md_read("../data/g.csv", g);
    md_read("../data/c.csv", c);
    md_read("../data/h.csv", h);
    printf("Print the lineshape function for check");
    printf("C: t(fs)\tC(cm^-2)\n");
    md_plot(c, md_C);
    printf("H: t(fs)\tH(cm^-1)\n");
    md_plot(h, md_H);
    printf("G: t(fs)\tG(unitless)\n");
    md_plot(g, md_G);

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
