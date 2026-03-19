
/* mymath.c */
#include "mymath.h"
#include <math.h>     // M_PI (si no, usamos acos(-1.0))
#include <complex.h>

void dft(const double complex *x, double complex *X, size_t N) {
    if (N == 0) return;

    // Factor común: -2π/N
    const double two_pi_over_N = 2.0 * acos(-1.0) / (double)N;

    for (size_t k = 0; k < N; ++k) {
        double complex sum = 0.0 + 0.0*I;
        for (size_t n = 0; n < N; ++n) {
            double angle = -two_pi_over_N * (double)(k * n);
            sum += x[n] * cexp(I * angle);
        }
        X[k] = sum;
    }
}

void vec_mat(const double *x, const double *A, double *y, size_t m, size_t n) {
    // y <- 0
    for (size_t j = 0; j < n; ++j) y[j] = 0.0;

    // y[j] += x[i] * A(i,j)
    for (size_t i = 0; i < m; ++i) {
        double xi = x[i];
        const double *Ai = A + i * n; // fila i de A
        for (size_t j = 0; j < n; ++j) {
            y[j] += xi * Ai[j];
        }
    }
}
