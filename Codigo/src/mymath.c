
/* mymath.c */
#include "mymath.h"
#include <math.h>     // M_PI (si no, usamos acos(-1.0))
#include <complex.h>

/* Convierte vector binario a complejo (-1/+1) */
void binary_to_complex(const int *binary, double complex *complex_arr, size_t N) {
    for (size_t i = 0; i < N; i++) {
        complex_arr[i] = (binary[i] == 0) ? 1.0 + 0.0*I : 0 + 0.0*I;
    }
}

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

void psd(const double complex *X, uint8_t *psd_out, size_t N)
{
    double complex dft_temp[N];
    dft(X, dft_temp, N);
    for (size_t k = 0; k < N; ++k) {
        psd_out[k] = rint(creal(dft_temp[k]) * creal(dft_temp[k]) + 
                          cimag(dft_temp[k]) * cimag(dft_temp[k]));
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


void print_bits(const uint8_t *v, size_t n, const char *name) {
    printf("%s = [", name);
    for (size_t i = 0; i < n; ++i) {
        printf("%u", (unsigned)v[i]);
        if (i + 1 < n) printf(", ");
    }
    printf("]\n");
}

void print_vector(const uint8_t *v, size_t n) {
    printf("[");
    for (size_t i = 0; i < n; i++) {
        printf("%u", (unsigned)v[i]);
        if (i + 1 < n) printf(", ");
    }
    printf("]");
}

void legendre_sequence(int p, int q, int *sequence, int flag)
{
    sequence[0] = (q*q + 1) >> 1; // esto es dividir por 2, pero usando bit shift para enteros
    for (int n = 1; n < p; n++) {
        int legendre_symbol = flag;
        for (int k = 1; k <= (p - 1) / 2; k++) {
            if ((n % p) == (k * k % p)) {
                legendre_symbol = -legendre_symbol; // cambiar el signo
                break;
            }
        }
        sequence[n] = (q*q + legendre_symbol * q) >> 1;
    }
}


// Compara lexicográficamente dos vectores de enteros de longitud `dim`.
// Retorna -1, 0, o 1 (como strcmp).
int lex_cmp(const int *a, const int *b, size_t dim) {
    for (size_t i = 0; i < dim; i++) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return  1;
    }
    return 0;
}

// Busca `key` en una lista ordenada lexicográficamente (ascendente) de `num_rows`
// vectores, cada uno de longitud `dim`.
// Si `flag == 0`, retorna el índice exacto si se encuentra y -1 si no.
// Si `flag == 1`, retorna el índice donde `key` debe insertarse para mantener
// el orden; si `key` ya existe, retorna su posición.
int binary_search_sorted_pairs(int **sorted_list, size_t num_rows,
                                      size_t dim, const int *key, int flag) {
    if (num_rows == 0) {
        return -1;
    }

    int cmp_last = lex_cmp(key, sorted_list[num_rows - 1], dim);
    if (flag == 1 && cmp_last > 0) {
        return -1;
    }

    int left = 0, right = (int)num_rows;

    while (left < right) {
        int mid = left + (right - left) / 2;
        int cmp = lex_cmp(key, sorted_list[mid], dim);

        if (cmp <= 0) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }

    if (left < (int)num_rows && lex_cmp(key, sorted_list[left], dim) == 0) {
        return left;
    }

    if (flag == 1) {
        return left;
    }

    return -1;
}

void CompressSequence(int N, int p, const int *sequence, int *compressed) {
    if (N % p != 0) {
        fprintf(stderr, "Error: N debe ser múltiplo de p para la compresión.\n");
        return;
    }
    for (int i = 0; i < p; i++) {
        compressed[i] = 0;
        for (int j = 0; j < N / p; j++) {
            compressed[i] += sequence[i + j * p];
        }
    }
}
