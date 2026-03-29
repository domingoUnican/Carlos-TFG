
/* mymath.h */
#ifndef MYMATH_H
#define MYMATH_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <stdint.h>
#include <stddef.h>   // size_t
#include <stdint.h>   // uint8_t
#include <complex.h>  // double complex

#ifdef __cplusplus
extern "C" {
#endif


void binary_to_complex(const int *binary, double complex *complex_arr, size_t N);


/**
 * @brief Calcula la DFT (Transformada Discreta de Fourier) de N muestras complejas.
 *
 * Convención: X[k] = sum_{n=0}^{N-1} x[n] * exp(-j * 2π * k * n / N)
 *
 * @param x   Vector de entrada (longitud N), complejo.
 * @param X   Vector de salida (longitud N), complejo.
 * @param N   Número de muestras.
 */
void dft(const double complex *x, double complex *X, size_t N);


void psd(const double complex *X, uint8_t *psd_out, size_t N);  


/**
 * @brief Multiplica un vector fila por una matriz en orden por filas (row-major).
 *
 * Calcula y = x^T * A
 *  - x: longitud m
 *  - A: matriz m×n en row-major (A[i*n + j] = A(i,j))
 *  - y: longitud n
 *
 * @param x   Vector de entrada (m).
 * @param A   Matriz de entrada (m×n), row-major.
 * @param y   Vector de salida (n).
 * @param m   Número de filas de A (longitud de x).
 * @param n   Número de columnas de A (longitud de y).
 */
void vec_mat(const double *x, const double *A, double *y, size_t m, size_t n);

void print_bits(const uint8_t *v, size_t n, const char *name);
void print_vector(const uint8_t *v, size_t n);
void legendre_sequence(int p, int q, int *sequence, int flag);
int lex_cmp(const int *a, const int *b, size_t dim);
int binary_search_sorted_pairs(int **sorted_list, size_t num_rows,
                                      size_t dim, const int *key, int flag);
void CompressSequence(int N, int p, const int *sequence, int *compressed);
#ifdef __cplusplus
}
#endif

#endif /* MYMATH_H */
