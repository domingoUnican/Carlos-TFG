
/* mymath.h */
#ifndef MYMATH_H
#define MYMATH_H

#include <stddef.h>   // size_t
#include <complex.h>  // double complex

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* MYMATH_H */
