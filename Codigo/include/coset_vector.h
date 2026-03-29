#ifndef COSET_VECTOR_H
#define COSET_VECTOR_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include "cyclotomic_cosets.h"
#include <stdlib.h>
#include <stdbool.h>    
#include <complex.h>


/* Versión modificada de BinaryCombinations que genera vectores completos */
typedef struct {
    uint8_t **vectors;        // Array de vectores de longitud N
    uint8_t **combinations;   // Array de combinaciones de cosets (para referencia)
    size_t num_vectors;       // Número total de vectores (2^num_cosets)
    size_t vector_length;     // Longitud de cada vector (N)
    size_t num_cosets;        // Número de cosets
} CosetVectors;



void PSD(const CosetVectors * cv, size_t N, uint8_t resultado[cv->num_vectors][N]);

CosetVectors* generate_coset_vectors(const CosetList *cl, size_t N);

void free_coset_vectors(CosetVectors *cv);

void print_coset_vectors(const CosetVectors *cv, const CosetList *cl);

void print_first_k_vectors(const CosetVectors *cv, size_t k);

int save_vectors(const CosetVectors *cv, const char *filename); 

int save_dft_both(const CosetVectors *cv, const char *filename1,  const char *filename2, size_t N);
#endif /* COSET_VECTOR_H */