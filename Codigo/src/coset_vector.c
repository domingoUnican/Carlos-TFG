#include "coset_vector.h"

// Hagamos otra función
// Hazme una funcion que devuelva una matriz de enteros 
// donde cada fila sea un vector generado a partir de las combinaciones de cosets.

void PSD(const CosetVectors * cv, size_t N, uint8_t resultado[cv->num_vectors][N]) {
    double complex x_temp[N];

    for (size_t i = 0; i < cv->num_vectors; i++) {
        for (size_t j = 0; j < N; j++) {
            x_temp[j] = (double)cv->vectors[i][j] + 0.0 * I;
        }
        psd(x_temp, resultado[i], N);
    }
}

/* Genera todos los vectores basados en combinaciones de cosets */
CosetVectors* generate_coset_vectors(const CosetList *cl, size_t N) {
    if (!cl || cl->len == 0 || N == 0) {
        fprintf(stderr, "ERROR: lista de cosets inválida\n");
        return NULL;
    }
    
    if (cl->len > 30) {  // Límite razonable
        fprintf(stderr, "ERROR: demasiados cosets (%zu)\n", cl->len);
        return NULL;
    }
    
    CosetVectors *result = (CosetVectors*)malloc(sizeof(CosetVectors));
    if (!result) {
        fprintf(stderr, "ERROR: sin memoria para CosetVectors\n");
        return NULL;
    }
    
    result->num_cosets = cl->len;
    result->num_vectors = 1ULL << cl->len;  // 2^num_cosets
    result->vector_length = N;
    
    // Reservar memoria para combinaciones
    result->combinations = (uint8_t**)malloc(result->num_vectors * sizeof(uint8_t*));
    result->vectors = (uint8_t**)malloc(result->num_vectors * sizeof(uint8_t*));
    
    if (!result->combinations || !result->vectors) {
        fprintf(stderr, "ERROR: sin memoria para arrays\n");
        if (result->combinations) free(result->combinations);
        if (result->vectors) free(result->vectors);
        free(result);
        return NULL;
    }
    
    // Inicializar todo a NULL
    memset(result->combinations, 0, result->num_vectors * sizeof(uint8_t*));
    memset(result->vectors, 0, result->num_vectors * sizeof(uint8_t*));
    
    // Generar cada combinación y su vector correspondiente
    for (size_t i = 0; i < result->num_vectors; i++) {
        // 1. Generar combinación binaria para los cosets
        result->combinations[i] = (uint8_t*)malloc(cl->len * sizeof(uint8_t));
        if (!result->combinations[i]) {
            fprintf(stderr, "ERROR: sin memoria para combination %zu\n", i);
            goto error_cleanup;
        }
        
        for (size_t j = 0; j < cl->len; j++) {
            result->combinations[i][j] = (i >> j) & 1;
        }
        
        // 2. Generar vector completo para esta combinación
        result->vectors[i] = generate_vector_for_combination(cl, result->combinations[i], N);
        if (!result->vectors[i]) {
            fprintf(stderr, "ERROR: sin memoria para vector %zu\n", i);
            goto error_cleanup;
        }
    }
    
    return result;
    
error_cleanup:
    // Limpiar memoria parcialmente asignada
    for (size_t j = 0; j < result->num_vectors; j++) {
        if (result->combinations && result->combinations[j]) {
            free(result->combinations[j]);
        }
        if (result->vectors && result->vectors[j]) {
            free(result->vectors[j]);
        }
    }
    if (result->combinations) free(result->combinations);
    if (result->vectors) free(result->vectors);
    free(result);
    return NULL;
}

void free_coset_vectors(CosetVectors *cv) {
    if (!cv) return;
    
    if (cv->combinations) {
        for (size_t i = 0; i < cv->num_vectors; i++) {
            if (cv->combinations[i]) free(cv->combinations[i]);
        }
        free(cv->combinations);
    }
    
    if (cv->vectors) {
        for (size_t i = 0; i < cv->num_vectors; i++) {
            if (cv->vectors[i]) free(cv->vectors[i]);
        }
        free(cv->vectors);
    }
    
    free(cv);
}


/* Imprime todas las combinaciones y sus vectores */
void print_coset_vectors(const CosetVectors *cv, const CosetList *cl) {
    if (!cv || !cl) {
        printf("(datos nulos)\n");
        return;
    }
    
    printf("=== Vectores basados en cosets ===\n");
    printf("Cosets: %zu\n", cl->len);
    printf("Vectores: %zu (2^%zu)\n\n", cv->num_vectors, cl->len);
    
    // Imprimir mapeo coset → elementos
    printf("Mapeo cosets:\n ");
    for (size_t i = 0; i < cl->len; i++) {
        const Coset *c = &cl->data[i];
        printf(" C%zu = { ", i);
        for (size_t j = 0; j < c->len; j++) {
            printf("%zu", c->data[j]);
            if (j + 1 < c->len) printf(", ");
        }
        printf(" }\n ");
    }
    printf("\n");
    
    
    for (size_t i = 0; i < cv->num_vectors; i++) {
        printf("Vector %zu: ", i);
        printf("[");
        for (size_t j = 0; j < cv->num_cosets; j++) {
            printf("%u", cv->combinations[i][j]);
            if (j + 1 < cv->num_cosets) printf(", ");
        }
        printf("] -> ");
        print_vector(cv->vectors[i], cv->vector_length);
        printf("\n");
    }
    printf("\n");
}


void print_first_k_vectors(const CosetVectors *cv, size_t k) {
    if (!cv) {
        printf("(datos nulos)\n");
        return;
    }
    
    size_t limit = (k < cv->num_vectors) ? k : cv->num_vectors;
    
    printf("=== Primeros %zu vectores ===\n", limit);
    for (size_t i = 0; i < limit; i++) {
        printf("Vector %zu: ", i);
        print_vector(cv->vectors[i], cv->vector_length);
        printf("\n");
    }
}


int save_vectors(const CosetVectors *cv, const char *filename) {
    if (!cv || !filename) {
        fprintf(stderr, "ERROR: parámetros inválidos\n");
        return 0;
    }
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "ERROR: no se pudo abrir archivo '%s'\n", filename);
        return 0;
    }
    
    for (size_t i = 0; i < cv->num_vectors; i++) {
        for (size_t j = 0; j < cv->vector_length; j++) {
            fprintf(file, "%u", cv->vectors[i][j]);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
    printf("Vectores guardados\n");
    return 1;
}


/* Versión 2: Guarda módulo DFT redondeado Y segundo fichero con transformación */
int save_dft_both(const CosetVectors *cv, const char *filename1,  const char *filename2, size_t N) {
    if (!cv || !filename1 || !filename2) {
        fprintf(stderr, "ERROR: parámetros inválidos\n");
        return 0;
    }
    
    // Calcular valor constante (N-1)/2
    double constant = (N+1)/2;
    
    FILE *file1 = fopen(filename1, "w");
    FILE *file2 = fopen(filename2, "w");
    
    if (!file1 || !file2) {
        fprintf(stderr, "ERROR: no se pudo abrir archivos\n");
        if (file1) fclose(file1);
        if (file2) fclose(file2);
        return 0;
    }
    
    double complex *time_domain = (double complex*)malloc(N * sizeof(double complex));
    double complex *freq_domain = (double complex*)malloc(N * sizeof(double complex));
    
    if (!time_domain || !freq_domain) {
        fprintf(stderr, "ERROR: sin memoria para DFT\n");
        fclose(file1); fclose(file2);
        if (time_domain) free(time_domain);
        if (freq_domain) free(freq_domain);
        return 0;
    }
    
    for (size_t i = 0; i < cv->num_vectors; i++) {
        
        binary_to_complex(cv->vectors[i], time_domain, N);
        dft(time_domain, freq_domain, N);
        
        
        // Escribir en primer fich: módulo redondeado
        printf("********\n");
        for (size_t j = 1; j < N; j++) {

            double modulo = pow(cabs(freq_domain[j]),2);
            double rounded = rint(modulo);

            printf("%f : %f : %f : %f\n",creal(time_domain[j]), creal(freq_domain[j]), modulo, rounded);
            fprintf(file1, "%.0f", rounded);
            if (j + 1 < N) fprintf(file1, " ");
        }
        fprintf(file1, "\n");
        
        // Escribir en segundo fich: (N-3)/2 - módulo redondeado
        for (size_t j = 1; j < N; j++) {
            double modulo = pow(cabs(freq_domain[j]),2);
            double rounded = rint(modulo);
            double transformed = constant - rounded;
            
            fprintf(file2, "%.0f", rint(transformed));
            
            if (j + 1 < N) fprintf(file2, " ");
        }
        fprintf(file2, "\n");
    }
    
    free(time_domain);
    free(freq_domain);
    fclose(file1);
    fclose(file2);

    printf("Modulo DFT redondeado guardado\n");
    printf("(|(N-3)/2 - modulo)| redondeado guardado\n\n");
    
    return 1;
}
