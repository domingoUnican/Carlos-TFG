#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mymath.h"
#include <complex.h>
#include <math.h>
#include "pairs_reader.h"
#define ALFABET_SIZE 2



/* int* generate_vector_for_combination(const CosetList *cl,  */
/*                                         const int *combination,  */
/*                                         size_t N) { */
/*     int *vector = (int*)calloc(N,sizeof(int)); */
/*     if (!vector) { */
/*         fprintf(stderr, "ERROR: sin memoria para vector\n"); */
/*         return NULL; */
/*     } */
    
/*     // Para cada elemento 0..N-1, asignamos el valor según su coset */
/*     for (size_t elem = 0; elem < N; elem++) { */
/*         int coset_idx = find_element_in_cosets(cl, elem); */
/*         if (coset_idx >= 0) { */
/*             vector[elem] = combination[coset_idx]; */
/*         } */
/*         // Si no se encuentra (no debería pasar), queda en 0 */
/*     } */
    
/*     return vector; */
/* } */

/* static void binary_to_complex(const int *binary, double complex *complex_arr, size_t N) { */
/*     for (size_t i = 0; i < N; i++) { */
/*         complex_arr[i] = (binary[i] == 0) ? 1.0 : 0.0; */
/*     } */
/* } */


/* bool is_less_than_compression(int N, int p, const int *sequence, const int *compressed, const int *bound_sequence) */
/* { */
/*     //return true; */
/*     int *temp = calloc(p, sizeof(int)); */
/*     int *temp_bound = calloc(p, sizeof(int)); */
/*     for (int i = 0; i < p; i++) */
/*     { */
/*         for (int j = 0; p * j < N; j++) */
/*         { */
/*             temp[i] += sequence[(i + p * j) % N]; */
/*             temp_bound[i] += bound_sequence[(i + p * j) % N]; */
/*         } */
/*         if (temp[i] > compressed[i] || temp_bound[i] + temp[i] < compressed[i])  */
/*         { */
/*             free(temp); */
/*             free(temp_bound); */
/*             return false; */
/*         } */
/*     } */
/*     free(temp); */
/*     free(temp_bound); */
/*     return true; */
/* } */

/* typedef struct { */
/*     FILE *f_comb; */
/*     FILE *f_psd; */
/*     FILE *f_cte; */
/*     size_t matches_found; */
/*     double threshold; */
/*     //double constant; constant is the same as threshold.  */
/*     double complex **dft_matrix; */
/*     int **psd_matrix; */
/*     size_t N; */
/*     size_t num_cosets; */
/*     double complex *current_dft; */
/*     int *current_psd; */
/*     int *bound_psd; */
/*     int *compression_a; */
/*     int *compression_b; */
/*     int p; */
/*     int q; */
/*     int coset_idx; */
/*     CosetList *cl; */
/*     int *current_combination; */
/*     int **candidate_pairs1; */
/*     size_t num_candidate_pairs1; */
/*     size_t dimension_candidate_pairs1; */
/*     int **candidate_pairs2; */
/*     size_t num_candidate_pairs2; */
/*     size_t dimension_candidate_pairs2; */
/* } DFSContext; */

/* bool is_compression(int N, int p, int *sequence, int *compressed) */
/* { */
/*     int temp[p]; */
/*     for (int i = 0; i < p; i++) { */
/*         temp[i] = compressed[i]; */
/*     } */
/*     for (int i = 0; i < p; i++) */
/*     { */
/*         for (int j = 0; p * j < N; j++) */
/*         { */
/*             temp[i] -= sequence[(i + p * j) % N]; */
/*         } */
/*         if (temp[i] != 0) */
/*         { */
/*             return false; */
/*         } */
/*     } */
/*     return true; */

/* } */
/* bool is_less_than_candidates(const DFSContext *ctx, const int *sequence, const int *bound_bit_sequence, int flag) */
/* { */
/*     int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2; */
/*     int num_candidate_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2; */
/*     int copia[ctx->N]; */
/*     int bound_sequence[ctx->N]; */
/*     int **pairs_temp = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2; */
/*     /\*printf("****************************************************************\n"); */
/*     printf("Comprobando compresión de candidatos para flag=%d\n", flag); */
/*     printf("Secuencia original para el coset %d flag=%d: ", ctx->coset_idx, flag); */
/*     for (size_t i = 0; i < ctx->N; i++) { */
/*         printf("%d ", sequence[i]); */
/*     } */
/*     printf("\n");*\/ */
/*     int compresion_count =(int) ( ctx->N / dimension_candidate_pairs ); */
/*     for (size_t j = 0; j < dimension_candidate_pairs; j++) { */
/*         copia[j] = sequence[j]; */
/*         bound_sequence[j] = bound_bit_sequence[j]; */
/*     } */
/*     for (size_t i=0; i<dimension_candidate_pairs; i++) { */
/*         for (size_t j=1; j<compresion_count; j++) { */
/*             copia[i] += sequence[i + j * dimension_candidate_pairs]; */
/*             bound_sequence[i] += bound_bit_sequence[i + j * dimension_candidate_pairs]; */
/*         } */
/*     } */
/*     /\* */
/*     printf("Vector comprimido para flag=%d: ", flag); */
/*     for (size_t i = 0; i < dimension_candidate_pairs; i++) { */
/*         printf("%d ", copia[i]); */
/*     } */
/*     printf("\n"); */
/*     printf("Bound comprimido para flag=%d: ", flag); */
/*     for (size_t i = 0; i < dimension_candidate_pairs; i++) { */
/*         printf("%d ", bound_sequence[i]); */
/*     } */
/*     printf("\n"); */
/*     *\/ */
    
/*     // Búsqueda binaria aprovechando orden lexicográfico */
/*     // El array está ordenado: pairs_temp[pos][j] >= pairs_temp[pos+1][j] */
/*     bool matches_a = false; */
/*     size_t left = 0, right = num_candidate_pairs - 1; */
/*     int total_iterations = 0; */
    
/*     for (size_t pos = 0; (pos < dimension_candidate_pairs) && (left <= right); pos++)  */
/*     { */
/*         // Verificar primero si right ya cumple la condición */
/*         total_iterations = 0; // Reiniciar contador para la búsqueda del right */
/*         if (copia[pos] > pairs_temp[right][pos]) { */
/*             // Si no cumple, hacer búsqueda binaria para encontrar el mínimo right donde copia[pos] <= pairs_temp[right][pos] */
/*             int l = left, r = right; */
/*             int new_right = -1; */
/*             while (l <= r) { */
/*                 total_iterations++; */
/*                 if (total_iterations > 100) { */
/*                     return true;  // Exceso de iteraciones, retornar true como salvaguarda */
/*                 } */
/*                 int mid = (l + r) >> 1; */
/*                 if (copia[pos] <= pairs_temp[mid][pos]) { */
/*                     new_right = mid; */
/*                     r = mid - 1;  // Buscar más a la izquierda para encontrar el mínimo */
/*                 } else { */
/*                     l = mid + 1; */
/*                 } */
/*             } */
/*             // Si encontramos un right válido, actualizar; si no quedamos fuera del rango */
/*             if (new_right != -1) { */
/*                 right = new_right; */
/*             } else { */
/*                 // No hay ningún elemento que cumpla la condición, terminar este pos */
/*                 return true; */
/*             } */
/*         } */
        
/*         // Verificar primero si left ya cumple la condición */
/*         total_iterations = 0; // Reiniciar contador para la búsqueda del left */
/*         if (copia[pos] + bound_sequence[pos] < pairs_temp[left][pos]) { */
/*             // Si no cumple, hacer búsqueda binaria para encontrar el máximo left donde copia[pos] + bound_sequence[pos] >= pairs_temp[left][pos] */
/*             int l = left, r = right; */
/*             int new_left = -1; */
/*             while (l <= r) { */
/*                 total_iterations++; */
/*                 if (total_iterations > 100) { */
/*                     return true;  // Exceso de iteraciones, retornar true como salvaguarda */
/*                 } */
/*                 int mid = (l + (r - l + 1)) >> 1;  // Redondear hacia arriba para búsqueda del máximo */
/*                 if (copia[pos] + bound_sequence[pos] >= pairs_temp[mid][pos]) { */
/*                     new_left = mid; */
/*                     l = mid + 1;  // Buscar más a la derecha para encontrar el máximo */
/*                 } else { */
/*                     r = mid - 1; */
/*                 } */
/*             } */
/*             // Si encontramos un left válido, actualizar; si no quedamos fuera del rango */
/*             if (new_left != -1 && new_left >= left) { */
/*                 left = new_left; */
/*             } else { */
/*                 // No hay ningún elemento que cumpla la condición, terminar este pos */
/*                 return true; */
/*             } */
/*         } */
/*         /\* */
/*         for (size_t j = 0; (j < dimension_candidate_pairs) && matches_a; j++) { */
/*             if (copia[j] > pairs_temp[pos][j] || bound_sequence[j] + copia[j] < pairs_temp[pos][j]) { */
/*                 /\* */
/*                 printf("No coincide con el par %zu para flag=%d\n", pos, flag); */
/*                 printf("Par candidato: "); */
/*                 for (size_t k = 0; k < dimension_candidate_pairs; k++) { */
/*                     printf("%d ", pairs_temp[pos][k]); */
/*                 } */
/*                 printf("\n"); */
/*                 matches_a = false; */
/*             }*\/ */
/*     } */
/*     /\*printf("¿Coincide con algún par candidato para flag=%d? %s\n", flag, matches_a ? "Sí" : "No");*\/ */
/*     return left<=right; // No se encontró ningún par que coincida */
/* } */


/* bool is_compression_of_candidates(const DFSContext *ctx, const int *sequence, int flag) */
/* { */
/*     int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2; */
/*     int num_candidate_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2; */
/*     int copia[ctx->N]; */
/*     int **pairs_temp = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2; */
/*     /\*printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"); */
/*     printf("Secuencia original para el coset %d flag=%d: ", ctx->coset_idx, flag); */
/*     for (size_t i = 0; i < ctx->N; i++) { */
/*         printf("%d ", sequence[i]); */
/*     } */
/*     printf("\n");*\/ */
/*     for (size_t j = 0; j < dimension_candidate_pairs; j++) { */
/*         copia[j] = sequence[j]; */
/*     } */
/*     int compresion_count =(int) ( ctx->N / dimension_candidate_pairs ); */
/*     for (size_t i=0; i<dimension_candidate_pairs; i++) { */
/*         for (size_t j=1; j<compresion_count; j++) { */
/*             copia[i] += sequence[i + j * dimension_candidate_pairs]; */
/*         } */
/*     } */
/*     /\*printf("imprimiendo la copia:\n"); */
/*     for (size_t j = 0; j < dimension_candidate_pairs; j++) { */
/*         printf("%d ", copia[j]); */
/*     } */
/*     printf("\n");*\/ */
    
/*     // Búsqueda binaria para encontrar copia en pairs_temp (ordenado lexicográficamente) */
/*     bool matches_a = false; */
/*     int left = 0, right = (int)num_candidate_pairs - 1; */
/*     int iterations = 0; */
/*     int max_iterations = 100;  // Failsafe contra bucles infinitos */
    
/*     while (left <= right && !matches_a && iterations < max_iterations) { */
/*         iterations++; */
/*         int mid = (left + right) >> 1; */
/*         if (mid > right || mid < left) { */
/*             printf("Error: índice medio fuera de rango (mid=%d, left=%d, right=%d)\n", mid, left, right); */
/*             break; // Evitar desbordamiento */
/*         } */
        
/*         // Verificar si copia es exactamente igual a pairs_temp[mid] */
/*         bool es_igual = true; */
/*         for (size_t j = 0; j < dimension_candidate_pairs && es_igual; j++) { */
/*             if (copia[j] != pairs_temp[mid][j]) { */
/*                 es_igual = false; */
/*             } */
/*         } */
        
/*         if (es_igual) { */
/*             matches_a = true; */
/*             break; */
/*         } */
        
/*         // Comparación lexicográfica para decidir dirección de búsqueda */
/*         int cmp = 0;  // -1: copia<mid, 0: copia==mid, 1: copia>mid */
/*         for (size_t j = 0; j < dimension_candidate_pairs; j++) { */
/*             if (copia[j] < pairs_temp[mid][j]) { */
/*                 cmp = -1; */
/*                 break; */
/*             } */
/*             if (copia[j] > pairs_temp[mid][j]) { */
/*                 cmp = 1; */
/*                 break; */
/*             } */
/*         } */
        
/*         if (cmp > 0) { */
/*             // copia es lexicográficamente mayor, buscar en la izquierda (valores mayores) */
/*             right = mid - 1; */
/*         } else if (cmp < 0) { */
/*             // copia es lexicográficamente menor, buscar en la derecha (valores menores) */
/*             left = mid + 1; */
/*         } else { */
/*             // cmp == 0: caso imposible - arrays iguales pero es_igual fue false */
/*             // Esto indica un error lógico, terminar búsqueda */
/*             return true; */
/*         } */
/*     } */
/*     /\*printf("Vector comprimido para flag=%d: ", flag); */
/*     for (size_t i = 0; i < dimension_candidate_pairs; i++) { */
/*         printf("%d ", copia[i]); */
/*     } */
/*     printf("\n"); */
/*     printf("¿Coincide con algún par candidato para flag=%d? %s\n", flag, matches_a ? "Sí" : "No"); */
/*     printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");*\/ */
/*     return matches_a || iterations >= max_iterations; // No se encontró ningún par que coincida */
/* } */

/* bool check_bound(const DFSContext *ctx) { */
/*     int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N); */
/*     int *temp0 = malloc(ctx->cl->len * sizeof(int)); */
/*     for (size_t j = ctx->coset_idx; j < ctx->cl->len; j++) { */
/*         temp0[j] = 1; */
/*     } */
/*     int *temp1 = generate_vector_for_combination(ctx->cl, temp0, ctx->N); */
/*     bool result = is_less_than_candidates(ctx, vector_bits, temp1, 1) && is_less_than_candidates(ctx, vector_bits, temp1, 0); */
    
/*     /\* This is for using the compression check*\/ */
/*     /\* */
/*     if  ( !is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_a, temp1)) */
/*     { */
/*         if (!is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_b, temp1)) { */
/*             result = false; */
/*         } */
/*     } */
/*     *\/ */
/*     free(vector_bits); */
/*     free(temp1); */
/*     free(temp0); */
    
/*     if (!result) { */
/*         return false; */
/*     } */
/*     int max_diff = -1; */
/*     for (size_t j = 1; j < ctx->N; j++) { */
/*         int bound_remaining_j = 0; */
/*         size_t start_idx =  ctx->coset_idx ; */
/*         for (size_t i = start_idx; i < ctx->num_cosets; i++) { */
/*             bound_remaining_j += ctx->psd_matrix[i][j]; */
/*         } */
/*         int diff = ctx->current_psd[j] - bound_remaining_j; */
/*         if (diff > max_diff) { */
/*             max_diff = diff; */
/*         } */
/*     } */
/*     return max_diff <= (int)ctx->threshold; */
/* } */

/* bool is_valid_combination(const DFSContext *ctx) { */
/*     int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N); */
    
/*     int *temp= malloc(ctx->N * sizeof(int)); */
/*     for (size_t j = 0; j < ctx->N; j++) { */
/*         temp[j] = vector_bits[j]; */
/*     } */
/*     free(vector_bits); */
    
/*     bool result = false; */
/*     /\* This does not uses the candidates, now we are going to use them*\/ */
/*     //bool is_compressed = is_compression(ctx->N, ctx->p, temp, ctx->compression_a); */
/*     //if (is_compressed ||  is_compression(ctx->N, ctx->p, temp, ctx->compression_b)){ */
/*     if (is_compression_of_candidates(ctx,temp,0) && is_compression_of_candidates(ctx,temp,1)) */
/*     { */
/*         int max_psd = -1; */
/*         for (size_t j = 1; j < ctx->N; j++) { */
/*             if (ctx->current_psd[j] > max_psd) { */
/*                 max_psd = ctx->current_psd[j]; */
/*             } */
/*         } */
/*         result = (max_psd <= (int)ctx->threshold); */
/*     } */
    
/*     free(temp); */
/*     return result; */
/* } */

/* void dfs_explore_combinations( DFSContext *ctx)  */
/* { */
/*     // Caso base: hemos asignado todos los cosets */
/*     if (ctx->coset_idx == ctx->cl->len) { */
        
        
/*         // PODA: Si cumple la condición, guardar */
/*         if (is_valid_combination(ctx)) { */
/*             ctx->matches_found++; */
/*             int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);             */
/*             for (size_t j = 0; j < ctx->N; j++) { */
/*                 fprintf(ctx->f_comb, "%u", vector_bits[j]); */
/*             } */
/*             fprintf(ctx->f_comb, "\n"); */
            
/*             // Guardar PSD y transformación */
/*             for (size_t j = 1; 2 * j < ctx->N; j++) { */
/*                 int psd_val = (int)rint(pow(cabs(ctx->current_dft[j]), 2)); */
/*                 fprintf(ctx->f_psd, "%d%s", psd_val, (j + 1 < ctx->N) ? " " : ""); */
                
/*                 int transformed = (int)rint(ctx->threshold - psd_val); */
/*                 fprintf(ctx->f_cte, "%d%s", transformed, (j + 1 < ctx->N) ? " " : ""); */
/*             } */
/*             fprintf(ctx->f_psd, "\n"); */
/*             fprintf(ctx->f_cte, "\n"); */
            
/*             free(vector_bits); */
/*         } */
/*         return; */
/*     } */
/*     for (int alfabet_val = 0; alfabet_val < ALFABET_SIZE; alfabet_val++) { */
/*         ctx->current_combination[ctx->coset_idx] = alfabet_val; */
/*         for (size_t j = 1; j < ctx->N; j++) { */
/*            ctx->current_dft[j] = ctx->current_dft[j] + alfabet_val * ctx->dft_matrix[ctx->coset_idx][j]; */
/*             ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2)); */
/*         } */
/*         ctx->coset_idx++; */
/*         if (check_bound(ctx)) { */
/*             dfs_explore_combinations(ctx); */
/*         } */
/*         ctx->coset_idx--; */
/*         ctx->current_combination[ctx->coset_idx] = 0; */
/*         for (size_t j = 1; j < ctx->N; j++) { */
/*             ctx->current_dft[j] = ctx->current_dft[j] - alfabet_val * ctx->dft_matrix[ctx->coset_idx][j]; */
/*             ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2)); */
/*         } */
/*     } */
/* } */


/* /\* Versión DFS del procesamiento *\/ */
/* void process_and_filter_vectors_dfs(CosetList *cl, size_t N, int p, int q) { */
/*     if (!cl || cl->len == 0) return; */

/*     double threshold = ((double)N + 1.0) / 2.0; */
    
/*     // Asignar matrices dinámicamente (OPCIÓN 1: Recomendada) */
/*     double complex **dft_matrix = malloc(cl->len * sizeof(double complex *)); */
/*     int **psd_matrix = malloc(cl->len * sizeof(int *)); */
    
/*     for (size_t i = 0; i < cl->len; i++) { */
/*         dft_matrix[i] = malloc(N * sizeof(double complex)); */
/*         psd_matrix[i] = malloc(N * sizeof(int)); */
/*     } */
    
/*     double complex *time_domain = malloc(N * sizeof(double complex)); */
/*     double complex *freq_domain = malloc(N * sizeof(double complex)); */
/*     int *current_psd = calloc(N, sizeof(int)); // inicializar a 0 */
/*     int *bound_psd = calloc(N, sizeof(int)); // inicializar a 0 */
/*     int *compression_a = malloc(p * sizeof(int));  */
/*     int *compression_b = malloc(p * sizeof(int)); */
/*     double complex *current_dft = calloc(N, sizeof(double complex)); // inicializar a 0 */
/*     int *combination = calloc(cl->len, sizeof(int)); // inicializar a 0 */
/*     legendre_sequence(p, q, compression_a, 1); */
/*     legendre_sequence(p, q, compression_b, -1); */
/*     //printf("El valor de cl->len es: %zu\n", cl->len); */
/*     // Calcular DFT y PSD para cada coset */
/*     for (size_t i = 0; i < cl->len; i++) { */
/*         combination[i] = 1; */
/*         int *vector_bits = generate_vector_for_combination(cl, combination, N); */
/*         //printf("\n Vector generado para el coset %zu\n", i); */
/*         binary_to_complex(vector_bits, time_domain, N); */
/*         dft(time_domain, freq_domain, N); */
        
/*         for (size_t j = 0; j < N; j++) { */
/*             dft_matrix[i][j] = freq_domain[j]; */
/*             psd_matrix[i][j] = (int)rint(pow(cabs(freq_domain[j]), 2)); */
/*             bound_psd[j] += psd_matrix[i][j]; */
/*         } */
/*         free(vector_bits); */
/*         combination[i] = 0; */
/*     } */
    
/*     free(time_domain); */
/*     free(freq_domain); */
    
/*     // Abrir archivos */
/*     FILE *f_comb = fopen("combinations.txt", "a"); */
/*     FILE *f_psd = fopen("dft.txt", "a"); */
/*     FILE *f_cte = fopen("cte-dft.txt", "a"); */
/*     size_t rows1, cols1; */
/*     int **candidate1 = read_pairs_file("pairs_5", &rows1, &cols1); */
/*     size_t rows2, cols2; */
/*     int **candidate2 = read_pairs_file("pairs_9", &rows2, &cols2); */
/*     if (!f_comb || !f_psd || !f_cte || !candidate1 || !candidate2) { */
/*         printf("ERROR: No se pudieron abrir los archivos o leer los pares.\n"); */
/*         if (f_comb) fclose(f_comb); */
/*         if (f_psd) fclose(f_psd); */
/*         if (f_cte) fclose(f_cte); */
/*         if (candidate1) free_pairs(candidate1, rows1); */
/*         if (candidate2) free_pairs(candidate2, rows2); */
        
/*         // Limpiar memoria antes de salir */
/*         for (size_t i = 0; i < cl->len; i++) { */
/*             free(dft_matrix[i]); */
/*             free(psd_matrix[i]); */
/*         } */
/*         free(dft_matrix); */
/*         free(psd_matrix); */
/*         free(current_dft); */
/*         free(current_psd); */
/*         free(bound_psd); */
/*         free(combination); */
/*         free(compression_a); */
/*         free(compression_b); */
/*         return; */
/*     } */
/*     // Preparar contexto */
/*     DFSContext ctx = { */
/*         .f_comb = f_comb, */
/*         .f_psd = f_psd, */
/*         .f_cte = f_cte, */
/*         .matches_found = 0, */
/*         .threshold = threshold, */
/*         .N = N, */
/*         .p = p, */
/*         .q = q, */
/*         .num_cosets = cl->len, */
/*         .dft_matrix = dft_matrix, */
/*         .psd_matrix = psd_matrix, */
/*         .current_dft = current_dft, */
/*         .current_psd = current_psd, */
/*         .bound_psd = bound_psd, */
/*         .compression_a = compression_a, */
/*         .compression_b = compression_b, */
/*         .cl = cl,  */
/*         .current_combination = combination, */
/*         .coset_idx = 0, */
/*         .candidate_pairs1 = candidate1, */
/*         .num_candidate_pairs1 = rows1, */
/*         .dimension_candidate_pairs1 = cols1, */
/*         .candidate_pairs2 = candidate2, */
/*         .num_candidate_pairs2 = rows2, */
/*         .dimension_candidate_pairs2 = cols2 */
/*     }; */

/*     printf("=== EXPLORACIÓN DFS (N=%zu, cosets=%zu) ===\n", N, cl->len); */
/*     printf("Condicion: Max(PSD) < %.2f\n\n", threshold); */
/*     // Iniciar DFS desde coset 0 */
/*     dfs_explore_combinations(&ctx); */
    
/*     // Limpieza completa de memoria */
/*     for (size_t i = 0; i < cl->len; i++) { */
/*         free(dft_matrix[i]); */
/*         free(psd_matrix[i]); */
/*     } */
/*     free(dft_matrix); */
/*     free(psd_matrix); */
/*     free(current_dft); */
/*     free(current_psd); */
/*     free(bound_psd); */
/*     free(combination); */
/*     free(compression_a); */
/*     free(compression_b); */
/*     if (candidate1) free_pairs(candidate1, rows1); */
/*     if (candidate2) free_pairs(candidate2, rows2); */
/*     fclose(f_comb); */
/*     fclose(f_psd); */
/*     fclose(f_cte); */

/*     printf("\n=== EXPLORACIÓN FINALIZADA ===\n"); */
/*     printf("Vectores encontrados: %zu\n\n", ctx.matches_found); */
/* } */





int main(void) {
    int p = 11;
    int q = 3;
    size_t N = (size_t)(p*q*q);
    /* size_t k; */
    /* int tamanos[N]; */
    /* for (size_t i = 0; i < N; i++) { */
    /*     tamanos[i] = -1; */
    /* } */
    /* for (size_t k = 1; k<N; k++) { */
    /*     if (gcd(N, k) == 1) { */
    /*         CosetList cl = cyclotomic_cosets(k, N); */
    /*         tamanos[k] = cl.len; */
    /*         free_cosetlist(&cl); */
    /*     } */
    /* } */
    /* for (int minimo_size = 2; minimo_size <= N; minimo_size++) { */
    /*     for (int k=1; k< N; k++){ */
    /*         if (tamanos[k]==minimo_size){ */
    /*             printf("k=%d tiene %d cosets\n",k,tamanos[k]); */
    /*             printf("N=%d k=%d\n\n",N,k); */
    /*             CosetList cl = cyclotomic_cosets(k, N); */
    /*             // Usar versión DFS en lugar de la iterativa */
    /*             process_and_filter_vectors_dfs(&cl, N, p, q); */
    /*             free_cosetlist(&cl); */
    /*         } */
    /*     } */
    /* } */
    /* printf("===FIN==="); */
    /* getchar(); */
    return 0;
}
