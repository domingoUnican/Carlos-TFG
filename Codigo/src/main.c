#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mymath.h"
#include <complex.h>
#include <math.h>
#include "pairs_reader.h"
#include "cyclotomic_cosets.h"
#define ALFABET_SIZE 2
// New version of gcd that works for negative numbers as well



int gcd(int a, int b) {
    a = abs(a);
    b = abs(b);
    while (b != 0) {
        int t = a % b;
        a = b;
        b = t;
    }
    return a;
}




typedef struct {
    FILE *f_comb;
    FILE *f_psd;
    FILE *f_cte;
    size_t matches_found;
    double threshold;
    //double constant; constant is the same as threshold.
    double complex **dft_matrix;
    int **psd_matrix;
    int N;
    int num_cosets;
    double complex *current_dft;
    int *current_psd;
    int *bound_psd;
    int *compression_a;
    int *compression_b;
    int p;
    int q;
    int coset_idx;
    CosetList *cl;
    int *current_combination;
    int **candidate_pairs1;
    size_t num_candidate_pairs1;
    size_t dimension_candidate_pairs1;
    int **candidate_pairs2;
    size_t num_candidate_pairs2;
    size_t dimension_candidate_pairs2;
} DFSContext;


// Busca un par en pairs_temp tal que compression[j] < pairs_temp[k][j] para todo j.
// Reutiliza binary_search_sorted_pairs ya que pairs_temp está ordenado lexicográficamente.
static bool exists_dominating_pair(int **pairs_temp, int num_rows, int dim,
                                   const int *compression) {
    if (num_rows == 0) return false;
    
    // Encuentra el índice de inserción (primer elemento lex. >= compression)
    int idx = binary_search_sorted_pairs(pairs_temp, num_rows, dim, compression, 1);
    
    // Verifica si el elemento en idx domina componente a componente
    if (idx != -1 && idx < num_rows) {
        bool all_strictly_greater = true;
        for (int j = 0; j < dim; j++) {
            if (pairs_temp[idx][j] < compression[j]) {
                all_strictly_greater = false;
                break;
            }
        }
        if (all_strictly_greater) return true;
    }
    return false;
}

// Busca un par en pairs_temp en el rango: compression[j] <= pairs_temp[k][j] <= compresion_bound[j] para todo j.
// Reutiliza binary_search_sorted_pairs ya que pairs_temp está ordenado lexicográficamente.
static bool exists_pair_in_range(int **pairs_temp, int num_rows, int dim,
                                 const int *compression, const int *compresion_bound) {
    if (num_rows == 0) return false;
    
    // Encuentra el índice de inserción (primer elemento lex. >= compression)
    int idx = binary_search_sorted_pairs(pairs_temp, num_rows, dim, compression, 1);
    
    // Verifica si el elemento en idx está en el rango componente a componente
    if (idx != -1 && idx < num_rows) {
        bool in_range = true;
        for (int j = 0; j < dim; j++) {
            if (pairs_temp[idx][j] < compression[j] || pairs_temp[idx][j] > compresion_bound[j]) {
                in_range = false;
                break;
            }
        }
        if (in_range) return true;
    }
    return false;
}


bool is_less_than_candidates(const DFSContext *ctx, const int *sequence, const int *bound_bit_sequence, int flag)
{
    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    int num_candidate_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int **pairs_temp = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    int compression[dimension_candidate_pairs];
    int compression_bound[dimension_candidate_pairs];

    CompressSequence(ctx->N, dimension_candidate_pairs, sequence, compression);
    CompressSequence(ctx->N, dimension_candidate_pairs, bound_bit_sequence, compression_bound);
    
    // compression_bound = compression - bound; luego usaremos como upper bound
    for (int j = 0; j < dimension_candidate_pairs; j++) {
        compression_bound[j] += compression[j];
    }

    return exists_pair_in_range(pairs_temp, num_candidate_pairs, dimension_candidate_pairs, compression, compression_bound);
}


bool is_compression_of_candidates(const DFSContext *ctx, const int *sequence, int flag)
{
    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    int num_candidate_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int **pairs_temp = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    int compression[dimension_candidate_pairs];
    CompressSequence(ctx->N, dimension_candidate_pairs, sequence, compression);
    int pos = binary_search_sorted_pairs(pairs_temp, num_candidate_pairs, dimension_candidate_pairs, compression, 0);
    return (pos != -1);
}

bool check_bound(const DFSContext *ctx) {
    int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    //int debug_psd[63]= {0,0,1,0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,0,0,1,0,1,1,0,1,0,1,1,1,0,0,0,1,1,1,0,1,1,1};
    //bool is_interesting = true;
    /*for (int j=0; (j<ctx->N) && is_interesting; j++) {
        is_interesting = is_interesting && (vector_bits[j] <= debug_psd[j]);
    }
    */
    int *temp0 = calloc(ctx->cl->len, sizeof(int));
    for (int j = ctx->coset_idx; j < ctx->cl->len; j++) {
        temp0[j] = 1;
    }
    int *temp1 = generate_vector_for_combination(ctx->cl, temp0, ctx->N);
    bool result = is_less_than_candidates(ctx, vector_bits, temp1, 1) && is_less_than_candidates(ctx, vector_bits, temp1, 0);

    /* This is for using the compression check
    if  ( !is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_a, temp1))
    {
        if (!is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_b, temp1)) {
            result = false;
        }
    }
    
   for (int j = 0; j < ctx->N; j++) {
           if (vector_bits[j] < debug_psd[j]) {
               printf("pos: %d\n\n\n", j);
               if (j < ctx->coset_idx) {
                   printf("Coset idx: %d\n", ctx->coset_idx);
                   return false;
               }
               j = ctx->N; // break
           }
       }
   if (is_interesting && !result) {
       printf("Comb: [");
       for (int j = 0; j < ctx->cl->len; j++) {
           printf("%d", ctx->current_combination[j]);
           if (j + 1 < ctx->cl->len) printf(", ");
       }
       printf("]\n-vec: [");
       for (int j = 0; j < ctx->N; j++) {
           printf("%d", vector_bits[j]);
           if (j + 1 < ctx->N) printf(", ");
       }
       printf("]\n Coxind:%d\n - Result: %s\n", ctx->coset_idx, result ? "true" : "false");
       for (int j = 0; j < ctx->N; j++) {
           if (vector_bits[j] < debug_psd[j]) {
               printf("pos: %d\n\n\n", j);
               j = ctx->N; // break
           }
       }
    } */   
    free(vector_bits);
    free(temp1);
    free(temp0);

    if (!result) {
        return false;
    }
    double max_lower_bound = 0.0;
    for (int j = 1; j < ctx->N; j++) {
        double remaining_abs_sum = 0.0;
        int start_idx =  ctx->coset_idx ;
        for (int i = start_idx; i < ctx->num_cosets; i++) {
            remaining_abs_sum += cabs(ctx->dft_matrix[i][j]);
        }
        double current_abs = cabs(ctx->current_dft[j]);
        double lower_bound_psd = pow(fmax(0.0, current_abs - remaining_abs_sum), 2.0);
        if (lower_bound_psd > max_lower_bound) {
            max_lower_bound = lower_bound_psd;
        }
    }
    return max_lower_bound <= ctx->threshold;
}

bool is_valid_combination(const DFSContext *ctx) {
    int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);


    int *temp= malloc(ctx->N * sizeof(int));
    for (int j = 0; j < ctx->N; j++) {
        temp[j] = vector_bits[j];
    }
    free(vector_bits);

    bool result = false;
    /* This does not uses the candidates, now we are going to use them */
    //bool is_compressed = is_compression(ctx->N, ctx->p, temp, ctx->compression_a);
    //if (is_compressed ||  is_compression(ctx->N, ctx->p, temp, ctx->compression_b)){
    if (is_compression_of_candidates(ctx,temp,0) && is_compression_of_candidates(ctx,temp,1))
    {
        int max_psd = -1;
        for (int j = 1; j < ctx->N; j++) {
            if (ctx->current_psd[j] > max_psd) {
                max_psd = ctx->current_psd[j];
            }
        }
        result = (max_psd <= (int)ctx->threshold);
    }

    free(temp);
    return result;
}

void dfs_explore_combinations( DFSContext *ctx)
{
    // Caso base: hemos asignado todos los cosets
    if (ctx->coset_idx == ctx->cl->len) {


        // PODA: Si cumple la condición, guardar
        if (is_valid_combination(ctx)) {
            ctx->matches_found++;
            int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
            for (int j = 0; j < ctx->N; j++) {
                fprintf(ctx->f_comb, "%u", vector_bits[j]);
            }
            fprintf(ctx->f_comb, "\n");

            // Guardar PSD y transformación
            for (int j = 1; 2 * j < ctx->N + 1; j++) {
                int psd_val = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
                fprintf(ctx->f_psd, "%d%s", psd_val, (j + 1 < ctx->N) ? " " : "");

                int transformed = (int)rint(ctx->threshold - psd_val);
                fprintf(ctx->f_cte, "%d%s", transformed, (j + 1 < ctx->N) ? " " : "");
            }
            fprintf(ctx->f_psd, "\n");
            fprintf(ctx->f_cte, "\n");

            free(vector_bits);
        }
        return;
    }
    for (int alfabet_val = 0; alfabet_val < ALFABET_SIZE; alfabet_val++) {
        ctx->current_combination[ctx->coset_idx] = alfabet_val;
        for (int j = 1; j < ctx->N; j++) {
           ctx->current_dft[j] = ctx->current_dft[j] + alfabet_val * ctx->dft_matrix[ctx->coset_idx][j];
            ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
        }
        ctx->coset_idx++;
        if (check_bound(ctx)) {
            dfs_explore_combinations(ctx);
        }
        ctx->coset_idx--;
        ctx->current_combination[ctx->coset_idx] = 0;
        for (int j = 1; j < ctx->N; j++) {
            ctx->current_dft[j] = ctx->current_dft[j] - alfabet_val * ctx->dft_matrix[ctx->coset_idx][j];
            ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
        }
    }
}


/* Versión DFS del procesamiento */
void process_and_filter_vectors_dfs(CosetList *cl, int N, int p, int q,
                                    const char *pairs_file1, const char *pairs_file2) {
    if (!cl || cl->len == 0) return;

    double threshold = ((double)N + 1.0) / 2.0;

    // Asignar matrices dinámicamente (OPCIÓN 1: Recomendada)
    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    int **psd_matrix = malloc(cl->len * sizeof(int *));

    for (int i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        psd_matrix[i] = malloc(N * sizeof(int));
    }

    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    int *current_psd = calloc(N, sizeof(int)); // inicializar a 0
    int *bound_psd = calloc(N, sizeof(int)); // inicializar a 0
    int *compression_a = malloc(p * sizeof(int));
    int *compression_b = malloc(p * sizeof(int));
    double complex *current_dft = calloc(N, sizeof(double complex)); // inicializar a 0
    int *combination = calloc(cl->len, sizeof(int)); // inicializar a 0
    legendre_sequence(p, q, compression_a, 1);
    legendre_sequence(p, q, compression_b, -1);
    //printf("El valor de cl->len es: %zu\n", cl->len);
    // Calcular DFT y PSD para cada coset
    for (int i = 0; i < cl->len; i++) {
        combination[i] = 1;
        int *vector_bits = generate_vector_for_combination(cl, combination, N);
        //printf("\n Vector generado para el coset %zu\n", i);
        binary_to_complex(vector_bits, time_domain, N);
        dft(time_domain, freq_domain, N);

        for (int j = 0; j < N; j++) {
            dft_matrix[i][j] = freq_domain[j];
            psd_matrix[i][j] = (int)rint(pow(cabs(freq_domain[j]), 2));
            bound_psd[j] += psd_matrix[i][j];
        }
        free(vector_bits);
        combination[i] = 0;
    }

    free(time_domain);
    free(freq_domain);

    // Abrir archivos
    char comb_filename[512];
    char psd_filename[512];
    char cte_filename[512];
    snprintf(comb_filename, sizeof(comb_filename), "%s_%s_combinations.txt", pairs_file1, pairs_file2);
    snprintf(psd_filename, sizeof(psd_filename), "%s_%s_dft.txt", pairs_file1, pairs_file2);
    snprintf(cte_filename, sizeof(cte_filename), "%s_%s_cte-dft.txt", pairs_file1, pairs_file2);

    FILE *f_comb = fopen(comb_filename, "a");
    FILE *f_psd = fopen(psd_filename, "a");
    FILE *f_cte = fopen(cte_filename, "a");
    size_t rows1, cols1;
    int **candidate1 = read_pairs_file(pairs_file1, &rows1, &cols1);
    size_t rows2, cols2;
    int **candidate2 = read_pairs_file(pairs_file2, &rows2, &cols2);
    if (!f_comb || !f_psd || !f_cte || !candidate1 || !candidate2) {
        printf("ERROR: No se pudieron abrir los archivos o leer los pares de: %s y %s.\n", pairs_file1, pairs_file2);
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);

        // Limpiar memoria antes de salir
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
        }
        free(dft_matrix);
        free(psd_matrix);
        free(current_dft);
        free(current_psd);
        free(bound_psd);
        free(combination);
        free(compression_a);
        free(compression_b);
        return;
    }
    // Preparar contexto
    DFSContext ctx = {
        .f_comb = f_comb,
        .f_psd = f_psd,
        .f_cte = f_cte,
        .matches_found = 0,
        .threshold = threshold,
        .N = N,
        .p = p,
        .q = q,
        .num_cosets = cl->len,
        .dft_matrix = dft_matrix,
        .psd_matrix = psd_matrix,
        .current_dft = current_dft,
        .current_psd = current_psd,
        .bound_psd = bound_psd,
        .compression_a = compression_a,
        .compression_b = compression_b,
        .cl = cl,
        .current_combination = combination,
        .coset_idx = 0,
        .candidate_pairs1 = candidate1,
        .num_candidate_pairs1 = rows1,
        .dimension_candidate_pairs1 = cols1,
        .candidate_pairs2 = candidate2,
        .num_candidate_pairs2 = rows2,
        .dimension_candidate_pairs2 = cols2
    };

    printf("=== EXPLORACIÓN DFS (N=%d, cosets=%d) ===\n", N, cl->len);
    printf("Condicion: Max(PSD) < %.2f\n\n", threshold);
    // Iniciar DFS desde coset 0
    dfs_explore_combinations(&ctx);

    // Limpieza completa de memoria
    for (int i = 0; i < cl->len; i++) {
        free(dft_matrix[i]);
        free(psd_matrix[i]);
    }
    free(dft_matrix);
    free(psd_matrix);
    free(current_dft);
    free(current_psd);
    free(bound_psd);
    free(combination);
    free(compression_a);
    free(compression_b);
    if (candidate1) free_pairs(candidate1, rows1);
    if (candidate2) free_pairs(candidate2, rows2);
    fclose(f_comb);
    fclose(f_psd);
    fclose(f_cte);

    printf("\n=== EXPLORACIÓN FINALIZADA ===\n");
    printf("Vectores encontrados: %zu\n\n", ctx.matches_found);
}





int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <pairs_file1> <pairs_file2>\n", argv[0]);
        return 1;
    }

    const char *pairs_file1 = argv[1];
    const char *pairs_file2 = argv[2];

    int p = 3;
    int q = 5;
    int N = (p*q*q);
    int tamanos[N];
    for (int i = 0; i < N; i++) {
        tamanos[i] = -1;
    }
    for (int k = 1; k<N; k++) {
        if (gcd(N, k) == 1) {
            CosetList cl = cyclotomic_cosets(k, N);
            tamanos[k] = cl.len;
            free_cosetlist(&cl);
        }
    }
    for (int minimo_size = 2; minimo_size <= N; minimo_size++) {
        for (int k=1; k< N; k++){
            if (tamanos[k]==minimo_size){
                printf("k=%d tiene %d cosets\n",k,tamanos[k]);
                printf("N=%d k=%d\n\n",N,k);
                CosetList cl = cyclotomic_cosets(k, N);
                // Usar versión DFS en lugar de la iterativa
                process_and_filter_vectors_dfs(&cl, N, p, q, pairs_file1, pairs_file2);
                free_cosetlist(&cl);
            }
        }
    }
    printf("===FIN===");
    //getchar();
    return 0;
}
