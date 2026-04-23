#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mymath.h"
#include <complex.h>
#include <math.h>
#include "pairs_reader.h"
#include "cyclotomic_cosets.h"
#define ALFABET_SIZE 2
#define DEPTH 20
// ...existing code...



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

static const char *path_basename(const char *path) {
    if (!path) return "";
    const char *slash = strrchr(path, '/');
    return slash ? (slash + 1) : path;
}

static void path_dirname(const char *path, char *out, size_t out_size) {
    if (!out || out_size == 0) return;
    if (!path) {
        snprintf(out, out_size, ".");
        return;
    }

    const char *slash = strrchr(path, '/');
    if (!slash) {
        snprintf(out, out_size, ".");
        return;
    }

    size_t len = (size_t)(slash - path);
    if (len == 0) {
        snprintf(out, out_size, "/");
        return;
    }

    if (len >= out_size) len = out_size - 1;
    memcpy(out, path, len);
    out[len] = '\0';
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
    int *positions_candidates1;
    int *positions_candidates2;
    double *current_bound;
    double **depth_exploration_bound;
    size_t num_candidate_pairs2;
    size_t dimension_candidate_pairs2;
    int pos_depth;
} DFSContext;


// ...existing code...
static bool exists_dominating_pair(int **pairs_temp, int num_rows, int dim,
                                   const int *compression) {
    if (num_rows == 0) return false;
    
    int idx = binary_search_sorted_pairs(pairs_temp, num_rows, dim, compression, 1);
    
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

// ...existing code...
static bool filter_dominating_positions_linear(int **pairs_temp, int num_rows, int dim,
                                               const int *compression, const int *compression_bound,
                                               int *positions) {
    if (positions[0] == -1) {
        return false;
    }

    int write_count = 0;

    for (int k = 0; positions[k] != -1; k++) {
        int idx = positions[k];
        if (idx < 0 || idx >= num_rows) {
            continue;
        }

        bool in_range = true;
        for (int j = 0; j < dim; j++) {
            if (pairs_temp[idx][j] < compression[j] || pairs_temp[idx][j] > compression_bound[j] + compression[j]) {
                in_range = false;
                break;
            }
        }
        if (in_range) {
            positions[write_count++] = idx;
        }
    }
    positions[write_count] = -1;
    return write_count > 0;
}

// ...existing code...
static bool exists_pair_in_range(int **pairs_temp, int num_rows, int dim,
                                 const int *compression, const int *compresion_bound) {
    if (num_rows == 0) return false;
    
    int idx = binary_search_sorted_pairs(pairs_temp, num_rows, dim, compression, 1);
    
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


bool is_less_than_candidates(DFSContext *ctx, const int *sequence, const int *bound_bit_sequence, int flag)
{
    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    int num_candidate_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int **pairs_temp = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    int *positions = (flag == 1) ? ctx->positions_candidates1 : ctx->positions_candidates2;
    int *compression = malloc(dimension_candidate_pairs * sizeof(int));
    int *compression_bound = malloc(dimension_candidate_pairs * sizeof(int));

    CompressSequence(ctx->N, dimension_candidate_pairs, sequence, compression);
    CompressSequence(ctx->N, dimension_candidate_pairs, bound_bit_sequence, compression_bound);
    for (int j = 0; j < dimension_candidate_pairs; j++) {
        compression_bound[j] += compression[j];
    }
    printf("Compression:      ");
    for (int j = 0; j < dimension_candidate_pairs; j++) {
        printf("%d,", compression[j]);
    }
    printf("\n");
    printf("Compression bound: ");
    for (int j = 0; j < dimension_candidate_pairs; j++) {
        printf("%d,", compression_bound[j]);
    }
    printf("\n");


    bool resultado =  filter_dominating_positions_linear(pairs_temp, num_candidate_pairs,
                                              dimension_candidate_pairs,
                                              compression, compression_bound,
                                              positions);
    free(compression);
    free(compression_bound);
    return resultado;
}


bool is_compression_of_candidates(const DFSContext *ctx, const int *sequence, int flag)
{
    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    int num_candidate_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int **pairs_temp = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    int *compression = malloc(dimension_candidate_pairs * sizeof(int));
    CompressSequence(ctx->N, dimension_candidate_pairs, sequence, compression);
    int pos = binary_search_sorted_pairs(pairs_temp, num_candidate_pairs, dimension_candidate_pairs, compression, 0);
    free(compression);
    return (pos != -1);
}

bool check_bound(DFSContext *ctx) {
    int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    int *temp0 = calloc(ctx->cl->len, sizeof(int));
    int pos_temp=0;
    for (int j = ctx->coset_idx; j < ctx->cl->len; j++) {
        temp0[pos_temp++] = j;
    }
    temp0[pos_temp] = -1;
    int *temp1 = generate_vector_for_combination(ctx->cl, temp0, ctx->N);
    printf("temp1:");
    for (int j = 0; j < ctx->N; j++) {
        printf("%d,", temp1[j]);
    }
    printf("\n");
    bool result = is_less_than_candidates(ctx, vector_bits, temp1, 1) && is_less_than_candidates(ctx, vector_bits, temp1, 0);
    free(vector_bits);
    free(temp1);
    free(temp0);
    if (!result) {
        return false;
    }
    double max_lower_bound = 0.0;
    int pos = 1;
    for (int j = 1; j < ctx->N; j++) {
        double current_abs = cabs(ctx->current_dft[j]);
        double lower_bound_psd = pow(fmax(0.0, current_abs - ctx->current_bound[j]), 2.0);
        if (lower_bound_psd > max_lower_bound) {
            max_lower_bound = lower_bound_psd;
            pos = j;
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

void dfs_explore_combinations(DFSContext *ctx)
{
    /* Caso base: hemos procesado todos los cosets.
       current_combination es [coset_idx1, coset_idx2, ..., -1]
       Verificar si cumple condición. */
    if (ctx->coset_idx == ctx->cl->len) {
        if (is_valid_combination(ctx)) {
            ctx->matches_found++;
            int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
            for (int j = 0; j < ctx->N; j++) {
                fprintf(ctx->f_comb, "%u", vector_bits[j]);
            }
            fprintf(ctx->f_comb, "\n");

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
        double complex *current_dft_backup = malloc(ctx->N * sizeof(double complex));
        int temporal_nuevo[45] = {0,1,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,1};
        int *current_psd_backup = malloc(ctx->N * sizeof(int));
        memcpy(current_dft_backup, ctx->current_dft, ctx->N * sizeof(double complex));
        memcpy(current_psd_backup, ctx->current_psd, ctx->N * sizeof(int));

        if (alfabet_val == 1) {
            ctx->current_combination[ctx->pos_depth] = ctx->coset_idx;
            ctx->pos_depth++;
            ctx->current_combination[ctx->pos_depth] = -1;
            for (int j = 0; j < ctx->N; j++) {
                ctx->current_dft[j] = ctx->current_dft[j] + ctx->dft_matrix[ctx->coset_idx][j];
                ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
            }
        }

        bool is_interesting = (ctx->coset_idx < 45) ? (temporal_nuevo[ctx->coset_idx] == alfabet_val) : true;

        int *positions_backup1 = malloc((ctx->num_candidate_pairs1 + 1) * sizeof(int));
        int *positions_backup2 = malloc((ctx->num_candidate_pairs2 + 1) * sizeof(int));
        memcpy(positions_backup1, ctx->positions_candidates1, (ctx->num_candidate_pairs1 + 1) * sizeof(int));
        memcpy(positions_backup2, ctx->positions_candidates2, (ctx->num_candidate_pairs2 + 1) * sizeof(int));

        ctx->coset_idx++;
        memcpy(ctx->current_bound, ctx->depth_exploration_bound[ctx->coset_idx], ctx->N * sizeof(double));

        bool pass_bound = check_bound(ctx);

        if (pass_bound && is_interesting) {
            dfs_explore_combinations(ctx);
        }

        if (is_interesting && !pass_bound) {
            printf("Poda en coset_idx=%d, alfabet_val=%d\n", ctx->coset_idx - 1, alfabet_val);
            int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
            printf("Vector actual: ");
            for (int j = 0; j < ctx->N; j++) {
                printf("%u,", vector_bits[j]);
            }
            printf("\n");
            printf("temporal_nuevo:");
            for (int j = 0; j < 45; j++) {
                printf("%d,", temporal_nuevo[j]);
            }
            printf("\n");
            printf("current combination: ");
            for (int j = 0; j <= ctx->pos_depth; j++) {
                printf("%d,", ctx->current_combination[j]);
            }
            printf("\n");

            free(vector_bits);
        }

        memcpy(ctx->positions_candidates1, positions_backup1, (ctx->num_candidate_pairs1 + 1) * sizeof(int));
        memcpy(ctx->positions_candidates2, positions_backup2, (ctx->num_candidate_pairs2 + 1) * sizeof(int));
        memcpy(ctx->current_dft, current_dft_backup, ctx->N * sizeof(double complex));
        memcpy(ctx->current_psd, current_psd_backup, ctx->N * sizeof(int));

        free(positions_backup1);
        free(positions_backup2);
        free(current_dft_backup);
        free(current_psd_backup);

        ctx->coset_idx--;

        if (alfabet_val == 1) {
            ctx->pos_depth--;
            ctx->current_combination[ctx->pos_depth] = -1;
        }
    }
}


// ...existing code (compute_depth_exploration_bounds, free_depth_exploration_bounds)...

static double **compute_depth_exploration_bounds(DFSContext *ctx) {
    int num_cosets = ctx->num_cosets;
    int N = ctx->N;

    double **bounds = calloc(num_cosets + 1, sizeof(double *));
    if (!bounds) return NULL;

    int *combo_arr = calloc(ctx->cl->len, sizeof(int));
    int *bound_arr = calloc(ctx->cl->len, sizeof(int));
    int *pos1_backup = malloc((ctx->num_candidate_pairs1 + 1) * sizeof(int));
    int *pos2_backup = malloc((ctx->num_candidate_pairs2 + 1) * sizeof(int));
    if (!combo_arr || !bound_arr || !pos1_backup || !pos2_backup) {
        free(combo_arr); free(bound_arr); free(pos1_backup); free(pos2_backup);
        free(bounds);
        return NULL;
    }
    memcpy(pos1_backup, ctx->positions_candidates1, (ctx->num_candidate_pairs1 + 1) * sizeof(int));
    memcpy(pos2_backup, ctx->positions_candidates2, (ctx->num_candidate_pairs2 + 1) * sizeof(int));

    for (int i = 0; i <= num_cosets; i++) {
        bounds[i] = calloc(N, sizeof(double));
        if (!bounds[i]) {
            for (int k = 0; k < i; k++) free(bounds[k]);
            free(bounds);
            free(combo_arr); free(bound_arr); free(pos1_backup); free(pos2_backup);
            return NULL;
        }
    }

    int direct_limit = (DEPTH - 1 < num_cosets) ? (DEPTH - 1) : num_cosets;
    for (int i = 1; i <= direct_limit; i++) {
        for (int j = 0; j < N; j++) {
            double acc = 0.0;
            for (int k = i; k < num_cosets; k++) {
                acc += cabs(ctx->dft_matrix[k][j]);
            }
            bounds[i][j] = acc;
        }
    }

    int combo_start_i = (DEPTH > 1) ? DEPTH : 1;
    if (combo_start_i <= num_cosets) {
        for (int i = combo_start_i; i <= num_cosets; i++) {
            int d = num_cosets - i;
            if (d <= 0) {
                continue;
            }

            if (d > DEPTH) {
                for (int j = 0; j < N; j++) {
                    double acc = 0.0;
                    for (int k = i; k < num_cosets; k++) {
                        acc += cabs(ctx->dft_matrix[k][j]);
                    }
                    bounds[i][j] = acc;
                }
                continue;
            }

            int start_coset = i;

            memset(bound_arr, 0, ctx->cl->len * sizeof(int));
            for (int k = 0; k < start_coset; k++) bound_arr[k] = k;
            bound_arr[start_coset] = -1;
            int *bound_seq = generate_vector_for_combination(ctx->cl, bound_arr, N);

            int num_combinations = 1 << d;
            for (int combo = 0; combo < num_combinations; combo++) {
                memset(combo_arr, 0, ctx->cl->len * sizeof(int));
                int pos = 0;
                for (int k = 0; k < d; k++) {
                    if ((combo >> k) & 1) {
                        combo_arr[pos++] = start_coset + k;
                    }
                }
                combo_arr[pos] = -1;
                int *seq = generate_vector_for_combination(ctx->cl, combo_arr, N);
                memcpy(ctx->positions_candidates1, pos1_backup, (ctx->num_candidate_pairs1 + 1) * sizeof(int));
                memcpy(ctx->positions_candidates2, pos2_backup, (ctx->num_candidate_pairs2 + 1) * sizeof(int));
                
                bool psd_below_threshold = true;
                for (int j = 1; j < N; j++) {
                    double complex val = 0.0;
                    for (int k = 0; k < d; k++) {
                        int bit = (combo >> k) & 1;
                        val += bit * ctx->dft_matrix[start_coset + k][j];
                    }
                    double psd_val = pow(cabs(val), 2.0);
                    if (psd_val > ctx->threshold) {
                        psd_below_threshold = false;
                        break;
                    }
                }
                bool passes = psd_below_threshold &&
                              is_less_than_candidates(ctx, seq, bound_seq, 1) &&
                              is_less_than_candidates(ctx, seq, bound_seq, 0);

                if (passes) {
                    for (int j = 0; j < N; j++) {
                        double complex val = 0.0;
                        for (int k = 0; k < d; k++) {
                            int bit = (combo >> k) & 1;
                            val += bit * ctx->dft_matrix[start_coset + k][j];
                        }
                        double abs_val = cabs(val);
                        if (abs_val > bounds[i][j]) {
                            bounds[i][j] = abs_val;
                        }
                    }
                }
                free(seq);
            }
            free(bound_seq);
        }
    }

    memcpy(ctx->positions_candidates1, pos1_backup, (ctx->num_candidate_pairs1 + 1) * sizeof(int));
    memcpy(ctx->positions_candidates2, pos2_backup, (ctx->num_candidate_pairs2 + 1) * sizeof(int));
    free(combo_arr);
    free(bound_arr);
    free(pos1_backup);
    free(pos2_backup);
    return bounds;
}

static void free_depth_exploration_bounds(double **bounds, int num_cosets) {
    if (!bounds) return;
    for (int i = 0; i <= num_cosets; i++) {
        free(bounds[i]);
    }
    free(bounds);
}

/* Versión DFS del procesamiento */
void process_and_filter_vectors_dfs(CosetList *cl, int N, int p, int q,
                                    const char *pairs_file1, const char *pairs_file2) {
    if (!cl || cl->len == 0) return;

    double threshold = ((double)N + 1.0) / 2.0;

    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    int **psd_matrix = malloc(cl->len * sizeof(int *));

    for (int i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        psd_matrix[i] = malloc(N * sizeof(int));
    }

    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    int *current_psd = calloc(N, sizeof(int));
    int *bound_psd = calloc(N, sizeof(int));
    int *compression_a = malloc(p * sizeof(int));
    int *compression_b = malloc(p * sizeof(int));
    double complex *current_dft = calloc(N, sizeof(double complex));
    int *combination = malloc(cl->len * sizeof(int));
    double *current_bound = malloc(N * sizeof(double));

    if (!current_bound || !combination) {
        printf("ERROR: No se pudo asignar memoria\n");
        free(current_dft);
        free(combination);
        free(compression_a);
        free(compression_b);
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
        }
        free(dft_matrix);
        free(psd_matrix);
        return;
    }

    /* Inicializar combination como lista vacía terminada en -1 */
    for (int i = 0; i < cl->len; i++) {
        combination[i] = -1;
    }

    legendre_sequence(p, q, compression_a, 1);
    legendre_sequence(p, q, compression_b, -1);

    /* Calcular DFT y PSD para cada coset usando buffer temporal */
    int *single_coset_comb = malloc(cl->len * sizeof(int));
    for (int i = 0; i < cl->len; i++) single_coset_comb[i] = -1;

    for (int i = 0; i < cl->len; i++) {
        single_coset_comb[0] = i;
        single_coset_comb[1] = -1;

        int *vector_bits = generate_vector_for_combination(cl, single_coset_comb, N);
        binary_to_complex(vector_bits, time_domain, N);
        dft(time_domain, freq_domain, N);

        for (int j = 0; j < N; j++) {
            dft_matrix[i][j] = freq_domain[j];
            psd_matrix[i][j] = (int)rint(pow(cabs(freq_domain[j]), 2));
            bound_psd[j] += psd_matrix[i][j];
        }
        free(vector_bits);
    }
    free(single_coset_comb);

    free(time_domain);
    free(freq_domain);

    /* Asegurar estado inicial DFS vacío */
    for (int i = 0; i < cl->len; i++) {
        combination[i] = -1;
    }

    // Inicializar current_bound con la suma de valores absolutos de la DFT para cada frecuencia
    for (int j = 0; j < N; j++) {
        current_bound[j] = 0.0;
        for (int i = 0; i < cl->len; i++) {
            current_bound[j] += cabs(dft_matrix[i][j]);
        }
    }

    // Abrir archivos de salida
    char comb_filename[512];
    char psd_filename[512];
    char cte_filename[512];
    char output_dir[512];
    const char *pairs_name1 = path_basename(pairs_file1);
    const char *pairs_name2 = path_basename(pairs_file2);
    path_dirname(pairs_file1, output_dir, sizeof(output_dir));

    snprintf(comb_filename, sizeof(comb_filename), "%s/%s_%s_combinations.txt", output_dir, pairs_name1, pairs_name2);
    snprintf(psd_filename, sizeof(psd_filename), "%s/%s_%s_dft.txt", output_dir, pairs_name1, pairs_name2);
    snprintf(cte_filename, sizeof(cte_filename), "%s/%s_%s_cte-dft.txt", output_dir, pairs_name1, pairs_name2);

    FILE *f_comb = fopen(comb_filename, "w");
    FILE *f_psd = fopen(psd_filename, "w");
    FILE *f_cte = fopen(cte_filename, "w");
    size_t rows1, cols1;
    int **candidate1 = read_pairs_file(pairs_file1, &rows1, &cols1);
    size_t rows2, cols2;
    int **candidate2 = read_pairs_file(pairs_file2, &rows2, &cols2);

    int *positions_candidates1 = malloc((rows1 + 1) * sizeof(int));
    int *positions_candidates2 = malloc((rows2 + 1) * sizeof(int));
    if (positions_candidates1) {
        for (size_t i = 0; i < rows1; i++) {
            positions_candidates1[i] = (int)i;
        }
        positions_candidates1[rows1] = -1;
    }
    if (positions_candidates2) {
        for (size_t i = 0; i < rows2; i++) {
            positions_candidates2[i] = (int)i;
        }
        positions_candidates2[rows2] = -1;
    }

    if (!f_comb || !f_psd || !f_cte || !candidate1 || !candidate2 ||
        !positions_candidates1 || !positions_candidates2) {
        printf("ERROR: No se pudieron abrir los archivos o leer los pares de: %s y %s.\n", pairs_file1, pairs_file2);
        if (!f_comb) perror("fopen comb_filename");
        if (!f_psd) perror("fopen psd_filename");
        if (!f_cte) perror("fopen cte_filename");
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);
        free(positions_candidates1);
        free(positions_candidates2);
        free(current_bound);

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
        .positions_candidates1 = positions_candidates1,
        .num_candidate_pairs1 = rows1,
        .dimension_candidate_pairs1 = cols1,
        .candidate_pairs2 = candidate2,
        .positions_candidates2 = positions_candidates2,
        .num_candidate_pairs2 = rows2,
        .dimension_candidate_pairs2 = cols2,
        .current_bound = current_bound,
        .depth_exploration_bound = NULL,
        .pos_depth = 0
    };

    ctx.depth_exploration_bound = compute_depth_exploration_bounds(&ctx);
    if (ctx.depth_exploration_bound == NULL) {
        printf("ERROR: No se pudieron calcular los bounds de exploración por profundidad.\n");
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
        free(positions_candidates1);
        free(positions_candidates2);
        free(current_bound);
        fclose(f_comb);
        fclose(f_psd);
        fclose(f_cte);
        return;
    }

    printf("\n=== EXPLORACIÓN DFS (N=%d, cosets=%d) ===\n", N, cl->len);
    printf("Condicion: Max(PSD) < %.2f\n\n", threshold);

    /* Asegurar que combination empieza vacío antes de DFS */
    for (int i = 0; i < cl->len; i++) combination[i] = -1;
    ctx.coset_idx = 0;
    ctx.pos_depth = 0;

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
    free(positions_candidates1);
    free(positions_candidates2);
    free(current_bound);
    if (candidate1) free_pairs(candidate1, rows1);
    if (candidate2) free_pairs(candidate2, rows2);
    fclose(f_comb);
    fclose(f_psd);
    fclose(f_cte);
    free_depth_exploration_bounds(ctx.depth_exploration_bound, ctx.num_cosets);
    printf("\n=== EXPLORACIÓN FINALIZADA ===\n");
    printf("Vectores encontrados: %zu\n\n", ctx.matches_found);
}


int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <pairs_file1> <pairs_file2> [p] [q] [k]\n", argv[0]);
        return 1;
    }

    const char *pairs_file1 = argv[1];
    const char *pairs_file2 = argv[2];

    int p = (argc >= 4) ? atoi(argv[3]) : 5;
    int q = (argc >= 5) ? atoi(argv[4]) : 3;
    int k_arg = (argc >= 6) ? atoi(argv[5]) : -1;

    if (p <= 0 || q <= 0) {
        fprintf(stderr, "Error: p y q deben ser positivos.\n");
        return 1;
    }

    int N = p * q * q;
    printf("N=%d, p=%d, q=%d\n", N, p, q);

    if (k_arg != -1) {
        /* Ejecutar con un k específico */
        if (k_arg <= 0 || k_arg >= N || gcd(N, k_arg) != 1) {
            fprintf(stderr, "Error: k=%d invalido. Debe cumplir 1 <= k < N y gcd(N,k)=1.\n", k_arg);
            return 1;
        }
        CosetList cl = cyclotomic_cosets(k_arg, N);
        printf("k=%d tiene %d cosets\n", k_arg, cl.len);
        process_and_filter_vectors_dfs(&cl, N, p, q, pairs_file1, pairs_file2);
        free_cosetlist(&cl);
    } else {
        /* Explorar todos los k coprimos con N, ordenados por número de cosets */
        int *tamanos = malloc(N * sizeof(int));
        for (int i = 0; i < N; i++) {
            tamanos[i] = -1;
        }
        for (int k = 1; k < N; k++) {
            if (gcd(N, k) == 1) {
                CosetList cl = cyclotomic_cosets(k, N);
                tamanos[k] = cl.len;
                free_cosetlist(&cl);
            }
        }
        for (int minimo_size = 1; minimo_size <= N; minimo_size++) {
            for (int k = 1; k < N; k++) {
                if (tamanos[k] == minimo_size) {
                    printf("k=%d tiene %d cosets\n", k, tamanos[k]);
                    printf("N=%d k=%d\n\n", N, k);
                    CosetList cl = cyclotomic_cosets(k, N);
                    process_and_filter_vectors_dfs(&cl, N, p, q, pairs_file1, pairs_file2);
                    free_cosetlist(&cl);
                }
            }
        }
        free(tamanos);
    }

    printf("===FIN===\n");
    return 0;
}