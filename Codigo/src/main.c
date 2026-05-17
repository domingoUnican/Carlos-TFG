#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mymath.h"
#include <complex.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include "pairs_reader.h"
#include "cyclotomic_cosets.h"
#define ALFABET_SIZE 2
#define DEPTH 0
#define DIMENSIONS 3
#ifndef USE_CHECK_BOUND
#define USE_CHECK_BOUND 1
#endif



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
    char comb_filename[512];
    char psd_filename[512];
    char cte_filename[512];
    char lp_filename[512];
    size_t matches_found;
    int threshold;
    //double constant; constant is the same as threshold.
    double complex **dft_matrix;
    double complex **dft_backup;
    int N;
    int num_cosets;
    double complex *current_dft;
    int *bound_psd;
    int *current_sequence;
    int *current_psd; 
    int *remaining_bits;
    int *compression_1;
    int *compression_2;
    int *bound_compression_1;
    int *bound_compression_2;
    int *compression_3;
    int *bound_compression_3;
    int *compression_a;
    int p;
    int q;
    int coset_idx;
    CosetList *cl;
    int *current_combination1;
    int *current_combination2;
    int **candidate_pairs1;
    size_t num_candidate_pairs1;
    size_t dimension_candidate_pairs1;
    int **candidate_pairs2;
    double *current_bound;
    double **depth_exploration_bound;
    size_t num_candidate_pairs2;
    size_t dimension_candidate_pairs2;
    int pos_depth_1;
    int pos_depth_2;
    int current_ones1;
    int current_ones2;
    int target_ones;
    int spectrum_size;
    int *suffix_ones1;
    int *suffix_ones2;
    bool *is_less_compression_a;
    bool *is_less_compression_b;
    bool **is_less_candidates_1;
    bool **is_less_candidates_2;
} DFSContext;

static int compare_int_asc(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    return (ia > ib) - (ia < ib);
}

static bool **allocate_bool_matrix(size_t rows, size_t cols) {
    bool **matrix = malloc(rows * sizeof(bool *));
    if (!matrix) return NULL;

    for (size_t i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(bool));
        if (!matrix[i]) {
            for (size_t j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            return NULL;
        }
        for (size_t j = 0; j < cols; j++) {
            matrix[i][j] = true;
        }
    }
    return matrix;
}

static void free_bool_matrix(bool **matrix, size_t rows) {
    if (!matrix) return;
    for (size_t i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void reorder_cosets_by_candidate_pairs(DFSContext *ctx) {
    return ; 
    if (!ctx || !ctx->cl || ctx->cl->len == 0) return;

    int num_cosets = ctx->cl->len;
    if (ctx->dimension_candidate_pairs1 == 0 || ctx->dimension_candidate_pairs2 == 0) return;

    int d1 = ctx->N / (int)ctx->dimension_candidate_pairs1;
    int d2 = ctx->N / (int)ctx->dimension_candidate_pairs2;
    int *weights1 = binomial_coefficients(d1);
    int *weights2 = binomial_coefficients(d2);
    if (!weights1 || !weights2) {
        free(weights1);
        free(weights2);
        return;
    }

    int values_count = (int)(ctx->num_candidate_pairs1 + ctx->num_candidate_pairs2);
    int *values = malloc((size_t)values_count * sizeof(int));
    int *element_weight = malloc((size_t)ctx->N * sizeof(int));
    int *coset_weight = malloc((size_t)num_cosets * sizeof(int));
    int *coset_indices = malloc(num_cosets * sizeof(int));
    if (!values || !element_weight || !coset_weight || !coset_indices) {
        free(values);
        free(element_weight);
        free(coset_weight);
        free(coset_indices);
        free(weights1);
        free(weights2);
        return;
    }

    for (int elem = 0; elem < ctx->N; elem++) {
        int col1 = elem % (int)ctx->dimension_candidate_pairs1;
        int col2 = elem % (int)ctx->dimension_candidate_pairs2;
        int nvals = 0;

        for (size_t i = 0; i < ctx->num_candidate_pairs1; i++) {
            int v = ctx->candidate_pairs1[i][col1];
            if (v >= 0 && v <= d1) {
                values[nvals++] = weights1[v];
            }
        }
        for (size_t i = 0; i < ctx->num_candidate_pairs2; i++) {
            int v = ctx->candidate_pairs2[i][col2];
            if (v >= 0 && v <= d2) {
                values[nvals++] = weights2[v];
            }
        }

        if (nvals == 0) {
            element_weight[elem] = INT_MAX;
        } else {
            qsort(values, (size_t)nvals, sizeof(int), compare_int_asc);
            element_weight[elem] = values[nvals / 2];
        }
    }

    for (int i = 0; i < num_cosets; i++) {
        coset_indices[i] = i;
        int min_w = INT_MAX;
        Coset *c = &ctx->cl->data[i];
        for (int j = 0; j < c->len; j++) {
            int elem = c->data[j];
            if (elem >= 0 && elem < ctx->N && element_weight[elem] < min_w) {
                min_w = element_weight[elem];
            }
        }
        coset_weight[i] = min_w;
    }

    for (int i = 0; i < num_cosets - 1; i++) {
        int best = i;
        for (int j = i + 1; j < num_cosets; j++) {
            int idx_best = coset_indices[best];
            int idx_j = coset_indices[j];
            if (coset_weight[idx_j] < coset_weight[idx_best] ||
                (coset_weight[idx_j] == coset_weight[idx_best] && idx_j < idx_best)) {
                best = j;
            }
        }
        if (best != i) {
            int tmp = coset_indices[i];
            coset_indices[i] = coset_indices[best];
            coset_indices[best] = tmp;
        }
    }

    Coset *reordered_cosets = NULL;
    double complex **reordered_dft = NULL;
    int **reordered_psd = NULL;

    reordered_cosets = malloc((size_t)num_cosets * sizeof(Coset));
    if (!reordered_cosets) {
        free(values);
        free(element_weight);
        free(coset_weight);
        free(coset_indices);
        free(weights1);
        free(weights2);
        return;
    }

    reordered_dft = malloc((size_t)num_cosets * sizeof(double complex *));
    if (!reordered_dft) {
        free(reordered_cosets);
        free(values);
        free(element_weight);
        free(coset_weight);
        free(coset_indices);
        free(weights1);
        free(weights2);
        return;
    }

    reordered_psd = malloc((size_t)num_cosets * sizeof(int *));
    if (!reordered_psd) {
        free(reordered_cosets);
        free(reordered_dft);
        free(values);
        free(element_weight);
        free(coset_weight);
        free(coset_indices);
        free(weights1);
        free(weights2);
        return;
    }

    for (int i = 0; i < num_cosets; i++) {
        int old_idx = coset_indices[i];
        reordered_cosets[i] = ctx->cl->data[old_idx];
        reordered_dft[i] = ctx->dft_matrix[old_idx];
    }

    free(ctx->cl->data);
    free(ctx->dft_matrix);
    ctx->cl->data = reordered_cosets;
    ctx->dft_matrix = reordered_dft;

    if (ctx->cl->positions) {
        for (int i = 0; i < num_cosets; i++) {
            Coset *c = &ctx->cl->data[i];
            for (int j = 0; j < c->len; j++) {
                int elem = c->data[j];
                if (elem >= 0 && elem < ctx->N) {
                    ctx->cl->positions[elem] = i;
                }
            }
        }
    }

    if (ctx->suffix_ones1 && ctx->suffix_ones2) {
        ctx->suffix_ones1[num_cosets] = 0;
        ctx->suffix_ones2[num_cosets] = 0;
        for (int i = num_cosets - 1; i >= 0; i--) {
            ctx->suffix_ones1[i] = ctx->suffix_ones1[i + 1] + ctx->cl->data[i].len;
            ctx->suffix_ones2[i] = ctx->suffix_ones2[i + 1] + ctx->cl->data[i].len;
        }
    }

    free(values);
    free(element_weight);
    free(coset_weight);
    free(coset_indices);
    free(weights1);
    free(weights2);
}

static bool is_less_than_candidates_for_vectors(DFSContext *ctx, const int *sequence, const int *bound_bit_sequence, int flag)
{
    int N = ctx->N;
    int full_dimension_1 = (int)ctx->dimension_candidate_pairs1;
    int full_dimension_2 = (int)ctx->dimension_candidate_pairs2;
    int half_dimension_1 = full_dimension_1 >> 1;
    int half_dimension_2 = full_dimension_2 >> 1;

    if (half_dimension_1 <= 0 || half_dimension_2 <= 0) {
        return false;
    }

    int *compression_1 = calloc((size_t)half_dimension_1, sizeof(int));
    int *compression_bound_1 = calloc((size_t)half_dimension_1, sizeof(int));
    int *compression_2 = calloc((size_t)half_dimension_2, sizeof(int));
    int *compression_bound_2 = calloc((size_t)half_dimension_2, sizeof(int));
    if (!compression_1 || !compression_bound_1 || !compression_2 || !compression_bound_2) {
        free(compression_1);
        free(compression_bound_1);
        free(compression_2);
        free(compression_bound_2);
        return false;
    }

    for (int i = 0; i < N; i++) {
        compression_1[i % half_dimension_1] += sequence[i];
        compression_bound_1[i % half_dimension_1] += bound_bit_sequence[i];
        compression_2[i % half_dimension_2] += sequence[i];
        compression_bound_2[i % half_dimension_2] += bound_bit_sequence[i];
    }

    bool feasible_1 = false;
    int offset_1 = (flag == 0) ? 0 : half_dimension_1;
    int other_offset_1 = (flag == 0) ? half_dimension_1 : 0;
    for (size_t r = 0; r < ctx->num_candidate_pairs1; r++) {
        bool ok = true;
        for (int c = 0; c < half_dimension_1; c++) {
            int v = ctx->candidate_pairs1[r][offset_1 + c];
            int l = compression_1[c];
            int u = l + compression_bound_1[c];
            if (!(l <= v && v <= u)) {
                ok = false;
                break;
            }
        }
        for (int c = 0; c < half_dimension_1 && ok; c++) {
            if (ctx->candidate_pairs1[r][other_offset_1 + c] != 0) {
                ok = false;
            }
        }
        if (ok) {
            feasible_1 = true;
            break;
        }
    }

    bool feasible_2 = false;
    int offset_2 = (flag == 0) ? 0 : half_dimension_2;
    int other_offset_2 = (flag == 0) ? half_dimension_2 : 0;
    for (size_t r = 0; r < ctx->num_candidate_pairs2; r++) {
        bool ok = true;
        for (int c = 0; c < half_dimension_2; c++) {
            int v = ctx->candidate_pairs2[r][offset_2 + c];
            int l = compression_2[c];
            int u = l + compression_bound_2[c];
            if (!(l <= v && v <= u)) {
                ok = false;
                break;
            }
        }
        for (int c = 0; c < half_dimension_2 && ok; c++) {
            if (ctx->candidate_pairs2[r][other_offset_2 + c] != 0) {
                ok = false;
            }
        }
        if (ok) {
            feasible_2 = true;
            break;
        }
    }

    free(compression_1);
    free(compression_bound_1);
    free(compression_2);
    free(compression_bound_2);
    return feasible_1 && feasible_2;
}

bool is_less_than_candidates(DFSContext *ctx, int flag)
{
    if (!ctx || !ctx->compression_3 || !ctx->bound_compression_3 || !ctx->compression_a || ctx->p <= 0) {
        return false;
    }
    if (flag == 1){
        bool ok_legendre = true;
        for (int i = 0; i < 2 * ctx->p; i++) {
            int l = ctx->compression_3[i];
            int u = l + ctx->bound_compression_3[i];
            int ta = ctx->compression_a[i];
            if (!(l <= ta && ta <= u)) {
                ok_legendre = false;
                break;
            }
        }
        ctx->is_less_compression_a[ctx->coset_idx] = ok_legendre;
        if (! ctx->is_less_compression_a[ctx->coset_idx]) {
            return false;
        }
    }
    

    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    size_t num_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int *compression = (flag == 1) ? ctx->compression_1 : ctx->compression_2;
    int *compression_bound = (flag == 1) ? ctx->bound_compression_1 : ctx->bound_compression_2;
    int ** candidate_pairs = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    bool ** is_less_cand = (flag == 1) ? ctx->is_less_candidates_1 : ctx->is_less_candidates_2;
    Coset *selected_coset = (ctx->coset_idx > 0) ? &ctx->cl->data[ctx->coset_idx - 1] : NULL;
    int half_dimension = dimension_candidate_pairs / 2;
    size_t prev_idx = (ctx->coset_idx > 0) ? (size_t)(ctx->coset_idx - 1) : 0;

    bool exists = false;
    for (size_t j = 0; j < num_pairs; j++) {
        bool ok = (ctx->coset_idx > 0) ? is_less_cand[prev_idx][j] : true;
        if (ok && selected_coset && half_dimension > 0) {
            for (int t = 0; t < selected_coset->len && ok; t++) {
                int elem = selected_coset->data[t];
                int c0 = elem % half_dimension;
                int c1 = c0 + half_dimension;

                int l0 = compression[c0];
                int u0 = l0 + compression_bound[c0];
                int v0 = candidate_pairs[j][c0];
                if (!(l0 <= v0 && v0 <= u0)) {
                    ok = false;
                    break;
                }

                int l1 = compression[c1];
                int u1 = l1 + compression_bound[c1];
                int v1 = candidate_pairs[j][c1];
                if (!(l1 <= v1 && v1 <= u1)) {
                    ok = false;
                }
            }
        }
        is_less_cand[ctx->coset_idx][j] = ok;
        exists = exists || ok;
    }
    return exists;
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

static bool compression_state_matches_candidates(const DFSContext *ctx, int flag) {
    const int *compression = (flag == 1) ? ctx->compression_1 : ctx->compression_2;
    int **candidate_pairs = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    size_t num_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int dimension = (flag == 1) ? (int)ctx->dimension_candidate_pairs1 : (int)ctx->dimension_candidate_pairs2;

    for (size_t r = 0; r < num_pairs; r++) {
        bool same = true;
        for (int c = 0; c < dimension; c++) {
            if (compression[c] != candidate_pairs[r][c]) {
                same = false;
                break;
            }
        }
        if (same) {
            return true;
        }
    }
    return false;
}

static bool check_weight_feasibility(const DFSContext *ctx) {
    if (ctx->current_ones1 > ctx->target_ones || ctx->current_ones2 > ctx->target_ones) {
        return false;
    }
    if (ctx->current_ones1 + ctx->suffix_ones1[ctx->coset_idx] < ctx->target_ones ||
        ctx->current_ones2 + ctx->suffix_ones2[ctx->coset_idx] < ctx->target_ones) {
        return false;
    }
    return true;
}

bool check_bound(DFSContext *ctx) {
    if (!check_weight_feasibility(ctx)) {
        return false;
    }

    bool result = is_less_than_candidates(ctx, 1) &&
                  is_less_than_candidates(ctx, 0);
    if (!result) {
        return false;
    }
    return result;
    double max_lower_bound_sum = 0.0;
    for (int j = 1; j < ctx->spectrum_size; j++) {
        double current_abs_0 = cabs(ctx->current_dft[j]);
        double current_abs_1 = cabs(ctx->current_dft[j + ctx->N]);
        ctx->current_psd[j] = (int)rint(pow(current_abs_0, 2));
        ctx->current_psd[j + ctx->N] = (int)rint(pow(current_abs_1, 2));

        double lower_bound_psd_0 = pow(fmax(0.0, current_abs_0 - ctx->current_bound[j]), 2.0);
        double lower_bound_psd_1 = pow(fmax(0.0, current_abs_1 - ctx->current_bound[j + ctx->N]), 2.0);
        double lower_bound_sum = lower_bound_psd_0 + lower_bound_psd_1;
        if (lower_bound_sum > max_lower_bound_sum) {
            max_lower_bound_sum = lower_bound_sum;
        }
    }
    return max_lower_bound_sum <= ctx->threshold;
}

bool is_valid_combination(const DFSContext *ctx) {
    bool result = false;
    if (compression_state_matches_candidates(ctx, 0) &&
        compression_state_matches_candidates(ctx, 1))
    {
        result = true;
        for (int j = 1; (j < ctx->spectrum_size) && result; j++) {
            double real_part_0 = creal(ctx->current_dft[j]);
            double imag_part_0 = cimag(ctx->current_dft[j]);
            int psd_0 = (int)rint(real_part_0 * real_part_0 + imag_part_0 * imag_part_0);

            double real_part_1 = creal(ctx->current_dft[j + ctx->N]);
            double imag_part_1 = cimag(ctx->current_dft[j + ctx->N]);
            int psd_1 = (int)rint(real_part_1 * real_part_1 + imag_part_1 * imag_part_1);

            ctx->current_psd[j] = psd_0;
            ctx->current_psd[j + ctx->N] = psd_1;
            result = ((psd_0 + psd_1) == ctx->threshold);
        }
    }
    return result;
}

void dfs_explore_combinations(DFSContext *ctx)
{
    if (ctx->current_ones1 > ctx->target_ones || ctx->current_ones2 > ctx->target_ones) {
        return;
    }

    if (ctx->current_ones1 + ctx->suffix_ones1[ctx->coset_idx] < ctx->target_ones ||
        ctx->current_ones2 + ctx->suffix_ones2[ctx->coset_idx] < ctx->target_ones ) {
        return;
    }

    /* Caso base: hemos procesado todos los cosets o ya tenemos el peso objetivo. */
    if (ctx->coset_idx == ctx->cl->len || 
        (ctx->current_ones1 == ctx->target_ones && ctx->current_ones2 == ctx->target_ones)) {
        if (is_valid_combination(ctx)) {
            ctx->matches_found++;
            for (int j = 0; j < ctx-> N * 2; j++) {
                fprintf(ctx->f_comb, "%u", ctx->current_sequence[j]);
            }
            fprintf(ctx->f_comb, "\n");

            // Recalcular DFT exacta desde la secuencia actual
            //double complex *time_domain_exact = malloc(ctx->N * sizeof(double complex));
            //binary_to_complex(ctx->current_sequence, time_domain_exact, ctx->N);
            //dft(time_domain_exact, ctx->current_dft, ctx->N);
            //free(time_domain_exact);
            for (int j = 1; j < ctx->spectrum_size; j++) {
                int psd_val = ctx->current_psd[j];
                int transformed = ctx->current_psd[j+ctx->N];

                if (j > 1) {
                    fprintf(ctx->f_psd, " ");
                    fprintf(ctx->f_cte, " ");
                }
                fprintf(ctx->f_psd, "%d", psd_val);
                fprintf(ctx->f_cte, "%d", transformed);
            }
            fprintf(ctx->f_psd, "\n");
            fprintf(ctx->f_cte, "\n");
            fflush(ctx->f_psd);
            fflush(ctx->f_cte);
            //find_matches_files(ctx->N, (size_t)(ctx->N - 1), ctx->psd_filename, ctx->cte_filename, ctx->lp_filename, ctx->comb_filename);

        }
        return;
    }
    printf("Exploring coset index %d with current ones (%d, %d) and target ones %d\n", ctx->coset_idx, 
                                                                                    ctx->current_ones1, 
                                                                                    ctx->current_ones2, 
                                                                                    ctx->target_ones);
    printf("Current sequence:\n");
    for (int i = 0; i < 2* ctx->N; i++) {
        printf("%d ", ctx->current_sequence[i]);
    }
    printf("\n");
    printf("compressions:\n");
    int dimension_candidate_pairs = ctx->dimension_candidate_pairs1;
    for (int i = 0; i < dimension_candidate_pairs; i++) {
        printf("%d ", ctx->compression_1[i]);
    }
    printf("\n");
    int num_pairs = ctx->num_candidate_pairs1;
    int ** candidate_pairs = ctx->candidate_pairs1;
    printf("candidates:\n");
    for (size_t j = 0; j < num_pairs; j++) {
        for (int c = 0; c < dimension_candidate_pairs; c++) {
            printf("%d ", candidate_pairs[j][c]);
        }
        printf("\n");
    }
    for (int alfabet_val_0 = 0; alfabet_val_0 < ALFABET_SIZE; alfabet_val_0++) {
        for (int alfabet_val_1 = 0; alfabet_val_1 < ALFABET_SIZE; alfabet_val_1++) {
            double complex *current_dft_backup = malloc((size_t)(2 * ctx->N) * sizeof(double complex));
            if (!current_dft_backup) {
                return;
            }
            memcpy(current_dft_backup, ctx->current_dft, (size_t)(2 * ctx->N) * sizeof(double complex));

            if (alfabet_val_0 == 1) {
                ctx->current_ones1 += ctx->cl->data[ctx->coset_idx].len;
                ctx->current_combination1[ctx->pos_depth_1] = ctx->coset_idx;
                ctx->pos_depth_1++;
                ctx->current_combination1[ctx->pos_depth_1] = -1;
            }

            if (alfabet_val_1 == 1) {
                ctx->current_ones2 += ctx->cl->data[ctx->coset_idx].len;
                ctx->current_combination2[ctx->pos_depth_2] = ctx->coset_idx;
                ctx->pos_depth_2++;
                ctx->current_combination2[ctx->pos_depth_2] = -1;
            }

            int dim1 = (int)ctx->dimension_candidate_pairs1 / 2;
            int dim2 = (int)ctx->dimension_candidate_pairs2 / 2;
            Coset *selected_coset = &ctx->cl->data[ctx->coset_idx];
            for (int t = 0; t < selected_coset->len; t++) {
                int elem = selected_coset->data[t];
                ctx->current_sequence[elem] = alfabet_val_0;
                ctx->current_sequence[elem + ctx->N] = alfabet_val_1;

                ctx->compression_1[elem % dim1] += alfabet_val_0;
                ctx->compression_1[(elem % dim1) + dim1] += alfabet_val_1;
                ctx->compression_2[elem % dim2] += alfabet_val_0;
                ctx->compression_2[(elem % dim2) + dim2] += alfabet_val_1;
                ctx->compression_3[elem % ctx->p] += alfabet_val_0;
                ctx->compression_3[(elem % ctx->p) + ctx->p] += alfabet_val_1;

                ctx->bound_compression_1[elem % dim1] -= 1;
                ctx->bound_compression_1[(elem % dim1) + dim1] -= 1;
                ctx->bound_compression_2[elem % dim2] -= 1;
                ctx->bound_compression_2[(elem % dim2) + dim2] -= 1;
                ctx->bound_compression_3[elem % ctx->p] -= 1;
                ctx->bound_compression_3[(elem % ctx->p) + ctx->p] -= 1;
            }

            for (int j = 0; j < ctx->spectrum_size; j++) {
                ctx->current_dft[j] += alfabet_val_0 * ctx->dft_matrix[ctx->coset_idx][j];
                ctx->current_dft[j + ctx->N] += alfabet_val_1 * ctx->dft_matrix[ctx->coset_idx][j];
            }

            ctx->coset_idx++;
            ctx->current_bound = ctx->depth_exploration_bound[ctx->coset_idx];
            bool can_recurse = check_bound(ctx);
            if (can_recurse) {
                dfs_explore_combinations(ctx);
            }

            memcpy(ctx->current_dft, current_dft_backup, (size_t)(2 * ctx->N) * sizeof(double complex));
            free(current_dft_backup);

            ctx->coset_idx--;
            for (int t = 0; t < selected_coset->len; t++) {
                int elem = selected_coset->data[t];
                ctx->current_sequence[elem] = 0;
                ctx->current_sequence[elem + ctx->N] = 0;

                ctx->compression_1[elem % dim1] -= alfabet_val_0;
                ctx->compression_1[(elem % dim1) + dim1] -= alfabet_val_1;
                ctx->compression_2[elem % dim2] -= alfabet_val_0;
                ctx->compression_2[(elem % dim2) + dim2] -= alfabet_val_1;
                ctx->compression_3[elem % ctx->p] -= alfabet_val_0;
                ctx->compression_3[(elem % ctx->p) + ctx->p] -= alfabet_val_1;

                ctx->bound_compression_1[elem % dim1] += 1;
                ctx->bound_compression_1[(elem % dim1) + dim1] += 1;
                ctx->bound_compression_2[elem % dim2] += 1;
                ctx->bound_compression_2[(elem % dim2) + dim2] += 1;
                ctx->bound_compression_3[elem % ctx->p] += 1;
                ctx->bound_compression_3[(elem % ctx->p) + ctx->p] += 1;
            }

            if (alfabet_val_0 == 1) {
                ctx->current_ones1 -= ctx->cl->data[ctx->coset_idx].len;
                ctx->pos_depth_1--;
                ctx->current_combination1[ctx->pos_depth_1] = -1;
            }
            if (alfabet_val_1 == 1) {
                ctx->current_ones2 -= ctx->cl->data[ctx->coset_idx].len;
                ctx->pos_depth_2--;
                ctx->current_combination2[ctx->pos_depth_2] = -1;
            }
        }
    }
}

/*
            if (alfabet_val_0 == 1 && alfabet_val_1 == 1) {
            ctx->current_combination[ctx->pos_depth] = ctx->coset_idx;
            ctx->current_ones += ctx->cl->data[ctx->coset_idx].len;
            ctx->pos_depth++;
            ctx->current_combination[ctx->pos_depth] = -1;

            Coset *selected_coset = &ctx->cl->data[ctx->coset_idx];
            for (int t = 0; t < selected_coset->len; t++) {
                int elem = selected_coset->data[t];
                if (elem >= 0 && elem < ctx->N) {
                    ctx->current_sequence[elem] = 1;
                    ctx->compresion_1[elem % (int)ctx->dimension_candidate_pairs1] += 1;
                    ctx->compresion_2[elem % (int)ctx->dimension_candidate_pairs2] += 1;
                    ctx->compresion_3[elem % ctx->p] += 1;
                }
            }

            for (int j = 0; j < ctx->spectrum_size; j++) {
                ctx->dft_backup[ctx->coset_idx][j] = ctx->current_dft[j];
                ctx->current_dft[j] = ctx->current_dft[j] + ctx->dft_matrix[ctx->coset_idx][j];
            }
        }

        Coset *processed_coset = &ctx->cl->data[ctx->coset_idx];
        for (int t = 0; t < processed_coset->len; t++) {
            int elem = processed_coset->data[t];
            if (elem >= 0 && elem < ctx->N) {
                ctx->remaining_bits[elem] = 0;
                ctx->bound_compression_1[elem % (int)ctx->dimension_candidate_pairs1] -= 1;
                ctx->bound_compression_2[elem % (int)ctx->dimension_candidate_pairs2] -= 1;
                ctx->bound_contrassion_3[elem % ctx->p] -= 1;
            }
        }

        ctx->coset_idx++;
        ctx->current_bound = ctx->depth_exploration_bound[ctx->coset_idx];

        if (ctx->current_ones <= ctx->target_ones &&
            ctx->current_ones + ctx->suffix_ones[ctx->coset_idx] >= ctx->target_ones &&
            check_bound(ctx)) {
            dfs_explore_combinations(ctx);
        }

        memcpy(ctx->current_dft, current_dft_backup, ctx->N * sizeof(double complex));

        free(current_dft_backup);

        ctx->coset_idx--;
        for (int t = 0; t < processed_coset->len; t++) {
            int elem = processed_coset->data[t];
            if (elem >= 0 && elem < ctx->N) {
                ctx->remaining_bits[elem] = 1;
                ctx->bound_compression_1[elem % (int)ctx->dimension_candidate_pairs1] += 1;
                ctx->bound_compression_2[elem % (int)ctx->dimension_candidate_pairs2] += 1;
                ctx->bound_contrassion_3[elem % ctx->p] += 1;
            }
        }

        if (alfabet_val == 1) {
            Coset *selected_coset = &ctx->cl->data[ctx->coset_idx];
            for (int t = 0; t < selected_coset->len; t++) {
                int elem = selected_coset->data[t];
                if (elem >= 0 && elem < ctx->N) {
                    ctx->current_sequence[elem] = 0;
                    ctx->compresion_1[elem % (int)ctx->dimension_candidate_pairs1] -= 1;
                    ctx->compresion_2[elem % (int)ctx->dimension_candidate_pairs2] -= 1;
                    ctx->compresion_3[elem % ctx->p] -= 1;
                }
            }

            ctx->pos_depth--;
            ctx->current_ones -= ctx->cl->data[ctx->coset_idx].len;
            ctx->current_combination[ctx->pos_depth] = -1;
        }
    }*/



// ...existing code (compute_depth_exploration_bounds, free_depth_exploration_bounds)...

static double **compute_depth_exploration_bounds(DFSContext *ctx) {
    int num_cosets = ctx->num_cosets;
    int N = ctx->N;

    double **bounds = calloc(num_cosets + 1, sizeof(double *));
    if (!bounds) return NULL;

    int *combo_arr = calloc(ctx->cl->len, sizeof(int));
    int *bound_arr = calloc(ctx->cl->len, sizeof(int));
    if (!combo_arr || !bound_arr) {
        free(combo_arr); free(bound_arr);
        free(bounds);
        return NULL;
    }

    for (int i = 0; i <= num_cosets; i++) {
        bounds[i] = calloc(2 * N, sizeof(double));
        if (!bounds[i]) {
            for (int k = 0; k < i; k++) free(bounds[k]);
            free(bounds);
            free(combo_arr); free(bound_arr);
            return NULL;
        }
    }

    int direct_limit = (DEPTH - 1 < num_cosets) ? (DEPTH - 1) : num_cosets;
    for (int i = 1; i <= direct_limit; i++) {
        for (int half = 0; half < 2; half++) {
            int out_offset = half * N;
            for (int j = 0; j < ctx->spectrum_size; j++) {
                double acc = 0.0;
                for (int k = i; k < num_cosets; k++) {
                    acc += cabs(ctx->dft_matrix[k][j]);
                }
                bounds[i][j + out_offset] = acc;
            }
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
                for (int half = 0; half < 2; half++) {
                    int out_offset = half * N;
                    for (int j = 0; j < ctx->spectrum_size; j++) {
                        double acc = 0.0;
                        for (int k = i; k < num_cosets; k++) {
                            acc += cabs(ctx->dft_matrix[k][j]);
                        }
                        bounds[i][j + out_offset] = acc;
                    }
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
                
                /*
                 * No podar por PSD del sufijo aislado: aunque prefijo y sufijo
                 * no compartan unos en tiempo, en frecuencia si puede haber
                 * cancelacion entre ambos. Este bound debe seguir siendo
                 * conservador respecto al prefijo actual.
                 */

                if (is_less_than_candidates_for_vectors(ctx, seq, bound_seq, 0)) {
                    for (int j = 0; j < ctx->spectrum_size; j++) {
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
                if (is_less_than_candidates_for_vectors(ctx, seq, bound_seq, 1)) {
                    for (int j = 0; j < ctx->spectrum_size; j++) {
                        double complex val = 0.0;
                        for (int k = 0; k < d; k++) {
                            int bit = (combo >> k) & 1;
                            val += bit * ctx->dft_matrix[start_coset + k][j];
                        }
                        double abs_val = cabs(val);
                        if (abs_val > bounds[i][j+N]) {
                            bounds[i][j+N] = abs_val;
                        }
                    }
                }

                free(seq);
            }
            free(bound_seq);
        }
    }

    free(combo_arr);
    free(bound_arr);
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

    int threshold = (N+1)/2;
    int spectrum_size = (N  + 1)/2;
    //int spectrum_size = 20;
    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    double complex **dft_backup = malloc((cl->len + 1) * sizeof(double complex *));

    for (int i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        dft_backup[i] = malloc(2 * N * sizeof(double complex));
    }
    dft_backup[cl->len] = malloc(2 * N * sizeof(double complex)); // Para backup temporal en DFS

    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    double complex *base_dft = malloc(N * sizeof(double complex));
    int *bound_psd = calloc(N, sizeof(int));
    int *compression_a = malloc(2 * p * sizeof(int));
    double complex *current_dft = calloc((size_t)(2 * N), sizeof(double complex));
    int *current_sequence = calloc((size_t)(2 * N), sizeof(int)); 
    int *current_psd = calloc((size_t)(2 * N), sizeof(int));
    int *remaining_bits = calloc((size_t)(2 * N), sizeof(int));
    int *compresion_1 = NULL;
    int *compresion_2 = NULL;
    int *bound_compression_1 = NULL;
    int *bound_compression_2 = NULL;
    int *compresion_3 = NULL;
    int *bound_compression_3 = NULL;
    int *combination1 = malloc(cl->len * sizeof(int)); // these are the cosets for the first sequence
    int *combination2 = malloc(cl->len * sizeof(int)); // These are the cosets for the second sequence
    int *suffix_ones1 = calloc(cl->len + 1, sizeof(int));
    int *suffix_ones2 = calloc(cl->len + 1, sizeof(int));
    bool *is_less_compression_a = malloc((size_t)(cl->len + 1) * sizeof(bool));
    double *current_bound = malloc((2 * N) * sizeof(double));

    if (!current_bound || !combination1 || !combination2 || !suffix_ones1 || !suffix_ones2 || !current_sequence || !current_psd || !remaining_bits ||
        !is_less_compression_a || !base_dft) {
        printf("ERROR: No se pudo asignar memoria\n");
        free(current_dft);
        free(current_sequence);
        free(current_psd);
        free(remaining_bits);
        free(combination1);
        free(combination2);
        free(suffix_ones1);
        free(suffix_ones2);
        free(is_less_compression_a);
        free(compression_a);
        free(compresion_3);
        free(bound_compression_3);
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(dft_backup[i]);    
        }
        free(dft_backup[cl->len]);
        free(dft_matrix);
        free(dft_backup);
        free(dft_backup);
        return;
    }

    /* Inicializar combination como lista vacía terminada en -1 */
    for (int i = 0; i < cl->len; i++) {
        combination1[i] = -1;
        combination2[i] = -1;
    }
    for (int i = 0; i <= cl->len; i++) {
        is_less_compression_a[i] = true;
    }

    legendre_sequence(p, q, compression_a);

    /* DFT base para la secuencia nula (todo 0 en binario -> todo 1 en dominio temporal). */
    int *zero_bits = calloc((size_t)N, sizeof(int));
    if (!zero_bits) {
        printf("ERROR: No se pudo asignar memoria para zero_bits\n");
        free(base_dft);
        return;
    }
    binary_to_complex(zero_bits, time_domain, N);
    dft(time_domain, base_dft, N);
    free(zero_bits);

    /* Calcular contribucion DFT por coset respecto a la base. */
    int *single_coset_comb = malloc(cl->len * sizeof(int));
    for (int i = 0; i < cl->len; i++) single_coset_comb[i] = -1;

    for (int i = 0; i < cl->len; i++) {
        single_coset_comb[0] = i;
        single_coset_comb[1] = -1;

        int *vector_bits = generate_vector_for_combination(cl, single_coset_comb, N);
        binary_to_complex(vector_bits, time_domain, N);
        dft(time_domain, freq_domain, N);

        for (int j = 0; j < N; j++) {
            dft_matrix[i][j] = freq_domain[j] - base_dft[j];
            dft_backup[i][j] = 0.0 + 0.0 * I;
            dft_backup[i][j+N] = 0.0 + 0.0 * I;
        }
        for (int j = 0; j < spectrum_size; j++) {
            bound_psd[j] += (int)rint(pow(cabs(freq_domain[j]), 2));
        }
        free(vector_bits);
    }
    free(single_coset_comb);

    free(time_domain);
    free(freq_domain);

    /* Asegurar estado inicial DFS vacío */
    for (int i = 0; i < cl->len; i++) {
        combination1[i] = -1;
        combination2[i] = -1;
    }
    for (int i = 0; i < 2 * N; i++) {
        remaining_bits[i] = 0;
    }
    for (int i = 0; i < cl->len; i++) {
        Coset *c = &cl->data[i];
        for (int j = 0; j < c->len; j++) {
            int elem = c->data[j];
            if (elem >= 0 && elem < N) {
                remaining_bits[elem] = 1;
                remaining_bits[elem + N] = 1; // Para la segunda secuencia
            }
        }
    }

    suffix_ones1[cl->len] = 0;
    suffix_ones2[cl->len] = 0;
    for (int i = cl->len - 1; i >= 0; i--) {
        suffix_ones1[i] = suffix_ones1[i + 1] + cl->data[i].len;
        suffix_ones2[i] = suffix_ones2[i + 1] + cl->data[i].len;
    }

    // Inicializar current_bound con la suma de valores absolutos de la DFT para cada frecuencia
    /*
     * current_dft[0..N-1]  -> DFT incremental de current_sequence[0..N-1]
     * current_dft[N..2N-1] -> DFT incremental de current_sequence[N..2N-1]
     * Ambas parten de la DFT base (secuencia todo 0 en binario) y se actualizan
     * sumando/restando deltas por coset durante el DFS.
     */
    for (int j = 0; j < spectrum_size; j++) {
        current_bound[j] = 0.0;
        current_bound[j + N] = 0.0;
        current_dft[j] = base_dft[j];
        current_dft[j + N] = base_dft[j];
        current_psd[j] = 0;
        current_psd[j + N] = 0;
        for (int i = 0; i < cl->len; i++) {
            double abs_delta = cabs(dft_matrix[i][j]);
            current_bound[j] += abs_delta;
            current_bound[j + N] += abs_delta;
        }
    }

    // Abrir archivos de salida
    char comb_filename[512];
    char psd_filename[512];
    char cte_filename[512];
    char lp_filename[512];
    char output_dir[512];
    const char *pairs_name1 = path_basename(pairs_file1);
    const char *pairs_name2 = path_basename(pairs_file2);
    path_dirname(pairs_file1, output_dir, sizeof(output_dir));

    snprintf(comb_filename, sizeof(comb_filename), "%s/%s_%s_combinations.txt", output_dir, pairs_name1, pairs_name2);
    snprintf(psd_filename, sizeof(psd_filename), "%s/%s_%s_dft.txt", output_dir, pairs_name1, pairs_name2);
    snprintf(cte_filename, sizeof(cte_filename), "%s/%s_%s_cte-dft.txt", output_dir, pairs_name1, pairs_name2);
    snprintf(lp_filename, sizeof(lp_filename), "%s/%s_%s_lp.txt", output_dir, pairs_name1, pairs_name2);

    FILE *f_comb = fopen(comb_filename, "w");
    FILE *f_psd = fopen(psd_filename, "w");
    FILE *f_cte = fopen(cte_filename, "w");
    size_t rows1, cols1;
    int **candidate1 = read_pairs_file(pairs_file1, &rows1, &cols1);
    size_t rows2, cols2;
    int **candidate2 = read_pairs_file(pairs_file2, &rows2, &cols2);

    compresion_1 = calloc(cols1 ? cols1 : 1, sizeof(int));
    compresion_2 = calloc(cols2 ? cols2 : 1, sizeof(int));
    bound_compression_1 = calloc(cols1 ? cols1 : 1, sizeof(int));
    bound_compression_2 = calloc(cols2 ? cols2 : 1, sizeof(int));
    compresion_3 = calloc((size_t)(2 * p), sizeof(int));
    bound_compression_3 = calloc((size_t)(2 * p ), sizeof(int));

    int max_value1 = (cols1 > 0) ? (2 * N / (int)cols1) : 0;
    int max_value2 = (cols2 > 0) ? (2 * N / (int)cols2) : 0;
    bool **is_less_candidates_1 = NULL;
    bool **is_less_candidates_2 = NULL;
    if (rows1 > 0) {
        is_less_candidates_1 = allocate_bool_matrix((size_t)(cl->len + 1), rows1);
    }
    if (rows2 > 0) {
        is_less_candidates_2 = allocate_bool_matrix((size_t)(cl->len + 1), rows2);
    }
    if (!f_comb || !f_psd || !f_cte || !candidate1 || !candidate2 ||
        !compresion_1 || !compresion_2 || !bound_compression_1 || !bound_compression_2 ||
        !compresion_3 || !bound_compression_3 ||
        cols1 == 0 || cols2 == 0) {
        printf("ERROR: No se pudieron abrir los archivos o leer los pares de: %s y %s.\n", pairs_file1, pairs_file2);
        if (!f_comb) perror("fopen comb_filename");
        if (!f_psd) perror("fopen psd_filename");
        if (!f_cte) perror("fopen cte_filename");
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);
        free(current_bound);
        free(base_dft);

        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
        }
        free(dft_matrix);
        free(current_dft);
        free(current_sequence);
        free(remaining_bits);
        free(compresion_1);
        free(compresion_2);
        free(bound_compression_1);
        free(bound_compression_2);
        free(bound_psd);
        free(combination1);
        free(combination2);
        free(suffix_ones1);
        free(suffix_ones2);
        free(is_less_compression_a);
        free(compression_a);
        free(compresion_3);
        free(bound_compression_3);
        return;
    }
    for (int i = 0; i < N; i++) {
        if (remaining_bits[i]) {
            int half_cols1 = (int)(cols1 / 2);
            int half_cols2 = (int)(cols2 / 2);
            bound_compression_1[i % half_cols1] += 1;
            bound_compression_1[(i % half_cols1) + half_cols1] += 1;
            bound_compression_2[i % half_cols2] += 1;
            bound_compression_2[(i % half_cols2) + half_cols2] += 1;
            bound_compression_3[i % p] += 1;
            bound_compression_3[(i % p) + p] += 1;
        }
    }

    combination1[0] = -1;
    combination2[0] = -1;
    // Preparar contexto
    DFSContext ctx = {
        .f_comb = f_comb,
        .f_psd = f_psd,
        .f_cte = f_cte,
        .comb_filename = "",
        .psd_filename = "",
        .cte_filename = "",
        .lp_filename = "",
        .matches_found = 0,
        .threshold = threshold,
        .N = N,
        .p = p,
        .q = q,
        .num_cosets = cl->len,
        .dft_matrix = dft_matrix,
        .dft_backup = dft_backup,
        .current_dft = current_dft,
        .bound_psd = bound_psd,
        .current_sequence = current_sequence,
        .remaining_bits = remaining_bits,
        .compression_1 = compresion_1,
        .compression_2 = compresion_2,
        .bound_compression_1 = bound_compression_1,
        .bound_compression_2 = bound_compression_2,
        .compression_3 = compresion_3,
        .bound_compression_3 = bound_compression_3,
        .compression_a = compression_a,
        .cl = cl,
        .current_combination1 = combination1,
        .current_combination2 = combination2,
        .coset_idx = 0,
        .candidate_pairs1 = candidate1,
        .num_candidate_pairs1 = rows1,
        .dimension_candidate_pairs1 = cols1,
        .candidate_pairs2 = candidate2,
        .num_candidate_pairs2 = rows2,
        .dimension_candidate_pairs2 = cols2,
        .current_bound = current_bound,
        .depth_exploration_bound = NULL,
        .pos_depth_1 = 0,
        .pos_depth_2 = 0,
        .current_ones1 = 0,
        .current_ones2 = 0,
        .target_ones = (N + 1) / 2,
        .spectrum_size = spectrum_size,
        .suffix_ones1 = suffix_ones1,
        .suffix_ones2 = suffix_ones2,
        .current_psd = current_psd,
        .is_less_compression_a = is_less_compression_a,
        .is_less_candidates_1 = is_less_candidates_1,
        .is_less_candidates_2 = is_less_candidates_2,
    };
    snprintf(ctx.comb_filename, sizeof(ctx.comb_filename), "%s", comb_filename);
    snprintf(ctx.psd_filename, sizeof(ctx.psd_filename), "%s", psd_filename);
    snprintf(ctx.cte_filename, sizeof(ctx.cte_filename), "%s", cte_filename);
    snprintf(ctx.lp_filename, sizeof(ctx.lp_filename), "%s", lp_filename);
    reorder_cosets_by_candidate_pairs(&ctx);
    
    // Actualizar punteros locales después de reordenación
    dft_matrix = ctx.dft_matrix;    
    
    ctx.depth_exploration_bound = compute_depth_exploration_bounds(&ctx);
    if (ctx.depth_exploration_bound == NULL) {
        printf("ERROR: No se pudieron calcular los bounds de exploración por profundidad.\n");
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
        }
        free(dft_matrix);
        free(current_dft);
        free(current_sequence);
        free(remaining_bits);
        free(compresion_1);
        free(compresion_2);
        free(bound_compression_1);
        free(bound_compression_2);
        free(bound_psd);
        free(combination1);
        free(combination2);
        free(suffix_ones1);
        free(suffix_ones2);
        free(is_less_compression_a);
        free_bool_matrix(is_less_candidates_1, (size_t)(cl->len + 1));
        free_bool_matrix(is_less_candidates_2, (size_t)(cl->len + 1));
        free(compression_a);
        free(compresion_3);
        free(bound_compression_3);
        free(base_dft);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);
        free(current_bound);
        fclose(f_comb);
        fclose(f_psd);
        fclose(f_cte);
        return;
    }

    printf("\n=== EXPLORACIÓN DFS (N=%d, cosets=%d) ===\n", N, cl->len);
    printf("Condicion: Max(PSD) < %d\n\n", threshold);

    /* Asegurar que combination empieza vacío antes de DFS */
    for (int i = 0; i < cl->len; i++) 
    {
        combination1[i] = -1;
        combination2[i] = -1;
    }
    ctx.coset_idx = 0;
    ctx.pos_depth_1 = 0;
    ctx.pos_depth_2 = 0;
    ctx.current_ones1 = 0;
    ctx.current_ones2 = 0;

    dfs_explore_combinations(&ctx);

    // Limpieza completa de memoria
    for (int i = 0; i < cl->len; i++) {
        free(dft_matrix[i]);
        free(dft_backup[i]);
    }
    free(dft_matrix);
    free(current_dft);
    free(current_sequence);
    free(remaining_bits);
    free(compresion_1);
    free(compresion_2);
    free(bound_compression_1);
    free(bound_compression_2);
    free(bound_psd);
    free(combination1);
    free(combination2);
    free(suffix_ones1);
    free(suffix_ones2);
    free(is_less_compression_a);    
    free_bool_matrix(is_less_candidates_1, (size_t)(cl->len + 1));
    free_bool_matrix(is_less_candidates_2, (size_t)(cl->len + 1));
    free(compression_a);
    free(compresion_3);
    free(bound_compression_3);
    free(base_dft);
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
