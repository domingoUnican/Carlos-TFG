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
#define DEPTH 15
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
    int *bound_psd;
    int *current_sequence;
    int *remaining_bits;
    int *compresion_1;
    int *compresion_2;
    int *bound_compression_1;
    int *bound_compression_2;
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
    double *current_bound;
    double **depth_exploration_bound;
    size_t num_candidate_pairs2;
    size_t dimension_candidate_pairs2;
    int pos_depth;
    int current_ones;
    int target_ones;
    int spectrum_size;
    int *suffix_ones;
    uint64_t ***ge_bitsets1;
    uint64_t ***le_bitsets1;
    size_t bitset_words1;
    int bitset_max_value1;
    uint64_t *bitset_work1;
    int *bitset_order1;
    int *bitset_score1;
    uint64_t ***ge_bitsets2;
    uint64_t ***le_bitsets2;
    size_t bitset_words2;
    int bitset_max_value2;
    uint64_t *bitset_work2;
    int *bitset_order2;
    int *bitset_score2;
} DFSContext;

static void free_candidate_bitsets(uint64_t ***ge, uint64_t ***le, size_t dim, int max_value) {
    if (ge) {
        for (size_t j = 0; j < dim; j++) {
            if (ge[j]) {
                for (int v = 0; v <= max_value; v++) {
                    free(ge[j][v]);
                }
                free(ge[j]);
            }
        }
        free(ge);
    }
    if (le) {
        for (size_t j = 0; j < dim; j++) {
            if (le[j]) {
                for (int v = 0; v <= max_value; v++) {
                    free(le[j][v]);
                }
                free(le[j]);
            }
        }
        free(le);
    }
}

static bool build_candidate_bitsets(int **pairs, size_t num_pairs, size_t dim, int max_value,
                                    uint64_t ****ge_out, uint64_t ****le_out, size_t *words_len_out) {
    if (!pairs || dim == 0 || !ge_out || !le_out || !words_len_out || max_value < 0) {
        return false;
    }

    size_t words_len = (num_pairs + 63u) / 64u;
    if (words_len == 0) {
        words_len = 1;
    }

    uint64_t ***ge = calloc(dim, sizeof(uint64_t **));
    uint64_t ***le = calloc(dim, sizeof(uint64_t **));
    if (!ge || !le) {
        free(ge);
        free(le);
        return false;
    }

    for (size_t j = 0; j < dim; j++) {
        ge[j] = calloc((size_t)max_value + 1u, sizeof(uint64_t *));
        le[j] = calloc((size_t)max_value + 1u, sizeof(uint64_t *));
        if (!ge[j] || !le[j]) {
            free_candidate_bitsets(ge, le, dim, max_value);
            return false;
        }
        for (int v = 0; v <= max_value; v++) {
            ge[j][v] = calloc(words_len, sizeof(uint64_t));
            le[j][v] = calloc(words_len, sizeof(uint64_t));
            if (!ge[j][v] || !le[j][v]) {
                free_candidate_bitsets(ge, le, dim, max_value);
                return false;
            }
        }
    }

    for (size_t i = 0; i < num_pairs; i++) {
        size_t word = i >> 6;
        uint64_t bit = 1ULL << (i & 63u);
        for (size_t j = 0; j < dim; j++) {
            int value = pairs[i][j];
            if (value < 0) value = 0;
            if (value > max_value) value = max_value;

            for (int l = 0; l <= value; l++) {
                ge[j][l][word] |= bit;
            }
            for (int u = value; u <= max_value; u++) {
                le[j][u][word] |= bit;
            }
        }
    }

    *ge_out = ge;
    *le_out = le;
    *words_len_out = words_len;
    return true;
}

static bool exists_candidate_in_range_bitset(const int *compression, const int *compression_bound,
                                               size_t dim,
                                               uint64_t ***ge, uint64_t ***le,
                                               size_t words_len, size_t num_pairs, int max_value,
                                               uint64_t *work,
                                               int *order_buf,
                                               int *score_buf) {
    if (!compression || !compression_bound || !ge || !le || !work || !order_buf || !score_buf) {
        return false;
    }
    if (num_pairs == 0) {
        return false;
    }

    for (size_t w = 0; w < words_len; w++) {
        work[w] = ~0ULL;
    }
    size_t last_word = (num_pairs - 1u) >> 6;
    size_t used_bits = ((num_pairs - 1u) & 63u) + 1u;
    if (used_bits < 64u) {
        uint64_t tail_mask = (1ULL << used_bits) - 1ULL;
        work[last_word] &= tail_mask;
    }

    for (size_t j = 0; j < dim; j++) {
        int l = compression[j];
        int u = compression[j] + compression_bound[j];

        if (u < 0 || l > max_value) {
            return false;
        }
        if (l < 0) l = 0;
        if (u > max_value) u = max_value;
        if (l > u) {
            return false;
        }

        int count = 0;
        for (size_t w = 0; w <= last_word; w++) {
            uint64_t bits = ge[j][l][w] & le[j][u][w];
            count += __builtin_popcountll((unsigned long long)bits);
        }
        if (count == 0) {
            return false;
        }

        order_buf[j] = (int)j;
        score_buf[j] = count;
    }

    for (size_t i = 1; i < dim; i++) {
        int key_order = order_buf[i];
        int key_score = score_buf[i];
        size_t p = i;
        while (p > 0 && score_buf[p - 1] > key_score) {
            score_buf[p] = score_buf[p - 1];
            order_buf[p] = order_buf[p - 1];
            p--;
        }
        score_buf[p] = key_score;
        order_buf[p] = key_order;
    }

    for (size_t s = 0; s < dim; s++) {
        int j = order_buf[s];
        int l = compression[j];
        int u = compression[j] + compression_bound[j];
        if (l < 0) l = 0;
        if (u > max_value) u = max_value;

        bool any = false;
        for (size_t w = 0; w <= last_word; w++) {
            work[w] &= ge[j][l][w] & le[j][u][w];
            if (work[w] != 0ULL) {
                any = true;
            }
        }
        if (!any) {
            return false;
        }
    }
    return true;
}

static int compare_int_asc(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    return (ia > ib) - (ia < ib);
}

void reorder_cosets_by_candidate_pairs(DFSContext *ctx) {
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
        reordered_psd[i] = ctx->psd_matrix[old_idx];
    }

    free(ctx->cl->data);
    free(ctx->dft_matrix);
    free(ctx->psd_matrix);
    ctx->cl->data = reordered_cosets;
    ctx->dft_matrix = reordered_dft;
    ctx->psd_matrix = reordered_psd;

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

    if (ctx->suffix_ones) {
        ctx->suffix_ones[num_cosets] = 0;
        for (int i = num_cosets - 1; i >= 0; i--) {
            ctx->suffix_ones[i] = ctx->suffix_ones[i + 1] + ctx->cl->data[i].len;
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
    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    int *compression = malloc((size_t)dimension_candidate_pairs * sizeof(int));
    int *compression_bound = malloc((size_t)dimension_candidate_pairs * sizeof(int));
    uint64_t ***ge = (flag == 1) ? ctx->ge_bitsets1 : ctx->ge_bitsets2;
    uint64_t ***le = (flag == 1) ? ctx->le_bitsets1 : ctx->le_bitsets2;
    size_t words_len = (flag == 1) ? ctx->bitset_words1 : ctx->bitset_words2;
    size_t num_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int max_value = (flag == 1) ? ctx->bitset_max_value1 : ctx->bitset_max_value2;
    uint64_t *work = (flag == 1) ? ctx->bitset_work1 : ctx->bitset_work2;
    int *order_buf = (flag == 1) ? ctx->bitset_order1 : ctx->bitset_order2;
    int *score_buf = (flag == 1) ? ctx->bitset_score1 : ctx->bitset_score2;
    if (!compression || !compression_bound) {
        free(compression);
        free(compression_bound);
        return false;
    }

    CompressSequence(ctx->N, dimension_candidate_pairs, sequence, compression);
    CompressSequence(ctx->N, dimension_candidate_pairs, bound_bit_sequence, compression_bound);

    bool resultado = exists_candidate_in_range_bitset(compression, compression_bound,
                                                      (size_t)dimension_candidate_pairs,
                                                      ge, le, words_len, num_pairs, max_value,
                                                      work, order_buf, score_buf);
    free(compression);
    free(compression_bound);
    return resultado;
}

bool is_less_than_candidates(DFSContext *ctx, int flag)
{
    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    size_t num_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int *compression = (flag == 1) ? ctx->compresion_1 : ctx->compresion_2;
    int *compression_bound = (flag == 1) ? ctx->bound_compression_1 : ctx->bound_compression_2;
    uint64_t ***ge = (flag == 1) ? ctx->ge_bitsets1 : ctx->ge_bitsets2;
    uint64_t ***le = (flag == 1) ? ctx->le_bitsets1 : ctx->le_bitsets2;
    size_t words_len = (flag == 1) ? ctx->bitset_words1 : ctx->bitset_words2;
    int max_value = (flag == 1) ? ctx->bitset_max_value1 : ctx->bitset_max_value2;
    uint64_t *work = (flag == 1) ? ctx->bitset_work1 : ctx->bitset_work2;
    int *order_buf = (flag == 1) ? ctx->bitset_order1 : ctx->bitset_order2;
    int *score_buf = (flag == 1) ? ctx->bitset_score1 : ctx->bitset_score2;

    bool resultado = exists_candidate_in_range_bitset(compression, compression_bound,
                                                      (size_t)dimension_candidate_pairs,
                                                      ge, le, words_len, num_pairs, max_value,
                                                      work, order_buf, score_buf);
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
    bool result = is_less_than_candidates(ctx, 1) &&
                  is_less_than_candidates(ctx, 0);
    if (!result) {
        return false;
    }
    double max_lower_bound = 0.0;
    for (int j = 1; j < ctx->spectrum_size; j++) {
        double current_abs = cabs(ctx->current_dft[j]);
        double lower_bound_psd = pow(fmax(0.0, current_abs - ctx->current_bound[j]), 2.0);
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

    bool result = true;
    if (is_compression_of_candidates(ctx,temp,0) && is_compression_of_candidates(ctx,temp,1))
    {
        int max_psd = -1;
        for (int j = 1; (j < ctx->spectrum_size) && result; j++) {
            if ((int)rint(pow(cabs(ctx->current_dft[j]), 2)) > max_psd) {
                max_psd = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
            }
        }
        result = (max_psd <= (int)ctx->threshold);
    }

    free(temp);
    return result;
}

void dfs_explore_combinations(DFSContext *ctx)
{
    if (ctx->current_ones > ctx->target_ones) {
        return;
    }

    if (ctx->current_ones + ctx->suffix_ones[ctx->coset_idx] < ctx->target_ones) {
        return;
    }

    /* Caso base: hemos procesado todos los cosets o ya tenemos el peso objetivo. */
    if (ctx->coset_idx == ctx->cl->len || ctx->current_ones == ctx->target_ones) {
        if (is_valid_combination(ctx)) {
            ctx->matches_found++;
            int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
            for (int j = 0; j < ctx->N; j++) {
                fprintf(ctx->f_comb, "%u", vector_bits[j]);
            }
            fprintf(ctx->f_comb, "\n");

            for (int j = 1; j < ctx->spectrum_size; j++) {
                int psd_val = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
                fprintf(ctx->f_psd, "%d%s", psd_val, (j + 1 < ctx->spectrum_size) ? " " : "");

                int transformed = (int)rint(ctx->threshold - psd_val);
                fprintf(ctx->f_cte, "%d%s", transformed, (j + 1 < ctx->spectrum_size) ? " " : "");
            }
            fprintf(ctx->f_psd, "\n");
            fprintf(ctx->f_cte, "\n");

            free(vector_bits);
        }
        return;
    }

    for (int alfabet_val = 0; alfabet_val < ALFABET_SIZE; alfabet_val++) {
        double complex *current_dft_backup = malloc(ctx->N * sizeof(double complex));
        memcpy(current_dft_backup, ctx->current_dft, ctx->N * sizeof(double complex));

        if (alfabet_val == 1) {
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
                }
            }

            for (int j = 0; j < ctx->spectrum_size; j++) {
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
                }
            }

            ctx->pos_depth--;
            ctx->current_ones -= ctx->cl->data[ctx->coset_idx].len;
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
    if (!combo_arr || !bound_arr) {
        free(combo_arr); free(bound_arr);
        free(bounds);
        return NULL;
    }

    for (int i = 0; i <= num_cosets; i++) {
        bounds[i] = calloc(N, sizeof(double));
        if (!bounds[i]) {
            for (int k = 0; k < i; k++) free(bounds[k]);
            free(bounds);
            free(combo_arr); free(bound_arr);
            return NULL;
        }
    }

    int direct_limit = (DEPTH - 1 < num_cosets) ? (DEPTH - 1) : num_cosets;
    for (int i = 1; i <= direct_limit; i++) {
        for (int j = 0; j < ctx->spectrum_size; j++) {
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
                for (int j = 0; j < ctx->spectrum_size; j++) {
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
                
                bool psd_below_threshold = true;
                for (int j = 1; j < ctx->spectrum_size; j++) {
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
                              is_less_than_candidates_for_vectors(ctx, seq, bound_seq, 1) &&
                              is_less_than_candidates_for_vectors(ctx, seq, bound_seq, 0);

                if (passes) {
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

    double threshold = ((double)N + 1.0) / 2.0;
    int spectrum_size = (N / 2) + 1;
    //int spectrum_size = 20;
    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    int **psd_matrix = malloc(cl->len * sizeof(int *));

    for (int i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        psd_matrix[i] = malloc(N * sizeof(int));
    }

    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    int *bound_psd = calloc(N, sizeof(int));
    int *compression_a = malloc(p * sizeof(int));
    int *compression_b = malloc(p * sizeof(int));
    double complex *current_dft = calloc(N, sizeof(double complex));
    int *current_sequence = calloc((size_t)N, sizeof(int));
    int *remaining_bits = calloc((size_t)N, sizeof(int));
    int *compresion_1 = NULL;
    int *compresion_2 = NULL;
    int *bound_compression_1 = NULL;
    int *bound_compression_2 = NULL;
    int *combination = malloc(cl->len * sizeof(int));
    int *suffix_ones = calloc(cl->len + 1, sizeof(int));
    double *current_bound = malloc(N * sizeof(double));

    if (!current_bound || !combination || !suffix_ones || !current_sequence || !remaining_bits) {
        printf("ERROR: No se pudo asignar memoria\n");
        free(current_dft);
        free(current_sequence);
        free(remaining_bits);
        free(combination);
        free(suffix_ones);
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
        }
        for (int j = 0; j < spectrum_size; j++) {
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
    for (int i = 0; i < N; i++) {
        remaining_bits[i] = 0;
    }
    for (int i = 0; i < cl->len; i++) {
        Coset *c = &cl->data[i];
        for (int j = 0; j < c->len; j++) {
            int elem = c->data[j];
            if (elem >= 0 && elem < N) {
                remaining_bits[elem] = 1;
            }
        }
    }

    suffix_ones[cl->len] = 0;
    for (int i = cl->len - 1; i >= 0; i--) {
        suffix_ones[i] = suffix_ones[i + 1] + cl->data[i].len;
    }

    // Inicializar current_bound con la suma de valores absolutos de la DFT para cada frecuencia
    for (int j = 0; j < spectrum_size; j++) {
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

    compresion_1 = calloc(cols1 ? cols1 : 1, sizeof(int));
    compresion_2 = calloc(cols2 ? cols2 : 1, sizeof(int));
    bound_compression_1 = calloc(cols1 ? cols1 : 1, sizeof(int));
    bound_compression_2 = calloc(cols2 ? cols2 : 1, sizeof(int));

    int max_value1 = (cols1 > 0) ? (N / (int)cols1) : 0;
    int max_value2 = (cols2 > 0) ? (N / (int)cols2) : 0;
    uint64_t ***ge_bitsets1 = NULL;
    uint64_t ***le_bitsets1 = NULL;
    uint64_t ***ge_bitsets2 = NULL;
    uint64_t ***le_bitsets2 = NULL;
    size_t bitset_words1 = 0;
    size_t bitset_words2 = 0;
    uint64_t *bitset_work1 = NULL;
    uint64_t *bitset_work2 = NULL;
    int *bitset_order1 = NULL;
    int *bitset_order2 = NULL;
    int *bitset_score1 = NULL;
    int *bitset_score2 = NULL;
    if (!f_comb || !f_psd || !f_cte || !candidate1 || !candidate2 ||
        !compresion_1 || !compresion_2 || !bound_compression_1 || !bound_compression_2 ||
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

        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
        }
        free(dft_matrix);
        free(psd_matrix);
        free(current_dft);
        free(current_sequence);
        free(remaining_bits);
        free(compresion_1);
        free(compresion_2);
        free(bound_compression_1);
        free(bound_compression_2);
        free_candidate_bitsets(ge_bitsets1, le_bitsets1, cols1, max_value1);
        free_candidate_bitsets(ge_bitsets2, le_bitsets2, cols2, max_value2);
        free(bitset_work1);
        free(bitset_work2);
        free(bitset_order1);
        free(bitset_order2);
        free(bitset_score1);
        free(bitset_score2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(compression_a);
        free(compression_b);
        return;
    }

    if (!build_candidate_bitsets(candidate1, rows1, cols1, max_value1,
                                 &ge_bitsets1, &le_bitsets1, &bitset_words1) ||
        !build_candidate_bitsets(candidate2, rows2, cols2, max_value2,
                                 &ge_bitsets2, &le_bitsets2, &bitset_words2)) {
        printf("ERROR: No se pudo construir indice de bitsets para candidatos.\n");
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);
        free(current_bound);
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
        }
        free(dft_matrix);
        free(psd_matrix);
        free(current_dft);
        free(current_sequence);
        free(remaining_bits);
        free(compresion_1);
        free(compresion_2);
        free(bound_compression_1);
        free(bound_compression_2);
        free_candidate_bitsets(ge_bitsets1, le_bitsets1, cols1, max_value1);
        free_candidate_bitsets(ge_bitsets2, le_bitsets2, cols2, max_value2);
        free(bitset_work1);
        free(bitset_work2);
        free(bitset_order1);
        free(bitset_order2);
        free(bitset_score1);
        free(bitset_score2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(compression_a);
        free(compression_b);
        return;
    }

    bitset_work1 = calloc(bitset_words1, sizeof(uint64_t));
    bitset_work2 = calloc(bitset_words2, sizeof(uint64_t));
    bitset_order1 = malloc(cols1 * sizeof(int));
    bitset_order2 = malloc(cols2 * sizeof(int));
    bitset_score1 = malloc(cols1 * sizeof(int));
    bitset_score2 = malloc(cols2 * sizeof(int));
    if (!bitset_work1 || !bitset_work2 || !bitset_order1 || !bitset_order2 || !bitset_score1 || !bitset_score2) {
        printf("ERROR: No se pudo asignar memoria para bitset_work.\n");
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);
        free(current_bound);
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
        }
        free(dft_matrix);
        free(psd_matrix);
        free(current_dft);
        free(current_sequence);
        free(remaining_bits);
        free(compresion_1);
        free(compresion_2);
        free(bound_compression_1);
        free(bound_compression_2);
        free_candidate_bitsets(ge_bitsets1, le_bitsets1, cols1, max_value1);
        free_candidate_bitsets(ge_bitsets2, le_bitsets2, cols2, max_value2);
        free(bitset_work1);
        free(bitset_work2);
        free(bitset_order1);
        free(bitset_order2);
        free(bitset_score1);
        free(bitset_score2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(compression_a);
        free(compression_b);
        return;
    }

    for (int i = 0; i < N; i++) {
        if (remaining_bits[i]) {
            bound_compression_1[i % (int)cols1] += 1;
            bound_compression_2[i % (int)cols2] += 1;
        }
    }

    combination[0] = -1;
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
        .bound_psd = bound_psd,
        .current_sequence = current_sequence,
        .remaining_bits = remaining_bits,
        .compresion_1 = compresion_1,
        .compresion_2 = compresion_2,
        .bound_compression_1 = bound_compression_1,
        .bound_compression_2 = bound_compression_2,
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
        .dimension_candidate_pairs2 = cols2,
        .current_bound = current_bound,
        .depth_exploration_bound = NULL,
        .pos_depth = 0,
        .current_ones = 0,
        .target_ones = (N + 1) / 2,
        .spectrum_size = spectrum_size,
        .suffix_ones = suffix_ones,
        .ge_bitsets1 = ge_bitsets1,
        .le_bitsets1 = le_bitsets1,
        .bitset_words1 = bitset_words1,
        .bitset_max_value1 = max_value1,
        .bitset_work1 = bitset_work1,
        .bitset_order1 = bitset_order1,
        .bitset_score1 = bitset_score1,
        .ge_bitsets2 = ge_bitsets2,
        .le_bitsets2 = le_bitsets2,
        .bitset_words2 = bitset_words2,
        .bitset_max_value2 = max_value2,
        .bitset_work2 = bitset_work2,
        .bitset_order2 = bitset_order2,
        .bitset_score2 = bitset_score2
    };
    reorder_cosets_by_candidate_pairs(&ctx);
    
    // Actualizar punteros locales después de reordenación
    dft_matrix = ctx.dft_matrix;
    psd_matrix = ctx.psd_matrix;
    
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
        free(current_sequence);
        free(remaining_bits);
        free(compresion_1);
        free(compresion_2);
        free(bound_compression_1);
        free(bound_compression_2);
        free_candidate_bitsets(ge_bitsets1, le_bitsets1, cols1, max_value1);
        free_candidate_bitsets(ge_bitsets2, le_bitsets2, cols2, max_value2);
        free(bitset_work1);
        free(bitset_work2);
        free(bitset_order1);
        free(bitset_order2);
        free(bitset_score1);
        free(bitset_score2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(compression_a);
        free(compression_b);
        if (candidate1) free_pairs(candidate1, rows1);
        if (candidate2) free_pairs(candidate2, rows2);
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
    ctx.current_ones = 0;

    dfs_explore_combinations(&ctx);

    // Limpieza completa de memoria
    for (int i = 0; i < cl->len; i++) {
        free(dft_matrix[i]);
        free(psd_matrix[i]);
    }
    free(dft_matrix);
    free(psd_matrix);
    free(current_dft);
    free(current_sequence);
    free(remaining_bits);
    free(compresion_1);
    free(compresion_2);
    free(bound_compression_1);
    free(bound_compression_2);
    free_candidate_bitsets(ge_bitsets1, le_bitsets1, cols1, max_value1);
    free_candidate_bitsets(ge_bitsets2, le_bitsets2, cols2, max_value2);
    free(bitset_work1);
    free(bitset_work2);
    free(bitset_order1);
    free(bitset_order2);
    free(bitset_score1);
    free(bitset_score2);
    free(bound_psd);
    free(combination);
    free(suffix_ones);
    free(compression_a);
    free(compression_b);
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