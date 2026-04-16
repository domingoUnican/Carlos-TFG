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
#define DEPTH 22
#define DIMENSIONS 6
// ...existing code...

void show_specific_lines(size_t line1, size_t line2, FILE *lp_file, size_t N, const char *path_combinations) {
    FILE *file = fopen(path_combinations, "r");
    if (!file) {
        printf("ERROR: No se pudo abrir '%s'\n", path_combinations);
        return;
    }

    char *line = malloc(N + 8);
    if (!line) { fclose(file); return; }

    size_t line_idx = 0;
    bool got_line1 = false;
    bool got_line2 = false;
    while (fgets(line, (int)(N + 8), file) != NULL) {
        if (line_idx == line1 || line_idx == line2) {
            size_t len = strlen(line);
            while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
                line[--len] = '\0';
            }
            fprintf(lp_file, "%s\n", line);
            if (line_idx == line1) got_line1 = true;
            if (line_idx == line2) got_line2 = true;
            if (got_line1 && got_line2) {
                break;
            }
        }
        line_idx++;
    }
    fprintf(lp_file, "\n");

    free(line);
    fclose(file);
}


void find_matches_files(size_t N, size_t row_ints, char* path_psd, char* path_cte, char* path_lp, char* path_comb) {
    FILE *f1 = fopen(path_psd, "rb");
    FILE *f2 = fopen(path_cte, "rb");
    FILE *lp_file = fopen(path_lp, "a"); // acumulamos los LP encontrados en cada llamada
    
    if (!f1 || !f2 || !lp_file) {
        if (f1) fclose(f1); if (f2) fclose(f2);
        return;
    }

    size_t L_dft = row_ints; // Cantidad de enteros por bloque binario
    size_t block_size = L_dft * sizeof(int);
    // 1. Cargar fichero 2 en memoria para evitar accesos a disco
    size_t capacity = 1000;
    size_t total_lines2 = 0;

    int *all_lines2 = malloc(capacity * block_size);

    int *temp_buffer = malloc(block_size);

    while (fread(temp_buffer, sizeof(int), L_dft, f2) == L_dft) {
        if (total_lines2 >= capacity) {
            capacity *= 2;
            all_lines2 = realloc(all_lines2, capacity * block_size);
        }
        // Copiamos el bloque al gran buffer de memoria
        memcpy(&all_lines2[total_lines2 * L_dft], temp_buffer, block_size);
        total_lines2++;
    }

    printf("Searching for LPS...\n");
    bool LP_found = false; // flag para parar si encuentra el primer LP
    int *line1_data = malloc(block_size);
    size_t line1_num = 0;

    // 2. Leer F1 bloque a bloque y comparar contra todo F2 en memoria
    rewind(f1);
    while (!LP_found && fread(line1_data, sizeof(int), L_dft, f1) == L_dft) {
        
        size_t j = 0;
        while (j < total_lines2 && !LP_found) {
            
            // Condición: i < j para evitar duplicados y simetrías
            if (line1_num < j) {
                // Comparamos el bloque actual de f1 con el bloque j en RAM de f2
                if (memcmp(line1_data, &all_lines2[j * L_dft], block_size) == 0) {
                    
                    // Si coinciden los espectros, es un Legendre Pair
                    show_specific_lines(line1_num, j, lp_file, N, path_comb);
                    
                    printf("LEGENDRE PAIR! Indices: %zu y %zu\n", line1_num, j);
                    LP_found = true; 
                    // El !LP_found en la condición del while romperá el bucle j automáticamente
                }
            }
            j++;
        }
        line1_num++;
    }

    // Limpieza
    free(all_lines2);
    free(temp_buffer);
    free(line1_data);
    fclose(f1); fclose(f2); fclose(lp_file);
}

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
    int dims_used;
    int max_value;
    size_t axis_len;
    size_t total_cells;
    size_t strides[DIMENSIONS];
    uint32_t *prefix;
} PrefixRangeIndex;

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
    int **psd_matrix;
    int N;
    int num_cosets;
    double complex *current_dft;
    int *bound_psd;
    int *current_sequence;
    int *current_psd; 
    int *remaining_bits;
    int *compresion_1;
    int *compresion_2;
    int *bound_compression_1;
    int *bound_compression_2;
    int *compresion_3;
    int *bound_contrassion_3;
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
    bool *is_less_compression_a;
    bool *is_less_compression_b;
    bool **is_less_candidates_1;
    bool **is_less_candidates_2;
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
    PrefixRangeIndex prefix_index1;
    PrefixRangeIndex prefix_index2;
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

static void free_prefix_range_index(PrefixRangeIndex *index) {
    if (!index) return;
    free(index->prefix);
    index->prefix = NULL;
    index->dims_used = 0;
    index->max_value = 0;
    index->axis_len = 0;
    index->total_cells = 0;
    for (int i = 0; i < DIMENSIONS; i++) {
        index->strides[i] = 0;
    }
}

static bool **alloc_bool_matrix_true(size_t rows, size_t cols) {
    bool **matrix = malloc(rows * sizeof(bool *));
    if (!matrix) {
        return NULL;
    }

    for (size_t i = 0; i < rows; i++) {
        matrix[i] = malloc(cols * sizeof(bool));
        if (!matrix[i]) {
            for (size_t k = 0; k < i; k++) {
                free(matrix[k]);
            }
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

static bool build_prefix_range_index(const int *const *pairs, size_t num_pairs,
                                     size_t dim, int max_value,
                                     PrefixRangeIndex *out) {
    if (!out) return false;

    out->dims_used = 0;
    out->max_value = 0;
    out->axis_len = 0;
    out->total_cells = 0;
    out->prefix = NULL;
    for (int i = 0; i < DIMENSIONS; i++) {
        out->strides[i] = 0;
    }

    if (!pairs || dim == 0 || max_value < 0) {
        fprintf(stderr, "DEBUG: build_prefix_range_index falló: pairs=%p, dim=%zu, max_value=%d\n", 
                (void*)pairs, dim, max_value);
        return false;
    }

    int dims_used = (int)dim;
    if (dims_used > DIMENSIONS) {
        dims_used = DIMENSIONS;
    }
    if (dims_used <= 0) {
        fprintf(stderr, "DEBUG: dims_used=%d es inválido\n", dims_used);
        return false;
    }

    size_t axis_len = (size_t)max_value + 1u;
    size_t total_cells = 1u;
    const size_t max_cells = 20000000u;
    for (int d = 0; d < dims_used; d++) {
        if (axis_len != 0 && total_cells > SIZE_MAX / axis_len) {
            fprintf(stderr, "DEBUG: overflow en stride computation (d=%d, axis_len=%zu, total_cells=%zu)\n",
                    d, axis_len, total_cells);
            return false;
        }
        total_cells *= axis_len;
        if (total_cells > max_cells) {
            fprintf(stderr, "DEBUG: total_cells=%zu excede max_cells=%zu (dims_used=%d, axis_len=%zu, max_value=%d)\n",
                    total_cells, max_cells, dims_used, axis_len, max_value);
            return false;
        }
    }

    uint32_t *grid = calloc(total_cells, sizeof(uint32_t));
    if (!grid) {
        fprintf(stderr, "DEBUG: malloc falló para total_cells=%zu\n", total_cells);
        return false;
    }

    out->strides[0] = 1u;
    for (int d = 1; d < dims_used; d++) {
        out->strides[d] = out->strides[d - 1] * axis_len;
    }

    for (size_t i = 0; i < num_pairs; i++) {
        size_t idx = 0u;
        for (int d = 0; d < dims_used; d++) {
            int value = pairs[i][d];
            if (value < 0) value = 0;
            if (value > max_value) value = max_value;
            idx += (size_t)value * out->strides[d];
        }
        grid[idx] += 1u;
    }

    for (int d = 0; d < dims_used; d++) {
        size_t step = out->strides[d];
        size_t block = step * axis_len;
        for (size_t base = 0; base < total_cells; base += block) {
            for (size_t off = 0; off < step; off++) {
                uint32_t running = 0u;
                for (size_t pos = 0; pos < axis_len; pos++) {
                    size_t idx = base + off + pos * step;
                    running += grid[idx];
                    grid[idx] = running;
                }
            }
        }
    }

    out->dims_used = dims_used;
    out->max_value = max_value;
    out->axis_len = axis_len;
    out->total_cells = total_cells;
    out->prefix = grid;
    return true;
}

static bool exists_candidate_in_prefix_range(const PrefixRangeIndex *index,
                                             const int *compression,
                                             const int *compression_bound) {
    if (!index || !index->prefix || !compression || !compression_bound || index->dims_used <= 0) {
        return true;
    }

    int lo[DIMENSIONS];
    int hi[DIMENSIONS];
    for (int d = 0; d < index->dims_used; d++) {
        int l = compression[d];
        int u = compression[d] + compression_bound[d];
        if (u < 0 || l > index->max_value) {
            return false;
        }
        if (l < 0) l = 0;
        if (u > index->max_value) u = index->max_value;
        if (l > u) {
            return false;
        }
        lo[d] = l;
        hi[d] = u;
    }

    long long count = 0;
    int masks = 1 << index->dims_used;
    for (int mask = 0; mask < masks; mask++) {
        size_t idx = 0u;
        bool valid = true;
        for (int d = 0; d < index->dims_used; d++) {
            int coord = (mask & (1 << d)) ? (lo[d] - 1) : hi[d];
            if (coord < 0) {
                valid = false;
                break;
            }
            idx += (size_t)coord * index->strides[d];
        }
        if (!valid) {
            continue;
        }

        long long term = (long long)index->prefix[idx];
        if (__builtin_popcount((unsigned int)mask) & 1) {
            count -= term;
        } else {
            count += term;
        }
    }

    return count > 0;
}

static bool exists_candidate_in_range_bitset_precomputed(const PrefixRangeIndex *prefix_index,
                                                         const int *compression,
                                                         const int *compression_bound,
                                                         size_t dim,
                                                         uint64_t ***ge, uint64_t ***le,
                                                         size_t words_len, size_t num_pairs, int max_value,
                                                         uint64_t *work,
                                                         int *order_buf,
                                                         int *score_buf) {
    if (!exists_candidate_in_prefix_range(prefix_index, compression, compression_bound)){
       return false;
    }

    return exists_candidate_in_range_bitset(compression, compression_bound,
                                            dim,
                                            ge, le,
                                            words_len, num_pairs, max_value,
                                            work,
                                            order_buf,
                                            score_buf);
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
    PrefixRangeIndex *prefix_index = (flag == 1) ? &ctx->prefix_index1 : &ctx->prefix_index2;
    if (!compression || !compression_bound) {
        free(compression);
        free(compression_bound);
        return false;
    }

    CompressSequence(ctx->N, dimension_candidate_pairs, sequence, compression);
    CompressSequence(ctx->N, dimension_candidate_pairs, bound_bit_sequence, compression_bound);

    bool resultado = exists_candidate_in_range_bitset_precomputed(prefix_index,
                                                                  compression, compression_bound,
                                                                  (size_t)dimension_candidate_pairs,
                                                                  ge, le, words_len, num_pairs, max_value,
                                                                  work, order_buf, score_buf);
    free(compression);
    free(compression_bound);
    return resultado;
}

bool is_less_than_candidates(DFSContext *ctx, int flag)
{
    if (!ctx || !ctx->compresion_3 || !ctx->bound_contrassion_3 || !ctx->compression_a || !ctx->compression_b || ctx->p <= 0) {
        return false;
    }
    Coset *selected_coset = &ctx->cl->data[ctx->coset_idx - 1];
    if (flag == 1){
        for (int t = 0; t < selected_coset->len; t++) {
            int elem = selected_coset->data[t];
            int l = ctx->compresion_3[elem % ctx->p];
            int u = l + ctx->bound_contrassion_3[elem % ctx->p];
            int ta = ctx->compression_a[elem % ctx->p];
            int tb = ctx->compression_b[elem % ctx->p];
            // Miramos la compresión anterior a la actual, que es la que se va a comparar con los candidatos
            bool p_a = ctx->is_less_compression_a[ctx->coset_idx - 1]; 
            bool p_b = ctx->is_less_compression_b[ctx->coset_idx - 1];
            ctx->is_less_compression_a[ctx->coset_idx] = p_a && (l <= ta && ta <= u);
            ctx->is_less_compression_b[ctx->coset_idx] = p_b && (l <= tb && tb <= u);

        }
        if (! ctx->is_less_compression_a[ctx->coset_idx] && ! ctx->is_less_compression_b[ctx->coset_idx]) {
            return false;
        }
    }

    int dimension_candidate_pairs = (flag == 1) ? ctx->dimension_candidate_pairs1 : ctx->dimension_candidate_pairs2;
    size_t num_pairs = (flag == 1) ? ctx->num_candidate_pairs1 : ctx->num_candidate_pairs2;
    int *compression = (flag == 1) ? ctx->compresion_1 : ctx->compresion_2;
    int *compression_bound = (flag == 1) ? ctx->bound_compression_1 : ctx->bound_compression_2;
    int ** candidate_pairs = (flag == 1) ? ctx->candidate_pairs1 : ctx->candidate_pairs2;
    bool ** is_less_cand = (flag == 1) ? ctx->is_less_candidates_1 : ctx->is_less_candidates_2;
    bool exists = false;
    for (int t = 0; t < selected_coset->len; t++) {
        int elem = selected_coset->data[t];
        int l = compression[elem % dimension_candidate_pairs];
        int u = l + compression_bound[elem % dimension_candidate_pairs];
        for (size_t j = 0; j < num_pairs; j++) {
            int v = candidate_pairs[j][elem % dimension_candidate_pairs];
            is_less_cand[ctx->coset_idx][j] = is_less_cand[ctx->coset_idx - 1][j] && (l <= v && v <= u);
            exists = exists || is_less_cand[ctx->coset_idx][j];
        }

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

bool check_bound(DFSContext *ctx) {
    bool result = is_less_than_candidates(ctx, 1) &&
                  is_less_than_candidates(ctx, 0);
    if (!result) {
        return false;
    }
    double max_lower_bound = 0.0;
    for (int j = 1; j < ctx->spectrum_size; j++) {
        double current_abs = cabs(ctx->current_dft[j]);
        ctx->current_psd[j] = (int)rint(pow(current_abs, 2));
        double lower_bound_psd = pow(fmax(0.0, current_abs - ctx->current_bound[j]), 2.0);
        if (lower_bound_psd > max_lower_bound) {
            max_lower_bound = lower_bound_psd;
        }
    }
    return max_lower_bound <= ctx->threshold;
}

bool is_valid_combination(const DFSContext *ctx) {
    //int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);

    int *temp= ctx->current_sequence; //convert_to_binary_vector(vector_bits, ctx->N);
    //free(vector_bits);
    int temp_psd = -1;
    bool result = false;
    if (is_compression_of_candidates(ctx,temp,0) && is_compression_of_candidates(ctx,temp,1))
    {
        result = true;
        for (int j = 1; (j < ctx->spectrum_size) && result; j++) {
            double real_part = creal(ctx->current_dft[j]);
            double imag_part = cimag(ctx->current_dft[j]);
            temp_psd = (int)(real_part * real_part + imag_part * imag_part);
            ctx->current_psd[j] = temp_psd;
            result = (temp_psd <= ctx->threshold);
        }
        
    }
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
            for (int j = 0; j < ctx->N; j++) {
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
                int transformed = ctx->threshold - psd_val;

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
                
                /*
                 * No podar por PSD del sufijo aislado: aunque prefijo y sufijo
                 * no compartan unos en tiempo, en frecuencia si puede haber
                 * cancelacion entre ambos. Este bound debe seguir siendo
                 * conservador respecto al prefijo actual.
                 */
                bool passes = is_less_than_candidates_for_vectors(ctx, seq, bound_seq, 1) &&
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

    int threshold = (N+1)/2;
    int spectrum_size = (N  + 1)/2;
    //int spectrum_size = 20;
    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    double complex **dft_backup = malloc(cl->len * sizeof(double complex *));
    int **psd_matrix = malloc(cl->len * sizeof(int *));

    for (int i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        psd_matrix[i] = malloc(N * sizeof(int));
        dft_backup[i] = malloc(N * sizeof(double complex));
    }

    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    int *bound_psd = calloc(N, sizeof(int));
    int *compression_a = malloc(p * sizeof(int));
    int *compression_b = malloc(p * sizeof(int));
    double complex *current_dft = calloc(N, sizeof(double complex));
    int *current_sequence = calloc((size_t)N, sizeof(int));
    int *current_psd = calloc((size_t)N, sizeof(int));
    int *remaining_bits = calloc((size_t)N, sizeof(int));
    int *compresion_1 = NULL;
    int *compresion_2 = NULL;
    int *bound_compression_1 = NULL;
    int *bound_compression_2 = NULL;
    int *compresion_3 = NULL;
    int *bound_contrassion_3 = NULL;
    int *combination = malloc(cl->len * sizeof(int));
    int *suffix_ones = calloc(cl->len + 1, sizeof(int));
    bool *is_less_compression_a = malloc((size_t)(cl->len + 1) * sizeof(bool));
    bool *is_less_compression_b = malloc((size_t)(cl->len + 1) * sizeof(bool));
    double *current_bound = malloc(N * sizeof(double));

    if (!current_bound || !combination || !suffix_ones || !current_sequence || !current_psd || !remaining_bits ||
        !is_less_compression_a || !is_less_compression_b) {
        printf("ERROR: No se pudo asignar memoria\n");
        free(current_dft);
        free(current_sequence);
        free(current_psd);
        free(remaining_bits);
        free(combination);
        free(suffix_ones);
        free(is_less_compression_a);
        free(is_less_compression_b);
        free(compression_a);
        free(compression_b);
        free(compresion_3);
        free(bound_contrassion_3);
        for (int i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
            free(dft_backup[i]);    
        }
        free(dft_matrix);
        free(psd_matrix);
        free(dft_backup);
        return;
    }

    /* Inicializar combination como lista vacía terminada en -1 */
    for (int i = 0; i < cl->len; i++) {
        combination[i] = -1;
    }
    for (int i = 0; i <= cl->len; i++) {
        is_less_compression_a[i] = true;
        is_less_compression_b[i] = true;
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
            dft_backup[i][j] = 0.0 + 0.0 * I;
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
        current_psd[j] = 0;
        for (int i = 0; i < cl->len; i++) {
            current_bound[j] += cabs(dft_matrix[i][j]);
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
    compresion_3 = calloc((size_t)p, sizeof(int));
    bound_contrassion_3 = calloc((size_t)p, sizeof(int));

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
    bool **is_less_candidates_1 = NULL;
    bool **is_less_candidates_2 = NULL;
    PrefixRangeIndex prefix_index1 = {0};
    PrefixRangeIndex prefix_index2 = {0};
    if (!f_comb || !f_psd || !f_cte || !candidate1 || !candidate2 ||
        !compresion_1 || !compresion_2 || !bound_compression_1 || !bound_compression_2 ||
        !compresion_3 || !bound_contrassion_3 ||
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
        free_prefix_range_index(&prefix_index1);
        free_prefix_range_index(&prefix_index2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(is_less_compression_a);
        free(is_less_compression_b);
        free(compression_a);
        free(compression_b);
        free(compresion_3);
        free(bound_contrassion_3);
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
        free_prefix_range_index(&prefix_index1);
        free_prefix_range_index(&prefix_index2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(is_less_compression_a);
        free(is_less_compression_b);
        free(compression_a);
        free(compression_b);
        free(compresion_3);
        free(bound_contrassion_3);
        return;
    }

    is_less_candidates_1 = alloc_bool_matrix_true((size_t)(cl->len + 1), rows1);
    is_less_candidates_2 = alloc_bool_matrix_true((size_t)(cl->len + 1), rows2);
    if (!is_less_candidates_1 || !is_less_candidates_2) {
        printf("ERROR: No se pudo asignar memoria para is_less_candidates_1/2.\n");
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
        free_prefix_range_index(&prefix_index1);
        free_prefix_range_index(&prefix_index2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(is_less_compression_a);
        free(is_less_compression_b);
        free(compression_a);
        free(compression_b);
        free(compresion_3);
        free(bound_contrassion_3);
        free_bool_matrix(is_less_candidates_1, (size_t)(cl->len + 1));
        free_bool_matrix(is_less_candidates_2, (size_t)(cl->len + 1));
        return;
    }

    if (!build_prefix_range_index((const int *const *)candidate1, rows1, cols1, max_value1, &prefix_index1) ||
        !build_prefix_range_index((const int *const *)candidate2, rows2, cols2, max_value2, &prefix_index2)) {
        printf("ADVERTENCIA: No se pudo construir prefiltro para primeras dimensiones; se usa chequeo completo.\n");
        free_prefix_range_index(&prefix_index1);
        free_prefix_range_index(&prefix_index2);
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
        free_prefix_range_index(&prefix_index1);
        free_prefix_range_index(&prefix_index2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(is_less_compression_a);
        free(is_less_compression_b);
        free_bool_matrix(is_less_candidates_1, (size_t)(cl->len + 1));
        free_bool_matrix(is_less_candidates_2, (size_t)(cl->len + 1));
        free(compression_a);
        free(compression_b);
        free(compresion_3);
        free(bound_contrassion_3);
        return;
    }

    for (int i = 0; i < N; i++) {
        if (remaining_bits[i]) {
            bound_compression_1[i % (int)cols1] += 1;
            bound_compression_2[i % (int)cols2] += 1;
            bound_contrassion_3[i % p] += 1;
        }
    }

    combination[0] = -1;
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
        .psd_matrix = psd_matrix,
        .current_dft = current_dft,
        .bound_psd = bound_psd,
        .current_sequence = current_sequence,
        .remaining_bits = remaining_bits,
        .compresion_1 = compresion_1,
        .compresion_2 = compresion_2,
        .bound_compression_1 = bound_compression_1,
        .bound_compression_2 = bound_compression_2,
        .compresion_3 = compresion_3,
        .bound_contrassion_3 = bound_contrassion_3,
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
        .current_psd = current_psd,
        .is_less_compression_a = is_less_compression_a,
        .is_less_compression_b = is_less_compression_b,
        .is_less_candidates_1 = is_less_candidates_1,
        .is_less_candidates_2 = is_less_candidates_2,
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
        .bitset_score2 = bitset_score2,
        .prefix_index1 = prefix_index1,
        .prefix_index2 = prefix_index2
    };
    snprintf(ctx.comb_filename, sizeof(ctx.comb_filename), "%s", comb_filename);
    snprintf(ctx.psd_filename, sizeof(ctx.psd_filename), "%s", psd_filename);
    snprintf(ctx.cte_filename, sizeof(ctx.cte_filename), "%s", cte_filename);
    snprintf(ctx.lp_filename, sizeof(ctx.lp_filename), "%s", lp_filename);
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
        free_prefix_range_index(&ctx.prefix_index1);
        free_prefix_range_index(&ctx.prefix_index2);
        free(bound_psd);
        free(combination);
        free(suffix_ones);
        free(is_less_compression_a);
        free(is_less_compression_b);
        free_bool_matrix(is_less_candidates_1, (size_t)(cl->len + 1));
        free_bool_matrix(is_less_candidates_2, (size_t)(cl->len + 1));
        free(compression_a);
        free(compression_b);
        free(compresion_3);
        free(bound_contrassion_3);
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
    for (int i = 0; i < cl->len; i++) combination[i] = -1;
    ctx.coset_idx = 0;
    ctx.pos_depth = 0;
    ctx.current_ones = 0;

    dfs_explore_combinations(&ctx);

    // Limpieza completa de memoria
    for (int i = 0; i < cl->len; i++) {
        free(dft_matrix[i]);
        free(psd_matrix[i]);
        free(dft_backup[i]);
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
    free_prefix_range_index(&ctx.prefix_index1);
    free_prefix_range_index(&ctx.prefix_index2);
    free(bound_psd);
    free(combination);
    free(suffix_ones);
    free(is_less_compression_a);
    free(is_less_compression_b);
    free_bool_matrix(is_less_candidates_1, (size_t)(cl->len + 1));
    free_bool_matrix(is_less_candidates_2, (size_t)(cl->len + 1));
    free(compression_a);
    free(compression_b);
    free(compresion_3);
    free(bound_contrassion_3);
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
