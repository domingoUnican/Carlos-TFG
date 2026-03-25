
/* cyclotomic_cosets.c */
#include "cyclotomic_cosets.h"
#include <stdio.h>
#include <stddef.h>  /* para size_t */
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/* --------- utilidades internas ya presentes --------- */

static size_t gcd_size(size_t a, size_t b) {
    while (b != 0) {
        size_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static void coset_init(Coset *c) {
    c->data = NULL;
    c->len = 0;
    c->cap = 0;
}

static bool coset_push(Coset *c, size_t x) {
    if (c->len == c->cap) {
        size_t new_cap = (c->cap ? c->cap * 2 : 8);
        size_t *new_data = (size_t*)realloc(c->data, new_cap * sizeof(size_t));
        if (!new_data) return false;
        c->data = new_data;
        c->cap = new_cap;
    }
    c->data[c->len++] = x;
    return true;
}

static void cosetlist_init(CosetList *cl, size_t N) {
    cl->data = NULL;
    cl->len = 0;
    cl->cap = 0;
    cl->positions = (size_t*)calloc(N, sizeof(size_t));
    if (cl->positions) {
        for (size_t i = 0; i < N; ++i) {
            cl->positions[i] = (size_t)-1;
        }
    }
}

static bool cosetlist_push(CosetList *cl, const Coset *c_in) {
    if (cl->len == cl->cap) {
        size_t new_cap = (cl->cap ? cl->cap * 2 : 8);
        Coset *new_data = (Coset*)realloc(cl->data, new_cap * sizeof(Coset));
        if (!new_data) return false;
        cl->data = new_data;
        cl->cap = new_cap;
    }
    cl->data[cl->len++] = *c_in;
    return true;
}

/* --------- API pública ya presente --------- */

CosetList cyclotomic_cosets(size_t k, size_t N) {
    CosetList res;
    cosetlist_init(&res, N);

    if (N == 0) {
        return res;
    }

    size_t k_mod = k % N;
    size_t g = gcd_size(k_mod, N);
    if (g != 1) {
        fprintf(stderr,
                "Aviso: gcd(k,N) = %zu != 1. Las órbitas pueden ser cadenas que terminan en ciclos.\n",
                g);
    }

    uint8_t *visited = (uint8_t*)calloc(N, sizeof(uint8_t));
    if (!visited) {
        fprintf(stderr, "Error: sin memoria para 'visited'.\n");
        return res;
    }

    for (size_t start = 0; start < N; ++start) {
        if (visited[start]) continue;

        Coset c;
        coset_init(&c);
        size_t coset_idx = res.len;

        size_t x = start;
        while (!visited[x]) {
            if (!coset_push(&c, x)) {
                fprintf(stderr, "Error: sin memoria al ampliar coset.\n");
                free(c.data);
                free(visited);
                free_cosetlist(&res);
                return res;
            }
            visited[x] = 1;
            res.positions[x] = coset_idx;
            x = (x * k_mod) % N;
        }

        if (!cosetlist_push(&res, &c)) {
            fprintf(stderr, "Error: sin memoria al ampliar lista de cosets.\n");
            free(c.data);
            free(visited);
            free_cosetlist(&res);
            return res;
        }
    }

    free(visited);
    return res;
}

void free_cosetlist(CosetList *cl) {
    if (!cl || !cl->data) {
        if (cl) { cl->len = cl->cap = 0; }
        return;
    }
    for (size_t i = 0; i < cl->len; ++i) {
        free(cl->data[i].data);
        cl->data[i].data = NULL;
        cl->data[i].len = cl->data[i].cap = 0;
    }
    free(cl->data);
    free(cl->positions);
    cl->data = NULL;
    cl->positions = NULL;
    cl->len = cl->cap = 0;
}

void cosetlist_print(const CosetList *cl) {
    if (!cl) { printf("(lista nula)\n"); return; }
    printf("Se encontraron %zu cosets:\n", cl->len);
    for (size_t i = 0; i < cl->len; ++i) {
        const Coset *c = &cl->data[i];
        printf("C%zu = { ", i);
        for (size_t j = 0; j < c->len; ++j) {
            printf("%zu", c->data[j]);
            if (j + 1 < c->len) printf(", ");
        }
        printf(" }\n");
    }
}

int find_element_in_cosets(const CosetList *cl, size_t element)
{
    if (!cl || !cl->data) return -1;
    return (cl->positions[element] == (size_t)-1) ? -1 : (int)cl->positions[element];
}

/* --------- NUEVAS FUNCIONES --------- */

/* Comprueba si 'a' está contenido en 'b' (a ⊆ b). Complejidad O(|a|·|b|). */
bool coset_is_subset_of(const Coset *a, const Coset *b) {
    if (!a || a->len == 0) {
        /* El conjunto vacío está contenido en cualquier conjunto. */
        return true;
    }
    if (!b || b->len == 0) {
        /* a no vacío no puede estar contenido en b vacío. */
        return false;
    }
    for (size_t i = 0; i < a->len; ++i) {
        size_t x = a->data[i];
        bool found = false;
        for (size_t j = 0; j < b->len; ++j) {
            if (b->data[j] == x) { found = true; break; }
        }
        if (!found) return false;
    }
    return true;
}

/* Devuelve vector de 0/1 indicando si 'c' está contenido en cada coset de 'cl'. */
uint8_t* coset_membership_vector(const Coset *c, const CosetList *cl) {
    if (!cl || cl->len == 0) return NULL;

    uint8_t *vec = (uint8_t*)malloc(cl->len * sizeof(uint8_t));
    if (!vec) {
        fprintf(stderr, "Error: sin memoria para vector de pertenencia.\n");
        return NULL;
    }

    if (!c || c->len == 0) {
        /* Si c es vacío, está contenido en todos: rellena con 1 */
        for (size_t i = 0; i < cl->len; ++i) vec[i] = 1;
        return vec;
    }

    for (size_t i = 0; i < cl->len; ++i) {
        vec[i] = coset_is_subset_of(c, &cl->data[i]) ? 1 : 0;
    }
    return vec;
}
