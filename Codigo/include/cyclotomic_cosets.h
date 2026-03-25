
/* cyclotomic_cosets.h */
#ifndef CYCLOTOMIC_COSETS_H
#define CYCLOTOMIC_COSETS_H

#include <stddef.h>  // size_t
#include <stdint.h>  // uint8_t
#include <stdbool.h>  // bool

#ifdef __cplusplus
extern "C" {
#endif

/** Un coset como array dinámico de elementos. */
typedef struct {
    size_t *data;  ///< elementos del coset
    size_t len;    ///< número de elementos
    size_t cap;    ///< capacidad interna (para realloc)
} Coset;

/** Lista de cosets (array dinámico de Coset). */
typedef struct {
    Coset *data;   ///< array de cosets
    size_t len;    ///< número de cosets
    size_t cap;    ///< capacidad interna
    size_t *positions;
} CosetList;

/* API principal que ya tenías */
CosetList cyclotomic_cosets(size_t k, size_t N);
void free_cosetlist(CosetList *cl);
void cosetlist_print(const CosetList *cl);

/**
 * @brief Comprueba si el coset 'a' está contenido en el coset 'b' (como conjuntos).
 *        Ignora el orden, compara por pertenencia de elementos.
 * @param a coset candidato a estar contenido
 * @param b coset supuestamente contenedor
 * @return true si a ⊆ b, false en caso contrario.
 */
bool coset_is_subset_of(const Coset *a, const Coset *b);

/**
 * @brief Devuelve un vector {0,1} de longitud cl->len indicando, para cada coset
 *        de la lista 'cl', si 'c' está contenido en él.
 * 
 * @param c  coset a comprobar (si es NULL o c->len==0, se considera conjunto vacío
 *           y devolverá todo '1')
 * @param cl lista de cosets (si es NULL o cl->len==0, devuelve NULL)
 * @return   puntero a array de uint8_t (0/1) de longitud cl->len. Debes liberar con free().
 */
uint8_t* coset_membership_vector(const Coset *c, const CosetList *cl);

/**
 * @brief Busca el coset que contiene un elemento dado.
 * @param cl lista de cosets donde se realiza la búsqueda.
 * @param element elemento a localizar.
 * @return índice del coset que contiene 'element'; devuelve -1 si no aparece.
 */
int find_element_in_cosets(const CosetList *cl, size_t element);


#ifdef __cplusplus
}
#endif

#endif /* CYCLOTOMIC_COSETS_H */
