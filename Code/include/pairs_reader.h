#ifndef PAIRS_READER_H
#define PAIRS_READER_H

#include <stddef.h>

/**
 * Lee un archivo de pares y devuelve un array de arrays
 * 
 * @param filename Nombre del archivo a leer (ej: "99_pairs.txt")
 * @param rows Puntero donde se almacenará el número de filas (pares)
 * @param cols Puntero donde se almacenará el número de columnas
 * 
 * @return Puntero a array bidimensional (int**), o NULL si hay error
 *         Debe ser liberado con free_pairs() cuando ya no se necesite
 */
int** read_pairs_file(const char *filename, size_t *rows, size_t *cols);

/**
 * Libera la memoria asignada por read_pairs_file()
 * 
 * @param pairs Puntero al array bidimensional
 * @param rows Número de filas
 */
void free_pairs(int **pairs, size_t rows);

/**
 * Imprime el contenido del array de pares
 * 
 * @param pairs Puntero al array bidimensional
 * @param rows Número de filas
 * @param cols Número de columnas
 */
void print_pairs(const int **pairs, size_t rows, size_t cols);

#endif // PAIRS_READER_H
