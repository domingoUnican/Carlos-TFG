#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pairs_reader.h"

int** read_pairs_file(const char *filename, size_t *rows, size_t *cols) {
    if (!filename || !rows || !cols) {
        fprintf(stderr, "ERROR: parámetros inválidos\n");
        return NULL;
    }

    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "ERROR: no se pudo abrir archivo '%s'\n", filename);
        return NULL;
    }

    // Primer paso: contar líneas y columnas
    size_t line_count = 0;
    size_t max_cols = 0;
    char buffer[4096];
    
    while (fgets(buffer, sizeof(buffer), file)) {
        // Ignorar líneas vacías
        if (buffer[0] == '\n' || buffer[0] == '\0') {
            continue;
        }
        
        // Contar espacios en blanco para determinar columnas
        size_t col_count = 0;
        int num;
        int offset = 0;
        int n = 0;
        
        // Parsear la línea para contar números
        while (sscanf(buffer + offset, "%d%n", &num, &n) == 1) {
            col_count++;
            offset += n;
        }
        
        if (col_count > 0) {
            line_count++;
            if (col_count > max_cols) {
                max_cols = col_count;
            }
        }
    }

    if (line_count == 0 || max_cols == 0) {
        fprintf(stderr, "ERROR: archivo %s vacío o formato inválido\n", filename);
        fclose(file);
        return NULL;
    }

    *rows = line_count;
    *cols = max_cols;

    // Segundo paso: asignar memoria para el array bidimensional
    int **pairs = (int**)malloc(line_count * sizeof(int*));
    if (!pairs) {
        fprintf(stderr, "ERROR: sin memoria para array de punteros\n");
        fclose(file);
        return NULL;
    }

    for (size_t i = 0; i < line_count; i++) {
        pairs[i] = (int*)malloc(max_cols * sizeof(int));
        if (!pairs[i]) {
            fprintf(stderr, "ERROR: sin memoria para fila %zu\n", i);
            // Limpiar memoria parcialmente asignada
            for (size_t j = 0; j < i; j++) {
                free(pairs[j]);
            }
            free(pairs);
            fclose(file);
            return NULL;
        }
        
        // Inicializar con ceros
        memset(pairs[i], 0, max_cols * sizeof(int));
    }

    // Reiniciar archivo
    rewind(file);

    // Tercer paso: leer los datos
    size_t current_row = 0;
    while (fgets(buffer, sizeof(buffer), file) && current_row < line_count) {
        // Ignorar líneas vacías
        if (buffer[0] == '\n' || buffer[0] == '\0') {
            continue;
        }

        size_t col_idx = 0;
        int num;
        int offset = 0;
        int n = 0;

        // Parsear números de la línea
        while (sscanf(buffer + offset, "%d%n", &num, &n) == 1 && col_idx < max_cols) {
            pairs[current_row][col_idx] = num;
            col_idx++;
            offset += n;
        }

        current_row++;
    }

    fclose(file);
    
    printf("Archivo '%s' leído exitosamente: %zu filas x %zu columnas\n", filename, line_count, max_cols);
    return pairs;
}

void free_pairs(int **pairs, size_t rows) {
    if (!pairs) return;
    
    for (size_t i = 0; i < rows; i++) {
        if (pairs[i]) {
            free(pairs[i]);
        }
    }
    free(pairs);
}

void print_pairs(const int **pairs, size_t rows, size_t cols) {
    if (!pairs) {
        printf("(datos nulos)\n");
        return;
    }

    printf("=== Array de Pares ===\n");
    printf("Dimensiones: %zu x %zu\n\n", rows, cols);
    
    for (size_t i = 0; i < rows; i++) {
        printf("Fila %zu: [", i);
        for (size_t j = 0; j < cols; j++) {
            printf("%d", pairs[i][j]);
            if (j + 1 < cols) printf(", ");
        }
        printf("]\n");
    }
}
