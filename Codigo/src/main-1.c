
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cyclotomic_cosets.h"
#include "coset_vector.h"
#include "mymath.h"
#include <complex.h>
#include <math.h>
#include <time.h>

const char *PATH_COMP_COSETS            = "comp_cte_cosets.bin";            // Secuencias originales filtradas
const char *PATH_COMP_DFT_COSETS        = "comp_dft_cosets.bin";            // Magnitudes de la DFT
const char *PATH_COMP_LP                = "comp_lp.bin";                    // Pares encontrados comprimidos

const char *PATH_COMBINATIONS           = "combinations.bin";               // Combinaciones descomprimidas
const char *PATH_DFT                    = "dft.bin";                        // DFT de las combinaciones
const char *PATH_CTEDFT                 = "cte-dft.bin";                    // Constante - DFT de las combinaciones

const char *PATH_LP                     = "lp.txt";                         // LPS




void show_specific_lines(size_t line1, size_t line2, FILE *lp_file, size_t N) {
    FILE *file = fopen(PATH_COMBINATIONS, "rb");
    if (!file) {
        printf("ERROR: No se pudo abrir '%s'\n", PATH_COMBINATIONS);
        return;
    }

    // Buffer para una secuencia (N bytes)
    unsigned char *seq = malloc(N);
    if (!seq) { fclose(file); return; }

    fseek(file, (long)(line1 * N), SEEK_SET);
    if (fread(seq, 1, N, file) == N) {
        for (size_t i = 0; i < N; i++) {
            fprintf(lp_file, "%d", seq[i]);
        }
        fprintf(lp_file, "\n");
    }

    // --- EXTRAER LÍNEA 2 ---
    fseek(file, (long)(line2 * N), SEEK_SET);
    if (fread(seq, 1, N, file) == N) {
        for (size_t i = 0; i < N; i++) {
            fprintf(lp_file, "%d", seq[i]);
        }
        fprintf(lp_file, "\n\n");
    }

    free(seq);
    fclose(file);
}

typedef struct {
    char *data;
    size_t original_index;
} DFTLine;

void find_matches_files(size_t N) {
    FILE *f1 = fopen(PATH_DFT, "rb");
    FILE *f2 = fopen(PATH_CTEDFT, "rb");
    FILE *lp_file = fopen(PATH_LP, "w"); //este si a texto para que podamos leerlo
    
    if (!f1 || !f2 || !lp_file) {
        if (f1) fclose(f1); if (f2) fclose(f2);
        return;
    }

    size_t L_dft = N - 1; // Cantidad de enteros por bloque
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
                    show_specific_lines(line1_num, j, lp_file, N);
                    
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



bool is_constant_on_cosets(int *x, const CosetList *cl) {
    for (size_t i = 0; i < cl->len; ++i) {
        const Coset *c = &cl->data[i];
        if (c->len <= 1) continue;

        // Tomamos el valor del primer índice del coset como referencia
        int reference_value = x[c->data[0]];
        
        for (size_t j = 1; j < c->len; ++j) {
            // Usamos un pequeño margen de error (epsilon) para comparaciones double
            if ((x[c->data[j]] != reference_value) ) {
                return false; 
            }
        }
    }
    return true;
}

void generar_opciones_comprimidas(int N, int C, const CosetList *cl) {
    if (N % C != 0) {
        printf("Error: N=%d no es divisible por C=%d\n", N, C);
        return;
    }

    int L = N / C;
    FILE *f = fopen(PATH_COMP_COSETS, "wb");
    if (!f) return;

    long long total_posibles = 1;
    for (int i = 0; i < L; i++) total_posibles *= (C + 1);

    int *secuencia = (int *)calloc(L, sizeof(int));
    size_t guardados = 0;
    // Convertimos temporalmente a double para usar la función de verificación
    double complex *temp_v = malloc(L * sizeof(double complex));

    printf("Possible compressed sequences: %zu \n", total_posibles);
    for (long long i = 0; i < total_posibles; i++) {

        for(int j=0; j<L; j++) temp_v[j] = (double)secuencia[j] +0.0*I; // Convertir a double complex para la función de verificación
        uint8_t temp_v_uint8[L];
        double complex temp_v_complex[L];
        /*
        dft(temp_v, temp_v_complex, L);
        printf("Generated dft %lld: ", i);
        for (int j = 0; j < L; j++) {
            printf("%.6f%+.6fi| ", creal(temp_v_complex[j]), cimag(temp_v_complex[j]));
        }
        printf("\n");
        */
        psd(temp_v, temp_v_uint8, L);
        bool is_less_than_half = true;
        /*
        printf("Checking the psd sequence %lld: ", i);
        for (int j = 0; j < L; j++) {
            printf("%d ", temp_v_uint8[j]);
        }
        printf("\n");
        */
        for (int j = 1; j < L && is_less_than_half; j++) {
            is_less_than_half = is_less_than_half && ((temp_v_uint8[j] << 1) < (N + 1));
        }
        if (is_constant_on_cosets(secuencia, cl) && is_less_than_half) {
            fwrite(secuencia, sizeof(int), L, f);
            guardados++;
        }
        /*
        // imprime el vector secuencia.
        printf("Generated: [");
        for (int j = 0; j < L; j++) {
            printf("%d", secuencia[j]);
            if (j + 1 < L) printf(", ");
        }
        printf("] - %s\n", (is_constant_on_cosets(secuencia, cl))? "ACCEPTED" : "REJECTED");
        printf(" %s\n", is_less_than_half ? "VERDADERO" : "FALSO");
        */
        for (int j = L - 1; j >= 0; j--) {
            if (secuencia[j] < C) {
                secuencia[j]++;
                break;
            } else {
                secuencia[j] = 0;
            }
        }
    }
    fflush(f);
    fclose(f);
    free(secuencia);
    free(temp_v);
}

void analizar_cosets_comprimidos(size_t N, size_t k, size_t C) {
    size_t L = N / C; 
    CosetList cl_comprimida = cyclotomic_cosets(k, L);

    FILE *f_in = fopen(PATH_COMP_COSETS, "rb");
    if (!f_in) {
        printf("ERROR: No existe el archivo %s\n", PATH_COMP_COSETS);
        return;
    }

    fseek(f_in, 0, SEEK_END);
    long tam = ftell(f_in);
    rewind(f_in);

    FILE *f_out = fopen(PATH_COMP_DFT_COSETS, "wb");
    int *v_int = malloc(L * sizeof(int));
    double complex *x_in = malloc(L * sizeof(double complex));
    double complex *X_out = malloc(L * sizeof(double complex));
    int *powers_int = malloc(L * sizeof(int));
    size_t count = 0;

    while (fread(v_int, sizeof(int), L, f_in) == L) {
        
        for (size_t i = 0; i < L; i++) {
            x_in[i] = (double)v_int[i] + 0.0 * I;
        }

        dft(x_in, X_out, L);

        for (size_t i = 0; i < L; i++) {
            double power = creal(X_out[i]) * creal(X_out[i]) + 
                           cimag(X_out[i]) * cimag(X_out[i]);
            powers_int[i] = (int)round(power);
        }
        
        fwrite(powers_int, sizeof(int), L, f_out);
        count++;
    }

    printf("Processed %zu secuences for DFT calculations\n", count);

    fclose(f_in); fclose(f_out);
    free(v_int); free(x_in); free(X_out); free(powers_int);
    free_cosetlist(&cl_comprimida);
}


void process_compressed_cosets(const double *comprimido, size_t N, const CosetList *cl) {
    
    // 1. Verificar si es constante en los cosets
    if (!is_constant_on_cosets(comprimido, cl)) {
        return;
    }

    // 2. Preparar datos para la DFT (de double a double complex)
    double complex *x_complex = malloc(N * sizeof(double complex));
    double complex *X_output = malloc(N * sizeof(double complex));
    
    if (!x_complex || !X_output) {
        fprintf(stderr, "Error de memoria en DFT\n");
        free(x_complex); free(X_output);
        return;
    }

    for (size_t i = 0; i < N; i++) {
        x_complex[i] = comprimido[i] + 0.0 * I;
    }

    // 3. Calcular la DFT usando tu función de mymath.c
    dft(x_complex, X_output, N);

    // 4. Guardar abs(DFT)^2 en el fichero
    FILE *f = fopen(PATH_COMP_DFT_COSETS, "a");
    if (f) {
        for (size_t k = 0; k < N; k++) {
            // Magnitud al cuadrado: real^2 + imag^2
            double power = creal(X_output[k]) * creal(X_output[k]) + 
                           cimag(X_output[k]) * cimag(X_output[k]);
            
            fprintf(f, "%f%s", power, (k == N - 1) ? "" : " ");
        }
        fprintf(f, "\n");
        fclose(f);
    }

    // Limpieza
    free(x_complex);
    free(X_output);
}

void buscar_pares_complementarios(size_t N, size_t C) {
    FILE *f_dft = fopen(PATH_COMP_DFT_COSETS, "rb");
    FILE *f_orig = fopen(PATH_COMP_COSETS, "rb");
    
    if (!f_dft || !f_orig) {
        printf("ERROR BUSCAR_PARES_COMPLEMENTARIOS: No se pudo abrir %s o %s\n", PATH_COMP_DFT_COSETS, PATH_COMP_COSETS);
        if (f_dft) fclose(f_dft);
        if (f_orig) fclose(f_orig);
        return;
    }

    double objetivo_double = (double) (N + 1) / 2.0;
    int objetivo = (int)round(objetivo_double);
    size_t L = N / C;

    // 1. Cargamos las magnitudes (Columna 1) y las secuencias originales
    size_t capacidad = 2000;
    size_t total = 0;
    int *col1 = malloc(capacidad * sizeof(int));
    int *secuencias_raw = malloc(capacidad * L * sizeof(int));

    // Buffers temporales para leer cada bloque
    int *temp_dft = malloc(L * sizeof(int));
    int *temp_orig = malloc(L * sizeof(int));

    char buffer[4096];
    // Leemos ambos archivos en paralelo para mantener la correspondencia de líneas
    while (fread(temp_dft, sizeof(int), L, f_dft) == L &&  fread(temp_orig, sizeof(int), L, f_orig) == L) {
        if (total >= capacidad) {
            capacidad *= 2;
            col1 = realloc(col1, capacidad * sizeof(int));
            secuencias_raw = realloc(secuencias_raw, capacidad * L * sizeof(int));
        }
        
        // Guardamos el valor en el índice 1 de la DFT (que es lo que comparamos)
        col1[total] = temp_dft[1]; 
        /*for(int k=0; k<L; k++) printf("%d ", temp_dft[k]);
        printf("\nComparando: %d + %d = %d (Busco %d)\n", col1[total], col1[total], col1[total]+col1[total], objetivo);
        */
        // Copiamos la secuencia original completa al bloque contiguo
        memcpy(&secuencias_raw[total * L], temp_orig, L * sizeof(int));
        
        total++;
    }
    fclose(f_dft);
    fclose(f_orig);
    free(temp_dft);
    free(temp_orig);

    // 2. Comparación y guardado de pares
    FILE *f_out = fopen(PATH_COMP_LP, "wb");
    if (!f_out) {
        printf("ERROR: No se pudo crear %s\n", PATH_COMP_LP);
        return;
    }

    size_t encontrados = 0;
    for (size_t i = 0; i < total; i++) {
        for (size_t j = i + 1; j < total; j++) {
            if ((col1[i] + col1[j]) == objetivo) {
                
                fwrite(&secuencias_raw[i * L], sizeof(int), L, f_out);
                fwrite(&secuencias_raw[j * L], sizeof(int), L, f_out);
                
                encontrados++;
            }
        }
    }

    printf("\nCompressed candidates: %d \n", encontrados);

    
    free(secuencias_raw);
    free(col1);
    fclose(f_out);
}

void generar_sub_combs(int n, int k, uint8_t **res, int *count) {
    for (int i = 0; i < (1 << n); i++) {
        int ones = 0;
        for (int j = 0; j < n; j++) {
            if ((i >> j) & 1) ones++;
        }
        if (ones == k) {
            for (int j = 0; j < n; j++) res[*count][j] = (i >> j) & 1;
            (*count)++;
        }
    }
}

// Coeficiente binomial
int nCr(int n, int r) {
    if (r > n || r < 0) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n / 2) r = n - r;
    long res = 1;
    for (int i = 1; i <= r; ++i) res = res * (n - i + 1) / i;
    return (int)res;
}

// Recursión para el producto cartesiano de todos los bloques
void expandir_recursivo(int bloque, int L, int C, int *pesos, uint8_t ***tablas, int *indices, FILE *f) {
    if (bloque == L) {
        // Escribimos la combinación completa bit a bit
        // IMPORTANTE: Sin fprintf y SIN \n para mantener el desplazamiento N exacto
        for (int m = 0; m < C; m++) {
            for (int b = 0; b < L; b++) {
                // Escribimos el valor numérico 0 o 1 (1 byte por bit)
                fputc(tablas[b][indices[b]][m], f);
            }
        }
        return;
    }
    
    int num = nCr(C, pesos[bloque]);
    for (int i = 0; i < num; i++) {
        indices[bloque] = i;
        expandir_recursivo(bloque + 1, L, C, pesos, tablas, indices, f);
    }
}




// Modificamos ligeramente procesar_linea para que devuelva el número de combinaciones generadas
long procesar_linea(int *pesos, int L, int C, FILE *f_out) {
    uint8_t ***tablas = malloc(L * sizeof(uint8_t **));
    long combinaciones_esta_linea = 1;

    for (int b = 0; b < L; b++) {
        int num = nCr(C, pesos[b]);
        combinaciones_esta_linea *= num; // Multiplicamos las posibilidades de cada bloque
        tablas[b] = malloc(num * sizeof(uint8_t *));
        for (int i = 0; i < num; i++) tablas[b][i] = malloc(C);
        int cnt = 0;
        generar_sub_combs(C, pesos[b], tablas[b], &cnt);
    }

    int *indices = malloc(L * sizeof(int));
    expandir_recursivo(0, L, C, pesos, tablas, indices, f_out);
    
    // Limpieza
    for (int b = 0; b < L; b++) {
        for (int i = 0; i < nCr(C, pesos[b]); i++) free(tablas[b][i]);
        free(tablas[b]);
    }
    free(tablas); 
    free(indices);

    return combinaciones_esta_linea;
}
// Función auxiliar para verificar si ya procesamos estos pesos
bool ya_procesado(int *pesos, int L, int **vistos, int *num_vistos) {
    for (int i = 0; i < *num_vistos; i++) {
        bool coinciden = true;
        for (int j = 0; j < L; j++) {
            if (vistos[i][j] != pesos[j]) {
                coinciden = false;
                break;
            }
        }
        if (coinciden) return true;
    }
    return false;
}

void descomprimir(int N, int C) {
    FILE *f_in = fopen(PATH_COMP_LP, "rb");
    FILE *f_out = fopen(PATH_COMBINATIONS, "wb");
    int L = (N / C);
    
    if (!f_in || !f_out) {
        printf("ERROR DESCOMPRIMIR: No se pudieron abrir %s o %s\n", PATH_COMP_LP, PATH_COMBINATIONS);
        if (f_in) fclose(f_in);
        if (f_out) fclose(f_out);
        return;
    }

    int *pesos = malloc(L * sizeof(int));
    long total_descomprimidos = 0;
    
    // Registro de pesos vistos (ajusta la capacidad según necesites)
    int capacidad_vistos = 1000;
    int num_vistos = 0;
    int **vistos = malloc(capacidad_vistos * sizeof(int *));

    // El archivo comp_lp tiene pares, pero fscanf saltará los espacios
    // leeremos número a número.
    while (fread(pesos, sizeof(int), L, f_in) == (size_t)L) {
        if (!ya_procesado(pesos, L, vistos, &num_vistos)) {
            
            // 1. Guardar en el registro de vistos
            if (num_vistos >= capacidad_vistos) {
                capacidad_vistos *= 2;
                vistos = realloc(vistos, capacidad_vistos * sizeof(int *));
            }
            vistos[num_vistos] = malloc(L * sizeof(int));
            memcpy(vistos[num_vistos], pesos, L * sizeof(int));
            num_vistos++;

            // 2. Descomprimir solo si es nuevo
            total_descomprimidos += procesar_linea(pesos, L, C, f_out);
        }
    }

    // Limpieza de memoria local
    for (int i = 0; i < num_vistos; i++) free(vistos[i]);
    free(vistos);
    free(pesos);
    
    fclose(f_in);
    fclose(f_out);

    printf("Decompression complete. Only processed: %d. Total vectors in  %s: %ld\n", num_vistos, PATH_COMBINATIONS, total_descomprimidos);
}


/* Versión optimizada para memoria: Lee de archivo en lugar de estructura */
int save_dft_from_file(size_t N) {
    
    double constant = (double)(N + 1) / 2.0;
    FILE *f_in = fopen(PATH_COMBINATIONS, "rb");
    FILE *file1 = fopen(PATH_DFT, "wb");
    FILE *file2 = fopen(PATH_CTEDFT, "wb");
    
    if (!f_in || !file1 || !file2) {
        fprintf(stderr, "ERROR SAVE_DFT_FROM_FILE: no se pudo abrir los archivos\n");
        if (f_in) fclose(f_in);
        if (file1) fclose(file1);
        if (file2) fclose(file2);
        return 0;
    }
    
    double complex *time_domain = (double complex*)malloc(N * sizeof(double complex));
    double complex *freq_domain = (double complex*)malloc(N * sizeof(double complex));
    // Buffer para leer la secuencia binaria (N bytes)
    unsigned char *seq_buffer = (unsigned char*)malloc(N);
    // Buffers para los resultados (N-1 enteros)
    int *buffer_dft = (int*)malloc((N - 1) * sizeof(int));
    int *buffer_ctedft = (int*)malloc((N - 1) * sizeof(int));

    if (!time_domain || !freq_domain || !seq_buffer || !buffer_dft || !buffer_ctedft) {
        fprintf(stderr, "ERROR: fallo de memoria en save_dft_from_file\n");
        return 0;
    }
    
    // Leemos bloques de N bytes directamente
    while (fread(seq_buffer, sizeof(unsigned char), N, f_in) == N) {
        
        // Convertir los bytes leídos a complejo
        for (size_t i = 0; i < N; i++) {
            // Asumiendo que se guardaron como caracteres '0' o '1'
            time_domain[i] = (seq_buffer[i] == 0) ? 1.0 + 0.0*I : 0.0 + 0.0*I;
        }
        // Calcular DFT
        dft(time_domain, freq_domain, N);
        
        // Escribir en primer fich: módulo redondeado
        for (size_t j = 1; j < N; j++) {
            double modulo = pow(cabs(freq_domain[j]), 2);
            double rounded = rint(modulo);
            double transformed = constant - rounded;

            buffer_dft[j-1] = (int)rounded;
            buffer_ctedft[j-1] = (int)rint(transformed);
        }
        
        fwrite(buffer_dft, sizeof(int), N - 1, file1);
        fwrite(buffer_ctedft, sizeof(int), N - 1, file2);
    }
    
    free(time_domain);
    free(freq_domain);
    free(seq_buffer);
    free(buffer_dft);
    free(buffer_ctedft);
    fclose(f_in);
    fclose(file1);
    fclose(file2);

    printf("DFT and cte-DFT processing completed\n");
    return 1;
}

int main(void) {

    clock_t start_time = clock();

    int N = 25;//15
    int k = 1;//2
    int C = 5;//3
    int L =N/C;
    printf("N=%d k=%d C=%d L=%d\n\n", N, k, C, L);

    // 1. Generar los cosets para el espacio comprimido
    CosetList cl_comprimida = cyclotomic_cosets(k, L);

    // 2. Generar todas las combinaciones de pesos (0 a C) que cumplen la simetría de cosets
    generar_opciones_comprimidas(N, C, &cl_comprimida);

    // 3. Calcular la DFT de esas opciones comprimidas y guardarlas
    analizar_cosets_comprimidos(N, k, C);

    // 4. Buscar qué pares de magnitudes suman el objetivo y guardarlos en comp_lp.txt
    buscar_pares_complementarios(N, C);

    // 5. Expandir esos pares de pesos a sus combinaciones binarias finales (0s y 1s)
    descomprimir(N, C);

    // 6. Guardar DFTS de las combinaciones y cte-DFT
    save_dft_from_file(N);
    
    // 7. Buscar pares de las combinaciones y cte-DFT
    find_matches_files(N);


    //PARA MEJORAR HABRIA QUE HACER QUE LOS PARES COMPLEMENTARIOS NO SE REPITIESEN

    // Limpieza final
    free_cosetlist(&cl_comprimida);


    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("TOTAL TIME: %.3f seconds\n", cpu_time_used);
    printf("\nThat's all, folks!");
    return 0;
}