#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

// Función auxiliar para escribir la DFT desenrollada automáticamente
string get_unrolled_dft(int s, int v, const string& array_name) {
    ostringstream oss;
    oss << "std::complex<double>(" << array_name << "[0])*omegaPowers[0]";
    for (int k = 1; k < v; ++k) {
        int exp = (s * k) % v;
        oss << " + std::complex<double>(" << array_name << "[" << k << "])*omegaPowers[" << exp << "]";
    }
    return oss.str();
}

void generarCodigo(std::ofstream& out, int longitud) {
    int i, v = longitud; 
    
    // Cálculo de raíces
    double pi = 3.14159265358979323846;
    double angle = -2.0 * pi / (double)longitud;
    std::complex<double> omega(cos(angle), sin(angle)); 
    
    // --- CABECERAS Y FUNCION HASH ---
    out << "#include <stdio.h>\n"
        << "#include <stdlib.h>\n"
        << "#include <complex>\n"
        << "#include <math.h>\n"
        << "#include <errno.h>\n"
        << "#include <string.h>\n"
        << "#include <vector>\n"
        << "#include <atomic>\n"
        << "#include <cstdint>\n"
        << "#include <omp.h>\n\n"
        << "// Compilar con: g++ -fopenmp -O3 -march=native -o compress." << longitud << " compress." << longitud << ".cpp\n\n";

    // Función Hash FNV-1a integrada en el código generado
    out << "inline uint64_t hash_array(int* arr, int len, uint64_t seed) {\n"
        << "    uint64_t hash = 14695981039346656037ULL + seed;\n"
        << "    for(int i = 0; i < len; i++) {\n"
        << "        hash ^= (uint64_t)arr[i];\n"
        << "        hash *= 1099511628211ULL;\n"
        << "    }\n"
        << "    return hash;\n"
        << "}\n\n";

    out << "int main(int argc, char *argv[])\n{\n";

    // --- GESTIÓN DE LOTES (Aplicado a B) ---
    out << "    int lote_id = 0;\n"
        << "    int inicio_bucle_b = 1;\n"
        << "    int fin_bucle_b = 10;\n\n"
        << "    if (argc >= 2) {\n"
        << "        lote_id = atoi(argv[1]);\n"
        << "        if (lote_id >= 1 && lote_id <= 10) {\n"
        << "            printf(\"Modo Lote: Ejecutando Fase 2 (B) para lote %d\\n\", lote_id);\n"
        << "            inicio_bucle_b = lote_id;\n"
        << "            fin_bucle_b = lote_id;\n"
        << "        }\n"
        << "    }\n\n";

    out << "    omp_set_num_threads(omp_get_max_threads());\n"
        << "    printf(\"Usando %d hilos OpenMP\\n\", omp_get_max_threads());\n\n";

    // --- INICIALIZACIÓN DEL BLOOM FILTER ---
    // Usamos 1000 millones de bits (~120 MB de RAM), seguro para Altamira
    out << "    uint64_t bf_size = 1000000007ULL;\n"
        << "    std::vector<std::atomic<uint8_t>> bloom_filter(bf_size);\n"
        << "    for(size_t k=0; k<bf_size; k++) bloom_filter[k].store(0, std::memory_order_relaxed);\n\n";

    out << "    int v = " << longitud << ";\n"
        << "    int objetivo = (2*v)+2;\n\n";

    out << "    char filename[64];\n"
        << "    if (lote_id > 0) sprintf(filename, \"Candidatos_B_lote_%d.txt\", lote_id);\n"
        << "    else sprintf(filename, \"Candidatos_B_completos.txt\");\n"
        << "    FILE *f = fopen(filename, \"w\");\n\n";

    out << "    std::complex<double> omegaPowers[" << longitud << "];\n";
    for(int k = 0; k < longitud; k++) {
        std::complex<double> p = pow(omega, k);
        out << "    omegaPowers[" << k << "] = std::complex<double>(" << p.real() << ", " << p.imag() << ");\n";
    }

    out << "\n"
        << "    int aOp1s[] = {0,3,7,9,13}; \n"
        << "    int aOm1s[] = {0,6,8,10}; \n"
        << "    int aOp3s[] = {0,2,4,11,14}; \n"
        << "    int aOm3s[] = {0,1,5,12,15}; \n"
        << "    int axm1[5], axp1[5], axm3[5], axp3[5], A[76]; \n"
        << "    int bOp1s[] = {0,6,7,12,14,15}; \n"
        << "    int bOm1s[] = {0,1,3,4,13}; \n"
        << "    int bOp3s[] = {0,8,9,10}; \n"
        << "    int bOm3s[] = {0,2,5,11}; \n"
        << "    int bxm1[5], bxp1[6], bxm3[4], bxp3[4], B[76]; \n\n";

    // ========================================================================
    // FASE 1: RECORRER SECUENCIAS A Y LLENAR BLOOM FILTER
    // ========================================================================
    out << "    printf(\"\\n--- FASE 1: Llenando Bloom Filter con combinaciones de A ---\\n\");\n";
    out << "    #pragma omp parallel for collapse(4)\n";
    out << "    for (int axp1_1 = 1; axp1_1 <= 10; ++axp1_1) {\n";
    out << "    for (int axp1_2 = 1; axp1_2 <= 10; ++axp1_2) {\n";
    out << "    for (int axp1_3 = 1; axp1_3 <= 10; ++axp1_3) {\n";
    out << "    for (int axp1_4 = 1; axp1_4 <= 10; ++axp1_4) {\n";
    out << "        axp1[1] = axp1_1; axp1[2] = axp1_2; axp1[3] = axp1_3; axp1[4] = axp1_4;\n";
    out << "        for(int i=1; i<=4; i++){\n"
        << "            switch (axp1[i]) {\n"
        << "                case 1 : A[aOp1s[i]] = -1 ; A[aOp1s[i]+15] = -1 ; A[aOp1s[i]+30] = 1 ; A[aOp1s[i]+45] = 1 ; A[aOp1s[i]+60] = 1 ; break;\n"
        << "                case 2 : A[aOp1s[i]] = -1 ; A[aOp1s[i]+15] = 1 ; A[aOp1s[i]+30] = -1 ; A[aOp1s[i]+45] = 1 ; A[aOp1s[i]+60] = 1 ; break;\n"
        << "                case 3 : A[aOp1s[i]] = -1 ; A[aOp1s[i]+15] = 1 ; A[aOp1s[i]+30] = 1 ; A[aOp1s[i]+45] = -1 ; A[aOp1s[i]+60] = 1 ; break;\n"
        << "                case 4 : A[aOp1s[i]] = -1 ; A[aOp1s[i]+15] = 1 ; A[aOp1s[i]+30] = 1 ; A[aOp1s[i]+45] = 1 ; A[aOp1s[i]+60] = -1 ; break;\n"
        << "                case 5 : A[aOp1s[i]] = 1 ; A[aOp1s[i]+15] = -1 ; A[aOp1s[i]+30] = -1 ; A[aOp1s[i]+45] = 1 ; A[aOp1s[i]+60] = 1 ; break;\n"
        << "                case 6 : A[aOp1s[i]] = 1 ; A[aOp1s[i]+15] = -1 ; A[aOp1s[i]+30] = 1 ; A[aOp1s[i]+45] = -1 ; A[aOp1s[i]+60] = 1 ; break;\n"
        << "                case 7 : A[aOp1s[i]] = 1 ; A[aOp1s[i]+15] = -1 ; A[aOp1s[i]+30] = 1 ; A[aOp1s[i]+45] = 1 ; A[aOp1s[i]+60] = -1 ; break;\n"
        << "                case 8 : A[aOp1s[i]] = 1 ; A[aOp1s[i]+15] = 1 ; A[aOp1s[i]+30] = -1 ; A[aOp1s[i]+45] = -1 ; A[aOp1s[i]+60] = 1 ; break;\n"
        << "                case 9 : A[aOp1s[i]] = 1 ; A[aOp1s[i]+15] = 1 ; A[aOp1s[i]+30] = -1 ; A[aOp1s[i]+45] = 1 ; A[aOp1s[i]+60] = -1 ; break;\n"
        << "                case 10 : A[aOp1s[i]] = 1 ; A[aOp1s[i]+15] = 1 ; A[aOp1s[i]+30] = 1 ; A[aOp1s[i]+45] = -1 ; A[aOp1s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    out << "    for (int axm1_1 = 1; axm1_1 <= 10; ++axm1_1) {\n"
        << "    for (int axm1_2 = 1; axm1_2 <= 10; ++axm1_2) {\n"
        << "    for (int axm1_3 = 1; axm1_3 <= 10; ++axm1_3) {\n"
        << "        axm1[1] = axm1_1; axm1[2] = axm1_2; axm1[3] = axm1_3;\n"
        << "        for(int i=1; i<=3; i++){\n"
        << "            switch (axm1[i]) {\n"
        << "                case 1 : A[aOm1s[i]] = -1 ; A[aOm1s[i]+15] = -1 ; A[aOm1s[i]+30] = -1 ; A[aOm1s[i]+45] = 1 ; A[aOm1s[i]+60] = 1 ; break;\n"
        << "                case 2 : A[aOm1s[i]] = -1 ; A[aOm1s[i]+15] = -1 ; A[aOm1s[i]+30] = 1 ; A[aOm1s[i]+45] = -1 ; A[aOm1s[i]+60] = 1 ; break;\n"
        << "                case 3 : A[aOm1s[i]] = -1 ; A[aOm1s[i]+15] = -1 ; A[aOm1s[i]+30] = 1 ; A[aOm1s[i]+45] = 1 ; A[aOm1s[i]+60] = -1 ; break;\n"
        << "                case 4 : A[aOm1s[i]] = -1 ; A[aOm1s[i]+15] = 1 ; A[aOm1s[i]+30] = -1 ; A[aOm1s[i]+45] = -1 ; A[aOm1s[i]+60] = 1 ; break;\n"
        << "                case 5 : A[aOm1s[i]] = -1 ; A[aOm1s[i]+15] = 1 ; A[aOm1s[i]+30] = -1 ; A[aOm1s[i]+45] = 1 ; A[aOm1s[i]+60] = -1 ; break;\n"
        << "                case 6 : A[aOm1s[i]] = -1 ; A[aOm1s[i]+15] = 1 ; A[aOm1s[i]+30] = 1 ; A[aOm1s[i]+45] = -1 ; A[aOm1s[i]+60] = -1 ; break;\n"
        << "                case 7 : A[aOm1s[i]] = 1 ; A[aOm1s[i]+15] = -1 ; A[aOm1s[i]+30] = -1 ; A[aOm1s[i]+45] = -1 ; A[aOm1s[i]+60] = 1 ; break;\n"
        << "                case 8 : A[aOm1s[i]] = 1 ; A[aOm1s[i]+15] = -1 ; A[aOm1s[i]+30] = -1 ; A[aOm1s[i]+45] = 1 ; A[aOm1s[i]+60] = -1 ; break;\n"
        << "                case 9 : A[aOm1s[i]] = 1 ; A[aOm1s[i]+15] = -1 ; A[aOm1s[i]+30] = 1 ; A[aOm1s[i]+45] = -1 ; A[aOm1s[i]+60] = -1 ; break;\n"
        << "                case 10 : A[aOm1s[i]] = 1 ; A[aOm1s[i]+15] = 1 ; A[aOm1s[i]+30] = -1 ; A[aOm1s[i]+45] = -1 ; A[aOm1s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    out << "    for (int axp3_1 = 1; axp3_1 <= 5; ++axp3_1) {\n"
        << "    for (int axp3_2 = 1; axp3_2 <= 5; ++axp3_2) {\n"
        << "    for (int axp3_3 = 1; axp3_3 <= 5; ++axp3_3) {\n"
        << "    for (int axp3_4 = 1; axp3_4 <= 5; ++axp3_4) {\n"
        << "        axp3[1] = axp3_1; axp3[2] = axp3_2; axp3[3] = axp3_3; axp3[4] = axp3_4;\n"
        << "        for(int i=1; i<=4; i++){\n"
        << "            switch (axp3[i]) {\n"
        << "                case 1 : A[aOp3s[i]] = -1 ; A[aOp3s[i]+15] = 1 ; A[aOp3s[i]+30] = 1 ; A[aOp3s[i]+45] = 1 ; A[aOp3s[i]+60] = 1 ; break;\n"
        << "                case 2 : A[aOp3s[i]] = 1 ; A[aOp3s[i]+15] = -1 ; A[aOp3s[i]+30] = 1 ; A[aOp3s[i]+45] = 1 ; A[aOp3s[i]+60] = 1 ; break;\n"
        << "                case 3 : A[aOp3s[i]] = 1 ; A[aOp3s[i]+15] = 1 ; A[aOp3s[i]+30] = -1 ; A[aOp3s[i]+45] = 1 ; A[aOp3s[i]+60] = 1 ; break;\n"
        << "                case 4 : A[aOp3s[i]] = 1 ; A[aOp3s[i]+15] = 1 ; A[aOp3s[i]+30] = 1 ; A[aOp3s[i]+45] = -1 ; A[aOp3s[i]+60] = 1 ; break;\n"
        << "                case 5 : A[aOp3s[i]] = 1 ; A[aOp3s[i]+15] = 1 ; A[aOp3s[i]+30] = 1 ; A[aOp3s[i]+45] = 1 ; A[aOp3s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    out << "    for (int axm3_1 = 1; axm3_1 <= 5; ++axm3_1) {\n"
        << "    for (int axm3_2 = 1; axm3_2 <= 5; ++axm3_2) {\n"
        << "    for (int axm3_3 = 1; axm3_3 <= 5; ++axm3_3) {\n"
        << "    for (int axm3_4 = 1; axm3_4 <= 5; ++axm3_4) {\n"
        << "        axm3[1] = axm3_1; axm3[2] = axm3_2; axm3[3] = axm3_3; axm3[4] = axm3_4;\n"
        << "        for(int i=1; i<=4; i++){\n"
        << "            switch (axm3[i]) {\n"
        << "                case 1 : A[aOm3s[i]] = -1 ; A[aOm3s[i]+15] = -1 ; A[aOm3s[i]+30] = -1 ; A[aOm3s[i]+45] = -1 ; A[aOm3s[i]+60] = 1 ; break;\n"
        << "                case 2 : A[aOm3s[i]] = -1 ; A[aOm3s[i]+15] = -1 ; A[aOm3s[i]+30] = -1 ; A[aOm3s[i]+45] = 1 ; A[aOm3s[i]+60] = -1 ; break;\n"
        << "                case 3 : A[aOm3s[i]] = -1 ; A[aOm3s[i]+15] = -1 ; A[aOm3s[i]+30] = 1 ; A[aOm3s[i]+45] = -1 ; A[aOm3s[i]+60] = -1 ; break;\n"
        << "                case 4 : A[aOm3s[i]] = -1 ; A[aOm3s[i]+15] = 1 ; A[aOm3s[i]+30] = -1 ; A[aOm3s[i]+45] = -1 ; A[aOm3s[i]+60] = -1 ; break;\n"
        << "                case 5 : A[aOm3s[i]] = 1 ; A[aOm3s[i]+15] = -1 ; A[aOm3s[i]+30] = -1 ; A[aOm3s[i]+45] = -1 ; A[aOm3s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    // Cálculo PSD para A y llenado del Bloom Filter
    out << "        int target_B[" << v << "];\n"
        << "        bool valid = true;\n"
        << "        for (int s = 1; s < v; s++) {\n"
        << "            std::complex<double> valA;\n";
    
    // Switch de frecuencias para desenrollar
    out << "            switch(s) {\n";
    for(int s=1; s<v; s++) {
        out << "                case " << s << ": valA = " << get_unrolled_dft(s, v, "A") << "; break;\n";
    }
    out << "            }\n";

    out << "            int psd_a = (int)round(std::real(valA)*std::real(valA) + std::imag(valA)*std::imag(valA));\n"
        << "            target_B[s] = objetivo - psd_a;\n"
        << "            if (target_B[s] < 0) { valid = false; break; }\n"
        << "        }\n"
        << "        if (valid) {\n"
        << "            uint64_t h1 = hash_array(&target_B[1], v-1, 1) % bf_size;\n"
        << "            uint64_t h2 = hash_array(&target_B[1], v-1, 2) % bf_size;\n"
        << "            uint64_t h3 = hash_array(&target_B[1], v-1, 3) % bf_size;\n"
        << "            bloom_filter[h1].store(1, std::memory_order_relaxed);\n"
        << "            bloom_filter[h2].store(1, std::memory_order_relaxed);\n"
        << "            bloom_filter[h3].store(1, std::memory_order_relaxed);\n"
        << "        }\n";

    // Cierres de los 15 bucles de A
    for(int j=0; j<15; j++) out << "    }\n";
    out << "    printf(\"Fase 1 completada. Bloom Filter poblado.\\n\");\n\n";

    // ========================================================================
    // FASE 2: RECORRER SECUENCIAS B (AQUÍ APLICAN LOS LOTES)
    // ========================================================================
    out << "    printf(\"\\n--- FASE 2: Buscando match para B (Lote %d) ---\\n\", lote_id);\n";
    out << "    #pragma omp parallel for collapse(4)\n";
    
    // APLICACIÓN DEL LOTE AL PRIMER BUCLE DE B
    out << "    for (int bxp1_1 = inicio_bucle_b; bxp1_1 <= fin_bucle_b; ++bxp1_1) {\n";
    out << "    for (int bxp1_2 = 1; bxp1_2 <= 10; ++bxp1_2) {\n";
    out << "    for (int bxp1_3 = 1; bxp1_3 <= 10; ++bxp1_3) {\n";
    out << "    for (int bxp1_4 = 1; bxp1_4 <= 10; ++bxp1_4) {\n";
    out << "    for (int bxp1_5 = 1; bxp1_5 <= 10; ++bxp1_5) {\n";
    out << "        bxp1[1] = bxp1_1; bxp1[2] = bxp1_2; bxp1[3] = bxp1_3; bxp1[4] = bxp1_4; bxp1[5] = bxp1_5;\n";
    out << "        for(int i=1; i<=5; i++){\n"
        << "            switch (bxp1[i]) {\n"
        << "                case 1 : B[bOp1s[i]] = -1 ; B[bOp1s[i]+15] = -1 ; B[bOp1s[i]+30] = 1 ; B[bOp1s[i]+45] = 1 ; B[bOp1s[i]+60] = 1 ; break;\n"
        << "                case 2 : B[bOp1s[i]] = -1 ; B[bOp1s[i]+15] = 1 ; B[bOp1s[i]+30] = -1 ; B[bOp1s[i]+45] = 1 ; B[bOp1s[i]+60] = 1 ; break;\n"
        << "                case 3 : B[bOp1s[i]] = -1 ; B[bOp1s[i]+15] = 1 ; B[bOp1s[i]+30] = 1 ; B[bOp1s[i]+45] = -1 ; B[bOp1s[i]+60] = 1 ; break;\n"
        << "                case 4 : B[bOp1s[i]] = -1 ; B[bOp1s[i]+15] = 1 ; B[bOp1s[i]+30] = 1 ; B[bOp1s[i]+45] = 1 ; B[bOp1s[i]+60] = -1 ; break;\n"
        << "                case 5 : B[bOp1s[i]] = 1 ; B[bOp1s[i]+15] = -1 ; B[bOp1s[i]+30] = -1 ; B[bOp1s[i]+45] = 1 ; B[bOp1s[i]+60] = 1 ; break;\n"
        << "                case 6 : B[bOp1s[i]] = 1 ; B[bOp1s[i]+15] = -1 ; B[bOp1s[i]+30] = 1 ; B[bOp1s[i]+45] = -1 ; B[bOp1s[i]+60] = 1 ; break;\n"
        << "                case 7 : B[bOp1s[i]] = 1 ; B[bOp1s[i]+15] = -1 ; B[bOp1s[i]+30] = 1 ; B[bOp1s[i]+45] = 1 ; B[bOp1s[i]+60] = -1 ; break;\n"
        << "                case 8 : B[bOp1s[i]] = 1 ; B[bOp1s[i]+15] = 1 ; B[bOp1s[i]+30] = -1 ; B[bOp1s[i]+45] = -1 ; B[bOp1s[i]+60] = 1 ; break;\n"
        << "                case 9 : B[bOp1s[i]] = 1 ; B[bOp1s[i]+15] = 1 ; B[bOp1s[i]+30] = -1 ; B[bOp1s[i]+45] = 1 ; B[bOp1s[i]+60] = -1 ; break;\n"
        << "                case 10: B[bOp1s[i]] = 1 ; B[bOp1s[i]+15] = 1 ; B[bOp1s[i]+30] = 1 ; B[bOp1s[i]+45] = -1 ; B[bOp1s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    out << "    for (int bxm1_1 = 1; bxm1_1 <= 10; ++bxm1_1) {\n"
        << "    for (int bxm1_2 = 1; bxm1_2 <= 10; ++bxm1_2) {\n"
        << "    for (int bxm1_3 = 1; bxm1_3 <= 10; ++bxm1_3) {\n"
        << "    for (int bxm1_4 = 1; bxm1_4 <= 10; ++bxm1_4) {\n"
        << "        bxm1[1] = bxm1_1; bxm1[2] = bxm1_2; bxm1[3] = bxm1_3; bxm1[4] = bxm1_4;\n"
        << "        for(int i=1; i<=4; i++){\n"
        << "            switch (bxm1[i]) {\n"
        << "                case 1 : B[bOm1s[i]] = -1 ; B[bOm1s[i]+15] = -1 ; B[bOm1s[i]+30] = -1 ; B[bOm1s[i]+45] = 1 ; B[bOm1s[i]+60] = 1 ; break;\n"
        << "                case 2 : B[bOm1s[i]] = -1 ; B[bOm1s[i]+15] = -1 ; B[bOm1s[i]+30] = 1 ; B[bOm1s[i]+45] = -1 ; B[bOm1s[i]+60] = 1 ; break;\n"
        << "                case 3 : B[bOm1s[i]] = -1 ; B[bOm1s[i]+15] = -1 ; B[bOm1s[i]+30] = 1 ; B[bOm1s[i]+45] = 1 ; B[bOm1s[i]+60] = -1 ; break;\n"
        << "                case 4 : B[bOm1s[i]] = -1 ; B[bOm1s[i]+15] = 1 ; B[bOm1s[i]+30] = -1 ; B[bOm1s[i]+45] = -1 ; B[bOm1s[i]+60] = 1 ; break;\n"
        << "                case 5 : B[bOm1s[i]] = -1 ; B[bOm1s[i]+15] = 1 ; B[bOm1s[i]+30] = -1 ; B[bOm1s[i]+45] = 1 ; B[bOm1s[i]+60] = -1 ; break;\n"
        << "                case 6 : B[bOm1s[i]] = -1 ; B[bOm1s[i]+15] = 1 ; B[bOm1s[i]+30] = 1 ; B[bOm1s[i]+45] = -1 ; B[bOm1s[i]+60] = -1 ; break;\n"
        << "                case 7 : B[bOm1s[i]] = 1 ; B[bOm1s[i]+15] = -1 ; B[bOm1s[i]+30] = -1 ; B[bOm1s[i]+45] = -1 ; B[bOm1s[i]+60] = 1 ; break;\n"
        << "                case 8 : B[bOm1s[i]] = 1 ; B[bOm1s[i]+15] = -1 ; B[bOm1s[i]+30] = -1 ; B[bOm1s[i]+45] = 1 ; B[bOm1s[i]+60] = -1 ; break;\n"
        << "                case 9 : B[bOm1s[i]] = 1 ; B[bOm1s[i]+15] = -1 ; B[bOm1s[i]+30] = 1 ; B[bOm1s[i]+45] = -1 ; B[bOm1s[i]+60] = -1 ; break;\n"
        << "                case 10: B[bOm1s[i]] = 1 ; B[bOm1s[i]+15] = 1 ; B[bOm1s[i]+30] = -1 ; B[bOm1s[i]+45] = -1 ; B[bOm1s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    out << "    for (int bxp3_1 = 1; bxp3_1 <= 5; bxp3_1++){\n"
        << "    for (int bxp3_2 = 1; bxp3_2 <= 5; bxp3_2++){\n"
        << "    for (int bxp3_3 = 1; bxp3_3 <= 5; bxp3_3++){\n"
        << "        bxp3[1] = bxp3_1; bxp3[2] = bxp3_2; bxp3[3] = bxp3_3;\n"
        << "        for(int i=1; i<=3; i++){\n"
        << "            switch (bxp3[i]) {\n"
        << "                case 1 : B[bOp3s[i]] = -1 ; B[bOp3s[i]+15] = 1 ; B[bOp3s[i]+30] = 1 ; B[bOp3s[i]+45] = 1 ; B[bOp3s[i]+60] = 1 ; break;\n"
        << "                case 2 : B[bOp3s[i]] = 1 ; B[bOp3s[i]+15] = -1 ; B[bOp3s[i]+30] = 1 ; B[bOp3s[i]+45] = 1 ; B[bOp3s[i]+60] = 1 ; break;\n"
        << "                case 3 : B[bOp3s[i]] = 1 ; B[bOp3s[i]+15] = 1 ; B[bOp3s[i]+30] = -1 ; B[bOp3s[i]+45] = 1 ; B[bOp3s[i]+60] = 1 ; break;\n"
        << "                case 4 : B[bOp3s[i]] = 1 ; B[bOp3s[i]+15] = 1 ; B[bOp3s[i]+30] = 1 ; B[bOp3s[i]+45] = -1 ; B[bOp3s[i]+60] = 1 ; break;\n"
        << "                case 5 : B[bOp3s[i]] = 1 ; B[bOp3s[i]+15] = 1 ; B[bOp3s[i]+30] = 1 ; B[bOp3s[i]+45] = 1 ; B[bOp3s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    out << "    for (int bxm3_1 = 1; bxm3_1 <= 5; bxm3_1++){\n"
        << "    for (int bxm3_2 = 1; bxm3_2 <= 5; bxm3_2++){\n"
        << "    for (int bxm3_3 = 1; bxm3_3 <= 5; bxm3_3++){\n"
        << "        bxm3[1] = bxm3_1; bxm3[2] = bxm3_2; bxm3[3] = bxm3_3;\n"
        << "        for(int i=1; i<=3; i++){\n"
        << "            switch (bxm3[i]) {\n"
        << "                case 1 : B[bOm3s[i]] = -1 ; B[bOm3s[i]+15] = -1 ; B[bOm3s[i]+30] = -1 ; B[bOm3s[i]+45] = -1 ; B[bOm3s[i]+60] = 1 ; break;\n"
        << "                case 2 : B[bOm3s[i]] = -1 ; B[bOm3s[i]+15] = -1 ; B[bOm3s[i]+30] = -1 ; B[bOm3s[i]+45] = 1 ; B[bOm3s[i]+60] = -1 ; break;\n"
        << "                case 3 : B[bOm3s[i]] = -1 ; B[bOm3s[i]+15] = -1 ; B[bOm3s[i]+30] = 1 ; B[bOm3s[i]+45] = -1 ; B[bOm3s[i]+60] = -1 ; break;\n"
        << "                case 4 : B[bOm3s[i]] = -1 ; B[bOm3s[i]+15] = 1 ; B[bOm3s[i]+30] = -1 ; B[bOm3s[i]+45] = -1 ; B[bOm3s[i]+60] = -1 ; break;\n"
        << "                case 5 : B[bOm3s[i]] = 1 ; B[bOm3s[i]+15] = -1 ; B[bOm3s[i]+30] = -1 ; B[bOm3s[i]+45] = -1 ; B[bOm3s[i]+60] = -1 ; break;\n"
        << "            }\n"
        << "        }\n";

    // Cálculo de PSD de B y verificación contra el Bloom Filter
    out << "        int psd_b[" << v << "];\n"
        << "        for (int s = 1; s < v; s++) {\n"
        << "            std::complex<double> valB;\n";
        
    out << "            switch(s) {\n";
    for(int s=1; s<v; s++) {
        out << "                case " << s << ": valB = " << get_unrolled_dft(s, v, "B") << "; break;\n";
    }
    out << "            }\n";

    out << "            psd_b[s] = (int)round(std::real(valB)*std::real(valB) + std::imag(valB)*std::imag(valB));\n"
        << "        }\n\n"
        << "        uint64_t h1 = hash_array(&psd_b[1], v-1, 1) % bf_size;\n"
        << "        uint64_t h2 = hash_array(&psd_b[1], v-1, 2) % bf_size;\n"
        << "        uint64_t h3 = hash_array(&psd_b[1], v-1, 3) % bf_size;\n\n"
        << "        if (bloom_filter[h1].load(std::memory_order_relaxed) &&\n"
        << "            bloom_filter[h2].load(std::memory_order_relaxed) &&\n"
        << "            bloom_filter[h3].load(std::memory_order_relaxed)) {\n"
        << "            \n"
        << "            #pragma omp critical\n"
        << "            {\n"
        << "                fprintf(f, \"CANDIDATO B ENCONTRADO: \");\n"
        << "                for(int x=0; x<v; x++) fprintf(f, \"%d \", B[x]);\n"
        << "                fprintf(f, \"\\n\");\n"
        << "                fflush(f);\n"
        << "            }\n"
        << "        }\n";

    // Cierres de los 15 bucles de B
    for(int j=0; j<15; j++) out << "    }\n";

    out << "    fclose(f);\n";
    out << "    printf(\"Lote %d finalizado exitosamente.\\n\", lote_id);\n";
    out << "    return 0;\n";
    out << "}\n";
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <longitud>\n";
        return 1;
    }
    int longitud = std::stoi(argv[1]);
    
    std::string filename = "compress." + std::to_string(longitud) + ".cpp";
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "No se pudo crear " << filename << "\n";
        return 1;
    }

    generarCodigo(out, longitud);

    std::cout << "Generador ejecutado. Codigo de busqueda guardado en: " << filename << "\n";
    return 0;
}