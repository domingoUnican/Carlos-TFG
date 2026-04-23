#!/bin/bash

# Asegurarnos de estar en la carpeta correcta
if [ ! -d "Tests" ] || [ ! -f "build/app" ]; then
    echo "Error: Ejecuta este script desde la carpeta Code/ o Codigo/ donde está build/app"
    exit 1
fi

# Definir los parámetros para cada longitud (Matemáticas: N = p * q^2)
# Formato del array: "N file1 file2 p q"
CONFIGS=(
    "27 9.txt 9.txt 3 3"
    "45 15.txt 9.txt 5 3"
    "63 21.txt 9.txt 7 3"
    "75 15.txt 25.txt 3 5"
)

for CONFIG in "${CONFIGS[@]}"; do
    read -r N f1 f2 p q <<< "$CONFIG"
    
    echo -e "\n========================================"
    echo "Procesando longitud N=$N (p=$p, q=$q)"
    echo "========================================"

    # Ir a la carpeta del test
    cd Tests/$N || { echo "⚠️ No existe la carpeta Tests/$N. Omitiendo..."; continue; }

    # Ejecutar la aplicación de búsqueda (silenciamos la salida para no manchar la terminal)
    echo "[1/3] Ejecutando el algoritmo DFS..."
    ../../build/app $f1 $f2 $p $q > run_log.txt

    # Archivos que genera tu programa C
    FILE_DFT="${f1}_${f2}_dft.txt"
    FILE_CTE="${f1}_${f2}_cte-dft.txt"
    FILE_COMB="${f1}_${f2}_combinations.txt"

    # Archivo de salida final
    OUT_FILE="../../Resultados_Pares_$N.txt"
    > "$OUT_FILE" # Limpiar/crear archivo

    if [ ! -f "$FILE_DFT" ] || [ ! -f "$FILE_CTE" ]; then
        echo "⚠️ No se generaron los archivos de espectro para N=$N."
        cd ../..
        continue
    fi

    echo "[2/3] Cruzando espectros (Meet-in-the-Middle)..."
    sort "$FILE_DFT" > temp_dft.txt
    sort "$FILE_CTE" > temp_cte.txt
    comm -12 temp_dft.txt temp_cte.txt > comunes.txt

    # Borramos líneas en blanco si las hay
    sed -i '/^$/d' comunes.txt

    PARES_ENCONTRADOS=0

    echo "[3/3] Extrayendo vectores binarios..."
    while IFS= read -r comun; do
        if [ -z "$comun" ]; then continue; fi

        PARES_ENCONTRADOS=$((PARES_ENCONTRADOS + 1))

        # Encontrar en qué línea de los .txt están esos valores PSD
        lineA=$(grep -n -F "$comun" "$FILE_DFT" | head -n 1 | cut -d: -f1)
        lineB=$(grep -n -F "$comun" "$FILE_CTE" | head -n 1 | cut -d: -f1)

        # Sacar el vector binario exacto de la línea correspondiente
        vecA=$(sed -n "${lineA}p" "$FILE_COMB")
        vecB=$(sed -n "${lineB}p" "$FILE_COMB")

        # Escribir en el archivo de resultados
        echo "Par $PARES_ENCONTRADOS:" >> "$OUT_FILE"
        echo "VectorA = $vecA" >> "$OUT_FILE"
        echo "VectorB = $vecB" >> "$OUT_FILE"
        echo "--------------------------------------------------------------------------------" >> "$OUT_FILE"

    done < comunes.txt

    if [ "$PARES_ENCONTRADOS" -gt 0 ]; then
        echo "¡Éxito! Se encontraron $PARES_ENCONTRADOS pares válidos."
        echo "Guardados en: Resultados_Pares_$N.txt"
    else
        echo "No se encontraron pares válidos para N=$N."
    fi

    # Limpiar basura temporal
    rm -f temp_dft.txt temp_cte.txt comunes.txt run_log.txt
    
    # Volver a la raíz del código
    cd ../..
done

echo -e "\n¡Proceso de automatización completado al 100%!"
