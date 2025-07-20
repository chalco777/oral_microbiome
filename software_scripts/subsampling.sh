#!/bin/bash
# Porcentaje de muestreo (por ejemplo, 0.1 para el 10%)
SAMPLE_NUMBER=909628
# Semilla para el muestreo aleatorio
SEED=100
OUTPUT_DIR=/media/crowfoot2/DATOS/CHALCO/met_mon/downsampling
INPUT_DIR=/media/crowfoot2/DATOS/CHALCO/met_mon/human_read_filtered
mkdir -p $OUTPUT_DIR
# Iterar sobre todos los archivos R1
for R1 in $INPUT_DIR/*_R1.fastq.gz; do
    # Obtener el nombre base y el archivo R2 correspondiente
    BASE_NAME="${R1##*/}"  # Eliminar la ruta para obtener solo el nombre de archivo
    ###*/: Esta parte indica que se elimine todo lo que está antes del último /, incluido el mismo /. Básicamente, conserva solo el nombre del archivo.
    BASE_NAME="${BASE_NAME%_R1.fastq.gz}"
    #%_R1.fastq.gz: Indica que se elimine la subcadena _R1.fastq.gz al final del valor de BASE_NAME
    R2="${INPUT_DIR}/${BASE_NAME}_R2.fastq.gz"

    # Verificar si el archivo R2 existe
    if [ -e "$R2" ]; then
        echo "Procesando $R1 y $R2"

        # Aplicar el muestreo a las lecturas intercaladas
        seqtk sample -s"$SEED" $R1 "$SAMPLE_NUMBER" | gzip > "${OUTPUT_DIR}/${BASE_NAME}_R1_sampled_fastq.gz"
        seqtk sample -s"$SEED" $R2 "$SAMPLE_NUMBER" | gzip > "${OUTPUT_DIR}/${BASE_NAME}_R2_sampled_fastq.gz"

    else
        echo "No se encontró el archivo emparejado para $R1"
    fi
done
