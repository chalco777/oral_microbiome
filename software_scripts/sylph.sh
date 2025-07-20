#!/bin/bash

# Activar el entorno de conda para Sylph
echo "Activando el entorno de conda para Sylph..."
source ~/anaconda3/etc/profile.d/conda.sh
conda activate sylph || { echo "ERROR: No se pudo activar el entorno Conda"; exit 1; }

# Definir el directorio de trabajo y el número de hilos
WORK_DIR="/media/crowfoot2/DATOS/CHALCO/met_mon/downsampling"  # Ajusta esta ruta según sea necesario
THREADS=15

# Moverse al directorio de trabajo
cd "$WORK_DIR" || { echo "ERROR: No se pudo acceder al directorio de trabajo $WORK_DIR"; exit 1; }

# Definir la ruta a la base de datos GTDB
GTDB_DB="/home/crowfoot2/databases/gtdb-r220/gtdb-r220-c200-dbv1.syldb"
GTDB_MET="/home/crowfoot2/databases/gtdb-r220/gtdb_r220_metadata.tsv.gz"

# Crear directorios de salida
OUTPUT_DIR="$WORK_DIR/sylph_results"
mkdir -p "$OUTPUT_DIR"
 
# Realizar el perfilado contra la base de datos GTDB-R220
echo "Perfilando contra la base de datos GTDB-R220..."
sylph profile "$GTDB_DB" -1 *_R1_sampled_fastq.gz -2 *_R2_sampled_fastq.gz -t "$THREADS" -o "$OUTPUT_DIR/results.tsv"

sylph profile -u "$GTDB_DB" -1 *_R1_sampled_fastq.gz -2 *_R2_sampled_fastq.gz -t "$THREADS" -o "$OUTPUT_DIR/results_unknown.tsv"

# Generar perfiles taxonómicos
echo "Generando perfiles taxonómicos..."
mkdir -p $OUTPUT_DIR/taxon
python /media/crowfoot2/DATOS/CHALCO/met_mon/programs/sylph-utils/sylph_to_taxprof.py -s "$OUTPUT_DIR/results.tsv" -m $GTDB_MET -o "$OUTPUT_DIR/taxon/down_"
python /media/crowfoot2/DATOS/CHALCO/met_mon/programs/sylph-utils/sylph_to_taxprof.py -s "$OUTPUT_DIR/results_unknown.tsv" -m $GTDB_MET -o "$OUTPUT_DIR/taxon/unknown_"

echo "Análisis completo. Resultados guardados en $OUTPUT_DIR"

# Contar la cantidad de lecturas para cada par de archivos FASTQ y guardarlas en un archivo TSV
echo "Contando lecturas para cada archivo FASTQ..."

echo -e "Muestra\tLecturas_R1\tLecturas_R2" > "$OUTPUT_DIR/countreads.tsv"

SAMPLES=($(ls *_R1_sampled_fastq.gz | sed 's/_R1_sampled_fastq.gz//'))

for SAMPLE in "${SAMPLES[@]}"; do
    R1_FILE="${SAMPLE}_R1_sampled_fastq.gz"
    R2_FILE="${SAMPLE}_R2_sampled_fastq.gz"
    
    if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
        # Contar lecturas en R1
        READ_COUNT_R1=$(gunzip -c "$R1_FILE" | wc -l)
        READ_COUNT_R1=$((READ_COUNT_R1 / 4))
        
        # Contar lecturas en R2
        READ_COUNT_R2=$(gunzip -c "$R2_FILE" | wc -l)
        READ_COUNT_R2=$((READ_COUNT_R2 / 4))
        
        echo -e "$SAMPLE\t$READ_COUNT_R1\t$READ_COUNT_R2" >> "$OUTPUT_DIR/countreads.tsv"
    else
        echo "No se encontraron archivos para la muestra $SAMPLE"
    fi
done

echo "Conteo de lecturas guardado en $OUTPUT_DIR/countreads.tsv"

# Desactivar el entorno de conda
conda deactivate

echo "Procesamiento con Sylph completado."
