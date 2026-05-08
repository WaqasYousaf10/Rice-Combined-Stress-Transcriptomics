#!/bin/bash
# =============================================================================
# RNA-seq Data Processing Pipeline for Rice Combined Stress Experiment
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: bash 01_linux_pipeline.sh
# =============================================================================

set -e

# =============================================================================
# CONFIGURATION - UPDATE THESE PATHS
# =============================================================================

WORKING_DIR="/path/to/your/RNAseq_data"
REF_GENOME="/path/to/rice_genome/IRGSP-1.0.fa"
REF_GTF="/path/to/rice_genome/IRGSP-1.0.gtf"
HISAT2_INDEX="/path/to/hisat2_index/rice_index"
THREADS=8
ADAPTERS="/path/to/adapters/TruSeq3-PE.fa"

# =============================================================================
# CREATE DIRECTORIES
# =============================================================================

mkdir -p ${WORKING_DIR}/{01_raw_data,02_trimmed,03_aligned,04_counts,05_reports}

# =============================================================================
# STEP 1: QUALITY CONTROL (FASTQC)
# =============================================================================

echo "Running FASTQC..."
for file in ${WORKING_DIR}/01_raw_data/*.fastq.gz; do
    fastqc $file -o ${WORKING_DIR}/05_reports/
done

# =============================================================================
# STEP 2: TRIMMING (TRIMMOMATIC)
# =============================================================================

echo "Running Trimmomatic..."
for R1 in ${WORKING_DIR}/01_raw_data/*_R1.fastq.gz; do
    SAMPLE=$(basename ${R1} _R1.fastq.gz)
    R2=${WORKING_DIR}/01_raw_data/${SAMPLE}_R2.fastq.gz
    
    TrimmomaticPE \
        -threads ${THREADS} \
        ${R1} ${R2} \
        ${WORKING_DIR}/02_trimmed/${SAMPLE}_R1_trimmed.fastq.gz \
        ${WORKING_DIR}/02_trimmed/${SAMPLE}_R1_unpaired.fastq.gz \
        ${WORKING_DIR}/02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz \
        ${WORKING_DIR}/02_trimmed/${SAMPLE}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

# =============================================================================
# STEP 3: ALIGNMENT (HISAT2)
# =============================================================================

echo "Running HISAT2 alignment..."
for R1 in ${WORKING_DIR}/02_trimmed/*_R1_trimmed.fastq.gz; do
    SAMPLE=$(basename ${R1} _R1_trimmed.fastq.gz)
    R2=${WORKING_DIR}/02_trimmed/${SAMPLE}_R2_trimmed.fastq.gz
    
    hisat2 -p ${THREADS} \
        -x ${HISAT2_INDEX} \
        -1 ${R1} -2 ${R2} \
        -S ${WORKING_DIR}/03_aligned/${SAMPLE}.sam
done

# =============================================================================
# STEP 4: CONVERT SAM TO BAM
# =============================================================================

echo "Converting SAM to BAM..."
for SAM in ${WORKING_DIR}/03_aligned/*.sam; do
    SAMPLE=$(basename ${SAM} .sam)
    samtools view -S -b ${SAM} > ${WORKING_DIR}/03_aligned/${SAMPLE}.bam
    samtools sort ${WORKING_DIR}/03_aligned/${SAMPLE}.bam -o ${WORKING_DIR}/03_aligned/${SAMPLE}_sorted.bam
    samtools index ${WORKING_DIR}/03_aligned/${SAMPLE}_sorted.bam
done

# =============================================================================
# STEP 5: QUANTIFICATION (FEATURECOUNTS)
# =============================================================================

echo "Running featureCounts..."
featureCounts -p -T ${THREADS} -a ${REF_GTF} -o ${WORKING_DIR}/04_counts/count_matrix.txt \
    ${WORKING_DIR}/03_aligned/*_sorted.bam

# Create count matrix for R
cut -f 1,7- ${WORKING_DIR}/04_counts/count_matrix.txt > ${WORKING_DIR}/04_counts/count_matrix_for_R.csv

echo "Pipeline completed successfully!"
