#!/bin/bash

# Define the input and output directories
path=$(pwd)
INPUT_DIR="${path}/VCT_OUT"  # Directory containing VCF files
OUTPUT_DIR="${path}/Combine_VCT"  # Directory to save combined VCF file
GENOME_REFERENCE="${path}/Genomic_Data/Sesamum_indicum.fasta"  # Reference genome file

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define the output combined VCF file
COMBINED_VCF_FILE="${OUTPUT_DIR}/combined.g.vcf.gz"

# Initialize an empty variable to hold the -V arguments for GATK
VCF_ARGUMENTS=""

# Loop through all .vcf.gz files in the input directory
for vcf_file in ${INPUT_DIR}/*.vcf.gz; do
    # Add each VCF file as a -V argument for GATK
    VCF_ARGUMENTS+="-V ${vcf_file} "
done

# Run GATK CombineGVCFs to combine all VCF files
gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CombineGVCFs \
    -R "$GENOME_REFERENCE" \
    $VCF_ARGUMENTS \
    -O "$COMBINED_VCF_FILE"

done
