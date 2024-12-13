#!/bin/bash
# Define the input and output directories
path=$(pwd)
INPUT_DIR="${path}/Combine_VCT"  # Directory containing BAM files
OUTPUT_DIR="${path}/GenotypeGVCFs_VCT"  # Directory to save VCF files
genomeName="${path}/Genomic_Data/Sesamum_indicum.fasta"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Read sample names from a file or list
filesList=$(cat fastQlist | sort | uniq)

# Loop through each sample name in the list
echo "$filesList" | while read file; do

    # Define input BAM file
    inputBamFile="${INPUT_DIR}/combined.g.vcf.gz"

    # Define output VCF file
    outputVcfFile="${OUTPUT_DIR}/joint_genotyped.vcf.gz"

    # Run GATK GenotypeGVCFs
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
        -R "$genomeName" \
        -V "$inputBamFile" \
        -O "$outputVcfFile"

done



