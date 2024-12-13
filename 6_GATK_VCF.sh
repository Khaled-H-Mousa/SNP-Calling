#!/bin/bash
# Define the input and output directories
path=$(pwd)
INPUT_DIR="${path}/Processed_BAM"  # Directory containing BAM files
OUTPUT_DIR="${path}/VCT_OUT"  # Directory to save VCF files
genomeName="${path}/Genomic_Data/Sesamum_indicum.fasta"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Read sample names from a file or list
filesList=$(cat fastQlist | sort | uniq)

# Loop through each sample name in the list
echo "$filesList" | while read file; do
    # Define input BAM file
    inputBamFile="${INPUT_DIR}/${file}_fixed.bam"

    # Define output VCF file
    outputVcfFile="${OUTPUT_DIR}/${file}.vcf.gz"

    #-ERC GVCF \
    # Run GATK HaplotypeCaller for variant calling
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true" HaplotypeCaller \
        -R "$genomeName" \
        -I "$inputBamFile" \
        -O "$outputVcfFile" \
        --emit-ref-confidence GVCF 

done

