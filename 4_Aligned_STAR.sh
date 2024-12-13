#!/bin/bash

# Define the output directory
path=$(pwd) 
INPUT_DIR="${path}/trimmed_filtered"  # Ensure there's no trailing slash
OUTPUT_DIR="${path}/Aligned_OUT"  # Replace with your actual output directory
GENOME_DIR="${path}/Genomic_Data"  # Ensure this points to your genome directory

# Make sure the output directory exists
mkdir -p "$OUTPUT_DIR"  # Create the output directory if it doesn't exist

# Read sample names from fastQlist file
filesList=$(cat fastQlist | sort | uniq)
#echo "Sample list: $filesList"

# Loop through each sample name in the list
echo "$filesList" | while read file; do
    # Construct input file names
    inputFileR1="${INPUT_DIR}/${file}_R1.trimmed.fastq.gz"
    inputFileR2="${INPUT_DIR}/${file}_R2.trimmed.fastq.gz"

        # Define output SAM file
        outputSamFile="${OUTPUT_DIR}/${file}.sam"

#bowtie2 --no-unal -p n -x index_name -1 reads_1.fastq -2 reads_2.fastq -S output.sam

        # Run Bowtie2 alignment
        bowtie2 -p 4 -x "${GENOME_DIR}/Sesamum_index" \
            -1 "$inputFileR1" -2 "$inputFileR2" -S "$outputSamFile"

done

