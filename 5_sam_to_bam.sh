#!/bin/bash

# Define the input and output directories
path=$(pwd) 
INPUT_DIR="${path}/Aligned_OUT"  # Directory containing SAM files
OUTPUT_DIR="${path}/Processed_BAM"  # Directory to save BAM, sorted BAM, and index files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"  

# Read sample names from a file or list
filesList=$(cat fastQlist | sort | uniq)

# Loop through each sample name in the list
echo "$filesList" | while read file; do
   
    
    # Step 1: Convert SAM to BAM
    samtools view -@ 8 -Sb -o "${OUTPUT_DIR}/${file}.bam" "${INPUT_DIR}/${file}.sam"
    
    # Step 2: Sort the BAM file
    samtools sort -@ 8 -o "${OUTPUT_DIR}/${file}_sorted.bam" "${OUTPUT_DIR}/${file}.bam"

    # Step 3: Index the sorted BAM file
    samtools index "${OUTPUT_DIR}/${file}_sorted.bam"
    
    # Step 4: Generate alignment statistics
    #samtools flagstat "${OUTPUT_DIR}/${file}_sorted.bam" > "${OUTPUT_DIR}/${file}_mapping_stats.txt"

    # Clean up intermediate files (optional)
    #rm "${OUTPUT_DIR}/${file}.bam"
    
    echo "Finished processing sample: $file"
done

