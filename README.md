# SNP Calling Analysis Workflow

This repository contains a detailed guide and scripts for performing SNP calling analysis. The workflow includes quality control, genome indexing, read alignment, variant calling, and filtering steps. Below are the tools used and the corresponding steps in the analysis.

---

## Tools Used
- **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):** Quality control of raw sequencing reads.
- **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic):** Trimming and filtering of raw reads to remove low-quality bases and adapter sequences.
- **[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml):** Genome indexing and alignment of reads to the reference genome.
- **[Samtools](http://www.htslib.org/):** Processing and manipulation of alignment files.
- **[GATK](https://gatk.broadinstitute.org/hc/en-us):** Variant calling and SNP analysis.

---

## Workflow Overview

### 1. Quality Control with FastQC
- **Objective**: Assess the quality of raw sequencing data.
- **Command**:
  ```bash
  fastqc -o qc_reports raw_data/*.fastq.gz
  ```

### 2. Read Trimming with Trimmomatic
- **Objective**: Remove low-quality bases and adapter sequences.
- **Command**:
  ```bash
  trimmomatic PE -threads 4 \
    raw_data/sample_R1.fastq.gz raw_data/sample_R2.fastq.gz \
    trimmed_data/sample_R1_paired.fastq.gz trimmed_data/sample_R1_unpaired.fastq.gz \
    trimmed_data/sample_R2_paired.fastq.gz trimmed_data/sample_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
  ```

### 3. Genome Indexing with Bowtie2
- **Objective**: Build a genome index for the reference genome.
- **Command**:
  ```bash
  bowtie2-build reference/genome.fa reference/genome_index
  ```

### 4. Read Alignment with Bowtie2
- **Objective**: Align reads to the reference genome.
- **Command**:
  ```bash
  bowtie2 -x reference/genome_index \
    -1 trimmed_data/sample_R1_paired.fastq.gz \
    -2 trimmed_data/sample_R2_paired.fastq.gz \
    -S alignments/sample.sam
  ```

### 5. Conversion and Sorting with Samtools
- **Objective**: Convert SAM to BAM, sort, and index the alignments.
- **Commands**:
  ```bash
  samtools view -Sb alignments/sample.sam > alignments/sample.bam
  samtools sort alignments/sample.bam -o alignments/sample_sorted.bam
  samtools index alignments/sample_sorted.bam
  ```

### 6. Variant Calling with GATK
- **Objective**: Call variants and identify SNPs.
- **Steps**:
  1. **HaplotypeCaller**:
     ```bash
     gatk HaplotypeCaller \
         -R reference/genome.fa \
         -I alignments/sample_sorted.bam \
         -O variants/sample.g.vcf.gz \
         -ERC GVCF
     ```
  2. **Combine GVCFs**:
     ```bash
     gatk CombineGVCFs \
         -R reference/genome.fa \
         --variant variants/sample1.g.vcf.gz \
         --variant variants/sample2.g.vcf.gz \
         -O variants/combined.g.vcf.gz
     ```
  3. **GenotypeGVCFs**:
     ```bash
     gatk GenotypeGVCFs \
         -R reference/genome.fa \
         -V variants/combined.g.vcf.gz \
         -O variants/cohort.vcf.gz
     ```

---

## Repository Structure
- `raw_data/`: Contains raw sequencing files.
- `trimmed_data/`: Contains trimmed and filtered sequencing files.
- `qc_reports/`: Contains quality control reports.
- `genome_index/`: Reference genome and index files.
- `alignments/`: Contains alignment files in BAM/SAM format.
- `variants/`: Contains GVCF and VCF files with variant data.
- `scripts/`: Contains bash scripts for each step of the workflow.

---

## Folder Structure
```plaintext
SNP-Calling-Analysis/
├── raw_data/         # Raw FASTQ files
├── trimmed_data/     # Trimmed FASTQ files
├── qc_reports/       # Quality control reports
├── genome_index/     # Reference genome and index files
├── alignments/       # BAM files and indexes
├── variants/         # GVCF and VCF files
├── scripts/          # Scripts for analysis steps
└── results/          # Final filtered SNP data
```

---

## Usage
Clone this repository and follow the steps in the workflow:
```bash
git clone https://github.com/Khaled-H-Mousa/snp-calling-workflow.git
cd snp-calling-workflow
bash scripts/run_workflow.sh
```

---

## References
1. [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us)
2. [Bowtie2 Documentation](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
3. [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

---

## Author
[Khaled-H-Mousa](https://github.com/Khaled-H-Mousa)


