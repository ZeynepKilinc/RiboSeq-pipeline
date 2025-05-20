#!/bin/bash
#SBATCH --job-name=riboseq_pipeline
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zeynep_kilinc@brown.edu
#SBATCH --output=logs/riboseq_%j.log


mkdir -p logs

module load star/2.7.10b-fj6kao2 \
	     bowtie2/2.5.3-qgscc2u \
	     subread/2.0.2-5agghnd \
	      fastqc/0.12.1-sk2rb3a 
	      
set -euo pipefail

source $HOME/scratch/ribo_seq/bin/activate
echo "Using cutadapt at $(which cutadapt)"
cutadapt --version
echo "Using STAR at $(which STAR)"
STAR --version | head -1

# ---- user‑specific paths ----
PROJECT_DIR=$HOME/scratch/ribo_seq

RAW_DIR=$PROJECT_DIR/raw_fastq          
TRIMMED_DIR=$PROJECT_DIR/trimmed_fastq  
NON_RRNA_DIR=$PROJECT_DIR/non_rrna_fastq
FASTQC_DIR=$PROJECT_DIR/fastqc_reports
GENOME_DIR=$PROJECT_DIR/reference_genome
STAR_INDEX=$PROJECT_DIR/star_index
STAR_OUT=$PROJECT_DIR/star_alignments
COUNTS_DIR=$PROJECT_DIR/counts
export RAW_DIR TRIMMED_DIR NON_RRNA_DIR FASTQC_DIR \
       GENOME_DIR STAR_INDEX STAR_OUT COUNTS_DIR
# --------------------------------

# Step 1: Trim adapters
mkdir -p $TRIMMED_DIR
echo "Found the following files in $RAW_DIR:"
ls -1 $RAW_DIR/*.fastq
for f in $RAW_DIR/*.fastq; do
    echo "Processing $f"
    base=$(basename "$f" .fastq)
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
             -m 18 \
             --trim-n \
             -o $TRIMMED_DIR/${base}_trimmed.fastq "$f"
done

# Step 2: FastQC on trimmed files
mkdir -p $FASTQC_DIR/trimmed
fastqc $TRIMMED_DIR/*.fastq.gz -o $FASTQC_DIR/trimmed

# Step 3: Extract rRNA reference
#cd $GENOME_DIR
#grep -A 1 -i "^>.* ribosomal RNA" GCF_000001405.26_GRCh38_rna.fna > human_rRNA_GRCh38.fa
#bowtie2-build human_rRNA_GRCh38.fa human_rRNA_GRCh38
#cd -

# Step 4: rRNA depletion
#mkdir -p $NON_RRNA_DIR
#for f in $TRIMMED_DIR/*.fastq.gz; do
#    base=$(basename "$f" _trimmed.fastq.gz)
#    bowtie2 -p 4 -x $GENOME_DIR/human_rRNA_GRCh38 \
#            -U "$f" \
#            --un $NON_RRNA_DIR/${base}_clean.fastq \
#            -S /dev/null
#done

# Step 5: FastQC on rRNA-depleted files
#mkdir -p $FASTQC_DIR/non_rRNA
#fastqc $NON_RRNA_DIR/*.fastq -o $FASTQC_DIR/non_rRNA

# Skip rRNA depletion; use trimmed reads directly
NON_RRNA_DIR=$TRIMMED_DIR   # Just use the trimmed reads

# Step 6: Build STAR genome index (if not already built)
mkdir -p $STAR_INDEX
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir $STAR_INDEX \
     --genomeFastaFiles $GENOME_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile $GENOME_DIR/Homo_sapiens.GRCh38.109.gtf \
     --sjdbOverhang 49 \
     --limitGenomeGenerateRAM 25000000000

# Step 7: STAR alignment (now directly from trimmed reads)
mkdir -p $STAR_OUT
for f in $NON_RRNA_DIR/*.fastq.gz; do
    base=$(basename "$f" _trimmed.fastq.gz)
    STAR --runThreadN 8 \
         --genomeDir $STAR_INDEX \
         --readFilesIn "$f" \
         --readFilesCommand zcat \
         --outFileNamePrefix $STAR_OUT/${base}_ \
         --outSAMtype BAM SortedByCoordinate
done

# Step 8: Quantify gene expression with featureCounts
mkdir -p $COUNTS_DIR
featureCounts -T 8 -s 1 \
  -a $GENOME_DIR/Homo_sapiens.GRCh38.109.gtf \
  -o $COUNTS_DIR/gene_counts.txt \
  $STAR_OUT/*Aligned.sortedByCoord.out.bam
