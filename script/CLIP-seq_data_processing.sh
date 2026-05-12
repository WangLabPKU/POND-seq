#!/bin/bash
# ==============================================================================
# CLIP-seq Analysis Pipeline: UMI Extraction, Alignment, and Peak Calling
# ==============================================================================

echo "Starting CLIP-seq data processing pipeline..."

# ------------------------------------------------------------------------------
# Step 1: Extract Unique Molecular Identifiers (UMIs)
# ------------------------------------------------------------------------------
# Extract a 10 bp random UMI from the reads and append it to the read name.
# This is crucial for distinguishing biological duplicates from PCR duplicates later.
# Note: If your data does NOT contain UMIs (standard CLIP-seq), skip this step.

# Process IP replicate 1
umi_tools extract --random-seed 1 \
    --bc-pattern NNNNNNNNNN \
    --log IP1.metrics \
    -I fastq/XXXXXXXXXXX.fastq.gz \
    -S umifq/IP1.fq.gz 

# Process IP replicate 2
umi_tools extract --random-seed 1 \
    --bc-pattern NNNNNNNNNN \
    --log IP2.metrics \
    -I fastq/XXXXXXXXXXX.fastq.gz \
    -S umifq/IP2.fq.gz 

# Process Input (Control)
umi_tools extract --random-seed 1 \
    --bc-pattern NNNNNNNNNN \
    --log input.metrics \
    -I fastq/XXXXXXXXXXX.fastq.gz \
    -S umifq/input.fq.gz 


# ------------------------------------------------------------------------------
# Step 2: Quality Control and Adapter Trimming
# ------------------------------------------------------------------------------
# Trim low-quality bases and adapters from the UMI-extracted reads.

trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
    umifq/IP1.fq.gz --gzip -o ./clean/ --basename IP1

trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
    umifq/IP2.fq.gz --gzip -o ./clean/ --basename IP2

trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 \
    umifq/input.fq.gz --gzip -o ./clean/ --basename input


# ------------------------------------------------------------------------------
# Step 3: Read Alignment to the Reference Genome
# ------------------------------------------------------------------------------
# Map the trimmed reads to the human genome. 
# Note: --outFilterMultimapNmax 1 ensures only uniquely mapped reads are kept, 
# which is standard practice for precise binding site identification.

for sample in IP1 IP2 input
do
    echo "Aligning sample: ${sample}..."
    STAR --runThreadN 20 \
        --runMode alignReads \
        --genomeDir Reference/GRCh38/index_v40.pri \
        --readFilesCommand zcat \
        --outFilterMultimapNmax 1 \
        --outSAMattrRGline ID:${sample} SM:${sample} PU:Illumina \
        --outSAMattributes All \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix align/${sample}_ \
        --readFilesIn clean/${sample}_trimmed.fq.gz 
done


# ------------------------------------------------------------------------------
# Step 4: PCR Deduplication
# ------------------------------------------------------------------------------
# Note for non-UMI data: Use Picard MarkDuplicates instead of umi_tools.
# picard MarkDuplicates -I align/XXX_Aligned.sortedByCoord.out.bam \
#     -O align/XXX_Aligned_dedupped.bam --REMOVE_DUPLICATES true \
#     --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
#     -M align/XXX_Aligned_MarkDuplicates_output.metrics


# Index the sorted BAM files before deduplication
samtools index -@ 20 align/input_Aligned.sortedByCoord.out.bam
samtools index -@ 20 align/IP1_Aligned.sortedByCoord.out.bam
samtools index -@ 20 align/IP2_Aligned.sortedByCoord.out.bam

# Perform deduplication
echo "Performing UMI-based deduplication..."
umi_tools dedup -I align/input_Aligned.sortedByCoord.out.bam \
    -S align/input_Aligned.sortedByCoord.umidedup.bam

umi_tools dedup -I align/IP1_Aligned.sortedByCoord.out.bam \
    -S align/IP1_Aligned.sortedByCoord.umidedup.bam

umi_tools dedup -I align/IP2_Aligned.sortedByCoord.out.bam \
    -S align/IP2_Aligned.sortedByCoord.umidedup.bam

# Index the deduplicated BAM files
samtools index -@ 20 align/input_Aligned.sortedByCoord.umidedup.bam
samtools index -@ 20 align/IP1_Aligned.sortedByCoord.umidedup.bam
samtools index -@ 20 align/IP2_Aligned.sortedByCoord.umidedup.bam


# ------------------------------------------------------------------------------
# Step 5: Peak Calling and Crosslink Site Identification (PureCLIP)
# ------------------------------------------------------------------------------
# Use PureCLIP to identify RBP crosslink sites by incorporating IP and Input data.

echo "Running PureCLIP for peak calling..."
pureclip -i align/IP1_Aligned.sortedByCoord.umidedup.bam \
         --bai align/IP1_Aligned.sortedByCoord.umidedup.bam.bai \
         -i align/IP2_Aligned.sortedByCoord.umidedup.bam \
         --bai align/IP2_Aligned.sortedByCoord.umidedup.bam.bai \
         -ibam align/input_Aligned.sortedByCoord.umidedup.bam \
         --ibai align/input_Aligned.sortedByCoord.umidedup.bam.bai \
         -fk -g Reference/GRCh38/GRCh38.p13.genome.pri.fa \
         -nt 20 \
         -nta 24 \
         -iv "chr1;chr2;chr3;" \
         -oa -o peak/PureCLIP.crosslink_sites.bed \
         -or peak/PureCLIP.binding_regions.bed 

echo "Pipeline execution completed!"

# Annotate crosslink sites with gene and exon information
bedtools intersect -a peak/PureCLIP.crosslink_sites.bed -b Reference/GRCh38/gencode.v40.pri.annotation_gene.bed -s -wa -wb > peak/PureCLIP.crosslink_sites_in_gene.bed
bedtools intersect -a peak/PureCLIP.crosslink_sites.bed -b Reference/GRCh38/gencode.v40.pri.annotation_exon.bed -s -wa -wb > peak/PureCLIP.crosslink_sites_in_exon.bed