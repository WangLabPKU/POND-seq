#!/bin/bash
# ==============================================================================
# Pipeline for Sequencing Data Preprocessing: QC, Alignment, and Quantification
# ==============================================================================

# Loop through all sample directories in the 'fastq' folder
for sample in fastq/*
do 
    # Extract the base name of the sample directory
    sample1=`basename $sample`
    
    echo "Processing sample: ${sample1} ..."

    # --------------------------------------------------------------------------
    # Step 1: Quality Control and Adapter Trimming using Trim Galore
    # --------------------------------------------------------------------------
    trim_galore -j 20 \
        -q 20 \
        --phred33 \
        --stringency 3 \
        --length 20 \
        -e 0.1 \
        --paired fastq/${sample1}/${sample1}_raw_1.fq.gz fastq/${sample1}/${sample1}_raw_2.fq.gz \
        --gzip \
        -o ./clean/ \
        --basename ${sample1}
        
    # Parameter details for Trim Galore:
    # -j 20: Use 20 cores for processing.
    # -q 20: Trim low-quality ends with a Phred score below 20.
    # --phred33: Instructs the tool to use ASCII+33 quality scores.
    # --stringency 3: Overlap with adapter sequence required to trim a read.
    # --length 20: Discard reads that become shorter than 20 bp after trimming.
    # -e 0.1: Maximum allowed error rate.
    # --paired: Perform paired-end validation (discards the pair if one read is too short).
    # --gzip: Compress the output files to save storage.

    # --------------------------------------------------------------------------
    # Step 2: Read Alignment using STAR
    # --------------------------------------------------------------------------
    STAR --runThreadN 24 \
        --runMode alignReads \
        --genomeDir Reference/GRCh38/index_ERCC \
        --readFilesCommand zcat \
        --outFilterMismatchNoverLmax 0.04 \
        --outFilterType BySJout \
        --outSAMattrRGline ID:${sample1} SM:${sample1} PU:Illumina \
        --outSAMattributes All \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix align/${sample1}_ \
        --readFilesIn clean/${sample1}_val_1.fq.gz clean/${sample1}_val_2.fq.gz \
        --outReadsUnmapped Fastx 
        
    # Parameter details for STAR:
    # --runThreadN 24: Use 24 threads for alignment.
    # --genomeDir: Path to the STAR reference genome index (GRCh38 + ERCC spike-ins).
    # --readFilesCommand zcat: Uncompress gzipped FASTQ files on the fly.
    # --outFilterMismatchNoverLmax 0.04: Alignment is output only if the ratio of mismatches to mapped length is <= 4%.
    # --outFilterType BySJout: Reduces spurious junctions.
    # --outSAMattrRGline: Add Read Group (RG) headers, essential for downstream variant calling or specific tools.
    # --outSAMtype BAM SortedByCoordinate: Output a coordinate-sorted BAM file directly, saving an extra samtools sort step.
    # --outReadsUnmapped Fastx: Output unmapped reads into separate FASTQ files for potential downstream troubleshooting.

done

echo "Alignment completed. Starting feature quantification..."

# --------------------------------------------------------------------------
# Step 3: Read Quantification using featureCounts (Subread package)
# --------------------------------------------------------------------------
featureCounts -a GRCh38.p13.gencode.v40.chr_ERCC.gtf \
    -g gene_id \
    -t exon \
    -T 20 \
    -O \
    -o rawCounts.tsv \
    -p \
    --fraction \
    --primary \
    -Q 30 \
    align/*_Aligned.sortedByCoord.out.bam
    
# Parameter details for featureCounts:
# -a: Path to the annotation file (GTF format, matching the GRCh38 and ERCC references).
# -g gene_id: Specify the attribute type used to group features (summarizing at the gene level).
# -t exon: Specify the feature type to count reads against (exons).
# -T 20: Use 20 threads.
# -O: Assign reads to all their overlapping meta-features (allows multi-overlap).
# -o rawCounts.tsv: Name of the output count matrix file.
# -p: Specify that the sequencing data is paired-end.
# --fraction: Assign fractional counts to multi-mapping reads.
# --primary: Count primary alignments only (ignore secondary alignments).
# -Q 30: Minimum mapping quality (MAPQ) score required for a read to be counted.