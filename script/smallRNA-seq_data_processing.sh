#!/bin/bash
# ==============================================================================
# Advanced RNA-Seq Pipeline: Custom ncRNA/miRNA Profiling and Sequential Alignment
# ==============================================================================

# ------------------------------------------------------------------------------
# Module 1: Build Custom Reference Sequences for Small ncRNA
# ------------------------------------------------------------------------------
# Download the human ncRNA FASTA file from Ensembl
wget https://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
gunzip Homo_sapiens.GRCh38.ncrna.fa.gz

# Filter out long non-coding RNAs (lncRNA, lincRNA) and non-standard scaffolds
# to create a specific reference index for small ncRNAs
seqkit seq Homo_sapiens.GRCh38.ncrna.fa -n | \
    grep -v lncRNA | \
    grep -v lincRNA | \
    grep -v "CHR_" | \
    grep -v "scaffold:" > filter_list

# Extract the filtered sequences based on the list
seqkit grep -n -f filter_list Homo_sapiens.GRCh38.ncrna.fa -o Homo_sapiens.GRCh38.small.ncRNA.fa

# Build the Bowtie index for the filtered small ncRNA reference
bowtie-build Homo_sapiens.GRCh38.small.ncRNA.fa Homo_sapiens.GRCh38.small.ncRNA.fa


# ------------------------------------------------------------------------------
# Module 2: Adapter Trimming and Alignment to Small ncRNA (e.g., snoRNA)
# ------------------------------------------------------------------------------
for sample in fastq/*
do
    samplename=`basename $sample`
    echo "Processing ${samplename} for small ncRNA alignment..."

    # Step 2.1: Trim adapters and low-quality bases
    # Specific adapters (-a and -a2) are provided for the library preparation
    trim_galore -j 20 -q 20 --phred33 --stringency 3 --length 15 -e 0.1 --paired \
        -a AGATCGGAAGAGC -a2 GATCGTCGGA \
        fastq/${samplename}/${samplename}_raw_1.fq.gz fastq/${samplename}/${samplename}_raw_2.fq.gz \
        -o ./clean --basename ${samplename} --dont_gzip
    
    # Step 2.2: Map trimmed reads to the custom small ncRNA reference using Bowtie
    # --un: Output reads that fail to align for downstream genome alignment
    # -a --best --strata: Report all valid alignments in the best stratum
    bowtie -a --best --strata -n 1 \
        -x Homo_sapiens.GRCh38.small.ncRNA.fa \
        -q -1 clean/${samplename}_val_1.fq -2 clean/${samplename}_val_2.fq \
        -S align/${samplename}.sam \
        --un align/${samplename}.unaligned.fq \
        --allow-contain -p 40 
    
    # Step 2.3: Sort BAM and calculate fractional counts for multi-mapping reads
    # The awk script parses the XM:i tag to distribute counts equally (1/N) among all mapped loci
    samtools view -hb -F 4 align/${samplename}.sam -@ 20 | \
        samtools sort -@ 20 - -o align/${samplename}.sorted.bam
        
    samtools view align/${samplename}.sorted.bam | \
        awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]} END {for(i in x) print i, x[i]}' | \
        sort -gr -k 2 > ${samplename}.counts.txt
done


# ------------------------------------------------------------------------------
# Module 3: Alignment and Quantification against Mature microRNAs (miRNAs)
# ------------------------------------------------------------------------------
for sample in fastq/*
do
    samplename=`basename $sample`
    echo "Quantifying mature miRNAs for ${samplename}..."

    # Map the original cleaned reads directly to the mature miRNA database using Bowtie
    bowtie -a --best --strata -n 1 \
        -x hsa_mature.fa \
        -q -1 clean/${samplename}_val_1.fq -2 clean/${samplename}_val_2.fq \
        -S align_to_mature/${samplename}.sam \
        --un align_to_mature/${samplename}.unaligned.fq \
        --allow-contain -p 40 
        
    # Sort BAM and calculate fractional counts (1/N) for mature miRNAs
    samtools view -hb -F 4 align_to_mature/${samplename}.sam -@ 20 | \
        samtools sort -@ 20 - -o align_to_mature/${samplename}.sorted.bam
        
    samtools view align_to_mature/${samplename}.sorted.bam | \
        awk -v OFS="\t" '{split($0,a,"XM:i:"); x[$3]+=1/a[2]} END {for(i in x) print i, x[i]}' | \
        sort -gr -k 2 > ${samplename}.maturecounts.txt
done

## For hairpin miRNAs, we can use the same approach as above but with a different reference index (e.g., hsa_hairpin.fa) to quantify hairpin miRNAs. 


# ------------------------------------------------------------------------------
# Module 4: Align Unmapped Reads to the Whole Genome
# ------------------------------------------------------------------------------
# Process the reads that did NOT map to the small ncRNA database
for sample in fastq/*
do 
    sample1=`basename $sample`
    echo "Aligning unmapped reads of ${sample1} to GRCh38 genome..."

    STAR --runThreadN 24 \
        --runMode alignReads \
        --genomeDir Reference/GRCh38/index_ERCC \
        --readFilesCommand cat \
        --outFilterMismatchNoverLmax 0.04 \
        --outFilterType BySJout \
        --outSAMattrRGline ID:${sample1} SM:${sample1} PU:Illumina \
        --outSAMattributes All \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix align_to_genome/${sample1}_ \
        --readFilesIn align/${sample1}.unaligned_1.fq align/${sample1}.unaligned_2.fq
done