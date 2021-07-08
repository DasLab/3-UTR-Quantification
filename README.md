 # Research Project RNA UTR Isoforms
Purpose of the project is to identify 3’ UTR isoforms from 3’ RNA sequencing data

## Installation

```bash
python, anaconda, cutadapt?
brew install samtools
brew install bedtools
brew install wget
```

## Comands for Pipeline

```bash
# Download [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

# Download [Fasta and GTF files](https://uswest.ensembl.org/Rattus_norvegicus/Info/Index)

# Make Folder Containing Data:
mkdir APA analysis

# Download Rnor6 reference genome (experimental data)

# Build bowtie index for genome
# Input rat genome, output bowtie index
bowtie2-build ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rnor_6

# Align FastQ to genome using bowtie
# input 
bowtie2 -x ref/Rnor_6 -U fastq/Proces_cleaned.fastq -b aligned/projections.bam
bowtie2 -x ref/Rnor_6 -U fastq/Proces_cleaned.fastq -S aligned/projections.sam

# Convert SAM to BAM File using Samtools:
samtools view -S -b aligned/projections.sam > aligned/projections.bam

# Convert BAM to BED File using Bedtools:
bedtools bamtobed -i aligned/projections.bam > aligned/projections.bed

# Sort BED file using Bedtools:
command: sortBed -i aligned/projections.bed > aligned/sortedprojections.bed

# Filter GTF file for stop codon annotations:
command: grep "stop_codon" ref/Rattus_norvegicus.Rnor_6.0.104.gtf > ref/rnor6_stopcodons.gtf

# Use [Bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html) closest command to find closest upstream stop codon:
command: bedtools closest -D ref -fu -a aligned/sortedprojections.bed -b ref/rnor6_stopcodons_sorted.gtf > closest/stopcodonsprojections.bed

# Sort GTF file:
command: sortBed -i ref/Rattus_norvegicus.Rnor_6.0.104.gtf > ref/rnor6sorted.gtf (bedtools accept gtf)

# Grep filtering stop codons for sorted:
grep "stop_codon" ref/rnor6sorted.gtf > ref/rnor6_stopcodons_sorted.gtf



```

