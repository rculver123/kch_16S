#!/bin/bash
#
#SBATCH --job-name=16Spipeline
#
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-user=<YOUR EMAIL>
#SBATCH --mail-type=ALL
#SBATCH -p kchuang

# Making folders
mkdir fwd_reads
mkdir fwd_reads/samples
mkdir rev_reads
mkdir rev_reads/samples
mv *R1* fwd_reads/samples/
mv *R2* rev_reads/samples/

# Optional:
#gunzip fwd_reads/samples/* 
#gunzip rev_reads/samples/*

# Rename files
cd fwd_reads/samples
rename '_R1_001.fastq' '.fastq' *
for file in *; do mv "$file" "${file/*_S/S}"; done
cd ../../rev_reads/samples
rename '_R2_001.fastq' '.fastq' *
for file in *; do mv "$file" "${file/*_S/S}"; done
cd ../..

# Run DADA2 pooled and make a tree using TreeDA
ml R
Rscript dada2_analysis.R

