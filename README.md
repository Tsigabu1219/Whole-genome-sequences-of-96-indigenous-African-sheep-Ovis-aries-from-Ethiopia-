# Whole Genome Sequencing Pipeline (96 Ethiopian Sheep)

This repository contains a complete WGS analysis pipeline including:

## Steps
1. FastQC quality control
2. MultiQC summary
3. Mapping with BWA
4. Sorting and duplicate marking (GATK)
5. Base Quality Score Recalibration (BQSR)
6. Variant calling (HaplotypeCaller)
7. Joint genotyping
8. VQSR filtering
9. SNP density analysis (CMplot)

## Tools
- BWA
- SAMtools
- GATK 4.3
- FastQC
- MultiQC
- R (CMplot, vcfR)


