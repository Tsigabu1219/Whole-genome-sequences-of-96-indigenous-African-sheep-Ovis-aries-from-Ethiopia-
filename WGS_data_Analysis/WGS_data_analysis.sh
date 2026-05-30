Code availability
###########################################################################
1. Raw sequence Quality check through FastQC and MultiQC
module load fastqc/0.11.9
fastqc -t 12 *.fastq.gz
2.MultiQC (summary of the fastqc results)
module load multiqc/1.8
multiqc *_fastqc.zip
Variant discovery
#!/bin/bash/
Create outputs folders
#^^^^^^^^^^^^^^^^^^^^^^^^^
mkdir -p ${output1}${BREED}/{BAMs,MDBAMs,RECTABs,RECBAMs,GVCFs,metrices,TMP{1,2},LOGs}
3. loading tools/modules
module load bwa/0.7.17
module load samtools/1.11
module load java/8
#2. define variables
REF="/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"
BREED="ABERGELLE"
output1="/home/output1/"
GATK="/home/apps/gatk/4.3.0.0/"
t="12"
TH="4"
dbSNP="/home/tsigabu/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz"
^^^^^^^^^^^^^^^^^^^^^^^^^
RGSM=(
ABETF0000003 
ABETF0000005 
ABETF0000007 
ABETF0000008 
ABETF0000009 
ABETF0000012 
ABETF0000013 
ABETF0000014 
ABETF0000015 
ABETF0000016 
ABETF0000020 
ABETF0000021 
ABETF0000023 
ABETF0000026 
ABETF0000033 
ABETF0000034 
ABETF0000035 
ABETF0000036 
ABETM0000006 
ABETM0000037 
BGETF0000002 
BGETF0000003 
BGETF0000004 
BGETM0000005 
BGETF0000006 
BGETF0000007 
BGETF0000008 
BGETM0000010 
BGETM0000012 
BGETF0000013 
BGETF0000014 
BGETF0000015 
BGETF0000016 
BGETF0000018 
BGETF0000019 
BGETM0000020 
BGETM0000021 
BGETF00000022 
BGETM0000023 
BGETM0000024 
BGETM0000025 
BGETF0000026 
BGETF0000027 
ELETF0000001 
ELETF0000002 
ELETF0000003 
ELETF0000004 
ELETF0000005 
ELETF0000006 
ELETF0000007 
ELETF0000008 
ELETM0000009 
ELETF0000011 
ELETF0000012 
ELETF0000013 
ELETF0000014 
ELETF0000016 
ELETF0000017 
ELETF0000018 
ELETF0000019 
ELETM0000020 
ELETM0000021 
ELETM0000022 
MEETM0000018 
MEETF0000022 
MEETF0000023 
MEETF0000024 
MEETF0000025 
MEETF0000026 
MEETF0000027 
MEETF0000031 
MEETF0000033 
MEETF0000035 
MEETF0000036 
MEETF0000037 
MEETF0000041 
MEETF0000042 
MEETF0000045 
NUETF0000001 
NUETF0000011 
NUETM0000002 
NUETM0000046 
NUETM0000048 
NUETF0000049 
NUETF0000050 
NUETF0000051 
NUETF0000052 
NUETF0000053 
NUETF0000056 
NUETF0000057 
NUETF0000059 
NUETF0000060 
NUETF0000061 
NUETF0000062 
NUETF0000066 
NUETM0000070) 
#############################################################################
4. #Mapping_script
bwa mem -t 12 -M -R '@RG\tID:\tSM:\tPL:ILLUMINA\tPU:\tLB:\tPI:' ${REF} ${INPUT}/R1_001.fastq.gz ${INPUT}/R2_001.fastq.gz | samtools view -bS - > ${output1}${BREED}/BAMs/sample1.pe.bam 
#Sorting BAM file
for i in ${!RGSM[@]}; do
log=${output1}${BREED}/LOGs/${RGSM[i]}_SortSam.log
	java -Xmx80G -jar ${GATK}gatk-package-4.3.0.0-local.jar SortSam \
                      -I ${output1}${BREED}/BAMs/${RGSM[i]}.pe.bam \
                      -O ${output1}${BREED}/BAMs/${RGSM[i]}.sorted.bam \
                -SORT_ORDER coordinate
                     --spark-master local[${TH}] \
                     --create-output-bam-index true \
          --tmp-dir ${output1}${BREED}/TMP1/ \
        2> >(tee -a ${log})

5.#Marking Duplicates
log=${output1}${BREED}/LOGs/${RGSM[i]}-MarkDuplicatesSpark
java -Xmx80G -jar ${GATK}gatk-package-4.3.0.0-local.jar MarkDuplicatesSpark \
		-I ${output1}${BREED}/BAMs/${RGSM[i]}.sorted.bam \
		-O ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam \
		--spark-master local[${TH}] \
 		--optical-duplicate-pixel-distance 2500 \
 		--read-validation-stringency LENIENT \
		--create-output-bam-index true \
                                   --remove-all-duplicates true \
		--tmp-dir ${output1}${BREED}/TMP1/ \
	 2> >(tee -a ${log})
6. #Estimate Library Complexity
java -Xmx80G -jar ${GATK}gatk-package-4.3.0.0-local.jar EstimateLibraryComplexity \
	-I ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam \
	-O ${output1}${BREED}/metrices/${RGSM[i]}-dedup-metrices.txt \
	--TMP_DIR ${output1}${BREED}/TMP2/ \
  2> >(tee -a ${log})
7.#Base quality score recalibration (BQSR)
java -Xmx80G -jar ${GATK}/gatk-package-4.3.0.0-local.jar BaseRecalibratorSpark \
                -I ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam \
                -R ${REF} \
                --known-sites ${dbSNP} \
               -O ${output1}${BREED}/RECTABs/${RGSM[i]}.before.recal.table \
              --spark-master local[${TH}] \
                --tmp-dir ${output1}${BREED}/TMP1/ \
        2> >(tee -a ${log})
8. #Apply BQSR
 java -Xmx80G -jar ${GATK}/gatk-package-4.3.0.0-local.jar ApplyBQSRSpark \
                -I ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam  \
               -R ${REF} \
               --bqsr-recal-file ${output1}${BREED}/RECTABs/${RGSM[i]}.before.recal.table \
                -O ${output1}${BREED}/RECBAMs/${RGSM[i]}-dedup.recal.bam \
               --create-output-bam-index true \
               --spark-master local[${TH}] \
               --tmp-dir ${output1}${BREED}/TMP2/ \
                2> >(tee -a ${log})
9.#Base Recalibration
 log=${output1}${BREED}/LOGs/${RGSM[i]} -BQSR.log
        java -Xmx80G -jar ${GATK}/gatk-package-4.3.0.0-local.jar BaseRecalibratorSpark \
                -I ${output1}${BREED}/RECBAMs/${RGSM[i]}-dedup.recal.bam \
                -R ${REF} \
                --known-sites ${dbSNP} \
               -O ${output1}${BREED}/RECTABs/${RGSM[i]}.after.recal.table \
              --spark-master local[${TH}] \
              --tmp-dir ${output1}${BREED}/TMP1/ \
        2> >(tee -a ${log})
10. #Variant (SNP)call
     java -Xmx80G -jar ${GATK}gatk-package-4.3.0.0-local.jar HaplotypeCaller \
                -R ${REF} \
                -I ${output1}${BREED}/RECBAMs/${RGSM[i]}-dedup.recal.bam \
                -O ${output1}${BREED}/GVCFs/${RGSM[i]}.g.vcf.gz \
                --tmp-dir ${output1}${BREED}/TMP1/ \
                -ERC GVCF \
                --dbsnp ${dbSNP} \
                -G StandardAnnotation 
                -G AS_StandardAnnotation
11.#Combine GVCFs 
#!/bin/bash/
gatk CombineGVCFs -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" --dbsnp "/home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz" --variant ABETF0000003.g.vcf.gz --variant ABETF0000005.g.vcf.gz --variant ABETF0000007.g.vcf.gz --variant ABETF0000008.g.vcf.gz --variant ABETF0000009.g.vcf.gz --variant ABETF0000012.g.vcf.gz --variant ABETF0000013.g.vcf.gz --variant ABETF0000014.g.vcf.gz --variant ABETF0000015.g.vcf.gz -variant ABETF0000016.g.vcf.gz --variant ABETF0000020.g.vcf.gz --variant ABETF0000021.g.vcf.gz --variant ABETF0000023.g.vcf.gz --variant ABETF0000026.g.vcf.gz --variant ABETF0000033.g.vcf.gz --variant ABETF0000034.g.vcf.gz --variant ABETF0000035.g.vcf.gz --variant ABETF0000036.g.vcf.gz --variant ABETM0000006.g.vcf.gz --variant ABETM0000037.g.vcf.gz --variant BGETF0000002.g.vcf.gz --variant BGETF0000003.g.vcf.gz --variant BGETF0000004.g.vcf.gz --variant BGETM0000005.g.vcf.gz --variant BGETF0000006.g.vcf.gz --variant BGETF0000007.g.vcf.gz --variant BGETF0000008.g.vcf.gz --variant BGETM0000010.g.vcf.gz --variant BGETM0000012.g.vcf.gz --variant BGETF0000013.g.vcf.gz --variant BGETF0000014.g.vcf.gz --variant BGETF0000015.g.vcf.gz --variant BGETF0000016.g.vcf.gz --variant BGETF0000018.g.vcf.gz --variant BGETF0000019.g.vcf.gz --variant BGETM0000020.g.vcf.gz --variant BGETM0000021.g.vcf.gz --variant BGETF0000022.g.vcf.gz --variant BGETM0000023.g.vcf.gz --variant BGETM0000024.g.vcf.gz --variant BGETM0000025.g.vcf.gz --variant BGETF0000026.g.vcf.gz --variant BGETF0000027.g.vcf.gz --variant ELETF0000001.g.vcf.gz --variant ELETF0000002.g.vcf.gz --variant ELETF0000003.g.vcf.gz --variant ELETF0000004.g.vcf.gz --variant ELETF0000005.g.vcf.gz --variant ELETF0000006.g.vcf.gz --variant ELETF0000007.g.vcf.gz --variant ELETF0000008.g.vcf.gz --variant ELETM0000009.g.vcf.gz --variant ELETF0000011.g.vcf.gz --variant ELETF0000012.g.vcf.gz --variant ELETF0000013.g.vcf.gz --variant ELETF0000014.g.vcf.gz --variant ELETF0000016.g.vcf.gz --variant ELETF0000017.g.vcf.gz --variant ELETF0000018.g.vcf.gz --variant ELETF0000019.g.vcf.gz --variant ELETM0000020.g.vcf.gz --variant ELETM0000021.g.vcf.gz --variant ELETM0000022.g.vcf.gz --variant MEETM0000018.g.vcf.gz --variant MEETF0000022.g.vcf.gz --variant MEETF0000023.g.vcf.gz --variant MEETF0000024.g.vcf.gz --variant MEETF0000025.g.vcf.gz --variant MEETF0000026.g.vcf.gz --variant MEETF0000027.g.vcf.gz --variant MEETF0000031.g.vcf.gz --variant MEETF0000033.g.vcf.gz --variant MEETF0000035.g.vcf.gz --variant MEETF0000036.g.vcf.gz --variant MEETF0000037.g.vcf.gz --variant MEETF0000041.g.vcf.gz --variant MEETF0000042.g.vcf.gz --variant MEETF0000045.g.vcf.gz --variant NUETF0000001.g.vcf.gz --variant NUETF0000011.g.vcf.gz --variant NUETM0000002.g.vcf.gz --variant NUETM0000046.g.vcf.gz --variant NUETM0000048.g.vcf.gz --variant NUETF0000049.g.vcf.gz --variant NUETF0000050.g.vcf.gz --variant NUETF0000051.g.vcf.gz --variant NUETF0000052.g.vcf.gz --variant NUETF0000053.g.vcf.gz --variant NUETF0000056.g.vcf.gz --variant NUETF0000057.g.vcf.gz --variant NUETF0000059.g.vcf.gz --variant NUETF0000060.g.vcf.gz --variant NUETF0000062.g.vcf.gz --variant NUETF0000061.g.vcf.gz --variant NUETF0000066.g.vcf.gz --variant NUETM0000070.g.vcf.gz --output 96_samples.g.vcf.gz

12.#Joint genotyping
#!bash
module load gatk/4.3.0.0 
gatk GenotypeGVCFs -R  "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" --dbsnp /home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz -V /home/tgebreselassie/output1/ABERGELLE/GVCFs/Final.Final2.119.g.vcf.gz -O /home/tgebreselassie/output1/ABERGELLE/VCF/96_samples.vcf.gz
13.#VQSR (Variant Quality Score Recalibration)
VQSR for SNP
Gatk VariantRecalibrator -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" -V 96_samples.vcf.gz -resource:ensembl,known=false,training=true,truth=true,prior=15 "/home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz" -resource:dbsnp,known=true,training=false,truth=false,prior=7 "/home/tgebreselassie/119_Omega/ovis_aries_rambouillet.vcf.gz" -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -tranche 90.0 -O 96_samples.SNP.recal --tranches-file Samples.SNP.tranches --rscript-file Samples_recal.plot.R
14.#VQSR
 Apply VQSR for SNP
 gatk ApplyVQSR -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" -V 96_samples.SNP.recal --truth-sensitivity-filter-level 99.0 -mode SNP --tranches-file Samples_recal.SNP.tranches --recal-file Samples_recal.SNP.recal -O 96_samples.vcf.gz
#VQSR for INDELs
gatk VariantRecalibrator 
-R "/home/REF/Reference.fa" -V 96_samples.vcf.gz --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP --max-gaussians 4 -mode INDEL -resource:ensembl,known=false,training=true,truth=true,prior=15 "/home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz" -resource:dbsnp,known=true,training=false,truth=false,prior=7 "/home/tgebreselassie/119_Omega/ovis_aries_rambouillet.vcf.gz" -O cohort_INDELs.recal --tranches-file cohort_INDELs.tranches 
Apply VQSR for SNP
 gatk ApplyVQSR -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" -V cohort_INDELs.tranches --truth-sensitivity-filter-level 99.0 -mode INDEL --tranches-file Samples_recal.INDEL.tranches --recal-file Samples_recal.INDEL.recal -O 96_samples.INDELs.VQSR.vcf.gz

15.#Selecting only biallelic and autosomal SNPs
gatk SelectVariants \
 -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" -V Sample.SNP.VQSR.vcf.gz --select-type-to-include SNP -O 96_samples.biallelic.vcf.gz --restrict-alleles-to BIALLELIC
16.#Selecting autosomal SNPs
gatk SelectVariants -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa" -V 96_samples.biallelic.vcf.gz --intervals 1 --intervals 2 --intervals 3 --intervals 4 --intervals 5 --intervals 6 --intervals 7 --intervals 8 --intervals 9 --intervals 10 --intervals 11 --intervals 12 --intervals 13 --intervals 14 --intervals 15 --intervals 16 --intervals 17 --intervals 18 --intervals 19 --intervals 20 --intervals 21 --intervals 22 --intervals 23 --intervals 24 --intervals 25 --intervals 26 -O 96_samples.vcf.gz

17.# BAM file summary statistics plots 
# Load libraries
module load R/4.2
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# Read Excel file
df <- read_excel("C:/Users/Tsigabu/Desktop/2026_analysis/Scientific_data/BAM_stat.xlsx")

# Set strict alphabetical breed order
df$Breed <- factor(df$Breed,
                   levels = sort(unique(df$Breed)))

# Define breed colors
breed_colors <- c("Abergelle" = "green",
                  "Begait" = "red",
                  "Elle" = "cyan",
                  "Majang" = "blue",
                  "Nuer" = "yellow",
                  "Simien" = "purple")

# Function to plot a metric
plot_metric <- function(metric) {
  ggplot(df, aes(x = Breed, y = .data[[metric]], fill = Breed)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(color = "black", width = 0.2, size = 1.5) +
    scale_fill_manual(values = breed_colors) +
    theme_bw(base_size = 18) +  # larger font for clarity
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_blank()) +
    labs(y = metric, x = "Breed")
}

# Generate plots
plot_a <- plot_metric("Depth of coverage")
plot_b <- plot_metric("MeanCoverage")
plot_c <- plot_metric("MeanBaseQ")
plot_d <- plot_metric("MeanMAPQ")

# Combine into a single figure with external labels
combined_plot <- (plot_a | plot_b) / (plot_c | plot_d) +
  plot_annotation(tag_levels = "a")

# Save high-resolution outputs
ggsave("C:/Users/Tsigabu/Desktop/2026_analysis/Scientific_data/BAM_stat_boxplots.pdf",
       combined_plot, width = 14, height = 12)   

ggsave("C:/Users/Tsigabu/Desktop/2026_analysis/Scientific_data/BAM_stat_boxplots.tiff",
       combined_plot, width = 14, height = 12, dpi = 600, compression = "lzw") 

18.#Variants density plot across autosomal chromosomes (Chr1-26)
module load R/4.2
 library(vcfR)
library(CMplot)
library(Cairo)

# =========================
# 1. READ VCF
# =========================
vcf <- read.vcfR("96_samples.vcf.gz")

dat <- data.frame(
  SNP = paste0("SNP_", 1:length(getPOS(vcf))),
  Chromosome = as.character(getCHROM(vcf)),
  Position = getPOS(vcf),
  P = runif(length(getPOS(vcf)), 0, 1)   # required for CMplot framework
)

# =========================
# 2. ORDER CHROMOSOMES
# =========================
dat$Chromosome <- factor(dat$Chromosome, levels = as.character(1:26))
dat <- dat[order(dat$Chromosome, dat$Position), ]

# =========================
# 3. SUBSAMPLE (CRITICAL for 39M SNPs)
# =========================
set.seed(123)
dat_small <- dat[sample(nrow(dat), 2000000), ]  # 2M SNPs for plotting

# =========================================================
# A. LINEAR MANHATTAN (PUBLICATION MAIN FIGURE)
# =========================================================
CairoPNG("Fig1_Manhattan_linear.png", width = 8000, height = 3000, res = 300)

CMplot(
  dat_small,
  plot.type = "m",
  col = c("darkgreen", "goldenrod"),
  bin.size = 1e6,
  file.output = FALSE,
  verbose = TRUE
)

dev.off()

# =========================================================
# B. CIRCULAR MANHATTAN + DENSITY 
# =========================================================
CairoPNG("Fig2_Manhattan_circular.png", width = 7000, height = 7000, res = 400)

CMplot(
  dat_small,
  plot.type = "c",
  type = "p",

  bin.size = 1e6,
  r = 0.4,

  chr.labels = paste0("Chr", 1:26),

  chr.den.col = c("darkgreen", "yellow", "red"),

  cir.axis = TRUE,
  cir.axis.col = "black",
  cir.chr.h = 1.3,
  outward = FALSE,

  file.output = FALSE,
  verbose = TRUE
)

dev.off()

# =========================================================
# C. SNP DENSITY PLOT (GENOME VIEW)
# =========================================================
CairoPNG("Fig3_SNP_density.png", width = 7000, height = 1800, res = 300)

CMplot(
  dat,
  plot.type = "d",
  bin.size = 1e6,
  col = c("darkgreen", "yellow", "red"),
  file.output = FALSE,
  verbose = TRUE
)

dev.off()