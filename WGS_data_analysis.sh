Code availability

1. Raw sequence Quality check through FastQC and MultiQC
module load fastqc/0.11.9
fastqc -t 12 *.fastq.gz
2.MultiQC (summary of the fastqc results)
module load multiqc/1.8
multiqc *_fastqc.zip
Variant discovery
#!/bin/bash/
#2. define variables
REF="/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"
BREED="ABERGELLE"
output1="/home/output1/"
GATK="/home/apps/gatk/4.3.0.0/"
t="12"
TH="4"
dbSNP="/home/tsigabu/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz"
Create outputs folders
#^^^^^^^^^^^^^^^^^^^^^^^^^
mkdir -p ${output1}${BREED}/{BAMs,MDBAMs,RECTABs,RECBAMs,GVCFs,metrices,TMP{1,2},LOGs}
3. loading tools/modules
module load bwa/0.7.17
module load samtools/1.11
module load java/8
^^^^^^^^^^^^^^^^^^^^^^^^^
RGSM=(sample1 sample2) 
#############################################################################
4. Mapping_script
bwa mem -t 12 -M -R '@RG\tID:\tSM:\tPL:ILLUMINA\tPU:\tLB:\tPI:' ${REF} ${INPUT}/R1_001.fastq.gz ${INPUT}/R2_001.fastq.gz | samtools view -bS - > ${output1}${BREED}/BAMs/sample1.pe.bam 

BAM sorting
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

5.Marking Duplicates
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
6. Estimate Library Complexity
java -Xmx80G -jar ${GATK}gatk-package-4.3.0.0-local.jar EstimateLibraryComplexity \
	-I ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam \
	-O ${output1}${BREED}/metrices/${RGSM[i]}-dedup-metrices.txt \
	--TMP_DIR ${output1}${BREED}/TMP2/ \
  2> >(tee -a ${log})
7.Base quality score recalibration (BQSR)
java -Xmx80G -jar ${GATK}/gatk-package-4.3.0.0-local.jar BaseRecalibratorSpark \
                -I ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam \
                -R ${REF} \
                --known-sites ${dbSNP} \
               -O ${output1}${BREED}/RECTABs/${RGSM[i]}.before.recal.table \
              --spark-master local[${TH}] \
                --tmp-dir ${output1}${BREED}/TMP1/ \
        2> >(tee -a ${log})
8. Apply BQSR
 java -Xmx80G -jar ${GATK}/gatk-package-4.3.0.0-local.jar ApplyBQSRSpark \
                -I ${output1}${BREED}/MDBAMs/${RGSM[i]}-dedup.bam  \
               -R ${REF} \
               --bqsr-recal-file ${output1}${BREED}/RECTABs/${RGSM[i]}.before.recal.table \
                -O ${output1}${BREED}/RECBAMs/${RGSM[i]}-dedup.recal.bam \
               --create-output-bam-index true \
               --spark-master local[${TH}] \
               --tmp-dir ${output1}${BREED}/TMP2/ \
                2> >(tee -a ${log})
9.Base Recalibration
 log=${output1}${BREED}/LOGs/${RGSM[i]} -BQSR.log
        java -Xmx80G -jar ${GATK}/gatk-package-4.3.0.0-local.jar BaseRecalibratorSpark \
                -I ${output1}${BREED}/RECBAMs/${RGSM[i]}-dedup.recal.bam \
                -R ${REF} \
                --known-sites ${dbSNP} \
               -O ${output1}${BREED}/RECTABs/${RGSM[i]}.after.recal.table \
              --spark-master local[${TH}] \
              --tmp-dir ${output1}${BREED}/TMP1/ \
        2> >(tee -a ${log})
10. Variant (SNP)call
     java -Xmx80G -jar ${GATK}gatk-package-4.3.0.0-local.jar HaplotypeCaller \
                -R ${REF} \
                -I ${output1}${BREED}/RECBAMs/${RGSM[i]}-dedup.recal.bam \
                -O ${output1}${BREED}/GVCFs/${RGSM[i]}.g.vcf.gz \
                --tmp-dir ${output1}${BREED}/TMP1/ \
                -ERC GVCF \
                --dbsnp ${dbSNP} \
                -G StandardAnnotation 
                -G AS_StandardAnnotation
11.Combine GVCFs 
#!/bin/bash/
gatk CombineGVCFs -R /home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa 
--dbsnp “/home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz” --variant ABETF0000003.g.vcf.gz --variant ABETF0000005.g.vcf.gz --variant ABETF0000007.g.vcf.gz --variant ABETF0000008.g.vcf.gz --variant ABETF0000009.g.vcf.gz --variant ABETF0000012.g.vcf.gz --variant ABETF0000013.g.vcf.gz --variant ABETF0000014.g.vcf.gz --variant ABETF0000015.g.vcf.gz --variant ABETF0000016.g.vcf.gz --variant ABETF0000020.g.vcf.gz --variant ABETF0000021.g.vcf.gz --variant ABETF0000023.g.vcf.gz --variant ABETF0000026.g.vcf.gz --variant ABETF0000033.g.vcf.gz --variant ABETF0000034.g.vcf.gz --variant ABETF0000035.g.vcf.gz --variant ABETF0000036.g.vcf.gz --variant ABETM0000006.g.vcf.gz --variant ABETM0000037.g.vcf.gz --variant BGETF0000002.g.vcf.gz --variant BGETF0000003.g.vcf.gz --variant BGETF0000004.g.vcf.gz --variant BGETM0000005.g.vcf.gz --variant BGETF0000006.g.vcf.gz --variant BGETF0000007.g.vcf.gz --variant BGETF0000008.g.vcf.gz --variant BGETM0000010.g.vcf.gz --variant BGETM0000012.g.vcf.gz --variant BGETF0000013.g.vcf.gz --variant BGETF0000014.g.vcf.gz --variant BGETF0000015.g.vcf.gz --variant BGETF0000016.g.vcf.gz --variant BGETF0000018.g.vcf.gz --variant BGETF0000019.g.vcf.gz --variant BGETM0000020.g.vcf.gz --variant BGETM0000021.g.vcf.gz --variant BGETF00000022.g.vcf.gz --variant BGETM0000023.g.vcf.gz --variant BGETM0000024.g.vcf.gz --variant BGETM0000025.g.vcf.gz --variant BGETF0000026.g.vcf.gz --variant BGETF0000027.g.vcf.gz --variant ELETF0000001.g.vcf.gz --variant ELETF0000002.g.vcf.gz --variant ELETF0000003.g.vcf.gz --variant ELETF0000004.g.vcf.gz --variant ELETF0000005.g.vcf.gz --variant ELETF0000006.g.vcf.gz --variant ELETF0000007.g.vcf.gz --variant ELETF0000008.g.vcf.gz --variant ELETM0000009.g.vcf.gz --variant ELETF0000011.g.vcf.gz --variant ELETF0000012.g.vcf.gz --variant ELETF0000013.g.vcf.gz --variant ELETF0000014.g.vcf.gz --variant ELETF0000016.g.vcf.gz --variant ELETF0000017.g.vcf.gz --variant ELETF0000018.g.vcf.gz --variant ELETF0000019.g.vcf.gz --variant ELETM0000020.g.vcf.gz --variant ELETM0000021.g.vcf.gz --variant ELETM0000022.g.vcf.gz --variant MEETM0000018.g.vcf.gz --variant MEETF0000022.g.vcf.gz --variant MEETF0000023.g.vcf.gz --variant MEETF0000024.g.vcf.gz --variant MEETF0000025.g.vcf.gz --variant MEETF0000026.g.vcf.gz --variant MEETF0000027.g.vcf.gz --variant MEETF0000031.g.vcf.gz --variant MEETF0000033.g.vcf.gz --variant MEETF0000035.g.vcf.gz --variant MEETF0000036.g.vcf.gz --variant MEETF0000037.g.vcf.gz --variant MEETF0000041.g.vcf.gz --variant MEETF0000042.g.vcf.gz --variant MEETF0000045.g.vcf.gz --variant NUETF0000001.g.vcf.gz --variant NUETF0000011.g.vcf.gz --variant NUETM0000002.g.vcf.gz --variant NUETM0000046.g.vcf.gz --variant NUETM0000048.g.vcf.gz --variant NUETF0000049.g.vcf.gz --variant NUETF0000050.g.vcf.gz --variant NUETF0000051.g.vcf.gz --variant NUETF0000052.g.vcf.gz --variant NUETF0000053.g.vcf.gz --variant NUETF0000056.g.vcf.gz --variant NUETF0000057.g.vcf.gz --variant NUETF0000059.g.vcf.gz --variant NUETF0000060.g.vcf.gz --variant NUETF0000062.g.vcf.gz --variant NUETF0000061.g.vcf.gz --variant NUETF0000066.g.vcf.gz --variant NUETM0000070.g.vcf.gz --variant AM31.g.vcf.gz --variant AM32.g.vcf.gz --variant AM42.g.vcf.gz --output 96_samples.g.vcf.gz

12.Joint genotyping
#!bash
module load gatk/4.3.0.0 
gatk GenotypeGVCFs \
-R  “/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa “
--dbsnp /home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz  
-V /home/tgebreselassie/output1/ABERGELLE/GVCFs/Final.Final2.119.g.vcf.gz  
-O /home/tgebreselassie/output1/ABERGELLE/VCF/96_samples.vcf.gz

13.VQSR (Variant Quality Score Recalibration)
VQSR for SNP
Gatk VariantRecalibrator \
-R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"  \
-V 96_samples.vcf.gz \
-resource:dbsnp,known=false,training=true,truth=true,prior=15.0${knownsites} \
-resource:dbSNP,known=true,training=false,truth=false,prior=2.0 ${dbSNP}  \
-mode SNP \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.0 \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
 -an FS \
-an MQ \
-an SOR \
-an DP \
-tranche 90.0 \
-O 96_samples.SNP.recal \
--tranches-file Samples.SNP.tranches \
--rscript-file Samples_recal.plot.R\

14.VQSR
 Apply VQSR for SNP
 gatk ApplyVQSR \ 
-R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"  \
-V 96_samples.SNP.recal \
--truth-sensitivity-filter-level 99.0 \
-mode SNP \
--tranches-file Samples_recal.SNP.tranches \
--recal-file Samples_recal.SNP.recal  \
-O 96_samples.vcf.gz\

VQSR for INDELs
gatk VariantRecalibrator 
-R  "/home/REF/Reference.fa" \
-V 96_samples.vcf.gz \
--trust-all-polymorphic -tranche 100.0 \
-tranche 99.95 -tranche 99.9 \
-tranche 99.8 -tranche 99.6 \
-tranche 99.5 -tranche 99.4 \
 -tranche 99.3 \
-tranche 99.0 \
-tranche 98.0 \
-tranche 97.0 -tranche 90.0 \
-an QD -an MQRankSum \
-an ReadPosRankSum \
-an FS  \
-an MQ \
-an SOR -an DP \
--max-gaussians 4 \
-mode INDEL \
-resource:ensembl,known=false,training=true,truth=true,prior=15 "/home/tgebreselassie/V2_dbSNP/indexed/GCA_016772045.1_current_ids.vcf.gz" \
-resource:dbsnp,known=true,training=false,truth=false,prior=7 "/home/tgebreselassie/119_Omega/ovis_aries_rambouillet.vcf.gz" \
-O cohort_INDELs.recal 
--tranches-file cohort_INDELs.tranches 

Apply VQSR for SNP
 gatk ApplyVQSR \ 
-R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"  \
-V cohort_INDELs.tranches \
--truth-sensitivity-filter-level 99.0 \
 -mode INDEL \
--tranches-file Samples_recal.INDEL.tranches \
--recal-file Samples_recal.INDEL.recal  \
-O 96_samples.INDELs.VQSR.vcf.gz \

15. Selecting only biallelic and autosomal SNPs
gatk SelectVariants \
 -R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"  \
 -V Sample.SNP.VQSR.vcf.gz --select-type-to-include SNP -O 96_samples.biallelic.vcf.gz --restrict-alleles-to BIALLELIC

16. Selecting autosomal SNPs
gatk SelectVariants  
-R "/home/tgebreselassie/REF/ARS-UI_Ramb_v2.0_genomic.fa"  \
-V 96_samples.biallelic.vcf.gz --intervals 1 --intervals 2 --intervals 3 --intervals 4 --intervals 5 --intervals 6 --intervals 7 --intervals 8 --intervals 9 --intervals 10 --intervals 11 --intervals 12 --intervals 13 --intervals 14 --intervals 15 --intervals 16 --intervals 17 --intervals 18 --intervals 19 --intervals 20 --intervals 21 --intervals 22 --intervals 23 --intervals 24 --intervals 25 --intervals 26 -O 96_samples.vcf.gz

17. Variants density plot
module load R/4.2
 library(vcfR)
library(CMplot)
library(Cairo)
# Read VCF
vcf <- read.vcfR("96_samples.vcf.gz")
# Prepare CMplot input
CMsnp <- data.frame(
  SNP = paste0("SNP_", 1:length(getPOS(vcf))),
  Chromosome = as.character(getCHROM(vcf)),
  Position = getPOS(vcf)
)
# Band range: 5 to 8
band_values <- 5:8
# Loop through bands
for (b in band_values) {

  CairoPNG(
    filename = paste0("SNP_density_band_", b, ".png"),
    width = 5000,
    height = 1600,
    res = 300
  )

  CMplot(
    CMsnp,
    plot.type = "d",
    bin.size = 1000000,   # 1 Mb windows
    col = c("darkgreen", "yellow", "red"),
    band = b,
    file.output = FALSE,
    verbose = FALSE
  )

  dev.off()
}
