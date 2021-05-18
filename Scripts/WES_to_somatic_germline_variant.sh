#!/bin/bash -l
#SBATCH -A snic2020-15-69
#SBATCH -p core
#SBATCH -n 9
#SBATCH -t 40:00:00
#SBATCH -J xiya
#SBATCH --mail-user xiya.song@scilifelab.se
#SBATCH --mail-type=ALL


#-----------------------normal.bam(xxxx) preprocessing----------------------

#----------1:transfer raw bam file into paired-end fastq files

#the raw data:analyzed by HiSeq 2000 with the 100-bp paired-end read option
#make a directory to save one sample's all results
#picard/2.23.4
module load bioinfo-tools
module load picard

mkdir $TMPDIR/xiya_tmp
mkdir $TMPDIR/xiya_tmp/raw_data_fastq

cd $TMPDIR/xiya_tmp/raw_data_fastq

java -jar -Xmx8g $PICARD_ROOT/picard.jar SamToFastq -INPUT /proj/sllstore2017024/nobackup/xiya/DNA-seq/raw_data/Original_bam/xxxx.bam \
                                                    -INCLUDE_NON_PF_READS true -VALIDATION_STRINGENCY SILENT \
                                                    -FASTQ xxxx_read1.fastq -SECOND_END_FASTQ xxxx_read2.fastq

#----------2:bwa_index, already done and generates index files

#----------3:bwa.mem: raw fastq to mapped sam using bwa/0.7.17

module load bwa
mkdir $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file
cd $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file
bwa mem -M -t 8 /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
                $TMPDIR/xiya_tmp/raw_data_fastq/xxxx_read1.fastq \
                $TMPDIR/xiya_tmp/raw_data_fastq/xxxx_read2.fastq > xxxx.sam

#----------4:Creating a BAM file and adding Read Group information

module load picard

cd $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file
java -Xmx16g -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \
                                          INPUT=xxxx.sam \
                                          OUTPUT=xxxx.bam \
                                          SORT_ORDER=coordinate RGID=xxxx-id \
                                          RGLB=xxxx-lib RGPL=ILLUMINA RGPU=xxxx-01 RGSM=xxxx
rm xxxx.sam

#--------------index bamfile--------------

cd $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file
java -Xmx16g -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=xxxx.bam 
#------------------------------------------

#----------start from here ,create a xxxx_results file and saved all following files into this foler: in the tmp folder of uppmax

#---------5:Local Realignment step 1 by older GATK

module load GATK
cd $TMPDIR/xiya_tmp 
mkdir $TMPDIR/xiya_tmp/xxxx_DNAseq_results
cd $TMPDIR/xiya_tmp/xxxx_DNAseq_results

java -Xmx16g -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -known /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    -I $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file/xxxx.bam \
    -o xxxx.intervals

#6: Local Realignment step 2

cd $TMPDIR/xiya_tmp/xxxx_DNAseq_results

java -Xmx16G -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -targetIntervals xxxx.intervals \
    -known /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    -I $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file/xxxx.bam \
    -o xxxx.realigned.bam


#7:markduplicates;
module load GATK/4.1.1.0
cd $TMPDIR/xiya_tmp/xxxx_DNAseq_results

gatk MarkDuplicates \
      -I xxxx.realigned.bam \
      -O xxxx.realigned.marked.bam \
      -M xxxx.realigned.marked.bam.metrics

rm xxxx.realigned.bam
 
#8:Base (Quality Score) Recalibration:step 1

cd $TMPDIR/xiya_tmp/xxxx_DNAseq_results

gatk BaseRecalibrator \
   -I xxxx.realigned.marked.bam \
   -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
   --known-sites /sw/data/GATK/hg38/dbsnp_138.hg38.vcf.gz \
   --known-sites /sw/data/GATK/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites /sw/data/GATK/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O xxxx.recal_data.table

#9: BQSR step 2

gatk ApplyBQSR \
   -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
   -I xxxx.realigned.marked.bam \
   --bqsr-recal-file xxxx.recal_data.table \
   -O xxxx.realigned.marked.BQSR.bam

#-------------------------------remove unneeded file--------------
cd $TMPDIR/xiya_tmp/xxxx_DNAseq_results
rm xxxx.realigned.marked.bam

#generated analysis ready bam: normal

#------------------------------------same procedure as tumor.bam-------------------------------------------

#----------1:transfer raw bam file into paired-end fastq files

#the raw data:analyzed by HiSeq 2000 with the 100-bp paired-end read option
#make a directory to save one sample's all results

#picard/2.23.4
module load bioinfo-tools
module load picard

mkdir $TMPDIR/xiya_tmp
mkdir $TMPDIR/xiya_tmp/raw_data_fastq

cd $TMPDIR/xiya_tmp/raw_data_fastq

java -jar -Xmx8g $PICARD_ROOT/picard.jar SamToFastq -INPUT /proj/sllstore2017024/nobackup/xiya/DNA-seq/raw_data/Original_bam/yyyy.bam \
                                                    -INCLUDE_NON_PF_READS true -VALIDATION_STRINGENCY SILENT \
                                                    -FASTQ yyyy_read1.fastq -SECOND_END_FASTQ yyyy_read2.fastq

#----------2:bwa_index, already done and generates index files

#----------3:bwa.mem: raw fastq to mapped sam using bwa/0.7.17

module load bwa
mkdir $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file

cd $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file

bwa mem -M -t 8 /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
                $TMPDIR/xiya_tmp/raw_data_fastq/yyyy_read1.fastq \
                $TMPDIR/xiya_tmp/raw_data_fastq/yyyy_read2.fastq > yyyy.sam

#----------4:Creating a BAM file and adding Read Group information

module load picard
cd $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file
java -Xmx16g -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups \
                                          INPUT=yyyy.sam \
                                          OUTPUT=yyyy.bam \
                                          SORT_ORDER=coordinate RGID=yyyy-id \
                                          RGLB=yyyy-lib RGPL=ILLUMINA RGPU=yyyy-01 RGSM=yyyy
rm yyyy.sam

#--------------index bamfile--------------

java -Xmx16g -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=yyyy.bam 
#------------------------------------------


#----------start from here ,create a xxxx_results file and saved all following files into this foler: in the tmp folder of uppmax

#---------5:Local Realignment step 1 by older GATK

module load GATK
cd $TMPDIR/xiya_tmp
mkdir $TMPDIR/xiya_tmp/yyyy_DNAseq_results
cd $TMPDIR/xiya_tmp/yyyy_DNAseq_results

java -Xmx16g -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -known /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    -I $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file/yyyy.bam \
    -o yyyy.intervals

#6: Local Realignment step 2

cd $TMPDIR/xiya_tmp/yyyy_DNAseq_results
java -Xmx16G -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    -targetIntervals yyyy.intervals \
    -known /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    -I $TMPDIR/xiya_tmp/bwa_mem_sam_bam_file/yyyy.bam \
    -o yyyy.realigned.bam

#7:markduplicates;

cd $TMPDIR/xiya_tmp/yyyy_DNAseq_results
module load GATK/4.1.1.0
gatk MarkDuplicates \
      -I yyyy.realigned.bam \
      -O yyyy.realigned.marked.bam \
      -M yyyy.realigned.marked.bam.metrics

rm yyyy.realigned.bam

#8:Base (Quality Score) Recalibration:step 1

cd $TMPDIR/xiya_tmp/yyyy_DNAseq_results
gatk BaseRecalibrator \
   -I yyyy.realigned.marked.bam \
   -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
   --known-sites /sw/data/GATK/hg38/dbsnp_138.hg38.vcf.gz \
   --known-sites /sw/data/GATK/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites /sw/data/GATK/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O yyyy.recal_data.table

#9: BQSR step 2

gatk ApplyBQSR \
   -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
   -I yyyy.realigned.marked.bam \
   --bqsr-recal-file yyyy.recal_data.table \
   -O yyyy.realigned.marked.BQSR.bam

#-------------------------------remove unneeded file--------------
cd $TMPDIR/xiya_tmp/yyyy_DNAseq_results
rm yyyy.realigned.marked.bam

#generated analysis ready bam: tumor

#----------------------------------Mutect2 calling
module load GATK/4.1.1.0
cd /proj/sllstore2017024/nobackup/xiya/DNA-seq/Mutect2-results
mkdir yyyy-xxxx
cd /proj/sllstore2017024/nobackup/xiya/DNA-seq/Mutect2-results/yyyy-xxxx

gatk  --java-options "-Xmx32G" Mutect2 \
   -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
   -I $TMPDIR/xiya_tmp/yyyy_DNAseq_results/yyyy.realigned.marked.BQSR.bam \
   -tumor yyyy \
   -I $TMPDIR/xiya_tmp/xxxx_DNAseq_results/xxxx.realigned.marked.BQSR.bam \
   -normal xxxx \
   --germline-resource /sw/data/GATK/Mutect2/af-only-gnomad.hg38.vcf.gz \
   --panel-of-normals /proj/sllstore2017024/nobackup/xiya/DNA-seq/somatic-hg38_1000g_pon.hg38.vcf.gz \
   --f1r2-tar-gz yyyy-xxxx.somatic.flr2.tar.gz \
   -O yyyy-xxxx.somatic.vcf.gz

#---------------------------------------------------------------------------------
#Filter somatic variation
#1: learn read orientation
module load GATK/4.1.1.0
cd /proj/sllstore2017024/nobackup/xiya/DNA-seq/Mutect2-results/yyyy-xxxx

gatk LearnReadOrientationModel -I yyyy-xxxx.somatic.flr2.tar.gz -O yyyy-read-orientation-model.tar.gz

#2: Run GetPileupSummaries to summarize read support for a set number of known variant sites.

gatk GetPileupSummaries \
    -I $TMPDIR/xiya_tmp/yyyy_DNAseq_results/yyyy.realigned.marked.BQSR.bam \
    -V /sw/data/GATK/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
    -L /sw/data/GATK/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
    -O yyyy-getpileupsummaries.table

#3: Estimate contamination with CalculateContamination.

gatk CalculateContamination \
        -I yyyy-getpileupsummaries.table \
        -tumor-segmentation yyyy-segments.table \
        -O yyyy-calculatecontamination.table

#4:Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:

gatk FilterMutectCalls -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
                       -V yyyy-xxxx.somatic.vcf.gz \
                       --tumor-segmentation yyyy-segments.table \
                       --contamination-table yyyy-calculatecontamination.table \
                       --ob-priors yyyy-read-orientation-model.tar.gz \
                       -O yyyy-xxxx.somatic.filtered.vcf.gz

#---------------------------------HaplotypeCaller calling---- (with allel specific model)

module load GATK/4.1.1.0
cd /proj/sllstore2017024/nobackup/xiya/DNA-seq/HaplotypeCaller-results
mkdir yyyy
cd /proj/sllstore2017024/nobackup/xiya/DNA-seq/HaplotypeCaller-results/yyyy

gatk --java-options "-Xmx8g" HaplotypeCaller \
  -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
  -I $TMPDIR/xiya_tmp/yyyy_DNAseq_results/yyyy.realigned.marked.BQSR.bam \
  -O yyyy.realigned.marked.BQSR.bam.g.vcf.gz \
  -ERC GVCF \
  -G StandardAnnotation \
  -G AS_StandardAnnotation \
  -G StandardHCAnnotation

wait
