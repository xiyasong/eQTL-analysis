#!/bin/bash -l
#SBATCH -A snic2020-15-69
##SBATCH -p core
#SBATCH -n 9
#SBATCH -t 10:00:00
#SBATCH -J xiya
#SBATCH --mail-user xiya.song@scilifelab.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load GATK/4.1.1.0

gatk --java-options "-Xmx4g" GenotypeGVCFs \
     -R /proj/sllstore2017024/nobackup/xiya/DNA-seq/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
     -V gendb://my_database_xxxx \
     -O ccRCC-xxxx.vcf.gz \
     --tmp-dir=$TMPDIR \
     -L xxxx \
     -G StandardAnnotation -G AS_StandardAnnotation

wait
