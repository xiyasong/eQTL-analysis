#!/bin/bash -l
#SBATCH -A snic2020-15-69
##SBATCH -p core
#SBATCH -n 9
#SBATCH -t 8:00:00
#SBATCH -J xiya
#SBATCH --mail-user xiya.song@scilifelab.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load GATK/4.1.1.0


gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path my_database_xxxx \
       --batch-size 50 \
       --sample-name-map sample_name_map \
       --tmp-dir=$TMPDIR \
       --reader-threads 5 \
       --intervals xxxx
wait



