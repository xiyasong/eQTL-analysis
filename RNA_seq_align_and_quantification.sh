#!/bin/bash -l
#SBATCH -A snic2020-15-69
#SBATCH -p core
#SBATCH -n 9
#SBATCH -t 30:00:00
#SBATCH -J xiya
#SBATCH --mail-user xiya.song@scilifelab.se
#SBATCH --mail-type=ALL

#RNA-seq analysis following by GTEx

#-----------1:# BAM to FASTQ conversion: use bamtofastq_new.sh

module load bioinfo-tools
module load picard

#----------(make a new folder to save xxxx.bam all the results)

cd $TMPDIR
mkdir $TMPDIR/xiya_RNA
mkdir $TMPDIR/xiya_RNA/xxxx_results
cd $TMPDIR/xiya_RNA/xxxx_results

java -jar -Xmx8g $PICARD_ROOT/picard.jar SamToFastq \
          -INPUT /proj/snic2020-16-69/nobackup/xiya/raw_bam_data/xxxx.bam \
          -INCLUDE_NON_PF_READS true \
          -INCLUDE_NON_PRIMARY_ALIGNMENTS false \
          -VALIDATION_STRINGENCY SILENT \
          -FASTQ xxxx_r1.fastq \
          -SECOND_END_FASTQ xxxx_r2.fastq

#----------2: # STAR alignment

module load bioinfo-tools
module load star/2.5.3a
cd $TMPDIR/xiya_RNA/xxxx_results

STAR --runMode alignReads --runThreadN 4 \
     --genomeDir /proj/snic2020-16-69/nobackup/xiya/star_index_oh100 \
     --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout \
     --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 \
     --readFilesIn $TMPDIR/xiya_RNA/xxxx_results/xxxx_r1.fastq \
                   $TMPDIR/xiya_RNA/xxxx_results/xxxx_r2.fastq \
     --outFileNamePrefix $TMPDIR/xiya_RNA/xxxx_results/xxxx. \
     --outSAMstrandField intronMotif --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMtype BAM Unsorted --outSAMunmapped Within \
     --genomeLoad NoSharedMemory --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimMainSegmentMultNmax 1 \
     --outSAMattributes NH HI AS nM NM ch --outSAMattrRGline ID:rg1 SM:sm1

#----------3: sort by samtools, by coordinate

module load bioinfo-tools
module load samtools

cd $TMPDIR/xiya_RNA/xxxx_results

#samtools -T: tmpprefix -o: outbam 
samtools sort -T xxxx.Aligned.sorted. \
              -o xxxx.Aligned.sorted.out.bam \
                 xxxx.Aligned.out.bam

#----------4: mark duplicates and also index the bam file

module load bioinfo-tools
module load samtools
module load picard

cd $TMPDIR/xiya_RNA/xxxx_results

java -jar $PICARD_ROOT/picard.jar MarkDuplicates \
          -I $TMPDIR/xiya_RNA/xxxx_results/xxxx.Aligned.sorted.out.bam \
          -O xxxx.Aligned.sorted.out.md.bam \
          -PROGRAM_RECORD_ID null -MAX_RECORDS_IN_RAM 500000 \
          -SORTING_COLLECTION_SIZE_RATIO 0.25 -TMP_DIR $TMPDIR/xiya_RNA/xxxx_results \
          -M xxxx.Aligned.sorted.out.md.marked_dup_metrics.txt \
          -ASSUME_SORT_ORDER coordinate -TAGGING_POLICY DontTag -OPTICAL_DUPLICATE_PIXEL_DISTANCE 100

#samtools index [-bc] [-m INT] aln.bam|aln.cram [out.index]

samtools index xxxx.Aligned.sorted.out.md.bam \
               xxxx.Aligned.sorted.out.md.bam.bai


#----------5: run the rnaseqc to get the gene-level quantification results and saved 

cd /proj/snic2020-16-69/nobackup/xiya

mkdir xxxx_rnaseqc_results
#rnaseqc [gtf] [bam] [output] {OPTIONS}

rnaseqc /proj/snic2020-16-69/nobackup/xiya/GENCODE_gencode.v26.GRCh38.ERCC.genes.gtf \
        $TMPDIR/xiya_RNA/xxxx_results/xxxx.Aligned.sorted.out.md.bam \
        /proj/snic2020-16-69/nobackup/xiya/xxxx_rnaseqc_results \
        -s xxxx \
        -vv

wait
