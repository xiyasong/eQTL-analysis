# Docker run of FastQTL for JP cohort eQTL
### Following Broad Institute eQTL Pipeline (https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl)

### 1 Prepare expression, normalize the expression data
```
docker run --rm -v /Users/mstermo/Degree_Project/data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash \
    -c "/src/eqtl_prepare_expression.py --output_dir /data /data/Total_tpm.gct.gz /data/Total_gene_counts.gct.gz \
        /data/GENCODE_gencode.v26.GRCh38.ERCC.genes.gtf /data/sample_participant_lookup.txt /data/ccRCC-chr1-22-X_snp.recalibrated.99.8_output.1000G_PPs.vcf.gz.filterby.Inbreed.GQ.annotated.PASS.split.vcf.bgz.chr_list eQTL-ccRCC- \
        --tpm_threshold 0.1 --count_threshold 6 --sample_frac_threshold 0.2 --normalization_method tmm"
```

### 2 Calculate PEER factors
```
docker run --rm -v /Users/mstermo/Degree_Project/data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash \-c "Rscript /src/run_PEER.R /data/eQTL-ccRCC-.expression.bed.gz eQTL-ccRCC 15 --output_dir /data"
```

### 3 Combine covariates
```
# adding sex variant
docker run --rm -v /Users/mstermo/Degree_Project/data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash \-c "/src/combine_covariates.py /data/eQTL-ccRCC.PEER_covariates.txt eQTL-ccRCC --add_covariates /data/Sex_cov.txt -o /data"
```

### 4 Run Fastqtl
```
### norminal pass
docker run --rm -v /Users/mstermo/Degree_Project/data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash -c "/opt/fastqtl/python/run_FastQTL_threaded.py /data/df1_for_eqtl.vcf.gz /data/eQTL-ccRCC.expression.bed.gz eQTL-ccRCC-final --covariates /data/Combined_covariates_final.txt --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 8 --output_dir /data"

### permute pass
docker run --rm -v /Users/mstermo/Degree_Project/data:/data -t broadinstitute/gtex_eqtl:V8 /bin/bash -c "/opt/fastqtl/python/run_FastQTL_threaded.py /data/df1_for_eqtl.vcf.gz /data/eQTL-ccRCC.expression.bed.gz eQTL-ccRCC-final --covariates /data/Combined_covariates_final.txt --window 1e6 --chunks 100 --threads 16 --permute 1000 10000 --output_dir /data"
```

### The following files will be generated:
```
eQTL-ccRCC-final.allpairs.txt.gz
eQTL-ccRCC-final.egenes.txt.gz
```
