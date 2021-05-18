#read genotype data and expression data 
expr<- read.csv2('eQTL-ccRCC.expression.bed.gz',header = TRUE,sep = '\t')
expr <- expr[,c(-1,-2,-3)]
expr[,-1]=as.data.frame(lapply(expr[,-1],as.numeric))
snps <- read.csv2('df1_for_eqtl.simple.vcf',header = TRUE,sep = '\t')
head(snps)
snps <- snps[,c(-1,-2,-4,-5,-6,-7,-8,-9)]
snps[snps == '0/0'] <- '0'
snps[snps == '0/1'] <- '1'
snps[snps == '1/1'] <- '2'
snps[snps == './.'] <- '0'
snps[,-1]=as.data.frame(lapply(snps[,-1],as.numeric))

head(expr)[1:12]
class(expr[2,3])
paste("Mean value of expression for gene ",expr$gene_id[1:10]," is ", rowMeans(expr[1:10, -1]))
paste("Mean value for SNP ",snps$snipid[1:10]," is ", rowMeans(snps[1:10, -1]))

#Choose interested genotype and gene,a simple plot showed expression variantion.
genotype = 'chr5_96895296_G_T'
genes = 'ENSG00000164308.16'
par(mfrow=c(1,length(genotype)))
plot(jitter(as.numeric(snps[1,-1]), factor = 0.5), as.numeric(expr[1,-1]),
     xlab = genotype[1], ylab = genes[1], col = "steelblue",
     main = paste(genes[1], "vs", genotype[1]), xlim= c(-0.5,2.5), xaxt="n")
axis(1, at =c (0,1,2), labels = c("0", "1", "2"))

#Transpose data frames so that we have variables, i.e. SNPs and expression levels as columns and samples as rows.
expr_trans = data.frame(t(expr[, -1]))
colnames(expr_trans)=t(expr[, 1])
expr_trans = tibble::rownames_to_column(expr_trans, "sample")
head(expr_trans)[1:10]

snps_trans = data.frame(t(snps[-1]))
colnames(snps_trans)=t(snps[, 1])
snps_trans = tibble::rownames_to_column(snps_trans, "sample")
head(snps_trans)[1:10]

#Another convenient way to display gene expression values by genotype is as box plots. 
#These provide a good, nonparametric, indication of the distributions. 
#To convey a sense of the frequency of each genotype in the sample it is useful to also add points for each individual to the plot. 
#Now use R’s ggplot2 library to generate visualization.
library(ggplot2)

#Reshape dataframes for use with ggplot2
dev.off()
#snp-1 vs gene-1, snp-2 vs gene-2
#This enable to draw boxplots the associations between first 10 genotypes and 10 genes.
snps_long = tidyr::gather(snps_trans[, 1:10], snp, genotype, -sample)
expr_long = tidyr::gather(expr_trans[, 1:10], gene, expression, -sample)
data_long <- cbind(snps_long, expr_long["expression"])
data_long$genotype <- as.factor(data_long$genotype)
head(data_long)

ggplot(data_long, aes(genotype, expression)) +
  geom_jitter(colour = "darkorange",alpha = 0.3, width = 0.02) +
  geom_boxplot(alpha = 0.5, fill = "steelblue") +
  facet_wrap(~snp) 

# This can draw snp-gene pair specific box plot __________________

snps_choose = 'chr11_74964851_T_G'
genes_choose = 'ENSG00000166435.15'
par(mfrow=c(1,length(snps_choose)))

dev.off()

for (index in seq(length(snps_choose))){
  genotype = snps_trans[[snps_choose[index]]]
  expression = expr_trans[[genes_choose[index]]]
  plot(jitter(genotype, factor = 0.4), expression,
       main=paste(snps_choose[index], "vs", genes_choose[index]), xlim= c(-0.5,2.5),
       xlab = "genotype", xaxt="n", col ="steelblue")
  axis(1, at=c(0,1,2), labels = c("0", "1", "2"))
}
#boxplot:

genoLong = tidyr::gather(snps_trans, snp, genotype,snps_choose)
#head(genoLong)
exprLong = tidyr::gather(expr_trans, gene, expression,genes_choose)
#head(exprLong)
dataLong = cbind(genoLong[,c("snp", "genotype")], exprLong[,c("gene", "expression")])
dataLong$comparison = paste(dataLong$snp, "vs", dataLong$gene)
dataLong$genotype = factor(dataLong$genotype)

ggplot(dataLong, aes(genotype, expression)) +
  geom_jitter(col="darkorange", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="steelblue") +
  facet_wrap(~comparison) 