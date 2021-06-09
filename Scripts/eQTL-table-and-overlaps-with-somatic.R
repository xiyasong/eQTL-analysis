rm(list=ls())
#install.packages("dplyr")
library(dplyr)
path_raw<-"/Users/mstermo/Degree_Project/data"

# Calculate the eGenes
setwd(path_raw)
somatic_mutation_gene <- read.table('only_somatic_gene.list.txt',header = FALSE,row.names = NULL)
#check the somatic mutation's frequency:
somatic_mutation_gene_vector <- as.vector(somatic_mutation_gene[['V1']])
somatic_mutation_gene_vector
somatic_mutation_gene_vector_more_5 <- somatic_mutation_gene_vector[table(somatic_mutation_gene_vector)>5]
#plot the frequencies
somatic_mutation_frequencies <- as.data.frame(table(somatic_mutation_gene_vector))
hist(somatic_mutation_frequencies$Freq, breaks = 50,xlim = c(0,15),plot = TRUE)
#filter the mutations that lower than 5
somatic_mutation_frequencies_5 <- filter(somatic_mutation_frequencies, Freq >=5)

#unique genes
somatic_mutation_unique_gene <- read.table('only_somatic_gene.unique.list.txt',header = FALSE,row.names = NULL)

#my japan cohort eqtl
setwd(path_raw) 

#Japan_eqtl= read.table('eQTL-ccRCC.genes.txt.annotated.txt',header = TRUE)
Japan_eqtl_final = read.table('eQTL-ccRCC-final.genes.annotated.txt',header = TRUE)
Japan_eqtl_sigpairs = read.table('eQTL-ccRCC-final.signifpairs.txt',header = TRUE)
#Japan_eqtl_0.05 <- Japan_eqtl %>% filter(qval <0.05)
Japan_eqtl_final_0.05 <- Japan_eqtl_final %>% filter(qval <0.05)
?write.csv
write.table(Japan_eqtl_final_0.05,file = 'Japan_eqtl_final_0.05',sep = '\t',quote = FALSE,row.names = FALSE)
#Japan_eqtl_0.2 <- Japan_eqtl  %>% filter(qval <0.2)

#my japan cohort egenes
#Japan_gene_list <- Japan_eqtl[,'gene_name']
#Japan_gene_0.05_list <- Japan_eqtl_0.05[,"gene_name"]
Japan_final_gene_0.05_list <- Japan_eqtl_final_0.05[,"gene_name"]

#Japan_gene_0.2_list <- Japan_eqtl_0.2[,"gene_name"]

#GTEx eqtl
GTEx_eqtl= read.table('Kidney_Cortex.v8.egenes.txt',header = TRUE, fill = TRUE)
GTEx_eqtl_0.05 <- GTEx_eqtl  %>% filter(qval <0.05)
#GTEx_eqtl_0.2 <- GTEx_eqtl  %>% filter(qval <0.2)

#GTEx egenes
GTEx_gene_list <- GTEx_eqtl[,'gene_name']

GTEx_gene_0.05_list <- GTEx_eqtl_0.05[,"gene_name"]
#GTEx_gene_0.2_list <- GTEx_eqtl_0.2[,"gene_name"]

#statistics overlap 
#Overlap_0.05 <- intersect(GTEx_gene_0.05_list,Japan_gene_0.05_list)

Overlap_0.05_final <- intersect(GTEx_gene_0.05_list,Japan_final_gene_0.05_list)

#Overlap_0.2 <- intersect(GTEx_gene_0.2_list,Japan_gene_0.2_list)

#Overlap_all <- intersect(GTEx_gene_list,Japan_gene_list)

#New_gene_all <-setdiff(Japan_gene_list,GTEx_gene_list)
#New_gene_0.05 <- setdiff(Japan_gene_0.05_list,GTEx_gene_0.05_list)

New_gene_final_0.05 <- setdiff(Japan_final_gene_0.05_list,GTEx_gene_0.05_list)

#New_gene_0.2 <- setdiff(Japan_gene_0.2_list,GTEx_gene_0.2_list)

#choose New_gene_0.05, see the overlap with somatic mutations?
class(New_gene_0.05)
class(somatic_mutation_unique_gene)

#convert somatic_mutation_unique_gene to vector | also later, compare the mutations which >5

somatic_mutation_unique_gene_vector <- as.vector(somatic_mutation_unique_gene[['V1']])
somatic_mutation_frequencies_5_vector <- as.vector(somatic_mutation_frequencies_5[['somatic_mutation_gene_vector']])
# check overlap , =995 genes
#Overlap_somatic_egene <- intersect(New_gene_0.05,somatic_mutation_unique_gene_vector)
Overlap_somatic_egene_final <- intersect(New_gene_final_0.05,somatic_mutation_unique_gene_vector)
write(Overlap_somatic_egene_final)
Overlap_somatic_egene_5 <- intersect(New_gene_final_0.05,somatic_mutation_frequencies_5_vector)

#Overlap_somatic_egene_5_mis <- intersect(New_gene_mis_0.05,somatic_mutation_frequencies_5_vector)

#Overlap_test <- intersect(Overlap_somatic_egene_5, Overlap_somatic_egene_5_mis)

Overlap_somatic_egene <- as.data.frame(Overlap_somatic_egene)
Overlap_somatic_egene_5 <- as.data.frame(Overlap_somatic_egene_5)
#write into txt
write.table(Overlap_somatic_egene_5, file = 'Overlap_93_genes.txt',quote = FALSE,row.names = FALSE)
write.table(Overlap_somatic_egene_final, file = 'Overlap_genes.txt',quote = FALSE,row.names = FALSE)
somatic_mutation_unique_gene


#Draw the plot for tss_distance

TSS_distance <- Japan_eqtl_final_0.05[,"tss_distance"]
TSS_distance
TSS_distance_0.05 <- d_0.05[,8]

f <- function( x ) x / 1000

TSS_distance_kb <- f(TSS_distance)
#TSS_distance_kb_0.05 <- f(TSS_distance_0.05)
hist(TSS_distance_kb,main = "Distance distribution between the most significant cis-eQTL and the gene TSS in Kb", xlab = "Distance(Kb)",col = NULL,breaks  = 30)
#hist(TSS_distance_kb_0.05,main = "distribution of the distance (TSS) in Kb", xlab = "Distance(Kb)",col = NULL,breaks  =30)
?hist
?write.csv
install.packages("qqman")
library(qqman)
?qqman
manhattan
str(gwasResults)
head(gwasResults)

#PCA plot

pca <- read.table("/Users/mstermo/Degree_Project/data/WES.eigenvec", header = T)
eigval <- read.table("/Users/mstermo/Degree_Project/data/WES.eigenval", header = F)
pcs <- paste0("PC", 1:nrow(eigval))
#eigval[nrow(eigval),1] <- 0
percentage <- eigval$V1/sum(eigval$V1)*100
eigval_df <- as.data.frame(cbind(pcs, eigval[,1], percentage), stringsAsFactors = F)
names(eigval_df) <- c("PCs", "variance", "proportion")
eigval_df$variance <- as.numeric(eigval_df$variance)
eigval_df$proportion <- as.numeric(eigval_df$proportion)

pc1_proportion <- paste0(round(eigval_df[1,3],2),"%")
pc2_proportion <- paste0(round(eigval_df[2,3],2),"%")

p_pca <- ggplot(pca,aes(PC1,PC2))+
  geom_point(size = 2.5)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(colour = "black", size=12, family = "Times New Roman"),
        axis.title = element_text(color="black",size = 15, family = "Times New Roman"),
        #aspect.ratio = 1,
        legend.text = element_text(colour = "black", size=15, family = "Times New Roman"),
        legend.position = c(0.15,0.15)
  )+
  labs(x=paste0("PC1(",pc1_proportion,")"),
       y=paste0("PC2(",pc2_proportion,")"))+
  coord_fixed()+
  scale_x_continuous(breaks=seq(-5,5,0.1))+
  scale_y_continuous(breaks=seq(-5,5,0.1))
  
p_pca

