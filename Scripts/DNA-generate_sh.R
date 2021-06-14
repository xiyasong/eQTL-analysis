install.packages("readr")
library(readr)
####################################
#sh bulk-DNA
rm(list=ls())

#read the temple sh file
path_raw<-"/Users/mstermo/Degree_Project"
setwd(path_raw)
script_test<-read_file("WGS_to_somatic_germline_variant.sh")

#get the filename list
filesname.normal.xxxx <- read.table("DNA-seq-normal-filelist_1.txt")
filesname.tumor.yyyy <- read.table("DNA-seq-tumor-filelist_1.txt")
#set the output folder
path_out<-paste0(path_raw,"/DNA_seq_sh_output")
setwd(path_out)
sh_file_namelist <- c()
times=28
for (i in 1:times){
  script_out<-gsub(pattern = "xxxx",replacement = filesname.normal.xxxx[i,],script_test)
  script_out_final<-gsub(pattern = "yyyy",replacement = filesname.tumor.yyyy[i,],script_out)
  write_file(script_out_final,paste0("DNA-seq-",filesname.tumor.yyyy[i,],".sh"))
  #get the sh file's namelist
  sh_file_namelist <- c(sh_file_namelist,paste0("DNA-seq-",filesname.tumor.yyyy[i,],".sh"))
  rm(script_out)
  rm(script_out_final)
}
##########################################################################
#perform sh bulk

library(readr)

path_raw<-"/Users/mstermo/Degree_Project"
setwd(path_raw)
script_test<-read_file("perform_sh_example.txt")

times=28
script_end<-NULL
for (i in 1:times){
  
  script_out<-gsub(pattern = "xxxx",replacement =sh_file_namelist[i],script_test)
  #script_end:最终合成的多个dos2unix---sbatch的命令文件
  script_end<-c(script_end,script_out)
  
}

write.table(script_end,file="Perform_sh_end.txt",sep="\t",row.names=F,col.names=F,quote=F)
