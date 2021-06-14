install.packages("readr")
library(readr)
####################################
#sh bulk-RNA
rm(list=ls())

#read the temple sh file
path_raw<-"/Users/mstermo/Degree_Project"
setwd(path_raw)
script_test<-read_file("RNA_seq_align_and_quantification.sh")

#get the filename list
filesname <- read.table("RNA-seq-filelist.txt")

#set the output folder
path_out<-paste0(path_raw,"/RNA_seq_sh_output")
setwd(path_out)
sh_file_namelist <- c()

times=25
for (i in 1:times){
  script_out<-gsub(pattern = "xxxx",replacement = filesname[i,],script_test)
  write_file(script_out,paste0("RNA-seq-",filesname[i,],".sh"))
  #get the sh file's namelist
  sh_file_namelist <- c(sh_file_namelist,paste0("RNA-seq-",filesname[i,],".sh"))
  rm(script_out)
}
#perform sh bulk
rm(list=ls())
rm(list=ls())
rm(list=ls())

library(readr)

path_raw<-"/Users/mstermo/Degree_Project"
setwd(path_raw)
script_test<-read_file("perform_sh_example.txt")

times=25
script_end<-NULL
for (i in 1:times){
  
  script_out<-gsub(pattern = "xxxx",replacement =sh_file_namelist[i],script_test)
  #script_end:最终合成的多个dos2unix---sbatch的命令文件
  script_end<-c(script_end,script_out)
  
}

write.table(script_end,file="Perform_sh_end.txt",sep="\t",row.names=F,col.names=F,quote=F)

