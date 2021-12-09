install.packages("readr")
library(readr)
####################################
#sh bulk-RNA
rm(list=ls())

#read the temple sh file
path_raw<-"/Users/mstermo/Degree_Project/data"
setwd(path_raw)
script_test<-read_file("GenomicsDBImport.sh")

#get the interval list
intervals<- read.table("intervals.list")

#set the output folder
path_out<-paste0(path_raw,"/GenomicsDBImport_sh_output")
setwd(path_out)
sh_file_namelist <- c()

times=23
for (i in 1:times){
  script_out<-gsub(pattern = "xxxx",replacement = intervals[i,],script_test)
  write_file(script_out,paste0("GenomicsDBImport-",intervals[i,],".sh"))
  #get the sh file's namelist
  sh_file_namelist <- c(sh_file_namelist,paste0("GenomicsDBImport-",intervals[i,],".sh"))
  rm(script_out)
}
#perform sh bulk
rm(list=ls())
rm(list=ls())
rm(list=ls())

library(readr)

path_raw<-"/Users/mstermo/Degree_Project/data"
setwd(path_raw)
script_test<-read_file("perform_sh_example.txt")

times=23
script_end<-NULL
for (i in 1:times){
  
  script_out<-gsub(pattern = "xxxx",replacement =sh_file_namelist[i],script_test)
  #script_end:最终合成的多个dos2unix---sbatch的命令文件
  script_end<-c(script_end,script_out)
  
}

write.table(script_end,file="Perform_sh_end.txt",sep="\t",row.names=F,col.names=F,quote=F)