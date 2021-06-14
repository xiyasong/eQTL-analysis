#---------create sample_participant_lookup file-----------
path_raw<-"/Users/mstermo/Degree_Project"
setwd(path_raw)
#get the filename list
sample_id <- read.table("RNA-seq-filelist-total.txt")
participant_id <- read.table("participant_id.txt")

sample_participant_lookup <- data.frame(sample_id,participant_id,row.names = NULL)
#colnames(sample_participant_lookup) <- c("sample_id","participant_id")
write.table(sample_participant_lookup, file = "sample_participant_lookup.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#---------------------------------------------------------

#create sample-name-map

path_raw<-"/Users/mstermo/Degree_Project"
setwd(path_raw)
#get the filename list
sample_id <- read.table("participant_id.txt")
path_to_sample <- read.table("path_to_sample.txt")

sample_name_map <- data.frame(sample_id,path_to_sample,row.names = NULL)
#colnames(sample_participant_lookup) <- c("sample_id","participant_id")
write.table(sample_name_map, file = "sample_name_map", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

