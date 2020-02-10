library(fastqcr)
library(magrittr)
args <- commandArgs(T)
#datapath <- paste0("/cluster/sharedata/langzhaobo/pscngs/",project,"/")
datapath <- argv[1]
part <- argv[2]
output <- paste0(datapath,"/",part,"_fastqc_stat.txt")
qc_raw_path <- paste0(datapath,"/qc_raw")
qc_clean_path <- paste0(datapath,"/qc_clean")

qc_raw <- qc_aggregate(qc_raw_path) %>% qc_stats()
qc_clean <- qc_aggregate(qc_clean_path) %>% qc_stats()
qc <- cbind(qc_raw, qc_clean)
qc <- qc[,c(1,4,9)]
colnames(qc) <- c('sample','tot.seq','tot.clean.seq')
qc$percentage <- round(as.numeric(qc$tot.clean.seq) / as.numeric(qc$tot.seq) * 100, 1)
  
write.table(qc,output,sep='\t',row.names = F, quote=F)
