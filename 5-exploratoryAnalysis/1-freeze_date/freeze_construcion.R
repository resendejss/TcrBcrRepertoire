################################################################################
## script to merge TRUST4 results
## jean resende
## big data
################################################################################
# -- pool TRUST4 results for all samples -- ####################################
## -- joining airr and report for all samples
dir <- "../../4-extraction/outputTrust4/"
samples <- gsub("_report.tsv","",list.files(dir)[
  stringr::str_detect(list.files(dir),"_report.tsv")])

for (i in 1:76) {
  report <- read.delim(paste(dir,samples[i],"_report.tsv",sep = ""))
  airr <- read.delim(paste(dir,samples[i],"_airr.tsv",sep = ""))
  report$by_id <- report$cid
  airr$by_id <- gsub("_0","",airr$sequence_id)
  airr$sample_id <- rep(samples[i],nrow(report))
  freeze.all <- merge(x=report, y=airr, by="by_id")
  write.csv(freeze.all, file = paste(samples[i],".csv",sep = ""))
}

## -- building a freeze (a table with all samples)
dirFreeze <- getwd()
path <- paste0(dirFreeze,"/",samples,".csv")

tcr.bcr_trust4 <- lapply(path, read.csv) 
tcr.bcr_trust4 <- plyr::ldply(tcr.bcr_trust4)

save(tcr.bcr_trust4, file = "freeze_tcrbcr_trust4.RData")
write.csv(tcr.bcr_trust4, file = "freeze_tcrbcr_trust4.csv")
################################################################################
