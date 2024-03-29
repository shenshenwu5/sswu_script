args=commandArgs(TRUE)
if (length(args) != 3) {
  print ("usage: <raw_count file> <gft file> <out.prefix> ")
  q()
}

raw_count <-  args[1]
gtf_file <-  args[2]
prefix <- args[3]
##########################################
library(GenomicFeatures)

Counts <- read.table(file = raw_count, header = TRUE, row.names = 1,sep='\t')
# Counts <- read.csv(file = "all_sample_raw_counts.csv", header = TRUE, row.names = 1)

# "../../dapseq/00_ref/Bipolaris_maydis_c5_gca_000338975.CocheC5_3.58.gtf2"
# txdb <- makeTxDbFromGFF(gtffile, format="gtf")
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
########################################
ebg <- exonsBy(txdb, by="gene")
head(Counts)

ebgList <- sum(width(reduce(ebg)))
genes <- intersect(rownames(Counts), names(ebgList))
Length <- as.vector(ebgList[genes])

FPKM <- t(t(Counts / t(Length)) * 1e9 / colSums(Counts))
head(FPKM)

TPM <- t(t(Counts / t(Length)) * 1e6 / colSums(Counts / t(Length)))
head(TPM)

CPM <- t(t(Counts) / colSums(Counts)) * 1e6
head(CPM)
###########################################
ot_tpm <- paste(prefix,'TPM.tsv',sep='_')
ot_fpkm <- paste(prefix,'FPKM.tsv',sep='_')
ot_cpm <- paste(prefix,'CPM.tsv',sep='_')
# write.csv(TPM,file=ot_tpm,quote=F,sep='\t')
# write.csv(FPKM,file=ot_fpkm,quote=F,'\t')

write.table(TPM, file=ot_tpm, append=FALSE, quote=FALSE, sep="\t")
write.table(FPKM, file=ot_fpkm, append=FALSE, quote=FALSE, sep="\t")
write.table(CPM, file=ot_cpm, append=FALSE, quote=FALSE, sep="\t")

