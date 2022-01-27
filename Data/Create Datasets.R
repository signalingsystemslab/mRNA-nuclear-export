setwd('~/Documents/Disk1/Documents/Ribo_Modeling_Reprocessed/Paper Final Run/')

# create directories
dir.create('Data', showWarnings = F)
dir.create('Data/raw_counts', showWarnings = F)
dir.create('Data/exons_rpkms', showWarnings = F)

# load caRNA data
caRNA_Naive_b1 <- read.table(file = 'merged_counts/chromatin.Naive.b1/D40000/all/all_counts.txt', header=T)
caRNA_Naive_b2 <- read.table(file = 'merged_counts/chromatin.Naive.b2/D40000/all/all_counts.txt', header=T)
caRNA_Naive_b3 <- read.table(file = 'merged_counts/chromatin.Naive.b3/D40000/all/all_counts.txt', header=T)
caRNA_Naive_all <- data.frame(caRNA_Naive_b1[, 1:6],
                              caRNA_Naive_b1[, -c(1:6)],
                              caRNA_Naive_b2[, -c(1:6)],
                              caRNA_Naive_b3[, -c(1:6)])
rm(caRNA_Naive_b1, caRNA_Naive_b2, caRNA_Naive_b3)

colnames(caRNA_Naive_all) <- gsub("merged_bams_mm10_vM14.chromatin.Naive.b..","", colnames(caRNA_Naive_all))
colnames(caRNA_Naive_all) <- gsub(".bam","", colnames(caRNA_Naive_all))


caRNA_Naive_b1 <- read.table(file = 'merged_counts/chromatin.Naive.b1/D40000/exons/exon_counts.txt', header=T)
caRNA_Naive_b2 <- read.table(file = 'merged_counts/chromatin.Naive.b2/D40000/exons/exon_counts.txt', header=T)
caRNA_Naive_b3 <- read.table(file = 'merged_counts/chromatin.Naive.b3/D40000/exons/exon_counts.txt', header=T)
caRNA_Naive_exons <- data.frame(caRNA_Naive_b1[,1:6],
                                caRNA_Naive_b1[, -c(1:6)],
                                caRNA_Naive_b2[, -c(1:6)],
                                caRNA_Naive_b3[, -c(1:6)])
rm(caRNA_Naive_b1, caRNA_Naive_b2, caRNA_Naive_b3)

colnames(caRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.chromatin.Naive.b..","", colnames(caRNA_Naive_exons))
colnames(caRNA_Naive_exons) <- gsub(".bam","", colnames(caRNA_Naive_exons))

save(caRNA_Naive_all, caRNA_Naive_exons, file='Data/raw_counts/caRNA_Naive.Rdata')

caRNA_LPA_b1 <- read.table(file = 'merged_counts/chromatin.LPA.b1/D40000/all/all_counts.txt', header=T)
caRNA_LPA_b2 <- read.table(file = 'merged_counts/chromatin.LPA.b2/D40000/all/all_counts.txt', header=T)
caRNA_LPA_b3 <- read.table(file = 'merged_counts/chromatin.LPA.b3/D40000/all/all_counts.txt', header=T)
caRNA_LPA_all <- data.frame(caRNA_LPA_b1[,1:6],
                            caRNA_LPA_b1[, -c(1:6)],
                            caRNA_LPA_b2[, -c(1:6)],
                            caRNA_LPA_b3[, -c(1:6)])
rm(caRNA_LPA_b1, caRNA_LPA_b2, caRNA_LPA_b3)

colnames(caRNA_LPA_all) <- gsub("merged_bams_mm10_vM14.chromatin.LPA.b..","", colnames(caRNA_LPA_all))
colnames(caRNA_LPA_all) <- gsub(".bam","", colnames(caRNA_LPA_all))


caRNA_LPA_b1 <- read.table(file = 'merged_counts/chromatin.LPA.b1/D40000/exons/exon_counts.txt', header=T)
caRNA_LPA_b2 <- read.table(file = 'merged_counts/chromatin.LPA.b2/D40000/exons/exon_counts.txt', header=T)
caRNA_LPA_b3 <- read.table(file = 'merged_counts/chromatin.LPA.b3/D40000/exons/exon_counts.txt', header=T)
caRNA_LPA_exons <- data.frame(caRNA_LPA_b1[,1:6],
                              caRNA_LPA_b1[, -c(1:6)],
                              caRNA_LPA_b2[, -c(1:6)],
                              caRNA_LPA_b3[, -c(1:6)])
rm(caRNA_LPA_b1, caRNA_LPA_b2, caRNA_LPA_b3)

colnames(caRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.chromatin.LPA.b..","", colnames(caRNA_LPA_exons))
colnames(caRNA_LPA_exons) <- gsub(".bam","", colnames(caRNA_LPA_exons))

save(caRNA_LPA_all, caRNA_LPA_exons, file='Data/raw_counts/caRNA_LPA.Rdata')

## load nucleoplasmic data

npRNA_Naive_b1 <- read.table(file = 'merged_counts/Nuc.Naive.b1/D40000/all/all_counts.txt', header=T)
npRNA_Naive_b3 <- read.table(file = 'merged_counts/Nuc.Naive.b3/D40000/all/all_counts.txt', header=T)
npRNA_Naive_all <- data.frame(npRNA_Naive_b1[,1:6],
                              npRNA_Naive_b1[, -c(1:6)],
                              npRNA_Naive_b3[, -c(1:6)])
rm(npRNA_Naive_b1, npRNA_Naive_b3)

colnames(npRNA_Naive_all) <- gsub("merged_bams_mm10_vM14.Nuc.Naive.b..","", colnames(npRNA_Naive_all))
colnames(npRNA_Naive_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.Naive.b..","", colnames(npRNA_Naive_all))
colnames(npRNA_Naive_all) <- gsub(".bam","", colnames(npRNA_Naive_all))


npRNA_Naive_b1 <- read.table(file = 'merged_counts/Nuc.Naive.b1/D40000/exons/exon_counts.txt', header=T)
npRNA_Naive_b3 <- read.table(file = 'merged_counts/Nuc.Naive.b3/D40000/exons/exon_counts.txt', header=T)
npRNA_Naive_exons <- data.frame(npRNA_Naive_b1[,1:6],
                                npRNA_Naive_b1[, -c(1:6)],
                                npRNA_Naive_b3[, -c(1:6)])
rm(npRNA_Naive_b1, npRNA_Naive_b3)

colnames(npRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.Nuc.Naive.b..","", colnames(npRNA_Naive_exons))
colnames(npRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14_2.Nuc.Naive.b..","", colnames(npRNA_Naive_exons))
colnames(npRNA_Naive_exons) <- gsub(".bam","", colnames(npRNA_Naive_exons))

save(npRNA_Naive_all, npRNA_Naive_exons, file='Data/raw_counts/npRNA_Naive.Rdata')

npRNA_LPA_b1 <- read.table(file = 'merged_counts/Nuc.LPA.b1/D40000/all/all_counts.txt', header=T)
npRNA_LPA_b3 <- read.table(file = 'merged_counts/Nuc.LPA.b3/D40000/all/all_counts.txt', header=T)
npRNA_LPA_all <- data.frame(npRNA_LPA_b1[,1:6],
                            npRNA_LPA_b1[, -c(1:6)],
                            npRNA_LPA_b3[, -c(1:6)])
rm(npRNA_LPA_b1, npRNA_LPA_b3)

colnames(npRNA_LPA_all) <- gsub("merged_bams_mm10_vM14.Nuc.LPA.b..","", colnames(npRNA_LPA_all))
colnames(npRNA_LPA_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.LPA.b..","", colnames(npRNA_LPA_all))
colnames(npRNA_LPA_all) <- gsub(".bam","", colnames(npRNA_LPA_all))


npRNA_LPA_b1 <- read.table(file = 'merged_counts/Nuc.LPA.b1/D40000/exons/exon_counts.txt', header=T)
npRNA_LPA_b3 <- read.table(file = 'merged_counts/Nuc.LPA.b3/D40000/exons/exon_counts.txt', header=T)
npRNA_LPA_exons <- data.frame(npRNA_LPA_b1[,1:6],
                              npRNA_LPA_b1[, -c(1:6)],
                              npRNA_LPA_b3[, -c(1:6)])
rm(npRNA_LPA_b1, npRNA_LPA_b3)

colnames(npRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.Nuc.LPA.b..","", colnames(npRNA_LPA_exons))
colnames(npRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14_2.Nuc.LPA.b..","", colnames(npRNA_LPA_exons))
colnames(npRNA_LPA_exons) <- gsub(".bam","", colnames(npRNA_LPA_exons))

save(npRNA_LPA_all, npRNA_LPA_exons, file='Data/raw_counts/npRNA_LPA.Rdata')

## load cytoplasmic data

cytoRNA_Naive_b1 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b1/D40000/all/all_counts.txt', header=T)
cytoRNA_Naive_b2 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b2/D40000/all/all_counts.txt', header=T)
cytoRNA_Naive_b3 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b3/D40000/all/all_counts.txt', header=T)
cytoRNA_Naive_all <- data.frame(cytoRNA_Naive_b1[,1:6],
                              cytoRNA_Naive_b1[, -c(1:6)],
                              cytoRNA_Naive_b2[, -c(1:6)],
                              cytoRNA_Naive_b3[, -c(1:6)])
rm(cytoRNA_Naive_b1, cytoRNA_Naive_b2, cytoRNA_Naive_b3)

colnames(cytoRNA_Naive_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.Naive.b..","", colnames(cytoRNA_Naive_all))
colnames(cytoRNA_Naive_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.Naive.b..","", colnames(cytoRNA_Naive_all))
colnames(cytoRNA_Naive_all) <- gsub(".bam","", colnames(cytoRNA_Naive_all))

# Switch b2 and b3
colnames(cytoRNA_Naive_all) <- gsub("b2","b4", colnames(cytoRNA_Naive_all))
colnames(cytoRNA_Naive_all) <- gsub("b3","b2", colnames(cytoRNA_Naive_all))
colnames(cytoRNA_Naive_all) <- gsub("b4","b3", colnames(cytoRNA_Naive_all))


cytoRNA_Naive_b1 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b1/D40000/exons/exon_counts.txt', header=T)
cytoRNA_Naive_b2 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b2/D40000/exons/exon_counts.txt', header=T)
cytoRNA_Naive_b3 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b3/D40000/exons/exon_counts.txt', header=T)
cytoRNA_Naive_exons <- data.frame(cytoRNA_Naive_b1[,1:6],
                                cytoRNA_Naive_b1[, -c(1:6)],
                                cytoRNA_Naive_b2[, -c(1:6)],
                                cytoRNA_Naive_b3[, -c(1:6)])
rm(cytoRNA_Naive_b1, cytoRNA_Naive_b2, cytoRNA_Naive_b3)

colnames(cytoRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.cytoplasmic.Naive.b..","", colnames(cytoRNA_Naive_exons))
colnames(cytoRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.Naive.b..","", colnames(cytoRNA_Naive_exons))
colnames(cytoRNA_Naive_exons) <- gsub(".bam","", colnames(cytoRNA_Naive_exons))

# Switch b2 and b3
colnames(cytoRNA_Naive_exons) <- gsub("b2","b4", colnames(cytoRNA_Naive_exons))
colnames(cytoRNA_Naive_exons) <- gsub("b3","b2", colnames(cytoRNA_Naive_exons))
colnames(cytoRNA_Naive_exons) <- gsub("b4","b3", colnames(cytoRNA_Naive_exons))

save(cytoRNA_Naive_all, cytoRNA_Naive_exons, file='Data/raw_counts/cytoRNA_Naive.Rdata')


cytoRNA_LPA_b1 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b1/D40000/all/all_counts.txt', header=T)
cytoRNA_LPA_b2 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b2/D40000/all/all_counts.txt', header=T)
cytoRNA_LPA_b3 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b3/D40000/all/all_counts.txt', header=T)
cytoRNA_LPA_all <- data.frame(cytoRNA_LPA_b1[,1:6],
                            cytoRNA_LPA_b1[, -c(1:6)],
                            cytoRNA_LPA_b2[, -c(1:6)],
                            cytoRNA_LPA_b3[, -c(1:6)])
rm(cytoRNA_LPA_b1, cytoRNA_LPA_b2, cytoRNA_LPA_b3)

colnames(cytoRNA_LPA_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.LPA.b..","", colnames(cytoRNA_LPA_all))
colnames(cytoRNA_LPA_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.LPA.b..","", colnames(cytoRNA_LPA_all))
colnames(cytoRNA_LPA_all) <- gsub(".bam","", colnames(cytoRNA_LPA_all))

# Switch b2 and b3
colnames(cytoRNA_LPA_all) <- gsub("b2","b4", colnames(cytoRNA_LPA_all))
colnames(cytoRNA_LPA_all) <- gsub("b3","b2", colnames(cytoRNA_LPA_all))
colnames(cytoRNA_LPA_all) <- gsub("b4","b3", colnames(cytoRNA_LPA_all))


cytoRNA_LPA_b1 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b1/D40000/exons/exon_counts.txt', header=T)
cytoRNA_LPA_b2 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b2/D40000/exons/exon_counts.txt', header=T)
cytoRNA_LPA_b3 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b3/D40000/exons/exon_counts.txt', header=T)
cytoRNA_LPA_exons <- data.frame(cytoRNA_LPA_b1[,1:6],
                              cytoRNA_LPA_b1[, -c(1:6)],
                              cytoRNA_LPA_b2[, -c(1:6)],
                              cytoRNA_LPA_b3[, -c(1:6)])
rm(cytoRNA_LPA_b1, cytoRNA_LPA_b2, cytoRNA_LPA_b3)

colnames(cytoRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.cytoplasmic.LPA.b..","", colnames(cytoRNA_LPA_exons))
colnames(cytoRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.LPA.b..","", colnames(cytoRNA_LPA_exons))
colnames(cytoRNA_LPA_exons) <- gsub(".bam","", colnames(cytoRNA_LPA_exons))

colnames(cytoRNA_LPA_exons) <- gsub("b2","b4", colnames(cytoRNA_LPA_exons))
colnames(cytoRNA_LPA_exons) <- gsub("b3","b2", colnames(cytoRNA_LPA_exons))
colnames(cytoRNA_LPA_exons) <- gsub("b4","b3", colnames(cytoRNA_LPA_exons))

save(cytoRNA_LPA_all, cytoRNA_LPA_exons, file='Data/raw_counts/cytoRNA_LPA.Rdata')


## load polyA data

polyARNA_Naive_b1 <- read.table(file = 'merged_counts/polyA.Naive.b1/D40000/all/all_counts.txt', header=T)
polyARNA_Naive_b2 <- read.table(file = 'merged_counts/polyA.Naive.b2/D40000/all/all_counts.txt', header=T)
polyARNA_Naive_b3 <- read.table(file = 'merged_counts/polyA.Naive.b3/D40000/all/all_counts.txt', header=T)
polyARNA_Naive_all <- data.frame(polyARNA_Naive_b1[,1:6],
                              polyARNA_Naive_b1[, -c(1:6)],
                              polyARNA_Naive_b2[, -c(1:6)],
                              polyARNA_Naive_b3[, -c(1:6)])
rm(polyARNA_Naive_b1, polyARNA_Naive_b2, polyARNA_Naive_b3)

colnames(polyARNA_Naive_all) <- gsub("merged_bams_mm10_vM14.polyA.Naive.b..","", colnames(polyARNA_Naive_all))
colnames(polyARNA_Naive_all) <- gsub(".bam","", colnames(polyARNA_Naive_all))


polyARNA_Naive_b1 <- read.table(file = 'merged_counts/polyA.Naive.b1/D40000/exons/exon_counts.txt', header=T)
polyARNA_Naive_b2 <- read.table(file = 'merged_counts/polyA.Naive.b2/D40000/exons/exon_counts.txt', header=T)
polyARNA_Naive_b3 <- read.table(file = 'merged_counts/polyA.Naive.b3/D40000/exons/exon_counts.txt', header=T)
polyARNA_Naive_exons <- data.frame(polyARNA_Naive_b1[,1:6],
                                polyARNA_Naive_b1[, -c(1:6)],
                                polyARNA_Naive_b2[, -c(1:6)],
                                polyARNA_Naive_b3[, -c(1:6)])
rm(polyARNA_Naive_b1, polyARNA_Naive_b2, polyARNA_Naive_b3)

colnames(polyARNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.polyA.Naive.b..","", colnames(polyARNA_Naive_exons))
colnames(polyARNA_Naive_exons) <- gsub(".bam","", colnames(polyARNA_Naive_exons))

save(polyARNA_Naive_all, polyARNA_Naive_exons, file='Data/raw_counts/polyARNA_Naive.Rdata')

polyARNA_LPA_b1 <- read.table(file = 'merged_counts/polyA.LPA.b1/D40000/all/all_counts.txt', header=T)
polyARNA_LPA_b2 <- read.table(file = 'merged_counts/polyA.LPA.b2/D40000/all/all_counts.txt', header=T)
polyARNA_LPA_b3 <- read.table(file = 'merged_counts/polyA.LPA.b3/D40000/all/all_counts.txt', header=T)
polyARNA_LPA_all <- data.frame(polyARNA_LPA_b1[,1:6],
                            polyARNA_LPA_b1[, -c(1:6)],
                            polyARNA_LPA_b2[, -c(1:6)],
                            polyARNA_LPA_b3[, -c(1:6)])
rm(polyARNA_LPA_b1, polyARNA_LPA_b2, polyARNA_LPA_b3)

colnames(polyARNA_LPA_all) <- gsub("merged_bams_mm10_vM14.polyA.LPA.b..","", colnames(polyARNA_LPA_all))
colnames(polyARNA_LPA_all) <- gsub(".bam","", colnames(polyARNA_LPA_all))


polyARNA_LPA_b1 <- read.table(file = 'merged_counts/polyA.LPA.b1/D40000/exons/exon_counts.txt', header=T)
polyARNA_LPA_b2 <- read.table(file = 'merged_counts/polyA.LPA.b2/D40000/exons/exon_counts.txt', header=T)
polyARNA_LPA_b3 <- read.table(file = 'merged_counts/polyA.LPA.b3/D40000/exons/exon_counts.txt', header=T)
polyARNA_LPA_exons <- data.frame(polyARNA_LPA_b1[,1:6],
                              polyARNA_LPA_b1[, -c(1:6)],
                              polyARNA_LPA_b2[, -c(1:6)],
                              polyARNA_LPA_b3[, -c(1:6)])
rm(polyARNA_LPA_b1, polyARNA_LPA_b2, polyARNA_LPA_b3)

colnames(polyARNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.polyA.LPA.b..","", colnames(polyARNA_LPA_exons))
colnames(polyARNA_LPA_exons) <- gsub(".bam","", colnames(polyARNA_LPA_exons))

save(polyARNA_LPA_all, polyARNA_LPA_exons, file='Data/raw_counts/polyARNA_LPA.Rdata')

## Create exons rpkms
load(file='Data/raw_counts/caRNA_Naive.Rdata')
Time_ca_Naive <- as.numeric(sapply(colnames(caRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_Naive <- sapply(colnames(caRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_Naive <- as.factor(Time_ca_Naive)

load(file='Data/raw_counts/npRNA_Naive.Rdata')
Time_np_Naive <- as.numeric(sapply(colnames(npRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_Naive <- sapply(colnames(npRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_Naive <- as.factor(Time_np_Naive)

load(file='Data/raw_counts/cytoRNA_Naive.Rdata')
Time_cyto_Naive <- as.numeric(sapply(colnames(cytoRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_Naive <- sapply(colnames(cytoRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_Naive <- as.factor(Time_cyto_Naive)

load(file='Data/raw_counts/polyARNA_Naive.Rdata')
Time_polyA_Naive <- as.numeric(sapply(colnames(polyARNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_polyA_Naive <- sapply(colnames(polyARNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_polyA_Naive <- as.factor(Time_polyA_Naive)

library(edgeR)
caRNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6])), log=T)
npRNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6])), log=T)
cytoRNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6])), log=T)
polyARNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = polyARNA_Naive_exons[,-c(1:6)], group = group_polyA_Naive, genes = polyARNA_Naive_exons[,1:6])), log=T)

row.names(caRNA_Naive_exons_cpm) <- row.names(npRNA_Naive_exons_cpm) <- row.names(cytoRNA_Naive_exons_cpm) <- row.names(polyARNA_Naive_exons_cpm) <- caRNA_Naive_exons$Geneid

save(caRNA_Naive_exons_cpm, file='Data/exons_cpms/caRNA_Naive_exons_cpm.Rdata')
save(npRNA_Naive_exons_cpm, file='Data/exons_cpms/npRNA_Naive_exons_cpm.Rdata')
save(cytoRNA_Naive_exons_cpm, file='Data/exons_cpms/cytoRNA_Naive_exons_cpm.Rdata')
save(polyARNA_Naive_exons_cpm, file='Data/exons_cpms/polyARNA_Naive_exons_cpm.Rdata')

caRNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6])), log=T)
npRNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6])), log=T)
cytoRNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6])), log=T)
polyARNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = polyARNA_Naive_exons[,-c(1:6)], group = group_polyA_Naive, genes = polyARNA_Naive_exons[,1:6])), log=T)

row.names(caRNA_Naive_exons_rpkm) <- row.names(npRNA_Naive_exons_rpkm) <- row.names(cytoRNA_Naive_exons_rpkm) <- row.names(polyARNA_Naive_exons_rpkm) <- caRNA_Naive_exons$Geneid

save(caRNA_Naive_exons_rpkm, file='Data/exons_rpkms/caRNA_Naive_exons_rpkm.Rdata')
save(npRNA_Naive_exons_rpkm, file='Data/exons_rpkms/npRNA_Naive_exons_rpkm.Rdata')
save(cytoRNA_Naive_exons_rpkm, file='Data/exons_rpkms/cytoRNA_Naive_exons_rpkm.Rdata')
save(polyARNA_Naive_exons_rpkm, file='Data/exons_rpkms/polyARNA_Naive_exons_rpkm.Rdata')

# get lib size and scaling factor
tmp <- calcNormFactors(DGEList(counts = caRNA_Naive_all[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_all[,1:6]))
lib_all <- data.frame( samples = row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors)
tmp <- calcNormFactors(DGEList(counts = npRNA_Naive_all[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_Naive_all[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = polyARNA_Naive_all[,-c(1:6)], group = group_polyA_Naive, genes = polyARNA_Naive_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

tmp <- calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6]))
lib_exons <- data.frame( samples = row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors)
tmp <- calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = polyARNA_Naive_exons[,-c(1:6)], group = group_polyA_Naive, genes = polyARNA_Naive_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))


load(file='Data/raw_counts/caRNA_LPA.Rdata')
Time_ca_LPA <- as.numeric(sapply(colnames(caRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_LPA <- sapply(colnames(caRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_LPA <- as.factor(Time_ca_LPA)

load(file='Data/raw_counts/npRNA_LPA.Rdata')
Time_np_LPA <- as.numeric(sapply(colnames(npRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_LPA <- sapply(colnames(npRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_LPA <- as.factor(Time_np_LPA)

load(file='Data/raw_counts/cytoRNA_LPA.Rdata')
Time_cyto_LPA <- as.numeric(sapply(colnames(cytoRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_LPA <- sapply(colnames(cytoRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_LPA <- as.factor(Time_cyto_LPA)

load(file='Data/raw_counts/polyARNA_LPA.Rdata')
Time_polyA_LPA <- as.numeric(sapply(colnames(polyARNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_polyA_LPA <- sapply(colnames(polyARNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_polyA_LPA <- as.factor(Time_polyA_LPA)

caRNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6])), log=T)
npRNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6])), log=T)
cytoRNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6])), log=T)
polyARNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = polyARNA_LPA_exons[,-c(1:6)], group = group_polyA_LPA, genes = polyARNA_LPA_exons[,1:6])), log=T)

row.names(caRNA_LPA_exons_cpm) <- row.names(npRNA_LPA_exons_cpm) <- row.names(cytoRNA_LPA_exons_cpm) <- row.names(polyARNA_LPA_exons_cpm) <- caRNA_LPA_exons$Geneid

save(caRNA_LPA_exons_cpm, file='Data/exons_cpms/caRNA_LPA_exons_cpm.Rdata')
save(npRNA_LPA_exons_cpm, file='Data/exons_cpms/npRNA_LPA_exons_cpm.Rdata')
save(cytoRNA_LPA_exons_cpm, file='Data/exons_cpms/cytoRNA_LPA_exons_cpm.Rdata')
save(polyARNA_LPA_exons_cpm, file='Data/exons_cpms/polyARNA_LPA_exons_cpm.Rdata')


caRNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6])), log=T)
npRNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6])), log=T)
cytoRNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6])), log=T)
polyARNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = polyARNA_LPA_exons[,-c(1:6)], group = group_polyA_LPA, genes = polyARNA_LPA_exons[,1:6])), log=T)

row.names(caRNA_LPA_exons_rpkm) <- row.names(npRNA_LPA_exons_rpkm) <- row.names(cytoRNA_LPA_exons_rpkm) <- row.names(polyARNA_LPA_exons_rpkm) <- caRNA_LPA_exons$Geneid

save(caRNA_LPA_exons_rpkm, file='Data/exons_rpkms/caRNA_LPA_exons_rpkm.Rdata')
save(npRNA_LPA_exons_rpkm, file='Data/exons_rpkms/npRNA_LPA_exons_rpkm.Rdata')
save(cytoRNA_LPA_exons_rpkm, file='Data/exons_rpkms/cytoRNA_LPA_exons_rpkm.Rdata')
save(polyARNA_LPA_exons_rpkm, file='Data/exons_rpkms/polyARNA_LPA_exons_rpkm.Rdata')

tmp <- calcNormFactors(DGEList(counts = caRNA_LPA_all[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = npRNA_LPA_all[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_LPA_all[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = polyARNA_LPA_all[,-c(1:6)], group = group_polyA_LPA, genes = polyARNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

tmp <- calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = polyARNA_LPA_exons[,-c(1:6)], group = group_polyA_LPA, genes = polyARNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

save(lib_all, lib_exons, file='Data/raw_counts/lib.size.Rdata')

########################################################################################################################
# Create dataset with last 5kb exons
caRNA_Naive_exons_5kb_b1 <- read.table(file = 'merged_counts/chromatin.Naive.b1/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
caRNA_Naive_exons_5kb_b2 <- read.table(file = 'merged_counts/chromatin.Naive.b2/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
caRNA_Naive_exons_5kb_b3 <- read.table(file = 'merged_counts/chromatin.Naive.b3/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)

caRNA_Naive_exons_5kb_all <- data.frame(caRNA_Naive_exons_5kb_b1[,1:6],
                                        caRNA_Naive_exons_5kb_b1[, -c(1:6)],
                                        caRNA_Naive_exons_5kb_b2[, -c(1:6)],
                                        caRNA_Naive_exons_5kb_b3[, -c(1:6)])
rm(caRNA_Naive_exons_5kb_b1, caRNA_Naive_exons_5kb_b2, caRNA_Naive_exons_5kb_b3)

colnames(caRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.chromatin.Naive.b..","", colnames(caRNA_Naive_exons_5kb_all))
colnames(caRNA_Naive_exons_5kb_all) <- gsub(".bam","", colnames(caRNA_Naive_exons_5kb_all))

save(caRNA_Naive_exons_5kb_all, file='Data/raw_counts/caRNA_Naive_5kb.Rdata')

npRNA_Naive_exons_5kb_b1 <- read.table(file = 'merged_counts/Nuc.Naive.b1/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
npRNA_Naive_exons_5kb_b3 <- read.table(file = 'merged_counts/Nuc.Naive.b3/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)

npRNA_Naive_exons_5kb_all <- data.frame(npRNA_Naive_exons_5kb_b1[,1:6],
                                        npRNA_Naive_exons_5kb_b1[, -c(1:6)],
                                        npRNA_Naive_exons_5kb_b3[, -c(1:6)])
rm(npRNA_Naive_exons_5kb_b1, npRNA_Naive_exons_5kb_b3)

colnames(npRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.Nuc.Naive.b..","", colnames(npRNA_Naive_exons_5kb_all))
colnames(npRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.Naive.b..","", colnames(npRNA_Naive_exons_5kb_all))
colnames(npRNA_Naive_exons_5kb_all) <- gsub(".bam","", colnames(npRNA_Naive_exons_5kb_all))

save(npRNA_Naive_exons_5kb_all, file='Data/raw_counts/npRNA_Naive_5kb.Rdata')

cytoRNA_Naive_exons_5kb_b1 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b1/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
cytoRNA_Naive_exons_5kb_b2 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b2/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
cytoRNA_Naive_exons_5kb_b3 <- read.table(file = 'merged_counts/cytoplasmic.Naive.b3/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)

cytoRNA_Naive_exons_5kb_all <- data.frame(cytoRNA_Naive_exons_5kb_b1[,1:6],
                                          cytoRNA_Naive_exons_5kb_b1[, -c(1:6)],
                                          cytoRNA_Naive_exons_5kb_b2[, -c(1:6)],
                                          cytoRNA_Naive_exons_5kb_b3[, -c(1:6)])
rm(cytoRNA_Naive_exons_5kb_b1, cytoRNA_Naive_exons_5kb_b2, cytoRNA_Naive_exons_5kb_b3)

colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.Naive.b..","", colnames(cytoRNA_Naive_exons_5kb_all))
colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.Naive.b..","", colnames(cytoRNA_Naive_exons_5kb_all))
colnames(cytoRNA_Naive_exons_5kb_all) <- gsub(".bam","", colnames(cytoRNA_Naive_exons_5kb_all))

colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("b2","b4", colnames(cytoRNA_Naive_exons_5kb_all))
colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("b3","b2", colnames(cytoRNA_Naive_exons_5kb_all))
colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("b4","b3", colnames(cytoRNA_Naive_exons_5kb_all))

save(cytoRNA_Naive_exons_5kb_all, file='Data/raw_counts/cytoRNA_Naive_5kb.Rdata')

## Create exons rpkms (normalised using total exons library size)
load(file='Data/raw_counts/caRNA_Naive_5kb.Rdata')
load(file='Data/raw_counts/caRNA_Naive.Rdata')

Time_ca_Naive_5kb <- as.numeric(sapply(colnames(caRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_Naive_5kb <- sapply(colnames(caRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_Naive_5kb <- as.factor(Time_ca_Naive_5kb)

load(file='Data/raw_counts/npRNA_Naive_5kb.Rdata')
load(file='Data/raw_counts/npRNA_Naive.Rdata')
Time_np_Naive_5kb <- as.numeric(sapply(colnames(npRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_Naive_5kb <- sapply(colnames(npRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_Naive_5kb <- as.factor(Time_np_Naive_5kb)

load(file='Data/raw_counts/cytoRNA_Naive_5kb.Rdata')
load(file='Data/raw_counts/cytoRNA_Naive.Rdata')
Time_cyto_Naive_5kb <- as.numeric(sapply(colnames(cytoRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_Naive_5kb <- sapply(colnames(cytoRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_Naive_5kb <- as.factor(Time_cyto_Naive_5kb)

library(edgeR)
y_ca <- calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6]))
y_ca_5kb <- DGEList(counts = caRNA_Naive_exons_5kb_all[,-c(1:6)], group = group_ca_Naive_5kb, genes = caRNA_Naive_exons_5kb_all[,1:6])
y_ca_5kb$samples <- y_ca$samples
caRNA_Naive_exons_rpkm_5kb <- rpkm(y_ca_5kb, log=T)
caRNA_Naive_exons_cpm_5kb <- cpm(y_ca_5kb, log=T)

y_np <- calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6]))
y_np_5kb <- DGEList(counts = npRNA_Naive_exons_5kb_all[,-c(1:6)], group = group_np_Naive_5kb, genes = npRNA_Naive_exons_5kb_all[,1:6])
y_np_5kb$samples <- y_np$samples
npRNA_Naive_exons_rpkm_5kb <- rpkm(y_np_5kb, log=T)
npRNA_Naive_exons_cpm_5kb <- cpm(y_np_5kb, log=T)

y_cyto <- calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6]))
y_cyto_5kb <- DGEList(counts = cytoRNA_Naive_exons_5kb_all[,-c(1:6)], group = group_cyto_Naive_5kb, genes = cytoRNA_Naive_exons_5kb_all[,1:6])
y_cyto_5kb$samples <- y_cyto$samples
cytoRNA_Naive_exons_rpkm_5kb <- rpkm(y_cyto_5kb, log=T)
cytoRNA_Naive_exons_cpm_5kb <- cpm(y_cyto_5kb, log=T)

row.names(caRNA_Naive_exons_rpkm_5kb) <- row.names(npRNA_Naive_exons_rpkm_5kb) <- row.names(cytoRNA_Naive_exons_rpkm_5kb) <- caRNA_Naive_exons_5kb_all$Geneid
row.names(caRNA_Naive_exons_cpm_5kb) <- row.names(npRNA_Naive_exons_cpm_5kb) <- row.names(cytoRNA_Naive_exons_cpm_5kb) <- caRNA_Naive_exons_5kb_all$Geneid

save(caRNA_Naive_exons_rpkm_5kb, file='Data/exons_rpkms/caRNA_Naive_exons_rpkm_5kb.Rdata')
save(npRNA_Naive_exons_rpkm_5kb, file='Data/exons_rpkms/npRNA_Naive_exons_rpkm_5kb.Rdata')
save(cytoRNA_Naive_exons_rpkm_5kb, file='Data/exons_rpkms/cytoRNA_Naive_exons_rpkm_5kb.Rdata')

save(caRNA_Naive_exons_cpm_5kb, file='Data/exons_cpms/caRNA_Naive_exons_cpm_5kb.Rdata')
save(npRNA_Naive_exons_cpm_5kb, file='Data/exons_cpms/npRNA_Naive_exons_cpm_5kb.Rdata')
save(cytoRNA_Naive_exons_cpm_5kb, file='Data/exons_cpms/cytoRNA_Naive_exons_cpm_5kb.Rdata')



## LPA Create dataset with last 5kb exons
caRNA_LPA_exons_5kb_b1 <- read.table(file = 'merged_counts/chromatin.LPA.b1/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
caRNA_LPA_exons_5kb_b2 <- read.table(file = 'merged_counts/chromatin.LPA.b2/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
caRNA_LPA_exons_5kb_b3 <- read.table(file = 'merged_counts/chromatin.LPA.b3/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
caRNA_LPA_exons_5kb_all <- data.frame(caRNA_LPA_exons_5kb_b1[,1:6],
                                      caRNA_LPA_exons_5kb_b1[, -c(1:6)],
                                      caRNA_LPA_exons_5kb_b2[, -c(1:6)],
                                      caRNA_LPA_exons_5kb_b3[, -c(1:6)])
rm(caRNA_LPA_exons_5kb_b1, caRNA_LPA_exons_5kb_b2, caRNA_LPA_exons_5kb_b3)

colnames(caRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.chromatin.LPA.b..","", colnames(caRNA_LPA_exons_5kb_all))
colnames(caRNA_LPA_exons_5kb_all) <- gsub(".bam","", colnames(caRNA_LPA_exons_5kb_all))

save(caRNA_LPA_exons_5kb_all, file='Data/raw_counts/caRNA_LPA_5kb.Rdata')

npRNA_LPA_exons_5kb_b1 <- read.table(file = 'merged_counts/Nuc.LPA.b1/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
npRNA_LPA_exons_5kb_b3 <- read.table(file = 'merged_counts/Nuc.LPA.b3/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
npRNA_LPA_exons_5kb_all <- data.frame(npRNA_LPA_exons_5kb_b1[,1:6],
                                      npRNA_LPA_exons_5kb_b1[, -c(1:6)],
                                      npRNA_LPA_exons_5kb_b3[, -c(1:6)])
rm(npRNA_LPA_exons_5kb_b1, npRNA_LPA_exons_5kb_b3)

colnames(npRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.Nuc.LPA.b..","", colnames(npRNA_LPA_exons_5kb_all))
colnames(npRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.LPA.b..","", colnames(npRNA_LPA_exons_5kb_all))
colnames(npRNA_LPA_exons_5kb_all) <- gsub(".bam","", colnames(npRNA_LPA_exons_5kb_all))

save(npRNA_LPA_exons_5kb_all, file='Data/raw_counts/npRNA_LPA_5kb.Rdata')

cytoRNA_LPA_exons_5kb_b1 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b1/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
cytoRNA_LPA_exons_5kb_b2 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b2/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
cytoRNA_LPA_exons_5kb_b3 <- read.table(file = 'merged_counts/cytoplasmic.LPA.b3/D40000_LPAgenes_last5kb/exon_counts.txt', header=T)
cytoRNA_LPA_exons_5kb_all <- data.frame(cytoRNA_LPA_exons_5kb_b1[,1:6],
                                        cytoRNA_LPA_exons_5kb_b1[, -c(1:6)],
                                        cytoRNA_LPA_exons_5kb_b2[, -c(1:6)],
                                        cytoRNA_LPA_exons_5kb_b3[, -c(1:6)])
rm(cytoRNA_LPA_exons_5kb_b1, cytoRNA_LPA_exons_5kb_b2, cytoRNA_LPA_exons_5kb_b3)

colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.LPA.b..","", colnames(cytoRNA_LPA_exons_5kb_all))
colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.LPA.b..","", colnames(cytoRNA_LPA_exons_5kb_all))
colnames(cytoRNA_LPA_exons_5kb_all) <- gsub(".bam","", colnames(cytoRNA_LPA_exons_5kb_all))

colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("b2","b4", colnames(cytoRNA_LPA_exons_5kb_all))
colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("b3","b2", colnames(cytoRNA_LPA_exons_5kb_all))
colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("b4","b3", colnames(cytoRNA_LPA_exons_5kb_all))

save(cytoRNA_LPA_exons_5kb_all, file='Data/raw_counts/cytoRNA_LPA_5kb.Rdata')

## Create exons rpkms (normalised using total exons library size)
load(file='Data/raw_counts/caRNA_LPA_5kb.Rdata')
load(file='Data/raw_counts/caRNA_LPA.Rdata')

Time_ca_LPA_5kb <- as.numeric(sapply(colnames(caRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_LPA_5kb <- sapply(colnames(caRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_LPA_5kb <- as.factor(Time_ca_LPA_5kb)

load(file='Data/raw_counts/npRNA_LPA_5kb.Rdata')
load(file='Data/raw_counts/npRNA_LPA.Rdata')
Time_np_LPA_5kb <- as.numeric(sapply(colnames(npRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_LPA_5kb <- sapply(colnames(npRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_LPA_5kb <- as.factor(Time_np_LPA_5kb)

load(file='Data/raw_counts/cytoRNA_LPA_5kb.Rdata')
load(file='Data/raw_counts/cytoRNA_LPA.Rdata')
Time_cyto_LPA_5kb <- as.numeric(sapply(colnames(cytoRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_LPA_5kb <- sapply(colnames(cytoRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_LPA_5kb <- as.factor(Time_cyto_LPA_5kb)

library(edgeR)
y_ca <- calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6]))
y_ca_5kb <- DGEList(counts = caRNA_LPA_exons_5kb_all[,-c(1:6)], group = group_ca_LPA_5kb, genes = caRNA_LPA_exons_5kb_all[,1:6])
y_ca_5kb$samples <- y_ca$samples
caRNA_LPA_exons_rpkm_5kb <- rpkm(y_ca_5kb, log=T)
caRNA_LPA_exons_cpm_5kb <- cpm(y_ca_5kb, log=T)

y_np <- calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6]))
y_np_5kb <- DGEList(counts = npRNA_LPA_exons_5kb_all[,-c(1:6)], group = group_np_LPA_5kb, genes = npRNA_LPA_exons_5kb_all[,1:6])
y_np_5kb$samples <- y_np$samples
npRNA_LPA_exons_rpkm_5kb <- rpkm(y_np_5kb, log=T)
npRNA_LPA_exons_cpm_5kb <- cpm(y_np_5kb, log=T)

y_cyto <- calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6]))
y_cyto_5kb <- DGEList(counts = cytoRNA_LPA_exons_5kb_all[,-c(1:6)], group = group_cyto_LPA_5kb, genes = cytoRNA_LPA_exons_5kb_all[,1:6])
y_cyto_5kb$samples <- y_cyto$samples
cytoRNA_LPA_exons_rpkm_5kb <- rpkm(y_cyto_5kb, log=T)
cytoRNA_LPA_exons_cpm_5kb <- cpm(y_cyto_5kb, log=T)

row.names(caRNA_LPA_exons_rpkm_5kb) <- row.names(npRNA_LPA_exons_rpkm_5kb) <- row.names(cytoRNA_LPA_exons_rpkm_5kb) <- caRNA_LPA_exons_5kb_all$Geneid
row.names(caRNA_LPA_exons_cpm_5kb) <- row.names(npRNA_LPA_exons_cpm_5kb) <- row.names(cytoRNA_LPA_exons_cpm_5kb) <- caRNA_LPA_exons_5kb_all$Geneid

save(caRNA_LPA_exons_rpkm_5kb, file='Data/exons_rpkms/caRNA_LPA_exons_rpkm_5kb.Rdata')
save(npRNA_LPA_exons_rpkm_5kb, file='Data/exons_rpkms/npRNA_LPA_exons_rpkm_5kb.Rdata')
save(cytoRNA_LPA_exons_rpkm_5kb, file='Data/exons_rpkms/cytoRNA_LPA_exons_rpkm_5kb.Rdata')

save(caRNA_LPA_exons_cpm_5kb, file='Data/exons_cpms/caRNA_LPA_exons_cpm_5kb.Rdata')
save(npRNA_LPA_exons_cpm_5kb, file='Data/exons_cpms/npRNA_LPA_exons_cpm_5kb.Rdata')
save(cytoRNA_LPA_exons_cpm_5kb, file='Data/exons_cpms/cytoRNA_LPA_exons_cpm_5kb.Rdata')

