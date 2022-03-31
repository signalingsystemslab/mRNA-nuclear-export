# Check working Directory -> must be Data folder

# create directories
dir.create('raw_counts', showWarnings = F)
dir.create('exons_rpkms', showWarnings = F)
dir.create('exons_cpms', showWarnings = F)

# load caRNA data
caRNA_Naive_rep1 <- read.table(file = 'merged_counts/chromatin.Naive.rep1/D40000/all/all_counts.txt.gz', header=T)
caRNA_Naive_rep2 <- read.table(file = 'merged_counts/chromatin.Naive.rep2/D40000/all/all_counts.txt.gz', header=T)
caRNA_Naive_rep3 <- read.table(file = 'merged_counts/chromatin.Naive.rep3/D40000/all/all_counts.txt.gz', header=T)
caRNA_Naive_all <- data.frame(caRNA_Naive_rep1[, 1:6],
                              caRNA_Naive_rep1[, -c(1:6)],
                              caRNA_Naive_rep2[, -c(1:6)],
                              caRNA_Naive_rep3[, -c(1:6)])
rm(caRNA_Naive_rep1, caRNA_Naive_rep2, caRNA_Naive_rep3)

colnames(caRNA_Naive_all) <- gsub("merged_bams_mm10_vM14.chromatin.Naive.rep..","", colnames(caRNA_Naive_all))
colnames(caRNA_Naive_all) <- gsub(".bam","", colnames(caRNA_Naive_all))

caRNA_Naive_rep1 <- read.table(file = 'merged_counts/chromatin.Naive.rep1/D40000/exons/exon_counts.txt.gz', header=T)
caRNA_Naive_rep2 <- read.table(file = 'merged_counts/chromatin.Naive.rep2/D40000/exons/exon_counts.txt.gz', header=T)
caRNA_Naive_rep3 <- read.table(file = 'merged_counts/chromatin.Naive.rep3/D40000/exons/exon_counts.txt.gz', header=T)
caRNA_Naive_exons <- data.frame(caRNA_Naive_rep1[,1:6],
                                caRNA_Naive_rep1[, -c(1:6)],
                                caRNA_Naive_rep2[, -c(1:6)],
                                caRNA_Naive_rep3[, -c(1:6)])
rm(caRNA_Naive_rep1, caRNA_Naive_rep2, caRNA_Naive_rep3)

colnames(caRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.chromatin.Naive.rep..","", colnames(caRNA_Naive_exons))
colnames(caRNA_Naive_exons) <- gsub(".bam","", colnames(caRNA_Naive_exons))

save(caRNA_Naive_all, caRNA_Naive_exons, file='raw_counts/caRNA_Naive.Rdata')

caRNA_LPA_rep1 <- read.table(file = 'merged_counts/chromatin.LPA.rep1/D40000/all/all_counts.txt.gz', header=T)
caRNA_LPA_rep2 <- read.table(file = 'merged_counts/chromatin.LPA.rep2/D40000/all/all_counts.txt.gz', header=T)
caRNA_LPA_rep3 <- read.table(file = 'merged_counts/chromatin.LPA.rep3/D40000/all/all_counts.txt.gz', header=T)
caRNA_LPA_all <- data.frame(caRNA_LPA_rep1[,1:6],
                            caRNA_LPA_rep1[, -c(1:6)],
                            caRNA_LPA_rep2[, -c(1:6)],
                            caRNA_LPA_rep3[, -c(1:6)])
rm(caRNA_LPA_rep1, caRNA_LPA_rep2, caRNA_LPA_rep3)

colnames(caRNA_LPA_all) <- gsub("merged_bams_mm10_vM14.chromatin.LPA.rep..","", colnames(caRNA_LPA_all))
colnames(caRNA_LPA_all) <- gsub(".bam","", colnames(caRNA_LPA_all))

caRNA_LPA_rep1 <- read.table(file = 'merged_counts/chromatin.LPA.rep1/D40000/exons/exon_counts.txt.gz', header=T)
caRNA_LPA_rep2 <- read.table(file = 'merged_counts/chromatin.LPA.rep2/D40000/exons/exon_counts.txt.gz', header=T)
caRNA_LPA_rep3 <- read.table(file = 'merged_counts/chromatin.LPA.rep3/D40000/exons/exon_counts.txt.gz', header=T)
caRNA_LPA_exons <- data.frame(caRNA_LPA_rep1[,1:6],
                              caRNA_LPA_rep1[, -c(1:6)],
                              caRNA_LPA_rep2[, -c(1:6)],
                              caRNA_LPA_rep3[, -c(1:6)])
rm(caRNA_LPA_rep1, caRNA_LPA_rep2, caRNA_LPA_rep3)

colnames(caRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.chromatin.LPA.rep..","", colnames(caRNA_LPA_exons))
colnames(caRNA_LPA_exons) <- gsub(".bam","", colnames(caRNA_LPA_exons))

save(caRNA_LPA_all, caRNA_LPA_exons, file='raw_counts/caRNA_LPA.Rdata')

## load nucleoplasmic data

npRNA_Naive_rep1 <- read.table(file = 'merged_counts/Nuc.Naive.rep1/D40000/all/all_counts.txt.gz', header=T)
npRNA_Naive_rep2 <- read.table(file = 'merged_counts/Nuc.Naive.rep2/D40000/all/all_counts.txt.gz', header=T)
npRNA_Naive_all <- data.frame(npRNA_Naive_rep1[,1:6],
                              npRNA_Naive_rep1[, -c(1:6)],
                              npRNA_Naive_rep2[, -c(1:6)])
rm(npRNA_Naive_rep1, npRNA_Naive_rep2)

colnames(npRNA_Naive_all) <- gsub("merged_bams_mm10_vM14.Nuc.Naive.rep..","", colnames(npRNA_Naive_all))
colnames(npRNA_Naive_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.Naive.rep..","", colnames(npRNA_Naive_all))
colnames(npRNA_Naive_all) <- gsub(".bam","", colnames(npRNA_Naive_all))

npRNA_Naive_rep1 <- read.table(file = 'merged_counts/Nuc.Naive.rep1/D40000/exons/exon_counts.txt.gz', header=T)
npRNA_Naive_rep2 <- read.table(file = 'merged_counts/Nuc.Naive.rep2/D40000/exons/exon_counts.txt.gz', header=T)
npRNA_Naive_exons <- data.frame(npRNA_Naive_rep1[,1:6],
                                npRNA_Naive_rep1[, -c(1:6)],
                                npRNA_Naive_rep2[, -c(1:6)])
rm(npRNA_Naive_rep1, npRNA_Naive_rep2)

colnames(npRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.Nuc.Naive.rep..","", colnames(npRNA_Naive_exons))
colnames(npRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14_2.Nuc.Naive.rep..","", colnames(npRNA_Naive_exons))
colnames(npRNA_Naive_exons) <- gsub(".bam","", colnames(npRNA_Naive_exons))

save(npRNA_Naive_all, npRNA_Naive_exons, file='raw_counts/npRNA_Naive.Rdata')

npRNA_LPA_rep1 <- read.table(file = 'merged_counts/Nuc.LPA.rep1/D40000/all/all_counts.txt.gz', header=T)
npRNA_LPA_rep2 <- read.table(file = 'merged_counts/Nuc.LPA.rep2/D40000/all/all_counts.txt.gz', header=T)
npRNA_LPA_all <- data.frame(npRNA_LPA_rep1[,1:6],
                            npRNA_LPA_rep1[, -c(1:6)],
                            npRNA_LPA_rep2[, -c(1:6)])
rm(npRNA_LPA_rep1, npRNA_LPA_rep2)

colnames(npRNA_LPA_all) <- gsub("merged_bams_mm10_vM14.Nuc.LPA.rep..","", colnames(npRNA_LPA_all))
colnames(npRNA_LPA_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.LPA.rep..","", colnames(npRNA_LPA_all))
colnames(npRNA_LPA_all) <- gsub(".bam","", colnames(npRNA_LPA_all))

npRNA_LPA_rep1 <- read.table(file = 'merged_counts/Nuc.LPA.rep1/D40000/exons/exon_counts.txt.gz', header=T)
npRNA_LPA_rep2 <- read.table(file = 'merged_counts/Nuc.LPA.rep2/D40000/exons/exon_counts.txt.gz', header=T)
npRNA_LPA_exons <- data.frame(npRNA_LPA_rep1[,1:6],
                              npRNA_LPA_rep1[, -c(1:6)],
                              npRNA_LPA_rep2[, -c(1:6)])
rm(npRNA_LPA_rep1, npRNA_LPA_rep2)

colnames(npRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.Nuc.LPA.rep..","", colnames(npRNA_LPA_exons))
colnames(npRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14_2.Nuc.LPA.rep..","", colnames(npRNA_LPA_exons))
colnames(npRNA_LPA_exons) <- gsub(".bam","", colnames(npRNA_LPA_exons))

save(npRNA_LPA_all, npRNA_LPA_exons, file='raw_counts/npRNA_LPA.Rdata')

## load cytoplasmic data

cytoRNA_Naive_rep1 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep1/D40000/all/all_counts.txt.gz', header=T)
cytoRNA_Naive_rep2 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep2/D40000/all/all_counts.txt.gz', header=T)
cytoRNA_Naive_rep3 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep3/D40000/all/all_counts.txt.gz', header=T)
cytoRNA_Naive_all <- data.frame(cytoRNA_Naive_rep1[,1:6],
                              cytoRNA_Naive_rep1[, -c(1:6)],
                              cytoRNA_Naive_rep2[, -c(1:6)],
                              cytoRNA_Naive_rep3[, -c(1:6)])
rm(cytoRNA_Naive_rep1, cytoRNA_Naive_rep2, cytoRNA_Naive_rep3)

colnames(cytoRNA_Naive_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.Naive.rep..","", colnames(cytoRNA_Naive_all))
colnames(cytoRNA_Naive_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.Naive.rep..","", colnames(cytoRNA_Naive_all))
colnames(cytoRNA_Naive_all) <- gsub(".bam","", colnames(cytoRNA_Naive_all))

cytoRNA_Naive_rep1 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep1/D40000/exons/exon_counts.txt.gz', header=T)
cytoRNA_Naive_rep2 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep2/D40000/exons/exon_counts.txt.gz', header=T)
cytoRNA_Naive_rep3 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep3/D40000/exons/exon_counts.txt.gz', header=T)
cytoRNA_Naive_exons <- data.frame(cytoRNA_Naive_rep1[,1:6],
                                cytoRNA_Naive_rep1[, -c(1:6)],
                                cytoRNA_Naive_rep2[, -c(1:6)],
                                cytoRNA_Naive_rep3[, -c(1:6)])
rm(cytoRNA_Naive_rep1, cytoRNA_Naive_rep2, cytoRNA_Naive_rep3)

colnames(cytoRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14.cytoplasmic.Naive.rep..","", colnames(cytoRNA_Naive_exons))
colnames(cytoRNA_Naive_exons) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.Naive.rep..","", colnames(cytoRNA_Naive_exons))
colnames(cytoRNA_Naive_exons) <- gsub(".bam","", colnames(cytoRNA_Naive_exons))

save(cytoRNA_Naive_all, cytoRNA_Naive_exons, file='raw_counts/cytoRNA_Naive.Rdata')


cytoRNA_LPA_rep1 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep1/D40000/all/all_counts.txt.gz', header=T)
cytoRNA_LPA_rep2 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep2/D40000/all/all_counts.txt.gz', header=T)
cytoRNA_LPA_rep3 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep3/D40000/all/all_counts.txt.gz', header=T)
cytoRNA_LPA_all <- data.frame(cytoRNA_LPA_rep1[,1:6],
                            cytoRNA_LPA_rep1[, -c(1:6)],
                            cytoRNA_LPA_rep2[, -c(1:6)],
                            cytoRNA_LPA_rep3[, -c(1:6)])
rm(cytoRNA_LPA_rep1, cytoRNA_LPA_rep2, cytoRNA_LPA_rep3)

colnames(cytoRNA_LPA_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.LPA.rep..","", colnames(cytoRNA_LPA_all))
colnames(cytoRNA_LPA_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.LPA.rep..","", colnames(cytoRNA_LPA_all))
colnames(cytoRNA_LPA_all) <- gsub(".bam","", colnames(cytoRNA_LPA_all))

cytoRNA_LPA_rep1 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep1/D40000/exons/exon_counts.txt.gz', header=T)
cytoRNA_LPA_rep2 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep2/D40000/exons/exon_counts.txt.gz', header=T)
cytoRNA_LPA_rep3 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep3/D40000/exons/exon_counts.txt.gz', header=T)
cytoRNA_LPA_exons <- data.frame(cytoRNA_LPA_rep1[,1:6],
                              cytoRNA_LPA_rep1[, -c(1:6)],
                              cytoRNA_LPA_rep2[, -c(1:6)],
                              cytoRNA_LPA_rep3[, -c(1:6)])
rm(cytoRNA_LPA_rep1, cytoRNA_LPA_rep2, cytoRNA_LPA_rep3)

colnames(cytoRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14.cytoplasmic.LPA.rep..","", colnames(cytoRNA_LPA_exons))
colnames(cytoRNA_LPA_exons) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.LPA.rep..","", colnames(cytoRNA_LPA_exons))
colnames(cytoRNA_LPA_exons) <- gsub(".bam","", colnames(cytoRNA_LPA_exons))

save(cytoRNA_LPA_all, cytoRNA_LPA_exons, file='raw_counts/cytoRNA_LPA.Rdata')

## Create exons rpkms
load(file='raw_counts/caRNA_Naive.Rdata')
Time_ca_Naive <- as.numeric(sapply(colnames(caRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_Naive <- sapply(colnames(caRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_Naive <- as.factor(Time_ca_Naive)

load(file='raw_counts/npRNA_Naive.Rdata')
Time_np_Naive <- as.numeric(sapply(colnames(npRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_Naive <- sapply(colnames(npRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_Naive <- as.factor(Time_np_Naive)

load(file='raw_counts/cytoRNA_Naive.Rdata')
Time_cyto_Naive <- as.numeric(sapply(colnames(cytoRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_Naive <- sapply(colnames(cytoRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_Naive <- as.factor(Time_cyto_Naive)

library(edgeR)
caRNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6])), log=T)
npRNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6])), log=T)
cytoRNA_Naive_exons_cpm <- cpm(calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6])), log=T)

row.names(caRNA_Naive_exons_cpm) <- row.names(npRNA_Naive_exons_cpm) <- row.names(cytoRNA_Naive_exons_cpm) <- caRNA_Naive_exons$Geneid

save(caRNA_Naive_exons_cpm, file='exons_cpms/caRNA_Naive_exons_cpm.Rdata')
save(npRNA_Naive_exons_cpm, file='exons_cpms/npRNA_Naive_exons_cpm.Rdata')
save(cytoRNA_Naive_exons_cpm, file='exons_cpms/cytoRNA_Naive_exons_cpm.Rdata')

caRNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6])), log=T)
npRNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6])), log=T)
cytoRNA_Naive_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6])), log=T)

row.names(caRNA_Naive_exons_rpkm) <- row.names(npRNA_Naive_exons_rpkm) <- row.names(cytoRNA_Naive_exons_rpkm) <- caRNA_Naive_exons$Geneid

save(caRNA_Naive_exons_rpkm, file='exons_rpkms/caRNA_Naive_exons_rpkm.Rdata')
save(npRNA_Naive_exons_rpkm, file='exons_rpkms/npRNA_Naive_exons_rpkm.Rdata')
save(cytoRNA_Naive_exons_rpkm, file='exons_rpkms/cytoRNA_Naive_exons_rpkm.Rdata')

# get lib size and scaling factor
tmp <- calcNormFactors(DGEList(counts = caRNA_Naive_all[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_all[,1:6]))
lib_all <- data.frame( samples = row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors)
tmp <- calcNormFactors(DGEList(counts = npRNA_Naive_all[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_Naive_all[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

tmp <- calcNormFactors(DGEList(counts = caRNA_Naive_exons[,-c(1:6)], group = group_ca_Naive, genes = caRNA_Naive_exons[,1:6]))
lib_exons <- data.frame( samples = row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors)
tmp <- calcNormFactors(DGEList(counts = npRNA_Naive_exons[,-c(1:6)], group = group_np_Naive, genes = npRNA_Naive_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_Naive_exons[,-c(1:6)], group = group_cyto_Naive, genes = cytoRNA_Naive_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

load(file='raw_counts/caRNA_LPA.Rdata')
Time_ca_LPA <- as.numeric(sapply(colnames(caRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_LPA <- sapply(colnames(caRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_LPA <- as.factor(Time_ca_LPA)

load(file='raw_counts/npRNA_LPA.Rdata')
Time_np_LPA <- as.numeric(sapply(colnames(npRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_LPA <- sapply(colnames(npRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_LPA <- as.factor(Time_np_LPA)

load(file='raw_counts/cytoRNA_LPA.Rdata')
Time_cyto_LPA <- as.numeric(sapply(colnames(cytoRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_LPA <- sapply(colnames(cytoRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_LPA <- as.factor(Time_cyto_LPA)

caRNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6])), log=T)
npRNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6])), log=T)
cytoRNA_LPA_exons_cpm <- cpm(calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6])), log=T)

row.names(caRNA_LPA_exons_cpm) <- row.names(npRNA_LPA_exons_cpm) <- row.names(cytoRNA_LPA_exons_cpm) <- caRNA_LPA_exons$Geneid

save(caRNA_LPA_exons_cpm, file='exons_cpms/caRNA_LPA_exons_cpm.Rdata')
save(npRNA_LPA_exons_cpm, file='exons_cpms/npRNA_LPA_exons_cpm.Rdata')
save(cytoRNA_LPA_exons_cpm, file='exons_cpms/cytoRNA_LPA_exons_cpm.Rdata')

caRNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6])), log=T)
npRNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6])), log=T)
cytoRNA_LPA_exons_rpkm <- rpkm(calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6])), log=T)

row.names(caRNA_LPA_exons_rpkm) <- row.names(npRNA_LPA_exons_rpkm) <- row.names(cytoRNA_LPA_exons_rpkm) <- caRNA_LPA_exons$Geneid

save(caRNA_LPA_exons_rpkm, file='exons_rpkms/caRNA_LPA_exons_rpkm.Rdata')
save(npRNA_LPA_exons_rpkm, file='exons_rpkms/npRNA_LPA_exons_rpkm.Rdata')
save(cytoRNA_LPA_exons_rpkm, file='exons_rpkms/cytoRNA_LPA_exons_rpkm.Rdata')

tmp <- calcNormFactors(DGEList(counts = caRNA_LPA_all[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = npRNA_LPA_all[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_LPA_all[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_all[,1:6]))
lib_all <- rbind(lib_all,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

tmp <- calcNormFactors(DGEList(counts = caRNA_LPA_exons[,-c(1:6)], group = group_ca_LPA, genes = caRNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = npRNA_LPA_exons[,-c(1:6)], group = group_np_LPA, genes = npRNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples=row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))
tmp <- calcNormFactors(DGEList(counts = cytoRNA_LPA_exons[,-c(1:6)], group = group_cyto_LPA, genes = cytoRNA_LPA_exons[,1:6]))
lib_exons <- rbind(lib_exons,data.frame(samples =row.names(tmp$samples), lib.size = tmp$samples$lib.size, norm.factors = tmp$samples$norm.factors))

save(lib_all, lib_exons, file='raw_counts/lib.size.Rdata')

