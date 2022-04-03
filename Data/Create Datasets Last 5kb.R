
########################################################################################################################
# Create dataset with last 5kb exons
caRNA_Naive_exons_5kb_rep1 <- read.table(file = 'merged_counts/chromatin.Naive.rep1/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
caRNA_Naive_exons_5kb_rep2 <- read.table(file = 'merged_counts/chromatin.Naive.rep2/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
caRNA_Naive_exons_5kb_rep3 <- read.table(file = 'merged_counts/chromatin.Naive.rep3/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)

caRNA_Naive_exons_5kb_all <- data.frame(caRNA_Naive_exons_5kb_rep1[,1:6],
                                        caRNA_Naive_exons_5kb_rep1[, -c(1:6)],
                                        caRNA_Naive_exons_5kb_rep2[, -c(1:6)],
                                        caRNA_Naive_exons_5kb_rep3[, -c(1:6)])
rm(caRNA_Naive_exons_5kb_rep1, caRNA_Naive_exons_5kb_rep2, caRNA_Naive_exons_5kb_rep3)

colnames(caRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.chromatin.Naive.rep..","", colnames(caRNA_Naive_exons_5kb_all))
colnames(caRNA_Naive_exons_5kb_all) <- gsub(".bam","", colnames(caRNA_Naive_exons_5kb_all))

save(caRNA_Naive_exons_5kb_all, file='raw_counts/caRNA_Naive_5kb.Rdata')

npRNA_Naive_exons_5kb_rep1 <- read.table(file = 'merged_counts/Nuc.Naive.rep1/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
npRNA_Naive_exons_5kb_rep2 <- read.table(file = 'merged_counts/Nuc.Naive.rep2/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)

npRNA_Naive_exons_5kb_all <- data.frame(npRNA_Naive_exons_5kb_rep1[,1:6],
                                        npRNA_Naive_exons_5kb_rep1[, -c(1:6)],
                                        npRNA_Naive_exons_5kb_rep2[, -c(1:6)])
rm(npRNA_Naive_exons_5kb_rep1, npRNA_Naive_exons_5kb_rep2)

colnames(npRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.Nuc.Naive.rep..","", colnames(npRNA_Naive_exons_5kb_all))
colnames(npRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.Naive.rep..","", colnames(npRNA_Naive_exons_5kb_all))
colnames(npRNA_Naive_exons_5kb_all) <- gsub(".bam","", colnames(npRNA_Naive_exons_5kb_all))

save(npRNA_Naive_exons_5kb_all, file='raw_counts/npRNA_Naive_5kb.Rdata')

cytoRNA_Naive_exons_5kb_rep1 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep1/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
cytoRNA_Naive_exons_5kb_rep2 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep2/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
cytoRNA_Naive_exons_5kb_rep3 <- read.table(file = 'merged_counts/cytoplasmic.Naive.rep3/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)

cytoRNA_Naive_exons_5kb_all <- data.frame(cytoRNA_Naive_exons_5kb_rep1[,1:6],
                                          cytoRNA_Naive_exons_5kb_rep1[, -c(1:6)],
                                          cytoRNA_Naive_exons_5kb_rep2[, -c(1:6)],
                                          cytoRNA_Naive_exons_5kb_rep3[, -c(1:6)])
rm(cytoRNA_Naive_exons_5kb_rep1, cytoRNA_Naive_exons_5kb_rep2, cytoRNA_Naive_exons_5kb_rep3)

colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.Naive.rep..","", colnames(cytoRNA_Naive_exons_5kb_all))
colnames(cytoRNA_Naive_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.Naive.rep..","", colnames(cytoRNA_Naive_exons_5kb_all))
colnames(cytoRNA_Naive_exons_5kb_all) <- gsub(".bam","", colnames(cytoRNA_Naive_exons_5kb_all))

save(cytoRNA_Naive_exons_5kb_all, file='raw_counts/cytoRNA_Naive_5kb.Rdata')

## Create exons rpkms (normalised using total exons library size)
load(file='raw_counts/caRNA_Naive_5kb.Rdata')
load(file='raw_counts/caRNA_Naive.Rdata')

Time_ca_Naive <- as.numeric(sapply(colnames(caRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_Naive <- sapply(colnames(caRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_Naive <- as.factor(Time_ca_Naive)

Time_ca_Naive_5kb <- as.numeric(sapply(colnames(caRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_Naive_5kb <- sapply(colnames(caRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_Naive_5kb <- as.factor(Time_ca_Naive_5kb)

load(file='raw_counts/npRNA_Naive_5kb.Rdata')
load(file='raw_counts/npRNA_Naive.Rdata')

Time_np_Naive <- as.numeric(sapply(colnames(npRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_Naive <- sapply(colnames(npRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_Naive <- as.factor(Time_np_Naive)

Time_np_Naive_5kb <- as.numeric(sapply(colnames(npRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_Naive_5kb <- sapply(colnames(npRNA_Naive_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_Naive_5kb <- as.factor(Time_np_Naive_5kb)

load(file='raw_counts/cytoRNA_Naive_5kb.Rdata')
load(file='raw_counts/cytoRNA_Naive.Rdata')

Time_cyto_Naive <- as.numeric(sapply(colnames(cytoRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_Naive <- sapply(colnames(cytoRNA_Naive_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_Naive <- as.factor(Time_cyto_Naive)

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

save(caRNA_Naive_exons_rpkm_5kb, file='exons_rpkms/caRNA_Naive_exons_rpkm_5kb.Rdata')
save(npRNA_Naive_exons_rpkm_5kb, file='exons_rpkms/npRNA_Naive_exons_rpkm_5kb.Rdata')
save(cytoRNA_Naive_exons_rpkm_5kb, file='exons_rpkms/cytoRNA_Naive_exons_rpkm_5kb.Rdata')

save(caRNA_Naive_exons_cpm_5kb, file='exons_cpms/caRNA_Naive_exons_cpm_5kb.Rdata')
save(npRNA_Naive_exons_cpm_5kb, file='exons_cpms/npRNA_Naive_exons_cpm_5kb.Rdata')
save(cytoRNA_Naive_exons_cpm_5kb, file='exons_cpms/cytoRNA_Naive_exons_cpm_5kb.Rdata')


## LPA Create dataset with last 5kb exons
caRNA_LPA_exons_5kb_rep1 <- read.table(file = 'merged_counts/chromatin.LPA.rep1/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
caRNA_LPA_exons_5kb_rep2 <- read.table(file = 'merged_counts/chromatin.LPA.rep2/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
caRNA_LPA_exons_5kb_rep3 <- read.table(file = 'merged_counts/chromatin.LPA.rep3/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
caRNA_LPA_exons_5kb_all <- data.frame(caRNA_LPA_exons_5kb_rep1[,1:6],
                                      caRNA_LPA_exons_5kb_rep1[, -c(1:6)],
                                      caRNA_LPA_exons_5kb_rep2[, -c(1:6)],
                                      caRNA_LPA_exons_5kb_rep3[, -c(1:6)])
rm(caRNA_LPA_exons_5kb_rep1, caRNA_LPA_exons_5kb_rep2, caRNA_LPA_exons_5kb_rep3)

colnames(caRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.chromatin.LPA.rep..","", colnames(caRNA_LPA_exons_5kb_all))
colnames(caRNA_LPA_exons_5kb_all) <- gsub(".bam","", colnames(caRNA_LPA_exons_5kb_all))

save(caRNA_LPA_exons_5kb_all, file='raw_counts/caRNA_LPA_5kb.Rdata')

npRNA_LPA_exons_5kb_rep1 <- read.table(file = 'merged_counts/Nuc.LPA.rep1/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
npRNA_LPA_exons_5kb_rep2 <- read.table(file = 'merged_counts/Nuc.LPA.rep2/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
npRNA_LPA_exons_5kb_all <- data.frame(npRNA_LPA_exons_5kb_rep1[,1:6],
                                      npRNA_LPA_exons_5kb_rep1[, -c(1:6)],
                                      npRNA_LPA_exons_5kb_rep2[, -c(1:6)])
rm(npRNA_LPA_exons_5kb_rep1, npRNA_LPA_exons_5kb_rep2)

colnames(npRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.Nuc.LPA.rep..","", colnames(npRNA_LPA_exons_5kb_all))
colnames(npRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.Nuc.LPA.rep..","", colnames(npRNA_LPA_exons_5kb_all))
colnames(npRNA_LPA_exons_5kb_all) <- gsub(".bam","", colnames(npRNA_LPA_exons_5kb_all))

save(npRNA_LPA_exons_5kb_all, file='raw_counts/npRNA_LPA_5kb.Rdata')

cytoRNA_LPA_exons_5kb_rep1 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep1/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
cytoRNA_LPA_exons_5kb_rep2 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep2/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
cytoRNA_LPA_exons_5kb_rep3 <- read.table(file = 'merged_counts/cytoplasmic.LPA.rep3/D40000_LPAgenes_last5kb/exon_counts.txt.gz', header=T)
cytoRNA_LPA_exons_5kb_all <- data.frame(cytoRNA_LPA_exons_5kb_rep1[,1:6],
                                        cytoRNA_LPA_exons_5kb_rep1[, -c(1:6)],
                                        cytoRNA_LPA_exons_5kb_rep2[, -c(1:6)],
                                        cytoRNA_LPA_exons_5kb_rep3[, -c(1:6)])
rm(cytoRNA_LPA_exons_5kb_rep1, cytoRNA_LPA_exons_5kb_rep2, cytoRNA_LPA_exons_5kb_rep3)

colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14.cytoplasmic.LPA.rep..","", colnames(cytoRNA_LPA_exons_5kb_all))
colnames(cytoRNA_LPA_exons_5kb_all) <- gsub("merged_bams_mm10_vM14_2.cytoplasmic.LPA.rep..","", colnames(cytoRNA_LPA_exons_5kb_all))
colnames(cytoRNA_LPA_exons_5kb_all) <- gsub(".bam","", colnames(cytoRNA_LPA_exons_5kb_all))

save(cytoRNA_LPA_exons_5kb_all, file='raw_counts/cytoRNA_LPA_5kb.Rdata')

## Create exons rpkms (normalised using total exons library size)
load(file='raw_counts/caRNA_LPA_5kb.Rdata')
load(file='raw_counts/caRNA_LPA.Rdata')

Time_ca_LPA <- as.numeric(sapply(colnames(caRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_LPA <- sapply(colnames(caRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_LPA <- as.factor(Time_ca_LPA)

Time_ca_LPA_5kb <- as.numeric(sapply(colnames(caRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_ca_LPA_5kb <- sapply(colnames(caRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_ca_LPA_5kb <- as.factor(Time_ca_LPA_5kb)

load(file='raw_counts/npRNA_LPA_5kb.Rdata')
load(file='raw_counts/npRNA_LPA.Rdata')

Time_np_LPA <- as.numeric(sapply(colnames(npRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_LPA <- sapply(colnames(npRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_LPA <- as.factor(Time_np_LPA)

Time_np_LPA_5kb <- as.numeric(sapply(colnames(npRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_np_LPA_5kb <- sapply(colnames(npRNA_LPA_exons_5kb_all)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_np_LPA_5kb <- as.factor(Time_np_LPA_5kb)

load(file='raw_counts/cytoRNA_LPA_5kb.Rdata')
load(file='raw_counts/cytoRNA_LPA.Rdata')

Time_cyto_LPA <- as.numeric(sapply(colnames(cytoRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][3]}))
Batch_cyto_LPA <- sapply(colnames(cytoRNA_LPA_exons)[-(1:6)], function(x){strsplit(x,fixed = T, split = ".")[[1]][4]})
group_cyto_LPA <- as.factor(Time_cyto_LPA)

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

save(caRNA_LPA_exons_rpkm_5kb, file='exons_rpkms/caRNA_LPA_exons_rpkm_5kb.Rdata')
save(npRNA_LPA_exons_rpkm_5kb, file='exons_rpkms/npRNA_LPA_exons_rpkm_5kb.Rdata')
save(cytoRNA_LPA_exons_rpkm_5kb, file='exons_rpkms/cytoRNA_LPA_exons_rpkm_5kb.Rdata')

save(caRNA_LPA_exons_cpm_5kb, file='exons_cpms/caRNA_LPA_exons_cpm_5kb.Rdata')
save(npRNA_LPA_exons_cpm_5kb, file='exons_cpms/npRNA_LPA_exons_cpm_5kb.Rdata')
save(cytoRNA_LPA_exons_cpm_5kb, file='exons_cpms/cytoRNA_LPA_exons_cpm_5kb.Rdata')

## Create gene_infos_with_length
library(biomaRt)
ensembl <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl', host = 'may2017.archive.ensembl.org')

gene_infos <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",'version'),
                 filters = "ensembl_gene_id",
                 values = unlist(lapply(strsplit(x = caRNA_Naive_exons_5kb_all$Geneid, split = ".", fixed = T),
                                        function(x){x[1]})),
                 mart = ensembl )

gene_infos$Length_5kb_naive = caRNA_Naive_exons_5kb_all$Length[match(gene_infos$ensembl_gene_id, sapply(caRNA_Naive_exons_5kb_all$Geneid, function(x){substr(x, 1, 18)}))]
gene_infos$Length_5kb_lpa = caRNA_LPA_exons_5kb_all$Length[match(gene_infos$ensembl_gene_id, sapply(caRNA_LPA_exons_5kb_all$Geneid, function(x){substr(x, 1, 18)}))]

save(gene_infos, file='gene_infos_with_length.Rdata')
