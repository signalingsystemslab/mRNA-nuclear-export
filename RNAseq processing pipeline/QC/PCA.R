## Load packages
library(edgeR)
library(ggplot2)

## Color scheme
# col1 <- data.frame(red = 0, green = 60/255, blue = 1)
# col2 <- data.frame(red = 0, green = 1, blue = 160/255)
col3 <- data.frame(red = 1, green = 165/255, blue = 0)
col4 <- data.frame(red = 1, green = 0, blue = 137/255)

colors_scheme <- c(
                   rgb(red = col3$red + 0.5*(1-col3$red), blue = col3$blue + 0.5*(1-col3$blue), green = col3$green + 0.5*(1-col3$green)),
                   rgb(red = 0.75*col3$red, green = 0.75*col3$green, blue = 0.75*col3$blue),
                   rgb(red = col3$red + 0.25*(1-col3$red), blue = col3$blue + 0.25*(1-col3$blue), green = col3$green),
                   # rgb(red = 0.5*col3$red, blue = 0.5*col3$blue, green = 0.5*col3$green), 
                   rgb(red = col4$red + 0.5*(1-col4$red), blue = col4$blue + 0.5*(1-col4$blue), green = col4$green + 0.5*(1-col4$green)),
                   rgb(red = 0.75*col4$red, green = 0.75*col4$green, blue = 0.75*col4$blue),
                   rgb(red = col4$red + 0.25*(1-col4$red), blue = col4$blue + 0.25*(1-col4$blue), green = col4$green)
                   )

interactions <- interaction(rep(c("Naive", "LPA"), each=3), rep(c("Rep1","Rep2", "Rep3"),2))

## Get Read Counts For all lanes 
plot_pca_many_exp <- function(data, barcodes, type, interactive=F, high_threshold = 5, topnumber=200, ...){
  lane_barcode <- unlist(lapply(strsplit(colnames(data[,-c(1:6)]),split = '.', fixed = T), function(x){x[[3]]}))
  lane <- unlist(lapply(strsplit(lane_barcode,split = '_'), function(x){x[[1]]}))
  barcode <- unlist(lapply(strsplit(lane_barcode,split = '_'), function(x){x[[2]]}))
  
  data_dglist <- DGEList(counts = data[, -c(1:6)][,match(paste(barcodes[,5], barcodes[,6], sep='_'),lane_barcode)], genes = data[, 1:6])
  non_zeros_indexes <- which(rowSums(cpm(data_dglist) > 0) >= 1)
  
  # all
  names <- data_dglist$genes$Geneid[non_zeros_indexes]
  all_cpm <- cpm(y = calcNormFactors( data_dglist[non_zeros_indexes, ]), log = T)
  all_rpkm <- rpkm(y = calcNormFactors( data_dglist[non_zeros_indexes, ]), log = T)
  
  lane_barcode <- unlist(lapply(strsplit(colnames(data_dglist$counts),split = '.', fixed = T), function(x){x[[3]]}))
  lane <- unlist(lapply(strsplit(lane_barcode,split = '_'), function(x){x[[1]]}))
  barcode <- unlist(lapply(strsplit(lane_barcode,split = '_'), function(x){x[[2]]}))
  
  if(type == 'all'){
    dat_pca <- t(scale(t(all_rpkm)))
    
  }else if (type == 'high'){
    names_filtered <- names[apply(all_cpm, 1, function(x){any(x > high_threshold)})]
    all_rpkm_filtered <- all_rpkm[apply(all_cpm, 1, function(x){any(x > high_threshold)}), ] 
    
    dat_pca <- t(scale(t(all_rpkm_filtered)))
    
  }else if (type == 'top'){
    names_filtered <- names[apply(all_cpm, 1, function(x){any(x > high_threshold)})]
    all_rpkm_filtered <-all_rpkm[apply(all_cpm, 1, function(x){any(x > high_threshold)}), ]
    
    var <- apply(all_rpkm_filtered, 1, function(x){var(x)})
    top_rpkm <- all_rpkm_filtered[order(var, decreasing = T)[1:topnumber],]
    
    dat_pca <- t(scale(t(top_rpkm)))
  }else if ( type == 'lpa'){
    lpa_genes_rpkm <- all_rpkm[match(LPA_induced_genes$GeneId_vM14, as.character(names)),][!is.na(match(LPA_induced_genes$GeneId_vM14, as.character(names))),]
    
    dat_pca <- t(scale(t(lpa_genes_rpkm)))
  }
  
  pca <- prcomp(t(dat_pca), scale.=F)
  
  explained_var <- pca$sdev^2/sum(pca$sdev^2) * 100
  dat <- data.frame(pca$x, Lane = lane, Condition = barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),2], Time = as.factor(barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),3]), Batch = as.factor(barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),4])) 
  
  
  # p <- ggplot(dat, aes(x=PC1, y=PC2, col=Condition, size=Time, shape=Batch, group=interaction(Condition, Time, Lane))) 
  p <- ggplot(dat, aes(x=PC1, y=PC2, col=interaction(Condition, Batch), group=interaction(Condition, Time, Lane))) 
  # p <- p + geom_point(alpha = 0.3) 
  p <- p + geom_text(aes(label=Time)) + scale_color_manual(breaks = interactions[interactions %in% levels(interaction(dat$Condition, dat$Batch, drop=T))], values = colors_scheme[match(interactions[interactions %in% levels(interaction(dat$Condition, dat$Batch, drop=T))], as.character(interactions))])
  p <- p + theme_bw() + labs(x=paste0('PC1 (', signif(explained_var[1],4),'%)'), y=paste0('PC2 (', signif(explained_var[2],4),'%)'))
  return(p)
}

# chromatin rep1 #####
Xap114L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap114L1/08_counts/all/all_counts.txt.gz', header = T)
Xap114L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap114L2/08_counts/all/all_counts.txt.gz', header = T)
Xbp051L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXbp051L1/08_counts/all/all_counts.txt.gz', header = T)
WAP033L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsWAP033L1/08_counts/all/all_counts.txt.gz', header = T)
WAP033L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsWAP033L2/08_counts/all/all_counts.txt.gz', header = T)
WAP033L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsWAP033L3/08_counts/all/all_counts.txt.gz', header = T)
WAP033L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsWAP033L4/08_counts/all/all_counts.txt.gz', header = T)
WAP033L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsWAP033L5/08_counts/all/all_counts.txt.gz', header = T)
WAP033L6 <- read.table( file = '../../Data/individual_counts/SxaQSEQsWAP033L6/08_counts/all/all_counts.txt.gz', header = T)

chromatin_rep1_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/chromatin_rep1_barcodes.txt')

chromatin_rep1_LPA_Naive <- data.frame(Xap114L1[, 1:6], # Gene descriptions
                                       Xap114L1[, -c(1:6)], 
                                       Xap114L2[, -c(1:6)], 
                                       Xbp051L1[, -c(1:6)], 
                                       WAP033L1[, -c(1:6)], 
                                       WAP033L2[, -c(1:6)], 
                                       WAP033L3[, -c(1:6)], 
                                       WAP033L4[, -c(1:6)], 
                                       WAP033L5[, -c(1:6)], 
                                       WAP033L6[, -c(1:6)])

rm(Xap114L1,Xap114L2,Xbp051L1,WAP033L1,WAP033L2,WAP033L3, WAP033L4, WAP033L5, WAP033L6)

# chromatin rep2 #### 
YAP006L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP006L5/08_counts/all/all_counts.txt.gz', header = T)
YAP006L6 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP006L6/08_counts/all/all_counts.txt.gz', header = T)
YAP012L6 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L6/08_counts/all/all_counts.txt.gz', header = T)
YAP012L7 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L7/08_counts/all/all_counts.txt.gz', header = T)
YAP012L8 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L8/08_counts/all/all_counts.txt.gz', header = T)
YBP005L6 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP005L6/08_counts/all/all_counts.txt.gz', header = T)

chromatin_rep2_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/chromatin_rep2_barcodes.txt')

chromatin_rep2_LPA_Naive <- data.frame(YAP006L5[, 1:6], # Gene descriptions
                                       YAP006L5[, -c(1:6)], 
                                       YAP006L6[, -c(1:6)],
                                       YAP012L6[, -c(1:6)], 
                                       YAP012L7[, -c(1:6)], 
                                       YAP012L8[, -c(1:6)], 
                                       YBP005L6[, -c(1:6)])

rm(YAP006L5, YAP006L6, YAP012L6, YAP012L7, YAP012L8, YBP005L6)

# chromatin rep3 ####
YAP006L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP006L3/08_counts/all/all_counts.txt.gz', header = T)
YAP006L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP006L4/08_counts/all/all_counts.txt.gz', header = T)
YAP012L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L4/08_counts/all/all_counts.txt.gz', header = T)
YAP012L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L5/08_counts/all/all_counts.txt.gz', header = T)
YBP005L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP005L4/08_counts/all/all_counts.txt.gz', header = T)
YBP005L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP005L5/08_counts/all/all_counts.txt.gz', header = T)

chromatin_rep3_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/chromatin_rep3_barcodes.txt')

chromatin_rep3_LPA_Naive <- data.frame(YAP006L3[, 1:6], # Gene descriptions
                                       YAP006L3[, -c(1:6)], 
                                       YAP006L4[, -c(1:6)], 
                                       YAP012L4[, -c(1:6)], 
                                       YAP012L5[, -c(1:6)], 
                                       YBP005L4[, -c(1:6)], 
                                       YBP005L5[, -c(1:6)])

rm(YAP006L3, YAP006L4, YAP012L4, YAP012L5, YBP005L4, YBP005L5)

## all chromatin ####
chromatin_all_LPA_Naive <- data.frame(chromatin_rep1_LPA_Naive[, 1:6], 
                            chromatin_rep1_LPA_Naive[, -c(1:6)], 
                            chromatin_rep2_LPA_Naive[, -c(1:6)],
                            chromatin_rep3_LPA_Naive[, -c(1:6)])

chromatin_all_LPA_Naive_barcodes <- rbind(chromatin_rep1_LPA_Naive_barcodes, chromatin_rep2_LPA_Naive_barcodes, chromatin_rep3_LPA_Naive_barcodes)

p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='all')
p + labs(title='Chromatin_LPA_Naive (All)')

p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='high')
p + labs(title='Chromatin_LPA_Naive (Highly expressed genes)')

p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='top')
p + labs(title='Chromatin_LPA_Naive (Highly expressed genes and high variance)')

# nucleoplasmic Rep1 ####
Xap097L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap097L2/08_counts/all/all_counts.txt.gz', header = T)
Xap103L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap103L1/08_counts/all/all_counts.txt.gz', header = T)
Xap103L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap103L2/08_counts/all/all_counts.txt.gz', header = T)
Xap115L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap115L1/08_counts/all/all_counts.txt.gz', header = T)
Xap115L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap115L2/08_counts/all/all_counts.txt.gz', header = T)
VBP026L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L1/08_counts/all/all_counts.txt.gz', header = T)
VBP026L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L4/08_counts/all/all_counts.txt.gz', header = T)
VBP026L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L5/08_counts/all/all_counts.txt.gz', header = T)
VBP026L6 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L6/08_counts/all/all_counts.txt.gz', header = T)
VBP026L7 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L7/08_counts/all/all_counts.txt.gz', header = T)
YBP008L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP008L3/08_counts/all/all_counts.txt.gz', header = T)
YBP030L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP030L5/08_counts/all/all_counts.txt.gz', header = T)

nucleoplasmic_rep1_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/nucleoplasm_rep1_barcodes.txt')

nucleoplasmic_rep1_LPA_Naive <- data.frame(Xap097L2[, 1:6], # Gene descriptions
                                         Xap097L2[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xap097L2",6], function(x) {grep(x, names(Xap097L2))})], 
                                         Xap103L1[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xap103L1",6], function(x) {grep(x, names(Xap103L1))})], 
                                         Xap103L2[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xap103L2",6], function(x) {grep(x, names(Xap103L2))})], 
                                         Xap115L1[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xap115L1",6], function(x) {grep(x, names(Xap115L1))})], 
                                         Xap115L2[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xap115L2",6], function(x) {grep(x, names(Xap115L2))})], 
                                         VBP026L1[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L1",6], function(x) {grep(x, names(VBP026L1))})], 
                                         VBP026L4[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L4",6], function(x) {grep(x, names(VBP026L4))})], 
                                         VBP026L5[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L5",6], function(x) {grep(x, names(VBP026L5))})], 
                                         VBP026L6[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L6",6], function(x) {grep(x, names(VBP026L6))})], 
                                         VBP026L7[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L7",6], function(x) {grep(x, names(VBP026L7))})], 
                                         YBP008L3[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "YBP008L3",6], function(x) {grep(x, names(YBP008L3))})],
                                         YBP030L5[, sapply(nucleoplasmic_rep1_LPA_Naive_barcodes[nucleoplasmic_rep1_LPA_Naive_barcodes[,5] == "YBP030L5",6], function(x) {grep(x, names(YBP030L5))})])

rm(Xap097L2,Xap103L1,Xap103L2, Xap115L1,Xap115L2,VBP026L1, VBP026L4, VBP026L5, VBP026L6, VBP026L7, YBP008L3, YBP030L5)
colnames(nucleoplasmic_rep1_LPA_Naive) <- gsub(pattern = ".3.bam", replacement= ".bam", colnames(nucleoplasmic_rep1_LPA_Naive))
colnames(nucleoplasmic_rep1_LPA_Naive) <- gsub(pattern = "YBP008L3.", replacement= "YBP008L3_", colnames(nucleoplasmic_rep1_LPA_Naive))
colnames(nucleoplasmic_rep1_LPA_Naive) <- gsub(pattern = "YBP030L5.", replacement= "YBP030L5_", colnames(nucleoplasmic_rep1_LPA_Naive))


# nucleoplasmic rep2 ####
YAP006L8 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP006L8/08_counts/all/all_counts.txt.gz', header = T)
YBP005L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP005L3/08_counts/all/all_counts.txt.gz', header = T)
YAP012L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L1/08_counts/all/all_counts.txt.gz', header = T)
YAP012L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L2/08_counts/all/all_counts.txt.gz', header = T)
YAP012L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP012L3/08_counts/all/all_counts.txt.gz', header = T)

nucleoplasmic_rep2_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/nucleoplasm_rep2_barcodes.txt')

nucleoplasmic_rep2_LPA_Naive <- data.frame(YAP006L8[, 1:6], # Gene descriptions
                                         YAP006L8[, sapply(nucleoplasmic_rep2_LPA_Naive_barcodes[nucleoplasmic_rep2_LPA_Naive_barcodes[,5] == "YAP006L8",6], function(x) {grep(x, names(YAP006L8))})], 
                                         YBP005L3[, sapply(nucleoplasmic_rep2_LPA_Naive_barcodes[nucleoplasmic_rep2_LPA_Naive_barcodes[,5] == "YBP005L3",6], function(x) {grep(x, names(YBP005L3))})], 
                                         YAP012L1[, sapply(nucleoplasmic_rep2_LPA_Naive_barcodes[nucleoplasmic_rep2_LPA_Naive_barcodes[,5] == "YAP012L1",6], function(x) {grep(x, names(YAP012L1))})], 
                                         YAP012L2[, sapply(nucleoplasmic_rep2_LPA_Naive_barcodes[nucleoplasmic_rep2_LPA_Naive_barcodes[,5] == "YAP012L2",6], function(x) {grep(x, names(YAP012L2))})], 
                                         YAP012L3[, sapply(nucleoplasmic_rep2_LPA_Naive_barcodes[nucleoplasmic_rep2_LPA_Naive_barcodes[,5] == "YAP012L3",6], function(x) {grep(x, names(YAP012L3))})])

rm(YAP006L8, YBP005L3, YAP012L1, YAP012L2, YAP012L3)

## all nucleoplasmic #####
nucleoplasmic_all_LPA_Naive <- data.frame(nucleoplasmic_rep1_LPA_Naive[, 1:6], 
                                          nucleoplasmic_rep1_LPA_Naive[, -c(1:6)], 
                                          nucleoplasmic_rep2_LPA_Naive[, -c(1:6)])

nucleoplasmic_all_LPA_Naive_barcodes <- rbind(nucleoplasmic_rep1_LPA_Naive_barcodes, nucleoplasmic_rep2_LPA_Naive_barcodes)

p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='all')
p + labs(title='nucleoplasmic_LPA_Naive (All)') # few genes highly expressed Naive.Rep2 at 25minutes but not in other samples

p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='high')
p + labs(title='nucleoplasmic_LPA_Naive (Highly expressed genes)')

p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='top')
p + labs(title='nucleoplasmic_LPA_Naive (Highly expressed genes and high variance)')


# cytoplasmic rep1 #####
Xap103L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap103L1/08_counts/all/all_counts.txt.gz', header = T)
VBP026L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L2/08_counts/all/all_counts.txt.gz', header = T)
VBP026L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L3/08_counts/all/all_counts.txt.gz', header = T)
Xbp081L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXbp081L1/08_counts/all/all_counts.txt.gz', header = T)
Xap097L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap097L2/08_counts/all/all_counts.txt.gz', header = T)
Xap115L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap115L1/08_counts/all/all_counts.txt.gz', header = T)
Xap115L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsXap115L2/08_counts/all/all_counts.txt.gz', header = T)
VBP026L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsVBP026L1/08_counts/all/all_counts.txt.gz', header = T)
YBP008L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP008L3/08_counts/all/all_counts.txt.gz', header = T)
YBP030L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP030L4/08_counts/all/all_counts.txt.gz', header = T)

cytoplasmic_rep1_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/cytoplasmic_rep1_barcodes.txt')

cytoplasmic_rep1_LPA_Naive <- data.frame(Xap103L1[, 1:6], # Gene descriptions
                                      Xap103L1[, sapply(cytoplasmic_rep1_LPA_Naive_barcodes[cytoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xap103L1",6], function(x) {grep(x, names(Xap103L1))})], 
                                      VBP026L2[, sapply(cytoplasmic_rep1_LPA_Naive_barcodes[cytoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L2",6], function(x) {grep(x, names(VBP026L2))})], 
                                      VBP026L3[, sapply(cytoplasmic_rep1_LPA_Naive_barcodes[cytoplasmic_rep1_LPA_Naive_barcodes[,5] == "VBP026L3",6], function(x) {grep(x, names(VBP026L3))})], 
                                      Xbp081L1[, sapply(cytoplasmic_rep1_LPA_Naive_barcodes[cytoplasmic_rep1_LPA_Naive_barcodes[,5] == "Xbp081L1",6], function(x) {grep(x, names(Xbp081L1))})], 
                                      Xap097L2[, -c(1:6)], 
                                      Xap115L1[, -c(1:6)], 
                                      Xap115L2[, -c(1:6)], 
                                      VBP026L1[, -c(1:6)], 
                                      YBP008L3[, -c(1:6)],
                                      YBP030L4[, -c(1:6)])

rm(Xap103L1,VBP026L2,VBP026L3,Xbp081L1,Xap097L2,Xap115L1, Xap115L2, VBP026L1, YBP008L3, YBP030L4)

colnames(cytoplasmic_rep1_LPA_Naive) <- gsub('YBP008L3.', 'YBP008L3_', colnames(cytoplasmic_rep1_LPA_Naive))
colnames(cytoplasmic_rep1_LPA_Naive) <- gsub('YBP030L4.', 'YBP030L4_', colnames(cytoplasmic_rep1_LPA_Naive))

# cytoplasmic rep2
YAP015L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L1/08_counts/all/all_counts.txt.gz', header = T)
YAP015L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L2/08_counts/all/all_counts.txt.gz', header = T)
YAP015L3 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L3/08_counts/all/all_counts.txt.gz', header = T)
YAP015L4 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L4/08_counts/all/all_counts.txt.gz', header = T)
YBP005L2 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP005L2/08_counts/all/all_counts.txt.gz', header = T)

cytoplasmic_rep2_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/cytoplasmic_rep2_barcodes.txt')

cytoplasmic_rep2_LPA_Naive <- data.frame(YAP015L1[, 1:6], # Gene descriptions
                                       YAP015L1[, -c(1:6)], 
                                       YAP015L2[, -c(1:6)], 
                                       YAP015L3[, -c(1:6)], 
                                       YAP015L4[, -c(1:6)], 
                                       YBP005L2[, -c(1:6)])

rm(YAP015L1,YAP015L2,YAP015L3,YAP015L4,YBP005L2)

# cytoplasmic rep3
YAP015L5 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L5/08_counts/all/all_counts.txt.gz', header = T)
YAP015L6 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L6/08_counts/all/all_counts.txt.gz', header = T)
YAP015L7 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L7/08_counts/all/all_counts.txt.gz', header = T)
YAP015L8 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYAP015L8/08_counts/all/all_counts.txt.gz', header = T)
YBP005L1 <- read.table( file = '../../Data/individual_counts/SxaQSEQsYBP005L1/08_counts/all/all_counts.txt.gz', header = T)

cytoplasmic_rep3_LPA_Naive_barcodes <- read.table(file = '../../Data/Barcodes/LPA_vs_Naive/cytoplasmic_rep3_barcodes.txt')

cytoplasmic_rep3_LPA_Naive <- data.frame(YAP015L5[, 1:6], # Gene descriptions
                                       YAP015L5[, -c(1:6)], 
                                       YAP015L6[, -c(1:6)], 
                                       YAP015L7[, -c(1:6)], 
                                       YAP015L8[, -c(1:6)], 
                                       YBP005L1[, -c(1:6)])

rm(YAP015L5,YAP015L6,YAP015L7,YAP015L8,YBP005L1)

## all cytoplamic #####
cytoplasmic_all_LPA_Naive <- data.frame(cytoplasmic_rep1_LPA_Naive[, 1:6], 
                                  cytoplasmic_rep1_LPA_Naive[, -c(1:6)],
                                  cytoplasmic_rep2_LPA_Naive[,-c(1:6)],
                                  cytoplasmic_rep3_LPA_Naive[, -c(1:6)])

cytoplasmic_all_LPA_Naive_barcodes <- rbind(cytoplasmic_rep1_LPA_Naive_barcodes,cytoplasmic_rep2_LPA_Naive_barcodes, cytoplasmic_rep3_LPA_Naive_barcodes)

p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='all')
p + labs(title='Cytoplasmic_LPA_Naive (All)')

p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='high')
p + labs(title='Cytoplasmic_LPA_Naive (Highly expressed genes)')

p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='top')
p + labs(title='Cytoplasmic_LPA_Naive (Highly expressed genes and high variance)')
