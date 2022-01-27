## Load packages
library(edgeR)
library(ggplot2)
library(plotly)
library(biomaRt)

setwd('/home/diane/Documents/Disk1/Documents/Ribo_Modeling_Reprocessed/Paper Final Run/QC')
## Color scheme
# col1 <- data.frame(red = 0, green = 60/255,blue = 1)
# col2 <- data.frame(red = 0, green = 1, blue = 160/255)
col3 <- data.frame(red = 1, green = 165/255, blue = 0)
col4 <- data.frame(red = 1, green = 0, blue = 137/255)

colors_scheme <- c(# rgb(red = col1$red + 0.5*(1-col1$red), blue = col1$blue + 0.5*(1-col1$blue), green = col1$green + 0.5*(1-col1$green)),
                   # rgb(red = col1$red + 0.25*(1-col1$red), blue = col1$blue + 0.25*(1-col1$blue), green = col1$green),
                   # rgb(red = 0.75*col1$red, green = 0.75*col1$green, blue = 0.75*col1$blue),
                   # rgb(red = 0.5*col1$red, blue = 0.5*col1$blue, green = 0.5*col1$green),
                   # rgb(red = col2$red + 0.5*(1-col2$red), blue = col2$blue + 0.5*(1-col2$blue), green = col2$green + 0.5*(1-col2$green)),
                   # rgb(red = col2$red + 0.25*(1-col2$red), blue = col2$blue + 0.25*(1-col2$blue), green = col2$green),
                   # rgb(red = 0.75*col2$red, green = 0.75*col2$green, blue = 0.75*col2$blue),
                   # rgb(red = 0.5*col2$red, blue = 0.5*col2$blue, green = 0.5*col2$green) # ,
                   rgb(red = col3$red + 0.5*(1-col3$red), blue = col3$blue + 0.5*(1-col3$blue), green = col3$green + 0.5*(1-col3$green)),
                   rgb(red = col3$red + 0.25*(1-col3$red), blue = col3$blue + 0.25*(1-col3$blue), green = col3$green),
                   rgb(red = 0.75*col3$red, green = 0.75*col3$green, blue = 0.75*col3$blue),
                   # rgb(red = 0.5*col3$red, blue = 0.5*col3$blue, green = 0.5*col3$green), 
                   rgb(red = col4$red + 0.5*(1-col4$red), blue = col4$blue + 0.5*(1-col4$blue), green = col4$green + 0.5*(1-col4$green)),
                   rgb(red = col4$red + 0.25*(1-col4$red), blue = col4$blue + 0.25*(1-col4$blue), green = col4$green),
                   rgb(red = 0.75*col4$red, green = 0.75*col4$green, blue = 0.75*col4$blue) # ,
                   # rgb(red = 0.5*col4$red, blue = 0.5*col4$blue, green = 0.5*col4$green)
)

interactions <- interaction(rep(c("Naive", "LPA"), each=3), rep(c("b1","b2", "b3"),2))

## Get Read Counts For all lanes 
plot_pca_few_exp <- function(data, barcodes, type, interactive=F, high_threshold = 5, topnumber=200, ...){
  #extract correct columns
  
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
  
  p <- ggplot(dat, aes(x=PC1, y=PC2, col=Time, size=Time, shape=Condition, group=interaction(Condition, Time, Lane))) + geom_point(alpha = 0.3) + theme_bw() + labs(x=paste0('PC1 (', signif(explained_var[1],4),'%)'), y=paste0('PC2 (', signif(explained_var[2],4),'%)'))
  return(p)
}
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
  p <- p + geom_text(aes(label=Time)) + scale_color_manual(breaks = interactions[interactions %in% levels(interaction(dat$Condition, dat$Batch, drop=T))], values = colors_scheme[match(levels(interaction(dat$Condition, dat$Batch, drop=T)), interactions)])
  p <- p + theme_bw() + labs(x=paste0('PC1 (', signif(explained_var[1],4),'%)'), y=paste0('PC2 (', signif(explained_var[2],4),'%)'))
  return(p)
}
plot_pca_many_exp_cpt <- function(data, barcodes, type, interactive=F, high_threshold = 5, topnumber=200, ...){
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
  dat <- data.frame(pca$x, Lane = lane, Compartment = barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),1], Condition = barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),2], Time = as.factor(barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),3]), Batch = as.factor(barcodes[match(lane_barcode, paste(barcodes[,5], barcodes[,6], sep='_')),4])) 
  
  p <- ggplot(dat, aes(x=PC1, y=PC2, col=Condition, size=Time, shape=Compartment, group=interaction(Condition, Time, Lane))) + geom_point(alpha = 0.3) + theme_bw() + labs(x=paste0('PC1 (', signif(explained_var[1],4),'%)'), y=paste0('PC2 (', signif(explained_var[2],4),'%)'))
  return(p)
}


## Conversion LPA induced genes to vM14
LPA_induced_genes <- read.csv( file = '../DEGs/top_table_summary.csv', header = T, as.is = T) # need to add new annotation
LPA_induced_genes <- LPA_induced_genes[LPA_induced_genes$Naive_non_bleeding,]
colnames(LPA_induced_genes)[1] <- "GeneId_vM14"

# listMarts(host='may2017.archive.ensembl.org') # Gencode VM14
# ensembl_mm_89 <- useMart(host='may2017.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
# attributes <- listAttributes(ensembl_mm_89)
# gene_attributes_mgi <- getBM(attributes = c('ensembl_gene_id', 'version', 'start_position', 'end_position', 'external_gene_name'), mart= ensembl_mm_89, filters = 'external_gene_name', values=LPA_induced_genes$Name)
#
# gene_attributes_mgi <- gene_attributes_mgi[!duplicated(gene_attributes_mgi$external_gene_name),]
#LPA_induced_genes$GeneId_vM14 <- paste0(gene_attributes_mgi$ensembl_gene_id, '.', gene_attributes_mgi$version)[match(LPA_induced_genes$Name, gene_attributes_mgi$external_gene_name)]

# chromatin b1 #####
Xap114L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap114L1/08_counts/all/all_counts.txt', header = T)
Xap114L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap114L2/08_counts/all/all_counts.txt', header = T)
Xbp051L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXbp051L1/08_counts/all/all_counts.txt', header = T)
WAP033L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsWAP033L1/08_counts/all/all_counts.txt', header = T)
WAP033L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsWAP033L2/08_counts/all/all_counts.txt', header = T)
WAP033L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsWAP033L3/08_counts/all/all_counts.txt', header = T)
WAP033L4 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsWAP033L4/08_counts/all/all_counts.txt', header = T)
WAP033L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsWAP033L5/08_counts/all/all_counts.txt', header = T)
WAP033L6 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsWAP033L6/08_counts/all/all_counts.txt', header = T)

chromatin_b1_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/chromatin_b1_barcodes.txt')

chromatin_b1_LPA_Naive <- data.frame(Xap114L1[, 1:6], # Gene descriptions
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

# chromatin b2 ####
YAP006L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L3/08_counts/all/all_counts.txt', header = T)
YAP006L4 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L4/08_counts/all/all_counts.txt', header = T)
YAP012L4 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L4/08_counts/all/all_counts.txt', header = T)
YAP012L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L5/08_counts/all/all_counts.txt', header = T)
YBP005L4 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L4/08_counts/all/all_counts.txt', header = T)
YBP005L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L5/08_counts/all/all_counts.txt', header = T)

chromatin_b2_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/chromatin_b2_barcodes.txt')

chromatin_b2_LPA_Naive <- data.frame(YAP006L3[, 1:6], # Gene descriptions
                           YAP006L3[, -c(1:6)], 
                           YAP006L4[, -c(1:6)], 
                           YAP012L4[, -c(1:6)], 
                           YAP012L5[, -c(1:6)], 
                           YBP005L4[, -c(1:6)], 
                           YBP005L5[, -c(1:6)])

rm(YAP006L3, YAP006L4, YAP012L4, YAP012L5, YBP005L4, YBP005L5)

# chromatin b3 #### 
YAP006L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L5/08_counts/all/all_counts.txt', header = T)
YAP006L6 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L6/08_counts/all/all_counts.txt', header = T)
YAP012L6 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L6/08_counts/all/all_counts.txt', header = T)
YAP012L7 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L7/08_counts/all/all_counts.txt', header = T)
YAP012L8 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L8/08_counts/all/all_counts.txt', header = T)
YBP005L6 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L6/08_counts/all/all_counts.txt', header = T)

chromatin_b3_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/chromatin_b3_barcodes.txt')

chromatin_b3_LPA_Naive <- data.frame(YAP006L5[, 1:6], # Gene descriptions
                           YAP006L5[, -c(1:6)], 
                           YAP006L6[, -c(1:6)],
                           YAP012L6[, -c(1:6)], 
                           YAP012L7[, -c(1:6)], 
                           YAP012L8[, -c(1:6)], 
                           YBP005L6[, -c(1:6)])

rm(YAP006L5, YAP006L6, YAP012L6, YAP012L7, YAP012L8, YBP005L6)

## all chromatin ####
chromatin_all_LPA_Naive <- data.frame(chromatin_b1_LPA_Naive[, 1:6], 
                            chromatin_b1_LPA_Naive[, -c(1:6)], 
                            chromatin_b2_LPA_Naive[, -c(1:6)],
                            chromatin_b3_LPA_Naive[, -c(1:6)])

chromatin_all_LPA_Naive_barcodes <- rbind(chromatin_b1_LPA_Naive_barcodes, chromatin_b2_LPA_Naive_barcodes, chromatin_b3_LPA_Naive_barcodes)

# p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='all')
# p + labs(title='Chromatin_LPA_Naive (All)')

# p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='high')
# p + labs(title='Chromatin_LPA_Naive (Highly expressed genes)')

# p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='top')
# p + labs(title='Chromatin_LPA_Naive (Highly expressed genes and high variance)')

p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes, type='lpa')
p + labs(title='Chromatin (LPA induced genes)') + theme(aspect.ratio=1)


p <- plot_pca_many_exp(chromatin_all_LPA_Naive, chromatin_all_LPA_Naive_barcodes[chromatin_all_LPA_Naive_barcodes$V2 == "Naive",], type='lpa')
p + labs(title='Chromatin (LPA induced genes)') + theme(aspect.ratio=1)


# nucleoplasmic b0 ####
# Xap103L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_Naive_RNAseq/SxaQSEQsXap103L2/08_counts/all/all_counts.txt', header = T)
# VBP026L4 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_Naive_RNAseq/SxaQSEQsVBP026L4/08_counts/all/all_counts.txt', header = T)
# VBP026L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_Naive_RNAseq/SxaQSEQsVBP026L5/08_counts/all/all_counts.txt', header = T)
# VBP026L6 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_Naive_RNAseq/SxaQSEQsVBP026L6/08_counts/all/all_counts.txt', header = T)
# VBP026L7 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_Naive_RNAseq/SxaQSEQsVBP026L7/08_counts/all/all_counts.txt', header = T)

# nucleoplasmic_b0_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/nucleoplasm_b0_barcodes.txt')

# nucleoplasmic_b0_LPA_Naive <- data.frame(Xap103L2[, 1:6], # Gene descriptions
#                                          Xap103L2[, sapply(nucleoplasmic_b0_LPA_Naive_barcodes[nucleoplasmic_b0_LPA_Naive_barcodes[,5] == "Xap103L2",6], function(x) {grep(x, names(Xap103L2))})], 
#                                          VBP026L4[, sapply(nucleoplasmic_b0_LPA_Naive_barcodes[nucleoplasmic_b0_LPA_Naive_barcodes[,5] == "VBP026L4",6], function(x) {grep(x, names(VBP026L4))})], 
#                                          VBP026L5[, sapply(nucleoplasmic_b0_LPA_Naive_barcodes[nucleoplasmic_b0_LPA_Naive_barcodes[,5] == "VBP026L5",6], function(x) {grep(x, names(VBP026L5))})], 
#                                          VBP026L6[, sapply(nucleoplasmic_b0_LPA_Naive_barcodes[nucleoplasmic_b0_LPA_Naive_barcodes[,5] == "VBP026L6",6], function(x) {grep(x, names(VBP026L6))})], 
#                                          VBP026L7[, sapply(nucleoplasmic_b0_LPA_Naive_barcodes[nucleoplasmic_b0_LPA_Naive_barcodes[,5] == "VBP026L7",6], function(x) {grep(x, names(VBP026L7))})])
# 
# rm(Xap103L2, VBP026L4, VBP026L5, VBP026L6, VBP026L7)
# 


# nucleoplasmic b1 ####
Xap097L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap097L2/08_counts/all/all_counts.txt', header = T)
Xap103L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap103L1/08_counts/all/all_counts.txt', header = T)
Xap103L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap103L2/08_counts/all/all_counts.txt', header = T)
Xap115L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap115L1/08_counts/all/all_counts.txt', header = T)
Xap115L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap115L2/08_counts/all/all_counts.txt', header = T)
VBP026L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L1/08_counts/all/all_counts.txt', header = T)
VBP026L4 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L4/08_counts/all/all_counts.txt', header = T)
VBP026L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L5/08_counts/all/all_counts.txt', header = T)
VBP026L6 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L6/08_counts/all/all_counts.txt', header = T)
VBP026L7 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L7/08_counts/all/all_counts.txt', header = T)
YBP008L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1_vs_M2/lanes/SxaQSEQsYBP008L3/08_counts/all/all_counts.txt', header = T)
YBP030L5 <- read.table( file = '~/Documents/Biggie_Expansion/Diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP030L5/08_counts/all/all_counts.txt', header = T)

nucleoplasmic_b1_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/nucleoplasm_b1_barcodes.txt')

nucleoplasmic_b1_LPA_Naive <- data.frame(Xap097L2[, 1:6], # Gene descriptions
                                         Xap097L2[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "Xap097L2",6], function(x) {grep(x, names(Xap097L2))})], 
                                         Xap103L1[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "Xap103L1",6], function(x) {grep(x, names(Xap103L1))})], 
                                         Xap103L2[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "Xap103L2",6], function(x) {grep(x, names(Xap103L2))})], 
                                         Xap115L1[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "Xap115L1",6], function(x) {grep(x, names(Xap115L1))})], 
                                         Xap115L2[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "Xap115L2",6], function(x) {grep(x, names(Xap115L2))})], 
                                         VBP026L1[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L1",6], function(x) {grep(x, names(VBP026L1))})], 
                                         VBP026L4[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L4",6], function(x) {grep(x, names(VBP026L4))})], 
                                         VBP026L5[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L5",6], function(x) {grep(x, names(VBP026L5))})], 
                                         VBP026L6[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L6",6], function(x) {grep(x, names(VBP026L6))})], 
                                         VBP026L7[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L7",6], function(x) {grep(x, names(VBP026L7))})], 
                                         YBP008L3[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "YBP008L3",6], function(x) {grep(x, names(YBP008L3))})],
                                         YBP030L5[, sapply(nucleoplasmic_b1_LPA_Naive_barcodes[nucleoplasmic_b1_LPA_Naive_barcodes[,5] == "YBP030L5",6], function(x) {grep(x, names(YBP030L5))})])

rm(Xap097L2,Xap103L1,Xap103L2, Xap115L1,Xap115L2,VBP026L1, VBP026L4, VBP026L5, VBP026L6, VBP026L7, YBP008L3, YBP030L5)
colnames(nucleoplasmic_b1_LPA_Naive) <- gsub(pattern = ".3.bam", replacement= ".bam", colnames(nucleoplasmic_b1_LPA_Naive))
colnames(nucleoplasmic_b1_LPA_Naive) <- gsub(pattern = "YBP008L3.", replacement= "YBP008L3_", colnames(nucleoplasmic_b1_LPA_Naive))
colnames(nucleoplasmic_b1_LPA_Naive) <- gsub(pattern = "YBP030L5.", replacement= "YBP030L5_", colnames(nucleoplasmic_b1_LPA_Naive))


# nucleoplasmic b3 ####
YAP006L8 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L8/08_counts/all/all_counts.txt', header = T)
YBP005L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L3/08_counts/all/all_counts.txt', header = T)
YAP012L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L1/08_counts/all/all_counts.txt', header = T)
YAP012L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L2/08_counts/all/all_counts.txt', header = T)
YAP012L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP012L3/08_counts/all/all_counts.txt', header = T)

nucleoplasmic_b3_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/nucleoplasm_b3_barcodes.txt')

nucleoplasmic_b3_LPA_Naive <- data.frame(YAP006L8[, 1:6], # Gene descriptions
                                         YAP006L8[, sapply(nucleoplasmic_b3_LPA_Naive_barcodes[nucleoplasmic_b3_LPA_Naive_barcodes[,5] == "YAP006L8",6], function(x) {grep(x, names(YAP006L8))})], 
                                         YBP005L3[, sapply(nucleoplasmic_b3_LPA_Naive_barcodes[nucleoplasmic_b3_LPA_Naive_barcodes[,5] == "YBP005L3",6], function(x) {grep(x, names(YBP005L3))})], 
                                         YAP012L1[, sapply(nucleoplasmic_b3_LPA_Naive_barcodes[nucleoplasmic_b3_LPA_Naive_barcodes[,5] == "YAP012L1",6], function(x) {grep(x, names(YAP012L1))})], 
                                         YAP012L2[, sapply(nucleoplasmic_b3_LPA_Naive_barcodes[nucleoplasmic_b3_LPA_Naive_barcodes[,5] == "YAP012L2",6], function(x) {grep(x, names(YAP012L2))})], 
                                         YAP012L3[, sapply(nucleoplasmic_b3_LPA_Naive_barcodes[nucleoplasmic_b3_LPA_Naive_barcodes[,5] == "YAP012L3",6], function(x) {grep(x, names(YAP012L3))})])

rm(YAP006L8, YBP005L3, YAP012L1, YAP012L2, YAP012L3)

## all nucleoplasmic #####
nucleoplasmic_all_LPA_Naive <- data.frame(nucleoplasmic_b1_LPA_Naive[, 1:6], 
                                          nucleoplasmic_b1_LPA_Naive[, -c(1:6)], 
                                          nucleoplasmic_b3_LPA_Naive[, -c(1:6)])

nucleoplasmic_all_LPA_Naive_barcodes <- rbind(nucleoplasmic_b1_LPA_Naive_barcodes, nucleoplasmic_b3_LPA_Naive_barcodes)

# p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='all')
# p + labs(title='nucleoplasmic_LPA_Naive (All)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='high')
# p + labs(title='nucleoplasmic_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='top')
# p + labs(title='nucleoplasmic_LPA_Naive (Highly expressed genes and high variance)')

p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes, type='lpa')
p + labs(title='Nucleoplasmic (LPA induced genes)') + theme(aspect.ratio=1)

p <- plot_pca_many_exp(nucleoplasmic_all_LPA_Naive, nucleoplasmic_all_LPA_Naive_barcodes[nucleoplasmic_all_LPA_Naive_barcodes$V2 == "Naive",], type='lpa')
p + labs(title='Nucleoplasmic (LPA induced genes)') + theme(aspect.ratio=1)

# cytoplasmic b0 #####
# Xap103L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap103L1/08_counts/all/all_counts.txt', header = T)
# VBP026L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L2/08_counts/all/all_counts.txt', header = T)
# VBP026L3 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L3/08_counts/all/all_counts.txt', header = T)
# Xbp081L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXbp081L1/08_counts/all/all_counts.txt', header = T)
# 
# cytoplasmic_b0_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/cytoplasmic_b0_barcodes.txt')
# 
# cytoplasmic_b0_LPA_Naive <- data.frame(Xap103L1[, 1:6], # Gene descriptions
#                                        Xap103L1[, sapply(cytoplasmic_b0_LPA_Naive_barcodes[cytoplasmic_b0_LPA_Naive_barcodes[,5] == "Xap103L1",6], function(x) {grep(x, names(Xap103L1))})], 
#                                        VBP026L2[, sapply(cytoplasmic_b0_LPA_Naive_barcodes[cytoplasmic_b0_LPA_Naive_barcodes[,5] == "VBP026L2",6], function(x) {grep(x, names(VBP026L2))})], 
#                                        VBP026L3[, sapply(cytoplasmic_b0_LPA_Naive_barcodes[cytoplasmic_b0_LPA_Naive_barcodes[,5] == "VBP026L3",6], function(x) {grep(x, names(VBP026L3))})], 
#                                        Xbp081L1[, sapply(cytoplasmic_b0_LPA_Naive_barcodes[cytoplasmic_b0_LPA_Naive_barcodes[,5] == "Xbp081L1",6], function(x) {grep(x, names(Xbp081L1))})])
# 
# rm(Xap103L1,VBP026L2,VBP026L3,Xbp081L1)


# Xap103L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap103L1/08_counts/exons/exon_counts.txt', header = T)
# VBP026L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L2/08_counts/exons/exon_counts.txt', header = T)
# VBP026L3 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L3/08_counts/exons/exon_counts.txt', header = T)
# Xbp081L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXbp081L1/08_counts/exons/exon_counts.txt', header = T)
# 
# cytoplasmic_b0_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/cytoplasmic_b0_barcodes.txt')
# 
# cytoplasmic_b0_LPA_Naive_exons <- data.frame(Xap103L1[, 1:6], # Gene descriptions
#                                              Xap103L1[, -c(1:6)], 
#                                              VBP026L2[, -c(1:6)], 
#                                              VBP026L3[, -c(1:6)], 
#                                              Xbp081L1[, -c(1:6)])
# 
# rm(Xap103L1,VBP026L2,VBP026L3,Xbp081L1)


# cytoplasmic b1 #####
Xap103L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap103L1/08_counts/all/all_counts.txt', header = T)
VBP026L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L2/08_counts/all/all_counts.txt', header = T)
VBP026L3 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L3/08_counts/all/all_counts.txt', header = T)
Xbp081L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXbp081L1/08_counts/all/all_counts.txt', header = T)
Xap097L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap097L2/08_counts/all/all_counts.txt', header = T)
Xap115L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap115L1/08_counts/all/all_counts.txt', header = T)
Xap115L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXap115L2/08_counts/all/all_counts.txt', header = T)
VBP026L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsVBP026L1/08_counts/all/all_counts.txt', header = T)
YBP008L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1_vs_M2/lanes/SxaQSEQsYBP008L3/08_counts/all/all_counts.txt', header = T)
YBP030L4 <- read.table( file = '~/Documents/Biggie_Expansion/Diane/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP030L4/08_counts/all/all_counts.txt', header = T)

cytoplasmic_b1_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/cytoplasmic_b1_barcodes_corrected.txt')

cytoplasmic_b1_LPA_Naive <- data.frame(Xap103L1[, 1:6], # Gene descriptions
                                      Xap103L1[, sapply(cytoplasmic_b1_LPA_Naive_barcodes[cytoplasmic_b1_LPA_Naive_barcodes[,5] == "Xap103L1",6], function(x) {grep(x, names(Xap103L1))})], 
                                      VBP026L2[, sapply(cytoplasmic_b1_LPA_Naive_barcodes[cytoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L2",6], function(x) {grep(x, names(VBP026L2))})], 
                                      VBP026L3[, sapply(cytoplasmic_b1_LPA_Naive_barcodes[cytoplasmic_b1_LPA_Naive_barcodes[,5] == "VBP026L3",6], function(x) {grep(x, names(VBP026L3))})], 
                                      Xbp081L1[, sapply(cytoplasmic_b1_LPA_Naive_barcodes[cytoplasmic_b1_LPA_Naive_barcodes[,5] == "Xbp081L1",6], function(x) {grep(x, names(Xbp081L1))})], 
                                      Xap097L2[, -c(1:6)], 
                                      Xap115L1[, -c(1:6)], 
                                      Xap115L2[, -c(1:6)], 
                                      VBP026L1[, -c(1:6)], 
                                      YBP008L3[, -c(1:6)],
                                      YBP030L4[, -c(1:6)])

rm(Xap103L1,VBP026L2,VBP026L3,Xbp081L1,Xap097L2,Xap115L1, Xap115L2, VBP026L1, YBP008L3, YBP030L4)

colnames(cytoplasmic_b1_LPA_Naive) <- gsub('YBP008L3.', 'YBP008L3_', colnames(cytoplasmic_b1_LPA_Naive))
colnames(cytoplasmic_b1_LPA_Naive) <- gsub('YBP030L4.', 'YBP030L4_', colnames(cytoplasmic_b1_LPA_Naive))

# 
# p <- plot_pca_few_exp(cytoplasmic_b1_LPA_Naive_exons, cytoplasmic_b1_LPA_Naive_barcodes, type='all')
# p + labs(title='Cytoplasmic_b1_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_LPA_Naive, cytoplasmic_b1_LPA_Naive_barcodes, type='high')
# p + labs(title='Cytoplasmic_b1_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_LPA_Naive, cytoplasmic_b1_LPA_Naive_barcodes, type='top')
# p + labs(title='Cytoplasmic_b1_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_LPA_Naive_exons, cytoplasmic_b1_LPA_Naive_barcodes, type='lpa')
# p + labs(title='Cytoplasmic_b1_LPA_Naive (LPA induced genes)')


# cytoplasmic b2 ##### => actually b3
YAP015L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L1/08_counts/all/all_counts.txt', header = T)
YAP015L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L2/08_counts/all/all_counts.txt', header = T)
YAP015L3 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L3/08_counts/all/all_counts.txt', header = T)
YAP015L4 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L4/08_counts/all/all_counts.txt', header = T)
YBP005L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L2/08_counts/all/all_counts.txt', header = T)

cytoplasmic_b3_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/cytoplasmic_b3_barcodes_corrected.txt')
cytoplasmic_b3_LPA_Naive_barcodes$V4 <- "b3"

cytoplasmic_b3_LPA_Naive <- data.frame(YAP015L1[, 1:6], # Gene descriptions
                                       YAP015L1[, -c(1:6)], 
                                       YAP015L2[, -c(1:6)], 
                                       YAP015L3[, -c(1:6)], 
                                       YAP015L4[, -c(1:6)], 
                                       YBP005L2[, -c(1:6)])

rm(YAP015L1,YAP015L2,YAP015L3,YAP015L4,YBP005L2)

# p <- plot_pca_few_exp(cytoplasmic_b2_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='all')
# p + labs(title='Cytoplasmic_b2_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='high')
# p + labs(title='Cytoplasmic_b2_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='top')
# p + labs(title='Cytoplasmic_b2_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='lpa')
# p + labs(title='Cytoplasmic_b2_LPA_Naive (LPA induced genes)')

# cytoplasmic b3 ##### => actually b2
YAP015L5 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L5/08_counts/all/all_counts.txt', header = T)
YAP015L6 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L6/08_counts/all/all_counts.txt', header = T)
YAP015L7 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L7/08_counts/all/all_counts.txt', header = T)
YAP015L8 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP015L8/08_counts/all/all_counts.txt', header = T)
YBP005L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L1/08_counts/all/all_counts.txt', header = T)

cytoplasmic_b2_LPA_Naive_barcodes <- read.table(file = '../Barcodes/LPA_vs_Naive/cytoplasmic_b2_barcodes_corrected.txt')
cytoplasmic_b2_LPA_Naive_barcodes$V4 <- "b2"

cytoplasmic_b2_LPA_Naive <- data.frame(YAP015L5[, 1:6], # Gene descriptions
                                       YAP015L5[, -c(1:6)], 
                                       YAP015L6[, -c(1:6)], 
                                       YAP015L7[, -c(1:6)], 
                                       YAP015L8[, -c(1:6)], 
                                       YBP005L1[, -c(1:6)])

rm(YAP015L5,YAP015L6,YAP015L7,YAP015L8,YBP005L1)

p <- plot_pca_few_exp(cytoplasmic_b3_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='all')
p + labs(title='Cytoplasmic_b3_LPA_Naive (All)')

p <- plot_pca_few_exp(cytoplasmic_b3_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='high')
p + labs(title='Cytoplasmic_b3_LPA_Naive (Highly expressed genes)')

p <- plot_pca_few_exp(cytoplasmic_b3_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='top')
p + labs(title='Cytoplasmic_b3_LPA_Naive (Highly expressed genes and high variance)')

p <- plot_pca_few_exp(cytoplasmic_b3_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='lpa')
p + labs(title='Cytoplasmic_b3_LPA_Naive (LPA induced genes)')

p <- plot_pca_few_exp(cytoplasmic_b3_LPA_Naive, cytoplasmic_b2_LPA_Naive_barcodes, type='lpa')
p + labs(title='Cytoplasmic_b3_LPA_Naive (LPA induced genes)')

## all cytoplamic #####
cytoplasmic_all_LPA_Naive <- data.frame(cytoplasmic_b1_LPA_Naive[, 1:6], 
                                  cytoplasmic_b1_LPA_Naive[, -c(1:6)],
                                  cytoplasmic_b3_LPA_Naive[,-c(1:6)],
                                  cytoplasmic_b2_LPA_Naive[, -c(1:6)])

colnames(cytoplasmic_all_LPA_Naive)

cytoplasmic_all_LPA_Naive_barcodes <- rbind(cytoplasmic_b1_LPA_Naive_barcodes,cytoplasmic_b3_LPA_Naive_barcodes, cytoplasmic_b2_LPA_Naive_barcodes)
head(cytoplasmic_all_LPA_Naive_barcodes)

# cytoplasmic_all_LPA_Naive_barcodes[,4] <- as.character(cytoplasmic_all_LPA_Naive_barcodes[,4])
# cytoplasmic_all_LPA_Naive_barcodes[cytoplasmic_all_LPA_Naive_barcodes[,4] == "b2", 4] <- "b4"
# cytoplasmic_all_LPA_Naive_barcodes[cytoplasmic_all_LPA_Naive_barcodes[,4] == "b3", 4] <- "b2"
# cytoplasmic_all_LPA_Naive_barcodes[cytoplasmic_all_LPA_Naive_barcodes[,4] == "b4", 4] <- "b3"

# cytoplasmic_all_LPA_Naive_barcodes[,4] <- as.factor(cytoplasmic_all_LPA_Naive_barcodes[,4])


# p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='all')
# p + labs(title='Cytoplasmic_LPA_Naive (All)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='high')
# p + labs(title='Cytoplasmic_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='top')
# p + labs(title='Cytoplasmic_LPA_Naive (Highly expressed genes and high variance)')

p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes, type='lpa')
p + labs(title='Cytoplasmic (LPA induced genes)') + theme(aspect.ratio=1)

p <- plot_pca_many_exp(cytoplasmic_all_LPA_Naive, cytoplasmic_all_LPA_Naive_barcodes[cytoplasmic_all_LPA_Naive_barcodes$V2=='Naive',], type='lpa')
p + labs(title='Cytoplasmic (LPA induced genes)') + theme(aspect.ratio=1)


# # polyA b0 #####
# VBP026L8 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/polyA_b1/VBP026L8/all_counts.txt', header = T)
# Xap099L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/polyA_b1/Xapp099L1/all_counts.txt', header = T)
# Xbp081L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXbp081L2/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b0_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/polyA_b0_barcodes.txt')
# 
# polyA_b0_LPA_Naive <- data.frame(VBP026L8[, 1:6], # Gene descriptions
#                                  VBP026L8[, -c(1:6)], 
#                                  Xap099L1[, -c(1:6)], 
#                                  Xbp081L2[, -c(1:6)])
# 
# rm(VBP026L8,Xap099L1,Xbp081L2)
# 
# p <- plot_pca_few_exp(polyA_b0_LPA_Naive, polyA_b0_LPA_Naive_barcodes, type='all')
# p + labs(title='polyA_b0_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(polyA_b0_LPA_Naive, polyA_b0_LPA_Naive_barcodes, type='high')
# p + labs(title='polyA_b0_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b0_LPA_Naive, polyA_b0_LPA_Naive_barcodes, type='top')
# p + labs(title='polyA_b0_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b0_LPA_Naive, polyA_b0_LPA_Naive_barcodes, type='lpa')
# p + labs(title='polyA_b0_LPA_Naive (LPA induced genes)')
# 
# 
# # polyA b1 #####
# VBP026L8 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/polyA_b1/VBP026L8/all_counts.txt', header = T)
# Xap099L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/polyA_b1/Xapp099L1/all_counts.txt', header = T)
# Xbp081L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsXbp081L2/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b1_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/polyA_b1_barcodes.txt')
# 
# polyA_b1_LPA_Naive <- data.frame(VBP026L8[, 1:6], # Gene descriptions
#                                  VBP026L8[, -c(1:6)], 
#                                  Xap099L1[, -c(1:6)], 
#                                  Xbp081L2[, -c(1:6)])
# 
# rm(VBP026L8,Xap099L1,Xbp081L2)
# 
# p <- plot_pca_few_exp(polyA_b1_LPA_Naive, polyA_b1_LPA_Naive_barcodes, type='all')
# p + labs(title='polyA_b1_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(polyA_b1_LPA_Naive, polyA_b1_LPA_Naive_barcodes, type='high')
# p + labs(title='polyA_b1_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b1_LPA_Naive, polyA_b1_LPA_Naive_barcodes, type='top')
# p + labs(title='polyA_b1_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b1_LPA_Naive, polyA_b1_LPA_Naive_barcodes, type='lpa')
# p + labs(title='polyA_b1_LPA_Naive (LPA induced genes)')
# 
# # polyA b2 #####
# YAP006L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L1/08_counts/all/all_counts.txt', header = T)
# YBP005L7 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L7/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b2_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/polyA_b2_barcodes.txt')
# 
# polyA_b2_LPA_Naive <- data.frame(YAP006L1[, 1:6], # Gene descriptions
#                                  YAP006L1[, -c(1:6)], 
#                                  YBP005L7[, -c(1:6)])
# 
# rm(YAP006L1,YBP005L7)
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='all')
# p + labs(title='polyA_b2_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='high')
# p + labs(title='polyA_b2_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='top')
# p + labs(title='polyA_b2_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='lpa')
# p + labs(title='polyA_b2_LPA_Naive (LPA induced genes)')
# 
# # polyA b2 #####
# YAP006L1 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYAP006L1/08_counts/all/all_counts.txt', header = T)
# YBP005L7 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/SxaQSEQsYBP005L7/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b2_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/polyA_b2_barcodes.txt')
# 
# polyA_b2_LPA_Naive <- data.frame(YAP006L1[, 1:6], # Gene descriptions
#                                  YAP006L1[, -c(1:6)], 
#                                  YBP005L7[, -c(1:6)])
# 
# rm(YAP006L1,YBP005L7)
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='all')
# p + labs(title='polyA_b2_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='high')
# p + labs(title='polyA_b2_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='top')
# p + labs(title='polyA_b2_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b2_LPA_Naive, polyA_b2_LPA_Naive_barcodes, type='lpa')
# p + labs(title='polyA_b2_LPA_Naive (LPA induced genes)')
# 
# # polyA b3 #####
# YAP006L2 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/polyA_b3/YAP006L2/all_counts.txt', header = T)
# YBP005L8 <- read.table( file = '~/Documents/Disk1/data/ribonomics/LPA_vs_Naive/lanes/polyA_b3/YBP005L8/all_counts.txt', header = T)
# 
# polyA_b3_LPA_Naive_barcodes <- read.table(file = '../LPA vs Naive/polyA_b3_barcodes.txt')
# 
# polyA_b3_LPA_Naive <- data.frame(YAP006L2[, 1:6], # Gene descriptions
#                                  YAP006L2[, -c(1:6)], 
#                                  YBP005L8[, -c(1:6)])
# 
# rm(YAP006L2,YBP005L8)
# 
# p <- plot_pca_few_exp(polyA_b3_LPA_Naive, polyA_b3_LPA_Naive_barcodes, type='all')
# p + labs(title='polyA_b3_LPA_Naive (All)')
# 
# p <- plot_pca_few_exp(polyA_b3_LPA_Naive, polyA_b3_LPA_Naive_barcodes, type='high')
# p + labs(title='polyA_b3_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b3_LPA_Naive, polyA_b3_LPA_Naive_barcodes, type='top')
# p + labs(title='polyA_b3_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b3_LPA_Naive, polyA_b3_LPA_Naive_barcodes, type='lpa')
# p + labs(title='polyA_b3_LPA_Naive (LPA induced genes)')
# 
# ## all polyA #####
# polyA_all_LPA_Naive <- data.frame(polyA_b1_LPA_Naive[, 1:6], 
#                                   polyA_b0_LPA_Naive[, -c(1:6)],
#                                           polyA_b1_LPA_Naive[, -c(1:6)],
#                                   polyA_b2_LPA_Naive[,-c(1:6)],
#                                           polyA_b3_LPA_Naive[, -c(1:6)])
# 
# polyA_all_LPA_Naive_barcodes <- rbind(polyA_b0_LPA_Naive_barcodes,polyA_b1_LPA_Naive_barcodes,polyA_b2_LPA_Naive_barcodes, polyA_b3_LPA_Naive_barcodes)
# 
# p <- plot_pca_many_exp(polyA_all_LPA_Naive, polyA_all_LPA_Naive_barcodes, type='all')
# p + labs(title='polyA_LPA_Naive (All)')
# 
# p <- plot_pca_many_exp(polyA_all_LPA_Naive, polyA_all_LPA_Naive_barcodes, type='high')
# p + labs(title='polyA_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(polyA_all_LPA_Naive, polyA_all_LPA_Naive_barcodes, type='top')
# p + labs(title='polyA_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(polyA_all_LPA_Naive, polyA_all_LPA_Naive_barcodes, type='lpa')
# p + labs(title='polyA_LPA_Naive (LPA induced genes)')
# 
# 
# 
# ## all b1 ----
# b1_all_LPA_Naive <- data.frame(chromatin_b1_LPA_Naive[, 1:6], 
#                                           chromatin_b1_LPA_Naive[, -c(1:6)], 
#                                           nucleoplasmic_b1_LPA_Naive[, -c(1:6)], 
#                                           cytoplasmic_b1_LPA_Naive[, -c(1:6)], 
#                                           polyA_b1_LPA_Naive[,-c(1:6)])
# 
# b1_all_LPA_Naive_barcodes <- rbind(chromatin_b1_LPA_Naive_barcodes, nucleoplasmic_b1_LPA_Naive_barcodes, cytoplasmic_b1_LPA_Naive_barcodes, polyA_b1_LPA_Naive_barcodes)
# 
# p <- plot_pca_many_exp_cpt(b1_all_LPA_Naive, b1_all_LPA_Naive_barcodes, type='all')
# p + labs(title='b1_LPA_Naive (All)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_LPA_Naive, b1_all_LPA_Naive_barcodes, type='high')
# p + labs(title='b1_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_LPA_Naive, b1_all_LPA_Naive_barcodes, type='top')
# p + labs(title='b1_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_LPA_Naive, b1_all_LPA_Naive_barcodes, type='lpa')
# p + labs(title='b1_LPA_Naive (LPA induced genes)')
# 
# 
# ## all b2 ----
# b2_all_LPA_Naive <- data.frame(chromatin_b2_LPA_Naive[, 1:6], 
#                                chromatin_b2_LPA_Naive[, -c(1:6)],
#                                cytoplasmic_b2_LPA_Naive[, -c(1:6)],
#                                polyA_b2_LPA_Naive[, -c(1:6)])
# 
# b2_all_LPA_Naive_barcodes <- rbind(chromatin_b2_LPA_Naive_barcodes, cytoplasmic_b2_LPA_Naive_barcodes, polyA_b2_LPA_Naive_barcodes)
# 
# p <- plot_pca_many_exp_cpt(b2_all_LPA_Naive, b2_all_LPA_Naive_barcodes, type='all')
# p + labs(title='b2_LPA_Naive (All)')
# 
# p <- plot_pca_many_exp_cpt(b2_all_LPA_Naive, b2_all_LPA_Naive_barcodes, type='high')
# p + labs(title='b2_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b2_all_LPA_Naive, b2_all_LPA_Naive_barcodes, type='top')
# p + labs(title='b2_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b2_all_LPA_Naive, b2_all_LPA_Naive_barcodes, type='lpa')
# p + labs(title='b2_LPA_Naive (LPA induced genes)')
# 
# 
# 
# ## all b3 ----
# b3_all_LPA_Naive <- data.frame(chromatin_b3_LPA_Naive[, 1:6], 
#                                chromatin_b3_LPA_Naive[, -c(1:6)], 
#                                nucleoplasmic_b3_LPA_Naive[, -c(1:6)], 
#                                cytoplasmic_b3_LPA_Naive[, -c(1:6)],
#                                polyA_b3_LPA_Naive[, -c(1:6)])
# 
# b3_all_LPA_Naive_barcodes <- rbind(chromatin_b3_LPA_Naive_barcodes, nucleoplasmic_b3_LPA_Naive_barcodes, cytoplasmic_b3_LPA_Naive_barcodes, polyA_b3_LPA_Naive_barcodes)
# 
# p <- plot_pca_many_exp_cpt(b3_all_LPA_Naive, b3_all_LPA_Naive_barcodes, type='all')
# p + labs(title='b3_LPA_Naive (All)')
# 
# p <- plot_pca_many_exp_cpt(b3_all_LPA_Naive, b3_all_LPA_Naive_barcodes, type='high')
# p + labs(title='b3_LPA_Naive (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b3_all_LPA_Naive, b3_all_LPA_Naive_barcodes, type='top')
# p + labs(title='b3_LPA_Naive (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b3_all_LPA_Naive, b3_all_LPA_Naive_barcodes, type='lpa')
# p + labs(title='b3_LPA_Naive (LPA induced genes)')
# 
# 
# 
# ############
# ## M1 M2 ###
# ############
# 
# ## chromatin b1 ####
# YBP008L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsYBP008L2/08_counts/all/all_counts.txt', header = T)
# 
# chromatin_b1_M1_M2_barcodes <- read.table(file = '../M1 vs M2/chromatin_b1_barcodes.txt')
# 
# chromatin_b1_M1_M2 <- data.frame(YBP008L2[, 1:6], # Gene descriptions
#                            YBP008L2[, -c(1:6)])
# 
# rm(YBP008L2)
# colnames(chromatin_b1_M1_M2) <- gsub('YBP008L2.', 'YBP008L2_', colnames(chromatin_b1_M1_M2))
# 
# p <- plot_pca_few_exp(chromatin_b1_M1_M2, chromatin_b1_M1_M2_barcodes, type='all')
# p + labs(title='Chromatin_b1_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(chromatin_b1_M1_M2, chromatin_b1_M1_M2_barcodes, type='high')
# p + labs(title='Chromatin_b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(chromatin_b1_M1_M2, chromatin_b1_M1_M2_barcodes, type='top')
# p + labs(title='Chromatin_b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(chromatin_b1_M1_M2, chromatin_b1_M1_M2_barcodes, type='lpa')
# p + labs(title='Chromatin_b1_M1_M2 (LPA induced genes)')
# 
# 
# ## nucleoplasmic b1 ####
# YBP008L5 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsYBP008L5/08_counts/all/all_counts.txt', header = T)
# 
# nucleoplasmic_b1_M1_M2_barcodes <- read.table(file = '../M1 vs M2/nucleoplasmic_b1_barcodes.txt')
# 
# nucleoplasmic_b1_M1_M2 <- data.frame(YBP008L5[, 1:6], # Gene descriptions
#                                YBP008L5[, -c(1:6)])
# 
# rm(YBP008L5)
# 
# colnames(nucleoplasmic_b1_M1_M2) <- gsub('YBP008L5.', 'YBP008L5_', colnames(nucleoplasmic_b1_M1_M2))
# 
# p <- plot_pca_few_exp(nucleoplasmic_b1_M1_M2, nucleoplasmic_b1_M1_M2_barcodes, type='all')
# p + labs(title='Nucleoplasmic_b1_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b1_M1_M2, nucleoplasmic_b1_M1_M2_barcodes, type='high')
# p + labs(title='Nucleoplasmic_b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b1_M1_M2, nucleoplasmic_b1_M1_M2_barcodes, type='top')
# p + labs(title='Nucleoplasmic_b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b1_M1_M2, nucleoplasmic_b1_M1_M2_barcodes, type='lpa')
# p + labs(title='Nucleoplasmic_b1_M1_M2 (LPA induced genes)')
# 
# 
# 
# 
# 
# ## cytoplasmic b1 #####
# YBP008L4 <- read.table(file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsYBP008L4/08_counts/all/all_counts.txt', header = T)
# 
# cytoplasmic_b1_M1_M2_barcodes <- read.table(file = '../M1 vs M2/cytoplasmic_b1_barcodes.txt')
# 
# cytoplasmic_b1_M1_M2 <- data.frame(YBP008L4[, 1:6], # Gene descriptions
#                                    YBP008L4[, -c(1:6)])
# 
# rm(YBP008L4)
# 
# colnames(cytoplasmic_b1_M1_M2) <- gsub('YBP008L4.', 'YBP008L4_', colnames(cytoplasmic_b1_M1_M2))
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2, cytoplasmic_b1_M1_M2_barcodes, type='all')
# p + labs(title='Cytoplasmic_b1_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2[,-7], cytoplasmic_b1_M1_M2_barcodes, type='all')
# p + labs(title='Cytoplasmic_b1_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2, cytoplasmic_b1_M1_M2_barcodes, type='high')
# p + labs(title='Cytoplasmic_b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2[,-7], cytoplasmic_b1_M1_M2_barcodes, type='high')
# p + labs(title='Cytoplasmic_b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2, cytoplasmic_b1_M1_M2_barcodes, type='top')
# p + labs(title='Cytoplasmic_b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2[,-7], cytoplasmic_b1_M1_M2_barcodes, type='top')
# p + labs(title='Cytoplasmic_b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2, cytoplasmic_b1_M1_M2_barcodes, type='lpa')
# p + labs(title='Cytoplasmic_b1_M1_M2 (LPA induced genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b1_M1_M2[,-7], cytoplasmic_b1_M1_M2_barcodes, type='lpa')
# p + labs(title='Cytoplasmic_b1_M1_M2 (LPA induced genes)')
# 
# 
# 
# 
# ## polyA b1 ####
# YBP008L3 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsYBP008L3/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b1_M1_M2_barcodes <- read.table(file = '../M1 vs M2/polyA_b1_barcodes.txt')
# 
# barcode <- unlist(lapply(strsplit(colnames(YBP008L3)[-(1:6)], '.', fixed = T), function(x){x[4]}))
# polyA_b1_M1_M2 <- data.frame(YBP008L3[, 1:6], # Gene descriptions
#                              YBP008L3[, -c(1:6)][,barcode %in% polyA_b1_M1_M2_barcodes$V6])
# 
# rm(YBP008L3)
# colnames(polyA_b1_M1_M2) <- gsub('YBP008L3.', 'YBP008L3_', colnames(polyA_b1_M1_M2))
# 
# p <- plot_pca_few_exp(polyA_b1_M1_M2, polyA_b1_M1_M2_barcodes, type='all')
# p + labs(title='PolyA_b1_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(polyA_b1_M1_M2, polyA_b1_M1_M2_barcodes, type='high')
# p + labs(title='PolyA_b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b1_M1_M2, polyA_b1_M1_M2_barcodes, type='top')
# p + labs(title='PolyA_b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b1_M1_M2, polyA_b1_M1_M2_barcodes, type='lpa')
# p + labs(title='PolyA_b1_M1_M2 (LPA induced genes)')
# 
# ## chromatin b2 ####
# Xbp115L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp115L1/08_counts/all/all_counts.txt', header = T)
# 
# chromatin_b2_M1_M2_barcodes <- read.table(file = '../M1 vs M2/chromatin_b2_barcodes.txt')
# 
# chromatin_b2_M1_M2 <- data.frame(Xbp115L1[, 1:6], # Gene descriptions
#                            Xbp115L1[, -c(1:6)])
# 
# rm(Xbp115L1)
# colnames(chromatin_b2_M1_M2) <- gsub('Xbp115L1.', 'Xbp115L1_', colnames(chromatin_b2_M1_M2))
# 
# p <- plot_pca_few_exp(chromatin_b2_M1_M2, chromatin_b2_M1_M2_barcodes, type='all')
# p + labs(title='Chromatin_b2_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(chromatin_b2_M1_M2, chromatin_b2_M1_M2_barcodes, type='high')
# p + labs(title='Chromatin_b2_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(chromatin_b2_M1_M2, chromatin_b2_M1_M2_barcodes, type='top')
# p + labs(title='Chromatin_b2_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(chromatin_b2_M1_M2, chromatin_b2_M1_M2_barcodes, type='lpa')
# p + labs(title='Chromatin_b2_M1_M2 (LPA induced genes)')
# 
# 
# 
# ## nucleoplasmic b2 ####
# Xbp124L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp124L1/08_counts/all/all_counts.txt', header = T)
# 
# nucleoplasmic_b2_M1_M2_barcodes <- read.table(file = '../M1 vs M2/nucleoplasmic_b2_barcodes.txt')
# 
# nucleoplasmic_b2_M1_M2 <- data.frame(Xbp124L1[, 1:6], # Gene descriptions
#                                Xbp124L1[, -c(1:6)])
# 
# rm(Xbp124L1)
# 
# colnames(nucleoplasmic_b2_M1_M2) <- gsub('Xbp124L1.', 'Xbp124L1_', colnames(nucleoplasmic_b2_M1_M2))
# 
# p <- plot_pca_few_exp(nucleoplasmic_b2_M1_M2, nucleoplasmic_b2_M1_M2_barcodes, type='all')
# p + labs(title='Nucleoplasmic_b2_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b2_M1_M2, nucleoplasmic_b2_M1_M2_barcodes, type='high')
# p + labs(title='Nucleoplasmic_b2_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b2_M1_M2, nucleoplasmic_b2_M1_M2_barcodes, type='top')
# p + labs(title='Nucleoplasmic_b2_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b2_M1_M2, nucleoplasmic_b2_M1_M2_barcodes, type='lpa')
# p + labs(title='Nucleoplasmic_b2_M1_M2 (LPA induced genes)')
# 
# 
# ## polyA b2 ####
# Xbp107L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp107L1/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b2_M1_M2_barcodes <- read.table(file = '../M1 vs M2/polyA_b2_barcodes.txt')
# 
# polyA_b2_M1_M2 <- data.frame(Xbp107L1[, 1:6], # Gene descriptions
#                              Xbp107L1[, -c(1:6)])
# 
# rm(Xbp107L1)
# colnames(polyA_b2_M1_M2) <- gsub('XBP107L1.', 'Xbp107L1_', colnames(polyA_b2_M1_M2))
# 
# p <- plot_pca_few_exp(polyA_b2_M1_M2, polyA_b2_M1_M2_barcodes, type='all')
# p + labs(title='PolyA_b2_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(polyA_b2_M1_M2, polyA_b2_M1_M2_barcodes, type='high')
# p + labs(title='PolyA_b2_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b2_M1_M2, polyA_b2_M1_M2_barcodes, type='top')
# p + labs(title='PolyA_b2_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b2_M1_M2, polyA_b2_M1_M2_barcodes, type='lpa')
# p + labs(title='PolyA_b2_M1_M2 (LPA induced genes)')
# 
# ## cytoplasmic b2 #####
# Xbp114L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp114L1/08_counts/all/all_counts.txt', header = T)
# 
# cytoplasmic_b2_M1_M2_barcodes <- read.table(file = '../M1 vs M2/cytoplasmic_b2_barcodes.txt')
# 
# cytoplasmic_b2_M1_M2 <- data.frame(Xbp114L1[, 1:6], # Gene descriptions
#                              Xbp114L1[, -c(1:6)])
# 
# rm(Xbp114L1)
# 
# colnames(cytoplasmic_b2_M1_M2) <- gsub('Xbp114L1.', 'Xbp114L1_', colnames(cytoplasmic_b2_M1_M2))
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_M1_M2, cytoplasmic_b2_M1_M2_barcodes, type='all')
# p + labs(title='Cytoplasmic_b2_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_M1_M2, cytoplasmic_b2_M1_M2_barcodes, type='high')
# p + labs(title='Cytoplasmic_b2_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_M1_M2, cytoplasmic_b2_M1_M2_barcodes, type='top')
# p + labs(title='Cytoplasmic_b2_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b2_M1_M2, cytoplasmic_b2_M1_M2_barcodes, type='lpa')
# p + labs(title='Cytoplasmic_b2_M1_M2 (LPA induced genes)')
# 
# 
# 
# 
# 
# 
# ## chromatin b3 ####
# Xbp126L1 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp126L1/08_counts/all/all_counts.txt', header = T)
# 
# chromatin_b3_M1_M2_barcodes <- read.table(file = '../M1 vs M2/chromatin_b3_barcodes.txt')
# 
# chromatin_b3_M1_M2 <- data.frame(Xbp126L1[, 1:6], # Gene descriptions
#                                      Xbp126L1[, -c(1:6)])
# 
# rm(Xbp126L1)
# 
# p <- plot_pca_few_exp(chromatin_b3_M1_M2, chromatin_b3_M1_M2_barcodes, type='all')
# p + labs(title='Chromatin_b3_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(chromatin_b3_M1_M2, chromatin_b3_M1_M2_barcodes, type='high')
# p + labs(title='Chromatin_b3_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(chromatin_b3_M1_M2, chromatin_b3_M1_M2_barcodes, type='top')
# p + labs(title='Chromatin_b3_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(chromatin_b3_M1_M2, chromatin_b3_M1_M2_barcodes, type='lpa')
# p + labs(title='Chromatin_b3_M1_M2 (LPA induced genes)')
# 
# 
# 
# ## nucleoplasmic b3 ####
# Xbp126L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp126L2/08_counts/all/all_counts.txt', header = T)
# 
# nucleoplasmic_b3_M1_M2_barcodes <- read.table(file = '../M1 vs M2/nucleoplasmic_b3_barcodes.txt')
# 
# nucleoplasmic_b3_M1_M2 <- data.frame(Xbp126L2[, 1:6], # Gene descriptions
#                                      Xbp126L2[, -c(1:6)])
# 
# rm(Xbp126L2)
# 
# p <- plot_pca_few_exp(nucleoplasmic_b3_M1_M2, nucleoplasmic_b3_M1_M2_barcodes, type='all')
# p + labs(title='Nucleoplasmic_b3_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b3_M1_M2, nucleoplasmic_b3_M1_M2_barcodes, type='high')
# p + labs(title='Nucleoplasmic_b3_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b3_M1_M2, nucleoplasmic_b3_M1_M2_barcodes, type='top')
# p + labs(title='Nucleoplasmic_b3_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(nucleoplasmic_b3_M1_M2, nucleoplasmic_b3_M1_M2_barcodes, type='lpa')
# p + labs(title='Nucleoplasmic_b3_M1_M2 (LPA induced genes)')
# 
# 
# ## cytoplasmic b3 #####
# Xbp124L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp124L2/08_counts/all/all_counts.txt', header = T)
# 
# cytoplasmic_b3_M1_M2_barcodes <- read.table(file = '../M1 vs M2/cytoplasmic_b3_barcodes.txt')
# 
# cytoplasmic_b3_M1_M2 <- data.frame(Xbp124L2[, 1:6], # Gene descriptions
#                              Xbp124L2[, -c(1:6)])
# 
# rm(Xbp124L2)
# 
# colnames(cytoplasmic_b3_M1_M2) <- gsub('Xbp124L2.', 'Xbp124L2_', colnames(cytoplasmic_b3_M1_M2))
# 
# p <- plot_pca_few_exp(cytoplasmic_b3_M1_M2, cytoplasmic_b3_M1_M2_barcodes, type='all')
# p + labs(title='Cytoplasmic_b3_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b3_M1_M2, cytoplasmic_b3_M1_M2_barcodes, type='high')
# p + labs(title='Cytoplasmic_b3_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b3_M1_M2, cytoplasmic_b3_M1_M2_barcodes, type='top')
# p + labs(title='Cytoplasmic_b3_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(cytoplasmic_b3_M1_M2, cytoplasmic_b3_M1_M2_barcodes, type='lpa')
# p + labs(title='Cytoplasmic_b3_M1_M2 (LPA induced genes)')
# 
# 
# 
# ## PolyA b3 #####
# 
# Xbp114L2 <- read.table( file = '~/Documents/Biggie/data2/diane/ribonomics/M1M2_RNAseq/SxaQSEQsXbp114L2/08_counts/all/all_counts.txt', header = T)
# 
# polyA_b3_M1_M2_barcodes <- read.table(file = '../M1 vs M2/polyA_b3_barcodes.txt')
# 
# polyA_b3_M1_M2 <- data.frame(Xbp114L2[, 1:6], # Gene descriptions
#                        Xbp114L2[, -c(1:6)])
# 
# rm(Xbp114L2)
# 
# colnames(polyA_b3_M1_M2) <- gsub('Xbp114L2.', 'Xbp114L2_', colnames(polyA_b3_M1_M2))
# 
# p <- plot_pca_few_exp(polyA_b3_M1_M2, polyA_b3_M1_M2_barcodes, type='all')
# p + labs(title='PolyA_b3_M1_M2 (All)')
# 
# p <- plot_pca_few_exp(polyA_b3_M1_M2, polyA_b3_M1_M2_barcodes, type='high')
# p + labs(title='PolyA_b3_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_few_exp(polyA_b3_M1_M2, polyA_b3_M1_M2_barcodes, type='top')
# p + labs(title='PolyA_b3_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_few_exp(polyA_b3_M1_M2, polyA_b3_M1_M2_barcodes, type='lpa')
# p + labs(title='PolyA_b3_M1_M2 (LPA induced genes)')
# 
# # compare with roberto => Not the same annotations => Hard to compare
# Xap114L1_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsXap114L1/06_counts/counts.txt', header = T)
# Xap114L2_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsXap114L2/06_counts/counts.txt', header = T)
# Xbp051L1_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsXbp051L1/06_counts/counts.txt', header = T)
# WAP033L1_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsWAP033L1/06_counts/counts.txt', header = T)
# WAP033L2_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsWAP033L2/06_counts/counts.txt', header = T)
# WAP033L3_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsWAP033L3/06_counts/counts.txt', header = T)
# WAP033L4_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsWAP033L4/06_counts/counts.txt', header = T)
# WAP033L5_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsWAP033L5/06_counts/counts.txt', header = T)
# WAP033L6_roberto <- read.table( file = '~/Documents/Biggie/data1/roberto/ribonomics/SxaQSEQsWAP033L6/06_counts/counts.txt', header = T)
# 
# rm(Xap114L1_roberto,Xap114L2_roberto,Xbp051L1_roberto,WAP033L1_roberto,WAP033L2_roberto,WAP033L3_roberto,WAP033L4_roberto, WAP033L5_roberto, WAP033L6_roberto)
# 
# 
# 
# 
# ## chromatin all M1 M2 ####
# chromatin_all_M1_M2 <- data.frame(chromatin_b1_M1_M2[, 1:6], 
#                                   chromatin_b1_M1_M2[, -c(1:6)], 
#                                   chromatin_b2_M1_M2[, -c(1:6)], 
#                                   chromatin_b3_M1_M2[, -c(1:6)])
# 
# chromatin_all_M1_M2_barcodes <- rbind(chromatin_b1_M1_M2_barcodes, chromatin_b2_M1_M2_barcodes, chromatin_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(chromatin_all_M1_M2, chromatin_all_M1_M2_barcodes, type='all')
# p + labs(title='Chromatin_M1_M2 (All)')
# 
# p <- plot_pca_many_exp(chromatin_all_M1_M2, chromatin_all_M1_M2_barcodes, type='high')
# p + labs(title='Chromatin_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(chromatin_all_M1_M2, chromatin_all_M1_M2_barcodes, type='top')
# p + labs(title='Chromatin_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(chromatin_all_M1_M2, chromatin_all_M1_M2_barcodes, type='lpa')
# p + labs(title='Chromatin_M1_M2 (LPA induced genes)')
# 
# 
# 
# ## nucleoplasmic all M1 M2 ####
# nucleoplasmic_all_M1_M2 <- data.frame(nucleoplasmic_b1_M1_M2[, 1:6], 
#                                   nucleoplasmic_b1_M1_M2[, -c(1:6)], 
#                                   nucleoplasmic_b2_M1_M2[, -c(1:6)], 
#                                   nucleoplasmic_b3_M1_M2[, -c(1:6)])
# 
# nucleoplasmic_all_M1_M2_barcodes <- rbind(nucleoplasmic_b1_M1_M2_barcodes, nucleoplasmic_b2_M1_M2_barcodes, nucleoplasmic_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(nucleoplasmic_all_M1_M2, nucleoplasmic_all_M1_M2_barcodes, type='all')
# p + labs(title='Nucleoplasmic_M1_M2 (All)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all_M1_M2, nucleoplasmic_all_M1_M2_barcodes, type='high')
# p + labs(title='Nucleoplasmic_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all_M1_M2, nucleoplasmic_all_M1_M2_barcodes, type='top')
# p + labs(title='Nucleoplasmic_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all_M1_M2, nucleoplasmic_all_M1_M2_barcodes, type='lpa')
# p + labs(title='Nucleoplasmic_M1_M2 (LPA induced genes)')
# 
# 
# 
# ## cytoplasmic all M1 M2 ####
# cytoplasmic_all_M1_M2 <- data.frame(cytoplasmic_b1_M1_M2[, 1:6], 
#                                       cytoplasmic_b1_M1_M2[, -c(1:6)], 
#                                       cytoplasmic_b2_M1_M2[, -c(1:6)], 
#                                       cytoplasmic_b3_M1_M2[, -c(1:6)])
# 
# cytoplasmic_all_M1_M2_barcodes <- rbind(cytoplasmic_b1_M1_M2_barcodes, cytoplasmic_b2_M1_M2_barcodes, cytoplasmic_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2, cytoplasmic_all_M1_M2_barcodes, type='all')
# p + labs(title='cytoplasmic_M1_M2 (All)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2[,-7], cytoplasmic_all_M1_M2_barcodes[-15,], type='all')
# p + labs(title='cytoplasmic_M1_M2 (All)')
# ggplotly(p)
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2, cytoplasmic_all_M1_M2_barcodes, type='high')
# p + labs(title='cytoplasmic_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2[,-7], cytoplasmic_all_M1_M2_barcodes[-15,], type='high')
# p + labs(title='cytoplasmic_M1_M2 (Highly expressed genes)')
# ggplotly(p)
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2, cytoplasmic_all_M1_M2_barcodes, type='top')
# p + labs(title='cytoplasmic_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2[,-7], cytoplasmic_all_M1_M2_barcodes[-15,], type='top')
# p + labs(title='cytoplasmic_M1_M2 (Highly expressed genes and high variance)')
# ggplotly(p)
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2, cytoplasmic_all_M1_M2_barcodes, type='lpa')
# p + labs(title='cytoplasmic_M1_M2 (LPA induced genes)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all_M1_M2[,-7], cytoplasmic_all_M1_M2_barcodes[-15,], type='lpa')
# p + labs(title='cytoplasmic_M1_M2 (LPA induced genes)')
# ggplotly(p)
# 
# 
# 
# ## polyA all M1 M2 #####
# polyA_all_M1_M2 <- data.frame(polyA_b1_M1_M2[, 1:6],
#                               polyA_b1_M1_M2[, -c(1:6)],
#                               polyA_b2_M1_M2[, -c(1:6)],
#                               polyA_b3_M1_M2[, -c(1:6)])
# 
# polyA_all_M1_M2_barcodes <- rbind(polyA_b1_M1_M2_barcodes, polyA_b2_M1_M2_barcodes, polyA_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(polyA_all_M1_M2, polyA_all_M1_M2_barcodes, type='all')
# p + labs(title='PolyA M1 M2 (All)')
# 
# p <- plot_pca_many_exp(polyA_all_M1_M2, polyA_all_M1_M2_barcodes, type='high')
# p + labs(title='PolyA M1 M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(polyA_all_M1_M2, polyA_all_M1_M2_barcodes, type='top')
# p + labs(title='PolyA M1 M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(polyA_all_M1_M2, polyA_all_M1_M2_barcodes, type='lpa')
# p + labs(title='PolyA M1 M2 (LPA induced genes)')
# 
# 
# 
# 
# ## all b1 ----
# b1_all_M1_M2 <- data.frame(chromatin_b1_M1_M2[, 1:6], 
#                                chromatin_b1_M1_M2[, -c(1:6)], 
#                                nucleoplasmic_b1_M1_M2[, -c(1:6)],
#                            cytoplasmic_b1_M1_M2[, -c(1:6)],
#                            polyA_b1_M1_M2[, -c(1:6)])
# 
# b1_all_M1_M2_barcodes <- rbind(chromatin_b1_M1_M2_barcodes, nucleoplasmic_b1_M1_M2_barcodes, cytoplasmic_b1_M1_M2_barcodes, polyA_b1_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2, b1_all_M1_M2_barcodes, type='all')
# p + labs(title='b1_M1_M2 (All)')
# ggplotly(p)
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2[,-55], b1_all_M1_M2_barcodes, type='all')
# p + labs(title='b1_M1_M2 (All)')
# ggplotly(p)
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2, b1_all_M1_M2_barcodes, type='high')
# p + labs(title='b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2[,-55], b1_all_M1_M2_barcodes, type='high')
# p + labs(title='b1_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2, b1_all_M1_M2_barcodes, type='top')
# p + labs(title='b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2[,-55], b1_all_M1_M2_barcodes, type='top')
# p + labs(title='b1_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2, b1_all_M1_M2_barcodes, type='lpa')
# p + labs(title='b1_M1_M2 (LPA induced genes)')
# 
# p <- plot_pca_many_exp_cpt(b1_all_M1_M2[,-55], b1_all_M1_M2_barcodes, type='lpa')
# p + labs(title='b1_M1_M2 (LPA induced genes)')
# 
# ## all b2 ----
# b2_all_M1_M2 <- data.frame(chromatin_b2_M1_M2[, 1:6], 
#                            chromatin_b2_M1_M2[, -c(1:6)], 
#                            nucleoplasmic_b2_M1_M2[, -c(1:6)],
#                            cytoplasmic_b2_M1_M2[, -c(1:6)],
#                            polyA_b2_M1_M2[, -c(1:6)])
# 
# b2_all_M1_M2_barcodes <- rbind(chromatin_b2_M1_M2_barcodes, nucleoplasmic_b2_M1_M2_barcodes, cytoplasmic_b2_M1_M2_barcodes, polyA_b2_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp_cpt(b2_all_M1_M2, b2_all_M1_M2_barcodes, type='all')
# p + labs(title='b2_M1_M2 (All)')
# 
# p <- plot_pca_many_exp_cpt(b2_all_M1_M2, b2_all_M1_M2_barcodes, type='high')
# p + labs(title='b2_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b2_all_M1_M2, b2_all_M1_M2_barcodes, type='top')
# p + labs(title='b2_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b2_all_M1_M2, b2_all_M1_M2_barcodes, type='lpa')
# p + labs(title='b2_M1_M2 (LPA induced genes)')
# 
# 
# 
# 
# ## all b3 ----
# b3_all_M1_M2 <- data.frame(chromatin_b3_M1_M2[, 1:6], 
#                                chromatin_b3_M1_M2[, -c(1:6)], 
#                                nucleoplasmic_b3_M1_M2[, -c(1:6)],
#                            cytoplasmic_b3_M1_M2[, -c(1:6)], 
#                            polyA_b3_M1_M2[, -c(1:6)])
# 
# b3_all_M1_M2_barcodes <- rbind(chromatin_b3_M1_M2_barcodes, nucleoplasmic_b3_M1_M2_barcodes, cytoplasmic_b3_M1_M2_barcodes, polyA_b3_M1_M2_barcodes )
# 
# p <- plot_pca_many_exp_cpt(b3_all_M1_M2, b3_all_M1_M2_barcodes, type='all')
# p + labs(title='b3_M1_M2 (All)')
# 
# p <- plot_pca_many_exp_cpt(b3_all_M1_M2, b3_all_M1_M2_barcodes, type='high')
# p + labs(title='b3_M1_M2 (Highly expressed genes)')
# 
# p <- plot_pca_many_exp_cpt(b3_all_M1_M2, b3_all_M1_M2_barcodes, type='top')
# p + labs(title='b3_M1_M2 (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp_cpt(b3_all_M1_M2, b3_all_M1_M2_barcodes, type='lpa')
# p + labs(title='b3_M1_M2 (LPA induced genes)')
# 
# 
# 
# 
# 
# ########################
# ### PCA all together ###
# ########################
# # chromatin ----
# chromatin_all <- data.frame(chromatin_b1_LPA_Naive[, 1:6], 
#                             chromatin_b1_LPA_Naive[, -c(1:6)], 
#                             chromatin_b2_LPA_Naive[, -c(1:6)],
#                             chromatin_b3_LPA_Naive[, -c(1:6)], 
#                             chromatin_b1_M1_M2[, -c(1:6)],
#                             chromatin_b2_M1_M2[, -c(1:6)],
#                             chromatin_b3_M1_M2[, -c(1:6)])
# 
# chromatin_all_barcodes <- rbind(chromatin_b1_LPA_Naive_barcodes, chromatin_b2_LPA_Naive_barcodes, chromatin_b3_LPA_Naive_barcodes, chromatin_b1_M1_M2_barcodes, chromatin_b2_M1_M2_barcodes, chromatin_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(chromatin_all, chromatin_all_barcodes, type='all')
# p + labs(title='Chromatin (All)')
# 
# p <- plot_pca_many_exp(chromatin_all, chromatin_all_barcodes, type='high')
# p + labs(title='Chromatin (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(chromatin_all, chromatin_all_barcodes, type='top')
# p + labs(title='Chromatin (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(chromatin_all, chromatin_all_barcodes, type='lpa')
# p + labs(title='Chromatin (LPA induced genes)')
# 
# 
# # nucleoplasmic ----
# nucleoplasmic_all <- data.frame(nucleoplasmic_b1_LPA_Naive[, 1:6], 
#                                 nucleoplasmic_b0_LPA_Naive[, -c(1:6)],
#                             nucleoplasmic_b1_LPA_Naive[, -c(1:6)], 
#                             nucleoplasmic_b3_LPA_Naive[, -c(1:6)], 
#                             nucleoplasmic_b1_M1_M2[, -c(1:6)],
#                             nucleoplasmic_b2_M1_M2[, -c(1:6)],
#                             nucleoplasmic_b3_M1_M2[, -c(1:6)])
# 
# nucleoplasmic_all_barcodes <- rbind(nucleoplasmic_b0_LPA_Naive_barcodes,nucleoplasmic_b1_LPA_Naive_barcodes, nucleoplasmic_b3_LPA_Naive_barcodes, nucleoplasmic_b1_M1_M2_barcodes, nucleoplasmic_b2_M1_M2_barcodes, nucleoplasmic_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(nucleoplasmic_all, nucleoplasmic_all_barcodes, type='all')
# p + labs(title='Nucleoplasmic (All)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all, nucleoplasmic_all_barcodes, type='high')
# p + labs(title='Nucleoplasmic (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all, nucleoplasmic_all_barcodes, type='top')
# p + labs(title='Nucleoplasmic (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(nucleoplasmic_all, nucleoplasmic_all_barcodes, type='lpa')
# p + labs(title='Nucleoplasmic (LPA induced genes)')
# 
# # cytoplasmic ----
# cytoplasmic_all <- data.frame(cytoplasmic_b1_LPA_Naive[, 1:6], 
#                               cytoplasmic_b0_LPA_Naive[, -c(1:6)], 
#                                 cytoplasmic_b1_LPA_Naive[, -c(1:6)], 
#                                 cytoplasmic_b2_LPA_Naive[, -c(1:6)], 
#                                 cytoplasmic_b3_LPA_Naive[, -c(1:6)], 
#                                 cytoplasmic_b1_M1_M2[, -c(1:6)],
#                                 cytoplasmic_b2_M1_M2[, -c(1:6)],
#                                 cytoplasmic_b3_M1_M2[, -c(1:6)])
# 
# cytoplasmic_all_barcodes <- rbind(cytoplasmic_b0_LPA_Naive_barcodes,cytoplasmic_b1_LPA_Naive_barcodes, cytoplasmic_b2_LPA_Naive_barcodes, cytoplasmic_b3_LPA_Naive_barcodes, cytoplasmic_b1_M1_M2_barcodes, cytoplasmic_b2_M1_M2_barcodes, cytoplasmic_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(cytoplasmic_all[,-374], cytoplasmic_all_barcodes[-342, ], type='all')
# p + labs(title='Cytoplasmic (All)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all[,-374], cytoplasmic_all_barcodes[-342, ], type='high')
# p + labs(title='Cytoplasmic (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all[,-374], cytoplasmic_all_barcodes[-342, ], type='top')
# p + labs(title='Cytoplasmic (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(cytoplasmic_all[,-374], cytoplasmic_all_barcodes[-342, ], type='lpa')
# p + labs(title='Cytoplasmic (LPA induced genes)')
# 
# # polyA ----
# polyA_all <- data.frame(polyA_b1_LPA_Naive[, 1:6],                                 
#                         polyA_b0_LPA_Naive[, -c(1:6)], 
#                                 polyA_b1_LPA_Naive[, -c(1:6)], 
#                         polyA_b2_LPA_Naive[, -c(1:6)], 
#                                 polyA_b3_LPA_Naive[, -c(1:6)], 
#                                 polyA_b1_M1_M2[, -c(1:6)],
#                                 polyA_b2_M1_M2[, -c(1:6)],
#                                 polyA_b3_M1_M2[, -c(1:6)])
# 
# polyA_all_barcodes <- rbind(polyA_b0_LPA_Naive_barcodes,polyA_b1_LPA_Naive_barcodes,polyA_b2_LPA_Naive_barcodes, polyA_b3_LPA_Naive_barcodes, polyA_b1_M1_M2_barcodes, polyA_b2_M1_M2_barcodes, polyA_b3_M1_M2_barcodes)
# 
# p <- plot_pca_many_exp(polyA_all, polyA_all_barcodes, type='all')
# p + labs(title='polyA (All)')
# 
# p <- plot_pca_many_exp(polyA_all, polyA_all_barcodes, type='high')
# p + labs(title='polyA (Highly expressed genes)')
# 
# p <- plot_pca_many_exp(polyA_all, polyA_all_barcodes, type='top')
# p + labs(title='polyA (Highly expressed genes and high variance)')
# 
# p <- plot_pca_many_exp(polyA_all, polyA_all_barcodes, type='lpa')
# p + labs(title='polyA (LPA induced genes)')
# 
# 
# 
