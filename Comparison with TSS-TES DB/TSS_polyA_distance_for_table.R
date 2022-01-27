# Set wd
setwd('Documents/Disk1/Documents/Ribo_Modeling_Reprocessed/Paper Final Run/5kb/')

# read sheet for polyA distance
gene_table <- read.table('Gene_for_polyA_TSS_distance.txt', header = T, sep = "\t", na.strings = c('?'))

# read bed file of CAGE peaks from: http://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz
Famtom5_CAGE_peak <- read.table(file = 'mm10_fair+new_CAGE_peaks_phase1and2.bed')
colnames(Famtom5_CAGE_peak) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'TSS start', 'TSS end', 'Color')

# read TSS from: http://reftss.clst.riken.jp/datafiles/current/mouse/refTSS_v3.1_mouse_coordinate.mm10.bed.gz
refTSS <- read.table(file = 'refTSS_v3.1_mouse_coordinate.mm10.bed')
colnames(refTSS) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'TSS start', 'TSS end', 'Color')

# read polyA from polyAsite: https://www.polyasite.unibas.ch/download/clusters/GRCm38-96/2-0/atlas.clusters.mm10.2-0.bed.gz
polyAsite <- read.table(file = 'polyAsite_atlas.clusters.mm10.2-0.bed')
colnames(polyAsite) <- c('chr', 'start', 'end', 'name', 'tpm', 'strand', 'percent samples', 'number of protocol', 'tpm2', 'cluster annotation', 'signal')
polyAsite$site <- as.numeric(sapply(polyAsite$name, function(x){unlist(strsplit(as.character(x),split = ':')[[1]][2])}))

# read polyA from Gencode
polyA_gencode <- read.table(file = 'gencode.vM14.polyAs.gtf')
polyA_gencode <- polyA_gencode[,c(1,3:5,7)]
polyA_gencode <- polyA_gencode[polyA_gencode$type == "polyA_site",]
colnames(polyA_gencode) <- c('chr', 'type', 'start', 'end', 'strand')

# Compare distance to TSS for annotated and observed => add colum to data table (NA if > 10kb)
position_of_closest_famtom5_annotated <- apply(gene_table, 1, function(x){peak_subset <- Famtom5_CAGE_peak[Famtom5_CAGE_peak$chr == paste0('chr',x['Chr']) & Famtom5_CAGE_peak$strand == x['Strand'],];
                                                                          distance <- abs(peak_subset$'TSS start'-as.numeric(x['Annotations...TSS']));
                                                                          closest <- peak_subset$'TSS start'[which.min(distance)]})

position_of_closest_famtom5_annotated <- unlist(lapply(position_of_closest_famtom5_annotated, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_famtom5_annotated <- abs(position_of_closest_famtom5_annotated - gene_table[,'Annotations...TSS'])

position_of_closest_famtom5_observed <- apply(gene_table, 1, function(x){peak_subset <- Famtom5_CAGE_peak[Famtom5_CAGE_peak$chr == paste0('chr',x['Chr']) & Famtom5_CAGE_peak$strand == x['Strand'],];
                                                                         distance <- abs(peak_subset$'TSS start'-as.numeric(x['Observed...TSS']));
                                                                         closest <- peak_subset$'TSS start'[which.min(distance)]})

position_of_closest_famtom5_observed <- unlist(lapply(position_of_closest_famtom5_observed, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_famtom5_observed <- abs(position_of_closest_famtom5_observed - gene_table[,'Observed...TSS'])




position_of_closest_refTSS_annotated <- apply(gene_table, 1, function(x){peak_subset <- refTSS[refTSS$chr == paste0('chr',x['Chr']) & refTSS$strand == x['Strand'],];
                                                                          distance <- abs(peak_subset$'TSS start'-as.numeric(x['Annotations...TSS']));
                                                                          closest <- peak_subset$'TSS start'[which.min(distance)]})

position_of_closest_refTSS_annotated <- unlist(lapply(position_of_closest_refTSS_annotated, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_refTSS_annotated <- abs(position_of_closest_refTSS_annotated - gene_table[,'Annotations...TSS'])

position_of_closest_refTSS_observed <- apply(gene_table, 1, function(x){peak_subset <- refTSS[refTSS$chr == paste0('chr',x['Chr']) & refTSS$strand == x['Strand'],];
                                                                        distance <- abs(peak_subset$'TSS start'-as.numeric(x['Observed...TSS']));
                                                                        closest <- peak_subset$'TSS start'[which.min(distance)]})

position_of_closest_refTSS_observed <- unlist(lapply(position_of_closest_refTSS_observed, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_refTSS_observed <- abs(position_of_closest_refTSS_observed - gene_table[,'Observed...TSS'])





position_of_closest_pA_annotated <- apply(gene_table, 1, function(x){peak_subset <- polyAsite[polyAsite$chr == x['Chr'] & polyAsite$strand == x['Strand'],];
                                                                      distance <- abs(peak_subset$site-as.numeric(x['Annotations...TES']));
                                                                      closest <- peak_subset$site[which.min(distance)]})

position_of_closest_pA_annotated <- unlist(lapply(position_of_closest_pA_annotated, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_pA_annotated <- abs(position_of_closest_pA_annotated - gene_table[,'Annotations...TES'])

position_of_closest_pA_observed <- apply(gene_table, 1, function(x){peak_subset <- polyAsite[polyAsite$chr == x['Chr'] & polyAsite$strand == x['Strand'],];
                                   distance <- abs(peak_subset$site-as.numeric(x['Observed...TES']));
                                   closest <- peak_subset$site[which.min(distance)]})

position_of_closest_pA_observed <- unlist(lapply(position_of_closest_pA_observed, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_pA_observed <- abs(position_of_closest_pA_observed - gene_table[,'Observed...TES'])




position_of_closest_pA_gencode_annotated <- apply(gene_table, 1, function(x){peak_subset <- polyA_gencode[polyA_gencode$chr == paste0('chr',x['Chr']) & polyA_gencode$strand == x['Strand'],];
                                                                            distance <- abs(peak_subset$start-as.numeric(x['Annotations...TES']));
                                                                            closest <- peak_subset$start[which.min(distance)]})

position_of_closest_pA_gencode_annotated <- unlist(lapply(position_of_closest_pA_gencode_annotated, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_pA_gencode_annotated <- abs(position_of_closest_pA_gencode_annotated - gene_table[,'Annotations...TES'])

position_of_closest_pA_gencode_observed <- apply(gene_table, 1, function(x){peak_subset <- polyA_gencode[polyA_gencode$chr == paste0('chr',x['Chr']) & polyA_gencode$strand == x['Strand'],];
                                                                            distance <- abs(peak_subset$start-as.numeric(x['Observed...TES']));
                                                                            closest <- peak_subset$start[which.min(distance)]})

position_of_closest_pA_gencode_observed <- unlist(lapply(position_of_closest_pA_gencode_observed, function(x){if(length(x) == 0){NA}else{x}}))
distance_of_closest_pA_gencode_observed <- abs(position_of_closest_pA_gencode_observed - gene_table[,'Observed...TES'])




# Make plots

library(ggplot2)
library(gridExtra)
dat <- data.frame(Position_TSS_Annotated = gene_table$'Annotations...TSS', 
                  Position_TSS_observed = gene_table$'Observed...TSS', 
                  Distance_TSS_between_observed_and_annotated = abs(gene_table$'Annotations...TSS' - gene_table$'Observed...TSS'), 
                  Distance_TES_between_observed_and_annotated = abs(gene_table$'Annotations...TES' - gene_table$'Observed...TES'), 
                  distance_of_closest_famtom5_annotated = 1+distance_of_closest_famtom5_annotated, 
                  distance_of_closest_famtom5_observed = 1+distance_of_closest_famtom5_observed,
                  distance_of_closest_refTSS_annotated = 1+distance_of_closest_refTSS_annotated, 
                  distance_of_closest_refTSS_observed = 1+distance_of_closest_refTSS_observed,
                  distance_of_closest_pA_annotated = 1+distance_of_closest_pA_annotated, 
                  distance_of_closest_pA_observed = 1+distance_of_closest_pA_observed,
                  distance_of_closest_pA_gencode_annotated = 1+distance_of_closest_pA_gencode_annotated, 
                  distance_of_closest_pA_gencode_observed = 1+distance_of_closest_pA_gencode_observed
)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p <- ggplot(dat, aes(x = distance_of_closest_famtom5_annotated, y = distance_of_closest_famtom5_observed, col= Distance_TSS_between_observed_and_annotated)) 
p <- p + geom_point(data = dat[dat$Distance_TSS_between_observed_and_annotated <= 50,], col='green3')
p <- p + geom_point(data = dat[dat$Distance_TSS_between_observed_and_annotated > 50,])
p <- p + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_y_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_color_continuous(trans='log10')
p <- p + theme_bw() + theme(aspect.ratio = 1)+labs(color= 'D(obs,annot)')

top <- ggplot(dat, aes(x = distance_of_closest_famtom5_annotated)) 
top <- top + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
top <- top + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
top <- top + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
top <- top + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
top

right <- ggplot(dat, aes(x = distance_of_closest_famtom5_observed)) 
right <- right + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
right <- right + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
right <- right + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
right <- right + theme_bw() + coord_flip() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
right

empty <- ggplot()+ theme(axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

legend <- g_legend(p)

g_main <- ggplotGrob(p + theme(legend.position="none"))
g_top <- ggplotGrob(top)
g_top$heights <- 0.25*g_main$heights
g_right <- ggplotGrob(right)
g_right$widths <- 0.25*g_main$widths

g <- rbind(g_top, g_main, size = "last")
g$widths <- g_main$widths 
  
g2 <- rbind(g_empty, g_right, size="last")
g2$widths <- g_right$widths

g3 <- cbind(g,g2, size="first")
g3$height <- unit.pmax(g$heights, g2$heights)
g3$widths

grid.newpage()
grid.draw(g3)

grid.arrange(arrangeGrob(g3, legend, nrow=1, widths=c(10,2), heights=10,respect = T))



# g <- ggplotGrob(p)
# 
# panel_id <- g$layout[g$layout$name == "panel",c("t","l")]
# g <- gtable_add_cols(g, unit(2,"cm"))
# 
# g <- gtable_add_grob(g, ggplotGrob(right), t = panel_id$t, l = ncol(g))
# g <- gtable_add_rows(g, unit(2,"cm"), 0)
# g <- gtable_add_grob(g, ggplotGrob(top),t = 1, l = panel_id$l)
# 
# grid.newpage()
# grid.draw(g)
# 
# ggMarginal(p, type = "histogram", binwidth=0.2)


p <- ggplot(dat, aes(x = distance_of_closest_refTSS_annotated, y = distance_of_closest_refTSS_observed, col= Distance_TSS_between_observed_and_annotated)) 
p <- p + geom_point(data = dat[dat$Distance_TSS_between_observed_and_annotated <= 50,], col='green3')
p <- p + geom_point(data = dat[dat$Distance_TSS_between_observed_and_annotated > 50,])
p <- p + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_y_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_color_continuous(trans='log10')
p <- p + theme_bw() + theme(aspect.ratio = 1)+labs(color= 'D(obs,annot)')

top <- ggplot(dat, aes(x = distance_of_closest_refTSS_annotated)) 
top <- top + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
top <- top + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
top <- top + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
top <- top + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
top

right <- ggplot(dat, aes(x = distance_of_closest_refTSS_observed)) 
right <- right + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
right <- right + geom_density(data = dat[dat$Distance_TSS_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
right <- right + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
right <- right + theme_bw() + coord_flip() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
right

empty <- ggplot()+ theme(axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

legend <- g_legend(p)

g_main <- ggplotGrob(p + theme(legend.position="none"))
g_top <- ggplotGrob(top)
g_top$heights <- 0.25*g_main$heights
g_right <- ggplotGrob(right)
g_right$widths <- 0.25*g_main$widths

g <- rbind(g_top, g_main, size = "last")
g$widths <- g_main$widths 

g2 <- rbind(g_empty, g_right, size="last")
g2$widths <- g_right$widths

g3 <- cbind(g,g2, size="first")
g3$height <- unit.pmax(g$heights, g2$heights)
g3$widths

grid.newpage()
grid.draw(g3)

grid.arrange(arrangeGrob(g3, legend, nrow=1, widths=c(10,2), heights=10,respect = T))



dat <- dat[!replace(gene_table$Disagree.with.Cufflinks.Strongly, is.na(gene_table$Disagree.with.Cufflinks.Strongly),FALSE),]

p <- ggplot(dat, aes(x = distance_of_closest_pA_annotated, y = distance_of_closest_pA_observed, col= Distance_TES_between_observed_and_annotated)) 
p <- p + geom_point(data = dat[dat$Distance_TES_between_observed_and_annotated <= 50,], col='green3')
p <- p + geom_point(data = dat[dat$Distance_TES_between_observed_and_annotated > 50,])
p <- p + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_y_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_color_continuous(trans='log10')
p <- p + theme_bw() + theme(aspect.ratio = 1)+labs(color= 'D(obs,annot)')

top <- ggplot(dat, aes(x = distance_of_closest_pA_annotated)) 
top <- top + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
top <- top + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
top <- top + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
top <- top + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
top

right <- ggplot(dat, aes(x = distance_of_closest_pA_observed)) 
right <- right + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
right <- right + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
right <- right + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
right <- right + theme_bw() + coord_flip() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
right

empty <- ggplot()+ theme(axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

legend <- g_legend(p)

g_main <- ggplotGrob(p + theme(legend.position="none"))
g_top <- ggplotGrob(top)
g_top$heights <- 0.25*g_main$heights
g_right <- ggplotGrob(right)
g_right$widths <- 0.25*g_main$widths

g <- rbind(g_top, g_main, size = "last")
g$widths <- g_main$widths 

g2 <- rbind(g_empty, g_right, size="last")
g2$widths <- g_right$widths

g3 <- cbind(g,g2, size="first")
g3$height <- unit.pmax(g$heights, g2$heights)
g3$widths

grid.newpage()
grid.draw(g3)

grid.arrange(arrangeGrob(g3, legend, nrow=1, widths=c(10,2), heights=10,respect = T))










p <- ggplot(dat, aes(x = distance_of_closest_pA_gencode_annotated, y = distance_of_closest_pA_gencode_observed, col= Distance_TES_between_observed_and_annotated)) 
p <- p + geom_point(data = dat[dat$Distance_TES_between_observed_and_annotated <= 50,], col='green3')
p <- p + geom_point(data = dat[dat$Distance_TES_between_observed_and_annotated > 50,])
p <- p + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_y_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_color_continuous(trans='log10')
p <- p + theme_bw() + theme(aspect.ratio = 1)+labs(color= 'D(obs,annot)')

top <- ggplot(dat, aes(x = distance_of_closest_pA_gencode_annotated)) 
top <- top + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
top <- top + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
top <- top + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
top <- top + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
top

right <- ggplot(dat, aes(x = distance_of_closest_pA_gencode_observed)) 
right <- right + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated <= 50,], fill='green3', alpha=0.3, bw=0.15)
right <- right + geom_density(data = dat[dat$Distance_TES_between_observed_and_annotated > 50,], fill = 'black', alpha=0.3, bw=0.15)
right <- right + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
right <- right + theme_bw() + coord_flip() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
right

empty <- ggplot()+ theme(axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

legend <- g_legend(p)

g_main <- ggplotGrob(p + theme(legend.position="none"))
g_top <- ggplotGrob(top)
g_top$heights <- 0.25*g_main$heights
g_right <- ggplotGrob(right)
g_right$widths <- 0.25*g_main$widths

g <- rbind(g_top, g_main, size = "last")
g$widths <- g_main$widths 

g2 <- rbind(g_empty, g_right, size="last")
g2$widths <- g_right$widths

g3 <- cbind(g,g2, size="first")
g3$height <- unit.pmax(g$heights, g2$heights)
g3$widths

grid.newpage()
grid.draw(g3)

grid.arrange(arrangeGrob(g3, legend, nrow=1, widths=c(10,2), heights=10,respect = T))

# Add to table
gene_table$position_of_closest_famtom5_annotated <- position_of_closest_famtom5_annotated
gene_table$position_of_closest_famtom5_observed <- position_of_closest_famtom5_observed
gene_table$position_of_closest_refTSS_annotated <- position_of_closest_refTSS_annotated
gene_table$position_of_closest_refTSS_observed <- position_of_closest_refTSS_observed
gene_table$position_of_closest_polyAsite_annotated <- position_of_closest_pA_annotated
gene_table$position_of_closest_polyAsite_observed <- position_of_closest_pA_observed
gene_table$position_of_closest_polyAsite_annotated_gencode <- position_of_closest_pA_gencode_annotated
gene_table$position_of_closest_polyAsite_observed_gencode <- position_of_closest_pA_gencode_observed

write.csv(gene_table,file = 'table_with_distances.csv')

# gene_length
length <- read.table('exon_counts_caRNAb1_200bp.txt', header = T)
p <- ggplot(length, aes(x=Length)) + geom_density(fill='gray', alpha=0.3)
p <- p + theme_bw() + theme(aspect.ratio = 0.25)
p

p <- ggplot(length, aes(x=Length)) + geom_histogram(fill='gray')
p <- p + theme_bw() + theme(aspect.ratio = 0.25)
p

length2 <- read.table('exon_counts_caRNAb1_50bp.txt', header = T)
p <- ggplot(length, aes(x=Length)) + geom_density(fill='gray', alpha=0.3)
p <- p + theme_bw() + theme(aspect.ratio = 0.25)
p

p <- ggplot(data.frame(length_200bp=length$Length,length_50bp=length2$Length[match(length$Geneid, length2$Geneid)]), aes(x=length_200bp, y=length_50bp)) + geom_point()
p <- p + theme_bw() + theme(aspect.ratio = 0.25)
p

