# read sheet for polyA distance
gene_table <- read.table('Gene_for_polyA_TSS_distance.txt', header = T, sep = "\t", na.strings = c('?'))

# read TSS from: http://reftss.clst.riken.jp/datafiles/current/mouse/refTSS_v3.1_mouse_coordinate.mm10.bed.gz
refTSS <- read.table(file = 'refTSS_v3.1_mouse_coordinate.mm10.bed.gz')
colnames(refTSS) <- c('chr', 'start', 'end', 'name', 'score', 'strand', 'TSS start', 'TSS end', 'Color')

# read polyA from polyAsite: https://www.polyasite.unibas.ch/download/clusters/GRCm38-96/2-0/atlas.clusters.mm10.2-0.bed.gz
polyAsite <- read.table(file = 'polyAsite_atlas.clusters.mm10.2-0.bed')
colnames(polyAsite) <- c('chr', 'start', 'end', 'name', 'tpm', 'strand', 'percent samples', 'number of protocol', 'tpm2', 'cluster annotation', 'signal')
polyAsite$site <- as.numeric(sapply(polyAsite$name, function(x){unlist(strsplit(as.character(x),split = ':')[[1]][2])}))

# Compare distance to TSS for annotated and observed
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

# Compare distance to TES for annotated and observed
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

# Make plots
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
dat <- data.frame(Position_TSS_Annotated = gene_table$'Annotations...TSS', 
                  Position_TSS_observed = gene_table$'Observed...TSS', 
                  Distance_TSS_between_observed_and_annotated = 1+abs(gene_table$'Annotations...TSS' - gene_table$'Observed...TSS'), 
                  Distance_TES_between_observed_and_annotated =1+abs(gene_table$'Annotations...TES' - gene_table$'Observed...TES'), 
                  distance_of_closest_refTSS_annotated = 1+distance_of_closest_refTSS_annotated, 
                  distance_of_closest_refTSS_observed = 1+distance_of_closest_refTSS_observed,
                  distance_of_closest_pA_annotated = 1+distance_of_closest_pA_annotated, 
                  distance_of_closest_pA_observed = 1+distance_of_closest_pA_observed,
                  disagree_with_cufflinks = gene_table$Disagree.with.Cufflinks.Strongly
)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

log10p1_trans <- function(){ trans_new(name = "log10p1",
                              transform = function(x){log10(1+x)},
                              inverse = function(x){10^x-1})
}

# refTSS
dat_ggplot <- dat[!replace(dat$disagree_with_cufflinks, is.na(dat$disagree_with_cufflinks),FALSE),]
p <- ggplot(dat_ggplot, aes(x = distance_of_closest_refTSS_annotated, y = distance_of_closest_refTSS_observed, col= Distance_TSS_between_observed_and_annotated)) 
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_y_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_color_continuous(trans='log10',  breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000), limits=c(1,2001), oob=squish)
p <- p + theme_bw() + theme(aspect.ratio = 1)
p <- p + labs(x="Distance of closest refTSS to annotated TSS", y="Distance of closest refTSS to observed TSS", col="Distance between\nobserved and\nannotated TSS (bp)")

top <- ggplot(dat_ggplot, aes(x = distance_of_closest_refTSS_annotated)) 
top <- top + geom_density(fill="black", alpha=0.3, bw=0.15)
top <- top + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
top <- top + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
top

right <- ggplot(dat_ggplot, aes(x = distance_of_closest_refTSS_observed)) 
right <- right + geom_density(fill = 'black', alpha=0.3, bw=0.15)
right <- right + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
right <- right + theme_bw() + coord_flip() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
right

empty <- ggplot()+ theme(axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
g_empty <- ggplotGrob(empty)

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

grid.arrange(arrangeGrob(g3, legend, nrow=1, widths=c(10,3), heights=10,respect = T))


# PolyASite
dat_ggplot <- dat[!replace(dat$disagree_with_cufflinks, is.na(dat$disagree_with_cufflinks),FALSE),]
p <- ggplot(dat_ggplot, aes(x = distance_of_closest_pA_annotated, y = distance_of_closest_pA_observed, col= Distance_TES_between_observed_and_annotated)) 
p <- p + geom_point()
p <- p + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_y_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
p <- p + scale_color_continuous(trans='log10',  breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000), limits=c(1,2001), oob=squish)
p <- p + theme_bw() + theme(aspect.ratio = 1)+labs(color= 'Distance between\nobserved and\nannotated TES (bp)', y='Distance of closest PolyASite to observed TES', x='Distance of closest PolyAsite to annotated TES')

top <- ggplot(dat_ggplot, aes(x = distance_of_closest_pA_annotated)) 
top <- top + geom_density(fill="black", alpha=0.3, bw=0.15)
top <- top + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
top <- top + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
top

right <- ggplot(dat, aes(x = distance_of_closest_pA_observed)) 
right <- right + geom_density(fill="black", alpha=0.3, bw=0.15)
right <- right + scale_x_continuous(trans='log10', breaks = 1+c(0,10,50,100,200,500,1000),labels = c(0,10,50,100,200,500,1000))
right <- right + theme_bw() + coord_flip() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
right

empty <- ggplot()+ theme(axis.ticks=element_blank(), panel.background=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
g_empty <- ggplotGrob(empty)

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

grid.arrange(arrangeGrob(g3, legend, nrow=1, widths=c(10,3), heights=10,respect = T))
