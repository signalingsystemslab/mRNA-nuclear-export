# check RBP
setwd('/home/diane/Documents/UCLA/Ribo_New/RBP')
manual_start_end <- read.csv('Curated list_naive.csv', header = T, sep = '\t')
transcripts <- unname(unlist(sapply(as.character(manual_start_end$Isoforms), function(x){strsplit(x, ', ')[[1]]})))

library(biomaRt)
ensembl89 <- useMart(biomart = 'ensembl', host="may2017.archive.ensembl.org")
mm10_89 <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl89)

gene_UTR_start_end <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 
                                           'chromosome_name', 'strand', 'transcript_start', 'transcript_end', 
                                           '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end', 
                                           'genomic_coding_start', 'genomic_coding_end'),
                            filters = 'ensembl_transcript_id', 
                            values = substr(transcripts,1,18), 
                            mart = mm10_89
)
head(gene_UTR_start_end) # multiple rows for 5'utr or 3'utr because it also has introns

server <- "https://may2017.rest.ensembl.org"
library(httr)

# get relative transport for all genes
library(compiler)
source('../Modeling/All_genes_last_5kbs/Model1/Hoffman2/functions_weighting_for_smoothing_include_negbinom_error_in_model_v13.R')

param_best_all <- c()
rel_transport_all <- c()
for ( gene in manual_start_end$Geneid ){
  gene_name <- manual_start_end$Name[manual_start_end$Geneid == gene]
  file <- paste0('../Modeling/All_genes_last_5kbs/Model1/Hoffman2/Optimisation/Modeling_v13/Optimisation/BFGS/naive/param_bfgs_1000ri_merged_', gene_name, '_Finished.Rdata')
  if(file.exists(file)){
    load(file)
    
    # get best parameter
    par <- bfgs_merged[['all']]$parameter[which.min(bfgs_merged[['all']]$value),]
    
    # convert parameter to original
    par_old <- parameters_convert_cp(par[-c(1:4)], log_in = T, log_out = T, direction = -1)

    # calculate relative transport
    rel_transport_all <- rbind(rel_transport_all, data.frame(GeneId=gene, Name=gene_name, Transport = par_old[1]+par_old[3]-par_old[2]))
    param_best_all <- rbind(param_best_all, data.frame(GeneID=gene, Name=gene_name, "k1' - ca -> " = par_old[1], "k2 - -> np" = par_old[2], "k2' - np -> " = par_old[3], "kdeg" = par_old[4], Transport  = par_old[1]+par_old[3]-par_old[2]))
  }
}

# dna_to_rna <- function(dna){
#   rna <- gsub("a", "u", dna)
#   rna <- gsub("A", "U", rna)
#   rna <- gsub("t", "a", rna)
#   rna <- gsub("T", "A", rna)
#   rna <- gsub("c", "@", rna)
#   rna <- gsub("g", "c", rna)
#   rna <- gsub("@", "g", rna)
#   rna <- gsub("C", "@", rna)
#   rna <- gsub("G", "C", rna)
#   rna <- gsub("@", "G", rna)
#  
#   return(rna) 
# }
# 
# 
# # order gene by transport value for all (more -> less)
# gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]
# 
# # for Meme get sequence of all transcripts and add weight in fasta
# MEME_fasta_5utr <- MEME_fasta_3utr <- c()
# MEME_fasta_weight <- c()
# 
# for (gene in gene_id_ordered){
#   dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
#   strand <- unique(dat_g$strand)
#   chr <- unique(dat_g$chromosome_name)
#   
#   rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
#   
#   gene_name <- unique(dat_g$external_gene_name)
#   print(gene_name)
#   
#   seq_5utr <- seq_3utr <- c()
#   
#   TSS_manual <- manual_start_end$TSS[manual_start_end$Geneid == gene ]
#   TES_manual <- manual_start_end$TES[manual_start_end$Geneid == gene ]
#   
#   MEME_fasta_weight <- c(MEME_fasta_weight, rep(1/length(unique(dat_g$ensembl_transcript_id)), length(unique(dat_g$ensembl_transcript_id))))
#   
#   # get utr start and end for each transcript
#   for (transcript in unique(dat_g$ensembl_transcript_id)){
#     MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0('> ', transcript , ' | ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
#     MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0('> ', transcript , ' | ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
#     
#     dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
#     
#     start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
#     end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
#     start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
#     end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
#     
#     # print(length(start_5utr_tmp))
#     # print(length(start_3utr_tmp))
#     
#     # all_sequences <- c(all_sequence, list())
#     seq_5utr <- c()
#     seq_3utr <- c()
#     
#     # In GET take - strand because RNA is reversed, to get real motifs
#     
#     if(strand == 1){
#       if( length(start_5utr_tmp) > 0 & !is.na(TSS_manual)){
#         for (i in 1:length(start_5utr_tmp)){
#           if(end_5utr_tmp[i] >= TSS_manual){
#             ext5 <- paste0("/sequence/region/mouse/", chr, ":",max(TSS_manual,start_5utr_tmp[i]),'..',end_5utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_5utr <- c(seq_5utr, content(r)$seq)
#           }
#         }
#       }
#       if( length(start_3utr_tmp) >0 & !is.na(TES_manual)){
#         for (i in 1:length(start_3utr_tmp)){
#           if(start_3utr_tmp[i] <= TES_manual){
#             ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',min(end_3utr_tmp[i], TES_manual),':', strand,"?mask=soft")
#             r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_3utr <- c(seq_3utr, content(r)$seq)
#           }
#         }
#       }
#     }else if (strand == -1){
#       if( length(start_5utr_tmp) >0 & !is.na(TSS_manual)){
#         for (i in 1:length(start_5utr_tmp)){
#           if(start_5utr_tmp[i] <= TSS_manual){
#             ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',min(TSS_manual,end_5utr_tmp[i]),':', strand,"?mask=soft")
#             r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_5utr <- c(seq_5utr, content(r)$seq)
#           }
#         }
#       }
#       if( length(start_3utr_tmp) > 0 & !is.na(TES_manual) ){
#         for (i in 1:length(start_3utr_tmp)){
#           if(end_3utr_tmp[i] >= TES_manual){
#             ext3 <- paste0("/sequence/region/mouse/", chr, ":", max(start_3utr_tmp[i], TES_manual),'..',end_3utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_3utr <- c(seq_3utr, content(r)$seq)
#           }
#         }
#       }
#     }
#     
#     # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
#     # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
#     MEME_fasta_3utr <- c(MEME_fasta_3utr, dna_to_rna(paste0(seq_3utr, collapse = '')))
#     MEME_fasta_5utr <- c(MEME_fasta_5utr, dna_to_rna(paste0(seq_5utr, collapse = '')))
#   }
# }
# 
# # remove zero length from file and weight
# zero_length_3utr <- which(sapply(MEME_fasta_3utr, function(x){nchar(x)==0}))
# MEME_fasta_3utr_to_write <- MEME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
# MEME_fasta_weight_3utr <- MEME_fasta_weight[-c(zero_length_3utr/2)]
#   
# zero_length_5utr <- which(sapply(MEME_fasta_5utr, function(x){nchar(x)==0}))
# MEME_fasta_5utr_to_write <- MEME_fasta_3utr[-c(zero_length_5utr-1,zero_length_5utr)]
# MEME_fasta_weight_5utr <- MEME_fasta_weight[-c(zero_length_5utr/2)]
# 
# ## Add weigh to MEME fasta
# MEME_fasta_3utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_3utr), collapse = ' '), MEME_fasta_3utr_to_write)
# MEME_fasta_5utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_5utr), collapse = ' '), MEME_fasta_5utr_to_write)
# 
# ## write to file
# write.table(x = MEME_fasta_3utr_to_write, file = 'meme/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# write.table(x = MEME_fasta_5utr_to_write, file = 'meme/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# 
# for (i in 1:(length(gene_id_ordered)-1)){
#   first_gene_less <- gene_id_ordered[i+1]
#   
#   first_gene_less_idx <- which.min(sapply(MEME_fasta_3utr, function(x){grep(first_gene_less, x)}))
#   
#   MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr[1:(first_gene_less_idx-1)]
#   MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_weight[1:(first_gene_less_idx-1)/2]
#   
#   zero_length_3utr_more <- which(sapply(MEME_fasta_3utr_to_write_more, function(x){nchar(x)==0}))
#   MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
#   MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_3utr_to_write_more_weight[-c(zero_length_3utr_more/2)]
# 
#   MEME_fasta_3utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_more_weight), collapse = ' '), MEME_fasta_3utr_to_write_more)
#   
#   
#   MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr[(first_gene_less_idx):length(MEME_fasta_3utr)]
#   MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
#   
#   zero_length_3utr_less <- which(sapply(MEME_fasta_3utr_to_write_less, function(x){nchar(x)==0}))
#   MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
#   MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_3utr_to_write_less_weight[-c(zero_length_3utr_less/2)]
#   
#   MEME_fasta_3utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_less_weight), collapse = ' '), MEME_fasta_3utr_to_write_less)
#   
#   MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr[1:(first_gene_less_idx-1)]
#   MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_weight[1:(first_gene_less_idx-1)/2]
#   
#   zero_length_5utr_more <- which(sapply(MEME_fasta_5utr_to_write_more, function(x){nchar(x)==0}))
#   MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
#   MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_5utr_to_write_more_weight[-c(zero_length_5utr_more/2)]
#   
#   MEME_fasta_5utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_more_weight), collapse = ' '), MEME_fasta_5utr_to_write_more)
#   
#   
#   MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr[(first_gene_less_idx):length(MEME_fasta_5utr)]
#   MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
#   
#   zero_length_5utr_less <- which(sapply(MEME_fasta_5utr_to_write_less, function(x){nchar(x)==0}))
#   MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
#   MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_5utr_to_write_less_weight[-c(zero_length_5utr_less/2)]
#   
#   MEME_fasta_5utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_less_weight), collapse = ' '), MEME_fasta_5utr_to_write_less)
#   
#   write.table(x = MEME_fasta_3utr_to_write_more, file = paste0('meme/3utr/meme_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_3utr_to_write_less, file = paste0('meme/3utr/meme_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_5utr_to_write_more, file = paste0('meme/5utr/meme_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_5utr_to_write_less, file = paste0('meme/5utr/meme_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   
# }
# 
# #### merge interval function for ame
# merge_intervals <- function(intervals, strand){
#   intervals_merged <- intervals[intervals$transcript == 1, c("start","end")]
#   merge_delim <- rep("",nrow(intervals_merged))
#   
#   for (i in 2:length(unique(intervals$transcript))){
#     intervals_tmp <- intervals[intervals$transcript == i, ]
#     
#     if(strand == 1){
#       # remove exons of 2nd transcript that ends before 1st exons of 1st transcript starts
#       intervals_tmp <- intervals_tmp[!(intervals_tmp$end < intervals_merged$start[1]),]
#       
#       # remove exons of 1st transcript that ends before 1st exons of 2nd transcript starts
#       merge_delim <- merge_delim[!(intervals_merged$end < intervals_tmp$start[1])]
#       intervals_merged <- intervals_merged[!(intervals_merged$end < intervals_tmp$start[1]),]
#       
#       # remove exons of 2nd transcript that start after last exons of 1st transcript ends
#       intervals_tmp <- intervals_tmp[!(intervals_tmp$start > intervals_merged$end[nrow(intervals_merged)]),]
#       
#       # remove exons of 1st transcript that start after last exons of 2nd transcript ends
#       merge_delim <- merge_delim[!(intervals_merged$start > intervals_tmp$end[nrow(intervals_tmp)])]
#       intervals_merged <- intervals_merged[!(intervals_merged$start > intervals_tmp$end[nrow(intervals_tmp)]),]
#     }else if (strand == -1){
#       # remove exons of 2nd transcript that ends before 1st exons of 1st transcript starts
#       intervals_tmp <- intervals_tmp[!(intervals_tmp$end < intervals_merged$start[nrow(intervals_merged)]),]
#       
#       # remove exons of 1st transcript that ends before 1st exons of 2nd transcript starts
#       merge_delim <- merge_delim[!(intervals_merged$end < intervals_tmp$start[nrow(intervals_tmp)])]
#       intervals_merged <- intervals_merged[!(intervals_merged$end < intervals_tmp$start[nrow(intervals_tmp)]),]
#       
#       # remove exons of 2nd transcript that start after last exons of 1st transcript ends
#       intervals_tmp <- intervals_tmp[!(intervals_tmp$start > intervals_merged$end[1]),]
#       
#       # remove exons of 1st transcript that start after last exons of 2nd transcript ends
#       merge_delim <- merge_delim[!(intervals_merged$start > intervals_tmp$end[1])]
#       intervals_merged <- intervals_merged[!(intervals_merged$start > intervals_tmp$end[1]),]
#     }
#     
#     # remove exons of 2nd transcript that are not overlapping with any exons of 1st transcript
#     overlapping <- c()
#     for ( j in 1:nrow(intervals_tmp)){
#       overlapping <- c(overlapping, any(intervals_tmp$start[j] >=  intervals_merged$start & intervals_tmp$start[j] <=  intervals_merged$end) | any(intervals_tmp$end[j] >=  intervals_merged$start & intervals_tmp$end[j] <=  intervals_merged$end))
#     } 
#     intervals_tmp <- intervals_tmp[overlapping,]
#     
#     # remove exons of 1st transcript that are not overlapping with any exons of 2nd transcript
#     overlapping <- c()
#     for ( j in 1:nrow(intervals_merged)){
#       overlapping <- c(overlapping, any(intervals_merged$start[j] >=  intervals_tmp$start & intervals_merged$start[j] <=  intervals_tmp$end) | any(intervals_merged$end[j] >=  intervals_tmp$start & intervals_merged$end[j] <=  intervals_tmp$end))
#     } 
#     intervals_merged <- intervals_merged[overlapping,]
#     merge_delim <- merge_delim[overlapping]
#     
#     # merge overlapping exons
#     
#     for (j in 1:nrow(intervals_tmp)){
#       start_in <- intervals_tmp$start[j] >=  intervals_merged$start & intervals_tmp$start[j] <=  intervals_merged$end
#       end_in <- intervals_tmp$end[j] >=  intervals_merged$start & intervals_tmp$end[j] <=  intervals_merged$end
#       
#       split <- FALSE
#       
#       if(j != nrow(intervals_tmp)){
#         start_in_next <- intervals_tmp$start[j+1] >=  intervals_merged$start & intervals_tmp$start[j+1] <=  intervals_merged$end
#         if(any(start_in_next == end_in)){ # => split 
#           split <- TRUE
#         }
#       }
# 
#       # if exon j start of 2nd transcript within other exon then modify the 1st transcript
#       if (any(start_in)){
#         if(which(start_in) != 1 & intervals_merged$start[which(start_in)] != intervals_tmp$start[j]){
#           merge_delim[which(start_in)-1] <- paste0(rep("N",10),collapse="")
#         }
#         intervals_merged$start[which(start_in)] <- intervals_tmp$start[j]
#       }
#         
#       
#       if (any(end_in)){
#         if(!split){
#           if(which(end_in) != nrow(intervals_merged) & intervals_merged$end[which(end_in)] != intervals_tmp$end[j]){
#             merge_delim[which(end_in)] <- paste0(rep("N",10),collapse="")
#           }
#           intervals_merged$end[which(end_in)] <- intervals_tmp$end[j]
#         }else{
#           intervals_merged <- rbind(intervals_merged[1:which(end_in),], intervals_merged[which(end_in):nrow(intervals_merged),])
#           if(which(end_in) !=1 ){
#             merge_delim <- c(merge_delim[1:(which(end_in)-1)],"", merge_delim[(which(end_in)-1):length(merge_delim)])
#           }else{
#             merge_delim <- c("", merge_delim[(which(end_in)-1):length(merge_delim)])
#           }
#           intervals_merged$end[which(end_in)] <- intervals_tmp$end[j]
#           intervals_merged$start[which(end_in)+1] <- intervals_tmp$end[j]+1
#         }
#       }
#     }
#   }
#   warnings()
#   res <- list(intervals_merged,merge_delim)
# }
# 
# 
# 
# 
# #### for Ame try to put parenthesis to merge trasncripts
# AME_fasta_5utr <- AME_fasta_3utr <- c()
# 
# for (gene in gene_id_ordered){
#   dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
#   strand <- unique(dat_g$strand)
#   chr <- unique(dat_g$chromosome_name)
#   
#   rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
#   
#   gene_name <- unique(dat_g$external_gene_name)
#   print(gene_name)
#   
#   AME_fasta_3utr <- c(AME_fasta_3utr, paste0('> ', gene , ' | ', gene_name, ' | ', strand, ' | ', rel_transport))
#   AME_fasta_5utr <- c(AME_fasta_5utr, paste0('> ', gene , ' | ', gene_name, ' | ', strand, ' | ', rel_transport))
#   
#   TSS_manual <- manual_start_end$TSS[manual_start_end$Geneid == gene ]
#   TES_manual <- manual_start_end$TES[manual_start_end$Geneid == gene ]
#   
#   start_5utr <- start_3utr <- c()
#   end_5utr <- end_3utr <- c()
#   
#   # get utr start and end for each transcript
#   for (transcript in unique(dat_g$ensembl_transcript_id)){
#     dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
#     
#     start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
#     end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
#     start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
#     end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
#     
#     start_5utr <- c(start_5utr, list(start_5utr_tmp))
#     end_5utr <- c(end_5utr, list(end_5utr_tmp))
#     start_3utr <- c(start_3utr, list(start_3utr_tmp))
#     end_3utr <- c(end_3utr, list(end_3utr_tmp))
#   }
#   
#   # get common intervals and get sequences and merge with 10 N (max motif width for RBP 8) when inconsistent between transcripts
#   if (length(unique(dat_g$ensembl_transcript_id)) > 1){
#     if(strand == 1){
#       dec <- F
#     }else if (strand == -1){
#       dec <- T
#     }
#     
#     if(any(unlist(lapply(start_3utr, function(l){length(l)})) == 0)){
#       start_3utr_merged <- c()
#       end_3utr_merged <- c()
#       merge_delim_3utr <- c()
#     }else{
#       intervals_3utr <- do.call("rbind", lapply(1:length(start_3utr), function(i){data.frame(start=start_3utr[[i]], end=end_3utr[[i]], "transcript"=i)}))
#       intervals_3utr <- intervals_3utr[order(intervals_3utr$start, decreasing = dec),]
#       
#       tmp <- merge_intervals(intervals_3utr, strand)
#       
#       start_3utr_merged <- tmp[[1]]$start
#       end_3utr_merged <- tmp[[1]]$end
#       merge_delim_3utr <- tmp[[2]]
#     }
#     
#     if(any(unlist(lapply(start_5utr, function(l){length(l)})) == 0)){
#       start_5utr_merged <- c()
#       end_5utr_merged <- c()
#       merge_delim_5utr <- c()
#     }else{
#       intervals_5utr <- do.call("rbind", lapply(1:length(start_5utr), function(i){data.frame(start=start_5utr[[i]], end=end_5utr[[i]], "transcript"=i)}))
#       intervals_5utr <- intervals_5utr[order(intervals_5utr$start, decreasing = dec),]
#       
#       tmp <- merge_intervals(intervals_5utr, strand)
#       
#       start_5utr_merged <- tmp[[1]]$start
#       end_5utr_merged <- tmp[[1]]$end
#       merge_delim_5utr <- tmp[[2]]
#     }
#       
#   }else{
#     start_5utr_merged <- unlist(start_5utr)
#     end_5utr_merged <- unlist(end_5utr)
#     start_3utr_merged <- unlist(start_3utr)
#     end_3utr_merged <- unlist(end_3utr)
#     
#     merge_delim_5utr <- rep("",length(start_5utr_merged))
# 
#     merge_delim_3utr <- rep("",length(start_3utr_merged))
#     
#   }
#   
# 
# 
#   seq_5utr <- c()
#   seq_3utr <- c()
#   # get sequences and merge with 10 N (max motif width for RBP 8)
#   # In GET take - strand because RNA is reversed, to get real motifs
#   if(strand == 1){
#     if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
#       for (i in 1:length(start_5utr_merged)){
#         if(start_5utr_merged[i] >= TSS_manual){
#           ext5 <- paste0("/sequence/region/mouse/", chr, ":",max(TSS_manual,start_5utr_merged[i]),'..',end_5utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_5utr <- c(seq_5utr, content(r)$seq)
#         }
#       }
#     }
#     if( length(start_3utr_merged) >0 & !is.na(TES_manual)){
#       for (i in 1:length(start_3utr_merged)){
#         if(start_3utr_merged[i] <= TES_manual){
#           ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..',min(end_3utr_merged[i], TES_manual),':', strand,"?mask=soft")
#           r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_3utr <- c(seq_3utr, content(r)$seq)
#         }
#       }
#     }
#   }else if (strand == -1){
#     if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
#       for (i in 1:length(start_5utr_merged)){
#         if(start_5utr_merged[i] <= TSS_manual){
#           ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',min(TSS_manual,end_5utr_merged[i]),':', strand,"?mask=soft")
#           r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_5utr <- c(seq_5utr, content(r)$seq)
#         }
#       }
#     }
#     if( length(start_3utr_merged) > 0 & !is.na(TES_manual) ){
#       for (i in 1:length(start_3utr_merged)){
#         if(end_3utr_merged[i] >= TES_manual){
#           ext3 <- paste0("/sequence/region/mouse/", chr, ":", max(start_3utr_merged[i], TES_manual),'..',end_3utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_3utr <- c(seq_3utr, content(r)$seq)
#         }
#       }
#     }
#   }
#   
#   # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
#   # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
#   AME_fasta_3utr <- c(AME_fasta_3utr, dna_to_rna(paste0(seq_3utr, merge_delim_3utr, sep="", collapse = '')))
#   AME_fasta_5utr <- c(AME_fasta_5utr, dna_to_rna(paste0(seq_5utr, merge_delim_5utr, sep="", collapse = '')))
#   
# }
# 
# # remove zero length from file and weight
# zero_length_3utr <- which(sapply(AME_fasta_3utr, function(x){nchar(x)==0}))
# AME_fasta_3utr_to_write <- AME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
# 
# zero_length_5utr <- which(sapply(AME_fasta_5utr, function(x){nchar(x)==0}))
# AME_fasta_5utr_to_write <- AME_fasta_3utr[-c(zero_length_5utr-1,zero_length_5utr)]
# 
# ## write to file
# write.table(x = AME_fasta_3utr_to_write, file = 'ame/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# write.table(x = AME_fasta_5utr_to_write, file = 'ame/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# 
# for (i in 1:(length(gene_id_ordered)-1)){
#   first_gene_less <- gene_id_ordered[i+1]
#   
#   first_gene_less_idx <- which.min(sapply(AME_fasta_3utr, function(x){grep(first_gene_less, x)}))
#   
#   AME_fasta_3utr_to_write_more <- AME_fasta_3utr[1:(first_gene_less_idx-1)]
# 
#   zero_length_3utr_more <- which(sapply(AME_fasta_3utr_to_write_more, function(x){nchar(x)==0}))
#   AME_fasta_3utr_to_write_more <- AME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
#   
#   AME_fasta_3utr_to_write_less <- AME_fasta_3utr[(first_gene_less_idx):length(AME_fasta_3utr)]
# 
#   zero_length_3utr_less <- which(sapply(AME_fasta_3utr_to_write_less, function(x){nchar(x)==0}))
#   AME_fasta_3utr_to_write_less <- AME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
# 
#   AME_fasta_5utr_to_write_more <- AME_fasta_5utr[1:(first_gene_less_idx-1)]
# 
#   zero_length_5utr_more <- which(sapply(AME_fasta_5utr_to_write_more, function(x){nchar(x)==0}))
#   AME_fasta_5utr_to_write_more <- AME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
# 
#   AME_fasta_5utr_to_write_less <- AME_fasta_5utr[(first_gene_less_idx):length(AME_fasta_5utr)]
# 
#   zero_length_5utr_less <- which(sapply(AME_fasta_5utr_to_write_less, function(x){nchar(x)==0}))
#   AME_fasta_5utr_to_write_less <- AME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
# 
#   write.table(x = AME_fasta_3utr_to_write_more, file = paste0('ame/3utr/ame_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_3utr_to_write_less, file = paste0('ame/3utr/ame_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_5utr_to_write_more, file = paste0('ame/5utr/ame_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_5utr_to_write_less, file = paste0('ame/5utr/ame_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   
# }
# 
# ### V2 make sure to extend if TES further away ###-------------
# 
# # order gene by transport value for all (more -> less)
# gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]
# 
# # for Meme get sequence of all transcripts and add weight in fasta
# MEME_fasta_5utr <- MEME_fasta_3utr <- c()
# MEME_fasta_weight <- c()
# 
# for (gene in gene_id_ordered){
#   dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
#   strand <- unique(dat_g$strand)
#   chr <- unique(dat_g$chromosome_name)
#   
#   rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
#   
#   gene_name <- unique(dat_g$external_gene_name)
#   print(gene_name)
#   
#   seq_5utr <- seq_3utr <- c()
#   
#   TSS_manual <- manual_start_end$TSS[manual_start_end$Geneid == gene ]
#   TES_manual <- manual_start_end$TES[manual_start_end$Geneid == gene ]
#   
#   MEME_fasta_weight <- c(MEME_fasta_weight, rep(1/length(unique(dat_g$ensembl_transcript_id)), length(unique(dat_g$ensembl_transcript_id))))
#   
#   # get utr start and end for each transcript
#   for (transcript in unique(dat_g$ensembl_transcript_id)){
#     MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0('>', which(gene_id_ordered == gene), '_', transcript , ' ', rel_transport, ' ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
#     MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0('>', which(gene_id_ordered == gene), '_', transcript , ' ', rel_transport, ' ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
#     
#     dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
#     
#     start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
#     end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
#     start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
#     end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
#     
#     # order if not ordered
#     order_5utr <- order(start_5utr_tmp, decreasing=ifelse(strand==-1, TRUE, FALSE))
#     start_5utr_tmp <- start_5utr_tmp[order_5utr]
#     end_5utr_tmp <- end_5utr_tmp[order_5utr]
#     
#     order_3utr <- order(start_3utr_tmp, decreasing=ifelse(strand==-1, TRUE, FALSE))
#     start_3utr_tmp <- start_3utr_tmp[order_3utr]
#     end_3utr_tmp <- end_3utr_tmp[order_3utr]
#     
#     # remove utr exons outside of TSS-TES range or extend range
#     
#     if(strand == 1){
#       
#       if ( length(start_5utr_tmp) > 0 & !is.na(TSS_manual)){
#         to_clip_or_remove <-  start_5utr_tmp <= TSS_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_5utr_tmp <= TSS_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_5utr_tmp[to_clip] <- TSS_manual
#           start_5utr_tmp <- start_5utr_tmp[!to_remove]
#           end_5utr_tmp <- end_5utr_tmp[!to_remove]
#         }else{ # extend
#           start_5utr_tmp[1] <- TSS_manual
#         }
#       }
#       
#       if ( length(start_3utr_tmp) > 0 &  !is.na(TES_manual)){
#         to_clip_or_remove <-  end_3utr_tmp >= TES_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_3utr_tmp >= TES_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_3utr_tmp[to_clip] <- TES_manual
#           end_3utr_tmp <- end_3utr_tmp[!to_remove]
#           start_3utr_tmp <- start_3utr_tmp[!to_remove]
#         }else{ # extend
#           end_3utr_tmp[length(end_3utr_tmp)] <- TES_manual
#         }
#       }
#       
#     }else if (strand == -1){
#       
#       if (length(start_5utr_tmp) > 0 & !is.na(TSS_manual)){
#         to_clip_or_remove <-  end_5utr_tmp >= TSS_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_5utr_tmp >= TSS_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_5utr_tmp[to_clip] <- TSS_manual
#           end_5utr_tmp <- end_5utr_tmp[!to_remove]
#           start_5utr_tmp <- start_5utr_tmp[!to_remove]
#         }else{ # extend
#           end_5utr_tmp[1] <- TSS_manual
#         }
#       }
#       
#       if (length(start_3utr_tmp) > 0 &  !is.na(TES_manual)){
#         to_clip_or_remove <-  start_3utr_tmp <= TES_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_3utr_tmp <= TES_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_3utr_tmp[to_clip] <- TES_manual
#           start_3utr_tmp <- start_3utr_tmp[!to_remove]
#           end_3utr_tmp <- end_3utr_tmp[!to_remove]
#         }else{ # extend
#           start_3utr_tmp[length(start_3utr_tmp)] <- TES_manual
#         }
#       }
#     }
#     
#     # print(length(start_5utr_tmp))
#     # print(length(start_3utr_tmp))
#     
#     # all_sequences <- c(all_sequence, list())
#     seq_5utr <- c()
#     seq_3utr <- c()
#     
#     # In GET take - strand because RNA is reversed, to get real motifs
#     
#     if(strand == 1){
#       if( length(start_5utr_tmp) > 0 & !is.na(TSS_manual)){
#         for (i in 1:length(start_5utr_tmp)){
#           if(end_5utr_tmp[i] >= TSS_manual){
#             ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',end_5utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_5utr <- c(seq_5utr, content(r)$seq)
#           }
#         }
#       }
#       if( length(start_3utr_tmp) >0 & !is.na(TES_manual)){
#         for (i in 1:length(start_3utr_tmp)){
#           if(start_3utr_tmp[i] <= TES_manual){
#             ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',end_3utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_3utr <- c(seq_3utr, content(r)$seq)
#           }
#         }
#       }
#     }else if (strand == -1){
#       if( length(start_5utr_tmp) >0 & !is.na(TSS_manual)){
#         for (i in 1:length(start_5utr_tmp)){
#           if(start_5utr_tmp[i] <= TSS_manual){
#             ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',end_5utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_5utr <- c(seq_5utr, content(r)$seq)
#           }
#         }
#       }
#       if( length(start_3utr_tmp) > 0 & !is.na(TES_manual) ){
#         for (i in 1:length(start_3utr_tmp)){
#           if(end_3utr_tmp[i] >= TES_manual){
#             ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',end_3utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_3utr <- c(seq_3utr, content(r)$seq)
#           }
#         }
#       }
#     }
#     
#     # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
#     # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
#     MEME_fasta_3utr <- c(MEME_fasta_3utr, dna_to_rna(paste0(seq_3utr, collapse = '')))
#     MEME_fasta_5utr <- c(MEME_fasta_5utr, dna_to_rna(paste0(seq_5utr, collapse = '')))
#   }
# }
# 
# # remove too short sequences length from file and weight
# len_th <- 8
# zero_length_3utr <- which(sapply(MEME_fasta_3utr, function(x){nchar(x) < len_th}))
# MEME_fasta_3utr_to_write <- MEME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
# MEME_fasta_weight_3utr <- MEME_fasta_weight[-c(zero_length_3utr/2)]
# 
# zero_length_5utr <- which(sapply(MEME_fasta_5utr, function(x){nchar(x) < len_th }))
# MEME_fasta_5utr_to_write <- MEME_fasta_3utr[-c(zero_length_5utr-1,zero_length_5utr)]
# MEME_fasta_weight_5utr <- MEME_fasta_weight[-c(zero_length_5utr/2)]
# 
# ## Add weigh to MEME fasta
# MEME_fasta_3utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_3utr), collapse = ' '), MEME_fasta_3utr_to_write)
# MEME_fasta_5utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_5utr), collapse = ' '), MEME_fasta_5utr_to_write)
# 
# ## write to file
# write.table(x = MEME_fasta_3utr_to_write, file = 'meme_v3/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# write.table(x = MEME_fasta_5utr_to_write, file = 'meme_v3/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# 
# for (i in 1:(length(gene_id_ordered)-1)){
#   first_gene_less <- gene_id_ordered[i+1]
#   
#   first_gene_less_idx <- which.min(sapply(MEME_fasta_3utr, function(x){grep(first_gene_less, x)}))
#   
#   MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr[1:(first_gene_less_idx-1)]
#   MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_weight[1:((first_gene_less_idx-1)/2)]
#   
#   zero_length_3utr_more <- which(sapply(MEME_fasta_3utr_to_write_more, function(x){nchar(x)<len_th}))
#   MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
#   MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_3utr_to_write_more_weight[-c(zero_length_3utr_more/2)]
#   
#   MEME_fasta_3utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_more_weight), collapse = ' '), MEME_fasta_3utr_to_write_more)
#   
#   
#   MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr[(first_gene_less_idx):length(MEME_fasta_3utr)]
#   MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
#   
#   zero_length_3utr_less <- which(sapply(MEME_fasta_3utr_to_write_less, function(x){nchar(x)< len_th}))
#   MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
#   MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_3utr_to_write_less_weight[-c(zero_length_3utr_less/2)]
#   
#   MEME_fasta_3utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_less_weight), collapse = ' '), MEME_fasta_3utr_to_write_less)
#   
#   MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr[1:(first_gene_less_idx-1)]
#   MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_weight[1:(first_gene_less_idx-1)/2]
#   
#   zero_length_5utr_more <- which(sapply(MEME_fasta_5utr_to_write_more, function(x){nchar(x) < len_th}))
#   MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
#   MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_5utr_to_write_more_weight[-c(zero_length_5utr_more/2)]
#   
#   MEME_fasta_5utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_more_weight), collapse = ' '), MEME_fasta_5utr_to_write_more)
#   
#   
#   MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr[(first_gene_less_idx):length(MEME_fasta_5utr)]
#   MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
#   
#   zero_length_5utr_less <- which(sapply(MEME_fasta_5utr_to_write_less, function(x){nchar(x) < len_th}))
#   MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
#   MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_5utr_to_write_less_weight[-c(zero_length_5utr_less/2)]
#   
#   MEME_fasta_5utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_less_weight), collapse = ' '), MEME_fasta_5utr_to_write_less)
#   
#   write.table(x = MEME_fasta_3utr_to_write_more, file = paste0('meme_v3/3utr/meme_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_3utr_to_write_less, file = paste0('meme_v3/3utr/meme_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_5utr_to_write_more, file = paste0('meme_v3/5utr/meme_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_5utr_to_write_less, file = paste0('meme_v3/5utr/meme_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   
# }
# 
# 
# ## AME_v2 common part only
# # order gene by transport value for all (more -> less)
# gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]
# 
# AME_fasta_5utr <- AME_fasta_3utr <- c()
# 
# for (gene in gene_id_ordered){
#   dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
#   strand <- unique(dat_g$strand)
#   chr <- unique(dat_g$chromosome_name)
#   
#   rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
#   
#   gene_name <- unique(dat_g$external_gene_name)
#   print(gene_name)
#   
#   AME_fasta_3utr <- c(AME_fasta_3utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
#   AME_fasta_5utr <- c(AME_fasta_5utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
#   
#   TSS_manual <- manual_start_end$TSS[manual_start_end$Geneid == gene ]
#   TES_manual <- manual_start_end$TES[manual_start_end$Geneid == gene ]
#   
#   start_5utr <- start_3utr <- c()
#   end_5utr <- end_3utr <- c()
#   
#   # seq_5utr <- seq_3utr <- c()
#   
#   # get utr start and end for each transcript
#   for (transcript in unique(dat_g$ensembl_transcript_id)){
#     dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
#     
#     start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
#     end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
#     start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
#     end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
#     
#     # order if not, weird case were it wasn't ordered
#     order_5utr <- order(start_5utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
#     start_5utr <- c(start_5utr, list(start_5utr_tmp[order_5utr]))
#     end_5utr <- c(end_5utr, list(end_5utr_tmp[order_5utr]))
#     order_3utr <- order(start_3utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
#     start_3utr <- c(start_3utr, list(start_3utr_tmp[order_3utr]))
#     end_3utr <- c(end_3utr, list(end_3utr_tmp[order_3utr]))
#   }
#   
#   # remove utr exons outside of TSS-TES range or extend range
#   if(strand == 1){
#     for (l in 1:length(start_5utr)){
#       if ( length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
#         to_clip_or_remove <-  start_5utr[[l]] <= TSS_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_5utr[[l]] <= TSS_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_5utr[[l]][to_clip] <- TSS_manual
#           start_5utr[[l]] <- start_5utr[[l]][!to_remove]
#           end_5utr[[l]] <- end_5utr[[l]][!to_remove]
#         }else{ # extend
#           start_5utr[[l]][1] <- TSS_manual
#         }
#       }
#     }
#     for (l in 1:length(start_3utr)){
#       if ( length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
#         to_clip_or_remove <-  end_3utr[[l]] >= TES_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_3utr[[l]] >= TES_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_3utr[[l]][to_clip] <- TES_manual
#           end_3utr[[l]] <- end_3utr[[l]][!to_remove]
#           start_3utr[[l]] <- start_3utr[[l]][!to_remove]      
#         }else{ # extend
#           end_3utr[[l]][length(end_3utr[[l]])] <- TES_manual
#         }
#       }
#     }
#   }else if (strand == -1){
#     for (l in 1:length(start_5utr)){
#       if (length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
#         to_clip_or_remove <-  end_5utr[[l]] >= TSS_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_5utr[[l]] >= TSS_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_5utr[[l]][to_clip] <- TSS_manual
#           end_5utr[[l]] <- end_5utr[[l]][!to_remove]
#           start_5utr[[l]] <- start_5utr[[l]][!to_remove]        
#         }else{ # extend
#           end_5utr[[l]][1] <- TSS_manual
#         }
#       }
#     }
#     
#     for (l in 1:length(start_3utr)){
#       if (length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
#         to_clip_or_remove <-  start_3utr[[l]] <= TES_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_3utr[[l]] <= TES_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_3utr[[l]][to_clip] <- TES_manual
#           start_3utr[[l]] <- start_3utr[[l]][!to_remove]
#           end_3utr[[l]] <- end_3utr[[l]][!to_remove]
#           merge_delim_3utr <- merge_delim_5utr[!to_remove]
#         }else{ # extend
#           start_3utr[[l]][length(start_3utr[[l]])] <- TES_manual
#         }
#       }
#     }
#   }
#   
# 
#   # get common intervals and get sequences and merge with 10 N (max motif width for RBP 8) when inconsistent between transcripts
#   if (length(unique(dat_g$ensembl_transcript_id)) > 1){
#     if(strand == 1){
#       dec <- F
#     }else if (strand == -1){
#       dec <- T
#     }
#     
#     if(any(unlist(lapply(start_3utr, function(l){length(l)})) == 0)){
#       start_3utr_merged <- c()
#       end_3utr_merged <- c()
#       merge_delim_3utr <- c()
#     }else{
#       intervals_3utr <- do.call("rbind", lapply(1:length(start_3utr), function(i){data.frame(start=start_3utr[[i]], end=end_3utr[[i]], "transcript"=i)}))
#       intervals_3utr <- intervals_3utr[order(intervals_3utr$start, decreasing = dec),]
#       
#       tmp <- merge_intervals(intervals_3utr, strand)
#       
#       start_3utr_merged <- tmp[[1]]$start
#       end_3utr_merged <- tmp[[1]]$end
#       merge_delim_3utr <- tmp[[2]]
#     }
#     
#     if(any(unlist(lapply(start_5utr, function(l){length(l)})) == 0)){
#       start_5utr_merged <- c()
#       end_5utr_merged <- c()
#       merge_delim_5utr <- c()
#     }else{
#       intervals_5utr <- do.call("rbind", lapply(1:length(start_5utr), function(i){data.frame(start=start_5utr[[i]], end=end_5utr[[i]], "transcript"=i)}))
#       intervals_5utr <- intervals_5utr[order(intervals_5utr$start, decreasing = dec),]
#       
#       tmp <- merge_intervals(intervals_5utr, strand)
#       
#       start_5utr_merged <- tmp[[1]]$start
#       end_5utr_merged <- tmp[[1]]$end
#       merge_delim_5utr <- tmp[[2]]
#     }
#     
#   }else{
#     start_5utr_merged <- unlist(start_5utr)
#     end_5utr_merged <- unlist(end_5utr)
#     start_3utr_merged <- unlist(start_3utr)
#     end_3utr_merged <- unlist(end_3utr)
#     
#     merge_delim_5utr <- rep("",length(start_5utr_merged))
#     
#     merge_delim_3utr <- rep("",length(start_3utr_merged))
#     
#   }
#   
#   seq_5utr <- c()
#   seq_3utr <- c()
#   # get sequences and merge with 10 N (max motif width for RBP 8)
#   # In GET take - strand because RNA is reversed, to get real motifs
#   if(strand == 1){
#     if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
#       for (i in 1:length(start_5utr_merged)){
#         if(start_5utr_merged[i] >= TSS_manual){
#           ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_5utr <- c(seq_5utr, content(r)$seq)
#         }
#       }
#     }
#     if( length(start_3utr_merged) >0 & !is.na(TES_manual)){
#       for (i in 1:length(start_3utr_merged)){
#         if(start_3utr_merged[i] <= TES_manual){
#           ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..', end_3utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_3utr <- c(seq_3utr, content(r)$seq)
#         }
#       }
#     }
#   }else if (strand == -1){
#     if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
#       for (i in 1:length(start_5utr_merged)){
#         if(start_5utr_merged[i] <= TSS_manual){
#           ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_5utr <- c(seq_5utr, content(r)$seq)
#         }
#       }
#     }
#     if( length(start_3utr_merged) > 0 & !is.na(TES_manual) ){
#       for (i in 1:length(start_3utr_merged)){
#         if(end_3utr_merged[i] >= TES_manual){
#           ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..',end_3utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_3utr <- c(seq_3utr, content(r)$seq)
#         }
#       }
#     }
#   }
#   
#   # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
#   # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
#   AME_fasta_3utr <- c(AME_fasta_3utr, dna_to_rna(paste0(seq_3utr, merge_delim_3utr, sep="", collapse = '')))
#   AME_fasta_5utr <- c(AME_fasta_5utr, dna_to_rna(paste0(seq_5utr, merge_delim_5utr, sep="", collapse = '')))
# 
#   warnings()
#   # tmp <- readline()
#   # if(tmp == "c"){next}else{eval(tmp)}
# }
# 
# # remove zero length from file and weight
# zero_length_3utr <- which(sapply(AME_fasta_3utr, function(x){nchar(x)==0}))
# AME_fasta_3utr_to_write <- AME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
# 
# zero_length_5utr <- which(sapply(AME_fasta_5utr, function(x){nchar(x)==0}))
# AME_fasta_5utr_to_write <- AME_fasta_5utr[-c(zero_length_5utr-1,zero_length_5utr)]
# 
# ## write to file reverse (increase order => low value = positive)
# AME_fasta_3utr_reverse_to_write <- AME_fasta_3utr_to_write[rep(2*((length(AME_fasta_3utr_to_write)/2):1),each=2) + c(-1,0)]
# AME_fasta_5utr_reverse_to_write <- AME_fasta_5utr_to_write[rep(2*((length(AME_fasta_5utr_to_write)/2):1),each=2) + c(-1,0)]
# 
# ## write to file
# write.table(x = AME_fasta_3utr_reverse_to_write, file = 'ame_v2/3utr/meme_3utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
# write.table(x = AME_fasta_5utr_reverse_to_write, file = 'ame_v2/5utr/meme_5utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
# 
# # Invert score to do (high value=positive)
# for (i in 1:(length(AME_fasta_3utr_to_write)/2)){
#   tmp <- strsplit(AME_fasta_3utr_to_write[2*i-1],split = " ")
#   tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
#   AME_fasta_3utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
# }
# for (i in 1:(length(AME_fasta_5utr_to_write)/2)){
#   tmp <- strsplit(AME_fasta_5utr_to_write[2*i-1],split = " ")
#   tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
#   AME_fasta_5utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
# }
# ## write to file
# write.table(x = AME_fasta_3utr_to_write, file = 'ame_v2/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# write.table(x = AME_fasta_5utr_to_write, file = 'ame_v2/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# 
# 
# # manual partition
# for (i in 1:(length(gene_id_ordered)-1)){
#   first_gene_less <- gene_id_ordered[i+1]
#   
#   first_gene_less_idx <- which.min(sapply(AME_fasta_3utr, function(x){grep(first_gene_less, x)}))
#   
#   AME_fasta_3utr_to_write_more <- AME_fasta_3utr[1:(first_gene_less_idx-1)]
#   
#   zero_length_3utr_more <- which(sapply(AME_fasta_3utr_to_write_more, function(x){nchar(x)==0}))
#   AME_fasta_3utr_to_write_more <- AME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
#   
#   AME_fasta_3utr_to_write_less <- AME_fasta_3utr[(first_gene_less_idx):length(AME_fasta_3utr)]
#   
#   zero_length_3utr_less <- which(sapply(AME_fasta_3utr_to_write_less, function(x){nchar(x)==0}))
#   AME_fasta_3utr_to_write_less <- AME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
#   
#   AME_fasta_5utr_to_write_more <- AME_fasta_5utr[1:(first_gene_less_idx-1)]
#   
#   zero_length_5utr_more <- which(sapply(AME_fasta_5utr_to_write_more, function(x){nchar(x)==0}))
#   AME_fasta_5utr_to_write_more <- AME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
#   
#   AME_fasta_5utr_to_write_less <- AME_fasta_5utr[(first_gene_less_idx):length(AME_fasta_5utr)]
#   
#   zero_length_5utr_less <- which(sapply(AME_fasta_5utr_to_write_less, function(x){nchar(x)==0}))
#   AME_fasta_5utr_to_write_less <- AME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
#   
#   write.table(x = AME_fasta_3utr_to_write_more, file = paste0('ame_v2/3utr/ame_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_3utr_to_write_less, file = paste0('ame_v2/3utr/ame_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_5utr_to_write_more, file = paste0('ame_v2/5utr/ame_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_5utr_to_write_less, file = paste0('ame_v2/5utr/ame_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   
# }
# 
# ### V4 use observe TES for each isoform ###-------------
# setwd('/home/diane/Documents/UCLA/Ribo_New/RBP')
# manual_start_end <- read.csv('Curated list_naive.csv', header = T, sep = '\t')
# transcripts <- unname(unlist(sapply(as.character(manual_start_end$Isoforms), function(x){strsplit(x, ', ')[[1]]})))
# 
# library(biomaRt)
# ensembl89 <- useMart(biomart = 'ensembl', host="may2017.archive.ensembl.org")
# mm10_89 <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl89)
# 
# gene_UTR_start_end <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 
#                                            'chromosome_name', 'strand', 'transcript_start', 'transcript_end', 
#                                            '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end', 
#                                            'genomic_coding_start', 'genomic_coding_end'),
#                             filters = 'ensembl_transcript_id', 
#                             values = substr(transcripts,1,18), 
#                             mart = mm10_89
# )
# head(gene_UTR_start_end) # multiple rows for 5'utr or 3'utr because it also has introns
# 
# server <- "https://may2017.rest.ensembl.org"
# library(httr)
# 
# # get relative transport for all genes
# library(compiler)
# source('../Modeling/All_genes_last_5kbs/Model1/Hoffman2/functions_weighting_for_smoothing_include_negbinom_error_in_model_v13.R')
# 
# param_best_all <- c()
# rel_transport_all <- c()
# for ( gene in manual_start_end$Geneid ){
#   gene_name <- manual_start_end$Name[manual_start_end$Geneid == gene]
#   file <- paste0('../Modeling/All_genes_last_5kbs/Model1/Hoffman2/Optimisation/Modeling_v13/Optimisation/BFGS/naive/param_bfgs_1000ri_merged_', gene_name, '_Finished.Rdata')
#   if(file.exists(file)){
#     load(file)
#     
#     # get best parameter
#     par <- bfgs_merged[['all']]$parameter[which.min(bfgs_merged[['all']]$value),]
#     
#     # convert parameter to original
#     par_old <- parameters_convert_cp(par[-c(1:4)], log_in = T, log_out = T, direction = -1)
#     
#     # calculate relative transport
#     rel_transport_all <- rbind(rel_transport_all, data.frame(GeneId=gene, Name=gene_name, Transport = par_old[1]+par_old[3]-par_old[2]))
#     param_best_all <- rbind(param_best_all, data.frame(GeneID=gene, Name=gene_name, "k1' - ca -> " = par_old[1], "k2 - -> np" = par_old[2], "k2' - np -> " = par_old[3], "kdeg" = par_old[4], Transport  = par_old[1]+par_old[3]-par_old[2]))
#   }
# }
# 
# # order gene by transport value for all (more -> less)
# gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]
# 
# # for Meme get sequence of all transcripts and add weight in fasta
# MEME_fasta_5utr <- MEME_fasta_3utr <- c()
# MEME_fasta_weight <- c()
# 
# for (gene in gene_id_ordered){
#   dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
#   strand <- unique(dat_g$strand)
#   chr <- unique(dat_g$chromosome_name)
#   
#   rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
#   
#   gene_name <- unique(dat_g$external_gene_name)
#   print(gene_name)
#   
#   seq_5utr <- seq_3utr <- c()
#   
#   TSS_manual <- as.numeric(strsplit(manual_start_end$TSSs[manual_start_end$Geneid == gene ], ", ")[[1]])
#   TES_manual <- as.numeric(strsplit(manual_start_end$TESs[manual_start_end$Geneid == gene ], ", ")[[1]])
#   transcript_manual <- sapply(strsplit(manual_start_end$Isoforms[manual_start_end$Geneid == gene ], ", ")[[1]], function(x){substr(x,1,18)})
#   
#   MEME_fasta_weight <- c(MEME_fasta_weight, rep(1/length(unique(dat_g$ensembl_transcript_id)), length(unique(dat_g$ensembl_transcript_id))))
#   
#   # get utr start and end for each transcript
#   for (transcript in unique(dat_g$ensembl_transcript_id)){
#     MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0('>', which(gene_id_ordered == gene), '_', transcript , ' ', rel_transport, ' ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
#     MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0('>', which(gene_id_ordered == gene), '_', transcript , ' ', rel_transport, ' ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
#     
#     dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
#     
#     start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
#     end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
#     start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
#     end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
#     
#     # order if not ordered
#     order_5utr <- order(start_5utr_tmp, decreasing=ifelse(strand==-1, TRUE, FALSE))
#     start_5utr_tmp <- start_5utr_tmp[order_5utr]
#     end_5utr_tmp <- end_5utr_tmp[order_5utr]
#     
#     order_3utr <- order(start_3utr_tmp, decreasing=ifelse(strand==-1, TRUE, FALSE))
#     start_3utr_tmp <- start_3utr_tmp[order_3utr]
#     end_3utr_tmp <- end_3utr_tmp[order_3utr]
#     
#     # remove utr exons outside of TSS-TES range or extend range
#     
#     if(strand == 1){
#       
#       if ( length(start_5utr_tmp) > 0 & !is.na(TSS_manual[transcript_manual == transcript])){
#         to_clip_or_remove <-  start_5utr_tmp <= TSS_manual[transcript_manual == transcript]
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_5utr_tmp <= TSS_manual[transcript_manual == transcript]
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_5utr_tmp[to_clip] <- TSS_manual[transcript_manual == transcript]
#           start_5utr_tmp <- start_5utr_tmp[!to_remove]
#           end_5utr_tmp <- end_5utr_tmp[!to_remove]
#         }else{ # extend
#           start_5utr_tmp[1] <- TSS_manual[transcript_manual == transcript]
#         }
#       }
#       
#       if ( length(start_3utr_tmp) > 0 &  !is.na(TES_manual[transcript_manual == transcript])){
#         to_clip_or_remove <-  end_3utr_tmp >= TES_manual[transcript_manual == transcript]
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_3utr_tmp >= TES_manual[transcript_manual == transcript]
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_3utr_tmp[to_clip] <- TES_manual[transcript_manual == transcript]
#           end_3utr_tmp <- end_3utr_tmp[!to_remove]
#           start_3utr_tmp <- start_3utr_tmp[!to_remove]
#         }else{ # extend
#           end_3utr_tmp[length(end_3utr_tmp)] <- TES_manual[transcript_manual == transcript]
#         }
#       }
#       
#     }else if (strand == -1){
#       
#       if (length(start_5utr_tmp) > 0 & !is.na(TSS_manual[transcript_manual == transcript])){
#         to_clip_or_remove <-  end_5utr_tmp >= TSS_manual[transcript_manual == transcript]
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_5utr_tmp >= TSS_manual[transcript_manual == transcript]
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_5utr_tmp[to_clip] <- TSS_manual[transcript_manual == transcript]
#           end_5utr_tmp <- end_5utr_tmp[!to_remove]
#           start_5utr_tmp <- start_5utr_tmp[!to_remove]
#         }else{ # extend
#           end_5utr_tmp[1] <- TSS_manual[transcript_manual == transcript]
#         }
#       }
#       
#       if (length(start_3utr_tmp) > 0 &  !is.na(TES_manual[transcript_manual == transcript])){
#         to_clip_or_remove <-  start_3utr_tmp <= TES_manual[transcript_manual == transcript]
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_3utr_tmp <= TES_manual[transcript_manual == transcript]
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_3utr_tmp[to_clip] <- TES_manual[transcript_manual == transcript]
#           start_3utr_tmp <- start_3utr_tmp[!to_remove]
#           end_3utr_tmp <- end_3utr_tmp[!to_remove]
#         }else{ # extend
#           start_3utr_tmp[length(start_3utr_tmp)] <- TES_manual[transcript_manual == transcript]
#         }
#       }
#     }
#     
#     # print(length(start_5utr_tmp))
#     # print(length(start_3utr_tmp))
#     
#     # all_sequences <- c(all_sequence, list())
#     seq_5utr <- c()
#     seq_3utr <- c()
#     
#     # In GET take - strand because RNA is reversed, to get real motifs
#     
#     if(strand == 1){
#       if( length(start_5utr_tmp) > 0 & !is.na(TSS_manual[transcript_manual == transcript])){
#         for (i in 1:length(start_5utr_tmp)){
#           if(end_5utr_tmp[i] >= TSS_manual[transcript_manual == transcript]){
#             ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',end_5utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_5utr <- c(seq_5utr, content(r)$seq)
#           }
#         }
#       }
#       if( length(start_3utr_tmp) >0 & !is.na(TES_manual[transcript_manual == transcript])){
#         for (i in 1:length(start_3utr_tmp)){
#           if(start_3utr_tmp[i] <= TES_manual[transcript_manual == transcript]){
#             ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',end_3utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_3utr <- c(seq_3utr, content(r)$seq)
#           }
#         }
#       }
#     }else if (strand == -1){
#       if( length(start_5utr_tmp) >0 & !is.na(TSS_manual[transcript_manual == transcript])){
#         for (i in 1:length(start_5utr_tmp)){
#           if(start_5utr_tmp[i] <= TSS_manual[transcript_manual == transcript]){
#             ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',end_5utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_5utr <- c(seq_5utr, content(r)$seq)
#           }
#         }
#       }
#       if( length(start_3utr_tmp) > 0 & !is.na(TES_manual[transcript_manual == transcript]) ){
#         for (i in 1:length(start_3utr_tmp)){
#           if(end_3utr_tmp[i] >= TES_manual[transcript_manual == transcript]){
#             ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',end_3utr_tmp[i],':', strand,"?mask=soft")
#             r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#             stop_for_status(r)
#             seq_3utr <- c(seq_3utr, content(r)$seq)
#           }
#         }
#       }
#     }
#     
#     # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
#     # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
#     MEME_fasta_3utr <- c(MEME_fasta_3utr, dna_to_rna(paste0(seq_3utr, collapse = '')))
#     MEME_fasta_5utr <- c(MEME_fasta_5utr, dna_to_rna(paste0(seq_5utr, collapse = '')))
#   }
# }
# 
# # remove too short sequences length from file and weight
# len_th <- 8
# zero_length_3utr <- which(sapply(MEME_fasta_3utr, function(x){nchar(x) < len_th}))
# MEME_fasta_3utr_to_write <- MEME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
# MEME_fasta_weight_3utr <- MEME_fasta_weight[-c(zero_length_3utr/2)]
# 
# zero_length_5utr <- which(sapply(MEME_fasta_5utr, function(x){nchar(x) < len_th }))
# MEME_fasta_5utr_to_write <- MEME_fasta_3utr[-c(zero_length_5utr-1,zero_length_5utr)]
# MEME_fasta_weight_5utr <- MEME_fasta_weight[-c(zero_length_5utr/2)]
# 
# ## Add weigh to MEME fasta
# MEME_fasta_3utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_3utr), collapse = ' '), MEME_fasta_3utr_to_write)
# MEME_fasta_5utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_5utr), collapse = ' '), MEME_fasta_5utr_to_write)
# 
# ## write to file
# write.table(x = MEME_fasta_3utr_to_write, file = 'meme_v4/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# write.table(x = MEME_fasta_5utr_to_write, file = 'meme_v4/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# 
# for (i in 1:(length(gene_id_ordered)-1)){
#   first_gene_less <- gene_id_ordered[i+1]
#   
#   first_gene_less_idx <- which.min(sapply(MEME_fasta_3utr, function(x){grep(first_gene_less, x)}))
#   
#   MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr[1:(first_gene_less_idx-1)]
#   MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_weight[1:((first_gene_less_idx-1)/2)]
#   
#   zero_length_3utr_more <- which(sapply(MEME_fasta_3utr_to_write_more, function(x){nchar(x)<len_th}))
#   MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
#   MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_3utr_to_write_more_weight[-c(zero_length_3utr_more/2)]
#   
#   MEME_fasta_3utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_more_weight), collapse = ' '), MEME_fasta_3utr_to_write_more)
#   
#   
#   MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr[(first_gene_less_idx):length(MEME_fasta_3utr)]
#   MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
#   
#   zero_length_3utr_less <- which(sapply(MEME_fasta_3utr_to_write_less, function(x){nchar(x)< len_th}))
#   MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
#   MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_3utr_to_write_less_weight[-c(zero_length_3utr_less/2)]
#   
#   MEME_fasta_3utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_less_weight), collapse = ' '), MEME_fasta_3utr_to_write_less)
#   
#   MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr[1:(first_gene_less_idx-1)]
#   MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_weight[1:(first_gene_less_idx-1)/2]
#   
#   zero_length_5utr_more <- which(sapply(MEME_fasta_5utr_to_write_more, function(x){nchar(x) < len_th}))
#   MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
#   MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_5utr_to_write_more_weight[-c(zero_length_5utr_more/2)]
#   
#   MEME_fasta_5utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_more_weight), collapse = ' '), MEME_fasta_5utr_to_write_more)
#   
#   
#   MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr[(first_gene_less_idx):length(MEME_fasta_5utr)]
#   MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
#   
#   zero_length_5utr_less <- which(sapply(MEME_fasta_5utr_to_write_less, function(x){nchar(x) < len_th}))
#   MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
#   MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_5utr_to_write_less_weight[-c(zero_length_5utr_less/2)]
#   
#   MEME_fasta_5utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_less_weight), collapse = ' '), MEME_fasta_5utr_to_write_less)
#   
#   write.table(x = MEME_fasta_3utr_to_write_more, file = paste0('meme_v4/3utr/meme_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_3utr_to_write_less, file = paste0('meme_v4/3utr/meme_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_5utr_to_write_more, file = paste0('meme_v4/5utr/meme_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = MEME_fasta_5utr_to_write_less, file = paste0('meme_v4/5utr/meme_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   
# }
# 
# 
# ## AME_v4 common part only
# # order gene by transport value for all (more -> less)
# gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]
# 
# AME_fasta_5utr <- AME_fasta_3utr <- c()
# 
# for (gene in gene_id_ordered){
#   dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
#   strand <- unique(dat_g$strand)
#   chr <- unique(dat_g$chromosome_name)
#   
#   rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
#   
#   gene_name <- unique(dat_g$external_gene_name)
#   print(gene_name)
#   
#   AME_fasta_3utr <- c(AME_fasta_3utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
#   AME_fasta_5utr <- c(AME_fasta_5utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
#   
#   
#   TSS_manual <- as.numeric(strsplit(manual_start_end$TSSs[manual_start_end$Geneid == gene ], ", ")[[1]])
#   TES_manual <- as.numeric(strsplit(manual_start_end$TESs[manual_start_end$Geneid == gene ], ", ")[[1]])
# 
#   if(strand == 1){
#     TSS_manual <- max(TSS_manual)
#     TES_manual <- min(TES_manual)
#   }else if (strand == -1){
#     TSS_manual <- min(TSS_manual)
#     TES_manual <- max(TES_manual)
#   }
#   
#   start_5utr <- start_3utr <- c()
#   end_5utr <- end_3utr <- c()
#   
#   # seq_5utr <- seq_3utr <- c()
#   
#   # get utr start and end for each transcript
#   for (transcript in unique(dat_g$ensembl_transcript_id)){
#     dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
#     
#     start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
#     end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
#     start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
#     end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
#     
#     # order if not, weird case were it wasn't ordered
#     order_5utr <- order(start_5utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
#     start_5utr <- c(start_5utr, list(start_5utr_tmp[order_5utr]))
#     end_5utr <- c(end_5utr, list(end_5utr_tmp[order_5utr]))
#     order_3utr <- order(start_3utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
#     start_3utr <- c(start_3utr, list(start_3utr_tmp[order_3utr]))
#     end_3utr <- c(end_3utr, list(end_3utr_tmp[order_3utr]))
#   }
#   
#   # remove utr exons outside of TSS-TES range or extend range
#   if(strand == 1){
#     for (l in 1:length(start_5utr)){
#       if ( length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
#         to_clip_or_remove <-  start_5utr[[l]] <= TSS_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_5utr[[l]] <= TSS_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_5utr[[l]][to_clip] <- TSS_manual
#           start_5utr[[l]] <- start_5utr[[l]][!to_remove]
#           end_5utr[[l]] <- end_5utr[[l]][!to_remove]
#         }else{ # extend
#           start_5utr[[l]][1] <- TSS_manual
#         }
#       }
#     }
#     for (l in 1:length(start_3utr)){
#       if ( length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
#         to_clip_or_remove <-  end_3utr[[l]] >= TES_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_3utr[[l]] >= TES_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_3utr[[l]][to_clip] <- TES_manual
#           end_3utr[[l]] <- end_3utr[[l]][!to_remove]
#           start_3utr[[l]] <- start_3utr[[l]][!to_remove]      
#         }else{ # extend
#           end_3utr[[l]][length(end_3utr[[l]])] <- TES_manual
#         }
#       }
#     }
#   }else if (strand == -1){
#     for (l in 1:length(start_5utr)){
#       if (length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
#         to_clip_or_remove <-  end_5utr[[l]] >= TSS_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & start_5utr[[l]] >= TSS_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           end_5utr[[l]][to_clip] <- TSS_manual
#           end_5utr[[l]] <- end_5utr[[l]][!to_remove]
#           start_5utr[[l]] <- start_5utr[[l]][!to_remove]        
#         }else{ # extend
#           end_5utr[[l]][1] <- TSS_manual
#         }
#       }
#     }
#     
#     for (l in 1:length(start_3utr)){
#       if (length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
#         to_clip_or_remove <-  start_3utr[[l]] <= TES_manual
#         if( any(to_clip_or_remove)) {
#           to_remove <- to_clip_or_remove & end_3utr[[l]] <= TES_manual
#           to_clip <- to_clip_or_remove & !to_remove
#           
#           start_3utr[[l]][to_clip] <- TES_manual
#           start_3utr[[l]] <- start_3utr[[l]][!to_remove]
#           end_3utr[[l]] <- end_3utr[[l]][!to_remove]
#           merge_delim_3utr <- merge_delim_5utr[!to_remove]
#         }else{ # extend
#           start_3utr[[l]][length(start_3utr[[l]])] <- TES_manual
#         }
#       }
#     }
#   }
#   
#   
#   # get common intervals and get sequences and merge with 10 N (max motif width for RBP 8) when inconsistent between transcripts
#   if (length(unique(dat_g$ensembl_transcript_id)) > 1){
#     if(strand == 1){
#       dec <- F
#     }else if (strand == -1){
#       dec <- T
#     }
#     
#     if(any(unlist(lapply(start_3utr, function(l){length(l)})) == 0)){
#       start_3utr_merged <- c()
#       end_3utr_merged <- c()
#       merge_delim_3utr <- c()
#     }else{
#       intervals_3utr <- do.call("rbind", lapply(1:length(start_3utr), function(i){data.frame(start=start_3utr[[i]], end=end_3utr[[i]], "transcript"=i)}))
#       intervals_3utr <- intervals_3utr[order(intervals_3utr$start, decreasing = dec),]
#       
#       tmp <- merge_intervals(intervals_3utr, strand)
#       
#       start_3utr_merged <- tmp[[1]]$start
#       end_3utr_merged <- tmp[[1]]$end
#       merge_delim_3utr <- tmp[[2]]
#     }
#     
#     if(any(unlist(lapply(start_5utr, function(l){length(l)})) == 0)){
#       start_5utr_merged <- c()
#       end_5utr_merged <- c()
#       merge_delim_5utr <- c()
#     }else{
#       intervals_5utr <- do.call("rbind", lapply(1:length(start_5utr), function(i){data.frame(start=start_5utr[[i]], end=end_5utr[[i]], "transcript"=i)}))
#       intervals_5utr <- intervals_5utr[order(intervals_5utr$start, decreasing = dec),]
#       
#       tmp <- merge_intervals(intervals_5utr, strand)
#       
#       start_5utr_merged <- tmp[[1]]$start
#       end_5utr_merged <- tmp[[1]]$end
#       merge_delim_5utr <- tmp[[2]]
#     }
#     
#   }else{
#     start_5utr_merged <- unlist(start_5utr)
#     end_5utr_merged <- unlist(end_5utr)
#     start_3utr_merged <- unlist(start_3utr)
#     end_3utr_merged <- unlist(end_3utr)
#     
#     merge_delim_5utr <- rep("",length(start_5utr_merged))
#     
#     merge_delim_3utr <- rep("",length(start_3utr_merged))
#     
#   }
#   
#   seq_5utr <- c()
#   seq_3utr <- c()
#   # get sequences and merge with 10 N (max motif width for RBP 8)
#   # In GET take - strand because RNA is reversed, to get real motifs
#   if(strand == 1){
#     if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
#       for (i in 1:length(start_5utr_merged)){
#         if(start_5utr_merged[i] >= TSS_manual){
#           ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_5utr <- c(seq_5utr, content(r)$seq)
#         }
#       }
#     }
#     if( length(start_3utr_merged) >0 & !is.na(TES_manual)){
#       for (i in 1:length(start_3utr_merged)){
#         if(start_3utr_merged[i] <= TES_manual){
#           ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..', end_3utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_3utr <- c(seq_3utr, content(r)$seq)
#         }
#       }
#     }
#   }else if (strand == -1){
#     if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
#       for (i in 1:length(start_5utr_merged)){
#         if(start_5utr_merged[i] <= TSS_manual){
#           ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_5utr <- c(seq_5utr, content(r)$seq)
#         }
#       }
#     }
#     if( length(start_3utr_merged) > 0 & !is.na(TES_manual) ){
#       for (i in 1:length(start_3utr_merged)){
#         if(end_3utr_merged[i] >= TES_manual){
#           ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..',end_3utr_merged[i],':', strand,"?mask=soft")
#           r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
#           stop_for_status(r)
#           seq_3utr <- c(seq_3utr, content(r)$seq)
#         }
#       }
#     }
#   }
#   
#   # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
#   # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
#   AME_fasta_3utr <- c(AME_fasta_3utr, dna_to_rna(paste0(seq_3utr, merge_delim_3utr, sep="", collapse = '')))
#   AME_fasta_5utr <- c(AME_fasta_5utr, dna_to_rna(paste0(seq_5utr, merge_delim_5utr, sep="", collapse = '')))
#   
#   warnings()
#   # tmp <- readline()
#   # if(tmp == "c"){next}else{eval(tmp)}
# }
# 
# # remove zero length from file and weight
# zero_length_3utr <- which(sapply(AME_fasta_3utr, function(x){nchar(x)==0}))
# AME_fasta_3utr_to_write <- AME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
# 
# zero_length_5utr <- which(sapply(AME_fasta_5utr, function(x){nchar(x)==0}))
# AME_fasta_5utr_to_write <- AME_fasta_5utr[-c(zero_length_5utr-1,zero_length_5utr)]
# 
# ## write to file reverse (increase order => low value = positive)
# AME_fasta_3utr_reverse_to_write <- AME_fasta_3utr_to_write[rep(2*((length(AME_fasta_3utr_to_write)/2):1),each=2) + c(-1,0)]
# AME_fasta_5utr_reverse_to_write <- AME_fasta_5utr_to_write[rep(2*((length(AME_fasta_5utr_to_write)/2):1),each=2) + c(-1,0)]
# 
# ## write to file
# write.table(x = AME_fasta_3utr_reverse_to_write, file = 'ame_v4/3utr/meme_3utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
# write.table(x = AME_fasta_5utr_reverse_to_write, file = 'ame_v4/5utr/meme_5utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
# 
# # Invert score to do (high value=positive)
# for (i in 1:(length(AME_fasta_3utr_to_write)/2)){
#   tmp <- strsplit(AME_fasta_3utr_to_write[2*i-1],split = " ")
#   tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
#   AME_fasta_3utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
# }
# for (i in 1:(length(AME_fasta_5utr_to_write)/2)){
#   tmp <- strsplit(AME_fasta_5utr_to_write[2*i-1],split = " ")
#   tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
#   AME_fasta_5utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
# }
# ## write to file
# write.table(x = AME_fasta_3utr_to_write, file = 'ame_v4/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# write.table(x = AME_fasta_5utr_to_write, file = 'ame_v4/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
# 
# 
# # manual partition
# for (i in 1:(length(gene_id_ordered)-1)){
#   first_gene_less <- gene_id_ordered[i+1]
#   
#   first_gene_less_idx <- which.min(sapply(AME_fasta_3utr, function(x){grep(first_gene_less, x)}))
#   
#   AME_fasta_3utr_to_write_more <- AME_fasta_3utr[1:(first_gene_less_idx-1)]
#   
#   zero_length_3utr_more <- which(sapply(AME_fasta_3utr_to_write_more, function(x){nchar(x)==0}))
#   AME_fasta_3utr_to_write_more <- AME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
#   
#   AME_fasta_3utr_to_write_less <- AME_fasta_3utr[(first_gene_less_idx):length(AME_fasta_3utr)]
#   
#   zero_length_3utr_less <- which(sapply(AME_fasta_3utr_to_write_less, function(x){nchar(x)==0}))
#   AME_fasta_3utr_to_write_less <- AME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
#   
#   AME_fasta_5utr_to_write_more <- AME_fasta_5utr[1:(first_gene_less_idx-1)]
#   
#   zero_length_5utr_more <- which(sapply(AME_fasta_5utr_to_write_more, function(x){nchar(x)==0}))
#   AME_fasta_5utr_to_write_more <- AME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
#   
#   AME_fasta_5utr_to_write_less <- AME_fasta_5utr[(first_gene_less_idx):length(AME_fasta_5utr)]
#   
#   zero_length_5utr_less <- which(sapply(AME_fasta_5utr_to_write_less, function(x){nchar(x)==0}))
#   AME_fasta_5utr_to_write_less <- AME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
#   
#   write.table(x = AME_fasta_3utr_to_write_more, file = paste0('ame_v4/3utr/ame_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_3utr_to_write_less, file = paste0('ame_v4/3utr/ame_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_5utr_to_write_more, file = paste0('ame_v4/5utr/ame_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   write.table(x = AME_fasta_5utr_to_write_less, file = paste0('ame_v4/5utr/ame_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
#   
# }
# 

### V5 and Final - use observe TES/TSS for each isoform - corrected ###-------------
setwd('/home/diane/Documents/UCLA/Ribo_New/RBP')
manual_start_end <- read.csv('Curated list_naive.csv', header = T, sep = '\t')
transcripts <- unname(unlist(sapply(as.character(manual_start_end$Isoforms), function(x){strsplit(x, ', ')[[1]]})))

library(biomaRt)
ensembl89 <- useMart(biomart = 'ensembl', host="may2017.archive.ensembl.org")
mm10_89 <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl89)

gene_UTR_start_end <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id', 
                                           'chromosome_name', 'strand', 'transcript_start', 'transcript_end', 
                                           '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end', 
                                           'genomic_coding_start', 'genomic_coding_end'),
                            filters = 'ensembl_transcript_id', 
                            values = substr(transcripts,1,18), 
                            mart = mm10_89
)
head(gene_UTR_start_end) # multiple rows for 5'utr or 3'utr because it also has introns

server <- "https://may2017.rest.ensembl.org"
library(httr)

dna_to_rna <- function(dna){
  rna <- gsub("t", "u", dna)
  rna <- gsub("T", "U", rna)
  return(rna) 
}

merge_intervals <- function(intervals, strand){
  intervals_merged <- intervals[intervals$transcript == 1, c("start","end")]
  merge_delim <- rep("",nrow(intervals_merged))

  for (i in 2:length(unique(intervals$transcript))){
    intervals_tmp <- intervals[intervals$transcript == i, ]

    if(strand == 1){
      # remove exons of 2nd transcript that ends before 1st exons of 1st transcript starts
      intervals_tmp <- intervals_tmp[!(intervals_tmp$end < intervals_merged$start[1]),]

      # remove exons of 1st transcript that ends before 1st exons of 2nd transcript starts
      merge_delim <- merge_delim[!(intervals_merged$end < intervals_tmp$start[1])]
      intervals_merged <- intervals_merged[!(intervals_merged$end < intervals_tmp$start[1]),]

      # remove exons of 2nd transcript that start after last exons of 1st transcript ends
      intervals_tmp <- intervals_tmp[!(intervals_tmp$start > intervals_merged$end[nrow(intervals_merged)]),]

      # remove exons of 1st transcript that start after last exons of 2nd transcript ends
      merge_delim <- merge_delim[!(intervals_merged$start > intervals_tmp$end[nrow(intervals_tmp)])]
      intervals_merged <- intervals_merged[!(intervals_merged$start > intervals_tmp$end[nrow(intervals_tmp)]),]
    }else if (strand == -1){
      # remove exons of 2nd transcript that ends before 1st exons of 1st transcript starts
      intervals_tmp <- intervals_tmp[!(intervals_tmp$end < intervals_merged$start[nrow(intervals_merged)]),]

      # remove exons of 1st transcript that ends before 1st exons of 2nd transcript starts
      merge_delim <- merge_delim[!(intervals_merged$end < intervals_tmp$start[nrow(intervals_tmp)])]
      intervals_merged <- intervals_merged[!(intervals_merged$end < intervals_tmp$start[nrow(intervals_tmp)]),]

      # remove exons of 2nd transcript that start after last exons of 1st transcript ends
      intervals_tmp <- intervals_tmp[!(intervals_tmp$start > intervals_merged$end[1]),]

      # remove exons of 1st transcript that start after last exons of 2nd transcript ends
      merge_delim <- merge_delim[!(intervals_merged$start > intervals_tmp$end[1])]
      intervals_merged <- intervals_merged[!(intervals_merged$start > intervals_tmp$end[1]),]
    }

    # remove exons of 2nd transcript that are not overlapping with any exons of 1st transcript
    overlapping <- c()
    for ( j in 1:nrow(intervals_tmp)){
      overlapping <- c(overlapping, any(intervals_tmp$start[j] >=  intervals_merged$start & intervals_tmp$start[j] <=  intervals_merged$end) | any(intervals_tmp$end[j] >=  intervals_merged$start & intervals_tmp$end[j] <=  intervals_merged$end))
    }
    intervals_tmp <- intervals_tmp[overlapping,]

    # remove exons of 1st transcript that are not overlapping with any exons of 2nd transcript
    overlapping <- c()
    for ( j in 1:nrow(intervals_merged)){
      overlapping <- c(overlapping, any(intervals_merged$start[j] >=  intervals_tmp$start & intervals_merged$start[j] <=  intervals_tmp$end) | any(intervals_merged$end[j] >=  intervals_tmp$start & intervals_merged$end[j] <=  intervals_tmp$end))
    }
    intervals_merged <- intervals_merged[overlapping,]
    merge_delim <- merge_delim[overlapping]

    # merge overlapping exons

    for (j in 1:nrow(intervals_tmp)){
      start_in <- intervals_tmp$start[j] >=  intervals_merged$start & intervals_tmp$start[j] <=  intervals_merged$end
      end_in <- intervals_tmp$end[j] >=  intervals_merged$start & intervals_tmp$end[j] <=  intervals_merged$end

      split <- FALSE

      if(j != nrow(intervals_tmp)){
        start_in_next <- intervals_tmp$start[j+1] >=  intervals_merged$start & intervals_tmp$start[j+1] <=  intervals_merged$end
        if(any(start_in_next == end_in)){ # => split
          split <- TRUE
        }
      }

      # if exon j start of 2nd transcript within other exon then modify the 1st transcript
      if (any(start_in)){
        if(which(start_in) != 1 & intervals_merged$start[which(start_in)] != intervals_tmp$start[j]){
          merge_delim[which(start_in)-1] <- paste0(rep("N",10),collapse="")
        }
        intervals_merged$start[which(start_in)] <- intervals_tmp$start[j]
      }


      if (any(end_in)){
        if(!split){
          if(which(end_in) != nrow(intervals_merged) & intervals_merged$end[which(end_in)] != intervals_tmp$end[j]){
            merge_delim[which(end_in)] <- paste0(rep("N",10),collapse="")
          }
          intervals_merged$end[which(end_in)] <- intervals_tmp$end[j]
        }else{
          intervals_merged <- rbind(intervals_merged[1:which(end_in),], intervals_merged[which(end_in):nrow(intervals_merged),])
          if(which(end_in) !=1 ){
            merge_delim <- c(merge_delim[1:(which(end_in)-1)],"", merge_delim[(which(end_in)-1):length(merge_delim)])
          }else{
            merge_delim <- c("", merge_delim[(which(end_in)-1):length(merge_delim)])
          }
          intervals_merged$end[which(end_in)] <- intervals_tmp$end[j]
          intervals_merged$start[which(end_in)+1] <- intervals_tmp$end[j]+1
        }
      }
    }
  }
  warnings()
  res <- list(intervals_merged,merge_delim)
}

# get relative transport for all genes
library(compiler)
source('../Modeling/All_genes_last_5kbs/Model1/Hoffman2/functions_weighting_for_smoothing_include_negbinom_error_in_model_v13.R')

param_best_all <- c()
rel_transport_all <- c()
for ( gene in manual_start_end$Geneid ){
  gene_name <- manual_start_end$Name[manual_start_end$Geneid == gene]
  file <- paste0('../Modeling/All_genes_last_5kbs/Model1/Hoffman2/Optimisation/Modeling_v13/Optimisation/BFGS/naive/param_bfgs_1000ri_merged_', gene_name, '_Finished.Rdata')
  if(file.exists(file)){
    load(file)

    # get best parameter
    par <- bfgs_merged[['all']]$parameter[which.min(bfgs_merged[['all']]$value),]

    # convert parameter to original
    par_old <- parameters_convert_cp(par[-c(1:4)], log_in = T, log_out = T, direction = -1)

    # calculate relative transport
    rel_transport_all <- rbind(rel_transport_all, data.frame(GeneId=gene, Name=gene_name, Transport = par_old[1]+par_old[3]-par_old[2]))
    param_best_all <- rbind(param_best_all, data.frame(GeneID=gene, Name=gene_name, "k1' - ca -> " = par_old[1], "k2 - -> np" = par_old[2], "k2' - np -> " = par_old[3], "kdeg" = par_old[4], Transport  = par_old[1]+par_old[3]-par_old[2]))
  }
}

# order gene by transport value for all (more -> less)
gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]

# for Meme get sequence of all transcripts and add weight in fasta
MEME_fasta_5utr <- MEME_fasta_3utr <- c()
MEME_fasta_weight <- c()

for (gene in gene_id_ordered){
  dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
  strand <- unique(dat_g$strand)
  chr <- unique(dat_g$chromosome_name)
  
  rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
  
  gene_name <- unique(dat_g$external_gene_name)
  print(gene_name)
  
  seq_5utr <- seq_3utr <- c()
  
  TSS_manual <- as.numeric(strsplit(manual_start_end$TSSs[manual_start_end$Geneid == gene ], ", ")[[1]])
  TES_manual <- as.numeric(strsplit(manual_start_end$TESs[manual_start_end$Geneid == gene ], ", ")[[1]])
  transcript_manual <- sapply(strsplit(manual_start_end$Isoforms[manual_start_end$Geneid == gene ], ", ")[[1]], function(x){substr(x,1,18)})
  
  MEME_fasta_weight <- c(MEME_fasta_weight, rep(1/length(unique(dat_g$ensembl_transcript_id)), length(unique(dat_g$ensembl_transcript_id))))
  
  # get utr start and end for each transcript
  for (transcript in unique(dat_g$ensembl_transcript_id)){
    MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0('>', which(gene_id_ordered == gene), '_', transcript , ' ', rel_transport, ' ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
    MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0('>', which(gene_id_ordered == gene), '_', transcript , ' ', rel_transport, ' ', gene_name, ' | ', gene, ' | ', strand, ' | ', rel_transport))
    
    dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
    
    start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
    end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
    start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
    end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
    
    # order if not ordered
    order_5utr <- order(start_5utr_tmp, decreasing=ifelse(strand==-1, TRUE, FALSE))
    start_5utr_tmp <- start_5utr_tmp[order_5utr]
    end_5utr_tmp <- end_5utr_tmp[order_5utr]
    
    order_3utr <- order(start_3utr_tmp, decreasing=ifelse(strand==-1, TRUE, FALSE))
    start_3utr_tmp <- start_3utr_tmp[order_3utr]
    end_3utr_tmp <- end_3utr_tmp[order_3utr]
    
    # remove utr exons outside of TSS-TES range or extend range
    
    if(strand == 1){
      
      if ( length(start_5utr_tmp) > 0 & !is.na(TSS_manual[transcript_manual == transcript])){
        to_clip_or_remove <-  start_5utr_tmp <= TSS_manual[transcript_manual == transcript]
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & end_5utr_tmp <= TSS_manual[transcript_manual == transcript]
          to_clip <- to_clip_or_remove & !to_remove
          
          start_5utr_tmp[to_clip] <- TSS_manual[transcript_manual == transcript]
          start_5utr_tmp <- start_5utr_tmp[!to_remove]
          end_5utr_tmp <- end_5utr_tmp[!to_remove]
        }else{ # extend
          start_5utr_tmp[1] <- TSS_manual[transcript_manual == transcript]
        }
      }
      
      if ( length(start_3utr_tmp) > 0 &  !is.na(TES_manual[transcript_manual == transcript])){
        to_clip_or_remove <-  end_3utr_tmp >= TES_manual[transcript_manual == transcript]
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & start_3utr_tmp >= TES_manual[transcript_manual == transcript]
          to_clip <- to_clip_or_remove & !to_remove
          
          end_3utr_tmp[to_clip] <- TES_manual[transcript_manual == transcript]
          end_3utr_tmp <- end_3utr_tmp[!to_remove]
          start_3utr_tmp <- start_3utr_tmp[!to_remove]
        }else{ # extend
          end_3utr_tmp[length(end_3utr_tmp)] <- TES_manual[transcript_manual == transcript]
        }
      }
      
    }else if (strand == -1){
      
      if (length(start_5utr_tmp) > 0 & !is.na(TSS_manual[transcript_manual == transcript])){
        to_clip_or_remove <-  end_5utr_tmp >= TSS_manual[transcript_manual == transcript]
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & start_5utr_tmp >= TSS_manual[transcript_manual == transcript]
          to_clip <- to_clip_or_remove & !to_remove
          
          end_5utr_tmp[to_clip] <- TSS_manual[transcript_manual == transcript]
          end_5utr_tmp <- end_5utr_tmp[!to_remove]
          start_5utr_tmp <- start_5utr_tmp[!to_remove]
        }else{ # extend
          end_5utr_tmp[1] <- TSS_manual[transcript_manual == transcript]
        }
      }
      
      if (length(start_3utr_tmp) > 0 &  !is.na(TES_manual[transcript_manual == transcript])){
        to_clip_or_remove <-  start_3utr_tmp <= TES_manual[transcript_manual == transcript]
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & end_3utr_tmp <= TES_manual[transcript_manual == transcript]
          to_clip <- to_clip_or_remove & !to_remove
          
          start_3utr_tmp[to_clip] <- TES_manual[transcript_manual == transcript]
          start_3utr_tmp <- start_3utr_tmp[!to_remove]
          end_3utr_tmp <- end_3utr_tmp[!to_remove]
        }else{ # extend
          start_3utr_tmp[length(start_3utr_tmp)] <- TES_manual[transcript_manual == transcript]
        }
      }
    }
    
    # print(length(start_5utr_tmp))
    # print(length(start_3utr_tmp))

    
    # all_sequences <- c(all_sequence, list())
    seq_5utr <- c()
    seq_3utr <- c()
    
    # In GET take - strand because RNA is reversed, to get real motifs
    
    if(strand == 1){
      if( length(start_5utr_tmp) > 0 & !is.na(TSS_manual[transcript_manual == transcript])){
        for (i in 1:length(start_5utr_tmp)){
          if(end_5utr_tmp[i] >= TSS_manual[transcript_manual == transcript]){
            ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',end_5utr_tmp[i],':', strand,"?mask=soft")
            r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
            stop_for_status(r)
            seq_5utr <- c(seq_5utr, content(r)$seq)
          }
        }
      }
      if( length(start_3utr_tmp) >0 & !is.na(TES_manual[transcript_manual == transcript])){
        for (i in 1:length(start_3utr_tmp)){
          if(start_3utr_tmp[i] <= TES_manual[transcript_manual == transcript]){
            ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',end_3utr_tmp[i],':', strand,"?mask=soft")
            r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
            stop_for_status(r)
            seq_3utr <- c(seq_3utr, content(r)$seq)
          }
        }
      }
    }else if (strand == -1){
      if( length(start_5utr_tmp) >0 & !is.na(TSS_manual[transcript_manual == transcript])){
        for (i in 1:length(start_5utr_tmp)){
          if(start_5utr_tmp[i] <= TSS_manual[transcript_manual == transcript]){
            ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_tmp[i],'..',end_5utr_tmp[i],':', strand,"?mask=soft")
            r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
            stop_for_status(r)
            seq_5utr <- c(seq_5utr, content(r)$seq)
          }
        }
      }
      if( length(start_3utr_tmp) > 0 & !is.na(TES_manual[transcript_manual == transcript]) ){
        for (i in 1:length(start_3utr_tmp)){
          if(end_3utr_tmp[i] >= TES_manual[transcript_manual == transcript]){
            ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_tmp[i],'..',end_3utr_tmp[i],':', strand,"?mask=soft")
            r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
            stop_for_status(r)
            seq_3utr <- c(seq_3utr, content(r)$seq)
          }
        }
      }
    }
    
    # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
    # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
    MEME_fasta_3utr <- c(MEME_fasta_3utr, dna_to_rna(paste0(seq_3utr, collapse = '')))
    MEME_fasta_5utr <- c(MEME_fasta_5utr, dna_to_rna(paste0(seq_5utr, collapse = '')))
  }
}

# remove too short sequences length from file and weight
len_th <- 8
zero_length_3utr <- which(sapply(MEME_fasta_3utr, function(x){nchar(x) < len_th}))
MEME_fasta_3utr_to_write <- MEME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
MEME_fasta_weight_3utr <- MEME_fasta_weight[-c(zero_length_3utr/2)]

zero_length_5utr <- which(sapply(MEME_fasta_5utr, function(x){nchar(x) < len_th }))
MEME_fasta_5utr_to_write <- MEME_fasta_5utr[-c(zero_length_5utr-1,zero_length_5utr)]
MEME_fasta_weight_5utr <- MEME_fasta_weight[-c(zero_length_5utr/2)]

## Add weigh to MEME fasta
MEME_fasta_3utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_3utr), collapse = ' '), MEME_fasta_3utr_to_write)
MEME_fasta_5utr_to_write <- c(paste0(c('>WEIGHTS', MEME_fasta_weight_5utr), collapse = ' '), MEME_fasta_5utr_to_write)

## write to file
write.table(x = MEME_fasta_3utr_to_write, file = 'meme_v5/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
write.table(x = MEME_fasta_5utr_to_write, file = 'meme_v5/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)

for (i in 1:(length(gene_id_ordered)-1)){
  first_gene_less <- gene_id_ordered[i+1]
  
  first_gene_less_idx <- which.min(sapply(MEME_fasta_3utr, function(x){grep(first_gene_less, x)}))
  
  MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr[1:(first_gene_less_idx-1)]
  MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_weight[1:((first_gene_less_idx-1)/2)]
  
  zero_length_3utr_more <- which(sapply(MEME_fasta_3utr_to_write_more, function(x){nchar(x)<len_th}))
  MEME_fasta_3utr_to_write_more <- MEME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
  MEME_fasta_3utr_to_write_more_weight <- MEME_fasta_3utr_to_write_more_weight[-c(zero_length_3utr_more/2)]
  
  MEME_fasta_3utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_more_weight), collapse = ' '), MEME_fasta_3utr_to_write_more)
  
  
  MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr[(first_gene_less_idx):length(MEME_fasta_3utr)]
  MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
  
  zero_length_3utr_less <- which(sapply(MEME_fasta_3utr_to_write_less, function(x){nchar(x)< len_th}))
  MEME_fasta_3utr_to_write_less <- MEME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
  MEME_fasta_3utr_to_write_less_weight <- MEME_fasta_3utr_to_write_less_weight[-c(zero_length_3utr_less/2)]
  
  MEME_fasta_3utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_3utr_to_write_less_weight), collapse = ' '), MEME_fasta_3utr_to_write_less)
  
  MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr[1:(first_gene_less_idx-1)]
  MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_weight[1:(first_gene_less_idx-1)/2]
  
  zero_length_5utr_more <- which(sapply(MEME_fasta_5utr_to_write_more, function(x){nchar(x) < len_th}))
  MEME_fasta_5utr_to_write_more <- MEME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
  MEME_fasta_5utr_to_write_more_weight <- MEME_fasta_5utr_to_write_more_weight[-c(zero_length_5utr_more/2)]
  
  MEME_fasta_5utr_to_write_more <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_more_weight), collapse = ' '), MEME_fasta_5utr_to_write_more)
  
  
  MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr[(first_gene_less_idx):length(MEME_fasta_5utr)]
  MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_weight[((first_gene_less_idx-1)/2+1):length(MEME_fasta_weight)]
  
  zero_length_5utr_less <- which(sapply(MEME_fasta_5utr_to_write_less, function(x){nchar(x) < len_th}))
  MEME_fasta_5utr_to_write_less <- MEME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
  MEME_fasta_5utr_to_write_less_weight <- MEME_fasta_5utr_to_write_less_weight[-c(zero_length_5utr_less/2)]
  
  MEME_fasta_5utr_to_write_less <- c(paste0(c('>WEIGHTS', MEME_fasta_5utr_to_write_less_weight), collapse = ' '), MEME_fasta_5utr_to_write_less)
  
  write.table(x = MEME_fasta_3utr_to_write_more, file = paste0('meme_v5/3utr/meme_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = MEME_fasta_3utr_to_write_less, file = paste0('meme_v5/3utr/meme_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = MEME_fasta_5utr_to_write_more, file = paste0('meme_v5/5utr/meme_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = MEME_fasta_5utr_to_write_less, file = paste0('meme_v5/5utr/meme_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
  
}


## AME_v5 main common part only
# order gene by transport value for all (more -> less)
gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]

AME_fasta_5utr <- AME_fasta_3utr <- c()

for (gene in gene_id_ordered){
  dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
  strand <- unique(dat_g$strand)
  chr <- unique(dat_g$chromosome_name)
  
  rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
  
  gene_name <- unique(dat_g$external_gene_name)
  print(gene_name)
  
  AME_fasta_3utr <- c(AME_fasta_3utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
  AME_fasta_5utr <- c(AME_fasta_5utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
  
  
  TSS_manual <- as.numeric(strsplit(manual_start_end$TSSs[manual_start_end$Geneid == gene ], ", ")[[1]])
  TES_manual <- as.numeric(strsplit(manual_start_end$TESs[manual_start_end$Geneid == gene ], ", ")[[1]])
  
  if(strand == 1){
    TSS_manual <- max(TSS_manual)
    TES_manual <- min(TES_manual)
  }else if (strand == -1){
    TSS_manual <- min(TSS_manual)
    TES_manual <- max(TES_manual)
  }
  
  start_5utr <- start_3utr <- c()
  end_5utr <- end_3utr <- c()
  
  # seq_5utr <- seq_3utr <- c()
  
  # get utr start and end for each transcript
  for (transcript in unique(dat_g$ensembl_transcript_id)){
    dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
    
    start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
    end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
    start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
    end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
    
    # order if not, weird case were it wasn't ordered
    order_5utr <- order(start_5utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
    start_5utr <- c(start_5utr, list(start_5utr_tmp[order_5utr]))
    end_5utr <- c(end_5utr, list(end_5utr_tmp[order_5utr]))
    order_3utr <- order(start_3utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
    start_3utr <- c(start_3utr, list(start_3utr_tmp[order_3utr]))
    end_3utr <- c(end_3utr, list(end_3utr_tmp[order_3utr]))
  }
  
  # remove utr exons outside of TSS-TES range or extend range
  if(strand == 1){
    for (l in 1:length(start_5utr)){
      if ( length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
        to_clip_or_remove <-  start_5utr[[l]] <= TSS_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & end_5utr[[l]] <= TSS_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          start_5utr[[l]][to_clip] <- TSS_manual
          start_5utr[[l]] <- start_5utr[[l]][!to_remove]
          end_5utr[[l]] <- end_5utr[[l]][!to_remove]
        }else{ # extend
          start_5utr[[l]][1] <- TSS_manual
        }
      }
    }
    for (l in 1:length(start_3utr)){
      if ( length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
        to_clip_or_remove <-  end_3utr[[l]] >= TES_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & start_3utr[[l]] >= TES_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          end_3utr[[l]][to_clip] <- TES_manual
          end_3utr[[l]] <- end_3utr[[l]][!to_remove]
          start_3utr[[l]] <- start_3utr[[l]][!to_remove]      
        }else{ # extend
          end_3utr[[l]][length(end_3utr[[l]])] <- TES_manual
        }
      }
    }
  }else if (strand == -1){
    for (l in 1:length(start_5utr)){
      if (length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
        to_clip_or_remove <-  end_5utr[[l]] >= TSS_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & start_5utr[[l]] >= TSS_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          end_5utr[[l]][to_clip] <- TSS_manual
          end_5utr[[l]] <- end_5utr[[l]][!to_remove]
          start_5utr[[l]] <- start_5utr[[l]][!to_remove]        
        }else{ # extend
          end_5utr[[l]][1] <- TSS_manual
        }
      }
    }
    
    for (l in 1:length(start_3utr)){
      if (length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
        to_clip_or_remove <-  start_3utr[[l]] <= TES_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & end_3utr[[l]] <= TES_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          start_3utr[[l]][to_clip] <- TES_manual
          start_3utr[[l]] <- start_3utr[[l]][!to_remove]
          end_3utr[[l]] <- end_3utr[[l]][!to_remove]
          merge_delim_3utr <- merge_delim_5utr[!to_remove]
        }else{ # extend
          start_3utr[[l]][length(start_3utr[[l]])] <- TES_manual
        }
      }
    }
  }
  
  
  # get common intervals and get sequences and merge with 10 N (max motif width for RBP 8) when inconsistent between transcripts
  if (length(unique(dat_g$ensembl_transcript_id)) > 1){
    if(strand == 1){
      dec <- F
    }else if (strand == -1){
      dec <- T
    }
    
    if(any(unlist(lapply(start_3utr, function(l){length(l)})) == 0)){
      start_3utr_merged <- c()
      end_3utr_merged <- c()
      merge_delim_3utr <- c()
    }else{
      intervals_3utr <- do.call("rbind", lapply(1:length(start_3utr), function(i){data.frame(start=start_3utr[[i]], end=end_3utr[[i]], "transcript"=i)}))
      intervals_3utr <- intervals_3utr[order(intervals_3utr$start, decreasing = dec),]
      
      tmp <- merge_intervals(intervals_3utr, strand)
      
      start_3utr_merged <- tmp[[1]]$start
      end_3utr_merged <- tmp[[1]]$end
      merge_delim_3utr <- tmp[[2]]
    }
    
    if(any(unlist(lapply(start_5utr, function(l){length(l)})) == 0)){
      start_5utr_merged <- c()
      end_5utr_merged <- c()
      merge_delim_5utr <- c()
    }else{
      intervals_5utr <- do.call("rbind", lapply(1:length(start_5utr), function(i){data.frame(start=start_5utr[[i]], end=end_5utr[[i]], "transcript"=i)}))
      intervals_5utr <- intervals_5utr[order(intervals_5utr$start, decreasing = dec),]
      
      tmp <- merge_intervals(intervals_5utr, strand)
      
      start_5utr_merged <- tmp[[1]]$start
      end_5utr_merged <- tmp[[1]]$end
      merge_delim_5utr <- tmp[[2]]
    }
    
  }else{
    start_5utr_merged <- unlist(start_5utr)
    end_5utr_merged <- unlist(end_5utr)
    start_3utr_merged <- unlist(start_3utr)
    end_3utr_merged <- unlist(end_3utr)
    
    merge_delim_5utr <- rep("",length(start_5utr_merged))
    
    merge_delim_3utr <- rep("",length(start_3utr_merged))
    
  }
  
  seq_5utr <- c()
  seq_3utr <- c()
  # get sequences and merge with 10 N (max motif width for RBP 8)
  # In GET take - strand because RNA is reversed, to get real motifs
  if(strand == 1){
    if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
      for (i in 1:length(start_5utr_merged)){
        if(start_5utr_merged[i] >= TSS_manual){
          ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_5utr <- c(seq_5utr, content(r)$seq)
        }
      }
    }
    if( length(start_3utr_merged) >0 & !is.na(TES_manual)){
      for (i in 1:length(start_3utr_merged)){
        if(start_3utr_merged[i] <= TES_manual){
          ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..', end_3utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_3utr <- c(seq_3utr, content(r)$seq)
        }
      }
    }
  }else if (strand == -1){
    if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
      for (i in 1:length(start_5utr_merged)){
        if(start_5utr_merged[i] <= TSS_manual){
          ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_5utr <- c(seq_5utr, content(r)$seq)
        }
      }
    }
    if( length(start_3utr_merged) > 0 & !is.na(TES_manual) ){
      for (i in 1:length(start_3utr_merged)){
        if(end_3utr_merged[i] >= TES_manual){
          ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..',end_3utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_3utr <- c(seq_3utr, content(r)$seq)
        }
      }
    }
  }
  
  # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
  # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
  AME_fasta_3utr <- c(AME_fasta_3utr, dna_to_rna(paste0(seq_3utr, merge_delim_3utr, sep="", collapse = '')))
  AME_fasta_5utr <- c(AME_fasta_5utr, dna_to_rna(paste0(seq_5utr, merge_delim_5utr, sep="", collapse = '')))
  
  warnings()
  # tmp <- readline()
  # if(tmp == "c"){next}else{eval(tmp)}
}

# remove zero length from file and weight
zero_length_3utr <- which(sapply(AME_fasta_3utr, function(x){nchar(x)==0}))
AME_fasta_3utr_to_write <- AME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]

zero_length_5utr <- which(sapply(AME_fasta_5utr, function(x){nchar(x)==0}))
AME_fasta_5utr_to_write <- AME_fasta_5utr[-c(zero_length_5utr-1,zero_length_5utr)]

## write to file reverse (increase order => low value = positive)
AME_fasta_3utr_reverse_to_write <- AME_fasta_3utr_to_write[rep(2*((length(AME_fasta_3utr_to_write)/2):1),each=2) + c(-1,0)]
AME_fasta_5utr_reverse_to_write <- AME_fasta_5utr_to_write[rep(2*((length(AME_fasta_5utr_to_write)/2):1),each=2) + c(-1,0)]

## write to file
write.table(x = AME_fasta_3utr_reverse_to_write, file = 'ame_v5/3utr/meme_3utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
write.table(x = AME_fasta_5utr_reverse_to_write, file = 'ame_v5/5utr/meme_5utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)

# Invert score to do (high value=positive)
for (i in 1:(length(AME_fasta_3utr_to_write)/2)){
  tmp <- strsplit(AME_fasta_3utr_to_write[2*i-1],split = " ")
  tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
  AME_fasta_3utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
}
for (i in 1:(length(AME_fasta_5utr_to_write)/2)){
  tmp <- strsplit(AME_fasta_5utr_to_write[2*i-1],split = " ")
  tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
  AME_fasta_5utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
}
## write to file
write.table(x = AME_fasta_3utr_to_write, file = 'ame_v5/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
write.table(x = AME_fasta_5utr_to_write, file = 'ame_v5/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)

# manual partition
for (i in 1:(length(gene_id_ordered)-1)){
  first_gene_less <- gene_id_ordered[i+1]
  
  first_gene_less_idx <- which.min(sapply(AME_fasta_3utr, function(x){grep(first_gene_less, x)}))
  
  AME_fasta_3utr_to_write_more <- AME_fasta_3utr[1:(first_gene_less_idx-1)]
  
  zero_length_3utr_more <- which(sapply(AME_fasta_3utr_to_write_more, function(x){nchar(x)==0}))
  AME_fasta_3utr_to_write_more <- AME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
  
  AME_fasta_3utr_to_write_less <- AME_fasta_3utr[(first_gene_less_idx):length(AME_fasta_3utr)]
  
  zero_length_3utr_less <- which(sapply(AME_fasta_3utr_to_write_less, function(x){nchar(x)==0}))
  AME_fasta_3utr_to_write_less <- AME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
  
  AME_fasta_5utr_to_write_more <- AME_fasta_5utr[1:(first_gene_less_idx-1)]
  
  zero_length_5utr_more <- which(sapply(AME_fasta_5utr_to_write_more, function(x){nchar(x)==0}))
  AME_fasta_5utr_to_write_more <- AME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
  
  AME_fasta_5utr_to_write_less <- AME_fasta_5utr[(first_gene_less_idx):length(AME_fasta_5utr)]
  
  zero_length_5utr_less <- which(sapply(AME_fasta_5utr_to_write_less, function(x){nchar(x)==0}))
  AME_fasta_5utr_to_write_less <- AME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
  
  write.table(x = AME_fasta_3utr_to_write_more, file = paste0('ame_v5/3utr/ame_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = AME_fasta_3utr_to_write_less, file = paste0('ame_v5/3utr/ame_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = AME_fasta_5utr_to_write_more, file = paste0('ame_v5/5utr/ame_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = AME_fasta_5utr_to_write_less, file = paste0('ame_v5/5utr/ame_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
  
}



## AME_v5 main isoform only
# order gene by transport value for all (more -> less)
gene_id_ordered <- rel_transport_all$GeneId[order(rel_transport_all$Transport, decreasing = T)]

AME_fasta_5utr <- AME_fasta_3utr <- c()

for (gene in gene_id_ordered){
  dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
  strand <- unique(dat_g$strand)
  chr <- unique(dat_g$chromosome_name)
  
  rel_transport <- rel_transport_all$Transport[rel_transport_all$GeneId == gene]
  
  gene_name <- unique(dat_g$external_gene_name)
  print(gene_name)
  
  AME_fasta_3utr <- c(AME_fasta_3utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
  AME_fasta_5utr <- c(AME_fasta_5utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
  
  
  TSS_manual <- as.numeric(strsplit(manual_start_end$TSSs[manual_start_end$Geneid == gene ], ", ")[[1]])[1]
  TES_manual <- as.numeric(strsplit(manual_start_end$TESs[manual_start_end$Geneid == gene ], ", ")[[1]])[1]
  
  start_5utr <- start_3utr <- c()
  end_5utr <- end_3utr <- c()
  
  transcript <- substr(strsplit(manual_start_end$Isoforms[manual_start_end$Geneid == gene ], ", ")[[1]][1],1,18)
  
  dat_t <- dat_g[dat_g$ensembl_transcript_id == transcript,]
  start_5utr_tmp <- dat_t$`5_utr_start`[!is.na(dat_t$`5_utr_start`)]
  end_5utr_tmp <- dat_t$`5_utr_end`[!is.na(dat_t$`5_utr_end`)]
  start_3utr_tmp <- dat_t$`3_utr_start`[!is.na(dat_t$`3_utr_start` )]
  end_3utr_tmp <- dat_t$`3_utr_end`[! is.na(dat_t$`3_utr_end`)]
    
  # order if not, weird case were it wasn't ordered
  order_5utr <- order(start_5utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
  start_5utr <- c(start_5utr, list(start_5utr_tmp[order_5utr]))
  end_5utr <- c(end_5utr, list(end_5utr_tmp[order_5utr]))
  order_3utr <- order(start_3utr_tmp, decreasing = ifelse(strand==-1,TRUE,FALSE))
  start_3utr <- c(start_3utr, list(start_3utr_tmp[order_3utr]))
  end_3utr <- c(end_3utr, list(end_3utr_tmp[order_3utr]))
  
  
  # remove utr exons outside of TSS-TES range or extend range
  if(strand == 1){
    for (l in 1:length(start_5utr)){
      if ( length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
        to_clip_or_remove <-  start_5utr[[l]] <= TSS_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & end_5utr[[l]] <= TSS_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          start_5utr[[l]][to_clip] <- TSS_manual
          start_5utr[[l]] <- start_5utr[[l]][!to_remove]
          end_5utr[[l]] <- end_5utr[[l]][!to_remove]
        }else{ # extend
          start_5utr[[l]][1] <- TSS_manual
        }
      }
    }
    for (l in 1:length(start_3utr)){
      if ( length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
        to_clip_or_remove <-  end_3utr[[l]] >= TES_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & start_3utr[[l]] >= TES_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          end_3utr[[l]][to_clip] <- TES_manual
          end_3utr[[l]] <- end_3utr[[l]][!to_remove]
          start_3utr[[l]] <- start_3utr[[l]][!to_remove]      
        }else{ # extend
          end_3utr[[l]][length(end_3utr[[l]])] <- TES_manual
        }
      }
    }
  }else if (strand == -1){
    for (l in 1:length(start_5utr)){
      if (length(start_5utr[[l]]) > 0 & !is.na(TSS_manual)){
        to_clip_or_remove <-  end_5utr[[l]] >= TSS_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & start_5utr[[l]] >= TSS_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          end_5utr[[l]][to_clip] <- TSS_manual
          end_5utr[[l]] <- end_5utr[[l]][!to_remove]
          start_5utr[[l]] <- start_5utr[[l]][!to_remove]        
        }else{ # extend
          end_5utr[[l]][1] <- TSS_manual
        }
      }
    }
    
    for (l in 1:length(start_3utr)){
      if (length(start_3utr[[l]]) > 0 &  !is.na(TES_manual)){
        to_clip_or_remove <-  start_3utr[[l]] <= TES_manual
        if( any(to_clip_or_remove)) {
          to_remove <- to_clip_or_remove & end_3utr[[l]] <= TES_manual
          to_clip <- to_clip_or_remove & !to_remove
          
          start_3utr[[l]][to_clip] <- TES_manual
          start_3utr[[l]] <- start_3utr[[l]][!to_remove]
          end_3utr[[l]] <- end_3utr[[l]][!to_remove]
          merge_delim_3utr <- merge_delim_5utr[!to_remove]
        }else{ # extend
          start_3utr[[l]][length(start_3utr[[l]])] <- TES_manual
        }
      }
    }
  }
  
  
  start_5utr_merged <- unlist(start_5utr)
  end_5utr_merged <- unlist(end_5utr)
  start_3utr_merged <- unlist(start_3utr)
  end_3utr_merged <- unlist(end_3utr)
    
  merge_delim_5utr <- rep("",length(start_5utr_merged))
  
  merge_delim_3utr <- rep("",length(start_3utr_merged))
  
  seq_5utr <- c()
  seq_3utr <- c()
  # get sequences and merge with 10 N (max motif width for RBP 8)
  # In GET take - strand because RNA is reversed, to get real motifs
  if(strand == 1){
    if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
      for (i in 1:length(start_5utr_merged)){
        if(start_5utr_merged[i] >= TSS_manual){
          ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_5utr <- c(seq_5utr, content(r)$seq)
        }
      }
    }
    if( length(start_3utr_merged) >0 & !is.na(TES_manual)){
      for (i in 1:length(start_3utr_merged)){
        if(start_3utr_merged[i] <= TES_manual){
          ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..', end_3utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_3utr <- c(seq_3utr, content(r)$seq)
        }
      }
    }
  }else if (strand == -1){
    if( length(start_5utr_merged) >0 & !is.na(TSS_manual)){
      for (i in 1:length(start_5utr_merged)){
        if(start_5utr_merged[i] <= TSS_manual){
          ext5 <- paste0("/sequence/region/mouse/", chr, ":", start_5utr_merged[i],'..',end_5utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext5, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_5utr <- c(seq_5utr, content(r)$seq)
        }
      }
    }
    if( length(start_3utr_merged) > 0 & !is.na(TES_manual) ){
      for (i in 1:length(start_3utr_merged)){
        if(end_3utr_merged[i] >= TES_manual){
          ext3 <- paste0("/sequence/region/mouse/", chr, ":", start_3utr_merged[i],'..',end_3utr_merged[i],':', strand,"?mask=soft")
          r <- GET(paste(server, ext3, sep = ""), content_type("application/json"))
          stop_for_status(r)
          seq_3utr <- c(seq_3utr, content(r)$seq)
        }
      }
    }
  }
  
  # MEME_fasta_3utr <- c(MEME_fasta_3utr, paste0(seq_3utr, collapse = '-'))
  # MEME_fasta_5utr <- c(MEME_fasta_5utr, paste0(seq_5utr, collapse = '-'))
  AME_fasta_3utr <- c(AME_fasta_3utr, dna_to_rna(paste0(seq_3utr, merge_delim_3utr, sep="", collapse = '')))
  AME_fasta_5utr <- c(AME_fasta_5utr, dna_to_rna(paste0(seq_5utr, merge_delim_5utr, sep="", collapse = '')))
  
  warnings()
  # tmp <- readline()
  # if(tmp == "c"){next}else{eval(tmp)}
}

# remove zero length from file and weight
zero_length_3utr <- which(sapply(AME_fasta_3utr, function(x){nchar(x)==0}))
if(length(zero_length_3utr)!=0){
  AME_fasta_3utr_to_write <- AME_fasta_3utr[-c(zero_length_3utr-1,zero_length_3utr)]
}else{
  AME_fasta_3utr_to_write <- AME_fasta_3utr
}

zero_length_5utr <- which(sapply(AME_fasta_5utr, function(x){nchar(x)==0}))
AME_fasta_5utr_to_write <- AME_fasta_5utr[-c(zero_length_5utr-1,zero_length_5utr)]

## write to file reverse (increase order => low value = positive)
AME_fasta_3utr_reverse_to_write <- AME_fasta_3utr_to_write[rep(2*((length(AME_fasta_3utr_to_write)/2):1),each=2) + c(-1,0)]
AME_fasta_5utr_reverse_to_write <- AME_fasta_5utr_to_write[rep(2*((length(AME_fasta_5utr_to_write)/2):1),each=2) + c(-1,0)]

## write to file
write.table(x = AME_fasta_3utr_reverse_to_write, file = 'ame_v5_main/3utr/meme_3utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
write.table(x = AME_fasta_5utr_reverse_to_write, file = 'ame_v5_main/5utr/meme_5utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)

# Invert score to do (high value=positive)
for (i in 1:(length(AME_fasta_3utr_to_write)/2)){
  tmp <- strsplit(AME_fasta_3utr_to_write[2*i-1],split = " ")
  tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
  AME_fasta_3utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
}
for (i in 1:(length(AME_fasta_5utr_to_write)/2)){
  tmp <- strsplit(AME_fasta_5utr_to_write[2*i-1],split = " ")
  tmp[[1]][2] <- -as.numeric(tmp[[1]][2])
  AME_fasta_5utr_to_write[2*i-1] <- paste(tmp[[1]], collapse=' ')
}
## write to file
write.table(x = AME_fasta_3utr_to_write, file = 'ame_v5_main/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
write.table(x = AME_fasta_5utr_to_write, file = 'ame_v5_main/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)


# manual partition
for (i in 1:(length(gene_id_ordered)-1)){
  first_gene_less <- gene_id_ordered[i+1]
  
  first_gene_less_idx <- which.min(sapply(AME_fasta_3utr, function(x){grep(first_gene_less, x)}))
  
  AME_fasta_3utr_to_write_more <- AME_fasta_3utr[1:(first_gene_less_idx-1)]
  zero_length_3utr_more <- which(sapply(AME_fasta_3utr_to_write_more, function(x){nchar(x)==0}))
  if(length(zero_length_3utr_more)!=0){
    AME_fasta_3utr_to_write_more <- AME_fasta_3utr_to_write_more[-c(zero_length_3utr_more-1,zero_length_3utr_more)]
  }
  
  AME_fasta_3utr_to_write_less <- AME_fasta_3utr[(first_gene_less_idx):length(AME_fasta_3utr)]
  zero_length_3utr_less <- which(sapply(AME_fasta_3utr_to_write_less, function(x){nchar(x)==0}))
  if(length(zero_length_3utr_less)!=0){
    AME_fasta_3utr_to_write_less <- AME_fasta_3utr_to_write_less[-c(zero_length_3utr_less-1,zero_length_3utr_less)]
  }
  
  AME_fasta_5utr_to_write_more <- AME_fasta_5utr[1:(first_gene_less_idx-1)]
  zero_length_5utr_more <- which(sapply(AME_fasta_5utr_to_write_more, function(x){nchar(x)==0}))
  if(length(zero_length_5utr_more)!=0){
    AME_fasta_5utr_to_write_more <- AME_fasta_5utr_to_write_more[-c(zero_length_5utr_more-1,zero_length_5utr_more)]
  }
  
  AME_fasta_5utr_to_write_less <- AME_fasta_5utr[(first_gene_less_idx):length(AME_fasta_5utr)]
  zero_length_5utr_less <- which(sapply(AME_fasta_5utr_to_write_less, function(x){nchar(x)==0}))
  if(length(zero_length_5utr_less)!=0){
    AME_fasta_5utr_to_write_less <- AME_fasta_5utr_to_write_less[-c(zero_length_5utr_less-1,zero_length_5utr_less)]
  }
  
  write.table(x = AME_fasta_3utr_to_write_more, file = paste0('ame_v5_main/3utr/ame_3utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = AME_fasta_3utr_to_write_less, file = paste0('ame_v5_main/3utr/ame_3utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = AME_fasta_5utr_to_write_more, file = paste0('ame_v5_main/5utr/ame_5utr_fasta_more_',i, '.txt'), row.names = F, quote = F, col.names = F)
  write.table(x = AME_fasta_5utr_to_write_less, file = paste0('ame_v5_main/5utr/ame_5utr_fasta_less_',i, '.txt'), row.names = F, quote = F, col.names = F)
  
}








