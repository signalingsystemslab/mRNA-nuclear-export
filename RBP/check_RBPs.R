library('httr')
server <- "https://may2017.rest.ensembl.org"

# check RBP
manual_start_end <- read.csv('../Manual Curation/Curated list_naive.csv', header = T, sep = '\t')
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

load(file = '../Modeling/all_results.Rdata')

param_best_all <- as.data.frame(Reduce(cbind,lapply(as.data.frame(all_best_param$naive$all), function(x){unlist(x)})))
colnames(param_best_all) <- colnames(all_best_param$naive$all)
param_best_all$GeneId <- unique(gene_UTR_start_end$ensembl_gene_id[match(row.names(param_best_all), gene_UTR_start_end$external_gene_name)])

# order gene by transport value for all (more -> less)
gene_id_ordered <- param_best_all$GeneId[order(param_best_all$`k1'k2'/k2`, decreasing = T)]

## AME_v5 main isoform only
# order gene by transport value for all (more -> less)
gene_id_ordered <- param_best_all$GeneId[order(param_best_all$`k1'k2'/k2`, decreasing = T)]

AME_fasta_5utr <- AME_fasta_3utr <- c()

for (gene in gene_id_ordered){
  dat_g <- gene_UTR_start_end[gene_UTR_start_end$ensembl_gene_id == substr(gene,1,18),]
  strand <- unique(dat_g$strand)
  chr <- unique(dat_g$chromosome_name)
  
  rel_transport <- param_best_all$`k1'k2'/k2`[param_best_all$GeneId == gene]
  
  gene_name <- unique(dat_g$external_gene_name)
  print(gene_name)
  
  AME_fasta_3utr <- c(AME_fasta_3utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
  AME_fasta_5utr <- c(AME_fasta_5utr, paste0('>', gene , ' ', rel_transport, ' ', gene_name, ' ', strand))
  
  
  TSS_manual <- as.numeric(strsplit(manual_start_end$TSSs[substr(manual_start_end$Geneid,1,18) == gene ], ", ")[[1]])[1]
  TES_manual <- as.numeric(strsplit(manual_start_end$TESs[substr(manual_start_end$Geneid,1,18) == gene ], ", ")[[1]])[1]
  
  start_5utr <- start_3utr <- c()
  end_5utr <- end_3utr <- c()
  
  transcript <- substr(strsplit(manual_start_end$Isoforms[substr(manual_start_end$Geneid,1,18) == gene ], ", ")[[1]][1],1,18)
  
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
dir.create('ame_main/3utr', showWarnings = F, recursive = T)
dir.create('ame_main/5utr', showWarnings = F, recursive = T)

write.table(x = AME_fasta_3utr_reverse_to_write, file = 'ame_main/3utr/meme_3utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)
write.table(x = AME_fasta_5utr_reverse_to_write, file = 'ame_main/5utr/meme_5utr_fasta_all_reverse.txt', row.names = F, quote = F, col.names = F)

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
write.table(x = AME_fasta_3utr_to_write, file = 'ame_main/3utr/meme_3utr_fasta_all.txt', row.names = F, quote = F, col.names = F)
write.table(x = AME_fasta_5utr_to_write, file = 'ame_main/5utr/meme_5utr_fasta_all.txt', row.names = F, quote = F, col.names = F)







