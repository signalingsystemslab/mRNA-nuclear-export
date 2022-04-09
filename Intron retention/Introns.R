library(biomaRt)

#### Process Intron data #####---------------------------------------------------------
PI_naive_rep1 <- read.table('SQUID/Naive/Nuc.Naive.rep1/Result/intron_PI.txt.gz', header=T, stringsAsFactors = F)
PI_naive_rep2 <- read.table('SQUID/Naive/Nuc.Naive.rep2/Result/intron_PI.txt.gz', header=T, stringsAsFactors = F)
PI_lpa_rep1 <- read.table('SQUID/LPA/Nuc.LPA.rep1/Result/intron_PI.txt.gz', header=T, stringsAsFactors = F)
PI_lpa_rep2 <- read.table('SQUID/LPA/Nuc.LPA.rep1/Result/intron_PI.txt.gz', header=T, stringsAsFactors = F)

ncol <- max(count.fields('Intron_transcript.txt', sep = "\t"))
intron_transcript_mapping <- read.table('Intron_transcript.txt', header=F, col.names=paste0('V', seq_len(ncol)), skip = 0, fill = T, stringsAsFactors = F)

PI_naive_rep1$Time <- paste0(c(10,15,20,25,30,35,40,50,60,90,120), collapse=",")
PI_naive_rep2$Time <- paste0(c(0,10,15,20,25,30,40,50,60,75,90,120), collapse=",")
PI_lpa_rep1$Time   <- paste0(c(0,10,15,20,25,30,35,40,50,60,90,120), collapse=",")
PI_lpa_rep2$Time   <- paste0(c(0,10,15,20,25,30,40,50,60,75,90,120), collapse=",")

# get annotations
ensembl89 <- useMart(biomart = 'ensembl', host="may2017.archive.ensembl.org", verbose=T)
mm_89 <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl89)

# select intron belonging to observed transcripts of induced genes
curated_list_naive <- read.table(file = "../Manual Curation/Curated list_naive.csv", header=T, sep ="\t")
all_transcript_naive <- unlist(strsplit(curated_list_naive$Isoforms,split = ", ", fixed=T))
naive_to_keep <- apply(intron_transcript_mapping, 1,function(x){tmp <- unlist(x[-1]); return(any(tmp %in% all_transcript_naive))})
intron_transcript_mapping_filtered_naive <- intron_transcript_mapping[naive_to_keep,]
intron_transcript_mapping_filtered_naive <- t(apply(intron_transcript_mapping_filtered_naive,1,function(x){tmp <- c(x[1], x[-1][x[-1] %in% all_transcript_naive]); tmp <- c(tmp, rep("", length(x)-length(tmp))); return(tmp)}))
intron_transcript_mapping_filtered_naive <- intron_transcript_mapping_filtered_naive[, !apply(intron_transcript_mapping_filtered_naive, 2, function(x){all(x=="")})]
# add geneid to introns:
intron_transcript_mapping_filtered_naive <- data.frame(Intrond_id = intron_transcript_mapping_filtered_naive[,1], 
                                                       Gene_id = curated_list_naive$Geneid[sapply(intron_transcript_mapping_filtered_naive[,2],function(x){grep(x,curated_list_naive$Isoforms)})],
                                                       Gene_name = curated_list_naive$Name[sapply(intron_transcript_mapping_filtered_naive[,2],function(x){grep(x,curated_list_naive$Isoforms)})],
                                                       intron_transcript_mapping_filtered_naive[,-1])

intron_annotated <- c()
# add intron attributes based on observed transcript: 
for (g in unique(intron_transcript_mapping_filtered_naive$Gene_id)){
  dat_g <- intron_transcript_mapping_filtered_naive[intron_transcript_mapping_filtered_naive$Gene_id == g,]
  transcript <- unique(unlist(dat_g[,-c(1:3)])[unlist(dat_g[,-c(1:3)])!=""])
  intron_start <- as.numeric(unlist(lapply(strsplit(dat_g$Intrond_id,split="_"),function(l){l[2]})))
  intron_end <- as.numeric(unlist(lapply(strsplit(dat_g$Intrond_id,split="_"),function(l){l[3]})))
  
  # get transcript_annotation and adjust with observed TSS and TES
  tmp <- getBM(mart = mm_89, attributes = c("ensembl_gene_id","ensembl_transcript_id", "strand", "ensembl_exon_id", "rank","exon_chrom_start", "exon_chrom_end"), filters = "ensembl_gene_id", values=substr(g,1,18))
  tmp <- tmp[tmp$ensembl_transcript_id %in% substr(transcript,1,18),]
  
  dat_t <- c()
  for( t_id in unique(tmp$ensembl_transcript_id) ){
    curated_index <- which(substr(unlist(strsplit(curated_list_naive$Isoforms[curated_list_naive$Geneid==g], ", ")),1,18) == t_id)
    TSS <- as.numeric(unlist(strsplit(curated_list_naive$TSSs[curated_list_naive$Geneid==g], ", "))[curated_index])
    TES <- as.numeric(unlist(strsplit(curated_list_naive$TESs[curated_list_naive$Geneid==g], ", "))[curated_index])
    
    dat_t_tmp <- tmp[tmp$ensembl_transcript_id == t_id, ]
    
    if(is.na(TSS)){
      TSS <- ifelse(unique(tmp$strand) == 1, dat_t_tmp[dat_t_tmp$rank == 1, "exon_chrom_start"], dat_t_tmp[dat_t_tmp$rank == 1, "exon_chrom_end"])
    }
    
    if( unique(tmp$strand) == 1 ){
      # remove exons if needed but it shouldn't (Lif does...)
      if(any(dat_t_tmp$exon_chrom_end < TSS)){
        print("removing exons")
        print(g)
        
        dat_t_tmp$rank <- dat_t_tmp$rank - sum(dat_t_tmp$exon_chrom_end < TSS)
        dat_t_tmp <- dat_t_tmp[dat_t_tmp$exon_chrom_end >= TSS,]
      }
      
      dat_t_tmp[dat_t_tmp$rank == 1, "exon_chrom_start"] <- TSS
      
      if(any(dat_t_tmp$exon_chrom_start > TES)){
        print("removing exons")
        print(g)
        dat_t_tmp$rank <- dat_t_tmp$rank - sum(dat_t_tmp$exon_chrom_start > TES)
        dat_t_tmp <- dat_t_tmp[dat_t_tmp$exon_chrom_end >= TSS,]
      }
      
      dat_t_tmp[dat_t_tmp$rank == max(dat_t_tmp$rank), "exon_chrom_end"] <- TES
      
    }else if ( unique(tmp$strand) == -1 ){
      # remove exons if needed but it shouldn't
      if(any(dat_t_tmp$exon_chrom_start > TSS)){
        print("removing exons")
        print(g)
        dat_t_tmp$rank <- dat_t_tmp$rank - sum(dat_t_tmp$exon_chrom_start > TSS)
        dat_t_tmp <- dat_t_tmp[dat_t_tmp$exon_chrom_start <= TSS,]
      }
      
      dat_t_tmp[dat_t_tmp$rank == 1, "exon_chrom_end"] <- TSS
      
      if(any(dat_t_tmp$exon_chrom_end < TES)){
        print("removing exons")
        print(g)
        
        dat_t_tmp$rank <- dat_t_tmp$rank - sum(dat_t_tmp$exon_chrom_end < TES)
        dat_t_tmp <- dat_t_tmp[dat_t_tmp$exon_chrom_end >= TES,]
      }
      
      dat_t_tmp[dat_t_tmp$rank == max(dat_t_tmp$rank), "exon_chrom_start"] <- TES
    } 
    dat_t <- rbind(dat_t, dat_t_tmp)
  }
  
  # check if junction is in an exon of another transcript, i.e. could lead to issue in PI measure
  attribute_E <- sapply(intron_start, function(x){any(x >= tmp$exon_chrom_start & x <= tmp$exon_chrom_end)}) | sapply(intron_end, function(x){any(x >= tmp$exon_chrom_start & x <= tmp$exon_chrom_end)})
  
  # check if junction is share by another intron i.e. could lead to issue in PI measure
  attribute_I <- sapply(1:length(intron_start), function(i){any(intron_start[i] == intron_start & intron_end[i] != intron_end)})
  
  dat_g$Attribute <- paste(sapply(1:nrow(dat_g), function(i){tmp <- paste0(ifelse(attribute_I[i], "I", ""), ifelse(attribute_E[i], "E", "")); if(tmp == ""){tmp <- "U"}; return(tmp)}))
  intron_annotated <- rbind(intron_annotated, dat_g)
  
}

head(intron_annotated)
unique(intron_annotated$Attribute)

# restrict to genes with only U introns
g_keep <- c()
for(g in unique(PI_naive_rep1$Gene_id)){
  if(all(intron_annotated$Attribute[intron_annotated$Gene_id == g] == "U")){
    g_keep <- c(g_keep, g)
  }
}
PI_naive_U <- list(rep1 = PI_naive_rep1[PI_naive_rep1$Intron_id %in% intron_annotated$Intrond_id[intron_annotated$Gene_id %in% g_keep], ],
                   rep2 = PI_naive_rep2[PI_naive_rep2$Intron_id %in% intron_annotated$Intrond_id[intron_annotated$Gene_id %in% g_keep], ])

avg_late_tp <- c(90,120)
avg_early_tp <- c(0,10)

for (b in names(PI_naive_U)){
  Time <- as.numeric(strsplit(PI_naive_U[[b]]$Time, split=',')[[1]])
  PI_naive_U[[b]]$Total_Junction_Counts <- sapply(1:nrow(PI_naive_U[[b]]), function(i){paste0(as.numeric(strsplit(PI_naive_U[[b]]$Inclusion_counts[i], split = ',')[[1]]) + as.numeric(strsplit(PI_naive_U[[b]]$Skipping_counts_EE[i], split = ',')[[1]]), collapse= ",")})
  
  # average late time points by batch
  PI_naive_U[[b]]$PI_Junction_avg_late <- sapply(PI_naive_U[[b]]$PI_Junction, function(x){mean(as.numeric(strsplit(x, split = ',')[[1]][Time %in% avg_late_tp]), na.rm=T)})
  PI_naive_U[[b]]$PI_Junction_avg_late_count <- sapply(PI_naive_U[[b]]$Total_Junction_Counts, function(x){mean(as.numeric(strsplit(x, split = ',')[[1]][Time %in% avg_late_tp]), na.rm=T)})
  PI_naive_U[[b]]$PI_Junction_wavg_late <- sapply(1:nrow(PI_naive_U[[b]]), function(i){sum(as.numeric(strsplit(PI_naive_U[[b]]$PI_Junction[i], split = ',')[[1]][Time %in% avg_late_tp])*(as.numeric(strsplit(PI_naive_U[[b]]$Total_Junction_Counts[i], split = ',')[[1]][Time %in% avg_late_tp])/sum(as.numeric(strsplit(PI_naive_U[[b]]$Total_Junction_Counts[i], split = ',')[[1]][Time %in% avg_late_tp]))), na.rm=T)})
  
  # average early time points by batch
  PI_naive_U[[b]]$PI_Junction_avg_early <- sapply(PI_naive_U[[b]]$PI_Junction, function(x){mean(as.numeric(strsplit(x, split = ',')[[1]][Time %in% avg_early_tp]), na.rm=T)})
  PI_naive_U[[b]]$PI_Junction_avg_early_count <- sapply(PI_naive_U[[b]]$Total_Junction_Counts, function(x){mean(as.numeric(strsplit(x, split = ',')[[1]][Time %in% avg_early_tp]), na.rm=T)})
  PI_naive_U[[b]]$PI_Junction_wavg_early <- sapply(1:nrow(PI_naive_U[[b]]), function(i){sum(as.numeric(strsplit(PI_naive_U[[b]]$PI_Junction[i], split = ',')[[1]][Time %in% avg_early_tp])*(as.numeric(strsplit(PI_naive_U[[b]]$Total_Junction_Counts[i], split = ',')[[1]][Time %in% avg_early_tp])/sum(as.numeric(strsplit(PI_naive_U[[b]]$Total_Junction_Counts[i], split = ',')[[1]][Time %in% avg_early_tp]))), na.rm=T)})
}

# check plot

# calculate max of intron PI per genes
all_intron_summary_measure <- list()

count_threshold <- 10
for(b in names(PI_naive_U)){
  all_intron_summary_measure[[b]] <- list()
  
  for (g in unique(PI_naive_U[[b]]$Gene_id) ){
    b_g <- PI_naive_U[[b]][PI_naive_U[[b]]$Gene_id==g,]
    
    if(g == "ENSMUSG00000021298.8"){
      print(b_g)
    }
    
    if(all(b_g$PI_Junction_avg_late_count >= count_threshold)){
      
      PI_naive_avg_avg_late <- mean(b_g$PI_Junction_avg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      PI_naive_avg_wavg_late <- mean(b_g$PI_Junction_wavg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      
      PI_naive_sd_avg_late <- sd(b_g$PI_Junction_avg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      PI_naive_sd_wavg_late <- sd(b_g$PI_Junction_wavg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      
      PI_naive_max_avg_late <- max(b_g$PI_Junction_avg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      PI_naive_max_wavg_late <- max(b_g$PI_Junction_wavg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      
      PI_naive_prod_avg_late <- prod(1-b_g$PI_Junction_avg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      PI_naive_prod_wavg_late <- prod(1-b_g$PI_Junction_wavg_late[b_g$PI_Junction_avg_late_count >= count_threshold])
      
      tmp <- data.frame(GeneId=g, 
                        PI_naive_avg_avg_late = PI_naive_avg_avg_late,
                        PI_naive_avg_wavg_late = PI_naive_avg_wavg_late,
                        
                        PI_naive_sd_avg_late = PI_naive_sd_avg_late,
                        PI_naive_sd_wavg_late = PI_naive_sd_wavg_late,
                        
                        PI_naive_max_avg_late = PI_naive_max_avg_late,
                        PI_naive_max_wavg_late = PI_naive_max_wavg_late,
                        
                        PI_naive_prod_avg_late = PI_naive_prod_avg_late,
                        PI_naive_prod_wavg_late = PI_naive_prod_wavg_late
      )
      tmp[,unlist(lapply(tmp, is.infinite))] <- NA
      all_intron_summary_measure[[b]] <- rbind(all_intron_summary_measure[[b]], tmp)
    }
  }
}


all_intron_summary_measure[['all']] <- list()
for (g in intersect(PI_naive_U[['rep1']]$Gene_id, PI_naive_U[['rep2']]$Gene_id)){
  rep1_g <- PI_naive_U[['rep1']][PI_naive_U[[b]]$Gene_id==g,]
  rep1_g$PI_Junction_avg_late[rep1_g$PI_Junction_avg_late_count < count_threshold] <- NA
  rep1_g$PI_Junction_wavg_late[rep1_g$PI_Junction_wavg_late_count < count_threshold] <- NA
  
  rep2_g <- PI_naive_U[['rep2']][PI_naive_U[[b]]$Gene_id==g,]
  rep2_g$PI_Junction_avg_late[rep2_g$PI_Junction_avg_late_count < count_threshold] <- NA
  rep2_g$PI_Junction_wavg_late[rep2_g$PI_Junction_wavg_late_count < count_threshold] <- NA
  
  rep2_g <- rep2_g[match(rep1_g$Intron_id,rep2_g$Intron_id),]
  
  PI_naive_avg_avg_late <- mean(rowMeans(cbind(rep1_g$PI_Junction_avg_late, rep2_g$PI_Junction_avg_late)))
  PI_naive_avg_wavg_late <- mean(rowMeans(cbind(rep1_g$PI_Junction_wavg_late, rep2_g$PI_Junction_wavg_late)))
  
  PI_naive_sd_avg_late <- sd(c(rep1_g$PI_Junction_avg_late, rep2_g$PI_Junction_avg_late))
  PI_naive_sd_wavg_late <- sd(c(rep1_g$PI_Junction_wavg_late, rep2_g$PI_Junction_wavg_late))
  
  PI_naive_max_avg_late <- max(rowMeans(cbind(rep1_g$PI_Junction_avg_late, rep2_g$PI_Junction_avg_late))) 
  PI_naive_max_wavg_late <- max(rowMeans(cbind(rep1_g$PI_Junction_wavg_late, rep2_g$PI_Junction_wavg_late))) 
  
  PI_naive_prod_avg_late <- prod(1- rowMeans(cbind(rep1_g$PI_Junction_avg_late, rep2_g$PI_Junction_avg_late)))
  PI_naive_prod_wavg_late <- prod(1- rowMeans(cbind(rep1_g$PI_Junction_wavg_late, rep2_g$PI_Junction_wavg_late)))
  
  
  tmp <- data.frame(GeneId=g, 
                    
                    PI_naive_avg_avg_late = PI_naive_avg_avg_late,
                    PI_naive_avg_wavg_late = PI_naive_avg_wavg_late,
                    
                    PI_naive_sd_avg_late = PI_naive_sd_avg_late,
                    PI_naive_sd_wavg_late = PI_naive_sd_wavg_late,
                    
                    PI_naive_max_avg_late = PI_naive_max_avg_late,
                    PI_naive_max_wavg_late = PI_naive_max_wavg_late,
                    
                    PI_naive_prod_avg_late = PI_naive_prod_avg_late,
                    PI_naive_prod_wavg_late = PI_naive_prod_wavg_late)
  
  tmp[,unlist(lapply(tmp, is.infinite))] <- NA
  all_intron_summary_measure[['all']] <- rbind(all_intron_summary_measure[['all']], tmp)
}

load('../Data/exons_rpkms/npRNA_Naive_exons_rpkm_5kb.Rdata')
load('../Modeling/gene_infos_with_length.Rdata')

# add PI=0 for genes without introns
for (g in row.names(npRNA_Naive_exons_rpkm_5kb)){
  if (! g %in% intron_transcript_mapping_filtered_naive$Gene_id ){ # not in the intron_annotation file i.e. no intron for the observed transcripts
    tmp <- data.frame(GeneId=g, 
                      
                      PI_naive_avg_avg_late = 0,
                      PI_naive_avg_wavg_late = 0,
                      
                      PI_naive_sd_avg_late = NA,
                      PI_naive_sd_wavg_late = NA,
                      
                      PI_naive_max_avg_late = 0,
                      PI_naive_max_wavg_late = 0,
                      
                      PI_naive_prod_avg_late = 1,
                      PI_naive_prod_wavg_late = 1
    )
    for(b in names(all_intron_summary_measure)){
      all_intron_summary_measure[[b]] <- rbind(all_intron_summary_measure[[b]], tmp)
    }
  }
}

for (b in names(all_intron_summary_measure)){
  all_intron_summary_measure[[b]]$Name <- gene_infos$external_gene_name[match(substr(all_intron_summary_measure[[b]]$GeneId,1,18), gene_infos$ensembl_gene_id)]
}

save(all_intron_summary_measure, file = 'all_intron_summary_measure.Rdata')

# average each intron then calculate prod - for kevin
splicing_prob <- data.frame()
for (g in unique(all_intron_summary_measure[['rep1']]$GeneId, all_intron_summary_measure[['rep2']]$GeneId)){
  tmp <- data.frame("GeneId" = g, 'Name'= unique(all_intron_summary_measure[['rep1']]$Name, all_intron_summary_measure[['rep2']]$Name)[match(g, unique(all_intron_summary_measure[['rep1']]$GeneId, all_intron_summary_measure[['rep2']]$GeneId))])
  for (b in c('rep1', 'rep2', 'all')){
    if(g %in% all_intron_summary_measure[[b]]$GeneId){
      tmp <- cbind(tmp, all_intron_summary_measure[[b]]$PI_naive_prod_avg_late[all_intron_summary_measure[[b]]$GeneId == g])
    }else{
      tmp <- cbind(tmp,NA)
    }
  }
  
  colnames(tmp) <- c('GeneId', 'Name', 'rep1', 'rep2', 'Avg')
  splicing_prob <- rbind(splicing_prob, tmp)
}

saveRDS(splicing_prob, file='Splicing_proba.Rds')
