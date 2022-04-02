## load datasets for chromatin

load("../Data/raw_counts/caRNA_Naive.Rdata")

## load libraries
library(edgeR)
library(biomaRt)

## load ensembl database
ensembl <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl', host = 'may2017.archive.ensembl.org')

# run DEG function
run_DEG <- function(dat, Time, Batch, group, RPKM_th, timeth, FC_th, pval, adjust, direction){
  y <- DGEList(counts = dat[ , -c(1:6) ], group = group, genes = dat[ , 1:6 ])
  
  # remove low gene for estimation of expression
  keep <- apply(y$counts, 1, function(x){ any( x > 32 ) }) 
  y <- DGEList(counts = dat[ keep, -c(1:6) ], group = group, genes = dat[ keep, 1:6 ])
  y <- calcNormFactors(y)
  
  design <- model.matrix(~ group + Batch)
  disp <- estimateDisp(y, design, robust = TRUE)
  
  # keep high express for DEG only
  rpkms <- rpkm(disp, log = T)
  hist(rpkms)
  hist(apply(rpkms,1,max))
  
  keep <- apply(rpkms, 1, function(x){ sum( x >= RPKM_th ) >= 3}) 
  
  
  # plotBCV(disp)
  fit <- glmQLFit(disp[ keep, ], design, robust = TRUE)
  # plotQLDisp(fit)
  
  top_table <- data.frame( matrix(nrow = 0, ncol = 12, 
                                  dimnames = list(c(), c("Geneid", "Chr", "Start", "End", "Strand", "Length", "logFC", "logCPM", "LR", "PValue", "FDR", "Time"))),
                           stringsAsFactors = FALSE)
  
  for(i in  which(as.numeric(levels(group)) <= timeth & as.numeric(levels(group)) != 0)){ # choose correct coefs
    topi <-  glmQLFTest(fit, coef = i )
    tt <- topTags(topi, n = sum(keep), p.value = pval, sort.by = "logFC", adjust.method = adjust)
    
    
    if(nrow(tt)>0){
      tt$table <- rapply(tt$table, as.character, classes = "factor", how = "replace")
      
      for(j in 1:nrow(tt$table)){
        up <- tt$table$logFC[j] >= log2(FC_th)
        down <- tt$table$logFC[j] <= -log2(FC_th)
        if (direction == "up"){
          selection <- up
        }else if (direction == "down"){
          selection <- down
        }else if (direction == "both"){
          selection <- up | down
        }else{
          stop("Wrong direction argument. Choose one of up, down or both.")
        }
        
        if( selection ){
          if( !(tt$table$Geneid[j] %in% top_table$Geneid)){
            top_table <- rbind(top_table, c(tt$table[j,], Time = levels(group)[i]), stringsAsFactors = FALSE)
          }else if (tt$table$FDR[j] < top_table$FDR[as.character(top_table$Geneid) == as.character(tt$table$Geneid[j])]){
            top_table[as.character(top_table$Geneid) == as.character(tt$table$Geneid[j]),] <- c(tt$table[j,], Time = levels(group)[i])
          }
        }
      }
    }
  }
  return(top_table)
}

# filter genes function
filter_genes <- function(top_table, mart){
  genemap <- getBM(attributes = c("ensembl_gene_id",'version', "external_gene_name", "description", "gene_biotype"),
                   filters = "ensembl_gene_id",
                   values = unlist(lapply(strsplit(x = top_table$Geneid, split = ".", fixed = T),
                                          function(x){x[1]})),
                   mart = mart )
  
  # keep only protein coding genes
  protein_coding_id <- paste(genemap$ensembl_gene_id[genemap$gene_biotype == 'protein_coding'], genemap$version[genemap$gene_biotype == 'protein_coding'], sep='.')
  top_table_protein_coding <- top_table_up[top_table_up$Geneid %in% protein_coding_id,]
  
  top_table_protein_coding$Name <- genemap$external_gene_name[match(top_table_protein_coding$Geneid, paste(genemap$ensembl_gene_id, genemap$version, sep='.'))]
  
  
  # Remove predicted genes == GmXXXX even if protein coding
  top_table_protein_coding <- top_table_protein_coding[-grep("^Gm[0-9]+$", top_table_protein_coding$Name),] 
  
  return(top_table_protein_coding)
}

## DEG for naive
dat <- caRNA_Naive_all

Time <- as.numeric(sapply(colnames(dat)[-(1:6)], function(x){strsplit(x, fixed = T, split = ".")[[1]][3]}))
Batch <- sapply(colnames(dat)[-(1:6)], function(x){strsplit(x, fixed = T, split = ".")[[1]][4]})
group <- as.factor(Time)

top_table_up <- run_DEG(dat, Time, Batch, group, RPKM_th = log2(1), timeth = 40,  FC_th = 10, pval = 0.01, adjust = 'fdr', direction = "up")
dim(top_table_up)

top_table_up_Naive_all_edgeR_FC10_RPKM1_FDR0.01_protein_coding <- filter_genes(top_table_up, ensembl)
dim(top_table_up_Naive_all_edgeR_FC10_RPKM1_FDR0.01_protein_coding)
sort(top_table_up_Naive_all_edgeR_FC10_RPKM1_FDR0.01_protein_coding$Name)

## Export list
gene_id <- top_table_up_Naive_all_edgeR_FC10_RPKM1_FDR0.01_protein_coding$Geneid
                 
genemap <- getBM(attributes = c("ensembl_gene_id",'version', "external_gene_name", "description", "gene_biotype"),
                 filters = "ensembl_gene_id",
                 values = unlist(lapply(strsplit(x = gene_id, split = ".", fixed = T),
                                        function(x){x[1]})),
                 mart = ensembl )

top_table_summary <- data.frame(Geneid = gene_id)
top_table_summary$Name <- genemap$external_gene_name[match(unlist(lapply(strsplit(x = gene_id, split = ".", fixed = T),
                                                                         function(x){x[1]})), genemap$ensembl_gene_id)]
top_table_summary$Description <- genemap$description[match(unlist(lapply(strsplit(x = gene_id, split = ".", fixed = T),
                                                                         function(x){x[1]})), genemap$ensembl_gene_id)]
top_table_summary$Gene_type <- genemap$gene_biotype[match(unlist(lapply(strsplit(x = gene_id, split = ".", fixed = T),
                                                                        function(x){x[1]})), genemap$ensembl_gene_id)]

top_table_summary <- top_table_summary[order(top_table_summary$Name),]

write.csv(top_table_summary, file = 'top_table_summary.csv', row.names = F)

