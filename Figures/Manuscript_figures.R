##### Manuscript Figures #####------------------------------------------------------
library(RColorBrewer)
library(grid)
library(ggplot2)
library(openxlsx)
library(biomaRt)
library(compiler)
library(scales)
library(deSolve)
library(ComplexHeatmap)
library(circlize)
library(MASS)
library(pheatmap)
library(gridExtra)
library(plotly)
# library(grImport2)

# DEG slightly more stringent
list_genes_cleaned <- read.csv('../Manual Curation/gene_list.csv', sep=',')

# load Data on exons whole genes, ActD data
load('../Data/exons_rpkms/caRNA_Naive_exons_rpkm.Rdata')
load('../Data/exons_rpkms/npRNA_Naive_exons_rpkm.Rdata')
load('../Data/exons_rpkms/cytoRNA_Naive_exons_rpkm.Rdata')

load('../Data/exons_rpkms/caRNA_LPA_exons_rpkm.Rdata')
load('../Data/exons_rpkms/npRNA_LPA_exons_rpkm.Rdata')
load('../Data/exons_rpkms/cytoRNA_LPA_exons_rpkm.Rdata')

load('../Data/exons_cpms/caRNA_Naive_exons_cpm.Rdata')
load('../Data/exons_cpms/npRNA_Naive_exons_cpm.Rdata')
load('../Data/exons_cpms/cytoRNA_Naive_exons_cpm.Rdata')

load('../ActD/BMDM_dataset.rda')

load(file = '../Data/exons_rpkms/caRNA_Naive_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/npRNA_Naive_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/cytoRNA_Naive_exons_rpkm_5kb.Rdata')

load(file = '../Data/exons_rpkms/caRNA_LPA_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/npRNA_LPA_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/cytoRNA_LPA_exons_rpkm_5kb.Rdata')

# load gene_infos
load('../Modeling/gene_infos_with_length.Rdata')
load('../Data/raw_counts/lib.size.Rdata')

# load modeling function
source(file = '../Modeling/functions_weighting_for_smoothing_include_negbinom_error_in_model.R')

# Color for figures
col_ca <- "#ee6565ff"
col_np <- "#3daf58ff"
col_cyto <- "#41a8efff"

col_naive_rep1 <- "#FFD280FF"
col_naive_rep2 <- "#BF7C00FF"
col_naive_rep3 <- "#FFA540FF"

col_lpa_rep1 <- "#FF80C4FF"
col_lpa_rep2 <- "#BF0667FF"
col_lpa_rep3 <- "#FF00A7FF"

# get gene annotation
ensembl89 <- useMart(biomart = 'ensembl', host="may2017.archive.ensembl.org", verbose=T)
mm_89 <- useDataset(dataset = 'mmusculus_gene_ensembl', mart = ensembl89)
gene_list_mm89 <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "version", "start_position", "end_position"), filters = "ensembl_gene_id", values = substr(rownames(caRNA_Naive_exons_rpkm), 1,18), mart = mm_89) 


#### Get all the data, fits and result ####---------------------------------------------
# => save data
all_data <- list()
all_fit <- list()
all_fit_CI <- list(LB=list(), UB=list())
all_best_param <- list()
all_param <- list()
all_annot <- list()
all_library <- list()

wcor <- function(x,y,w){
  mxw <- sum(w*x, na.rm=T)/sum(w[!is.na(x)])
  myw <- sum(w*y, na.rm=T)/sum(w[!is.na(y)])
  
  wcovxy <- sum(w*(x-mxw)*(y-myw), na.rm=T)/sum(w[!is.na(x) &  !is.na(y)])
  wcovxx <- sum(w*(x-mxw)^2, na.rm=T)/sum(w[!is.na(x)])
  wcovyy <- sum(w*(y-myw)^2, na.rm=T)/sum(w[!is.na(y)])
  wcor <- wcovxy/sqrt(wcovxx * wcovyy)
  
  return(wcor)
}

lib <- lib_exons[,c("lib.size", "norm.factors")]
row.names(lib) <- lib_exons$samples

for (condition in c('naive','lpa')){
  all_data[[condition]] <- list()
  all_fit[[condition]] <- list()
  all_fit_CI[["LB"]][[condition]] <- list()
  all_fit_CI[["UB"]][[condition]] <- list()
  all_best_param[[condition]] <- list()
  all_library[[condition]] <- list()
  all_annot[[condition]] <- list()
  
  for (g in gene_infos$ensembl_gene_id[order(gene_infos$external_gene_name)]){
    # get caRNA experimental data, to know the genes done
    if (condition == 'naive'){
      raw_data <- caRNA_Naive_exons_rpkm_5kb
    }else if(condition == 'lpa'){
      raw_data <- caRNA_LPA_exons_rpkm_5kb
    }
    
    ######## Do that before end and simplify!!!!
    # get library size for condition 
    lib_mat <- as.matrix(lib[grep(condition, row.names(lib), ignore.case = T),])
    row.names(lib_mat) <- gsub(pattern = condition, replacement = '', x = row.names(lib_mat), ignore.case = T)
    row.names(lib_mat) <- gsub(pattern = '..', replacement = '.', x = row.names(lib_mat), fixed = T)
    
    
    if (g %in% sapply(row.names(raw_data), function(s){substr(s,1,18)})){
      i <- unname(which(sapply(row.names(raw_data), function(s){substr(s,1,18)}) == g ))
      g_id <- rownames(raw_data)[i]
      
      gene_name <- unique(gene_infos$external_gene_name[gene_infos$ensembl_gene_id == substr(g_id,1,18)])
      
      # Load best param result if exists
      file <- paste0('../Modeling/Results/Optimisation/BFGS/',condition, '/param_bfgs_1000ri_merged_', gene_name, '_Finished.Rdata')
      if (file.exists(file)){
        load(file)
        bfgs <- bfgs_merged
      }else{
        next
      }
      
      print(paste0(gene_name,' - ',condition))
      
      # get full experimental data and associated gene length for neg-binom
      if (condition == "naive") {
        data_all <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_Naive_exons_rpkm_5kb,   g_id, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                          npRNA   = as.matrix(create_dat_g_rpkm(npRNA_Naive_exons_rpkm_5kb,   g_id, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                          cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_Naive_exons_rpkm_5kb, g_id, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
        gene_length <- gene_infos$Length_5kb_naive[gene_infos$ensembl_gene_id == substr(g_id,1,18)]
      }else if (condition == "lpa"){
        data_all <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_LPA_exons_rpkm_5kb,   g_id, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                          npRNA   = as.matrix(create_dat_g_rpkm(npRNA_LPA_exons_rpkm_5kb,   g_id, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                          cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_LPA_exons_rpkm_5kb, g_id, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
        gene_length <- gene_infos$Length_5kb_lpa[gene_infos$ensembl_gene_id == substr(g_id,1,18)]
      }
      time_data <- data_all$caRNA[!apply(data_all$caRNA, 1, function(x){all(is.na(x[-1]))}),1]
      
      # write experimental data to list
      all_data[[condition]] <- c(all_data[[condition]],list(data_all))
      names(all_data[[condition]])[length(all_data[[condition]])] <- gene_name
      
      # split data by cpt
      dat_g_exp_ca <- dat_g_exp_np <- dat_g_exp_cyto <- data.frame(RPKM=c(), Time = c(), batch = c(), Cpt = c(), Type=c(), Condition=c())
      for (bs in 2:ncol(data_all$caRNA)){
        dat_g_exp_ca <- rbind(dat_g_exp_ca, data.frame(RPKM=data_all$caRNA[,bs], Time = data_all$caRNA[,1], batch = c("rep1", "rep2", "rep3")[bs-1], Cpt = "Chromatin", Type="Experiment", Condition = condition))
        dat_g_exp_np <- rbind(dat_g_exp_np, data.frame(RPKM=data_all$npRNA[,bs], Time = data_all$npRNA[,1], batch = c("rep1", "rep2", "rep3")[bs-1], Cpt = "Nucleoplasm", Type="Experiment", Condition = condition))
        dat_g_exp_cyto <- rbind(dat_g_exp_cyto, data.frame(RPKM=data_all$cytoRNA[,bs], Time = data_all$cytoRNA[,1], batch = c("rep1", "rep2", "rep3")[bs-1], Cpt = "Cytoplasm", Type="Experiment", Condition = condition))
      }
      
      for (b in c('all', 'rep1', 'rep2', 'rep3')){   
        if(b == "all"){
          batch <- c("rep1", "rep2","rep3")
        }else{
          batch <- b
        }
        
        if (b == "rep3"){
          no_np <- TRUE
        }else{
          no_np <- FALSE
        }
        
        if (condition == "naive") {
          data <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_Naive_exons_rpkm_5kb,   g_id, replicate=batch, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                        npRNA   = as.matrix(create_dat_g_rpkm(npRNA_Naive_exons_rpkm_5kb,   g_id, replicate=batch, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                        cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_Naive_exons_rpkm_5kb, g_id, replicate=batch, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
        }else if (condition == "lpa"){
          data <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_LPA_exons_rpkm_5kb,   g_id, replicate=batch, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                        npRNA   = as.matrix(create_dat_g_rpkm(npRNA_LPA_exons_rpkm_5kb,   g_id, replicate=batch, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                        cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_LPA_exons_rpkm_5kb, g_id, replicate=batch, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
        }
        
        par <- bfgs[[b]]$parameter[which.min( bfgs[[b]]$value),]
        par_all <- bfgs[[b]]$parameter
        span <- par[1]
        sigmat <- par[2]
        sigmab <- par[3]
        
        # convert par and add param combination and write to list
        if(b %in% c('all', 'rep1') & condition == 'naive'){
          par_old <- c(par[1], abs(par[2:3]), par[4],parameters_convert_cp(par = par[-c(1:4)], log_in = T, log_out = T, direction = -1, parameter = 'row', no_np = no_np))
          par_old_all <- cbind(par_all[,1], abs(par_all[,2:3]), par_all[,4],parameters_convert_cp(par = par_all[,-c(1:4)], log_in = T, log_out = T, direction = -1, parameter = 'col', no_np = no_np))
          names(par_old) <- colnames(par_old_all) <- c('spar', 'sigmat', 'sigmab', 'ca0_rep1', "k1'", "k2", "k2'", "kdeg")
        }else{
          par_old <- c(par[1], abs(par[2:3]), parameters_convert_cp(par = par[-c(1:3)], log_in = T, log_out = T, direction = -1, parameter = 'row', no_np = no_np))
          par_old_all <- cbind(par_all[,1], abs(par_all[,2:3]), parameters_convert_cp(par = par_all[,-c(1:3)], log_in = T, log_out = T, direction = -1, parameter = 'col', no_np = no_np))
          if(no_np){
            names(par_old) <- colnames(par_old_all) <-c('spar', 'sigmat', 'sigmab', "k1'k2'", "max(k2,kdeg)", "min(k2,kdeg)")
          }else{
            names(par_old) <- colnames(par_old_all) <- c('spar', 'sigmat', 'sigmab', "k1'", "k2", "k2'", "kdeg")
          }
        }
        if(!no_np){
          par_old <- c(par_old, "k1'/k2" = unname(par_old["k1'"] - par_old["k2"]), "k1'k2'/k2" = unname(par_old["k1'"] + par_old["k2'"]- par_old["k2"]), "k2'/k2" = unname(par_old["k2'"]- par_old["k2"]))
          par_old_all <- cbind(par_old_all, "k1'/k2" = unname(par_old_all[,"k1'"] - par_old_all[,"k2"]), "k1'k2'/k2" = unname(par_old_all[,"k1'"] + par_old_all[,"k2'"]- par_old_all[,"k2"]), "k2'/k2" = unname(par_old_all[,"k2'"]- par_old_all[,"k2"]))
        }
        
        if(!(b %in% names(all_best_param[[condition]]))){
          all_best_param[[condition]][[b]] <- list()
          all_param[[condition]][[b]] <- list()
        }
        
        all_best_param[[condition]][[b]] <- rbind(all_best_param[[condition]][[b]], par_old)
        row.names(all_best_param[[condition]][[b]])[nrow(all_best_param[[condition]][[b]])] <- gene_name
        
        all_param[[condition]][[b]][[gene_name]] <- cbind(par_old_all,value=bfgs[[b]]$value)
        
        ## Extrapolate the chromatin fraction
        # times to extrapolate
        extra_t_end <- time_data[length(time_data)]
        extra_t_seq <- seq(time_data[1],extra_t_end + 1, by=1/60)
        
        # data to extrapolate
        data_extrapolation <- data$caRNA
        # replace missing t=0 by parameters 
        if(any(is.na(data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])]))){
          p <- par[-c(3+1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
          data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])] <- par[3+1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,]))]
        }else{
          p <- par
        }
        
        # remove time point if all na
        to_remove <- apply(data_extrapolation, 1, function(x){all(is.na(x[-1]))})
        data_extrapolation <- data_extrapolation[!to_remove,]
        
        dat_g_exp_pattern_ca <- paste0('chromatin.',sapply(nchar(as.character(dat_g_exp_ca$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp_ca$Time),'.', dat_g_exp_ca$batch)
        lib_indices_ca <- match(dat_g_exp_pattern_ca, row.names(lib_mat))
        libsize_ca <- as.integer(floor((lib_mat[,"lib.size"]*lib_mat[,"norm.factors"])[lib_indices_ca]))
        names(libsize_ca) <- dat_g_exp_pattern_ca
        
        dat_g_exp_pattern_np <- paste0('Nuc.',sapply(nchar(as.character(dat_g_exp_np$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp_np$Time),'.', dat_g_exp_np$batch)
        lib_indices_np <- match(dat_g_exp_pattern_np, row.names(lib_mat))
        libsize_np <- as.integer(floor((lib_mat[,"lib.size"]*lib_mat[,"norm.factors"])[lib_indices_np]))
        names(libsize_np) <- dat_g_exp_pattern_np
        
        dat_g_exp_pattern_cyto <- paste0('cytoplasmic.',sapply(nchar(as.character(dat_g_exp_cyto$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp_cyto$Time),'.', dat_g_exp_cyto$batch)
        lib_indices_cyto <- match(dat_g_exp_pattern_cyto, row.names(lib_mat))
        libsize_cyto <- as.integer(floor((lib_mat[,"lib.size"]*lib_mat[,"norm.factors"])[lib_indices_cyto]))
        names(libsize_cyto) <- dat_g_exp_pattern_cyto
        
        all_library[[condition]][[b]] <- list(caRNA = libsize_ca, npRNA = libsize_np, cytoRNA = libsize_cyto)
        
        SSEs_g_total <- TSE_reps_cp(par = par, fn_fit=fn_all_ode_cp, data = data_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=b, results = "sum", library_size = lib_mat, length = gene_length)
        SSEs_g_det <- TSE_reps_cp(par = par, fn_fit=fn_all_ode_cp, data = data_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=b, results = "detailed", library_size = lib_mat, length = gene_length)
        
        print(SSEs_g_total)
        print(min(bfgs[[b]]$value))
        
        par_theo <- par
        par_theo[2:3] <- c(0,1*10^-6)
        
        SSEs_g_total_theo <- TSE_reps_theo_cp(par = par_theo, fn_fit=fn_all_ode_cp, data = data_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=b, results = "sum", library_size = lib_mat, length = gene_length)
        
        # for each batch
        for (bs in 2:(ncol(data$caRNA))){
          proba_ca <- sapply(1:nrow(dat_g_exp_ca), function(i){ 2^dat_g_exp_ca$RPKM[i] * gene_length / 10^3 / 10^6 * dbeta( x = 2^dat_g_exp_ca$RPKM[i] * gene_length / 10^3 / 10^6, shape1 = 2^dat_g_exp_ca$RPKM[i] * gene_length/10^3 * libsize_ca[i]/10^6 + 1, shape2 = libsize_ca[i] - 2^dat_g_exp_ca$RPKM[i] * gene_length/10^3 * libsize_ca[i]/10^6 + 1)})
          proba_np <- sapply(1:nrow(dat_g_exp_np), function(i){ 2^dat_g_exp_np$RPKM[i] * gene_length / 10^3 / 10^6 * dbeta( x = 2^dat_g_exp_np$RPKM[i] * gene_length / 10^3 / 10^6, shape1 = 2^dat_g_exp_np$RPKM[i] * gene_length/10^3 * libsize_np[i]/10^6 + 1, shape2 = libsize_np[i] - 2^dat_g_exp_np$RPKM[i] * gene_length/10^3 * libsize_np[i]/10^6 + 1)})
          proba_cyto <- sapply(1:nrow(dat_g_exp_cyto), function(i){ 2^dat_g_exp_cyto$RPKM[i] * gene_length / 10^3 / 10^6 * dbeta( x = 2^dat_g_exp_cyto$RPKM[i] * gene_length / 10^3 / 10^6, shape1 = 2^dat_g_exp_cyto$RPKM[i] * gene_length/10^3 * libsize_cyto[i]/10^6 + 1, shape2 = libsize_cyto[i] - 2^dat_g_exp_cyto$RPKM[i] * gene_length/10^3 * libsize_cyto[i]/10^6 + 1)})
          
          if(ncol(data$caRNA) > 2){
            weight_ca <- proba_ca[dat_g_exp_ca$batch == paste0('rep',bs-1)]
            weight_np <- proba_np[dat_g_exp_np$batch == paste0('rep',bs-1)]
            weight_cyto <- proba_cyto[dat_g_exp_cyto$batch == paste0('rep',bs-1)]
          }else{
            weight_ca <- proba_ca[dat_g_exp_ca$batch == batch]
            weight_np <- proba_np[dat_g_exp_np$batch == batch]
            weight_cyto <- proba_cyto[dat_g_exp_cyto$batch == batch]
          }
          weight_ca <- weight_ca/max(weight_ca, na.rm=T)
          weight_np <- weight_np/max(weight_np, na.rm=T)
          weight_cyto <- weight_cyto/max(weight_cyto, na.rm=T)
          
          weight <- weight_ca
          # add max weigth to the missing t=0 values replaced by parameter (i.e. sure values)
          if(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))){
            weight[which(is.na(data$caRNA[data$caRNA[,1]==0,-1]))] <- 1
          }
          
          # extrapolate caRNA
          data_extra <- data.frame(x=as.vector(data_extrapolation[,1]), y=as.vector(data_extrapolation[,bs]))
          weight <- weight[!is.na(weight)]
          data_extra <- data_extra[!is.na(data_extra$y),]
          
          # order data to extrapolate and weight
          weight <- weight[order(data_extra$x)]
          data_extra <- data_extra[order(data_extra$x),]
          
          extrapolated_input <-  cbind(extra_t_seq, 2^predict(smooth.spline(x=data_extra, w = weight , spar = span), extra_t_seq)$y)
          
          # simulate data           
          fit_all <- fn_all_ode_cp(t = seq(min(time_data),max(time_data)), p[-c(1:3)], extrapolated_input = extrapolated_input, convert_param=T, log_in=T,u=input_fun_cp, no_np = no_np)
          fit <- fit_all
          if(any(colnames(fit_all) == "dcytoRNA")){
            fit <- fit_all[,-which(colnames(fit_all) == "dcytoRNA")]
          }
          fit[,-1] <- log2(fit[,-1])
          if (nrow(fit) != length(seq(min(time_data),max(time_data)))){
            next
          }
          fit <- cbind(fit, caRNA=log2(extrapolated_input[(seq(min(time_data),max(time_data)) - time_data[1])*60+1,2]))
          if(no_np){
            fit <- cbind(fit, npRNA=NA)
          }
          
          # check if exist because all is done one after each other
          if(!(b %in% names(all_fit[[condition]]))){
            all_fit[[condition]][[b]] <- list()
          }
          if(!(gene_name %in% names(all_fit[[condition]][[b]]))){
            all_fit[[condition]][[b]][[gene_name]] <- list(caRNA=fit[,c("time","caRNA")], npRNA=fit[,c("time","npRNA")], cytoRNA=fit[,c("time","cytoRNA")])
          }else{
            all_fit[[condition]][[b]][[gene_name]][["caRNA"]] <- cbind(all_fit[[condition]][[b]][[gene_name]][["caRNA"]], fit[,"caRNA"])
            all_fit[[condition]][[b]][[gene_name]][["npRNA"]] <- cbind(all_fit[[condition]][[b]][[gene_name]][["npRNA"]], fit[,"npRNA"])
            all_fit[[condition]][[b]][[gene_name]][["cytoRNA"]] <- cbind(all_fit[[condition]][[b]][[gene_name]][["cytoRNA"]], fit[,"cytoRNA"])
          }
          
          # calculate dmu (diff for caRNA, and analytic for np and cyto)
          dmu_ca <- data.frame(caRNA=c((extrapolated_input[2,2]-extrapolated_input[1,2])/(extrapolated_input[2,1]-extrapolated_input[1,1]), 
                                       0.5*((extrapolated_input[,2][-c(nrow(extrapolated_input)-1:0)]-extrapolated_input[,2][-c(1,nrow(extrapolated_input))])/(extrapolated_input[,1][-c(nrow(extrapolated_input)-1:0)]-extrapolated_input[,1][-c(1,nrow(extrapolated_input))])+
                                              (extrapolated_input[,2][-c(1:2)]-extrapolated_input[,2][-c(1,nrow(extrapolated_input))])/(extrapolated_input[,1][-c(1:2)]-extrapolated_input[,1][-c(1,nrow(extrapolated_input))])),
                                       (extrapolated_input[nrow(extrapolated_input),2]-extrapolated_input[nrow(extrapolated_input)-1,2])/(extrapolated_input[nrow(extrapolated_input),1]-extrapolated_input[nrow(extrapolated_input)-1,1])))
          dmu <- data.frame(caRNA = dmu_ca[(seq(min(time_data),max(time_data)) - time_data[1])*60+1,]/(2^fit[,"caRNA"]))
          
          for (col in colnames(fit)){
            if (col == "npRNA"){
              dmu_col <- log(2)*data.frame(npRNA=c((fit[2,col]-fit[1,col])/(fit[2,1]-fit[1,1]), 
                                                   0.5*((fit[,col][-c(nrow(fit)-1:0)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(nrow(fit)-1:0)]-fit[,1][-c(1,nrow(fit))])+
                                                          (fit[,col][-c(1:2)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(1:2)]-fit[,1][-c(1,nrow(fit))])),
                                                   (fit[nrow(fit),col]-fit[nrow(fit)-1,col])/(fit[nrow(fit),1]-fit[nrow(fit)-1,1])))
              dmu <- cbind(dmu, npRNA = dmu_col)
            }
            if (col == "cytoRNA"){
              dmu_col <- log(2)*data.frame(cytoRNA=c((fit[2,col]-fit[1,col])/(fit[2,1]-fit[1,1]), 
                                                     0.5*((fit[,col][-c(nrow(fit)-1:0)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(nrow(fit)-1:0)]-fit[,1][-c(1,nrow(fit))])+
                                                            (fit[,col][-c(1:2)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(1:2)]-fit[,1][-c(1,nrow(fit))])),
                                                     (fit[nrow(fit),col]-fit[nrow(fit)-1,col])/(fit[nrow(fit),1]-fit[nrow(fit)-1,1])))
              dmu <- cbind(dmu, cytoRNA = dmu_col)
            }
          }
          
          libsize <- data.frame(caRNA = libsize_ca[dat_g_exp_ca$batch == ifelse(length(batch)>1, colnames(data$caRNA)[bs],b)],
                                npRNA = if(!all(is.na(libsize_np[dat_g_exp_np$batch == ifelse(length(batch)>1, colnames(data$npRNA)[bs],b)]))){libsize_np[dat_g_exp_np$batch == ifelse(length(batch)>1, colnames(data$npRNA)[bs],b)]}else{NA},
                                cytoRNA = libsize_cyto[dat_g_exp_cyto$batch == ifelse(length(batch)>1, colnames(data$cytoRNA)[bs],b)])
          
          if(no_np){
            libsize <- libsize[, -which(colnames(libsize) == "npRNA")]
          }                      
          sigma <- sqrt(dmu^2*sigmat^2 + sigmab^2)[match(sort(unique(dat_g_exp_ca$Time)), fit[,"time"]),]
          m <- 1/(exp(sigma^2) - 1)
          m[!is.na(sigma) & exp(sigma^2) == 1] <- (1/sigma^2)[!is.na(sigma) & exp(sigma^2) == 1] # approx if sigma close to zero
          mu <- exp(log(2^fit[match(dat_g_exp_ca$Time[dat_g_exp_ca$batch == ifelse(length(batch)>1,colnames(data$caRNA)[bs],b)],fit[,1]),match(colnames(libsize), colnames(fit))]) + log(libsize/10^6*gene_length/10^3) + 1/2*sigma[,match(colnames(libsize), colnames(sigma))]^2)
          
          ub_ca <- log2(sapply(1:nrow(libsize), function(i){qnbinom(p = 0.025, size = m[i,"caRNA"], mu = mu[i, "caRNA"], lower.tail=F)})/(libsize[,"caRNA"]/10^6*gene_length/10^3))
          lb_ca <- log2(sapply(1:nrow(libsize), function(i){qnbinom(p = 0.025, size = m[i,"caRNA"], mu = mu[i, "caRNA"], lower.tail=T)})/(libsize[,"caRNA"]/10^6*gene_length/10^3))
          
          if(!no_np){
            ub_np <- log2(sapply(1:nrow(libsize), function(i){qnbinom(p = 0.025, size = m[i,"npRNA"], mu = mu[i, "npRNA"], lower.tail=F)})/(libsize[,"npRNA"]/10^6*gene_length/10^3))
            lb_np <- log2(sapply(1:nrow(libsize), function(i){qnbinom(p = 0.025, size = m[i,"npRNA"], mu = mu[i, "npRNA"], lower.tail=T)})/(libsize[,"npRNA"]/10^6*gene_length/10^3))
          }else{
            lb_np <- NA
            ub_np <- NA
          }
          
          ub_cyto <- log2(sapply(1:nrow(libsize), function(i){qnbinom(p = 0.025, size = m[i,"cytoRNA"], mu = mu[i, "cytoRNA"], lower.tail=F)})/(libsize[,"cytoRNA"]/10^6*gene_length/10^3))
          lb_cyto <- log2(sapply(1:nrow(libsize), function(i){qnbinom(p = 0.025, size = m[i,"cytoRNA"], mu = mu[i, "cytoRNA"], lower.tail=T)})/(libsize[,"cytoRNA"]/10^6*gene_length/10^3))
          
          if(!(b %in% names(all_fit_CI[["LB"]][[condition]]))){
            all_fit_CI[["LB"]][[condition]][[b]] <- list()
            all_fit_CI[["UB"]][[condition]][[b]] <- list()
          }
          if(!(gene_name %in% names(all_fit_CI[["LB"]][[condition]][[b]]))){
            all_fit_CI[["LB"]][[condition]][[b]][[gene_name]] <- list(caRNA = cbind(data$caRNA[,"t"], lb_ca),
                                                                      npRNA = cbind(data$npRNA[,"t"], lb_np),
                                                                      cytoRNA = cbind(data$cytoRNA[,"t"], lb_cyto))
            all_fit_CI[["UB"]][[condition]][[b]][[gene_name]] <- list(caRNA = cbind(data$caRNA[,"t"], ub_ca),
                                                                      npRNA= cbind(data$npRNA[,"t"], ub_np), 
                                                                      cytoRNA = cbind(data$cytoRNA[,"t"], ub_cyto))
          }else{
            all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["caRNA"]] <- cbind(all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["caRNA"]], lb_ca)
            all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["npRNA"]] <- cbind(all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["npRNA"]], lb_np)
            all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["cytoRNA"]] <- cbind(all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["cytoRNA"]], lb_cyto)
            all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["caRNA"]] <- cbind(all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["caRNA"]], ub_ca)
            all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["npRNA"]] <- cbind(all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["npRNA"]], ub_np)
            all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["cytoRNA"]] <- cbind(all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["cytoRNA"]], ub_cyto)
          }
          
          autocor_error_ca <- cor(c((fit[fit[,1] %in% data$caRNA[,1],'caRNA']-data$caRNA[,-1])[!is.na(data$caRNA[,-1])],NA),c(NA,(fit[fit[,1] %in% data$caRNA[,1],'caRNA']-data$caRNA[,-1])[!is.na(data$caRNA[,-1])]),use="pairwise.complete.obs")
          autocor_error_np <- cor(c((fit[fit[,1] %in% data$npRNA[,1],'npRNA']-data$npRNA[,-1])[!is.na(data$npRNA[,-1])],NA),c(NA,(fit[fit[,1] %in% data$npRNA[,1],'npRNA']-data$npRNA[,-1])[!is.na(data$npRNA[,-1])]),use="pairwise.complete.obs")
          autocor_error_cyto <- cor(c((fit[fit[,1] %in% data$cytoRNA[,1],'cytoRNA']-data$cytoRNA[,-1])[!is.na(data$cytoRNA[,-1])],NA),c(NA,(fit[fit[,1] %in% data$cytoRNA[,1],'cytoRNA']-data$cytoRNA[,-1])[!is.na(data$cytoRNA[,-1])]),use="pairwise.complete.obs")
          
          mean_error_over_range_ca <- mean(abs((fit[fit[,1] %in% data$caRNA[,1],'caRNA']-data$caRNA[,-1])),na.rm=T)/(max(data$caRNA[,-1], na.rm=T)-min(data$caRNA[,-1], na.rm=T))
          mean_error_over_range_np <- mean(abs((fit[fit[,1] %in% data$npRNA[,1],'npRNA']-data$npRNA[,-1])),na.rm=T)/(max(data$npRNA[,-1], na.rm=T)-min(data$npRNA[,-1], na.rm=T))
          mean_error_over_range_cyto <- mean(abs((fit[fit[,1] %in% data$cytoRNA[,1],'cytoRNA']-data$cytoRNA[,-1])),na.rm=T)/(max(data$cytoRNA[,-1], na.rm=T)-min(data$cytoRNA[,-1],na.rm=T))
          
          autocor_error_ca_times_mean_error_over_range <- autocor_error_ca*mean_error_over_range_ca
          autocor_error_np_times_mean_error_over_range <- autocor_error_np*mean_error_over_range_np
          autocor_error_cyto_times_mean_error_over_range <- autocor_error_cyto*mean_error_over_range_cyto
          
          wautocor_error_ca <- wcor(c((fit[fit[,1] %in% data$caRNA[,1],'caRNA']-data$caRNA[,-1])[!is.na(data$caRNA[,-1])],NA),c(NA,(fit[fit[,1] %in% data$caRNA[,1],'caRNA']-data$caRNA[,-1])[!is.na(data$caRNA[,-1])]),c(1,(weight_ca[!is.na(data$caRNA[,-1])][-1]+weight_ca[!is.na(data$caRNA[,-1])][-length(weight_ca[!is.na(data$caRNA[,-1])])])/2,1))
          wautocor_error_np <- wcor(c((fit[fit[,1] %in% data$npRNA[,1],'npRNA']-data$npRNA[,-1])[!is.na(data$npRNA[,-1])],NA),c(NA,(fit[fit[,1] %in% data$npRNA[,1],'npRNA']-data$npRNA[,-1])[!is.na(data$npRNA[,-1])]),c(1,(weight_np[!is.na(data$npRNA[,-1])][-length(weight_np[!is.na(data$npRNA[,-1])])]+weight_np[!is.na(data$npRNA[,-1])][-1])/2,1))
          wautocor_error_cyto <- wcor(c((fit[fit[,1] %in% data$cytoRNA[,1],'cytoRNA']-data$cytoRNA[,-1])[!is.na(data$cytoRNA[,-1])],NA),c(NA,(fit[fit[,1] %in% data$cytoRNA[,1],'cytoRNA']-data$cytoRNA[,-1])[!is.na(data$cytoRNA[,-1])]),c(1,(weight_cyto[!is.na(data$cytoRNA[,-1])][-length(weight_cyto[!is.na(data$cytoRNA[,-1])])]+weight_cyto[!is.na(data$cytoRNA[,-1])][-1])/2,1))
          
          wmean_error_ca_over_range <- sum(weight_ca * abs((fit[fit[,1] %in% data$caRNA[,1],'caRNA']-data$caRNA[,-1])),na.rm=T) /sum(weight_ca, na.rm=T)/(max(data$caRNA[,-1], na.rm=T)-min(data$caRNA[,-1], na.rm=T))
          wmean_error_np_over_range <- sum(weight_np * abs((fit[fit[,1] %in% data$npRNA[,1],'npRNA']-data$npRNA[,-1])),na.rm=T)/sum(weight_np, na.rm=T)/(max(data$npRNA[,-1], na.rm=T)-min(data$npRNA[,-1], na.rm=T))
          wmean_error_cyto_over_range <- sum(weight_cyto * abs((fit[fit[,1] %in% data$cytoRNA[,1],'cytoRNA']-data$cytoRNA[,-1])),na.rm=T)/sum(weight_cyto, na.rm=T)/(max(data$cytoRNA[,-1], na.rm=T)-min(data$cytoRNA[,-1],na.rm=T))
          
          wautocor_error_ca_times_wmean_error_over_range <- wautocor_error_ca * wmean_error_ca_over_range
          wautocor_error_np_times_wmean_error_over_range <- wautocor_error_np * wmean_error_np_over_range
          wautocor_error_cyto_times_wmean_error_over_range <- wautocor_error_cyto * wmean_error_cyto_over_range
          
          tmp <- c(max(abs(c(autocor_error_ca_times_mean_error_over_range,autocor_error_np_times_mean_error_over_range, autocor_error_cyto_times_mean_error_over_range)), na.rm=T), 
                   max((1+abs(c(autocor_error_ca,autocor_error_np, autocor_error_cyto)))*c(mean_error_over_range_ca,mean_error_over_range_np,mean_error_over_range_cyto), na.rm=T),
                   max(abs(c(wautocor_error_ca_times_wmean_error_over_range,wautocor_error_np_times_wmean_error_over_range, wautocor_error_cyto_times_wmean_error_over_range)), na.rm=T),
                   max((1+abs(c(wautocor_error_ca,wautocor_error_np, wautocor_error_cyto)))*c(wmean_error_ca_over_range,wmean_error_np_over_range,wmean_error_cyto_over_range), na.rm=T),
                   max(abs(c(wmean_error_ca_over_range, wmean_error_np_over_range, wmean_error_cyto_over_range)), na.rm=T), 
                   SSEs_g_det[[bs-1]])
          
          if(!(b %in% names(all_annot[[condition]]))){
            all_annot[[condition]][[b]] <- list()
          }
          
          if(!(gene_name %in% names(all_annot[[condition]][[b]]))){
            all_annot[[condition]][[b]][[gene_name]] <- list(length=gene_length, annot_bs = data.frame(tmp))
          }else{
            all_annot[[condition]][[b]][[gene_name]][["annot_bs"]] <- cbind(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]], tmp)
          }
        }
        
        if(!is.null(all_fit[[condition]][[b]][[gene_name]])){
          colnames(all_fit[[condition]][[b]][[gene_name]][["caRNA"]]) <- colnames(data$caRNA)
          colnames(all_fit[[condition]][[b]][[gene_name]][["npRNA"]]) <- colnames(data$npRNA)
          colnames(all_fit[[condition]][[b]][[gene_name]][["cytoRNA"]]) <- colnames(data$cytoRNA)
          
          colnames(all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["caRNA"]]) <- colnames(data$caRNA)
          colnames(all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["npRNA"]]) <- colnames(data$npRNA)
          colnames(all_fit_CI[["LB"]][[condition]][[b]][[gene_name]][["cytoRNA"]]) <- colnames(data$cytoRNA)
          
          colnames(all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["caRNA"]]) <- colnames(data$caRNA)
          colnames(all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["npRNA"]]) <- colnames(data$npRNA)
          colnames(all_fit_CI[["UB"]][[condition]][[b]][[gene_name]][["cytoRNA"]]) <- colnames(data$cytoRNA)
          
          if(is.null(ncol(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]]))){
            all_annot[[condition]][[b]][[gene_name]][["annot_bs"]] <- data.frame(all_annot[[condition]][[b]][[gene_name]][["annot"]])
          }
          colnames(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]]) <- colnames(data$caRNA)[-1]
          
          
          row.names(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]]) <- c("max(|autocor(error)|*mean(|error|)/range)", 
                                                                                 "max((1+|autocor(error)|)*mean(|error|)/range)", 
                                                                                 "max(|w_autocor(error)|*w_mean(|error|)/range)", 
                                                                                 "max((1+|w_autocor(error)|)*w_mean(|error|)/range)", 
                                                                                 "max(w_mean(|error|)/range)",
                                                                                 "NLL_ca", "NLL_np", "NLL_cyto") 
          
          all_annot[[condition]][[b]][[gene_name]][["annot"]] <- data.frame(total=c(as.vector(apply(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]],1, mean, na.rm=T))[-c(6:8)], SSEs_g_total, apply(all_annot[[condition]][[b]][[gene_name]]$annot_bs,1,sum, na.rm=T)[6:8], SSEs_g_total_theo, SSEs_g_total - SSEs_g_total_theo))
          
          row.names(all_annot[[condition]][[b]][[gene_name]][["annot"]]) <- c(row.names(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]])[-grep("NLL",row.names(all_annot[[condition]][[b]][[gene_name]][["annot_bs"]]))],"NLL", "NLL_ca", "NLL_np", "NLL_cyto", "NLL_theo", "NNL-NLL_theo")
        }
      }
    }
  }
}
stop()

# Get Profile and CI data
# => save data
merge_function_for_rep3 <- function(dat_g, th){
  dat_g[["k1'k2'/k2"]] <- NULL
  if(all(c("kdeg","k2") %in% names(dat_g))){
    if( ((max(dat_g$kdeg$Profile$kdeg) > min(dat_g$k2$Profile$k2)) & (min(dat_g$kdeg$Profile$kdeg) < min(dat_g$k2$Profile$k2))) || ((max(dat_g$k2$Profile$k2) > min(dat_g$kdeg$Profile$kdeg)) & (min(dat_g$k2$Profile$k2) < min(dat_g$kdeg$Profile$kdeg)))){
      tmp <- rbind(dat_g$kdeg$Profile, dat_g$k2$Profile)
      tmp$type <- replace(rep("max", nrow(tmp)), sapply(1:nrow(tmp),function(i){colnames(tmp[,c("k2","kdeg")])[which.min(tmp[i,c("k2","kdeg")])] == tmp$whichPar[i]}), "min")
      
      tmp_min <- tmp[tmp$type == "min",]
      tmp_min$kval <- apply(tmp_min[,c("k2","kdeg")],1,min)
      tmp_max <- tmp[tmp$type == "max",]
      tmp_max$kval <- apply(tmp_max[,c("k2","kdeg")],1,max)
      
      tmp_min <- tmp_min[order(tmp_min$kval),]
      tmp_max <- tmp_max[order(tmp_max$kval),]
      
      not_finished <- TRUE
      while(not_finished){
        to_remove <- c()
        for(i in 2:(nrow(tmp_min)-1)){
          if(tmp_min[i,"whichPar"] != tmp_min[i-1,"whichPar"]){
            if (tmp_min[i,"value"] > tmp_min[i-1,"value"] & tmp_min[i,"value"] > tmp_min[which.max(tmp_min$whichPar == unlist(tmp_min[i-1,"whichPar"]) & tmp_min$kval > unlist(tmp_min[i,"kval"])),"value"]){
              to_remove <- c(to_remove,i)
            }
          }
        }
        if(length(to_remove) == 0){
          not_finished <- FALSE
        }else{
          tmp_min <- tmp_min[-to_remove,]
        }
      }
      
      not_finished <- TRUE
      while(not_finished){
        to_remove <- c()
        for(i in 2:(nrow(tmp_max)-1)){
          if(tmp_max[i,"whichPar"] != tmp_max[i-1,"whichPar"]){
            if (tmp_max[i,"value"] > tmp_max[i-1,"value"] & tmp_max[i,"value"] > tmp_max[which.max(tmp_max$whichPar == unlist(tmp_max[i-1,"whichPar"]) & tmp_max$kval > unlist(tmp_max[i,"kval"])),"value"]){
              to_remove <- c(to_remove,i)
            }
          }
        }
        if(length(to_remove) == 0){
          not_finished <- FALSE
        }else{
          tmp_max <- tmp_max[-to_remove,]
        }
      }
      tmp <- rbind(tmp_min,tmp_max)
      
      i <- 2 
      to_remove <- c()
      while (i < nrow(tmp)){
        if(tmp$whichPar[i] != tmp$whichPar[i-1] && tmp$whichPar[i] != tmp$whichPar[i+1] && !((i-1) %in% to_remove)){
          to_remove <- c(to_remove,i)
        }
        i <- i+1
      }
      
      if(!is.null(to_remove)){
        tmp <- tmp[-to_remove,]
      }
      
      # check that it goes above the threshold
      if(tmp$value[1] >= th){
        if(tmp$value[nrow(tmp)] >= th){
          tmp <- tmp[(which(tmp$value <= th)[1]-1):(which(tmp$value <= th)[length(which(tmp$value <= th))]+1),]
        }else{
          tmp <- tmp[(which(tmp$value <= th)[1]-1):(which(tmp$value <= th)[length(which(tmp$value <= th))]),]
        }
      }else{
        if(tmp$value[nrow(tmp)] >= th){
          tmp <- tmp[(which(tmp$value <= th)[1]):(which(tmp$value <= th)[length(which(tmp$value <= th))]+1),]
        }else{
          tmp <- tmp[(which(tmp$value <= th)[1]):(which(tmp$value <= th)[length(which(tmp$value <= th))]),]
        }
      }
      
      tmp_kdeg <- tmp_k2 <- tmp
      
      tmp_kdeg[tmp$whichPar == "k2","k2"] <- tmp[tmp$whichPar == "k2","kdeg"]
      tmp_kdeg[tmp$whichPar == "k2","kdeg"] <- tmp[tmp$whichPar == "k2","kval"]
      tmp_kdeg$whichPar <- "kdeg"
      tmp_kdeg$constraint <- tmp_kdeg$kdeg - dat_g$kdeg$Profile$kdeg[dat_g$kdeg$Profile$constraint == 0]
      
      tmp_k2[tmp$whichPar == "kdeg","kdeg"] <- tmp[tmp$whichPar == "kdeg","k2"]
      tmp_k2[tmp$whichPar == "kdeg","k2"] <- tmp[tmp$whichPar == "kdeg","kval"]
      tmp_k2$whichPar <- "k2"
      tmp_k2$constraint <- tmp_k2$k2 - dat_g$k2$Profile$k2[dat_g$k2$Profile$constraint == 0]
      
      if(tmp$value[1] >= th){
        dat_g$kdeg$CI$lower <- dat_g$k2$CI$lower <- tmp_kdeg$kdeg[2] + (th-tmp_kdeg$value[2])/(tmp_kdeg$value[1]-tmp_kdeg$value[2])*(tmp_kdeg$kdeg[1]-tmp_min$kdeg[2])
      }else if( is.finite(dat_g$k2$CI$lower) || is.finite(dat_g$kdeg$CI$lower)){
        lower <- min(c(dat_g$k2$CI$lower,dat_g$kdeg$CI$lower)[is.finite(c(dat_g$k2$CI$lower,dat_g$kdeg$CI$lower))])
        if(lower <= min(tmp$kval)){
          dat_g$kdeg$CI$lower <- dat_g$k2$CI$lower <- lower
        }else{
          dat_g$kdeg$CI$lower <- dat_g$k2$CI$lower <- -Inf
        }
      }
      
      if(tmp$value[nrow(tmp)] >= th){
        dat_g$kdeg$CI$upper <- dat_g$k2$CI$upper <- tmp_kdeg$kdeg[nrow(tmp_kdeg)-1] + (th-tmp_kdeg$value[nrow(tmp_kdeg)-1])/(tmp_kdeg$value[nrow(tmp_kdeg)]-tmp_kdeg$value[nrow(tmp_kdeg)-1])*(tmp_kdeg$kdeg[nrow(tmp_kdeg)]-tmp_kdeg$kdeg[nrow(tmp_kdeg)-1])
      }else if( is.finite(dat_g$k2$CI$upper) || is.finite(dat_g$kdeg$CI$upper)){
        upper <- max(c(dat_g$k2$CI$upper,dat_g$kdeg$CI$upper)[is.finite(c(dat_g$k2$CI$upper,dat_g$kdeg$CI$upper))])
        if(upper >= max(tmp$kval)){
          dat_g$kdeg$CI$upper <- dat_g$k2$CI$upper <- upper
        }else{
          dat_g$kdeg$CI$upper <- dat_g$k2$CI$upper <- Inf
        }
      }
      
      plot(tmp$kval, tmp$value, pch=20, cex=0.5)
      points(dat_g$k2$Profile$k2, dat_g$k2$Profile$value, col="blue")
      points(dat_g$kdeg$Profile$kdeg, dat_g$kdeg$Profile$value, col="red")
      
      dat_g$k2$Profile <- tmp_k2[,!colnames(tmp_k2) %in% c("type", "kval")]
      dat_g$kdeg$Profile <- tmp_kdeg[,!colnames(tmp_kdeg) %in% c("type", "kval")]
    }else{
      plot(dat_g$kdeg$Profile$kdeg,dat_g$kdeg$Profile$value , col="red", pch=20, cex=0.5, xlim=range(c(dat_g$kdeg$Profile$kdeg, dat_g$k2$Profile$k2)))
      points(dat_g$k2$Profile$k2, dat_g$k2$Profile$value , col="blue", pch=20, cex=0.5)
    }
    
    names(dat_g)[names(dat_g)=="kdeg"] <- "min(k2,kdeg)"
    names(dat_g)[names(dat_g)=="k2"] <- "max(k2,kdeg)"
    
    row.names(dat_g[["min(k2,kdeg)"]][["CI"]]) <- dat_g[["min(k2,kdeg)"]][["CI"]][,"name"] <- as.character("min(k2,kdeg)")
    row.names(dat_g[["max(k2,kdeg)"]][["CI"]]) <- dat_g[["max(k2,kdeg)"]][["CI"]][,"name"] <- as.character("max(k2,kdeg)")
  }  
  return(dat_g)
}

all_param_CI <- list()

for (condition in c('naive', 'lpa')){
  all_param_CI[[condition]] <- list()
  if(condition == 'naive'){
    raw_data <- caRNA_Naive_exons_rpkm_5kb
  }else if (condition == 'lpa'){
    raw_data <- caRNA_LPA_exons_rpkm_5kb
  }
  
  for (g in gene_infos$ensembl_gene_id[order(gene_infos$external_gene_name)]){
    if (g %in% sapply(row.names(raw_data), function(s){substr(s,1,18)})){
      i <- unname(which(sapply(row.names(raw_data), function(s){substr(s,1,18)}) == g ))
      g_id <- rownames(raw_data)[i]
      
      gene_name <- unique(gene_infos$external_gene_name[gene_infos$ensembl_gene_id == substr(g_id,1,18)])
      
      print(paste0(g, ' - ', gene_name))
      # print(gene_name)
      
      for ( batch in c('all','rep1','rep2','rep3')){
        dat_CI_tmp <- c()
        dat_profile_tmp <- list()
        dat_profile_tmp_all <- list()
        
        if(!is.null(all_param[[condition]][[batch]][[gene_name]])){
          best <- which.min(all_param[[condition]][[batch]][[gene_name]][,"value"])
          min_val <- min(all_param[[condition]][[batch]][[gene_name]][,"value"])
          below_CI <- all_param[[condition]][[batch]][[gene_name]][,"value"] <= min_val + qchisq(0.99,1,lower.tail = T)/2
          param_best <- all_param[[condition]][[batch]][[gene_name]][best, -which(colnames(all_param[[condition]][[batch]][[gene_name]])=="value")]
          
          if(batch == 'rep3'){
            no_np <- T
          }else{
            no_np <- F
          }
          
          consistency_check <- consistency_check_approx <- consistency_check_exact <- c()
          consistency_check_lower <- consistency_check_lower_approx <- consistency_check_lower_exact <- c()
          consistency_check_upper <- consistency_check_upper_approx <- consistency_check_upper_exact <- c()
          
          # keep results to avoid reloading
          
          file_approx <- paste0('../Modeling/Results/Profile_dMod_integrate/', condition, '/profile_dMod_', batch, '_', gene_name,'_Finished.Rdata')
          file_approx_partial <- paste0('../Modeling/Results/Profile_dMod_integrate/', condition, '/sa_bfgs_', batch, '_', gene_name,'.Rdata')
          if(file.exists(file_approx) | file.exists(file_approx_partial)){
            if(file.exists(file_approx)){
              load(file_approx)
            }else{
              load(file_approx_partial)
            }
            
            if(batch == 'rep3'){
              print("approx")
              dat_g <- merge_function_for_rep3(dat_g, th = 2*(min_val + qchisq(0.95,1,lower.tail = T)/2))
            }
            
            for(i in 1:length(dat_g)){
              dat_g[[i]]$computing_time <- NULL
            }
            
            dat_profile_tmp_all[["approx"]] <- dat_g 
            
            for (par in colnames(all_param[[condition]][[batch]][[gene_name]])){
              if(par %in% c("value", "spar", "sigmat", "sigmab", "ca0_rep1")){
                next
              }
              if (par %in% names(dat_g)){
                l <- dat_g[[par]]
                if (batch == 'rep3' && l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)") && max(dat_g$"min(k2,kdeg)"$Profile$kdeg) >= min(dat_g$"max(k2,kdeg)"$Profile$k2)){ 
                  check_lower <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][, "max(k2,kdeg)"] < ifelse(is.na(l$CI$lower),-Inf,l$CI$lower) | all_param[[condition]][[batch]][[gene_name]][, "min(k2,kdeg)"] < ifelse(is.na(l$CI$lower),-Inf,l$CI$lower)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                  check_upper <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][, "max(k2,kdeg)"] > ifelse(is.na(l$CI$upper),Inf,l$CI$upper) | all_param[[condition]][[batch]][[gene_name]][, "min(k2,kdeg)"] > ifelse(is.na(l$CI$upper),Inf,l$CI$upper)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                }else{
                  check_lower <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][,as.character(l$CI$name)] < ifelse(is.na(l$CI$lower),-Inf,l$CI$lower)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                  check_upper <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][,as.character(l$CI$name)] > ifelse(is.na(l$CI$upper), Inf,l$CI$upper)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                }
                
                if(is.na(l$CI$lower) | is.infinite(l$CI$lower)){
                  check_lower <- FALSE
                }
                if(is.na(l$CI$upper) | is.infinite(l$CI$upper)){
                  check_upper <- FALSE
                }
                
                consistency_check_lower_approx <- c(consistency_check_lower_approx,check_lower)
                consistency_check_upper_approx <- c(consistency_check_upper_approx,check_upper)
                
                if(check_lower & check_upper){
                  consistency_check_approx <- c(consistency_check_approx,T)
                  dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower='approximate', method_upper="approximate"))
                  if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                    dat_profile_tmp[[as.character(l$CI$name)]] <- cbind(dat_g[[as.character(l$CI$name)]]$Profile[,c(1:5, which(colnames(dat_g[[as.character(l$CI$name)]]$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "approximate")
                  }else{
                    dat_profile_tmp[[as.character(l$CI$name)]] <- cbind(dat_g[[as.character(l$CI$name)]]$Profile[,c(1:5, which(colnames(dat_g[[as.character(l$CI$name)]]$Profile)==as.character(l$CI$name)))], method = "approximate")
                  }
                }else{
                  consistency_check_approx <- c(consistency_check_approx, F)
                }
              }else{
                consistency_check_lower_approx <- c(consistency_check_lower_approx,F)
                consistency_check_upper_approx <- c(consistency_check_upper_approx,F)
                consistency_check_approx <- c(consistency_check_approx, F)
              }
              names(consistency_check_approx)[length(consistency_check_approx)] <- par
              names(consistency_check_lower_approx)[length(consistency_check_lower_approx)] <- par
              names(consistency_check_upper_approx)[length(consistency_check_upper_approx)] <- par
            }
            consistency_check <- consistency_check_approx
            consistency_check_lower <- consistency_check_lower_approx
            consistency_check_upper <- consistency_check_upper_approx
          }  
          
          if( length(consistency_check) == 0 || !all(consistency_check) || any(is.na(consistency_check))){
            file_exact <- paste0('../Modeling/Results/Profile_dMod_optimize/', condition, '/profile_dMod_', batch, '_', gene_name,'_Finished.Rdata')
            file_exact_partial <- paste0('../Modeling/Results/Profile_dMod_optimize/', condition, '/sa_bfgs_', batch, '_', gene_name,'.Rdata')
            if(file.exists(file_exact) | file.exists(file_exact_partial)){
              if(file.exists(file_exact)){
                load(file_exact)
              }else{
                load(file_exact_partial)
              }
              
              if(batch == 'rep3'){
                print("exact")
                dat_g <- merge_function_for_rep3(dat_g, th = 2*(min_val + qchisq(0.95,1,lower.tail = T)/2))
              }
              
              for(i in 1:length(dat_g)){
                dat_g[[i]]$computing_time <- NULL
              }
              
              dat_profile_tmp_all[["exact"]] <- dat_g 
              
              for (par in colnames(all_param[[condition]][[batch]][[gene_name]])){
                if(par %in% c("value", "spar", "sigmat", "sigmab", "ca0_rep1")){
                  next
                }
                if (par %in% names(dat_g)){
                  l <- dat_g[[par]]
                  if (batch == 'rep3' && l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)") && max(dat_g$"min(k2,kdeg)"$Profile$kdeg) >= min(dat_g$"max(k2,kdeg)"$Profile$k2)){ 
                    check_lower <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][, "max(k2,kdeg)"] < ifelse(is.na(l$CI$lower),-Inf,l$CI$lower) | all_param[[condition]][[batch]][[gene_name]][, "min(k2,kdeg)"] < ifelse(is.na(l$CI$lower),-Inf,l$CI$lower)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                    check_upper <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][, "max(k2,kdeg)"] > ifelse(is.na(l$CI$upper),Inf,l$CI$upper) | all_param[[condition]][[batch]][[gene_name]][, "min(k2,kdeg)"] > ifelse(is.na(l$CI$upper),Inf,l$CI$upper)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                  }else{
                    check_lower <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][,as.character(l$CI$name)] < ifelse(is.na(l$CI$lower),-Inf,l$CI$lower)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                    check_upper <- all(all_param[[condition]][[batch]][[gene_name]][,"value"][all_param[[condition]][[batch]][[gene_name]][,as.character(l$CI$name)] > ifelse(is.na(l$CI$upper), Inf,l$CI$upper)] > min_val + qchisq(0.95,1,lower.tail = T)/2)
                  }
                  
                  if(is.na(l$CI$lower) | is.infinite(l$CI$lower)){
                    check_lower <- FALSE
                  }
                  if(is.na(l$CI$upper) | is.infinite(l$CI$upper)){
                    check_upper <- FALSE
                  }
                }else{
                  check_lower <- FALSE
                  check_upper <- FALSE
                }
                
                consistency_check_lower_exact <- c(consistency_check_lower_exact,check_lower)
                consistency_check_upper_exact <- c(consistency_check_upper_exact,check_upper)
                names(consistency_check_lower_exact)[length(consistency_check_lower_exact)] <- par
                names(consistency_check_upper_exact)[length(consistency_check_upper_exact)] <- par
                
                if(length(consistency_check_approx)>0){
                  if(!consistency_check_approx[par] & check_lower & check_upper){
                    consistency_check_exact <- c(consistency_check_exact,T)
                    dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower='exact', method_upper='exact'))
                    if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                      dat_profile_tmp[[as.character(l$CI$name)]] <- cbind(dat_g[[as.character(l$CI$name)]]$Profile[,c(1:5, which(colnames(dat_g[[as.character(l$CI$name)]]$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "exact")
                    }else{
                      dat_profile_tmp[[as.character(l$CI$name)]] <- cbind(dat_g[[as.character(l$CI$name)]]$Profile[,c(1:5, which(colnames(dat_g[[as.character(l$CI$name)]]$Profile)==as.character(l$CI$name)))], method = "exact")
                    }
                  }else{
                    consistency_check_exact <- c(consistency_check_exact, F)
                  }
                }else{
                  if(check_lower & check_upper){
                    consistency_check_exact <- c(consistency_check_exact,T)
                    dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower='exact', method_upper="exact"))
                    if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                      dat_profile_tmp[[as.character(l$CI$name)]] <- cbind(dat_g[[as.character(l$CI$name)]]$Profile[,c(1:5, which(colnames(dat_g[[as.character(l$CI$name)]]$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "exact")
                    }else{
                      dat_profile_tmp[[as.character(l$CI$name)]] <- cbind(dat_g[[as.character(l$CI$name)]]$Profile[,c(1:5, which(colnames(dat_g[[as.character(l$CI$name)]]$Profile)==as.character(l$CI$name)))], method = "exact")
                    }
                  }else{
                    consistency_check_exact <- c(consistency_check_exact, F)
                  }
                } 
                names(consistency_check_exact)[length(consistency_check_exact)] <- par
              }
              
              if(length(consistency_check_approx)>0){
                consistency_check <- consistency_check | consistency_check_exact[names(consistency_check)]
              }else{
                consistency_check <- consistency_check_exact
              }
            }
          }
          
          # check if one side was consistent with approx or exact method
          if(!is.null(ncol(dat_CI_tmp)) & !('lower_min' %in% colnames(dat_CI_tmp))){
            dat_CI_tmp <- cbind(dat_CI_tmp, lower_min=NA, upper_max=NA)
          }
          
          if(length(consistency_check) == 0 || !all(consistency_check) || any(is.na(consistency_check))){
            for (par in colnames(all_param[[condition]][[batch]][[gene_name]])){
              if(par %in% c("value", "spar", "sigmat", "sigmab", "ca0_rep1")){
                next
              }
              if(length(consistency_check) == 0 || !(par %in% names(consistency_check)) || !consistency_check[par] || is.na(consistency_check[par])){
                # if approx method lower bound worked use that, else if exact method lower bound worked use that, else use optim
                if (length(consistency_check_approx) !=0 && par %in% names(consistency_check_lower_approx) && consistency_check_lower_approx[par]){
                  l <- dat_profile_tmp_all[["approx"]][[par]]
                  dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower='approximate', method_upper=NA, lower_min = NA, upper_max = NA))
                  if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                    dat_profile_tmp[[par]] <- cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "approximate")
                  }else{
                    dat_profile_tmp[[par]] <- cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==as.character(l$CI$name)))], method = "approximate")
                  }
                }else if(length(consistency_check_exact) !=0 && par %in% names(consistency_check_lower_exact) && consistency_check_lower_exact[par]) {
                  l <- dat_profile_tmp_all[["exact"]][[par]]
                  dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower='exact', method_upper=NA, lower_min = NA, upper_max = NA))
                  if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                    dat_profile_tmp[[par]] <- cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "exact")
                  }else{
                    dat_profile_tmp[[par]] <- cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==as.character(l$CI$name)))], method = "exact")
                  }
                }else{
                  lower <- min(all_param[[condition]][[batch]][[gene_name]][which(all_param[[condition]][[batch]][[gene_name]][,"value"] <= min_val + 1/2 * qchisq(0.95,1)), par])
                  lower_min <- -Inf
                  dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(name=par, value=NA, lower = lower, upper = NA, method_lower = 'optim', method_upper = NA, lower_min = ifelse(length(lower_min)>0,lower_min,-Inf), upper_max = NA))
                }
                
                # if approx method upper bound worked use that, else if exact method lower bound worked use that, else use optim
                if (length(consistency_check_approx) != 0 && par %in% names(consistency_check_upper_approx) && consistency_check_upper_approx[par]){
                  l <- dat_profile_tmp_all[["approx"]][[par]]
                  if(par %in% dat_CI_tmp$name){
                    dat_CI_tmp$upper[dat_CI_tmp$name == par] <- l$CI$upper 
                    dat_CI_tmp$method_upper[dat_CI_tmp$name == par] <- 'approximate'
                    if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  rbind(dat_profile_tmp[[as.character(l$CI$name)]],cbind(l$Profile[l$Profile$constraint >= 0, c(1:5, which(colnames(l$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "approximate"))
                    }else{
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  rbind(dat_profile_tmp[[as.character(l$CI$name)]],cbind(l$Profile[l$Profile$constraint >= 0, c(1:5, which(colnames(l$Profile)==as.character(l$CI$name)))], method = "approximate"))
                    }
                  }else{
                    if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "approximate")
                    }else{
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==as.character(l$CI$name)))], method = "approximate")
                    }
                    dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower=NA, method_upper='approximate', lower_min = NA, upper_max = NA))
                  }
                }else if( length(consistency_check_exact) != 0 && par %in% names(consistency_check_upper_exact) && consistency_check_upper_exact[par]) {
                  l <- dat_profile_tmp_all[["exact"]][[par]]
                  if(par %in% dat_CI_tmp$name){
                    dat_CI_tmp$upper[dat_CI_tmp$name == par] <- l$CI$upper 
                    dat_CI_tmp$method_upper[dat_CI_tmp$name == par] <- 'exact'
                    if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  rbind(dat_profile_tmp[[as.character(l$CI$name)]],cbind(l$Profile[l$Profile$constraint >= 0, c(1:5, which(colnames(l$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "exact"))
                    }else{
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  rbind(dat_profile_tmp[[as.character(l$CI$name)]],cbind(l$Profile[l$Profile$constraint >= 0, c(1:5, which(colnames(l$Profile)==as.character(l$CI$name)))], method = "exact"))
                    }
                  }else{
                    if(l$CI$name %in% c("min(k2,kdeg)","max(k2,kdeg)")){
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  rbind(dat_profile_tmp[[as.character(l$CI$name)]],cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==unique(dat_g[[as.character(l$CI$name)]]$Profile$whichPar)))], method = "exact"))
                    }else{
                      dat_profile_tmp[[as.character(l$CI$name)]] <-  rbind(dat_profile_tmp[[as.character(l$CI$name)]],cbind(l$Profile[l$Profile$constraint <= 0, c(1:5, which(colnames(l$Profile)==as.character(l$CI$name)))], method = "exact"))
                    }
                    dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(l$CI, method_lower=NA, method_upper='exact', lower_min = NA, upper_max = NA))
                  }
                }else{
                  upper <- max(all_param[[condition]][[batch]][[gene_name]][which(all_param[[condition]][[batch]][[gene_name]][,"value"] <= min_val + 1/2 * qchisq(0.95,1)), par])
                  upper_max <- Inf
                  if(par %in% dat_CI_tmp$name){
                    dat_CI_tmp$upper[dat_CI_tmp$name == par] <- upper
                    dat_CI_tmp$upper_max[dat_CI_tmp$name == par] <- ifelse(length(upper_max)>0,upper_max,Inf)
                    dat_CI_tmp$method_upper[dat_CI_tmp$name == par] <- 'optim'
                  }else{
                    dat_CI_tmp <- rbind(dat_CI_tmp, data.frame(name=par, value=NA, lower = NA, upper = upper, method_lower = NA, method_upper = 'optim', lower_min = NA, upper_max = ifelse(length(upper_max)>0,upper_max,Inf)))
                  }
                }
                row.names(dat_CI_tmp)[nrow(dat_CI_tmp)] <- par
              }
            }
            
            # check both bound for b2 k2, kdeg if disjunct
            if(batch == 'rep3' && max(dat_profile_tmp[["min(k2,kdeg)"]]$k2) < min(dat_profile_tmp[["max(k2,kdeg)"]]$kdeg)){
              if( dat_CI_tmp["max(k2,kdeg)", "method_lower"] == "optim" & dat_CI_tmp["min(k2,kdeg)", "method_upper"] != "optim" ){ # if disjunct
                dat_CI_tmp["max(k2,kdeg)", "lower_min"] <- dat_CI_tmp["min(k2,kdeg)", "upper"]
              }else if(dat_CI_tmp["max(k2,kdeg)", "method_lower"] == "optim" & dat_CI_tmp["min(k2,kdeg)", "method_lower"] != "optim" ){
                dat_CI_tmp["max(k2,kdeg)", "lower_min"] <- dat_CI_tmp["min(k2,kdeg)", "lower"]
              }
              if( dat_CI_tmp["min(k2,kdeg)", "method_upper"] == "optim" & dat_CI_tmp["max(k2,kdeg)", "method_lower"] != "optim" ){
                dat_CI_tmp["min(k2,kdeg)", "upper_max"] <- dat_CI_tmp["max(k2,kdeg)", "lower"]
              }else if(dat_CI_tmp["min(k2,kdeg)", "method_upper"] == "optim" & dat_CI_tmp["max(k2,kdeg)", "method_upper"] != "optim" ){
                dat_CI_tmp["min(k2,kdeg)", "upper_max"] <- dat_CI_tmp["max(k2,kdeg)", "upper"]
              }
            }
          }
          if(!('lower_min' %in% colnames(dat_CI_tmp))){
            dat_CI_tmp <- cbind(dat_CI_tmp, lower_min=NA, upper_max=NA)
          }
        }
        
        # reorder dat_CI_tmp
        dat_CI_tmp <- dat_CI_tmp[names(param_best)[names(param_best) %in% dat_CI_tmp$name],]
        dat_CI_tmp$value <- param_best[names(param_best) %in% dat_CI_tmp$name]
        
        colnames(dat_CI_tmp)[1] <- 'whichPar'
        all_param_CI[[condition]][[batch]][[gene_name]] <- list(CI = dat_CI_tmp, Profile = dat_profile_tmp) 
      }
    }
  }
}

save(all_annot, all_best_param, all_param, all_data, all_fit, all_fit_CI, all_library, all_param_CI, file = 'all_results.Rdata')
stop()

load(file = 'all_results.Rdata')

#### Figure 1E - Heatmap of expression ####---------------------------------------------------------------
dat <- data.frame()[1:273,]

for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- caRNA_Naive_exons_rpkm[,grep(t, colnames(caRNA_Naive_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    tmp <- tmp[rownames(tmp) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],]
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    tmp <- tmp[names(tmp) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected]]
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Chromatin.Naive.", t)
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- npRNA_Naive_exons_rpkm[,grep(t, colnames(npRNA_Naive_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    tmp <- tmp[rownames(tmp) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],]
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    tmp <- tmp[names(tmp) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected]]
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Nuc.Naive.", t)
  
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- cytoRNA_Naive_exons_rpkm[,grep(t, colnames(cytoRNA_Naive_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    tmp <- tmp[rownames(tmp) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],]
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    tmp <- tmp[names(tmp) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected]]
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Cyto.Naive.", t)
}
row.names(dat) <- row.names(caRNA_Naive_exons_rpkm)[rownames(caRNA_Naive_exons_rpkm) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected]]

dat <- dat[row.names(dat) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],] 

colors_spectral=rev(brewer.pal(11,"Spectral"))
colors_spectral=colorRampPalette(colors_spectral, space="rgb")(100)

h <- hclust(dist(as.matrix(dat[,1:13])), method = 'ward.D2')

dat_2 <- dat
dat_2[dat > 10] <- 10
dat_2[dat < -5] <- -5
pheatmap(as.matrix(dat_2), color = colors_spectral, cluster_rows = h, cluster_cols = F, show_rownames = F,gaps_col = c(13,26), breaks = seq(-5,10,length.out = 101))

# svg('Fig1E.svg')
reorder <- c(h$order[cutree(h, k=3)[h$order] == 3],h$order[cutree(h, k=3)[h$order] == 1], rev(h$order[cutree(h, k=3)[h$order] == 2]))
pheatmap(as.matrix(dat_2[reorder,]), color = colors_spectral, cluster_rows = F, cluster_cols = F, show_rownames = F,gaps_col = c(13,26), breaks = seq(-5,10,length.out = 101))
# dev.off()

stop()


#### Figure 1G - Example genes ####---------------------------------------------------------------
for (gene_name in list(c("Cd74", "Btg2"),c("Arl5b", "Ccr3"))){
  g1 <- paste(gene_list_mm89$ensembl_gene_id[gene_list_mm89$external_gene_name == gene_name[1]], gene_list_mm89$version[gene_list_mm89$external_gene_name == gene_name[1]], sep='.')
  g2 <- paste(gene_list_mm89$ensembl_gene_id[gene_list_mm89$external_gene_name == gene_name[2]], gene_list_mm89$version[gene_list_mm89$external_gene_name == gene_name[2]], sep='.')
  
  dat_g <- data.frame()
  for(t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
    for (b in c("rep1","rep2","rep3")){
      indexes_ca <- intersect(grep(t, colnames(caRNA_Naive_exons_rpkm)), grep(b, colnames(caRNA_Naive_exons_rpkm)))
      if(length(indexes_ca) == 1){
        dat_g <- rbind(dat_g, data.frame(log2RPKM = caRNA_Naive_exons_rpkm[g1,indexes_ca], Time = as.numeric(t), Fraction = "Chromatin", Batch = b, Gene=gene_name[1]))
        dat_g <- rbind(dat_g, data.frame(log2RPKM = caRNA_Naive_exons_rpkm[g2,indexes_ca], Time = as.numeric(t), Fraction = "Chromatin", Batch = b, Gene=gene_name[2]))
      }
      indexes_np <- intersect(grep(t, colnames(npRNA_Naive_exons_rpkm)), grep(b, colnames(npRNA_Naive_exons_rpkm)))
      if(length(indexes_np) == 1){
        dat_g <- rbind(dat_g,data.frame(log2RPKM = npRNA_Naive_exons_rpkm[g1,indexes_np], Time = as.numeric(t), Fraction = "Nucleoplasm", Batch = b, Gene=gene_name[1]))
        dat_g <- rbind(dat_g,data.frame(log2RPKM = npRNA_Naive_exons_rpkm[g2,indexes_np], Time = as.numeric(t), Fraction = "Nucleoplasm", Batch = b, Gene=gene_name[2]))
      }
      indexes_cyt <- intersect(grep(t, colnames(cytoRNA_Naive_exons_rpkm)), grep(b, colnames(cytoRNA_Naive_exons_rpkm)))
      if(length(indexes_cyt) == 1){
        dat_g <- rbind(dat_g,data.frame(log2RPKM = cytoRNA_Naive_exons_rpkm[g1,indexes_cyt], Time = as.numeric(t), Fraction = "Cytoplasm", Batch = b, Gene=gene_name[1]))
        dat_g <- rbind(dat_g,data.frame(log2RPKM = cytoRNA_Naive_exons_rpkm[g2,indexes_cyt], Time = as.numeric(t), Fraction = "Cytoplasm", Batch = b, Gene=gene_name[2]))
      }
    }
  }
  dat_g$Fraction <- factor(dat_g$Fraction, levels = c("Chromatin", "Nucleoplasm", "Cytoplasm"))
  p <- ggplot(dat_g, aes(x = Time, y=log2RPKM, col = Gene, group=interaction(Batch, Fraction, Gene))) + geom_point() + geom_line(alpha=0.5)
  p <- p + geom_smooth(mapping = aes(group=interaction(Fraction, Gene)), se = F) + facet_grid(~Fraction)
  p <- p + theme_bw() + theme(aspect.ratio=1, panel.background = element_rect(fill = "transparent")) + labs(y='log2(FPKM)', title=paste0(gene_name, collapse = " and "))
  p <- p + scale_x_continuous(breaks=c(0,30,60,90,120))
  print(p)
}
stop()

#### Figure 1F - Correlation of Compartment ####---------------------------------------------------------------
dat <- data.frame()[1:nrow(caRNA_Naive_exons_cpm),]

for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- caRNA_Naive_exons_cpm[,grep(t, colnames(caRNA_Naive_exons_cpm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Chromatin.Naive.", t)
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- npRNA_Naive_exons_cpm[,grep(t, colnames(npRNA_Naive_exons_cpm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Nuc.Naive.", t)
  
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- cytoRNA_Naive_exons_cpm[,grep(t, colnames(cytoRNA_Naive_exons_cpm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Cyto.Naive.", t)
}
row.names(dat) <- row.names(caRNA_Naive_exons_cpm)


head(dat)
dat2 <- dat[row.names(dat) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],]

# Average
dat_cor <- cor(dat)
dat_cor2 <- cor(dat2)

dim(dat_cor)

cor_mat <- dat_cor
range(dat_cor)
range(dat_cor2)

diag(cor_mat) <- NA

dat_cor2[14:26,1:13] <- NA
dat_cor2[27:39,1:26] <- NA

pheatmap(dat_cor2,cluster_rows = F, cluster_cols = F, gaps_row=c(13,26), gaps_col=c(13,26), na_col = "white")
stop()

#### Figure S1A - Correlation between replicates ####--------------------------------------------------------------
dat <- cbind(caRNA_Naive_exons_cpm[,c(NA, 24, 12,1,25,13,2,26,14,3,27,15,4,28,16,5,29,17,6,7,30,18,8,31,19,9,32,20,33,21,10,34,22,11,35,23)], 
             npRNA_Naive_exons_cpm[,c(NA,12,NA,5,13,NA,1,14,NA,6,15,NA,2,16,NA,7,17,NA,3,8,18,NA,4,19,NA,9,20,NA,21,NA,10,22,NA,11,23,NA)], 
             cytoRNA_Naive_exons_cpm[,c(5,25,13,6,26,14,1,27,15,7,28,16,2,29,17,8,30,18,3,9,31,19,4,32,20,10,33,21,34,22,11,35,23,12,36,24)]) # switched b2 and b3 in the data
colnames(dat)
dat2 <- dat[row.names(dat) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],]

dat_cor <- cor(dat)
dat_cor2 <- cor(dat2)

dat_cor2[1:(12*3),-c(1:(12*3))] <- NA
dat_cor2[37:72,c(1:36,73:108)] <- NA
dat_cor2[73:108,-c(73:108)] <- NA

dat_cor2[1,] <- NA
dat_cor2[36 + c(1,3,6,9,12,15,18,22,25,28,30,33,36),] <- NA
dat_cor2[,36 + c(1,3,6,9,12,15,18,22,25,28,30,33,36)] <- NA

dim(dat_cor)

cor_mat <- dat_cor
range(dat_cor)
range(dat_cor2)

pheatmap(dat_cor2, cluster_rows = F, cluster_cols = F, gaps_row = c(36,72), gaps_col = c(36,72), na_col = "white")
stop()

#### Figure S1B - PCA of technical replicates ####--------------------------------------------------------------
# => See PCA_DEGs.R (500x396)

#### Figure S1C - Heatmap of Naive replicates ####--------------------------------------------------------------
# 500*591
dat_ca <- caRNA_Naive_exons_rpkm[match(list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],rownames(caRNA_Naive_exons_rpkm)),]
dat_ca <- dat_ca[,c(1:11,24:35,12:23)]
dat_np <- npRNA_Naive_exons_rpkm[match(list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],rownames(npRNA_Naive_exons_rpkm)),]
dat_np <- dat_np[,c(5,1,6,2,7,3,8,4,9,10,11,12:23)] 
dat_cyto <- cytoRNA_Naive_exons_rpkm[match(list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],rownames(cytoRNA_Naive_exons_rpkm)),]
dat_cyto <- dat_cyto[,c(5,6,1,7,2,8,3,9,4,10,11,12,13:36)] 
dat_cyto <- dat_cyto[,c(1:12,25:36,13:24)] #switch 2 and 3

dat <- cbind(dat_ca, dat_np, dat_cyto)

colors_spectral <- rev(brewer.pal(11,"Spectral"))
colors_spectral <- colorRampPalette(colors_spectral, space="rgb")(100)

dat_2 <- dat
dat_2[dat > 10] <- 10
dat_2[dat < -5] <- -5

pheatmap(as.matrix(dat_2[match(h$labels,rownames(dat_2)),][reorder,]), color = colors_spectral, cluster_rows = F, cluster_cols = F, show_rownames = F,gaps_col = c(35,35+23), breaks = seq(-5,10,length.out = 101))
stop()

#### Figure S5B - Heatmap of LPA replicates ####--------------------------------------------------------------
dat_ca <- caRNA_LPA_exons_rpkm[match(list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],rownames(caRNA_LPA_exons_rpkm)),]
dat_ca <- dat_ca[,c(1:12,25:36,13:24)]
dat_np <- npRNA_LPA_exons_rpkm[match(list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],rownames(npRNA_LPA_exons_rpkm)),]
dat_np <- dat_np[,c(5,6,1,7,2,8,3,9,4,10,11,12:24)] # reorder
dat_cyto <- cytoRNA_LPA_exons_rpkm[match(list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],rownames(cytoRNA_LPA_exons_rpkm)),]
dat_cyto <- dat_cyto[,c(5,6,1,7,2,8,3,9,4,10,11,12,13:36)] # reorder
dat_cyto <- dat_cyto[,c(1:12,25:36,13:24)] #switch 2 and 3

dat <- cbind(dat_ca, dat_np, dat_cyto)

colors_spectral=rev(brewer.pal(11,"Spectral"))
colors_spectral=colorRampPalette(colors_spectral, space="rgb")(100)

dat_2 <- dat
dat_2[dat > 10] <- 10
dat_2[dat < -5] <- -5

#remove the ones that will not be done for LPA due to low counts, and not high enough expression
dat_2[c("ENSMUSG00000001348.15", "ENSMUSG00000027358.6",  "ENSMUSG00000007655.16", "ENSMUSG00000035448.9",  "ENSMUSG00000032578.7",  "ENSMUSG00000040663.8", 
        "ENSMUSG00000029417.9",  "ENSMUSG00000018648.15", "ENSMUSG00000003545.3",  "ENSMUSG00000021303.14", "ENSMUSG00000024486.5",  "ENSMUSG00000034997.4",
        "ENSMUSG00000025383.2",  "ENSMUSG00000034394.14", "ENSMUSG00000028862.6",  "ENSMUSG00000057440.7",  "ENSMUSG00000005125.12", "ENSMUSG00000011179.8", 
        "ENSMUSG00000028583.14", "ENSMUSG00000020205.8",  "ENSMUSG00000054855.13", "ENSMUSG00000061878.15", "ENSMUSG00000026700.4",  "ENSMUSG00000045672.15", 
        "ENSMUSG00000025804.5",  "ENSMUSG00000003352.13", "ENSMUSG00000028341.9",  "ENSMUSG00000018916.5",  "ENSMUSG00000043421.8",  "ENSMUSG00000006205.13",
        "ENSMUSG00000029819.6",  "ENSMUSG00000058013.11"), ] <- NA

pheatmap(as.matrix(dat_2[match(h$labels,rownames(dat_2)),][reorder,]), color = colors_spectral, cluster_rows = F, cluster_cols = F, show_rownames = F,gaps_col = c(36,36+24), breaks = seq(-5,10,length.out = 101))
stop()

#### Figure S2G - Comparison of TSS and TES manual annotation to DB ####--------------------------------------------------------------
## See ../Comparison with TSS-TES DB/TSS_polyA_distance_for_table.R

#### Figure 6A - Heatmap Naive and LPA FC ####--------------------------------------------------------------
dat <- data.frame()[1:273,]
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- caRNA_Naive_exons_rpkm[rownames(caRNA_Naive_exons_rpkm) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected], grep(t, colnames(caRNA_Naive_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Chromatin.Naive.", t)
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- npRNA_Naive_exons_rpkm[rownames(caRNA_Naive_exons_rpkm) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected], ][,grep(t, colnames(npRNA_Naive_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  } else {
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Nuc.Naive.", t)
  
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- cytoRNA_Naive_exons_rpkm[rownames(caRNA_Naive_exons_rpkm) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected], ][,grep(t, colnames(cytoRNA_Naive_exons_rpkm))]
  if (!is.null(ncol(tmp))) {
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Cyto.Naive.", t)
}
row.names(dat) <- rownames(caRNA_Naive_exons_rpkm)[rownames(caRNA_Naive_exons_rpkm) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected]]

for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- caRNA_LPA_exons_rpkm[rownames(caRNA_LPA_exons_rpkm)  %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected] , grep(t, colnames(caRNA_LPA_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Chromatin.LPA.", t)
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- npRNA_LPA_exons_rpkm[rownames(caRNA_LPA_exons_rpkm)  %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],grep(t, colnames(npRNA_LPA_exons_rpkm))]
  if(!is.null(ncol(tmp))){
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  } else {
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Nuc.LPA.", t)
  
}
for (t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  tmp <- cytoRNA_LPA_exons_rpkm[rownames(caRNA_LPA_exons_rpkm)  %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],grep(t, colnames(cytoRNA_LPA_exons_rpkm))]
  if (!is.null(ncol(tmp))) {
    dat <- cbind(dat,log2(rowMeans(2^tmp, na.rm = T)))
  }else{
    dat <- cbind(dat,tmp)
  }
  colnames(dat)[ncol(dat)] <- paste0("Cyto.LPA.", t)
}

dat <- dat[row.names(dat) %in% list_genes_cleaned$Geneid[list_genes_cleaned$Naive.Selected],] 

dat_2 <- dat
dat_2[dat > 10] <- 10
dat_2[dat < -5] <- -5
dat_2[!(rownames(dat_2) %in% list_genes_cleaned$Geneid[list_genes_cleaned$LPA.Selected]),40:78] <- NA

## plot FC
head(dat)
datFC <- data.frame(dat[,1:13] - dat[,1],dat[,14:26]- dat[,14],dat[,27:39] - dat[,27], dat[,40:52]-dat[,40], dat[,53:65]-dat[,53], dat[,66:78]-dat[,66])
datFC[datFC > 10] <- 10

# Remove LPA genes that were not selected for modeling
datFC[!(rownames(dat_2) %in% list_genes_cleaned$Geneid[list_genes_cleaned$LPA.Selected]),40:78] <- NA

# Change column names
colnames(datFC) <- gsub("Nuc", "Nucleoplasm", colnames(datFC))
colnames(datFC) <- gsub("Cyto", "Cytoplasm", colnames(datFC))

colors_spectral <- rev(brewer.pal(11,"Spectral"))
colors_spectral <- colorRampPalette(colors_spectral, space="rgb")(101)

change_legend <- function(lgd){
  lgd@grob$children[[grep("legend.body",names(lgd@grob$children))]]$children[[3]] <- NULL
  lgd@grob$children[[grep("legend.body",names(lgd@grob$children))]]$children[[3]]$gp$col <- "black"
  segment_length <- abs(lgd@grob$children[[grep("legend.body",names(lgd@grob$children))]]$children[[3]]$x0 - lgd@grob$children[[grep("legend.body",names(lgd_Heatmap@grob$children))]]$children[[3]]$x1)
  lgd@grob$children[[grep("legend.body",names(lgd@grob$children))]]$children[[3]]$x0 <- segment_length/2 + lgd@grob$children[[grep("legend.body",names(lgd_Heatmap@grob$children))]]$children[[3]]$x0
  lgd@grob$children[[grep("legend.body",names(lgd@grob$children))]]$children[[3]]$x1 <- segment_length/2 + lgd@grob$children[[grep("legend.body",names(lgd_Heatmap@grob$children))]]$children[[3]]$x1
  return(lgd)
}

h_width <- unit(24,"mm")
h_height <- unit(90, "mm")
h_gap <- unit(1,"mm")
legend_width <- unit(2,"mm")
title_size <- 12
subtitle_size <- 8
annot_size <- 6

# Color functions
Heatmap_col_fun <- colorRamp2(breaks=seq(-max(datFC, na.rm=T), max(datFC, na.rm=T), length.out=101)[-c(1:(floor((min(datFC, na.rm=T)+max(datFC, na.rm=T))/(2*max(datFC, na.rm=T))*100)-1))], colors_spectral[-c(1:(floor((min(datFC, na.rm=T)+max(datFC, na.rm=T))/(2*max(datFC, na.rm=T))*100)-1))])


ht_opt(TITLE_PADDING = unit(1, "mm"))

Time_labels <- as.numeric(unlist(lapply(strsplit(colnames(datFC), '.', fixed=T), function(l){l[3]})))

Heat_naive <- Heatmap(matrix = as.matrix(datFC[match(h$labels, row.names(dat_2)),][reorder,])[,grep("Naive", colnames(datFC))], 
                      col = Heatmap_col_fun,
                      show_heatmap_legend = F,
                      height = h_height,
                      width = 3 * h_width + 2 * h_gap,
                      column_gap = h_gap,
                      column_split = factor(unlist(lapply(strsplit(colnames(datFC)[grep("Naive", colnames(datFC))], '.', fixed=T), function(l){l[1]})), levels=c("Chromatin", "Nucleoplasm","Cytoplasm")),
                      column_names_gp = gpar(fontsize=annot_size),
                      show_column_names = T,
                      show_row_names = F,
                      column_labels = structure(replace(Time_labels[grep("Naive", colnames(datFC))],!(Time_labels[grep("Naive", colnames(datFC))] %in% c(0,30,60,120)),""), names = Time_labels[grep("Naive", colnames(datFC))]),
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title_gp = gpar(fontsize = subtitle_size),
                      na_col = "grey",
                      name="Naive"
)

Heat_LPA <- Heatmap(matrix = as.matrix(datFC[match(h$labels, row.names(dat_2)),][reorder,])[,grep("LPA", colnames(datFC))], 
                    col = Heatmap_col_fun,
                    show_heatmap_legend = F,
                    height = h_height,
                    width = 3 * h_width + 2 * h_gap,
                    column_gap = h_gap,
                    column_split = factor(unlist(lapply(strsplit(colnames(datFC)[grep("LPA", colnames(datFC))], '.', fixed=T), function(l){l[1]})), levels=c("Chromatin", "Nucleoplasm","Cytoplasm")),
                    column_names_gp = gpar(fontsize=annot_size),
                    show_column_names = T,
                    show_row_names = F,
                    column_labels = structure(replace(Time_labels[grep("LPA", colnames(datFC))],!(Time_labels[grep("LPA", colnames(datFC))] %in% c(0,30,60,120)),""), names = Time_labels[grep("LPA", colnames(datFC))]),
                    cluster_rows = F,
                    cluster_columns = F,
                    column_title_gp = gpar(fontsize = subtitle_size),
                    na_col = "grey",
                    name="LPA"
)

lgd_Heatmap <- Legend(title=expression(log[2](FC)), col= Heatmap_col_fun, labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width, at = seq(0,10, by=2))
lgd_Heatmap <- change_legend(lgd_Heatmap)

gb_hnaive <- grid.grabExpr(draw(Heat_naive, column_title="Naive", column_title_gp = gpar(fontsize = title_size, fontface='bold')))
gb_hlpa <- grid.grabExpr(draw(Heat_LPA, column_title="Tolerized", column_title_gp = gpar(fontsize = title_size, fontface='bold')))
gb_legend <- grid.grabExpr(draw(lgd_Heatmap))

color_subtitle <- function(gb){
  for(vp in names(gb$children)){
    if(grepl("text", vp)){
      if(length(gb$children[[vp]]$label) == 1 ){
        if(gb$children[[vp]]$label == "Chromatin"){
          gp <- c(gb$children[[vp]]$gp, gpar(col=col_ca))
          class(gp) <- "gpar"
          gb$children[[vp]]$gp <- gp
        }else if(gb$children[[vp]]$label == "Nucleoplasm"){
          gp <- c(gb$children[[vp]]$gp, gpar(col=col_np))
          class(gp) <- "gpar"
          gb$children[[vp]]$gp <- gp
        }else if(gb$children[[vp]]$label == "Cytoplasm"){
          gp <- c(gb$children[[vp]]$gp, gpar(col=col_cyto))
          class(gp) <- "gpar"
          gb$children[[vp]]$gp <- gp
        }
      }
    }
  }
  return(gb)
}

gb_hnaive <- color_subtitle(gb_hnaive)
gb_hlpa <- color_subtitle(gb_hlpa)

vp_naive <- viewport(name="naive", x = 2*h_gap, y = 0, width= 3*h_width+2*h_gap, just = c("left", "bottom"))
vp_lpa <- viewport(name="LPA", x = 2*h_gap + vp_naive$x + vp_naive$width , y = vp_naive$y, width= 3*h_width + 2*h_gap, just = c("left", "bottom"))
vp_legend <- viewport(name="legend", x = 2 * h_gap + vp_lpa$x + vp_lpa$width, y=0, width = lgd_Heatmap@grob$vp$width, just=c(0,0))

#svg(file = 'Figures/Figure6A_Heatmap_FC_Naive_LPA.svg', width = 7, height = 7)
grid.newpage()
pushViewport(vp_naive)
grid.draw(gb_hnaive)
popViewport()

pushViewport(vp_lpa)
grid.draw(gb_hlpa)
popViewport()

pushViewport(vp_legend)
grid.draw(gb_legend)
popViewport()
#dev.off()
stop()

# #### Figure 6B - Comparison of basal ####------------------------------------------
cor_test <- cor.test(dat$Chromatin.Naive.000, dat$Chromatin.LPA.000)
# svg(filename = "Figures/Fig6_canaive_vs_catol.svg", width=5, height=5)
p <- ggplot(data=dat, aes(x=Chromatin.Naive.000, y=Chromatin.LPA.000)) + geom_point()
p <- p + coord_cartesian(xlim=c(-10,7),ylim=c(-10,7)) + theme_bw(base_size = 24) + theme(aspect.ratio = 1)
p <- p + geom_abline(aes(intercept=0, slope=1), col="green")
p <- p + labs(x=expression(log[2](FPKM) ~ " - Naive Basal") , y= expression(log[2](FPKM) ~ " - Tolerized Basal"), parse=T)
p
# dev.off()

datFC2 <- data.frame(naive_FC = apply(datFC[, grep("Chromatin.Naive", colnames(datFC))], 1 ,max),
                     lpa_FC = apply(datFC[, grep("Chromatin.LPA", colnames(datFC))], 1 ,max))

# svg(filename = "Figures/Fig6_FC_canaive_vs_catol.svg", width=5, height=5)
p <- ggplot(data=datFC2, aes(x=naive_FC, y=lpa_FC)) + geom_point()
p <- p + coord_cartesian(xlim=c(1,12.5),ylim=c(1,12.5)) + theme_bw(base_size = 24) + theme(aspect.ratio = 1)
p <- p + geom_abline(aes(intercept=0, slope=1), col="green")
p <- p + labs(x=expression(log[2](FC) ~ " - Naive") , y= expression(log[2](FC) ~ " - Tolerized"), parse=T)
p
# dev.off()
stop()

#### Figure 2B - Example plots for fitting workflow ####---------------------------------------------------------------
gene_name <- "Arl5b"
g1 <- paste(gene_list_mm89$ensembl_gene_id[gene_list_mm89$external_gene_name == gene_name], gene_list_mm89$version[gene_list_mm89$external_gene_name == gene_name], sep='.')

dat_g <- data.frame()
for(t in c("000","010","015","020","025","030","035","040","050","060","075","090","120")){
  for (b in c("rep1","rep2","rep3")){
    indexes_ca <- intersect(grep(t, colnames(caRNA_Naive_exons_rpkm_5kb)), grep(b, colnames(caRNA_Naive_exons_rpkm_5kb)))
    if(length(indexes_ca) == 1){
      dat_g <- rbind(dat_g, data.frame(log2RPKM = caRNA_Naive_exons_rpkm_5kb[g1,indexes_ca], Time = as.numeric(t), Fraction = "Chromatin", Batch = b, Gene=gene_name[1]))
    }
    indexes_np <- intersect(grep(t, colnames(npRNA_Naive_exons_rpkm_5kb)), grep(b, colnames(npRNA_Naive_exons_rpkm_5kb)))
    if(length(indexes_np) == 1){
      dat_g <- rbind(dat_g,data.frame(log2RPKM = npRNA_Naive_exons_rpkm_5kb[g1,indexes_np], Time = as.numeric(t), Fraction = "Nucleoplasm", Batch = b, Gene=gene_name[1]))
    }
    indexes_cyt <- intersect(grep(t, colnames(cytoRNA_Naive_exons_rpkm_5kb)), grep(b, colnames(cytoRNA_Naive_exons_rpkm_5kb)))
    if(length(indexes_cyt) == 1){
      dat_g <- rbind(dat_g,data.frame(log2RPKM = cytoRNA_Naive_exons_rpkm_5kb[g1,indexes_cyt], Time = as.numeric(t), Fraction = "Cytoplasm", Batch = b, Gene=gene_name[1]))
    }
  }
}
dat_g$Fraction <- factor(dat_g$Fraction, levels = c("Chromatin", "Nucleoplasm", "Cytoplasm"))


### Example for Arl5b b3

for (b in "rep2"){
  p <- ggplot(dat_g[dat_g$Batch == b & dat_g$Fraction == "Chromatin",], aes(x = Time, y=log2RPKM)) + geom_point()
  p <- p + theme_bw() + theme(aspect.ratio=1, panel.background = element_rect(fill = "transparent")) + labs(y='log2(RPKM)', title=paste0(gene_name, collapse = " and "))
  p <- p + scale_x_continuous(breaks=c(0,30,60,90,120))
  print(p)
}

ca_smooth <- predict(smooth.spline(x=dat_g$Time[ dat_g$Fraction == "Chromatin" & dat_g$Batch == "rep2"], y=dat_g$log2RPKM[ dat_g$Fraction == "Chromatin" & dat_g$Batch == "rep2"], spar = 0.45),seq(0,121,by=1/60))$y
p <- ggplot(dat_g[ dat_g$Fraction == "Chromatin" & dat_g$Batch == "b3",], aes(x = Time, y=log2RPKM)) + geom_point(alpha=0.3)
p <- p + geom_line(data = data.frame(x=seq(0,121,1/60), y=ca_smooth), aes(x=x, y=y),col="black", size=1) 
p <- p + theme_bw() + theme(aspect.ratio=1, panel.background = element_rect(fill = "transparent")) + labs(y='log2(RPKM)', title=paste0(gene_name, collapse = " and "))
p <- p + scale_x_continuous(breaks=c(0,30,60,90,120))
print(p)

for (b in "rep2"){
  p <- ggplot(dat_g[dat_g$Batch == b & dat_g$Fraction == "Nucleoplasm",], aes(x = Time, y=log2RPKM)) + geom_point()
  p <- p + theme_bw() + theme(aspect.ratio=1, panel.background = element_rect(fill = "transparent")) + labs(y='log2(RPKM)', title=paste0(gene_name, collapse = " and "))
  p <- p + scale_x_continuous(breaks=c(0,30,60,90,120))
  print(p)
}

for (b in "rep2"){
  p <- ggplot(dat_g[dat_g$Batch == b & dat_g$Fraction == "Cytoplasm",], aes(x = Time, y=log2RPKM)) + geom_point()
  p <- p + theme_bw() + theme(aspect.ratio=1, panel.background = element_rect(fill = "transparent")) + labs(y='log2(RPKM)', title=paste0(gene_name, collapse = " and "))
  p <- p + scale_x_continuous(breaks=c(0,30,60,90,120))
  print(p)
}

param <- c(20,10,0.02,0.02)

fit <- fn_all_ode_cp(0:120, param, convert_param = F, log_in = F, cbind(seq(0,121,1/60),2^ca_smooth), u = input_fun_cp, no_np=FALSE)

dat_g <- data.frame(Time=rep(fit[,"time"],2), y=c(log2(fit[,"npRNA"]),log2(fit[,"cytoRNA"])), cpt=factor(rep(c("Nucleoplasmic", "Cytoplasmic"), each=nrow(fit)),levels =  c("Nucleoplasmic", "Cytoplasmic"))) 
p <- ggplot(dat_g, aes(x=Time, y=y)) + geom_line(size=1) + facet_grid(.~cpt)
p + theme_bw() + theme(aspect.ratio = 1) + scale_x_continuous(breaks=c(0,30,60,90,120))
stop()

#### Figure 2C - Heatmap of data vs Fit Naive ####---------------------------------------------------------------
dat_sim <- list()
annot <- list()

for ( condition in c("naive", "lpa")){
  dat_sim_tmp <- data.frame()
  annot_tmp <- data.frame()
  for (g in names(all_data[[condition]])){
    print(g)
    
    batch <- 'rep2'
    
    data_all <- all_data[[condition]][[g]]  
    data <- list(caRNA = data_all$caRNA[,c("t", "rep2")],
                 npRNA = data_all$npRNA[,c("t", "rep2")],
                 cytoRNA = data_all$cytoRNA[,c("t", "rep2")])
    fit <- all_fit[[condition]][[batch]][[g]]
    
    annot_tmp <- rbind(annot_tmp, unlist(all_annot[[condition]][[batch]][[g]][["annot"]]))
    row.names(annot_tmp)[nrow(annot_tmp)] <- g
    
    dat_sim_tmp <- rbind(dat_sim_tmp, c(data$caRNA[!is.na(data$caRNA[,-1]),-1], 
                                        fit$caRNA[fit$caRNA[,1] %in% data$caRNA[!is.na(data$caRNA[,-1]),1],batch], 
                                        data$npRNA[!is.na(data$npRNA[,-1]),-1], 
                                        fit$npRNA[fit$npRNA[,1] %in% data$npRNA[!is.na(data$npRNA[,-1]),1],batch],
                                        data$cytoRNA[!is.na(data$cytoRNA[,-1]),-1], 
                                        fit$cytoRNA[fit$cytoRNA[,1] %in% data$cytoRNA[!is.na(data$cytoRNA[,-1]),1],batch]))
    
    row.names(dat_sim_tmp)[nrow(dat_sim_tmp)] <- g
  }
  
  colnames(dat_sim_tmp) <- c(rep(data$caRNA[!is.na(data$caRNA[,-1]),1],2),rep(data$npRNA[!is.na(data$caRNA[,-1]),1],2),rep(data$cytoRNA[!is.na(data$cytoRNA[,-1]),1],2))
  colnames(annot_tmp) <- c( "max(|autocor(error)|*mean(|error|)/range)", 
                            "max((1+|autocor(error)|)*mean(|error|)/range)", 
                            "max(|w_autocor(error)|*w_mean(|error|)/range)", 
                            "max((1+|w_autocor(error)|)*w_mean(|error|)/range)", 
                            "max(w_mean(|error|)/range)",
                            "NLL", "NLL_ca","NLL_np","NLL_cyto","NLL_theo", "NNL-NLL_theo") 
  
  dat_sim[[condition]] <- dat_sim_tmp
  annot[[condition]] <- annot_tmp
}

colors_spectral <- rev(brewer.pal(11,"Spectral"))
colors_spectral <- colorRampPalette(colors_spectral, space="rgb")(101)

h_width <- unit(14, "mm")
h_height <- unit(90, "mm")
h_gap <- unit(1,"mm")
annot_width <- unit(2, 'mm')
legend_width <- unit(2,"mm")
title_size <- 12
subtitle_size <- 8
annot_size <- 6

dat_2 <- as.matrix(dat_sim[['naive']])
dat_2[dat_2 > 10] <- 10
dat_2[dat_2 < -5] <- -5

dat_3 <- dat_2[match(h$labels,paste0(gene_infos$ensembl_gene_id,'.',gene_infos$version)[match(row.names(dat_2),gene_infos$external_gene_name)]),][reorder,]
dat_3 <- dat_3[!apply(dat_3,1,function(x){all(is.na(x))}),]

dat_annot <- annot[['naive']][match(h$labels,paste0(gene_infos$ensembl_gene_id,'.',gene_infos$version)[match(row.names(annot[['naive']]),gene_infos$external_gene_name)]),][reorder,]
dat_annot <- dat_annot[!apply(dat_annot,1,function(x){all(is.na(x))}),]

dat_annot2 <- data.frame(dat_annot[,1])
colnames(dat_annot2) <- "max(|autocor(error)|*mean(|error|)/range)"

# Color functions
Heatmap_col_fun <- colorRamp2(breaks=seq(min(dat_3), max(dat_3), length.out=101), colors_spectral)
Qual_fit_col_fun <- colorRamp2(breaks=seq(min(dat_annot[,1]), max(dat_annot[,1]),length.out=101),colorRampPalette(c("green3","gold2", "darkorange2","darkred"))(101))
NLL_col_fun <- colorRamp2(breaks=seq(min(dat_annot[,c(6,10)]), max(dat_annot[,c(6,10)]),length.out=101),colorRampPalette(c("white","darkred"))(101))
NLL_minus_NLLtheo_col_fun <- colorRamp2(breaks=seq(0, max(dat_annot[,c(6,10)])-min(dat_annot[,c(6,10)]),length.out=101),colorRampPalette(c("white","darkred"))(101))

ht_opt(TITLE_PADDING = unit(1, "mm"))
hca <- Heatmap(dat_3[,1:24], col = Heatmap_col_fun,
               cluster_rows=F, cluster_columns=F, 
               show_row_names=F, show_column_names=T, 
               column_names_gp = gpar(fontsize=annot_size),
               column_labels = structure(replace(colnames(dat_3)[1:24],!(colnames(dat_3)[1:24] %in% c(0,30,60,120)),""), names = colnames(dat_3)[1:24]), 
               column_split = rep(c("Data","Smoothed"),each=12), 
               column_title_gp = gpar(fontsize = subtitle_size), column_gap = h_gap, width=2*h_width+h_gap, height = h_height, show_heatmap_legend = F)
hnp <- Heatmap(dat_3[,25:48], col=Heatmap_col_fun, 
               cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, column_names_gp = gpar(fontsize=annot_size),
               column_labels = structure(replace(colnames(dat_3)[1:24],!(colnames(dat_3)[1:24] %in% c(0,30,60,120)),""), names = colnames(dat_3)[1:24]),
               column_split = rep(c("Data","Simulation"),each=12), column_title_gp = gpar(fontsize = subtitle_size), column_gap = h_gap, width=2*h_width+h_gap, height = h_height, show_heatmap_legend = F)
hcyto <- Heatmap(dat_3[,49:72], col=Heatmap_col_fun, 
                 cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, column_names_gp = gpar(fontsize=annot_size),
                 column_labels = structure(replace(colnames(dat_3)[1:24],!(colnames(dat_3)[1:24] %in% c(0,30,60,120)),""), names = colnames(dat_3)[1:24]), 
                 column_split = rep(c("Data","Simulation"),each=12), column_title_gp = gpar(fontsize = subtitle_size), column_gap = h_gap, width=2*h_width+h_gap, height = h_height, show_heatmap_legend = F)

annot2 <- Heatmap(matrix = as.matrix(dat_annot2), col = Qual_fit_col_fun, show_column_names = T, column_names_gp = gpar(fontsize=annot_size), column_labels = "Fit Quality",
                  cluster_rows=F, cluster_columns=F, show_row_names = F, show_heatmap_legend = F, width=annot_width, height=h_height)

annot51 <- Heatmap(as.matrix(dat_annot[,c(6,10)]),col = NLL_col_fun, show_column_names = T, column_names_gp = gpar(fontsize=annot_size), column_labels = c(expression(NLL),expression(NLL[theo])),
                   cluster_rows=F, cluster_columns=F, show_row_names = F, show_heatmap_legend = F, width=2*annot_width, height=h_height)
annot52 <- Heatmap(as.matrix(dat_annot[,11]),col = NLL_minus_NLLtheo_col_fun, show_column_names = T, column_names_gp = gpar(fontsize=annot_size), column_labels = c(expression(NLL~-~NLL[theo])),
                   cluster_rows=F, cluster_columns=F, show_row_names = F, show_heatmap_legend = F, width=annot_width, height=h_height, name="annot_52")

annot5 <- annot51 + annot52

lgd_Heatmap <- Legend(title=expression(log[2](FPKM)), col= Heatmap_col_fun, labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)
lgd_Fit_quality <- Legend(title = "Fit Quality", col = Qual_fit_col_fun, at=seq(0,0.2,by=0.05), labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)
lgd_NLL <- Legend(title = parse(text = "NLL~'&'~NLL[theo]"), col = NLL_col_fun, at=c(-100,-50,0,50), labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)
lgd_NLL_minus_NLLtheo <- Legend(title = parse(text = "NLL~-~NLL[theo]"), col = NLL_minus_NLLtheo_col_fun, at = c(0,50,100,150), labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)

lgd_Heatmap <- change_legend(lgd_Heatmap)
lgd_Fit_quality <- change_legend(lgd_Fit_quality)
lgd_NLL <- change_legend(lgd_NLL)
lgd_NLL_minus_NLLtheo <- change_legend(lgd_NLL_minus_NLLtheo)

lgd <- packLegend(lgd_Heatmap, lgd_Fit_quality, lgd_NLL, lgd_NLL_minus_NLLtheo)

gb_hca <- grid.grabExpr(draw(hca, column_title="Chromatin", column_title_gp = gpar(fontsize = title_size, fontface='bold', col = col_ca)))
gb_hnp <- grid.grabExpr(draw(hnp, column_title="Nucleoplasm", column_title_gp = gpar(fontsize = title_size, fontface='bold', col = col_np)))
gb_hcyto <- grid.grabExpr(draw(hcyto, column_title="Cytoplasm", column_title_gp = gpar(fontsize = title_size, fontface='bold', col = col_cyto)))
gb_annot5 <- grid.grabExpr(draw(annot5,gap=unit(0,"mm")))
gb_annot2 <- grid.grabExpr(draw(annot2))
gb_legend <- grid.grabExpr(draw(lgd))

vp_annot2 <- viewport(name="annot2", x = h_gap, y = 0, width=annot_width, just = c(0, 0))
vp_annot5 <- viewport(name="annot5", x = h_gap + vp_annot2$x + vp_annot2$width, y = (annot2@column_names_param$anno@height-annot5@ht_list$annot_52@column_names_param$anno@height)/2, width=3*annot_width, just = c(0, 0))
vp_ca <- viewport(name="ca", x = 2* h_gap + vp_annot5$x + vp_annot5$width, y = (annot2@column_names_param$anno@height-hca@column_names_param$anno@height+(subtitle_size+title_size)/12*unit(1,"lines"))/2, width= 2*h_width+h_gap, just = c("left", "bottom"))
vp_np <- viewport(name="np", x = 2* h_gap + vp_ca$x + vp_ca$width , y = vp_ca$y, width= 2*h_width + h_gap, just = c("left", "bottom"))
vp_cyto <- viewport(name="cyto", x = 2* h_gap + vp_np$x + vp_np$width, y = vp_ca$y, width= 2*h_width + h_gap, just = c("left", "bottom"))
vp_legend <- viewport(name="legend", x = 2 * h_gap + vp_cyto$x + vp_cyto$width, y=(gb_annot2$childrenvp$global$parent$height - lgd@grob$vp$height)/2 , width = lgd@grob$vp$width, just=c(0,0))

# svg(file = 'Figures/Figure2C_Heatmap_Fit_Naive_b3_v13.svg', width = 7, height = 7)

grid.newpage()
pushViewport(vp_annot2)
grid.draw(gb_annot2)
popViewport()

pushViewport(vp_annot5)
grid.draw(gb_annot5)
popViewport()

pushViewport(vp_ca)
grid.draw(gb_hca)
popViewport()

pushViewport(vp_np)
grid.draw(gb_hnp)
popViewport()

pushViewport(vp_cyto)
grid.draw(gb_hcyto)
popViewport()

pushViewport(vp_legend)
grid.draw(gb_legend)
popViewport()

# dev.off()

stop()

#### Figure S6A - Heatmap of data vs Fit LPA ####---------------------------------------------------------------
dat_2 <- as.matrix(dat_sim[['lpa']])
dat_2[dat_2 > 10] <- 10
dat_2[dat_2 < -5] <- -5

dat_3 <- dat_2[match(h$labels[reorder],paste0(gene_infos$ensembl_gene_id,'.',gene_infos$version)[match(row.names(dat_2),gene_infos$external_gene_name)]),]
dat_3 <- dat_3[!apply(dat_3,1,function(x){all(is.na(x))}),]

dat_annot <- annot[['lpa']][match(h$labels,paste0(gene_infos$ensembl_gene_id,'.',gene_infos$version)[match(row.names(annot[['lpa']]),gene_infos$external_gene_name)]),][reorder,]
dat_annot <- dat_annot[!apply(dat_annot,1,function(x){all(is.na(x))}),]

dat_annot2 <- data.frame(dat_annot[,1])
colnames(dat_annot2) <- "max(|autocor(error)|*mean(|error|)/range)"

# Color functions
Heatmap_col_fun <- colorRamp2(breaks=seq(min(dat_3), max(dat_3), length.out=101), colors_spectral)
Qual_fit_col_fun <- colorRamp2(breaks=seq(min(dat_annot[,1]), max(dat_annot[,1]),length.out=101),colorRampPalette(c("green3","gold2", "darkorange2","darkred"))(101))
NLL_col_fun <- colorRamp2(breaks=seq(min(dat_annot[,c(6,10)]), max(dat_annot[,c(6,10)]),length.out=101),colorRampPalette(c("white","darkred"))(101))
NLL_minus_NLLtheo_col_fun <- colorRamp2(breaks=seq(0, max(dat_annot[,c(6,10)])-min(dat_annot[,c(6,10)]),length.out=101),colorRampPalette(c("white","darkred"))(101))

ht_opt(TITLE_PADDING = unit(1, "mm"))
hca <- Heatmap(dat_3[,1:24], col = Heatmap_col_fun,
               cluster_rows=F, cluster_columns=F, 
               show_row_names=F, show_column_names=T, 
               column_names_gp = gpar(fontsize=annot_size),
               column_labels = structure(replace(colnames(dat_3)[1:24],!(colnames(dat_3)[1:24] %in% c(0,30,60,120)),""), names = colnames(dat_3)[1:24]), 
               column_split = rep(c("Data","Smoothed"),each=12), 
               column_title_gp = gpar(fontsize = subtitle_size), column_gap = h_gap, width=2*h_width+h_gap, height = h_height, show_heatmap_legend = F)
hnp <- Heatmap(dat_3[,25:48], col=Heatmap_col_fun, 
               cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, column_names_gp = gpar(fontsize=annot_size),
               column_labels = structure(replace(colnames(dat_3)[1:24],!(colnames(dat_3)[1:24] %in% c(0,30,60,120)),""), names = colnames(dat_3)[1:24]),
               column_split = rep(c("Data","Simulation"),each=12), column_title_gp = gpar(fontsize = subtitle_size), column_gap = h_gap, width=2*h_width+h_gap, height = h_height, show_heatmap_legend = F)
hcyto <- Heatmap(dat_3[,49:72], col=Heatmap_col_fun, 
                 cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=T, column_names_gp = gpar(fontsize=annot_size),
                 column_labels = structure(replace(colnames(dat_3)[1:24],!(colnames(dat_3)[1:24] %in% c(0,30,60,120)),""), names = colnames(dat_3)[1:24]), 
                 column_split = rep(c("Data","Simulation"),each=12), column_title_gp = gpar(fontsize = subtitle_size), column_gap = h_gap, width=2*h_width+h_gap, height = h_height, show_heatmap_legend = F)

annot2 <- Heatmap(matrix = as.matrix(dat_annot2), col = Qual_fit_col_fun, show_column_names = T, column_names_gp = gpar(fontsize=annot_size), column_labels = "Fit Quality",
                  cluster_rows=F, cluster_columns=F, show_row_names = F, show_heatmap_legend = F, width=annot_width, height=h_height)

annot51 <- Heatmap(as.matrix(dat_annot[,c(6,10)]),col = NLL_col_fun, show_column_names = T, column_names_gp = gpar(fontsize=annot_size), column_labels = c(expression(NLL),expression(NLL[theo])),
                   cluster_rows=F, cluster_columns=F, show_row_names = F, show_heatmap_legend = F, width=2*annot_width, height=h_height)
annot52 <- Heatmap(as.matrix(dat_annot[,11]),col = NLL_minus_NLLtheo_col_fun, show_column_names = T, column_names_gp = gpar(fontsize=annot_size), column_labels = c(expression(NLL~-~NLL[theo])),
                   cluster_rows=F, cluster_columns=F, show_row_names = F, show_heatmap_legend = F, width=annot_width, height=h_height, name="annot_52")

annot5 <- annot51 + annot52

lgd_Heatmap <- Legend(title=expression(log[2](FPKM)), col= Heatmap_col_fun, labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)
lgd_Fit_quality <- Legend(title = "Fit Quality", col = Qual_fit_col_fun, at=seq(0,0.2,by=0.05), labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)
lgd_NLL <- Legend(title = parse(text = "NLL~'&'~NLL[theo]"), col = NLL_col_fun, at=c(-100,-50,0,50), labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)
lgd_NLL_minus_NLLtheo <- Legend(title = parse(text = "NLL~-~NLL[theo]"), col = NLL_minus_NLLtheo_col_fun, at = c(0,50,100,150), labels_gp = gpar(fontsize=annot_size), title_gp=gpar(fontsize=annot_size), grid_width = legend_width)

lgd_Heatmap <- change_legend(lgd_Heatmap)
lgd_Fit_quality <- change_legend(lgd_Fit_quality)
lgd_NLL <- change_legend(lgd_NLL)
lgd_NLL_minus_NLLtheo <- change_legend(lgd_NLL_minus_NLLtheo)

lgd <- packLegend(lgd_Heatmap, lgd_Fit_quality, lgd_NLL, lgd_NLL_minus_NLLtheo)

gb_hca <- grid.grabExpr(draw(hca, column_title="Chromatin", column_title_gp = gpar(fontsize = title_size, fontface='bold', col = col_ca)))
gb_hnp <- grid.grabExpr(draw(hnp, column_title="Nucleoplasm", column_title_gp = gpar(fontsize = title_size, fontface='bold', col = col_np)))
gb_hcyto <- grid.grabExpr(draw(hcyto, column_title="Cytoplasm", column_title_gp = gpar(fontsize = title_size, fontface='bold', col = col_cyto)))
gb_annot5 <- grid.grabExpr(draw(annot5,gap=unit(0,"mm")))
gb_annot2 <- grid.grabExpr(draw(annot2))
gb_legend <- grid.grabExpr(draw(lgd))

vp_annot2 <- viewport(name="annot2", x = h_gap, y = 0, width=annot_width, just = c(0, 0))
vp_annot5 <- viewport(name="annot5", x = h_gap + vp_annot2$x + vp_annot2$width, y = (annot2@column_names_param$anno@height-annot5@ht_list$annot_52@column_names_param$anno@height)/2, width=3*annot_width, just = c(0, 0))
vp_ca <- viewport(name="ca", x = 2* h_gap + vp_annot5$x + vp_annot5$width, y = (annot2@column_names_param$anno@height-hca@column_names_param$anno@height+(subtitle_size+title_size)/12*unit(1,"lines"))/2, width= 2*h_width+h_gap, just = c("left", "bottom"))
vp_np <- viewport(name="np", x = 2* h_gap + vp_ca$x + vp_ca$width , y = vp_ca$y, width= 2*h_width + h_gap, just = c("left", "bottom"))
vp_cyto <- viewport(name="cyto", x = 2* h_gap + vp_np$x + vp_np$width, y = vp_ca$y, width= 2*h_width + h_gap, just = c("left", "bottom"))
vp_legend <- viewport(name="legend", x = 2 * h_gap + vp_cyto$x + vp_cyto$width, y=(gb_annot2$childrenvp$global$parent$height - lgd@grob$vp$height)/2 , width = lgd@grob$vp$width, just=c(0,0))

# svg(file = 'Figures/FigureS5A_Heatmap_Fit_LPA_b3_v13.svg', width = 7, height = 7)

grid.newpage()
pushViewport(vp_annot2)
grid.draw(gb_annot2)
popViewport()

pushViewport(vp_annot5)
grid.draw(gb_annot5)
popViewport()

pushViewport(vp_ca)
grid.draw(gb_hca)
popViewport()

pushViewport(vp_np)
grid.draw(gb_hnp)
popViewport()

pushViewport(vp_cyto)
grid.draw(gb_hcyto)
popViewport()

pushViewport(vp_legend)
grid.draw(gb_legend)
popViewport()

# dev.off()

stop()

#### Figure 2D - Fit Examples Malt1 and Egr1 ####---------------------------------------------------------------

## Line Plot examples:  Plot Malt 1 and Egr1 naive  
for (g in c("Egr1", "Malt1", 'Fos', 'Arl5c', 'Zhx2', 'Tnfaip3', 'Snx18', 'Nfkbie', 'Cd44', 'Tnf', 'Socs3')){
  dat_g <- data.frame(RPKM=c(), RPKM_ub = c(), RPKM_lb=c(), Type = c(), Time = c(), batch = c(), Cpt = c(), Condition=c())
  
  for (condition in c('naive', 'lpa')){
    if( g %in% names(data_all <- all_data[[condition]])){
      data_all <- all_data[[condition]][[g]]  
      
      for (b in c("rep1", "rep2", "rep3")){
        if(b == "rep3"){
          no_np <- TRUE
        }else{
          no_np <- FALSE
        }
        
        dat_g <- rbind(dat_g, data.frame(RPKM=data_all$caRNA[,b], RPKM_ub = NA, RPKM_lb=NA, Time = data_all$caRNA[,1], batch = b, Cpt = "Chromatin", Type="Experiment", Condition = condition))
        if(!no_np){
          dat_g <- rbind(dat_g, data.frame(RPKM=data_all$npRNA[,b], RPKM_ub = NA, RPKM_lb=NA, Time = data_all$npRNA[,1], batch = b, Cpt = "Nucleoplasm", Type="Experiment", Condition = condition))
        }
        dat_g <- rbind(dat_g, data.frame(RPKM=data_all$cytoRNA[,b], RPKM_ub = NA, RPKM_lb=NA, Time = data_all$cytoRNA[,1], batch = b, Cpt = "Cytoplasm", Type="Experiment", Condition = condition))
        
        if(g %in% row.names(all_best_param[[condition]][[b]])){
          par <- all_best_param[[condition]][[b]][g,]
          fit <- all_fit[[condition]][[b]][[g]]
          fit_ub <- all_fit_CI$UB[[condition]][[b]][[g]]
          fit_lb <- all_fit_CI$LB[[condition]][[b]][[g]]
          
          dat_g <- rbind(dat_g, data.frame(RPKM = fit$caRNA[,b], 
                                           RPKM_ub = fit_ub$caRNA[match(fit$caRNA[,"t"],fit_ub$caRNA[,"t"]),b], 
                                           RPKM_lb = fit_lb$caRNA[match(fit$caRNA[,"t"],fit_lb$caRNA[,"t"]),b], 
                                           Time = fit$caRNA[,"t"], batch = b, Cpt = "Chromatin", Type = "Fit", Condition = condition))
          if(!no_np){
            dat_g <- rbind(dat_g, data.frame(RPKM=fit$npRNA[,b], 
                                             RPKM_ub= fit_ub$npRNA[match(fit$npRNA[,"t"],fit_lb$npRNA[,"t"]),b], 
                                             RPKM_lb = fit_lb$npRNA[match(fit$npRNA[,"t"],fit_lb$npRNA[,"t"]),b], 
                                             Time = fit$npRNA[,"t"], batch = b, Cpt = "Nucleoplasm", Type = "Fit", Condition = condition))
          }
          
          dat_g <- rbind(dat_g, data.frame(RPKM=fit$cytoRNA[,b], 
                                           RPKM_ub= fit_ub$cytoRNA[match(fit$cytoRNA[,"t"],fit_lb$cytoRNA[,"t"]),b], 
                                           RPKM_lb = fit_lb$cytoRNA[match(fit$cytoRNA[,"t"],fit_lb$cytoRNA[,"t"]),b], 
                                           Time = fit$cytoRNA[,"t"], batch = b, Cpt = "Cytoplasm", Type = "Fit", Condition = condition))
        }
      }
    }
  }
  assign(paste0('dat_g_',g),dat_g)
}

dat_g <- dat_g_Malt1
gene_name <- 'Malt1'

dat_g <- dat_g_Egr1
gene_name <- 'Egr1'

for(gene_name in c("Egr1", "Malt1")){
  dat_g <- get(paste0("dat_g_", gene_name))
  dat_g$Cpt <- factor(dat_g$Cpt, levels = c("Chromatin", "Nucleoplasm", "Cytoplasm"), ordered = T)
  dat_g$batch <- factor(dat_g$batch, levels = c("rep1", "rep2", "rep3"), ordered = T)
  
  # svg("Figures/Figure2D_Malt1_v13.svg", width=4, height=3)
  # svg("Figures/Figure2D_Egr1_v13.svg", width=4, height=3)
  # plot Naive
  p_naive <- ggplot(dat_g[!is.na(dat_g$Type) & dat_g$Condition == "naive",], aes(x = Time, y = RPKM, group=batch)) + geom_line(data = dat_g[dat_g$Type == "Fit" & dat_g$Condition == "naive",], aes(col=batch), size=1, alpha=0.7) 
  p_naive <- p_naive + geom_point(data = dat_g[dat_g$Type == "Experiment" & dat_g$Condition == "naive",] , alpha = 0.5, aes(col=batch))
  p_naive <- p_naive + geom_errorbar(data = dat_g[!is.na(dat_g$Type) & dat_g$Condition == "naive" & !is.na(dat_g$RPKM_ub) & !is.na(dat_g$RPKM_lb),],mapping=aes(x=Time,ymin = RPKM_lb, ymax=RPKM_ub, color = batch, group=interaction(Type, batch)), size=0.5, alpha=0.3)
  p_naive <- p_naive + theme_bw() + theme(aspect.ratio=1) + scale_y_continuous(breaks = function(lims) {seq(floor(lims[1]), ceiling(lims[2]), 1)}) + ylab('log2(RPKM)')
  p_naive <- p_naive + labs(title=paste0(gene_name,  " - Naive"), color="Replicates") + scale_x_continuous(breaks = c(0,30,60,90,120))
  p_naive <- p_naive + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = rel(0.8)))
  p_naive <- p_naive + facet_wrap(~ Cpt)  + scale_color_manual(values = c("rep1" = "#FFd280FF", "rep2"="#BF7C00FF", "rep3"="#FFA540FF"), labels=c("Rep1","Rep2","Rep3"))
  p_naive <- p_naive + labs(y=expression(log[2](FPKM))) + theme(legend.position = "bottom")
  
  g1 <- ggplotGrob(p_naive)
  col <- c("Chromatin" = col_ca, "Nucleoplasm" = col_np, "Cytoplasm" = col_cyto)
  for(i in grep("strip-t",g1$layout$name)){
    Cpt <- g1$grobs[[i]]$grobs[[1]]$children[[grep("text", names(g1$grobs[[i]]$grobs[[1]]$children))]]$children[[1]]$label
    g1$grobs[[i]]$grobs[[1]]$children[[grep("rect", names(g1$grobs[[i]]$grobs[[1]]$children))]]$gp$fill <- col[Cpt]
  }
  grid.newpage()
  grid.draw(g1)
  # dev.off()
}
stop()

#### Figure 3 - Get Data for Plots 3AB, 4A, S5AB, 6B ####----------------------------------------------------------------

annot <- data.frame()

dat_ggplot <- data.frame(Gene_name = c(), Parameter_value = c(), Parameter_name=c(), 
                         CIL_0.95_Min = c(), CIL_0.95_Max = c(), CIL_0.95_Method = c(),
                         CIU_0.95_Min = c(), CIU_0.95_Max = c(), CIU_0.95_Method = c(),
                         Condition = c(), Batch=c())

for (condition in c('naive', 'lpa')){
  for(b in c('all', 'rep1', 'rep2', 'rep3')){
    for (g in names(all_param[[condition]][[b]])){
      if(!is.null(all_annot[[condition]][[b]][[g]])){
        par_best <- unlist(all_best_param[[condition]][[b]][g,])
        
        annot <- rbind(annot, t(all_annot[[condition]][[b]][[g]]$annot[c(1,6,10:11),]))
        row.names(annot)[nrow(annot)] <- paste0(g,'_',condition, '_', b)
        
        dat_ggplot <- rbind(dat_ggplot,
                            data.frame(Gene_name = g,
                                       Parameter_value = par_best, 
                                       Parameter_name = names(par_best),
                                       CIL_0.95_Min = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "lower_min"],
                                       CIL_0.95_Max = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "lower"],
                                       CIL_0.95_Method = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "method_lower"],
                                       CIU_0.95_Min = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "upper"],
                                       CIU_0.95_Max = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "upper_max"],
                                       CIU_0.95_Method = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "method_upper"],
                                       Condition = condition, Batch=b))
      }
    }
  }
}
colnames(annot) <- c( "max(|autocor(error)|*mean(|error|)/range)", 
                      "NLL", "NLL_theo", "NNL-NLL_theo")

#### Figure 3A plot number of identifiable #####
dat <- data.frame(batch = "rep1", 
                  do.call("rbind", lapply(1:length(all_param_CI$naive$rep1), function(i){cbind(Gene=names(all_param_CI$naive$rep1)[i],all_param_CI$naive$rep1[[i]]$CI[,c(1:4,7:8)])})))
dat <- rbind(dat, data.frame(batch = "rep2", 
                             do.call("rbind", lapply(1:length(all_param_CI$naive$rep2), function(i){cbind(Gene=names(all_param_CI$naive$rep2)[i],all_param_CI$naive$rep2[[i]]$CI[,c(1:4,7:8)])}))))

dat$batch <- factor(dat$batch, levels=c("rep1", "rep2"), labels=c("Rep1", "Rep2"))

head(dat)
dat$lower_min[is.na(dat$lower_min)] <- dat$lower[is.na(dat$lower_min)]
dat$upper_max[is.na(dat$upper_max)] <- dat$upper[is.na(dat$upper_max)]
head(dat)

dat$facet_label <- factor(dat$whichPar, levels=unique(dat$whichPar), labels=c("log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](k[1]*minute/k[2])", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])", "log[10](k[2]*minute/k[2])"))
dat$color_label <- factor(dat$whichPar, levels=unique(dat$whichPar), labels=c("k[1]*minute", "k[2]", "k[2]*minute", "k[deg]", "k[1]*minute/k[2]", "paste(k[1]*minute,k[2]*minute)/k[2]","k[2]*minute/k[2]"))

dat$CIrange_max <- dat$upper_max-dat$lower_min
dat$CIrange_min <- dat$upper-dat$lower

dat$identifiability <- NA
dat$identifiability[is.finite(dat$lower_min) & is.finite(dat$upper_max)] <- "Identifiable"
dat$identifiability[!is.finite(dat$lower_min) & !is.finite(dat$upper_max)] <- "Non-identifiable"
dat$identifiability[is.finite(dat$lower_min) & !is.finite(dat$upper_max)] <- "LB identifiable"
dat$identifiability[!is.finite(dat$lower_min) & is.finite(dat$upper_max)] <- "UB identifiable" 

dat$identifiability <- factor(dat$identifiability, levels=c("Non-identifiable", 'UB identifiable', "LB identifiable", "Identifiable"))
colors <- c("k1'"="red", "k2"="yellow", "k2'"="green", "kdeg"="blue", "k1'/k2"="black","k1'k2'/k2"="black", "k2'/k2"="black")
#### Identifiability ####-------------------------------------------------------
for (par in  c("k1'", "k2", "k2'", "kdeg", "k1'/k2", "k1'k2'/k2", "k2'/k2")){
  # svg(filename =  paste0('Figure3A_identifiability_',gsub('/','over',par),'.svg'), width=3)
  p <- ggplot(dat[dat$whichPar == par,]) 
  p <- p + geom_bar(aes(x=batch, fill=identifiability)) + facet_grid(~whichPar)
  p <- p + theme_bw() + labs(fill="Identifiability", y="Number", x="")
  p <- p + theme(aspect.ratio=2, strip.background = element_rect(color=colors[par], fill=alpha(colors[par],0.1), size=1.5, linetype="solid"))
  p <- p + scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
  print(p)
  # dev.off()
}

#### Figure 3A right plot of CIrange distribution of identifiable genes for k1', k2' ####

# svg(filename = 'Figure3B_distribution of CIrange.svg', width=4)
p <- ggplot(dat[is.finite(dat$CIrange_max) & dat$whichPar %in% c("k1'", "k2", "k2'", "kdeg"),]) 
bw <- bw.ucv(dat$CIrange_min)
p <- p + geom_density(mapping = aes(x=CIrange_min, fill=whichPar, color=whichPar, linetype=batch), alpha=0.1, bw="ucv")
p <- p + scale_color_manual(name="Parameter", values=c("k1'"="red", "k2"="yellow", "k2'" = "green", "kdeg"="blue"))
p <- p + scale_fill_manual(name="Parameter", values=c("k1'"="red", "k2"="yellow", "k2'" = "green", "kdeg"="blue"))
p <- p + scale_linetype_manual(values=c("Rep1"=1, "Rep2" = 2), labels = c("Rep1", "Rep2"))
p <- p + theme_bw() + theme(aspect.ratio=2/3) + labs(x=expression(log[10]("95% CIU/CIL")), linetype="Replicate")
p
#dev.off()


# Data
annot <- data.frame()

dat_ggplot <- data.frame(Gene_name = c(), Parameter_value = c(), Parameter_name=c(), 
                         CIL_0.95_Min = c(), CIL_0.95_Max = c(), CIL_0.95_Method = c(),
                         CIU_0.95_Min = c(), CIU_0.95_Max = c(), CIU_0.95_Method = c(),
                         Condition = c(), Batch=c())

for (condition in c('naive', 'lpa')){
  for(b in c('all', 'rep1', 'rep2', 'rep3')){
    for (g in names(all_param[[condition]][[b]])){
      if(!is.null(all_annot[[condition]][[b]][[g]])){
          
        par_best <- unlist(all_best_param[[condition]][[b]][g,])
        
        annot <- rbind(annot, t(all_annot[[condition]][[b]][[g]]$annot[c(1,6,10:11),]))
        row.names(annot)[nrow(annot)] <- paste0(g,'_',condition, '_', b)
        
        dat_ggplot <- rbind(dat_ggplot,
                            data.frame(Gene_name = g,
                                       Parameter_value = par_best, 
                                       Parameter_name = names(par_best),
                                       CIL_0.95_Min = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "lower_min"],
                                       CIL_0.95_Max = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "lower"],
                                       CIL_0.95_Method = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "method_lower"],
                                       CIU_0.95_Min = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "upper"],
                                       CIU_0.95_Max = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "upper_max"],
                                       CIU_0.95_Method = all_param_CI[[condition]][[b]][[g]]$CI[names(par_best), "method_upper"],
                                       Condition = condition, Batch=b))
      }
    }
  }
}
colnames(annot) <- c( "max(|autocor(error)|*mean(|error|)/range)", 
                      "NLL", "NLL_theo", "NNL-NLL_theo")



###### Fig 3B Distribution of Parameter for identifiable ones ####
palette <- colorRampPalette(colors = c("green", "yellow", "orange", "red"))(100)
make_palette <- function(value, palette){
  res <- palette[1+(value-min(value, na.rm = T))/(max(value, na.rm = T) -min(value, na.rm = T))*99]
}

dat_ggplot$facet_label <- factor(dat_ggplot$Parameter_name, levels=unique(dat_ggplot$Parameter_name), labels=c("spar","sigma[t]","sigma[b]","ca0[b1]", "log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](k[1]*minute/k[2])", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])","log[10](k[2]*minute/k[2])", "log[10](paste(k[1]*minute,k[2]*minute))", "log[10](max(k[2],k[deg]))", "log[10](min(k[2],k[deg]))"))
dat_ggplot$color_label <- factor(dat_ggplot$Parameter_name, levels=unique(dat_ggplot$Parameter_name), labels=c("spar","sigma[t]","sigma[b]","ca0[b1]", "k[1]*minute", "k[2]", "k[2]*minute", "k[deg]", "k[1]*minute/k[2]", "paste(k[1]*minute,k[2]*minute)/k[2]", "paste(k[2]*minute/k[2]", "paste(k[1]*minute,k[2]*minute)", "max(k[2],k[deg])", "min(k[2],k[deg])"))

# svg(filename = 'Fig3B_v13_ident.svg')
p <- ggplot(dat_ggplot[dat_ggplot$Parameter_name %in% c("k1'", "k2", "k2'", "kdeg") & dat_ggplot$Condition == "naive" & dat_ggplot$Batch %in% c("rep1", "rep2") & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max) ,]) 
p <- p + geom_density(mapping = aes(x=Parameter_value, fill=Parameter_name, color=Parameter_name, linetype=Batch), alpha=0.1)
p <- p + scale_color_manual(name="Parameter", values=c("red", "yellow","green", "blue"))
p <- p + scale_fill_manual(name="Parameter", values=c("red", "yellow", "green", "blue"))
p <- p + scale_linetype_manual(values=c("rep1"=1, "rep2" = 2), labels = c("Rep1", "Rep2"))
p <- p + theme_bw() + theme(aspect.ratio=1/4) + labs(x=expression(log[10]("Parameter value")), linetype="Replicate")
p
# dev.off()
stop()

# test similarity of distribution
ks.test(dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "k1'" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)],
        dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "k1'" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)])
ks.test(dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "k2'" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)],
        dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "k2'" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)])
ks.test(dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "k2" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)],
        dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "k2" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)])
ks.test(dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "kdeg" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)],
        dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "kdeg" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)])

# test JSD: 
x_rep1 <- dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "kdeg" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)]
x_rep2 <- dat_ggplot$Parameter_value[dat_ggplot$Parameter_name == "kdeg" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max)]
ecdf_rep1 <- ecdf(x_rep1)
ecdf_rep2 <- ecdf(x_rep2)
philentropy::distance(x = rbind(ecdf_rep1(sort(c(x_rep1,x_rep2))),ecdf_rep2(sort(c(x_rep1,x_rep2)))), method = "jensen-shannon")
philentropy::distance(x = rbind(ecdf_rep1(seq(min(x_rep1,x_rep2),max(x_rep1,x_rep2),by=0.0001)),ecdf_rep2(seq(min(x_rep1,x_rep2),max(x_rep1,x_rep2),by=0.0001))), method = "jensen-shannon")
philentropy::JSD(x = rbind(ecdf_rep1(seq(min(x_rep1,x_rep2),max(x_rep1,x_rep2),by=0.00001)),ecdf_rep2(seq(min(x_rep1,x_rep2),max(x_rep1,x_rep2),by=0.00001))))

epsilon <- 0.5 * min(c(diff(sort(x_rep1)),diff(sort(x_rep2)))) 
KLrep1 = (1/length(x_rep1))*sum(log((ecdf_rep1(x_rep1) - ecdf_rep1(x_rep1-epsilon))/(ecdf_rep2(x_rep1) - ecdf_rep2(x_rep1-epsilon)))) # need to linearly interpolate ecdf
KLrep2 = (1/length(x_rep2))*sum(log((ecdf_rep1(x_rep2) - ecdf_rep1(x_rep2-epsilon))/(ecdf_rep2(x_rep1) - ecdf_rep2(x_rep2-epsilon)))) # need to linearly interpolate ecdf
JSD <- 0.5 * (KLrep1 + KLrep2)

dat_ggplot_rep1 <- dat_ggplot[dat_ggplot$Batch == 'rep1',]
dat_ggplot_rep2 <- dat_ggplot[dat_ggplot$Batch == 'rep2',]

dat_ggplot2 <- dat_ggplot_rep1

# Add column for b3 value and CI
colnames(dat_ggplot2) <- gsub(pattern = "Parameter_value" , replacement = "Parameter_value_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$Parameter_value_rep2 <- dat_ggplot_rep2$Parameter_value[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                      paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIL_0.95_Min" , replacement = "CIL_0.95_Min_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIL_0.95_Min_rep2 <- dat_ggplot_rep2$CIL_0.95_Min[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIL_0.95_Max" , replacement = "CIL_0.95_Max_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIL_0.95_Max_rep2 <- dat_ggplot_rep2$CIL_0.95_Max[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIU_0.95_Min" , replacement = "CIU_0.95_Min_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIU_0.95_Min_rep2 <- dat_ggplot_rep2$CIU_0.95_Min[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIU_0.95_Max" , replacement = "CIU_0.95_Max_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIU_0.95_Max_rep2 <- dat_ggplot_rep2$CIU_0.95_Max[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

# add color for annot
dat_ggplot2$color_rep1 <- annot[paste(dat_ggplot2$Gene_name, dat_ggplot2$Condition, 'rep1', sep="_"),1]
dat_ggplot2$color_rep2 <- annot[paste(dat_ggplot2$Gene_name, dat_ggplot2$Condition, 'rep2', sep="_"),1]
dat_ggplot2$color_max <- apply(dat_ggplot2[,paste("color", c("rep1", "rep2"), sep="_")],1, max)

gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

#### Figure 3C plot reproducibility of identifiable genes foe k1', k2' ####
make_replicate_plot2 <- function(data, condition,parameter, xylim=NULL, color_facet=NULL, font_size=12, legend = "none"){
  
  p <- ggplot(data[data$Condition == condition & data$Parameter_name == parameter & apply(data[,c("CIL_0.95_Min_rep1", "CIU_0.95_Max_rep1", "CIL_0.95_Min_rep2", "CIU_0.95_Max_rep2")],1,function(x){all(is.na(x))}),],aes(x=Parameter_value_rep1,y=Parameter_value_rep2))
  p <- p + geom_point(mapping=aes(col=color_max), size=2) + facet_wrap(.~facet_label, scale='free', drop = T, labeller = "label_parsed")
  p <- p + geom_segment(mapping = aes(x=CIL_0.95_Min_rep1,y=Parameter_value_rep2,xend=CIL_0.95_Max_rep1,yend = Parameter_value_rep2, col=color_rep1), alpha=0.3, linetype=2)
  p <- p + geom_segment(mapping = aes(x=CIU_0.95_Min_rep1,y=Parameter_value_rep2,xend=CIU_0.95_Max_rep1,yend = Parameter_value_rep2, col=color_rep1), alpha=0.3, linetype=2)
  p <- p + geom_segment(mapping = aes(x=CIL_0.95_Max_rep1,y=Parameter_value_rep2,xend=CIU_0.95_Min_rep1,yend = Parameter_value_rep2, col=color_rep1), alpha=0.3)
  p <- p + geom_segment(mapping = aes(y=CIL_0.95_Min_rep2,x=Parameter_value_rep1,yend=CIL_0.95_Max_rep2,xend = Parameter_value_rep1, col=color_rep2), alpha=0.3, linetype = 2)
  p <- p + geom_segment(mapping = aes(y=CIU_0.95_Min_rep2,x=Parameter_value_rep1,yend=CIU_0.95_Max_rep2,xend = Parameter_value_rep1, col=color_rep2), alpha=0.3, linetype = 2)
  p <- p + geom_segment(mapping = aes(y=CIL_0.95_Max_rep2,x=Parameter_value_rep1,yend=CIU_0.95_Min_rep2,xend = Parameter_value_rep1, col=color_rep2), alpha=0.3)
  p <- p + theme_bw(base_size = font_size) + theme(aspect.ratio=1) + scale_color_gradientn(colours = c("green3", "gold2", "darkorange2", "darkred"), limits=range(c(data$color_rep1,data$color_rep2)))+labs(col='Fit Quality')
  p <- p + geom_abline(slope=1, intercept=0, linetype=1, col='black', alpha=0.5) + geom_abline(slope=1, intercept=log10(2), linetype=2, col='red', alpha=0.5) + geom_abline(slope=1, intercept=-log10(2), linetype=2, col='red', alpha=0.5)
  p <- p + labs(x="Replicate 1", y="Replicate 2")
  p <- p + theme(strip.text = element_text(margin=margin(0.2,0,0.2,0,unit = "line")))
  
  if(!is.null(xylim)){
    p <- p + coord_cartesian(xlim=xylim, ylim=xylim)
  }
  if(!is.null(color_facet)){
    p <- p + theme(strip.background = element_rect(color=color_facet, fill=alpha(color_facet,0.1), size=1.5, linetype="solid"))
  }
  
  if(legend == "none"){
    p <- p + theme(legend.position = "none")
  }else if(legend == "only"){
    p <- gglegend(p)
  }
  
  return(p)
}

# svg(filename = 'Figure3C_reproducibility_identifiable_k1p.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k1'", color_facet = "red", legend = "none", font_size=24, xylim = c(-1.75,0.75)) # xylim=c(-1.5,1.1)
p 
# dev.off()

# svg(filename = 'Figure3C_reproducibility_identifiable_k2.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k2", color_facet = "yellow", legend = "none", font_size=24, xylim=c(-3,0.5)) # xylim=c(-1.5,1.1)
p 
#dev.off()

#svg(filename = "Figure3C_reproducibility_identifiable_k2'.svg", width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k2'", color_facet = "green", xylim=c(-3,0), legend = "none", font_size=24)
p
#dev.off()

#svg(filename = 'Fig3C_reproducibility_identifiable_kdeg.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "kdeg", color_facet = "blue", xylim=c(-4,0), legend = "none", font_size=24)
p
#dev.off()

#svg(filename = 'Fig3C_reproducibikity_legend.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "kdeg", color_facet = "blue", xylim=c(-4,0), legend = "only", font_size=24)
grid.draw(p)
#dev.off()

#svg(filename = 'Fig3D_k1prime_over_k2_identifiable.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k1'/k2", xylim=c(-0.2,1.2), legend = "none", font_size=24)
p
#dev.off()

# svg(filename = 'Fig3D_k1primek2prime_over_k2_identifiable.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k1'k2'/k2", xylim=c(-2.7,0.5), legend = "none", font_size=24)
p
# dev.off()

#svg(filename = 'Fig3D_k2prime_over_k2_identifiable.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k2'/k2", xylim=c(-2,1), legend = "none", font_size=24)
p
#dev.off()
stop()

#### Suplementary - Fit quality ####--------------------------------------------
dat_ggplot$color <- annot[paste(dat_ggplot$Gene_name, dat_ggplot$Condition, dat_ggplot$Batch, sep="_"),1]

# svg(filename = "Fit_quality_histogram.svg", width=5, height=5)
p <- ggplot(dat_ggplot[dat_ggplot$Condition == "naive" & dat_ggplot$Batch %in% c("rep1", "rep2") & dat_ggplot$Parameter_name == "spar",]) + geom_density(aes(x=color, linetype=Batch), fill="grey", alpha=0.3)
p <- p + scale_linetype_discrete(name = "Replicates", labels = c("rep1"="Rep1", "rep2"="Rep2"))
for (g in c('Tnfaip2', 'Il10', 'Cpd', 'Cd44')){
  p <- p + geom_segment(data = dat_ggplot[dat_ggplot$Gene_name == g & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & dat_ggplot$Parameter_name == "spar",], aes(x = color, xend=color, y=0, yend=2.5, linetype=Batch, col=color))
  p <- p + geom_text(data = dat_ggplot[dat_ggplot$Gene_name == g & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & dat_ggplot$Parameter_name == "spar",], aes(x = color, label=Gene_name, col=color), y = 5, fontface="bold", size=rel(4))
}
p <- p + scale_color_gradientn(colours = c("green3", "gold2", "darkorange2", "darkred"), limits=range(dat_ggplot$color[dat_ggplot$Batch %in% c("rep1","rep2")]))
p <- p + theme_bw() + theme(aspect.ratio=0.5) + labs(x="Fit Quality", y="Density", col="Fit Quality")
p
# dev.off()

dat_ggplot[dat_ggplot$Gene_name %in% c("Tnfaip2", "Il10", "Cpd", "Cd44") & dat_ggplot$Parameter_name=="spar" & dat_ggplot$Condition == "naive" & dat_ggplot$Batch %in% c("rep1", "rep2"),]
sum(dat_ggplot$Gene_name[dat_ggplot$color >= 0.06 & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & dat_ggplot$Parameter_name == "spar"] %in% dat_ggplot$Gene_name[dat_ggplot$color >= 0.06 & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & dat_ggplot$Parameter_name == "spar"])
sum(dat_ggplot$Gene_name[dat_ggplot$color >= 0.1 & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep1" & dat_ggplot$Parameter_name == "spar"] %in% dat_ggplot$Gene_name[dat_ggplot$color >= 0.1 & dat_ggplot$Condition == "naive" & dat_ggplot$Batch == "rep2" & dat_ggplot$Parameter_name == "spar"])

sum(dat_ggplot$Gene_name[dat_ggplot$color >= 0.06 & dat_ggplot$Condition == "lpa" & dat_ggplot$Batch == "rep1" & dat_ggplot$Parameter_name == "spar"] %in% dat_ggplot$Gene_name[dat_ggplot$color >= 0.06 & dat_ggplot$Condition == "lpa" & dat_ggplot$Batch == "rep2" & dat_ggplot$Parameter_name == "spar"])
sum(dat_ggplot$Gene_name[dat_ggplot$color >= 0.1 & dat_ggplot$Condition == "lpa" & dat_ggplot$Batch == "rep1" & dat_ggplot$Parameter_name == "spar"] %in% dat_ggplot$Gene_name[dat_ggplot$color >= 0.1 & dat_ggplot$Condition == "lpa" & dat_ggplot$Batch == "rep2" & dat_ggplot$Parameter_name == "spar"])


for (g in c("Ccl5", "Il10", 'Acp5', 'Cybb', 'Cd44', 'Anxa5', 'Id3', 'Mtmr12', 'Cpd', 'Tnfaip2')){
  dat_g <- data.frame(RPKM=c(), RPKM_ub = c(), RPKM_lb=c(), Type = c(), Time = c(), batch = c(), Cpt = c(), Condition=c())
  
  for (condition in c('naive', "lpa")){
    if( g %in% names(data_all <- all_data[[condition]])){
      data_all <- all_data[[condition]][[g]]  
      
      for (b in c("rep1", "rep2", "rep3")){
        if(b == "rep3"){
          no_np <- TRUE
        }else{
          no_np <- FALSE
        }
        
        dat_g <- rbind(dat_g, data.frame(RPKM=data_all$caRNA[,b], RPKM_ub = NA, RPKM_lb=NA, Time = data_all$caRNA[,1], batch = b, Cpt = "Chromatin", Type="Experiment", Condition = condition))
        if(!no_np){
          dat_g <- rbind(dat_g, data.frame(RPKM=data_all$npRNA[,b], RPKM_ub = NA, RPKM_lb=NA, Time = data_all$npRNA[,1], batch = b, Cpt = "Nucleoplasm", Type="Experiment", Condition = condition))
        }
        dat_g <- rbind(dat_g, data.frame(RPKM=data_all$cytoRNA[,b], RPKM_ub = NA, RPKM_lb=NA, Time = data_all$cytoRNA[,1], batch = b, Cpt = "Cytoplasm", Type="Experiment", Condition = condition))
        
        if(g %in% row.names(all_best_param[[condition]][[b]])){
          par <- all_best_param[[condition]][[b]][g,]
          fit <- all_fit[[condition]][[b]][[g]]
          fit_ub <- all_fit_CI$UB[[condition]][[b]][[g]]
          fit_lb <- all_fit_CI$LB[[condition]][[b]][[g]]
          
          dat_g <- rbind(dat_g, data.frame(RPKM = fit$caRNA[,b], 
                                           RPKM_ub = fit_ub$caRNA[match(fit$caRNA[,"t"],fit_ub$caRNA[,"t"]),b], 
                                           RPKM_lb = fit_lb$caRNA[match(fit$caRNA[,"t"],fit_lb$caRNA[,"t"]),b], 
                                           Time = fit$caRNA[,"t"], batch = b, Cpt = "Chromatin", Type = "Fit", Condition = condition))
          if(!no_np){
            dat_g <- rbind(dat_g, data.frame(RPKM=fit$npRNA[,b], 
                                             RPKM_ub= fit_ub$npRNA[match(fit$npRNA[,"t"],fit_lb$npRNA[,"t"]),b], 
                                             RPKM_lb = fit_lb$npRNA[match(fit$npRNA[,"t"],fit_lb$npRNA[,"t"]),b], 
                                             Time = fit$npRNA[,"t"], batch = b, Cpt = "Nucleoplasm", Type = "Fit", Condition = condition))
          }
          
          dat_g <- rbind(dat_g, data.frame(RPKM=fit$cytoRNA[,b], 
                                           RPKM_ub= fit_ub$cytoRNA[match(fit$cytoRNA[,"t"],fit_lb$cytoRNA[,"t"]),b], 
                                           RPKM_lb = fit_lb$cytoRNA[match(fit$cytoRNA[,"t"],fit_lb$cytoRNA[,"t"]),b], 
                                           Time = fit$cytoRNA[,"t"], batch = b, Cpt = "Cytoplasm", Type = "Fit", Condition = condition))
        }
      }
    }
  }
  assign(paste0('dat_g_',g),dat_g)
}


plot_gene_fit <- function(gene_name, condition="naive", batch){
  dat_g <- get(paste0("dat_g_",gene_name))
  dat_g$Cpt <- factor(dat_g$Cpt, levels = c("Chromatin", "Nucleoplasm", "Cytoplasm"), ordered = T)
  
  p <- ggplot(dat_g[!is.na(dat_g$Type) & dat_g$Condition == condition & dat_g$batch %in% batch,], aes(x = Time, y = RPKM, group=batch)) + geom_line(data = dat_g[dat_g$Type == "Fit" & dat_g$Condition == condition & dat_g$batch %in% batch,], aes(col=batch), size=1, alpha=0.7) 
  p <- p + geom_point(data = dat_g[dat_g$Type == "Experiment" & dat_g$Condition == condition & dat_g$batch %in% batch,] , alpha = 0.5, aes(col=batch))
  # p <- p + geom_errorbar(data = dat_g[!is.na(dat_g$Type) & dat_g$Condition == condition & dat_g$batch %in% batch & !is.na(dat_g$RPKM_ub) & !is.na(dat_g$RPKM_lb),],mapping=aes(x=Time,ymin = RPKM_lb, ymax=RPKM_ub, color = batch, group=interaction(Type, batch)), size=0.5, alpha=0.3)
  p <- p + theme_bw() + theme(aspect.ratio=1) + scale_y_continuous(breaks = function(lims) {seq(floor(lims[1]), ceiling(lims[2]), 1)}) + ylab('log2(RPKM)')
  
  
  if (length(batch)==1){  
    fit <- unique(dat_ggplot$color[dat_ggplot$Gene_name==gene_name & dat_ggplot$Batch %in% batch & dat_ggplot$Condition == condition])
    p <- p + labs(title=paste0(gene_name,  " - ", toupper(substr(condition,1,1)), substring(condition,2)), subtitle = paste0("Fit quality: ", signif(fit,2)) , color="Replicates") + scale_x_continuous(breaks = c(0,30,60,90,120))
  }else{
    p <- p + labs(title=paste0(gene_name,  " - ", toupper(substr(condition,1,1)), substring(condition,2)), color="Replicates") + scale_x_continuous(breaks = c(0,30,60,90,120))
  }
  p <- p + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = rel(0.8)))
  p <- p + facet_wrap(~ Cpt)  + scale_color_manual(values = c("rep1" = "#FFd280FF", "rep3"="#FFA540FF", "rep2"="#BF7C00FF"), labels=c("rep1" = "Rep1", "rep3"="Rep3","rep2"="Rep2"))
  p <- p + labs(y=expression(log[2](FPKM))) + theme(legend.position = "bottom")
  
  g1 <- ggplotGrob(p)
  col <- c("Chromatin" = col_ca, "Nucleoplasm" = col_np, "Cytoplasm" = col_cyto)
  for(i in grep("strip-t",g1$layout$name)){
    Cpt <- g1$grobs[[i]]$grobs[[1]]$children[[grep("text", names(g1$grobs[[i]]$grobs[[1]]$children))]]$children[[1]]$label
    g1$grobs[[i]]$grobs[[1]]$children[[grep("rect", names(g1$grobs[[i]]$grobs[[1]]$children))]]$gp$fill <- col[Cpt]
  }
  return(g1)
}
# dev.off()


# Good gene example 
unique(unname(apply(dat_ggplot[dat_ggplot$Condition == "naive" &dat_ggplot$color < 0.015 & dat_ggplot$Batch == "rep2", c("Gene_name","Batch","color")], 1, function(x){paste(x,collapse=".")})))
# svg(filename = "Fit_quality_Tnfaip2.svg", width=5, height=5)
p <- plot_gene_fit("Tnfaip2", condition="naive",batch="rep2")
grid.newpage()
grid.draw(p)
# dev.off()

# Good gene example 
unique(unname(apply(dat_ggplot[dat_ggplot$Condition == "naive" & dat_ggplot$color > 0.035 & dat_ggplot$color < 0.055 & dat_ggplot$Batch == "rep2", c("Gene_name","Batch","color")], 1, function(x){paste(x,collapse=".")})))
# svg(filename = "Fit_quality_Il10.svg", width=5, height=5)
p <- plot_gene_fit("Il10", condition="naive",batch="rep2")
grid.newpage()
grid.draw(p)
# dev.off()

# Medium example
unique(unname(apply(dat_ggplot[dat_ggplot$Condition == "naive" &dat_ggplot$color > 0.09 & dat_ggplot$color < 0.15 & dat_ggplot$Batch == "rep2", c("Gene_name","Batch","color")], 1, function(x){paste(x,collapse=".")})))
# svg(filename = "Fit_quality_Cpd.svg", width=5, height=5)
p <- plot_gene_fit("Cpd", condition="naive",batch="rep2")
grid.newpage()
grid.draw(p)
# dev.off()

# bad gene example
unique(unname(apply(dat_ggplot[dat_ggplot$Condition == "naive" & dat_ggplot$color > 0.15 & dat_ggplot$Batch == "rep2", c("Gene_name","Batch","color")], 1, function(x){paste(x,collapse=".")})))
# svg(filename = "Fit_quality_Cd44.svg", width=5, height=5)
p <- plot_gene_fit("Cd44", condition="naive",batch="rep2")
grid.newpage()
grid.draw(p)
# dev.off()

###### Fig 3B Distribution of Parameter for identifiable ones ####
palette <- colorRampPalette(colors = c("green", "yellow", "orange", "red"))(100)
make_palette <- function(value, palette){
  res <- palette[1+(value-min(value, na.rm = T))/(max(value, na.rm = T) -min(value, na.rm = T))*99]
}

dat_ggplot$facet_label <- factor(dat_ggplot$Parameter_name, levels=unique(dat_ggplot$Parameter_name), labels=c("spar","sigma[t]","sigma[b]","ca0[rep1]", "log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](k[1]*minute/k[2])", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])","log[10](k[2]*minute/k[2])", "log[10](paste(k[1]*minute,k[2]*minute))", "log[10](max(k[2],k[deg]))", "log[10](min(k[2],k[deg]))"))
dat_ggplot$color_label <- factor(dat_ggplot$Parameter_name, levels=unique(dat_ggplot$Parameter_name), labels=c("spar","sigma[t]","sigma[b]","ca0[rep1]", "k[1]*minute", "k[2]", "k[2]*minute", "k[deg]", "k[1]*minute/k[2]", "paste(k[1]*minute,k[2]*minute)/k[2]", "paste(k[2]*minute/k[2]", "paste(k[1]*minute,k[2]*minute)", "max(k[2],k[deg])", "min(k[2],k[deg])"))

# svg(filename = 'Fig3B_v13_ident.svg')
p <- ggplot(dat_ggplot[dat_ggplot$Parameter_name %in% c("k1'", "k2", "k2'", "kdeg") & dat_ggplot$Condition == "naive" & dat_ggplot$Batch %in% c("rep1", "rep2") & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max) ,]) 
p <- p + geom_density(mapping = aes(x=Parameter_value, fill=Parameter_name, color=Parameter_name, linetype=Batch), alpha=0.1)
p <- p + scale_color_manual(name="Parameter", values=c("red", "yellow","green", "blue"))
p <- p + scale_fill_manual(name="Parameter", values=c("red", "yellow", "green", "blue"))
p <- p + scale_linetype_manual(values=c("rep1"=1, "rep2" = 2), labels = c("Rep1", "Rep2"))
p <- p + theme_bw() + theme(aspect.ratio=1/4) + labs(x=expression(log[10]("Parameter value")), linetype="Replicate")
p
# dev.off()
stop()

dat_ggplot_rep1 <- dat_ggplot[dat_ggplot$Batch == 'rep1',]
dat_ggplot_rep2 <- dat_ggplot[dat_ggplot$Batch == 'rep2',]

dat_ggplot2 <- dat_ggplot_rep1

# Add column for rep2 value and CI
colnames(dat_ggplot2) <- gsub(pattern = "Parameter_value" , replacement = "Parameter_value_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$Parameter_value_rep2 <- dat_ggplot_rep2$Parameter_value[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                      paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIL_0.95_Min" , replacement = "CIL_0.95_Min_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIL_0.95_Min_rep2 <- dat_ggplot_rep2$CIL_0.95_Min[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIL_0.95_Max" , replacement = "CIL_0.95_Max_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIL_0.95_Max_rep2 <- dat_ggplot_rep2$CIL_0.95_Max[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIU_0.95_Min" , replacement = "CIU_0.95_Min_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIU_0.95_Min_rep2 <- dat_ggplot_rep2$CIU_0.95_Min[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

colnames(dat_ggplot2) <- gsub(pattern = "CIU_0.95_Max" , replacement = "CIU_0.95_Max_rep1", x = colnames(dat_ggplot2))
dat_ggplot2$CIU_0.95_Max_rep2 <- dat_ggplot_rep2$CIU_0.95_Max[match(paste(dat_ggplot_rep1$Gene_name,dat_ggplot_rep1$Parameter_name,dat_ggplot_rep1$Condition,sep = '.'),
                                                                paste(dat_ggplot_rep2$Gene_name,dat_ggplot_rep2$Parameter_name,dat_ggplot_rep2$Condition,sep = '.'))]

# add color for annot
dat_ggplot2$color_rep1 <- annot[paste(dat_ggplot2$Gene_name, dat_ggplot2$Condition, 'rep1', sep="_"),1]
dat_ggplot2$color_rep2 <- annot[paste(dat_ggplot2$Gene_name, dat_ggplot2$Condition, 'rep2', sep="_"),1]
dat_ggplot2$color_max <- apply(dat_ggplot2[,paste("color", c("rep1", "rep2"), sep="_")],1, max)

gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}

#### Figure 3C plot reproducibility of identifiable genes foe k1', k2' ####
make_replicate_plot2 <- function(data, condition,parameter, xylim=NULL, color_facet=NULL, font_size=12, legend = "none"){
  
  p <- ggplot(data[data$Condition == condition & data$Parameter_name == parameter & apply(data[,c("CIL_0.95_Min_rep1", "CIU_0.95_Max_rep1", "CIL_0.95_Min_rep2", "CIU_0.95_Max_rep2")],1,function(x){all(is.na(x))}),],aes(x=Parameter_value_rep1,y=Parameter_value_rep2))
  p <- p + geom_point(mapping=aes(col=color_max), size=2) + facet_wrap(.~facet_label, scale='free', drop = T, labeller = "label_parsed")
  p <- p + geom_segment(mapping = aes(x=CIL_0.95_Min_rep1,y=Parameter_value_rep2,xend=CIL_0.95_Max_rep1,yend = Parameter_value_rep2, col=color_rep1), alpha=0.3, linetype=2)
  p <- p + geom_segment(mapping = aes(x=CIU_0.95_Min_rep1,y=Parameter_value_rep2,xend=CIU_0.95_Max_rep1,yend = Parameter_value_rep2, col=color_rep1), alpha=0.3, linetype=2)
  p <- p + geom_segment(mapping = aes(x=CIL_0.95_Max_rep1,y=Parameter_value_rep2,xend=CIU_0.95_Min_rep1,yend = Parameter_value_rep2, col=color_rep1), alpha=0.3)
  p <- p + geom_segment(mapping = aes(y=CIL_0.95_Min_rep2,x=Parameter_value_rep1,yend=CIL_0.95_Max_rep2,xend = Parameter_value_rep1, col=color_rep2), alpha=0.3, linetype = 2)
  p <- p + geom_segment(mapping = aes(y=CIU_0.95_Min_rep2,x=Parameter_value_rep1,yend=CIU_0.95_Max_rep2,xend = Parameter_value_rep1, col=color_rep2), alpha=0.3, linetype = 2)
  p <- p + geom_segment(mapping = aes(y=CIL_0.95_Max_rep2,x=Parameter_value_rep1,yend=CIU_0.95_Min_rep2,xend = Parameter_value_rep1, col=color_rep2), alpha=0.3)
  p <- p + theme_bw(base_size = font_size) + theme(aspect.ratio=1) + scale_color_gradientn(colours = c("green3", "gold2", "darkorange2", "darkred"), limits=range(c(data$color_rep1,data$color_rep2)))+labs(col='Fit Quality')
  p <- p + geom_abline(slope=1, intercept=0, linetype=1, col='black', alpha=0.5) + geom_abline(slope=1, intercept=log10(2), linetype=2, col='red', alpha=0.5) + geom_abline(slope=1, intercept=-log10(2), linetype=2, col='red', alpha=0.5)
  p <- p + labs(x="Replicate 1", y="Replicate 2")
  p <- p + theme(strip.text = element_text(margin=margin(0.2,0,0.2,0,unit = "line")))
  
  if(!is.null(xylim)){
    p <- p + coord_cartesian(xlim=xylim, ylim=xylim)
  }
  if(!is.null(color_facet)){
    p <- p + theme(strip.background = element_rect(color=color_facet, fill=alpha(color_facet,0.1), size=1.5, linetype="solid"))
  }
  
  if(legend == "none"){
    p <- p + theme(legend.position = "none")
  }else if(legend == "only"){
    p <- gglegend(p)
  }
  
  return(p)
}

# svg(filename = 'Figure3C_reproducibility_identifiable_k1p.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k1'", color_facet = "red", legend = "none", font_size=24, xylim = c(-1.75,0.75)) # xylim=c(-1.5,1.1)
p 
# dev.off()

# svg(filename = 'Figure3C_reproducibility_identifiable_k2.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k2", color_facet = "yellow", legend = "none", font_size=24, xylim=c(-3,0.5)) # xylim=c(-1.5,1.1)
p 
#dev.off()

#svg(filename = "Figure3C_reproducibility_identifiable_k2'.svg", width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k2'", color_facet = "green", xylim=c(-3,0), legend = "none", font_size=24)
p
#dev.off()

#svg(filename = 'Fig3C_reproducibility_identifiable_kdeg.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "kdeg", color_facet = "blue", xylim=c(-4,0), legend = "none", font_size=24)
p
#dev.off()

#svg(filename = 'Fig3C_reproducibikity_legend.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "kdeg", color_facet = "blue", xylim=c(-4,0), legend = "only", font_size=24)
grid.draw(p)
#dev.off()

#svg(filename = 'Fig3D_k1prime_over_k2_identifiable.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k1'/k2", xylim=c(-0.2,1.2), legend = "none", font_size=24)
p
#dev.off()

# svg(filename = 'Fig3D_k1primek2prime_over_k2_identifiable.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k1'k2'/k2", xylim=c(-2.7,0.5), legend = "none", font_size=24)
p
# dev.off()

#svg(filename = 'Fig3D_k2prime_over_k2_identifiable.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='naive', parameter =  "k2'/k2", xylim=c(-2,1), legend = "none", font_size=24)
p
#dev.off()
stop()

dat_ggplot$color <- annot[paste(dat_ggplot$Gene_name, dat_ggplot$Condition, dat_ggplot$Batch, sep="_"),1]

#### Figure 3E - Comparaison of decay with ActD HL ####---------------------------
dat_ggplot3 <- cbind(dat_ggplot[dat_ggplot$Condition == 'naive', ],
                     dat_ggplot[dat_ggplot$Condition == 'lpa',][match(apply(dat_ggplot[dat_ggplot$Condition == 'naive', c(1,3,11)],1,function(x){paste0(x, collapse=".")}), apply(dat_ggplot[dat_ggplot$Condition == 'lpa', c(1,3,11)],1,function(x){paste0(x, collapse=".")})),][,c(2,4:9,14)])

colnames(dat_ggplot3)[c(2,4:9,14)] <- paste0(colnames(dat_ggplot3)[c(2,4:9,14)],'_Naive')
colnames(dat_ggplot3)[15:22] <- paste0(colnames(dat_ggplot3)[15:22],'_LPA')
dat_ggplot3 <- dat_ggplot3[,-10]

dat_ggplot3$color_max <- apply(dat_ggplot3[,paste("color", c("Naive", "LPA"), sep="_")],1, max)

dat_ggplot_HL <- dat_ggplot3[dat_ggplot3$Parameter_name == 'kdeg',]

dat_ggplot_HL$ActD_HL <- BMDM_dataset$halflife$Naive_Basal$HL[match(dat_ggplot_HL$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL$ActD_CI95_L <- BMDM_dataset$halflife$Naive_Basal$CI_95_LB[match(dat_ggplot_HL$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL$ActD_CI95_U <- BMDM_dataset$halflife$Naive_Basal$CI_95_UB[match(dat_ggplot_HL$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL$HL_Naive <- log(2)/10^dat_ggplot_HL$Parameter_value_Naive
dat_ggplot_HL$HL_CIU_0.95_Min_Naive <- log(2)/10^dat_ggplot_HL$CIL_0.95_Max_Naive
dat_ggplot_HL$HL_CIU_0.95_Max_Naive <- log(2)/10^dat_ggplot_HL$CIL_0.95_Min_Naive
dat_ggplot_HL$HL_CIL_0.95_Min_Naive <- log(2)/10^dat_ggplot_HL$CIU_0.95_Max_Naive
dat_ggplot_HL$HL_CIL_0.95_Max_Naive <- log(2)/10^dat_ggplot_HL$CIU_0.95_Min_Naive


dat_ggplot_HL_2 <- cbind(dat_ggplot_HL[dat_ggplot_HL$Batch == 'rep1',], dat_ggplot_HL[dat_ggplot_HL$Batch == 'rep2',23:30])
colnames(dat_ggplot_HL_2)[31:38] <- paste(colnames(dat_ggplot_HL)[23:30], '_rep2', sep='')

dat_ggplot_HL$Batch <- factor(dat_ggplot_HL$Batch, levels=c('all','rep1','rep2','rep3'), labels=c('Together', 'Rep1', 'Rep2', 'Rep3'))

# svg(filename = 'Fig3E_ident.svg', width=5, height=5)
p <- ggplot(dat_ggplot_HL[dat_ggplot_HL$Batch %in% c('Rep1','Rep2') & is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive),], aes(x=HL_Naive,y=ActD_HL)) + geom_point(mapping = aes(shape=Batch))
p <- p + geom_segment(data = dat_ggplot_HL_2[is.na(dat_ggplot_HL_2$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL_2$HL_CIL_0.95_Min_Naive) & is.na(dat_ggplot_HL_2$HL_CIU_0.95_Max_Naive_rep2) & is.na(dat_ggplot_HL_2$HL_CIL_0.95_Min_Naive_rep2),], aes(x = HL_Naive, xend = HL_Naive_rep2, y = ActD_HL, yend = ActD_HL), alpha=0.1, linetype=2)
p <- p + scale_x_continuous(trans='log10', limits=c(5,10000)) + scale_y_continuous(trans='log10', limits = c(10,500))
p <- p + geom_abline(slope = 1, col='green') 
p <- p + geom_abline(slope = 1, intercept = log10(2), col='red', linetype=2)
p <- p + geom_abline(slope = 1, intercept = -log10(2), col='red', linetype=2)
p <- p + theme_bw(base_size = 24) + theme(aspect.ratio=1)
p <- p + labs(x='Model HL (min)', y='ActD HL (min)', shape="Replicate")
p1 <- p + theme(legend.position = "none")
p1
# dev.off()

cor(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)], use='pairwise.complete.obs', method="spearman")
cor(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)& dat_ggplot_HL$Batch=="Rep1"],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)& dat_ggplot_HL$Batch=="Rep1"], use='pairwise.complete.obs', method="spearman")
cor(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Rep2"],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Rep2"], use='pairwise.complete.obs', method="spearman")
cor(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Together"],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)&dat_ggplot_HL$Batch=="Together"], use='pairwise.complete.obs', method="spearman")

cor.test(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)], use='pairwise.complete.obs', method="spearman")$p.value
cor.test(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Rep1"],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Rep1"], use='pairwise.complete.obs', method="spearman")$p.value
cor.test(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Rep2"],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Rep2"], use='pairwise.complete.obs', method="spearman")$p.value
cor.test(dat_ggplot_HL$HL_Naive[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive)& dat_ggplot_HL$Batch=="Together"],dat_ggplot_HL$ActD_HL[is.na(dat_ggplot_HL$HL_CIU_0.95_Max_Naive) & is.na(dat_ggplot_HL$HL_CIL_0.95_Min_Naive) & dat_ggplot_HL$Batch=="Together"], use='pairwise.complete.obs', method="spearman")$p.value


# svg(filename = 'Figures/Fig3D_lgd.svg', width=5, height=5)
lgd <- gglegend(p)
grid.draw(lgd)
# dev.off()


quantile(dat_ggplot_HL$ActD_HL, c(0.05,0.25,0.5,0.75,0.95), na.rm = T)
quantile(dat_ggplot_HL$HL_Naive[dat_ggplot_HL$Batch == 'Rep1'], c(0.05,0.25,0.5,0.75,0.95), na.rm = T)
quantile(dat_ggplot_HL$HL_Naive[dat_ggplot_HL$Batch == 'Rep2'], c(0.05,0.25,0.5,0.75,0.95), na.rm = T)

cor(log(dat_ggplot_HL$ActD_HL[dat_ggplot_HL$Batch == 'Rep1']),log(dat_ggplot_HL$HL_Naive[dat_ggplot_HL$Batch == 'Rep1']),use = 'pairwise.complete.obs', method = 'spearman')
cor(log(dat_ggplot_HL$ActD_HL[dat_ggplot_HL$Batch == 'Rep2']),log(dat_ggplot_HL$HL_Naive[dat_ggplot_HL$Batch == 'Rep2']),use = 'pairwise.complete.obs', method='spearman')


#### Figure 4A - Example gene location on transport distribution ####---------------------------------------------------------------
# svg(filename = 'Figures/Fig4A_transport_dist_v13_2.svg', width=5, height=2.5)
p <- ggplot(dat_ggplot2[dat_ggplot2$Parameter_name %in% c("k1'k2'/k2") & dat_ggplot2$Condition == "naive" & apply(dat_ggplot2[,c("CIL_0.95_Min_rep1", "CIU_0.95_Max_rep1", "CIL_0.95_Min_rep2", "CIU_0.95_Max_rep2")],1,function(x){all(is.na(x))}),]) 
p <- p + geom_density(mapping = aes(x=Parameter_value_rep1, fill="k1'k2'/k2", color="k1'k2'/k2", linetype="Rep1"), alpha=0.1)
p <- p + geom_density(mapping = aes(x=Parameter_value_rep2, fill="k1'k2'/k2", color="k1'k2'/k2", linetype="Rep2"), alpha=0.1)
p <- p + scale_x_continuous(limits = c(-3,0.5)) 
p <- p + scale_color_manual(values=c("k1'k2'/k2"="black", "Egr1"="blue", "Malt1"="darkred"), guide=FALSE) 
p <- p + scale_fill_manual(values=c("k1'k2'/k2"="black"), guide=FALSE)
p <- p + scale_linetype_manual(values=c("Rep1"=1, "Rep2" = 2))
p <- p + geom_segment(data = dat_ggplot2[dat_ggplot2$Gene_name == "Egr1" & dat_ggplot2$Parameter_name %in% c("k1'k2'/k2") & dat_ggplot2$Condition == "naive",], aes(x = Parameter_value_rep1, xend=Parameter_value_rep1, y=0, yend=0.1, linetype="Rep1", col = "Egr1")) 
p <- p + geom_segment(data = dat_ggplot2[dat_ggplot2$Gene_name == "Egr1" & dat_ggplot2$Parameter_name %in% c("k1'k2'/k2") & dat_ggplot2$Condition == "naive",], aes(x = Parameter_value_rep2, xend=Parameter_value_rep2, y=0, yend=0.1, linetype="Rep2", col="Egr1")) 
p <- p + geom_segment(data = dat_ggplot2[dat_ggplot2$Gene_name == "Malt1" & dat_ggplot2$Parameter_name %in% c("k1'k2'/k2")& dat_ggplot2$Condition == "naive",], aes(x = Parameter_value_rep1, xend=Parameter_value_rep1, y=0, yend=0.1, linetype="Rep1", col="Malt1"))
p <- p + geom_segment(data = dat_ggplot2[dat_ggplot2$Gene_name == "Malt1" & dat_ggplot2$Parameter_name %in% c("k1'k2'/k2")& dat_ggplot2$Condition == "naive",], aes(x = Parameter_value_rep2, xend=Parameter_value_rep2, y=0, yend=0.1, linetype="Rep2", col="Malt1")) 
p <- p + annotate("text", x = mean(unlist(dat_ggplot2[dat_ggplot2$Gene_name == "Egr1" & dat_ggplot2$Parameter_name %in% c("k1'k2'/k2") & dat_ggplot2$Condition == "naive",c("Parameter_value_rep1","Parameter_value_rep2")]), na.rm=T), y = 0.15, label = "Egr1", color="blue", fontface="bold", size=rel(6))
p <- p + annotate("text", x = mean(unlist(dat_ggplot2[dat_ggplot2$Gene_name == "Malt1" & dat_ggplot2$Parameter_name %in% c("k1'k2'/k2") & dat_ggplot2$Condition == "naive",c("Parameter_value_rep1","Parameter_value_rep2")]), na.rm=T), y = 0.15, label = "Malt1", color="darkred", fontface="bold", size=rel(6))
p + theme_bw() + theme(aspect.ratio=0.5) + labs(x="Effective Transport (log10)", fill="Parameter", color="Gene", y="Density", linetype="Replicate")
# dev.off()

stop()



#### Figure 4B - distribution of transport efficiency ####
# svg(filename = 'Figures/Fig4B_efficiency_dist_v13_2.svg', width=5, height=2.5)
p <- ggplot(dat_ggplot2[dat_ggplot2$Parameter_name %in% c("k2'/k2") & dat_ggplot2$Condition == "naive" & apply(dat_ggplot2[,c("CIL_0.95_Min_rep1", "CIU_0.95_Max_rep1", "CIL_0.95_Min_rep2", "CIU_0.95_Max_rep2")],1,function(x){all(is.na(x))}),]) 
p <- p + geom_density(mapping = aes(x=Parameter_value_rep1, fill="k2'/k2", color="k2'/k2", linetype="Rep1"), alpha=0.1)
p <- p + geom_density(mapping = aes(x=Parameter_value_rep2, fill="k2'/k2", color="k2'/k2", linetype="Rep2"), alpha=0.1)
p <- p + scale_x_continuous(limits = c(-2,1)) 
p <- p + scale_color_manual(values=c("k2'/k2"="black", "Egr1"="blue", "Malt1"="darkred"), guide=FALSE) 
p <- p + scale_fill_manual(values=c("k2'/k2"="black"), guide=FALSE)
p <- p + scale_linetype_manual(values=c("Rep1"=1, "Rep2" = 2))
p + theme_bw() + theme(aspect.ratio=0.5) + labs(x="Transport Efficiency (log10)", fill="Parameter", color="Gene", y="Density", linetype="Replicate")
# dev.off()

stop()

#### Figure 4C - Comparison transport vs efficiency####
dat_ggplot4 <- rbind(cbind(dat_ggplot2[,c(1:5,7:8,10,12:13,20)], Batch="Rep1"),setNames(cbind(dat_ggplot2[,c(1,15,3,16:19,10,12:13, 21)], Batch="Rep2"),c(names(dat_ggplot2[,c(1:5,7:8,10,12:13,20)]),"Batch")))
colnames(dat_ggplot4) <- gsub("_rep1","", colnames(dat_ggplot4))
dat_ggplot4 <-  cbind(dat_ggplot4[dat_ggplot4$Parameter_name=="k1'k2'/k2",], dat_ggplot4[dat_ggplot4$Parameter_name=="k2'/k2",c(2,4:7)])
colnames(dat_ggplot4)[13:17] <- paste0(colnames(dat_ggplot4)[c(2,4:7)], "_eff")

# svg(filename = 'Figures/Fig4B_efficiency_vs transport.svg', width=5, height=5)
p <- ggplot(dat_ggplot4[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))}),]) 
p <- p + geom_point(mapping = aes(x=Parameter_value, y= Parameter_value_eff, color=color, shape=Batch))
p <- p + geom_segment(mapping = aes(x=CIL_0.95_Min,y=Parameter_value_eff,xend=CIL_0.95_Max,yend = Parameter_value_eff, col=color), alpha=0.3, linetype=2)
p <- p + geom_segment(mapping = aes(x=CIU_0.95_Min,y=Parameter_value_eff,xend=CIU_0.95_Max,yend = Parameter_value_eff, col=color), alpha=0.3, linetype=2)
p <- p + geom_segment(mapping = aes(x=CIL_0.95_Max,y=Parameter_value_eff,xend=CIU_0.95_Min,yend = Parameter_value_eff, col=color), alpha=0.3)
p <- p + geom_segment(mapping = aes(y=CIL_0.95_Min_eff,x=Parameter_value,yend=CIL_0.95_Max_eff,xend = Parameter_value, col=color), alpha=0.3, linetype = 2)
p <- p + geom_segment(mapping = aes(y=CIU_0.95_Min_eff,x=Parameter_value,yend=CIU_0.95_Max_eff,xend = Parameter_value, col=color), alpha=0.3, linetype = 2)
p <- p + geom_segment(mapping = aes(y=CIL_0.95_Max_eff,x=Parameter_value,yend=CIU_0.95_Min_eff,xend = Parameter_value, col=color), alpha=0.3)
p <- p + theme_bw(base_size = 24) + theme(aspect.ratio=1) + scale_color_gradientn(colours = c("green3", "gold2", "darkorange2", "darkred"), limits=range(c(dat_ggplot2$color_rep1,dat_ggplot2$color_rep2)))+labs(col='Fit Quality', shape="Replicate")
p <- p + coord_cartesian(xlim=c(-2.5,0.5), ylim=c(-2.5,1.5))
p <- p + labs(y="Transport Efficiency (log10)", x = "Effective transport (log10)")
p + theme(legend.position = "none")
# dev.off()

cor_test <- cor.test(dat_ggplot4$Parameter_value[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))})],dat_ggplot4$Parameter_value_eff[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))})], use = 'pairwise.complete.obs')
cor_test
cor_test <- cor.test(dat_ggplot4$Parameter_value[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))}) & dat_ggplot4$Batch == "Rep1"],dat_ggplot4$Parameter_value_eff[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))}) & dat_ggplot4$Batch == "Rep1"], use = 'pairwise.complete.obs')
cor_test
cor_test <- cor.test(dat_ggplot4$Parameter_value[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))}) & dat_ggplot4$Batch == "Rep2"],dat_ggplot4$Parameter_value_eff[ dat_ggplot4$Condition == "naive" & apply(dat_ggplot4[,c(4,7,14,17)],1,function(x){all(is.na(x))}) & dat_ggplot4$Batch == "Rep2"], use = 'pairwise.complete.obs')
cor_test

# svg(filename = 'Figures/Fig4B_efficiency_vs transport_leg.svg', width=5, height=5)
leg <- gglegend(p)
grid.draw(leg)
# dev.off()

##### Genes controlling Tolerance ######-------------------------
gene_tol <- getBM(attributes=c('external_gene_name', 'go_id','name_1006', 'namespace_1003'), filter='external_gene_name', values=row.names(all_best_param$naive$all),mart = mm_89)
gene_gobp <- gene_tol[gene_tol$namespace_1003 == "biological_process",]

gene_cat <- list("Cell growth genes" = c("Cd44", "Cdkn1a", "Cxcl16", "Dcbld2", "Edn1", "Gng4", "Hbegf", "Ifrd1", "Mmp14", "Sphk1"),
                 "Death or survival genes" = c("Bcl2l11", "Ppp1r15a", "Sod2", "Clec4e"),
                 "Cell adhesion genes" =  c("Adora2a", "Bmp2", "Cav1", "Ccl2", "Ccl5", "Cd44", "Cd74", "Cd83", "Dusp1", "Icam1", "Icosl", "Ifnrep1", "Il10", "Il12b", "Il1b", "Il1rn", "Il23a", "Il6", "Irf1", "Itga5", "Itgav", "Jag1", "Malt1", "Mmp14", "Nfkbid", "Pdpn", "Plau", "Rnd1", "S100a10", "Sdc4", "Serpine1", "Spp1", "Tgm2", "Tnf", "Tnfaip3", "Tnfsf4", "Tnfsf9", "Nlrp3", "Src"),
                 "Tissue remodeling genes" = c("Mmp13", "Mmp14", "Spp1", "Marcksl1"),
                 "Inflammatory response TF" = c("Junb", "Egr1", "Fos", "Egr2", "Fosb", "Fosl2", "Irf1", "Rel", "Relb"),
                 "Tolerance genes" = c("Nfkbia", "Tnfaip3", "Nfkrep3", "Nfkbid", "Nfkbib", "Tnfaip2", "Nfkbie", "Dusp1", "Dusp2", "Zfp36", "Socs3", "Dusp4", "Dusp8", "Dusp14", "Dusp16", "Dusp5", "Phlda1", "Il1rn"),
                 "Positive regulation of MAPK" = c("Bmp2", "Cd74",  "Gdf15", "Il1b", "Map2k3", "Tnf"),     
                 "Negative regulation of MAPK" = c("Cav1",  "Dusp1", "Dusp14", "Dusp16", "Dusp2", "Dusp4", "Dusp5", "Dusp8", "Spred1"),  
                 "Negative regulation of NFKB" = c("Nfkbia", "Nfkbib", "Nfkbid", "Tnfaip3", "Nfkrep3", "Tnfaip2", "Nfkbie"),
                 "Cytokines" = c("Tnf", "Ifnrep1", "Il1a", "Il10", "Il6", "Il1b", "Il23a", "Il12b"),
                 "Chemokine" = c("Cxcl10", "Ccl7", "Cxcl2", "Cxcl1", "Ccl2", "Ccl4", "Ccl12", "Ccl9", "Cxcl9", "Ccl5", "Cxcl16"))

dat_ggplot_6 <- cbind(dat_ggplot[dat_ggplot$Parameter_name == c("k1'k2'/k2") & dat_ggplot$Condition == "naive" & dat_ggplot$Batch %in% c("rep1","rep2") & apply(dat_ggplot[,c("CIL_0.95_Min", "CIU_0.95_Max")],1,function(x){all(is.na(x))}),], Cat = "All induced genes")
for ( go  in names(gene_cat)){
  for (gene in gene_cat[[go]]){
    tmp <- dat_ggplot[dat_ggplot$Gene_name == gene & dat_ggplot$Parameter_name == c("k1'k2'/k2") & dat_ggplot$Condition == "naive" & dat_ggplot$Batch %in% c("rep1","rep2") & apply(dat_ggplot[,c("CIL_0.95_Min", "CIU_0.95_Max")],1,function(x){all(is.na(x))}),]
    if(nrow(tmp) > 0){
      dat_ggplot_6 <- rbind(dat_ggplot_6, cbind(tmp, Cat = go))
    }
  }
}

dat_ggplot_6$Cat <- factor(dat_ggplot_6$Cat, levels = rev(c("All induced genes", names(gene_cat))))
dat_ggplot_6$Batch <- factor(dat_ggplot_6$Batch, levels = c("rep1","rep2"), labels=c("Rep1", "Rep2"))

library(ggridges)
# svg(filename = 'Fig4E_gene_cat.svg', width=7)
p <- ggplot(dat_ggplot_6, aes(y=Cat, x=Parameter_value, fill=Cat, linetype=Batch, shape=Batch)) 
p <- p + geom_density_ridges(mapping = aes(point_shape = Batch), panel_scaling=F, scale=0.95,
                             jittered_points = TRUE,
                             position = position_points_jitter(width = 0, height = 0), point_size = 2, point_alpha = 1, alpha = 0.5
)
p <- p + theme_bw()
p <- p + guides(fill=FALSE) + labs(y="", x="Effective transport rate (log10)", point_shape="Replicates", linetype="Replicates") 
print(p) 
# dev.off()

#### Figure S6B - Distri param lpa #####
dat_lpa <- data.frame(batch = "rep1", 
                      do.call("rbind", lapply(1:length(all_param_CI$lpa$rep1), function(i){cbind(Gene=names(all_param_CI$lpa$rep1)[i],all_param_CI$lpa$rep1[[i]]$CI[,c(1:4,7:8)])})))
dat_lpa <- rbind(dat_lpa, data.frame(batch = "rep2", 
                                     do.call("rbind", lapply(1:length(all_param_CI$lpa$rep2), function(i){cbind(Gene=names(all_param_CI$lpa$rep2)[i],all_param_CI$lpa$rep2[[i]]$CI[,c(1:4,7:8)])}))))

dat_lpa$batch <- factor(dat_lpa$batch, levels=c("rep1", "rep2"), labels=c("Rep1", "Rep2"))

head(dat_lpa)
dat_lpa$lower_min[is.na(dat_lpa$lower_min)] <- dat_lpa$lower[is.na(dat_lpa$lower_min)]
dat_lpa$upper_max[is.na(dat_lpa$upper_max)] <- dat_lpa$upper[is.na(dat_lpa$upper_max)]
head(dat_lpa)

dat_lpa$facet_label <- factor(dat_lpa$whichPar, levels=unique(dat_lpa$whichPar), labels=c("log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](k[1]*minute/k[2])", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])", "log[10](k[2]*minute/k[2])"))
dat_lpa$color_label <- factor(dat_lpa$whichPar, levels=unique(dat_lpa$whichPar), labels=c("k[1]*minute", "k[2]", "k[2]*minute", "k[deg]", "k[1]*minute/k[2]", "paste(k[1]*minute,k[2]*minute)/k[2]","k[2]*minute/k[2]"))

dat_lpa$CIrange_max <- dat_lpa$upper_max-dat_lpa$lower_min
dat_lpa$CIrange_min <- dat_lpa$upper-dat_lpa$lower

dat_lpa$identifiability <- NA
dat_lpa$identifiability[is.finite(dat_lpa$lower_min) & is.finite(dat_lpa$upper_max)] <- "Identifiable"
dat_lpa$identifiability[!is.finite(dat_lpa$lower_min) & !is.finite(dat_lpa$upper_max)] <- "Non-identifiable"
dat_lpa$identifiability[is.finite(dat_lpa$lower_min) & !is.finite(dat_lpa$upper_max)] <- "LB identifiable"
dat_lpa$identifiability[!is.finite(dat_lpa$lower_min) & is.finite(dat_lpa$upper_max)] <- "UB identifiable" 

dat_lpa$identifiability <- factor(dat_lpa$identifiability, levels=c("Non-identifiable", 'UB identifiable', "LB identifiable", "Identifiable"))
colors <- c("k1'"="red", "k2"="yellow", "k2'"="green", "kdeg"="blue", "k1'/k2"="black","k1'k2'/k2"="black", "k2'/k2"="black")
for (par in  c("k1'", "k2", "k2'", "kdeg", "k1'/k2", "k1'k2'/k2", "k2'/k2")){
  # svg(filename =  paste0('Figures/FigureS6B_identifiability_lpa_',gsub('/','over',par),'.svg'), width=3)
  p <- ggplot(dat_lpa[dat_lpa$whichPar == par,]) 
  p <- p + geom_bar(aes(x=batch, fill=identifiability)) + facet_grid(~whichPar)
  p <- p + theme_bw() + labs(fill="Identifiability", y="Number", x="")
  p <- p + theme(aspect.ratio=2, strip.background = element_rect(color=colors[par], fill=alpha(colors[par],0.1), size=1.5, linetype="solid"))
  p <- p + scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)
  print(p)
  # dev.off()
}

#### Figure S6A right plot of CIrange distribution of identifiable genes for k1', k2' ####
# svg(filename = 'Figures/FigureS6B_distribution of CIrange_lpa.svg', width=6)
p <- ggplot(dat_lpa[is.finite(dat_lpa$CIrange_max) & dat_lpa$whichPar %in% c("k1'", "k2", "k2'", "kdeg"),]) 
bw <- bw.ucv(dat_lpa$CIrange_min)
p <- p + geom_density(mapping = aes(x=CIrange_min, fill=whichPar, color=whichPar, linetype=batch), alpha=0.1, bw="ucv")
p <- p + scale_color_manual(name="Parameter", values=c("k1'"="red", "k2"="yellow", "k2'" = "green", "kdeg"="blue"))
p <- p + scale_fill_manual(name="Parameter", values=c("k1'"="red", "k2"="yellow", "k2'" = "green", "kdeg"="blue"))
p <- p + scale_linetype_manual(values=c("Rep1"=1, "Rep2" = 2), labels = c("Rep1", "Rep2"))
p <- p + theme_bw() + theme(aspect.ratio=1/4) + labs(x=expression(log[10]("95% CIU/CIL")), linetype="Replicate")
p
# dev.off()

#### Figure S6C - Distribution of parameter fits in Tolerized ####----------------
# svg(filename = 'Figures/FigS6C_distri_param_lpa.svg', width=7, height=4)
p <- ggplot(dat_ggplot[dat_ggplot$Parameter_name %in% c("k1'", "k2", "k2'", "kdeg") & dat_ggplot$Condition == "lpa" & dat_ggplot$Batch %in% c("rep1", "rep2") & is.na(dat_ggplot$CIL_0.95_Min) & is.na(dat_ggplot$CIU_0.95_Max),]) 
p <- p + geom_density(mapping = aes(x=Parameter_value, fill=Parameter_name, color=Parameter_name, linetype=Batch), alpha=0.1)
p <- p + scale_color_manual(values=c("red", "yellow","green", "blue")) 
p <- p + scale_fill_manual(values=c("red", "yellow", "green", "blue"))
p <- p + scale_linetype_manual(values=c("rep1"=1, "rep2" = 2), labels = c("Rep1", "Rep3"))
p <- p + theme_bw() + theme(aspect.ratio=1/4) + labs(x="Parameter value (log10)", fill="Parameter", color="Parameter", linetype="Replicate")
# p <- p + scale_x_continuous(limits=c(-3.5,8.5), breaks=c(seq(-5,8, by=2.5)))
p
# dev.off()
stop()

#### Figure S6C - Reproducibility of parameter fits in Tolerized ####----------------
# svg(filename = 'Figures/FigS6C_k1prime_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "k1'", color_facet = "red", legend = "none", font_size=24, xylim = c(-2,0.5)) # xylim=c(-1.5,1.1)
p 
# dev.off()

# svg(filename = 'Figures/FigS6C_k2_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "k2", color_facet = "yellow", legend = "none", font_size=24,xylim = c(-3,0)) # xylim=c(-1.5,1.1)
p
# dev.off()

# svg(filename = 'Figures/FigS6C_k2prime_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "k2'", color_facet = "green", xylim=c(-2.7,-0.2), legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/FigS6C_kdeg_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "kdeg", color_facet = "blue", xylim=c(-3.5,-0.3), legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/FigS6C_legend_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "kdeg", color_facet = "blue", xylim=c(-4,0), legend = "only", font_size=24)
grid.draw(p)
# dev.off()

# svg(filename = 'Figures/FigS6C_k1prime_over_k2_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "k1'/k2", legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/FigS6C_transport_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "k1'k2'/k2", xylim=c(-2.7,0.5), legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/FigS6C_transport_efficiency_lpa_rep.svg', width=5, height=5)
p <- make_replicate_plot2(data = dat_ggplot2, condition='lpa', parameter =  "k2'/k2", xylim=c(-2,1.5), legend = "none", font_size=24)
p
# dev.off()

stop()

#### Create Supplementary XLS with all parameters and CI ####---------------------
writeXLSwithType <- function(wb, sheet_name, data){
  # write header 
  openxlsx::writeData(wb, sheet_name, t(colnames(data)), startCol = 2, startRow = 1, rowNames = FALSE, colNames = FALSE)
  openxlsx::addStyle(wb, sheet_name, style = openxlsx::createStyle( textDecoration = "bold", border = "Bottom"), rows = 1, cols=1:(ncol(data)+1))
  # write row.names
  openxlsx::writeData(wb, sheet_name, row.names(data), startCol = 1, startRow = 2, rowNames = FALSE, colNames = FALSE)
  # write cell by cell to keep type
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      openxlsx::writeData(wb, sheet_name, data[i,j], startRow = i+1, startCol = j+1, rowNames = FALSE, colNames = FALSE)
    }
  }
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = 1:(ncol(tmp2)+1), widths = "auto")
}

wb <- openxlsx::createWorkbook()
openxlsx::modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Arial")

for (condition in c('naive', 'lpa')){
  if(condition == "naive"){
    condition_name <- "Naive"
  }else if(condition == "lpa"){
    condition_name <- "Tolerized"
  }
  
  for (batch in c('all','rep1','rep2','rep3')){
    sheet_name <- paste0(condition_name, "_", c("Together", "Rep1", "Rep2", "Rep3")[c('all','rep1','rep2','rep3') ==batch])
    openxlsx::addWorksheet(wb, sheet_name)
    
    tmp <- all_best_param[[condition]][[batch]][,colnames(all_best_param[[condition]][[batch]]) %in% c("k1'","k2","k2'","kdeg","k1'/k2", "k1'k2'/k2", "k2'/k2","k1'k2'", "max(k2,kdeg)", "min(k2,kdeg)")]
    tmp2 <- tmp
    for(col in colnames(tmp)){
      i <- which(colnames(tmp2) == col)
      if(i != ncol(tmp2)){
        tmp2 <- cbind(tmp2[,1:i], NA, tmp2[,(i+1):ncol(tmp2)]) 
      }else{
        tmp2 <- cbind(tmp2[,1:i], NA) 
      }
      if(i == 1){
        colnames(tmp2)[i] <- col
      }else if (i == ncol(tmp2)-2){
        colnames(tmp2)[ncol(tmp2)] <- colnames(tmp)[ncol(tmp)]
      }
      colnames(tmp2)[(i+1)] <- paste0(colnames(tmp2)[i],'_CI95')
    }
    
    for (gene in names(all_param_CI[[condition]][[batch]])){
      for(par in all_param_CI[[condition]][[batch]][[gene]]$CI$whichPar){
        tmp_CI <- all_param_CI[[condition]][[batch]][[gene]]$CI[all_param_CI[[condition]][[batch]][[gene]]$CI$whichPar == par,]
        tmp2[row.names(tmp2) == gene, paste0(par,'_CI95')] <- paste0("(", ifelse(!is.na(tmp_CI[,'lower_min']), paste0("< ", signif(tmp_CI[,"lower"],3)), signif(tmp_CI[,"lower"],3)),", ",
                                                                     ifelse(!is.na(tmp_CI[,'upper_max']), paste0("> ", signif(tmp_CI[,"upper"],3)), signif(tmp_CI[,"upper"],3)),")")
      }
    }
    
    for (par in c("k1'","k2","k2'","kdeg","k1'/k2", "k1'k2'/k2","k2'/k2", "k1'k2'", "max(k2,kdeg)", "min(k2,kdeg)")){
      if( par %in% colnames(tmp2)){
        column_i <- which(colnames(tmp2) == par)
        colnames(tmp2)[column_i + 0:1] <- gsub(par,paste0("log10(",par,")"), colnames(tmp2)[column_i + 0:1], fixed=T)
      }
    }
    
    if(batch == "rep3"){
      colnames(tmp2) <- gsub("max(k2,kdeg)", "k2 or kdeg", colnames(tmp2), fixed=T)
      colnames(tmp2) <- gsub("min(k2,kdeg)", "kdeg or k2", colnames(tmp2), fixed=T)
    }
    
    writeXLSwithType(wb, sheet_name, tmp2)
  }
}
openxlsx::saveWorkbook(wb, "Parameters_fit.xlsx", overwrite = TRUE)

#### Figure 6B - Naive vs LPA parameter fits ####---------------------------------
make_condition_plot2 <- function(data, batch, parameter, xylim=NULL, color_facet=NULL, font_size=12, legend = "none"){
  
  p <- ggplot(data[data$Batch %in% batch & data$Parameter_name == parameter & apply(data[,c(4,8,15,19)],1,function(x){all(is.na(x))}),], aes(x=Parameter_value_Naive,y=Parameter_value_LPA))
  p <- p + geom_segment(mapping = aes(x=CIL_0.95_Min_Naive,y=Parameter_value_LPA,xend=CIL_0.95_Max_Naive,yend = Parameter_value_LPA, col=color_Naive), alpha=0.3, linetype=2)
  p <- p + geom_segment(mapping = aes(x=CIU_0.95_Min_Naive,y=Parameter_value_LPA,xend=CIU_0.95_Max_Naive,yend = Parameter_value_LPA, col=color_Naive), alpha=0.3, linetype=2)
  p <- p + geom_segment(mapping = aes(x=CIL_0.95_Max_Naive,y=Parameter_value_LPA,xend=CIU_0.95_Min_Naive,yend = Parameter_value_LPA, col=color_Naive), alpha=0.3)
  p <- p + geom_segment(mapping = aes(y=CIL_0.95_Min_LPA,x=Parameter_value_Naive,yend=CIL_0.95_Max_LPA,xend = Parameter_value_Naive, col=color_LPA), alpha=0.3, linetype = 2)
  p <- p + geom_segment(mapping = aes(y=CIU_0.95_Min_LPA,x=Parameter_value_Naive,yend=CIU_0.95_Max_LPA,xend = Parameter_value_Naive, col=color_LPA), alpha=0.3, linetype = 2)
  p <- p + geom_segment(mapping = aes(y=CIL_0.95_Max_LPA,x=Parameter_value_Naive,yend=CIU_0.95_Min_LPA,xend = Parameter_value_Naive, col=color_LPA), alpha=0.3)
  p <- p + geom_point(mapping=aes(col=color_max, shape=Batch), size=2) + facet_wrap(.~facet_label, scale='free', drop = T, labeller = "label_parsed")
  p <- p + theme_bw(base_size = font_size) + theme(aspect.ratio=1) + scale_color_gradientn(colours = c("green3", "gold2", "darkorange2", "darkred"), limits=range(c(data$color_Naive,data$color_LPA)))+labs(col='Fit Quality')
  p <- p + geom_abline(slope=1, intercept=0, linetype=1, col='black', alpha=0.5) + geom_abline(slope=1, intercept=log10(2), linetype=2, col='red', alpha=0.5) + geom_abline(slope=1, intercept=-log10(2), linetype=2, col='red', alpha=0.5)
  p <- p + scale_shape_manual(values=c("all"=18, "rep1"=16,"rep3"=15,"rep2"=17), labels = function(x){gsub('b','Rep',x)})
  p <- p + labs(x="Naive", y="LPA", shape="Replicate")
  p <- p + theme(strip.text = element_text(margin=margin(0.2,0,0.2,0,unit = "line")))
  if(!is.null(xylim)){
    p <- p + coord_cartesian(xlim=xylim, ylim=xylim)
  }
  if(!is.null(color_facet)){
    p <- p + theme(strip.background = element_rect(color=color_facet, fill=alpha(color_facet,0.1), size=1.5, linetype="solid"))
  }
  
  if(legend == "none"){
    p <- p + theme(legend.position = "none")
  }else if(legend == "only"){
    p <- gglegend(p)
  }
  
  return(p)
}

# svg(filename = 'Figures/Fig6B_k1prime_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "k1'", color_facet = "red", xylim= c(-2,0.5),legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/Fig6B_k2_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "k2", color_facet = "yellow", xylim= c(-3.5,0),legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/Fig6B_k2prime_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "k2'", color_facet = "green", xylim= c(-3,0),legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/Fig6B_kdeg_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "kdeg", color_facet = "blue", xylim = c(-4,-0.5), legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/Fig6B_legend_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "kdeg", color_facet = "blue", xylim= c(-4,0.5),legend = "only", font_size=24)
grid.draw(p)
# dev.off()

# svg(filename = 'Figures/Fig6B_k1prime_over_k2_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "k1'/k2", xylim= c(-0.2,2),legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/Fig6B_k1primek2prime_over_k2_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "k1'k2'/k2", xylim= c(-3,0.5),legend = "none", font_size=24)
p
# dev.off()

# svg(filename = 'Figures/Fig6B_k2prime_over_k2_v2.svg', width=5, height=5)
p <- make_condition_plot2(data = dat_ggplot3, batch = c('rep1','rep2'), parameter =  "k2'/k2", xylim= c(-2,2),legend = "none", font_size=24)
p
# dev.off()

stop()


#### Figure S4D - Comparison ActD HL Basal vs 1h ####-------------------------------
dat_ggplot_HL2 <- dat_ggplot3[dat_ggplot3$Parameter_name == 'kdeg',]

dat_ggplot_HL2$ActD_HL_basal <- BMDM_dataset$halflife$Naive_Basal$HL[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_CI95_L_basal <- BMDM_dataset$halflife$Naive_Basal$CI_95_LB[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_CI95_U_basal <- BMDM_dataset$halflife$Naive_Basal$CI_95_UB[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_HL_1h <- BMDM_dataset$halflife$Naive_LPA_1h$HL[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_CI95_L_1h <- BMDM_dataset$halflife$Naive_LPA_1h$CI_95_LB[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_CI95_U_1h <- BMDM_dataset$halflife$Naive_LPA_1h$CI_95_UB[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_HL_3h <- BMDM_dataset$halflife$Naive_LPA_3h$HL[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_CI95_L_3h <- BMDM_dataset$halflife$Naive_LPA_3h$CI_95_LB[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]
dat_ggplot_HL2$ActD_CI95_U_3h <- BMDM_dataset$halflife$Naive_LPA_3h$CI_95_UB[match(dat_ggplot_HL2$Gene_name, BMDM_dataset$gene_infos$Symbol[match(row.names(BMDM_dataset$halflife$Naive_Basal),BMDM_dataset$gene_infos$Geneid)])]

dat_ggplot_HL2$HL_Naive <- log(2)/10^dat_ggplot_HL2$Parameter_value_Naive
dat_ggplot_HL2$HL_CIU_0.95_Min_Naive <- log(2)/10^dat_ggplot_HL2$CIL_0.95_Max_Naive
dat_ggplot_HL2$HL_CIU_0.95_Max_Naive <- log(2)/10^dat_ggplot_HL2$CIL_0.95_Min_Naive
dat_ggplot_HL2$HL_CIL_0.95_Min_Naive <- log(2)/10^dat_ggplot_HL2$CIU_0.95_Max_Naive
dat_ggplot_HL2$HL_CIL_0.95_Max_Naive <- log(2)/10^dat_ggplot_HL2$CIU_0.95_Min_Naive

dat_ggplot_HL_2 <- cbind(dat_ggplot_HL2[dat_ggplot_HL2$Batch == 'rep1',], dat_ggplot_HL2[dat_ggplot_HL2$Batch == 'rep2',c(13,23:36)])
colnames(dat_ggplot_HL_2)[37:51] <- paste(colnames(dat_ggplot_HL_2)[c(13,23:36)], '_rep2', sep='')
dat_ggplot_HL_2$col_naive_max <- apply(dat_ggplot_HL_2[,c(13,37)],1,max)

p <- ggplot(data = dat_ggplot_HL_2, aes(x=ActD_HL_basal, y=ActD_HL_1h)) + geom_point()
p <- p + geom_segment(mapping = aes(x= ActD_CI95_L_basal, xend=ActD_CI95_U_basal, y=ActD_HL_1h, yend=ActD_HL_1h), alpha=0.1, linetype=2)
p <- p + geom_segment(mapping = aes(x= ActD_HL_basal, xend=ActD_HL_basal, y=ActD_CI95_L_1h, yend=ActD_CI95_U_1h), alpha=0.1, linetype=2)
p <- p + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + coord_cartesian(xlim=c(5,500), ylim = c(10,500))
p <- p + geom_abline(slope = 1, col='green') 
p <- p + geom_abline(slope = 1, intercept = log10(2), col='red', linetype=2)
p <- p + geom_abline(slope = 1, intercept = -log10(2), col='red', linetype=2)
p <- p + theme_bw() + theme(aspect.ratio=1) + labs(x="ActD HL at basal (min)", y="ActD HL at 1h after stimulation (min)")
p

# svg(filename = 'Figures/FigS3D.svg', width=7, height=7)
p <- ggplot(data = dat_ggplot_HL_2, aes(x=ActD_HL_basal, y=ActD_HL_1h)) + geom_point()
p <- p + geom_segment(mapping = aes(x= ActD_CI95_L_basal, xend=ActD_CI95_U_basal, y=ActD_HL_1h, yend=ActD_HL_1h), alpha=0.1, linetype=2)
p <- p + geom_segment(mapping = aes(x= ActD_HL_basal, xend=ActD_HL_basal, y=ActD_CI95_L_1h, yend=ActD_CI95_U_1h), alpha=0.1, linetype=2)
p <- p + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + coord_cartesian(xlim=c(5,500), ylim = c(10,500))
p <- p + geom_abline(slope = 1, col='green') 
p <- p + geom_abline(slope = 1, intercept = log10(2), col='red', linetype=2)
p <- p + geom_abline(slope = 1, intercept = -log10(2), col='red', linetype=2)
p <- p + theme_bw(base_size = 24) + theme(aspect.ratio=1) + labs(x="ActD HL at basal (min)", y="ActD HL at 1h after stimulation (min)")
p1 <- p + theme(legend.position = "none")
p1
# dev.off()
stop()

#### Figure 4D - Malt1 and Egr same graph ####-------------------------------------
dat_malt1_egr1 <- rbind(cbind(dat_g_Malt1, Name="Malt1"),cbind(dat_g_Egr1, Name="Egr1"))
dat_malt1_egr1$Cpt <- factor(dat_malt1_egr1$Cpt, levels=c("Chromatin", "Nucleoplasm", "Cytoplasm"),ordered=T)

# svg(filename = 'Figure4C_Malt1_Egr1_v2.svg', width=6, height=7)
p <- ggplot(data = dat_malt1_egr1[dat_malt1_egr1$Condition == "naive",], mapping = aes(x=Time,y=RPKM, col=Name, group=interaction(Name,batch))) + geom_point(data = dat_malt1_egr1[dat_malt1_egr1$Type == "Experiment" & dat_malt1_egr1$Condition == "naive",]) + facet_grid(Cpt~.)
p <- p + geom_line(data=dat_malt1_egr1[dat_malt1_egr1$Type == "Fit" & dat_malt1_egr1$Condition == "naive",], alpha=0.5)
p <- p + scale_color_manual(values = c(Malt1="darkred", Egr1="blue"))
p <- p + theme_bw() + theme(aspect.ratio=0.5) + scale_x_continuous(breaks = c(0,30,60,90,120))
p <- p + labs(x="Time (min)", y=parse(text="log[2](FPKM)"))
g <- ggplotGrob(p)
col <- c("Chromatin" = col_ca, "Nucleoplasm" = col_np, "Cytoplasm" = col_cyto)
for(i in grep("strip-r",g$layout$name)){
  Cpt <- g$grobs[[i]]$grobs[[1]]$children[[grep("text", names(g$grobs[[i]]$grobs[[1]]$children))]]$children[[1]]$label
  g$grobs[[i]]$grobs[[1]]$children[[grep("rect", names(g$grobs[[i]]$grobs[[1]]$children))]]$gp$fill <- col[Cpt]
}
grid.draw(g)
# dev.off()

stop()

#### Process Intron data #####---------------------------------------------------------
load('../Intron retention/all_intron_summary_measure.Rdata')

#### Figure 5A - Plot Transport vs gene_length ####-------------------------------------------
# color cytokines : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3140102/
cytokines <-  c(paste0('Il', c(1:34, '1a','1b',paste0('1f', c(5:10)),'12a','12b', paste0('17', c('', 'a', 'b','c', 'd', 'f')),'23a')),
                paste0('Ifn',c(paste0('a',c(1:16,'b')),'b1','g')),
                'Tnf')
chemokines<-c(paste0('Ccl', c(1:28,'27a', '27b','21a', '21b')),
              paste0('Cxcl', 1:17),
              'Cx3cl1','Xcl1')

# get gene length from manual TSS and TES
curated_list_naive$Length <- apply(curated_list_naive, 1, function(x){
  TSS <- as.numeric(strsplit(x[6], split = ", ")[[1]])
  TES <- as.numeric(strsplit(x[7], split = ", ")[[1]])
  length <- abs(mean(TES-TSS)) + 1
  return(length)
})

dat_ggplot5 <- dat_ggplot[dat_ggplot$Parameter_name == "k1'k2'/k2",]
colnames(dat_ggplot5)[2] <- "transport"
colnames(dat_ggplot5)[4:9] <- paste0(colnames(dat_ggplot5)[4:9], '_transport')

dat_ggplot5$Length <- curated_list_naive$Length[match(dat_ggplot5$Gene_name, curated_list_naive$Name)]/1000
dat_ggplot5$type <- "other"
dat_ggplot5$type[dat_ggplot5$Gene_name %in% cytokines] <- "cytokine"
dat_ggplot5$type[dat_ggplot5$Gene_name %in% chemokines] <- "chemokine"

dat_ggplot5_naive <- dat_ggplot5[dat_ggplot5$Condition == "naive",]

dat_ggplot5_naive <- dat_ggplot5_naive[dat_ggplot5_naive$Batch %in% c('rep1','rep2'),]

cor.test(dat_ggplot5_naive$transport[dat_ggplot5_naive$Batch == "rep1" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], log10(dat_ggplot5_naive$Length[dat_ggplot5_naive$Batch == "rep1" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})]), use="pairwise.complete.obs")
cor.test(dat_ggplot5_naive$transport[dat_ggplot5_naive$Batch == "rep2" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], log10(dat_ggplot5_naive$Length[dat_ggplot5_naive$Batch == "rep2" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})]), use="pairwise.complete.obs")


cor_test <- cor.test(dat_ggplot5_naive$transport[apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], log10(dat_ggplot5_naive$Length[ apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})]), use="pairwise.complete.obs")

library(lmerTest)
fit_length <- lmer(transport ~ Length + (1 | Gene_name), data = dat_ggplot5_naive)
summary(fit_length)


# svg(filename = 'Figures/Figure5A_transport_vs_Length_v3.svg', width=4,height=4)
p <- ggplot(dat_ggplot5_naive[apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))}),], aes(y=transport, x=Length, col=type)) + geom_point(mapping=aes(shape=Batch))
p <- p + theme_bw() + theme(aspect.ratio=1) + scale_x_continuous(trans='log10') + scale_y_continuous(limits=c(-2.6,0.5))
p <- p + scale_shape(labels = c("Rep1","Rep2")) + coord_cartesian(ylim=c(-2.5,0.5), xlim=c(0.75,300))
p <- p + geom_segment(aes(x=Length,xend=Length, y= CIL_0.95_Max_transport, yend=CIU_0.95_Min_transport), alpha=0.05)
p <- p + geom_segment(aes(x=Length,xend=Length, y= CIL_0.95_Min_transport, yend=CIL_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=Length,xend=Length, y= CIU_0.95_Min_transport, yend=CIU_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + labs(y=parse(text = "`Relative Transport `(log[10])"), x='Length (kbp)', shape="Replicate") + scale_color_manual(values=c('other'='grey3', 'cytokine'='red3','chemokine'='skyblue3'))
p <- p + annotate(geom = "text", x=100, y=0.5, hjust=1, vjust=1,label=bquote("Pearson's" ~ rho == .(signif(cor_test$estimate,2))))
p + theme(legend.position = "none") + coord_cartesian(xlim=c(1.5,100))
# dev.off()
stop()

#### Figure 5B - Plot Transport vs Nb of introns ####-------------------------------------------
# get number of intron from annotation for observed isoforms
gene_list_mm89_with_intron_number <- getBM(attributes = c("ensembl_gene_id","external_gene_name", "version", "ensembl_transcript_id", "start_position", "end_position", "rank"), filters = "ensembl_gene_id", values = substr(rownames(caRNA_Naive_exons_rpkm), 1,18), mart = mm_89) 
curated_list_naive$NbIntron <- apply(curated_list_naive, 1, function(x){
  isoforms <- substr(strsplit(x[8], split = ", ")[[1]],1,18)
  intron_num <- c()
  for (t in isoforms){
    dat_t <- gene_list_mm89_with_intron_number[gene_list_mm89_with_intron_number$ensembl_transcript_id == t,]
    intron_num <- c(intron_num, max(dat_t$rank)-1)
  }
  intron_nb <- intron_num[1]
  return(intron_nb)
})

dat_ggplot5$Nb_introns <- curated_list_naive$NbIntron[match(dat_ggplot5$Gene_name, curated_list_naive$Name)]

log10_1p <- function() trans_new("log10_1p", function(x) log10(1+x), function(x) 10^x-1)


dat_ggplot5_naive <- dat_ggplot5[dat_ggplot5$Condition == "naive",]

dat_ggplot5_naive <- dat_ggplot5_naive[dat_ggplot5_naive$Batch %in% c('rep1','rep2'),]

cor.test(dat_ggplot5_naive$transport[dat_ggplot5_naive$Batch == "rep1" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], dat_ggplot5_naive$Nb_introns[dat_ggplot5_naive$Batch == "rep1" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor.test(dat_ggplot5_naive$transport[dat_ggplot5_naive$Batch == "rep2" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], dat_ggplot5_naive$Nb_introns[dat_ggplot5_naive$Batch == "rep2" & apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], use="pairwise.complete.obs")


cor_test <- cor.test(dat_ggplot5_naive$transport[apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], dat_ggplot5_naive$Nb_introns[apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))})], use="pairwise.complete.obs")

fit_intron <- lmer(transport ~ Nb_introns + (1 | Gene_name), data = dat_ggplot5_naive)
tmp <- summary(fit_intron)


# svg(filename = 'Figures/Figure5B_transport_vs_Intons_v3.svg', width=8,height=4)
p <- ggplot(dat_ggplot5_naive[apply(dat_ggplot5_naive[,c(4,8)],1,function(x){all(is.na(x))}),], aes(y=transport, x=Nb_introns, col=type)) + geom_point(mapping=aes(shape=Batch))
p <- p + theme_bw() + theme(aspect.ratio=1) + scale_y_continuous(limits=c(-2.6,0.5))
p <- p + scale_shape(labels = c("Rep1","Rep2")) + coord_cartesian(ylim=c(-2.5,0.5), xlim=c(0,40))
p <- p + geom_segment(aes(x=Nb_introns,xend=Nb_introns, y= CIL_0.95_Max_transport, yend=CIU_0.95_Min_transport), alpha=0.05)
p <- p + geom_segment(aes(x=Nb_introns,xend=Nb_introns, y= CIL_0.95_Min_transport, yend=CIL_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=Nb_introns,xend=Nb_introns, y= CIU_0.95_Min_transport, yend=CIU_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + labs(y=parse(text = "`Relative Transport `(log[10])"), x='Average number of introns', shape="Replicate")
p <- p + annotate(geom = "text", x=15, y=0.5, hjust=1, vjust=1,label=bquote("Pearson's" ~ rho == .(signif(cor_test$estimate,2))))
p <- p + scale_color_manual(values=c('other'='grey3', 'cytokine'='red3','chemokine'='skyblue3'))
p + theme(legend.position = "none") + coord_cartesian(xlim=c(0,15))
# dev.off()
stop()


#### Figure 5C - Plot introns vs Transport ####--------------------------------------
# Add parameters
dat3 <- c()
for (b in names(all_intron_summary_measure)){
  for (g in all_intron_summary_measure[[b]]$Name){
    dat3 <- rbind(dat3, cbind(all_intron_summary_measure[[b]][all_intron_summary_measure[[b]]$Name == g,], Batch=b, 
                              t(unlist(all_best_param$naive[[b]][g, c("kdeg", "k1'k2'/k2", "k2'/k2")])),
                              data.frame(CIL_0.95_Min_kdeg = all_param_CI$naive[[b]][[g]]$CI["kdeg", "lower_min"],
                                         CIL_0.95_Max_kdeg = all_param_CI$naive[[b]][[g]]$CI["kdeg", "lower"],
                                         CIL_0.95_Method_kdeg = all_param_CI$naive[[b]][[g]]$CI["kdeg", "method_lower"],
                                         CIU_0.95_Min_kdeg = all_param_CI$naive[[b]][[g]]$CI["kdeg", "upper"],
                                         CIU_0.95_Max_kdeg = all_param_CI$naive[[b]][[g]]$CI["kdeg", "upper_max"],
                                         CIU_0.95_Method_kdeg = all_param_CI$naive[[b]][[g]]$CI["kdeg", "method_upper"]),
                              data.frame(CIL_0.95_Min_transport = all_param_CI$naive[[b]][[g]]$CI["k1'k2'/k2", "lower_min"],
                                         CIL_0.95_Max_transport = all_param_CI$naive[[b]][[g]]$CI["k1'k2'/k2", "lower"],
                                         CIL_0.95_Method_transport = all_param_CI$naive[[b]][[g]]$CI["k1'k2'/k2", "method_lower"],
                                         CIU_0.95_Min_transport = all_param_CI$naive[[b]][[g]]$CI["k1'k2'/k2", "upper"],
                                         CIU_0.95_Max_transport = all_param_CI$naive[[b]][[g]]$CI["k1'k2'/k2", "upper_max"],
                                         CIU_0.95_Method_transport = all_param_CI$naive[[b]][[g]]$CI["k1'k2'/k2", "method_upper"]),
                              data.frame(CIL_0.95_Min_eff = all_param_CI$naive[[b]][[g]]$CI["k2'/k2", "lower_min"],
                                         CIL_0.95_Max_eff = all_param_CI$naive[[b]][[g]]$CI["k2'/k2", "lower"],
                                         CIL_0.95_Method_eff = all_param_CI$naive[[b]][[g]]$CI["k2'/k2", "method_lower"],
                                         CIU_0.95_Min_eff = all_param_CI$naive[[b]][[g]]$CI["k2'/k2", "upper"],
                                         CIU_0.95_Max_eff = all_param_CI$naive[[b]][[g]]$CI["k2'/k2", "upper_max"],
                                         CIU_0.95_Method_eff = all_param_CI$naive[[b]][[g]]$CI["k2'/k2", "method_upper"])
    )
    )
  }
}
dat3$HL <- log(2)/10^dat3$kdeg

cor(log10(dat3$PI_naive_prod_avg_late)[apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1,function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor(log10(dat3$PI_naive_prod_avg_late)[ dat3$Batch == 'rep1' & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[ dat3$Batch == 'rep1' & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor(log10(dat3$PI_naive_prod_avg_late)[dat3$Batch == 'rep2' & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[ dat3$Batch == 'rep2' & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")

cor(dat3$kdeg[apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg","CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg","CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor(dat3$kdeg[dat3$Batch == 'rep1' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg","CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[dat3$Batch == 'rep1' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg", "CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor(dat3$kdeg[dat3$Batch == 'rep2' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg", "CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[dat3$Batch == 'rep2' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg", "CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")

cor(log10(dat3$PI_naive_prod_avg_late)[apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg")],1, function(x){all(is.na(x))})], dat3$kdeg[apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor(log10(dat3$PI_naive_prod_avg_late)[dat3$Batch == 'rep1' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg")],1, function(x){all(is.na(x))})], dat3$kdeg[dat3$Batch == 'rep1' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor(log10(dat3$PI_naive_prod_avg_late)[dat3$Batch == 'rep2' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg")],1, function(x){all(is.na(x))})], dat3$kdeg[dat3$Batch == 'rep2' & apply(dat3[,c("CIL_0.95_Min_kdeg","CIU_0.95_Max_kdeg")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")

cor_test <- cor.test(log10(dat3$PI_naive_prod_avg_late)[apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor.test(log10(dat3$PI_naive_prod_avg_late)[dat3$Batch == "rep1" & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[dat3$Batch == "rep1" & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")
cor.test(log10(dat3$PI_naive_prod_avg_late)[dat3$Batch == "rep2" & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3$`k1'k2'/k2`[dat3$Batch == "rep2" & apply(dat3[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")



dat3$transport <- dat3$`k1'k2'/k2`
dat3$type <- "other"
dat3$type[dat3$Name %in% cytokines] <- "cytokine"
dat3$type[dat3$Name %in% chemokines] <- "chemokine"

dat_ggplot5_naive$PI_naive_prod_avg_late <- dat3$PI_naive_prod_avg_late[match(interaction(dat_ggplot5_naive$Gene_name, dat_ggplot5_naive$Batch),interaction(dat3$Name, dat3$Batch))]

fit_splicing <- rlm(transport ~ PI_naive_prod_avg_late, data = dat3)
tmp <- summary(fit_splicing)

dat3_rep1rep2 <- dat3[dat3$Batch %in% c("rep1","rep2"),]

library(ggrepel)
cor_test <- cor.test(log10(dat3_rep1rep2$PI_naive_prod_avg_late)[apply(dat3_rep1rep2[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3_rep1rep2$`k1'k2'/k2`[apply(dat3_rep1rep2[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")

# svg(filename = 'Figures/Figure5C_transport vs PI_v3.svg', width=4, height=4)
p <- ggplot(dat3_rep1rep2[apply(dat3_rep1rep2[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))}),], aes(x=PI_naive_prod_avg_late, y=`k1'k2'/k2`, col=type)) + geom_point(mapping = aes(shape=Batch)) 
p <- p + theme_bw() + theme(aspect.ratio=1)
p <- p + scale_x_continuous(trans="log10") + scale_shape_discrete(labels=c("rep1"="Rep1", "rep2"="Rep2")) 
p <- p + labs(x = parse(text="prod((1-PI[intron]), intron)"), y=bquote("Relative Transport "~(log[10])), color="", shape="Replicate")
p <- p + geom_segment(aes(x=PI_naive_prod_avg_late,xend=PI_naive_prod_avg_late, y= CIL_0.95_Max_transport, yend=CIU_0.95_Min_transport), alpha=0.05)
p <- p + geom_segment(aes(x=PI_naive_prod_avg_late,xend=PI_naive_prod_avg_late, y= CIL_0.95_Min_transport, yend=CIL_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=PI_naive_prod_avg_late,xend=PI_naive_prod_avg_late, y= CIU_0.95_Min_transport, yend=CIU_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + coord_cartesian(ylim=c(-2.5,0.5))
p <- p + annotate(geom = "text", x=0.001, y=0.5, hjust=0, vjust=1,label=bquote("Pearson's" ~ rho == .(signif(cor_test$estimate,2))))
p <- p + scale_color_manual(values=c('other'='grey3', 'cytokine'='red3','chemokine'='skyblue3'))
p + theme(legend.position = "none")
# dev.off()

# svg(filename = 'Figures/Figure5C_legend_v3.svg', width=4, height=4)
p
# dev.off()

dat4 <- dat3_rep1rep2[apply(dat3_rep1rep2[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))}),]
d <- highlight_key(dat4, ~Name )
p <- ggplot(d, aes(x=PI_naive_prod_avg_late, y=`k1'k2'/k2`, col=type)) + geom_point(mapping = aes(shape=Batch)) 
p <- p + theme_bw()
p <- p + scale_x_continuous(trans="log10") + scale_shape_discrete(labels=c("rep1"="Rep1", "rep2"="Rep2")) 
p <- p + geom_segment(aes(x=PI_naive_prod_avg_late,xend=PI_naive_prod_avg_late, y= CIL_0.95_Max_transport, yend=CIU_0.95_Min_transport), alpha=0.05)
p <- p + geom_segment(aes(x=PI_naive_prod_avg_late,xend=PI_naive_prod_avg_late, y= CIL_0.95_Min_transport, yend=CIL_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=PI_naive_prod_avg_late,xend=PI_naive_prod_avg_late, y= CIU_0.95_Min_transport, yend=CIU_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + coord_cartesian(ylim=c(-2.5,0.5))
p <- p + theme(legend.position = "none")
pp <- ggplotly(p, tooltip="all")
highlight( pp, on = "plotly_hover", off = "plotly_deselect", color = "red" )
#

cor_test <- cor.test(log10(dat3_rep1rep2$PI_naive_prod_avg_late)[apply(dat3_rep1rep2[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], dat3_rep1rep2$`k1'k2'/k2`[apply(dat3_rep1rep2[,c("CIL_0.95_Min_transport","CIU_0.95_Max_transport")],1, function(x){all(is.na(x))})], use="pairwise.complete.obs")

#### Figure 5D - HL vs Transport ####---------------------
head(dat_ggplot)

dat_ggplot6 <- cbind(dat_ggplot[dat_ggplot$Parameter_name=="kdeg" & dat_ggplot$Batch %in% c('rep1','rep2'),-c(3,12,13)],
                     dat_ggplot[dat_ggplot$Parameter_name=="k1'k2'/k2" & dat_ggplot$Batch %in% c('rep1','rep2'),-c(1,3,10,11,12,13,14)])
names(dat_ggplot6)[2] <- "kdeg"
names(dat_ggplot6)[12] <- "transport"

dat_ggplot6 <- cbind(dat_ggplot6,data.frame(log(2)/10^dat_ggplot6[,c(2:4,6:7)],dat_ggplot6[,c(5,8)]))

names(dat_ggplot6)[3:8] <- paste0(names(dat_ggplot6)[3:8],'_kdeg')
names(dat_ggplot6)[13:18] <- paste0(names(dat_ggplot6)[13:18],"_transport")
names(dat_ggplot6)[19] <- "HL"
names(dat_ggplot6)[20:25] <- paste0(names(dat_ggplot6)[20:25],"_HL")
names(dat_ggplot6)[20:25] <- gsub("CIL","@", names(dat_ggplot6)[20:25])
names(dat_ggplot6)[20:25] <- gsub("CIU","CIL", names(dat_ggplot6)[20:25])
names(dat_ggplot6)[20:25] <- gsub("@","CIU", names(dat_ggplot6)[20:25], fixed=T)
names(dat_ggplot6)[20:25] <- gsub("Min","@", names(dat_ggplot6)[20:25])
names(dat_ggplot6)[20:25] <- gsub("Max","Min", names(dat_ggplot6)[20:25])
names(dat_ggplot6)[20:25] <- gsub("@","Max", names(dat_ggplot6)[20:25], fixed=T)

dat_ggplot6$PI <- dat3$PI_naive_prod_avg_late[match(paste(dat_ggplot6$Gene_name,dat_ggplot6$Batch, sep='.'), paste(dat3$Name,dat3$Batch,sep = '.'))]
dat_ggplot6 <- dat_ggplot6[apply(dat_ggplot6[,c("CIL_0.95_Min_kdeg", "CIU_0.95_Max_kdeg", "CIL_0.95_Min_transport", "CIU_0.95_Max_transport")],1,function(x){all(is.na(x))}),]

cor_test <- cor.test(dat_ggplot6$transport, log10(dat_ggplot6$HL), use="pairwise.complete.obs")
# svg(filename = 'Figures/Figure5D_transport_vs_HL.svg', width=4,height=4)
p <- ggplot(dat_ggplot6, aes(y=transport, x=HL)) + geom_point(mapping=aes(shape=Batch, color=PI))
p <- p + theme_bw() + theme(aspect.ratio=1) + scale_x_continuous(trans='log10', limits = c(3,1000)) + scale_y_continuous(limits=c(-2.6,0.5))
p <- p + geom_segment(aes(x=CIL_0.95_Max_HL,xend=CIU_0.95_Min_HL, y= transport, yend=transport), alpha=0.05)
p <- p + geom_segment(aes(x=CIL_0.95_Min_HL,xend=CIL_0.95_Max_HL, y= transport, yend=transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=CIU_0.95_Min_HL,xend=CIU_0.95_Max_HL, y= transport, yend=transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=HL,xend=HL, y= CIL_0.95_Max_transport, yend=CIU_0.95_Min_transport), alpha=0.05)
p <- p + geom_segment(aes(x=HL,xend=HL, y= CIL_0.95_Min_transport, yend=CIL_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=HL,xend=HL, y= CIU_0.95_Min_transport, yend=CIU_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + scale_color_continuous(trans='log10')
p <- p + scale_shape(labels = c("Rep1","Rep2")) + coord_cartesian(ylim=c(-2.5,0.5))
p <- p + labs(y=parse(text = "`Relative Transport `(log[10])"), x='HL (min)', shape="Replicate")
# p <- p + annotate(geom = "text", x=1000, y=0.5, hjust=1, vjust=1,label=bquote("Pearson's" ~ rho == .(signif(cor_test$estimate,2)) ~ ', p-value' == .(ifelse(cor_test$p.value < 2.2e-16, "< 2.2e-16", signif(cor_test$p.value,2)))))
p <- p + annotate(geom = "text", x=1000, y=0.5, hjust=1, vjust=1,label=bquote("Pearson's" ~ rho == .(signif(cor_test$estimate,2))))
p + theme(legend.position = "none")
# dev.off()
stop()

# svg(filename = 'Figures/Figure5C_legend.svg', width=5, height=5)
grid.draw(gglegend(p))
# dev.off()

#### Figure 7A - Responsiveness vs HL ####-------------------------------
dat_points <- c()
for (b in c("rep1", "rep2")){
  dat_points <- rbind(dat_points,
                      data.frame(k1prime = unlist(all_best_param$naive[[b]][,"k1'"]), 
                                 k1prime_CIL = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'", "lower"]})),
                                 k1prime_CILmin = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'", "lower_min"]})), 
                                 k1prime_CIU = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'", "upper"]})),
                                 k1prime_CIUmax = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'", "upper_max"]})),
                                 k2 = unlist(all_best_param$naive[[b]][,"k2"]), 
                                 k2_CIL = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2", "lower"]})),
                                 k2_CILmin = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2", "lower_min"]})), 
                                 k2_CIU = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2", "upper"]})),
                                 k2_CIUmax = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2", "upper_max"]})),
                                 k2prime = unlist(all_best_param$naive[[b]][,"k2'"]), 
                                 k2prime_CIL = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2'", "lower"]})),
                                 k2prime_CILmin = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2'", "lower_min"]})), 
                                 k2prime_CIU = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2'", "upper"]})),
                                 k2prime_CIUmax = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k2'", "upper_max"]})),
                                 transport = unlist(all_best_param$naive[[b]][,"k1'k2'/k2"]), 
                                 transport_CIL = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'k2'/k2", "lower"]})),
                                 transport_CILmin = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'k2'/k2", "lower_min"]})), 
                                 transport_CIU = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'k2'/k2", "upper"]})),
                                 transport_CIUmax = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'k2'/k2", "upper_max"]})),
                                 kdeg = unlist(all_best_param$naive[[b]][,"kdeg"]), 
                                 kdeg_CIL = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["kdeg", "lower"]})),
                                 kdeg_CILmin = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["kdeg", "lower_min"]})), 
                                 kdeg_CIU = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["kdeg", "upper"]})),
                                 kdeg_CIUmax = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["kdeg", "upper_max"]})),
                                 k1prime_over_k2 = unlist(all_best_param$naive[[b]][,"k1'/k2"]), 
                                 k1prime_over_k2_CIL = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'/k2", "lower"]})),
                                 k1prime_over_k2_CILmin = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'/k2", "lower_min"]})), 
                                 k1prime_over_k2_CIU = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'/k2", "upper"]})),
                                 k1prime_over_k2_CIUmax = unlist(lapply(all_param_CI$naive[[b]],function(x){x$CI["k1'/k2", "upper_max"]})),
                                 Batch=b))
}



HL <- c(6,60,300)
kdeg <- log10(log(2)/HL)
k2 <- c(-1,-0.5,0)
k1prime <- k2 + mean(dat_points$k1prime_over_k2)
transport <- c(-1.5, -1, -0.5)
k2prime <- transport - mean(dat_points$k1prime_over_k2)
x0 <- 1
alpha <- 10

dat_sim <- t_half_max <- c()
for (i in 1:length(kdeg)){
  t_half <- c()
  for (j in 1:length(k2)){
    for ( k in 1:length(transport)){
      z0 <- 10^(k1prime[j] + k2prime[k] - k2[j] - kdeg[i]) * x0;
      znorm = c(rep(1,31), (1 - alpha)/(10^kdeg[i]-10^k2[j])*(10^kdeg[i]*exp(-10^k2[j]*seq(1,1000)) - 10^k2[j]*exp(-10^kdeg[i]*seq(1,1000))) + alpha)
      dat_sim <- rbind(dat_sim, 
                       data.frame(t=seq(-30,1000),
                                  x=c(rep(x0,31),rep(alpha * x0,1000)), 
                                  z0 = z0,
                                  znorm = znorm,
                                  z= z0 * znorm, 
                                  kdeg=kdeg[i], HL=HL[i],
                                  k2 = k2[j], 
                                  transport = transport[k]
                       )
      )
    }
    f_t <- function(t){
      if(kdeg[i] != k2[j]){
        res <- 10^kdeg[i]*exp(-10^k2[j]*t)-10^k2[j]*exp(-10^kdeg[i]*t)-1/2*(10^kdeg[i] - 10^k2[j])
      }else{
        res <- (1-10^kdeg[i]*t)*exp(-10^kdeg[i]*t)-1/2
      }
      return(res)
    };
    t_half <- c(t_half, uniroot(f_t,c(0,1000000))$root)
  }
  t_half_max <- rbind(t_half_max, data.frame(t = range(t_half),
                                             kdeg=kdeg[i], HL=HL[i]))
}

dat_sim$kdeg <- as.factor(dat_sim$kdeg)
dat_sim$HL <- as.factor(dat_sim$HL)
dat_sim$transport <- as.factor(dat_sim$transport)

t_half_max$ymin <- 0.8
t_half_max$ymax <- 1/2*(alpha+1)

ribbon_color <- c("red", "orange", "chartreuse4")


# svg("Figures/Fig7A_responsiveness_vs_HL.svg", width=5, height=5)
p <- ggplot(data = dat_sim[dat_sim$transport == -1, ], aes(x=t, y=znorm, group=interaction(k2, HL, transport), col=HL))
p <- p + geom_line(aes(y=znorm))
p <- p + theme_bw() + theme(aspect.ratio=1)
p <- p + scale_color_manual(values=c("6"="red","60"="orange", "300"="chartreuse4"))
p <- p + geom_segment(x=-Inf, xend=Inf, y=1/2*(alpha+1),yend=1/2*(alpha+1),linetype=2, col="black")
for(i in 1:length(HL)){
  p <- p + geom_ribbon(data = t_half_max[t_half_max$HL==HL[i],], mapping = aes(x=t, ymin=ymin, ymax=ymax), fill = ribbon_color[i], alpha=0.3,inherit.aes = F)
}
p <- p + labs(x="Time (min)", y=parse(text="z/z[0]"), color="HL (min)")
p <- p + geom_segment(y=1, yend=1/2*(alpha+1), x=-5, xend=-5, arrow = arrow(ends = "both", length = unit(0.025, "npc")), size=0.1, col="black")
p <- p + annotate(geom = "text" , x=-25, y=1/2*(1/2*(alpha+1)-1)+ 1, label=parse(text="Delta/2"), size=3)
p <- p + geom_segment(y=1/2*(alpha+1),yend=alpha, x=-5, xend=-5, arrow = arrow(ends = "both", length = unit(0.025, "npc")), size=0.1, col="black")
p <- p + annotate(geom = "text" , x=-25, y=alpha-1/2*(1/2*(alpha+1)-1), label=parse(text="Delta/2"), size=3)
p <- p + coord_cartesian(xlim=c(-30,480))
p
# dev.off()
stop()

#### Figure 7B - Heatmap Responsiveness vs k2 + points ####-------------------------------
kdeg <- seq(-3.5,0.5,by=0.05)
HL <- log(2)/10^kdeg
k2 <- rev(seq(-3.5,0.5,by=0.05))
k1prime <- -0.2
k2prime <- -1.5

responsiveness <- matrix(nrow=length(k2), ncol=length(kdeg))
for(j in 1:length(kdeg)){
  for(i in 1:length(k2)){
    f_tmax <- function(t){
      if( 10^(k2[i]-kdeg[j]) != 1){
        # res <- 10^kdeg[j]*exp(-10^k2[i]*t)*(1-10^(k2[i]-kdeg[j])*exp(-10^(kdeg[j]-k2[i])*t))-alpha/(2*(alpha-1))*10^kdeg[j]*(1 - 10^(k2[i]-kdeg[j]))
        res <- 10^kdeg[j]*exp(-10^k2[i]*t) - 10^k2[i]*exp(-10^kdeg[j]*t)-1/2*(10^kdeg[j] - 10^k2[i])
        
      }else{
        res <- (1+10^kdeg[j]*t)*exp(-10^kdeg[j]*t)-1/2
      }
      return(res)
    }
    responsiveness[i,j] <- uniroot(f_tmax,c(0,100000))$root
  }
}

dat_responsiveness <- data.frame(responsiveness = as.vector(responsiveness), 
                                 kdeg = rep(kdeg,each=length(k2)), 
                                 HL = rep(HL,each=length(k2)), 
                                 k1prime = k1prime, 
                                 k2 = rep(k2,length(kdeg)), 
                                 k2prime = k2prime, which="k2")

# plot the CI of k2 vs kdeg to show that the limits are higher than kdeg

hl_trans <- function(kdeg){log(2)/10^kdeg}

dat_points2 <- dat_points[apply(dat_points[,c("k2_CILmin","k2_CIUmax","kdeg_CILmin", "kdeg_CIUmax")],1, function(x){all(is.na(x))}) , ]

# svg('Fig7B_responsiveness_heatmap kdeg vs k2_v2.svg', width=5, height=5)
p <- ggplot(dat_responsiveness[dat_responsiveness$which == "k2",],aes(x=kdeg,y=k2)) + geom_tile(aes(fill=responsiveness)) 
p <- p + scale_fill_gradientn(trans="log10", colours = hm.palette(100))
p <- p + scale_shape_discrete(labels=c("Rep1","Rep2"))
p <- p + coord_cartesian(xlim = c(-3.5,0.5), ylim=c(-3.5,0.5), expand = F)
p <- p + theme_bw() + labs(fill = parse(text="t[1/2]~(min)"), y=parse(text="k[2]~(log[10])"), x=parse(text="k[deg]~(log[10])"), shape = "Replicate")
p <- p + geom_point(data=dat_points2, aes(shape=Batch), alpha=0.3)
p <- p + geom_segment(data=dat_points2, aes(x=kdeg, xend=kdeg, y = k2_CIL, yend=k2_CIU), alpha=0.1)
p <- p + geom_segment(data=dat_points2, aes(x=kdeg, xend=kdeg, y = k2_CILmin, yend=k2_CIL), lty=2, alpha=0.1)
p <- p + geom_segment(data=dat_points2, aes(x=kdeg, xend=kdeg, y = k2_CIU, yend= k2_CIUmax), lty=2, alpha=0.1)
p <- p + geom_segment(data=dat_points2, aes(x=kdeg_CIL, xend=kdeg_CIU, y = k2, yend= k2), alpha=0.1)
p <- p + geom_segment(data=dat_points2, aes(x=kdeg_CILmin, xend=kdeg_CIL, y = k2, yend=k2), lty=2, alpha=0.1)
p <- p + geom_segment(data=dat_points2, aes(x=kdeg_CIU, xend=kdeg_CIUmax, y = k2, yend=k2), lty=2, alpha=0.1)
p <- p + scale_x_continuous(sec.axis = sec_axis(trans = hl_trans, name= "HL (min)", breaks = c(600, 300, 120, 60, 30, 12, 6, 3,1)))
p + theme(aspect.ratio = 1, panel.border = element_blank())
# dev.off()
stop()

#### Figure 7C - Responsiveness vs HL and transport ####
dat_points$HL <- log(2)/10^dat_points$kdeg

dat_points2 <- dat_points[apply(dat_points[,c("k2_CILmin","k2_CIUmax","kdeg_CILmin", "kdeg_CIUmax")],1, function(x){all(is.na(x))}) , ]
dat_points2$HL_CIL <-  log(2)/10^dat_points2$kdeg_CIU
dat_points2$HL_CIU <-  log(2)/10^dat_points2$kdeg_CIL

dat_points2$thalfmax <- apply(dat_points2, 1, function(x){
  kdeg <- as.numeric(x["kdeg"])
  k2 <- as.numeric(x["k2"])
  
  f_t <- function(t){
    if(kdeg != k2){
      res <- 10^kdeg*exp(-10^k2*t)-10^k2*exp(-10^kdeg*t)-1/2*(10^kdeg - 10^k2)
    }else{
      res <- (1-10^kdeg*t)*exp(-10^kdeg*t)-1/2
    }
    return(res)
  }
  
  thalfmax <- uniroot(f_t,c(0,1000000))$root
})

dat_points2$thalfmax_max <- apply(dat_points2, 1, function(x){
  kdeg_max <- as.numeric(x["kdeg_CIU"])
  k2_max <- as.numeric(x["k2_CIU"])
  
  f_tmax <- function(t){
    if(kdeg_max != k2_max){
      res <- 10^kdeg_max*exp(-10^k2_max*t)-10^k2_max*exp(-10^kdeg_max*t)-1/2*(10^kdeg_max - 10^k2_max)
    }else{
      res <- (1-10^kdeg_max*t)*exp(-10^kdeg_max*t)-1/2
    }
    return(res)
  }
  thalfmax_max <- uniroot(f_tmax,c(0,1000000))$root
})

dat_points2$thalfmax_min <- apply(dat_points2, 1, function(x){
  kdeg_min <- as.numeric(x["kdeg_CIL"])
  k2_min <- as.numeric(x["k2_CIL"])
  
  f_tmin <- function(t){
    if(kdeg_min != k2_min){
      res <- 10^kdeg_min*exp(-10^k2_min*t)-10^k2_min*exp(-10^kdeg_min*t)-1/2*(10^kdeg_min - 10^k2_min)
    }else{
      res <- (1-10^kdeg_min*t)*exp(-10^kdeg_min*t)-1/2
    }
    return(res)
  }
  thalfmax_min <- uniroot(f_tmin,c(0,1000000))$root
})


p <- ggplot(dat_points2, aes(x=HL, y=thalfmax, shape=Batch)) + geom_point()
p <- p + coord_cartesian(xlim=c(0,400), ylim=c(0,400))
p <- p + theme_bw() + theme(aspect.ratio = 1)
p <- p + geom_segment(aes(x=HL_CIL, xend=HL_CIU, y=thalfmax, yend=thalfmax), alpha=0.05)
p <- p + geom_segment(aes(x=HL, xend=HL, y=thalfmax_min, yend=thalfmax_max), alpha=0.05)
p1 <- p + labs(x="HL (min)", y="Time to half max induction (min)")

p <- ggplot(dat_points2, aes(x=transport, y=thalfmax, shape=Batch)) + geom_point()
p <- p + coord_cartesian(ylim=c(0,400), xlim=c(-2.5,0.5))
p <- p + theme_bw() + theme(aspect.ratio = 1)
p <- p + geom_segment(aes(x=transport_CIL, xend=transport_CIU, y=thalfmax, yend=thalfmax), alpha=0.05)
p <- p + geom_segment(aes(x=transport, xend=transport, y=thalfmax_min, yend=thalfmax_max), alpha=0.05)
p2 <- p + labs(x="Relative Transport", y="Time to half max induction (min)")

p_legend <- gglegend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

# svg(filename = "Figures/Fig7C_Responsiveness of genes_v2.svg", width=10, height=5)
gridExtra::grid.arrange(p1,p2,p_legend, ncol=3, widths=c(5,5,1))
# dev.off()

stop()

#### Figure 7D - Correlation Transport HL ####------
# svg(filename = 'Figures/Fig7D_transport_vs_decay_v2.svg', width=5,height=5)
p <- ggplot(dat_ggplot6, aes(y=transport, x=kdeg)) + geom_point(mapping=aes(shape=Batch))
p <- p + theme_bw() + theme(aspect.ratio=1) + scale_x_continuous(limits = c(-3,-0.5), sec.axis = sec_axis(trans = hl_trans, name= "HL (min)", breaks = c(600, 300, 120, 60, 30, 12, 6, 3,1))) + scale_y_continuous(limits=c(-2.6,0.5))
p <- p + geom_segment(aes(x=CIL_0.95_Max_kdeg,xend=CIU_0.95_Min_kdeg, y= transport, yend=transport), alpha=0.05)
p <- p + geom_segment(aes(x=CIL_0.95_Min_kdeg,xend=CIL_0.95_Max_kdeg, y= transport, yend=transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=CIU_0.95_Min_kdeg,xend=CIU_0.95_Max_kdeg, y= transport, yend=transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=kdeg,xend=kdeg, y= CIL_0.95_Max_transport, yend=CIU_0.95_Min_transport), alpha=0.05)
p <- p + geom_segment(aes(x=kdeg,xend=kdeg, y= CIL_0.95_Min_transport, yend=CIL_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + geom_segment(aes(x=kdeg,xend=kdeg, y= CIU_0.95_Min_transport, yend=CIU_0.95_Max_transport), alpha=0.05, linetype=2)
p <- p + scale_shape(labels = c("Rep1","Rep2")) + coord_cartesian(ylim=c(-2.5,0.5))
p <- p + labs(y=parse(text = "`Relative Transport `(log[10])"), x=parse(text="k[deg]~(log[10])"), shape="Replicate")
# p <- p + annotate(geom = "text", x=1000, y=0.5, hjust=1, vjust=1,label=bquote("Pearson's" ~ rho == .(signif(cor_test$estimate,2)) ~ ', p-value' == .(ifelse(cor_test$p.value < 2.2e-16, "< 2.2e-16", signif(cor_test$p.value,2)))))
# p + theme(legend.position = "none")
p
# dev.off()

#### Figure 7E - Heatmap steady transport vs kdeg ####----------------
kdeg <- transport <-seq(-3.5,0.5,by=0.05)
transport <- rev(transport)

x0 <- 1 
alpha <- 10

zinf <- log2(10^transport%*%10^(-t(kdeg)))
dat4 <- data.frame(value = as.vector(zinf), kdeg=rep(kdeg, each=length(transport)), transport=rep(transport,length(kdeg)))

hm.palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'Spectral')))

# svg('Fig7E_sim_ss_transport_vs_kdeg.svg',height=5,width=5)
p <- ggplot(dat4, aes(x=kdeg, y=transport)) + geom_tile(aes(fill=value)) 
p <- p + scale_fill_gradientn(colours = hm.palette(100))
p <- p + coord_cartesian(xlim = c(-3.5,0.5), ylim=c(-3.5,0.5), expand = F)
p <- p + theme_bw() + labs(fill = parse(text="log[2](z[infinity]/x[infinity])"), y=parse(text="`Relative Transport `(log[10])"), x=parse(text="k[deg]~(log[10])"))
p <- p + geom_segment(y=-3, yend=-3, x=-2.5, xend=-0.5, arrow = arrow(length = unit(0.025, "npc")))
p <- p + annotate("text", x=-1.5, y=-3.1, label = "increase responsiveness")
p <- p + geom_segment(y=-3, yend=-1, x=-0.5, xend=-0.5, arrow = arrow(length = unit(0.025, "npc")))
p <- p + annotate("text", y=-2, x=-0.25, label = "compensate for\ndecrease in level", angle=90)
p <- p + scale_x_continuous(sec.axis = sec_axis(trans = hl_trans, name= "HL (min)", breaks = c(600, 300, 120, 60, 30, 12, 6, 3,1)))
p + theme(aspect.ratio = 1, panel.border = element_blank())
# dev.off()
stop()

#### Supplementary pdfs - Fit naive and lpa ####------------------------------------
height <- 2*3
width <- 2*8

dir.create('Supplementary Documents',showWarnings = F)
for (condition in names(all_data)){
  pdf(file = paste0('Supplementary Documents/all_genes_optim_',condition, '.pdf'), width = width, height=height)
  
  for (g in names(all_data[[condition]])){
    print(paste(g ," - ", condition))
    # get data
    dat_g <- data.frame(RPKM=c(), RPKM_ub = c(), RPKM_lb=c(), Time = c(), batch = c(), Cpt = c(), Type=c(), Condition=c())
    
    for (Cpt in names(all_data[[condition]][[g]])){
      for (b in colnames(all_data[[condition]][[g]][[Cpt]])[-1]){
        dat_g <- rbind(dat_g, data.frame(RPKM = all_data[[condition]][[g]][[Cpt]][,b], 
                                         RPKM_ub = NA, RPKM_lb = NA, 
                                         Time = all_data[[condition]][[g]][[Cpt]][,1], 
                                         batch = b,
                                         Cpt = Cpt,
                                         Type = "Experiment",
                                         Condition = condition ))
      }
    }
    
    # add fit
    for (batch in names(all_fit[[condition]])){
      for (Cpt in names(all_fit[[condition]][[batch]][[g]])){
        for (b in colnames(all_fit[[condition]][[batch]][[g]][[Cpt]])[-1]){
          dat_g <- rbind(dat_g, data.frame(RPKM = all_fit[[condition]][[batch]][[g]][[Cpt]][,b], 
                                           RPKM_ub = all_fit_CI[['UB']][[condition]][[batch]][[g]][[Cpt]][match(all_fit[[condition]][[batch]][[g]][[Cpt]][,1], all_fit_CI[['UB']][[condition]][[batch]][[g]][[Cpt]][,1]),b], 
                                           RPKM_lb = all_fit_CI[['LB']][[condition]][[batch]][[g]][[Cpt]][match(all_fit[[condition]][[batch]][[g]][[Cpt]][,1], all_fit_CI[['LB']][[condition]][[batch]][[g]][[Cpt]][,1]),b], 
                                           Time = all_fit[[condition]][[batch]][[g]][[Cpt]][,1], 
                                           batch = b,
                                           Cpt = Cpt,
                                           Type = ifelse(batch=="all","together","separate"),
                                           Condition = condition ))
        }
      }
    }
    
    # make SSE_table
    SSE_table <- as.data.frame(matrix(nrow=3*3+2, ncol=0))
    row.names(SSE_table) <- c(paste("-NLL",rep(c("rep1", "rep2","rep3"), each=3), rep(c("ca", "np","cyto"), 3), sep=" "), "Total", "Total with regul")
    for (batch in names(all_annot[[condition]])){
      if(!is.null(all_annot[[condition]][[batch]][[g]])){
          
        if(batch == "all"){
          tmp <- unlist(all_annot[[condition]][[batch]][[g]]$annot_bs[grep("NLL", row.names(all_annot[[condition]][[batch]][[g]]$annot_bs)),])
          SSE_table <- cbind(SSE_table, c(tmp, sum(tmp, na.rm=T), all_annot[[condition]][[batch]][[g]]$annot["NLL",]))
          colnames(SSE_table)[ncol(SSE_table)] <- "Together" 
        }else{
          tmp <- unlist(all_annot[[condition]][[batch]][[g]]$annot_bs[grep("NLL", row.names(all_annot[[condition]][[batch]][[g]]$annot_bs)),])
          SSE_table <- cbind(SSE_table, NA)
          colnames(SSE_table)[ncol(SSE_table)] <- batch 
          SSE_table[grep(batch, row.names(SSE_table)),batch] <- tmp
          SSE_table[grep("Total", row.names(SSE_table)), batch] <- c(sum(tmp, na.rm=T), all_annot[[condition]][[batch]][[g]]$annot["NLL",])
        }
      }
    }
    
    # make parameter_table
    par_table <- as.data.frame(matrix(nrow=3*3+2, ncol=0))
    row.names(par_table) <- c("spar", "sigmab", "sigmat", "ca0_rep1", "k1'", "k2", "k2'", "kdeg", "k1'k2'", "k1'/k2", "k1'k2'/k2")
    for (batch in names(all_best_param[[condition]])){
      tmp <- unlist(all_best_param[[condition]][[batch]][g,])
      par_table <- cbind(par_table, signif(tmp[match(row.names(par_table), names(tmp))],4))
      colnames(par_table)[ncol(par_table)] <- ifelse(batch == "all", "Together", batch)
      
      if(batch == 'rep3'){
        par_table['k2', batch] <- paste(signif(tmp["max(k2,kdeg)"],4), " or ", signif(tmp["min(k2,kdeg)"],4))
        par_table['kdeg', batch] <- paste(signif(tmp["min(k2,kdeg)"],4), " or ", signif(tmp["max(k2,kdeg)"],4))
        par_table["k1'k2'/k2", batch] <- paste(signif(tmp["k1'k2'"] - tmp["max(k2,kdeg)"],4), " or ", signif(tmp["k1'k2'"]  - tmp["min(k2,kdeg)"],4))
      }else{
        par_table["k1'k2'", ifelse(batch == "all", "Together", batch)] <- signif(tmp["k1'"]+tmp["k2'"],4)
      }
    }
    row.names(par_table) <-  c("spar", "sigma[b]", "sigma[t]", "ca[paste(0,',', rep[1])]", "log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](paste(k[1]*minute,k[2]*minute))", "log[10](k[1]*minute/k[2])", "transport == log[10](paste(k[1]*minute, k[2]*minute)/k[2])")
    
    dat_g$Cpt <- factor(dat_g$Cpt, levels = c("caRNA", "npRNA", "cytoRNA"), labels = c("Chromatin", "Nucleoplasm", "Cytoplasm"))
    dat_g$Type <- factor(dat_g$Type, levels = c("Experiment", "together", "separate"), labels = c("Experiment", "Fit - All together", "Fit - Separately"))
    dat_g$batch <- factor(dat_g$batch, levels = c("rep1","rep2","rep3"), labels = c("Rep1", "Rep2", "Rep3"))
    
    p1 <- ggplot(dat_g[dat_g$Type != "Experiment",], mapping = aes(x = Time, y = RPKM, col=batch)) + geom_line(data = dat_g[dat_g$Type != "Experiment",], aes(linetype=Type, group=interaction(Type, batch))) 
    p1 <- p1 + facet_grid(Type ~ Cpt, drop = T)
    p1 <- p1 + geom_errorbar(mapping=aes(x=Time,ymin = RPKM_lb, ymax=RPKM_ub, color = batch, group=interaction(Type, batch)))
    p1 <- p1 + geom_point(data = dat_g[dat_g$Type == "Experiment",-which(colnames(dat_g)=="Type")])
    p1 <- p1 + theme_bw() + theme(aspect.ratio=1) 
    p1 <- p1 + labs(title=paste0(g), color="Replicate", y="log2(FPKM)") 
    p1 <- p1 + theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = rel(0.8))) + scale_x_continuous(breaks = c(0,30,60,90,120))
    
    # change color of the strips
    g1 <- ggplotGrob(p1)
    col <- c("Chromatin" = col_ca, "Nucleoplasm" = col_np, "Cytoplasm" = col_cyto)
    for(i in grep("strip-t",g1$layout$name)){
      Cpt <- g1$grobs[[i]]$grobs[[1]]$children[[grep("text", names(g1$grobs[[i]]$grobs[[1]]$children))]]$children[[1]]$label
      g1$grobs[[i]]$grobs[[1]]$children[[grep("rect", names(g1$grobs[[i]]$grobs[[1]]$children))]]$gp$fill <- col[Cpt]
    }
    
    tt <- ttheme_default(base_size=6, rowhead=list(fg_params=list(parse=T, fontface = 2)), colhead=list(fg_params=list(fontface = 2)))
    tbl1 <- tableGrob(replace(as.matrix(signif(SSE_table,4)), is.na(SSE_table), ""), theme = tt)
    tbl2 <- tableGrob(replace(as.matrix(par_table), is.na(par_table), ""), theme = tt)
    
    p2 <- grid.arrange(g1, tbl1, tbl2, widths = width*c(4/8,2/8,2/8), heights = height, layout_matrix=matrix(c(1,2,3), ncol=3, byrow = T))
    grid.draw(p2)
  }
  dev.off()
}

stop()

#### Supplementary Pdfs - Profile Analysis ####---------------------------------
for (condition in c('naive', 'lpa')){
  pdf(file = paste0('Supplementary Documents/Profile_dMod_merged_', condition,'_with_optim_points.pdf'), width=19)
  
  for (g in names(all_param_CI[[condition]][[1]])){
    print(paste(g, " - ", condition))
    dat_ggplot_points <- dat_ggplot <- dat_CI <- dat_blank <- c()
    for (batch in names(all_param_CI[[condition]])){
      # add column for new param
      below_th <- all_param[[condition]][[batch]][[g]][all_param[[condition]][[batch]][[g]][,"value"] <= min(all_param[[condition]][[batch]][[g]][,"value"])+ qchisq(0.99,1,lower.tail = T)/2,]
      dat_ggplot_points <- rbind(dat_ggplot_points, data.frame(whichPar = rep(colnames(below_th)[-ncol(below_th)], each=nrow(below_th)), par_val = as.vector(below_th[,-ncol(below_th)]), value = rep(below_th[,"value"],ncol(below_th)-1), batch = ifelse(batch=="all", "Together",batch)))
      for( par in names(all_param_CI[[condition]][[batch]][[g]]$Profile) ){
        dat_ggplot <- rbind(dat_ggplot, data.frame(all_param_CI[[condition]][[batch]][[g]]$Profile[[par]][,c(1,2,5)], par_val=all_param_CI[[condition]][[batch]][[g]]$Profile[[par]][,unique(all_param_CI[[condition]][[batch]][[g]]$Profile[[par]]$whichPar)], method = all_param_CI[[condition]][[batch]][[g]]$Profile[[par]][,"method"],batch = ifelse(batch=="all", "Together",batch)))
      }
      
      dat_CI <- rbind(dat_CI, data.frame(all_param_CI[[condition]][[batch]][[g]]$CI, batch=ifelse(batch=="all", "Together",batch)))
      dat_blank <- rbind(dat_blank, 
                         data.frame(whichPar = all_param_CI[[condition]][[batch]][[g]]$CI$whichPar , x= all_param_CI[[condition]][[batch]][[g]]$CI$lower, y = min(below_th[,"value"]),batch = ifelse(batch=="all", "Together",batch)), 
                         data.frame(whichPar = all_param_CI[[condition]][[batch]][[g]]$CI$whichPar , x= all_param_CI[[condition]][[batch]][[g]]$CI$upper, y = min(below_th[,"value"])+ qchisq(0.99,1,lower.tail = T)/2,batch = ifelse(batch=="all", "Together",batch))) 
    }
    
    dat_ggplot$value <- dat_ggplot$value/2
    
    dat_ggplot$method <- factor(dat_ggplot$method, levels=c("approximate", "exact", "optim"))
    dat_CI$method_lower <- factor(dat_CI$method_lower, levels=c("approximate", "exact", "optim"))
    dat_CI$method_upper <- factor(dat_CI$method_upper, levels=c("approximate", "exact", "optim"))
    
    dat_ggplot$whichParOriginal <- dat_ggplot$whichPar
    dat_ggplot$whichPar[dat_ggplot$whichPar %in% c("k2", "kdeg") & dat_ggplot$batch == "rep3"] <- "k2 or kdeg"
    dat_ggplot$whichPar <- factor(dat_ggplot$whichPar, 
                                  levels = c("k1'", "k2", "k2'", "kdeg", "k1'k2'", "k2 or kdeg", "k1'/k2", "k1'k2'/k2", "k2'/k2"),
                                  labels = c("log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](paste(k[1]*minute, k[2]*minute))", "log[10](paste(k[2], ' or ', k[deg]))", "log[10](paste(k[1]*minute/k[2]))", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])", "log[10](k[2]*minute/k[2])"))
    
    dat_CI$whichPar <- as.character(dat_CI$whichPar)
    dat_CI$whichPar[dat_CI$whichPar %in% c("min(k2,kdeg)", "max(k2,kdeg)") & dat_CI$batch == "rep3"] <- "k2 or kdeg"
    dat_CI$whichPar <- factor(dat_CI$whichPar, 
                              levels = c("k1'", "k2", "k2'", "kdeg", "k1'k2'", "k2 or kdeg", "k1'/k2", "k1'k2'/k2", "k2'/k2"),
                              labels = c("log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](paste(k[1]*minute, k[2]*minute))", "log[10](paste(k[2], ' or ', k[deg]))", "log[10](paste(k[1]*minute/k[2]))", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])", "log[10](k[2]*minute/k[2])"))
    
    dat_ggplot_points$whichPar <- as.character(dat_ggplot_points$whichPar)
    dat_ggplot_points$whichPar[dat_ggplot_points$whichPar %in% c("min(k2,kdeg)", "max(k2,kdeg)") & dat_ggplot_points$batch == "rep3"] <- "k2 or kdeg"
    dat_ggplot_points$whichPar <- factor(dat_ggplot_points$whichPar, 
                                         levels = c("k1'", "k2", "k2'","kdeg", "k1'k2'", "k2 or kdeg", "k1'/k2", "k1'k2'/k2", "k2'/k2"),
                                         labels = c("log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](paste(k[1]*minute, k[2]*minute))", "log[10](paste(k[2], ' or ', k[deg]))", "log[10](paste(k[1]*minute/k[2]))", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])", "log[10](k[2]*minute/k[2])"))
    
    dat_blank$whichPar <- as.character(dat_blank$whichPar)
    dat_blank$whichPar[dat_blank$whichPar %in% c("min(k2,kdeg)", "max(k2,kdeg)") & dat_blank$batch == "rep3"] <- "k2 or kdeg"
    dat_blank$whichPar <- factor(dat_blank$whichPar, 
                                 levels = c("k1'", "k2", "k2'","kdeg", "k1'k2'", "k2 or kdeg", "k1'/k2", "k1'k2'/k2", "k2'/k2"),
                                 labels = c("log[10](k[1]*minute)", "log[10](k[2])", "log[10](k[2]*minute)", "log[10](k[deg])", "log[10](paste(k[1]*minute, k[2]*minute))", "log[10](paste(k[2], ' or ', k[deg]))", "log[10](paste(k[1]*minute/k[2]))", "log[10](paste(k[1]*minute,k[2]*minute)/k[2])", "log[10](k[2]*minute/k[2])"))
    
    # dat_blank <- dat_ggplot[dat_ggplot$constraint == 0,]
    # dat_blank <- rbind(cbind(dat_blank, y = dat_blank$value, x = apply(dat_blank, 1, function(x){min(dat_ggplot$par_val[dat_ggplot$whichPar == x[3] & dat_ggplot$batch == x[6] & dat_ggplot$method == x[5]])})) , # min
    #                    cbind(dat_blank, y = dat_blank$value + qchisq(0.99,1,lower.tail = T)/2, x = apply(dat_blank, 1, function(x){max(dat_ggplot$par_val[dat_ggplot$whichPar == x[3] & dat_ggplot$batch == x[6] & dat_ggplot$method == x[5]])}))) # max
    
    # dat_blank <- rbind(do.call(rbind, lapply(unique()data.frame(dat_CI, x=dat_CI$lower, y=min(dat_ggplot_points$value)),
    #                    data.frame(dat_CI, x=dat_CI$upper, y=min(below_th[,"value"]) + qchisq(0.99,1,lower.tail = T)/2))
    # 
    p_blank <- ggplot() + geom_blank(data=dat_blank, mapping = aes(y=y, x=x))
    p_blank <- p_blank + facet_grid( batch~whichPar , scales = 'free', labeller = "label_parsed")
    p_blank <- p_blank + geom_hline(data = do.call(rbind, lapply(unique(dat_blank$batch), function(b){tmp <- dat_blank[dat_blank$batch == b,]; cbind(tmp[1:(nrow(tmp)/2),],y_th=tmp[1:(nrow(tmp)/2),"y"]+qchisq(0.95,1,lower.tail = T)/2)})), mapping=aes(yintercept = y_th), linetype=2)
    p_blank_test <- ggplot_build(p_blank)
    
    p <- p_blank + theme_bw() + labs(title= g, y = 'NLL', x='Parameter value')
    p <- p + geom_vline(data=dat_CI, mapping=aes(xintercept=lower, col=method_lower), linetype=2, size=0.5, alpha=0.5)
    p <- p + geom_vline(data=dat_CI, mapping=aes(xintercept=upper, col=method_upper), linetype=2, size=0.5, alpha=0.5)
    if(any(!is.na(dat_CI$lower_min))){
      p <- p + geom_vline(data=dat_CI, mapping=aes(xintercept=lower_min, col=method_lower), linetype=2, size=0.5, alpha=0.5)
      p <- p + geom_rect(data=dat_CI, mapping=aes(xmin=lower_min, xmax=lower, ymin=-Inf, ymax=Inf, fill=method_lower), alpha=0.1)
    }
    if(any(!is.na(dat_CI$upper_max))){
      p <- p + geom_vline(data=dat_CI, mapping=aes(xintercept=upper_max, col=method_upper), linetype=2, size=0.5, alpha=0.5)
      p <- p + geom_rect(data=dat_CI, mapping=aes(xmin=upper, xmax=upper_max, ymin=-Inf, ymax=Inf, fill=method_upper), alpha=0.1)
    }
    p <- p + theme(aspect.ratio=1, axis.text = element_text(size = rel(0.5)))
    # p <- p + scale_x_continuous(limits = c(dat_blank$value, dat_blank$value + qchisq(0.99,1,lower.tail = T)/2))
    
    p <- p + geom_line(data = dat_ggplot, mapping = aes(col=method,x = par_val, y=value, group=interaction(whichParOriginal, method)))
    p <- p + geom_point(data = dat_ggplot[dat_ggplot$constraint == 0,], mapping = aes(col=method,x = par_val, y=value))
    
    p <- p + geom_point(data = dat_ggplot_points[interaction(dat_ggplot_points$whichPar, dat_ggplot_points$batch) %in% unique(interaction(dat_CI$whichPar, dat_CI$batch)),], mapping=aes(x=par_val, y=value), alpha=0.1)
    p <- p + scale_color_discrete(drop=FALSE) + scale_fill_discrete(drop=FALSE)
    
    p2 <- ggplot_build(p)
    
    p2$layout$panel_scales_y <- p_blank_test$layout$panel_scales_y
    p2$layout$setup_panel_params()
    
    # p2$layout$panel_params <- lapply( 1:length(p2$layout$panel_params), function(panel_idx){ p2$layout$panel_params[[panel_idx]]$y <- p_blank_test$layout$panel_params[[panel_idx]]$y; p2$layout$panel_params[[panel_idx]];})
    
    
    table_CI <- dat_CI
    table_CI$value <- as.character(signif(table_CI$value,4))
    # table_CI$whichPar <- paste0('log10(',table_CI$whichPar,')')
    table_CI$lower <- as.character(signif(table_CI$lower,4))
    table_CI$lower[is.infinite(table_CI$lower_min)] <- paste0('< ', table_CI$lower[is.infinite(table_CI$lower_min)] )
    table_CI$upper <- as.character(signif(table_CI$upper,4))
    table_CI$upper[is.infinite(table_CI$upper_max)] <- paste0('> ', table_CI$upper[is.infinite(table_CI$upper_max)] )
    
    table_CI <- table_CI[, c("batch", "whichPar", "value", "lower", "upper", "method_lower", "method_upper")]
    names(table_CI) <- c("Replicate", "Par","Best value", "CI95 LB", "CI95 UB", "Method LB", "Method UB")
    
    tt <- ttheme_default(base_size=6, parse=T, colhead=list(fg_params=list(fontface=2)))
    tbl1 <- tableGrob(table_CI, theme = tt,rows = NULL)
    
    p3 <- grid.arrange(ggplot_gtable(p2), tbl1, widths = c(14/19,5/19), heights = 1, layout_matrix=matrix(c(1,2), ncol=2, byrow = T))
    
    # grid.newpage()
    grid.draw(p3)
  }
  dev.off()
}


