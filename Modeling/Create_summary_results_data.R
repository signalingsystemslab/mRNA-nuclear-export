library(biomaRt)
library(compiler)
library(deSolve)
# library(grImport2)

# DEG slightly more stringent
list_genes_cleaned <- read.csv('../Manual Curation/gene_list.csv', sep=',')

# load Data on exons whole genes, ActD data
load(file = '../Data/exons_rpkms/caRNA_Naive_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/npRNA_Naive_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/cytoRNA_Naive_exons_rpkm_5kb.Rdata')

load(file = '../Data/exons_rpkms/caRNA_LPA_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/npRNA_LPA_exons_rpkm_5kb.Rdata')
load(file = '../Data/exons_rpkms/cytoRNA_LPA_exons_rpkm_5kb.Rdata')

# load gene_infos
load('../Data/gene_infos_with_length.Rdata')
load('../Data/raw_counts/lib.size.Rdata')

# load modeling function
source(file = 'functions_weighting_for_smoothing_include_negbinom_error_in_model.R')

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

all_param_naive <- all_param$naive
all_param_lpa <- all_param$lpa

save(all_param_naive, file = 'all_param_naive.Rdata')
save(all_param_lpa, file = 'all_param_lpa.Rdata')
save(all_annot, all_best_param, all_data, all_fit, all_fit_CI, all_library, all_param_CI, file = 'all_results.Rdata')
stop()
