## Get input arguments
print('Starting the Run')
args = commandArgs(TRUE)

print('Given Arguments:')
print(args)

i <- as.numeric(args[1])
condition <- as.character(args[2])
replicate <- as.character(args[3])
folder <- as.character(args[4])

print(paste0('i: ', i))
print(paste0('condition: ', condition))
print(paste0('replicate: ', replicate))
print(paste0('folder: ', folder))

## Load Packages
print('Loading packages')
for (pkg in c('compiler', 'deSolve', 'dMod', 'numDeriv')){
  if (!(pkg %in% installed.packages())){
    install.packages(pkg, repos = 'https://cran.cnr.berkeley.edu/')
  }
}

library(compiler)
library(deSolve)
library(dMod)
library(numDeriv)
source('functions_weighting_for_smoothing_include_negbinom_error_in_model.R')
source('profile_dMod_with_save.R')

# pre-compile hessian and grad function
grad_cp <- cmpfun(grad)
hessian_cp <- cmpfun(hessian)

TSE_reps_oldpar <- function(x,...){
  par_names <- names(x)
  pars_new <- c(x[ !(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg")) ], parameters_convert_cp(x[(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))], log_in = T, log_out = T, direction = 1, no_np = no_np))
  names(pars_new) <- c(par_names[!(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))], par_names[(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))])
  x_new <- pars_new[par_names]
  res <- TSE_reps_cp(par = x_new,...)
  return(res)
}
  
TSE_reps_cp_grad <- function(x, ...){grad_cp(TSE_reps_oldpar, x, ...)}
TSE_reps_cp_hessian <- function(x, ...){ hessian_cp(TSE_reps_oldpar, x, ...) }

################
## Load Rdata ##
################
print('Creating needed folders')

scratch_dir <- Sys.getenv("SCRATCH")
if(scratch_dir == ""){
  scratch_dir <- "."
}
print(scratch_dir)
out_dir <- paste0(scratch_dir,'/',folder)
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir, '/Profile_dMod_integrate')
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir,'/', condition)
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir,'/Run_log')
dir.create(file.path(out_dir), showWarnings = FALSE)

out_dir_path <- paste0(scratch_dir,'/',folder,'/Profile_dMod_integrate/', condition)

####################################
## Get data for all compartments ##
####################################
print('Loading data')

## Exons last 5kb ----
if (condition == 'naive'){
  load('../Data/exons_rpkms/caRNA_Naive_exons_rpkm_5kb.Rdata')
  load('../Data/exons_rpkms/npRNA_Naive_exons_rpkm_5kb.Rdata')
  load('../Data/exons_rpkms/cytoRNA_Naive_exons_rpkm_5kb.Rdata')
  g <- rownames(caRNA_Naive_exons_rpkm_5kb)[i]
}else if (condition == 'lpa'){
  load('../Data/exons_rpkms/caRNA_LPA_exons_rpkm_5kb.Rdata')
  load('../Data/exons_rpkms/npRNA_LPA_exons_rpkm_5kb.Rdata')
  load('../Data/exons_rpkms/cytoRNA_LPA_exons_rpkm_5kb.Rdata')
  g <- rownames(caRNA_LPA_exons_rpkm_5kb)[i]
}

load('../Data/gene_infos_with_length.Rdata')
gene_name <- gene_infos$external_gene_name[gene_infos$ensembl_gene_id == substr(g,1,18)]

load('../Data/raw_counts/lib.size.Rdata')
lib <- lib_exons[,c("lib.size", "norm.factors")]
row.names(lib) <- lib_exons$samples

# lib$Time <- sapply(row.names(lib),function(x){strsplit(x,'.', fixed=T)[[1]][3]})
# lib$Cpt <- sapply(row.names(lib),function(x){strsplit(x,'.', fixed=T)[[1]][1]})
# lib$replicate <- sapply(row.names(lib),function(x){strsplit(x,'.', fixed=T)[[1]][4]})

print(g)
print(gene_name)

if (condition == "naive"){
  data_rpkm_all <- list( caRNA   = as.matrix(create_dat_g_rpkm_cp(caRNA_Naive_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         npRNA   = as.matrix(create_dat_g_rpkm_cp(npRNA_Naive_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         cytoRNA = as.matrix(create_dat_g_rpkm_cp(cytoRNA_Naive_exons_rpkm_5kb, g, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
  length <- gene_infos$Length_5kb_naive[gene_infos$ensembl_gene_id == substr(g,1,18)]
}else{
  data_rpkm_all <- list( caRNA   = as.matrix(create_dat_g_rpkm_cp(caRNA_LPA_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         npRNA   = as.matrix(create_dat_g_rpkm_cp(npRNA_LPA_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         cytoRNA = as.matrix(create_dat_g_rpkm_cp(cytoRNA_LPA_exons_rpkm_5kb, g, replicate=c("rep1", "rep2", "rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
  length <- gene_infos$Length_5kb_lpa[gene_infos$ensembl_gene_id == substr(g,1,18)]
}

time_data <- sort(unique(unlist(lapply(data_rpkm_all, function(x){x[,"t"]}))))

# Do the optimisation
print('Starting the Profile Analysis')

# if not finished continue
if (!file.exists(paste0(out_dir_path,'/sa_bfgs_merged', gene_name, '_Finished.Rdata',sep=""))){
  # load_best
  load(paste0(scratch_dir, '/',folder, '/Optimisation/BFGS/',condition, '/param_bfgs_1000ri_',replicate, '_', gene_name, '_Finished.Rdata'))
  
  # for (replicate in c('rep1', 'b3', 'all', 'rep3')){
    if (!file.exists(paste0(out_dir_path,'/profile_dMod_', replicate, '_', gene_name, '_Finished.Rdata',sep=""))){
      sink(file = paste0(out_dir_path, '/Run_log/', gene_name, '_', replicate, '.txt'), append = TRUE, split = FALSE)
           
      if(replicate == "all"){
        b <- c("rep1", "rep2","rep3")
      }else{
        b <- replicate
      }
      
      if (condition == "naive") {
        data <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_Naive_exons_rpkm_5kb,   g, replicate=b, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      npRNA   = as.matrix(create_dat_g_rpkm(npRNA_Naive_exons_rpkm_5kb,   g, replicate=b, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_Naive_exons_rpkm_5kb, g, replicate=b, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
        data_all <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_Naive_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2","rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      npRNA   = as.matrix(create_dat_g_rpkm(npRNA_Naive_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2","rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_Naive_exons_rpkm_5kb, g, replicate=c("rep1", "rep2","rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
      }else if (condition == "lpa"){
        data <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_LPA_exons_rpkm_5kb,   g, replicate=b, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      npRNA   = as.matrix(create_dat_g_rpkm(npRNA_LPA_exons_rpkm_5kb,   g, replicate=b, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_LPA_exons_rpkm_5kb, g, replicate=b, time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
        data_all <- list( caRNA   = as.matrix(create_dat_g_rpkm(caRNA_LPA_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2","rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      npRNA   = as.matrix(create_dat_g_rpkm(npRNA_LPA_exons_rpkm_5kb,   g, replicate=c("rep1", "rep2","rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                      cytoRNA = as.matrix(create_dat_g_rpkm(cytoRNA_LPA_exons_rpkm_5kb, g, replicate=c("rep1", "rep2","rep3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
      }
      
      lib_mat <- as.matrix(lib[grep(condition, row.names(lib), ignore.case = T),])
      row.names(lib_mat) <- gsub(pattern = condition, replacement = '', x = row.names(lib_mat), ignore.case = T)
      row.names(lib_mat) <- gsub(pattern = '..', replacement = '.', x = row.names(lib_mat), fixed = T)
      
      tmp_file <- paste0(out_dir_path,'/sa_bfgs_', replicate, '_', gene_name, '.Rdata',sep="")
      bfgs <- bfgs[[replicate]]
      par_best <- bfgs$parameter[which.min( bfgs$value),]
    
      if(replicate == 'rep3'){
        no_np <- TRUE
      }else{
        no_np <- FALSE
      }
      
      if( replicate %in% c('rep1','all') & condition == "naive"){
        no_ca0 <- TRUE
      }else{ 
        no_ca0 <- FALSE 
      }
      
      # convert add names and combine param to the best
      if ( no_np ){
        par_names_best <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'k2'", "k2", "kdeg")
      }else{
        par_names_best <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'", "k2", "k2'", "kdeg")
      }
      
      if( !no_ca0 ){
        par_names_best <- par_names_best[par_names_best != "rep1_ca0"]
      }
      
      par_best_old <- c(par_best[!(par_names_best %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg")) ], parameters_convert_cp(par_best[(par_names_best %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg")) ], log_in = T, log_out = T, direction = -1, no_np = no_np))
      names(par_best_old) <- c(par_names_best[!(par_names_best %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))],par_names_best[(par_names_best %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))])
      
      ## add combined parameters
      if ( no_np ){
        par_best_old_transformed <- c(par_best_old, "k1'k2'/k2" = unname(par_best_old["k1'k2'"] - par_best_old["k2"])) 
        par_best_old_transformed <- par_best_old_transformed[-which(names(par_best_old_transformed) == "k1'k2'")]
      }else{
        par_best_old_transformed <- c(par_best_old,
                          "k1'/k2" = unname(par_best_old["k1'"] - par_best_old["k2"]),
                          "k1'k2'/k2" = unname(par_best_old["k1'"] + par_best_old["k2'"] - par_best_old["k2"]))
        par_best_old_transformed <- par_best_old_transformed[-which(names(par_best_old_transformed) %in% c("k1'","k2'"))]
      }

      
      # For dMod
      # obj_vgh <- function(pars, fixed = NULL, scale=2, no_np = F, no_ca0 = F ){
      #   if ( no_np ){
      #     par_names <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'k2'", "k2", "kdeg", "k1'k2'/k2",)
      #   }else{
      #     par_names <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'", "k2", "k2'", "kdeg", "k1'/k2","k1'k2'/k2")
      #   }
      #   
      #   if ( !no_ca0 ){
      #     par_names <- par_names[-which(par_names == "rep1_ca0")]
      #   }
      #   
      #   if ( !is.null(fixed) ){
      #     pars <- c(pars,fixed)
      #     fixed_index <- which( par_names == names(fixed) )
      #     
      #     if (names(fixed) == "k1'/k2"){
      #       tmp <- pars
      #       tmp["k1'"] <- tmp[fixed_index] + tmp["k2"]
      #     }else if( names(fixed) == "k1'k2'/k2" ){
      #       tmp <- pars
      #       if(!no_np){
      #         tmp["k1'"] <- tmp[fixed_index] + tmp["k2"] - tmp["k2'"]
      #       }else{
      #         tmp["k1'k2'"] <- tmp[fixed_index] + tmp["k2"]
      #       }
      #     }else{
      #       tmp <- pars
      #       if(!no_np){
      #         tmp["k1'/k2"] <- tmp["k1'"] - tmp["k2"]
      #         tmp["k1'k2/k2"] <- tmp["k1'"] + tmp["k2'"] - tmp["k2"]
      #       }else{
      #         tmp["k1'k2/k2"] <- tmp["k1'k2'"] - tmp["k2"]
      #       }
      #     }
      #   }else{
      #     tmp <- pars
      #     if(!no_np){
      #       tmp["k1'/k2"] <- tmp["k1'"] - tmp["k2"]
      #       tmp["k1'k2'/k2"] <- tmp["k1'"] + tmp["k2'"] - tmp["k2"]
      #     }else{
      #       tmp["k1'k2'/k2"] <- tmp["k1'k2'"] - tmp["k2"]
      #     }
      #   }
      #   pars <- tmp
      # 
      #   # reorder
      #   pars_ordered <- pars[par_names]
      #   
      #   # transform and reorder
      #   pars_new <- c(pars_ordered[ !(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg")) ], parameters_convert_cp(pars_ordered[(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))], log_in = T, log_out = T, direction = 1, no_np = no_np))
      #   names(pars_new) <- c(par_names[!(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))], par_names[(par_names %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))])
      #   pars_new_ordered <- pars_new[par_names]
      #   
      #   val <- unname(TSE_reps_cp(par = pars_new_ordered[-which(par_names %in% c("k1'/k2", "k1'k2'/k2"))], fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length))
      #   gradient_val <- TSE_reps_cp_grad(x = pars_new_ordered[-which(par_names %in% c("k1'/k2", "k1'k2'/k2"))], fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length)
      #   # add gradient for additional parnames
      #   for (parname_add in par_names[which(par_names %in% c("k1'/k2", "k1'k2'/k2"))]){
      #     if (parname_add == "k1'/k2"){
      #       gradient_val <- c(gradient_val, gradient_val[par_names == "k1'"] - gradient_val[par_names == "k2"])
      #     }else if(parname_add == "k1'k2'/k2"){
      #       if(no_np){
      #         gradient_val <- c(gradient_val, gradient_val[par_names == "k1'k2'"] - gradient_val[par_names == "k2"])
      #       }else{
      #         gradient_val <- c(gradient_val, gradient_val[par_names == "k1'"] + gradient_val[par_names == "k2'"] - gradient_val[par_names == "k2"])
      #       }
      #     }
      #   }
      #   
      #   hessian_val <- TSE_reps_cp_hessian(x = pars_new_ordered[-which(par_names %in% c("k1'/k2", "k1'k2'/k2"))],fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length)
      #   # add hessian for additional parnames
      #   for (parname_add in par_names[which(par_names %in% c("k1'/k2", "k1'k2'/k2"))]){
      #     if (parname_add == "k1'/k2"){
      #       vec <- rep(0,ncol(hessian_val))
      #       vec[which(par_names == "k1'")] <- 1
      #       vec[which(par_names == "k2")] <- -1
      #     }else if(parname_add == "k1'k2'/k2"){
      #       vec <- rep(0,ncol(hessian_val))
      #       if(no_np){
      #         vec[which(par_names == "k1'k2'")] <- 1
      #         vec[which(par_names == "k2")] <- -1
      #       }else{
      #         vec[which(par_names == "k1'")] <- 1
      #         vec[which(par_names == "k2'")] <- 1
      #         vec[which(par_names == "k2")] <- -1
      #       }
      #     }
      #     hessian_val_add <- hessian_val %*% vec
      #     hessian_val <- cbind(hessian_val, hessian_val_add)
      #     hessian_val <- rbind(hessian_val, c(hessian_val_add, t(hessian_val_add) %*% vec))
      #   }
      #   
      #   # hessian_val <- matrix(0, ncol=length(pars_new_ordered), nrow=length(pars_new_ordered))
      #   
      #   if ( !is.null(fixed) ){
      #     # 2* because not included in the chisq
      #     res <- objlist(
      #       value = scale * val,
      #       gradient = scale * gradient_val[-fixed_index],
      #       hessian = scale * hessian_val[-fixed_index,][,-fixed_index]
      #     )
      #   }else{
      #     # print(val)
      #     res <- objlist(
      #       value =  scale * val, 
      #       gradient = scale * gradient_val,
      #       hessian = scale * hessian_val
      #     )
      #   }
      #   
      #   #attr(res, controls$attr.name) <- res$value
      #   attr(res, "env") <- NULL
      #   return(res)
      # }
      # 
      # class(obj_vgh) <- c("objfn", "fn")
      # attr(obj_vgh, "conditions") <- NULL
      # attr(obj_vgh, "parameters") <- names(par_best_old)
      
      
      obj_vgh_main <- function(pars, fixed = NULL, scale=2, no_np = F, no_ca0 = F ){
        if ( no_np ){
          par_names <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'k2'", "k2", "kdeg")
        }else{
          par_names <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'", "k2", "k2'", "kdeg")
        }
        
        if ( !no_ca0 ){
          par_names <- par_names[-which(par_names == "rep1_ca0")]
        }
        
        if ( !is.null(fixed) ){
          pars <- c(pars,fixed)
          fixed_index <- which( par_names == names(fixed) )
        }
        
        # reorder
        pars_ordered <- pars[par_names]
        
        val <- unname(TSE_reps_oldpar(x = pars_ordered, fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length))
        gradient_val <- TSE_reps_cp_grad(x = pars_ordered, fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length)
        hessian_val <- TSE_reps_cp_hessian(x = pars_ordered,fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length)
        
        if ( val == 1*10^16){
          val <- NA
        }
        
        if ( !is.null(fixed) ){
          # 2* because not included in the chisq
          res <- objlist(
            value = scale * val,
            gradient = scale * gradient_val[-fixed_index],
            hessian = scale * hessian_val[-fixed_index,][,-fixed_index]
          )
        }else{
          # print(val)
          res <- objlist(
            value =  scale * val, 
            gradient = scale * gradient_val,
            hessian = scale * hessian_val
          )
        }
        
        #attr(res, controls$attr.name) <- res$value
        attr(res, "env") <- NULL
        return(res)
      }
      
      class(obj_vgh_main) <- c("objfn", "fn")
      attr(obj_vgh_main, "conditions") <- NULL
      attr(obj_vgh_main, "parameters") <- names(par_best_old)

      obj_vgh_transformed <- function(pars, fixed = NULL, scale=2, no_np = F, no_ca0 = F ){
        if ( no_np ){
          par_names <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k2", "kdeg", "k1'k2'/k2")
          par_names_old <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'k2'", "k2", "kdeg")
        }else{
          par_names <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k2", "kdeg", "k1'/k2","k1'k2'/k2")
          par_names_old <- c("span", "sigma_t", "sigma_b", "rep1_ca0", "k1'","k2", "k2'", "kdeg")
        }
        
        if ( !no_ca0 ){
          par_names <- par_names[-which(par_names == "rep1_ca0")]
          par_names_old <- par_names_old[-which(par_names_old == "rep1_ca0")]
        }
        
        if ( !is.null(fixed) ){
          pars <- c(pars,fixed)
          fixed_index <- which( par_names == names(fixed) )
        }  
        
        if(!no_np){
          pars_old <- c(pars, "k1'" = unname(pars["k1'/k2"] + pars["k2"]), "k2'" = unname(pars["k1'k2'/k2"] - pars["k1'/k2"]))
          pars_old <- pars_old[-which(names(pars_old) %in% c("k1'/k2", "k1'k2'/k2"))]
        }else{
          pars_old <- c(pars, "k1'k2'" = unname(pars["k1'k2'/k2"] + pars["k2"]))
          pars_old <- pars_old[-which(names(pars_old) %in% c("k1'k2'/k2"))]
        }
        

        # reorder
        pars_old_ordered <- pars_old[par_names_old]
        
        val <- unname(TSE_reps_oldpar(x = pars_old_ordered, fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length))
        
        if ( val == 1*10^16){
          val <- NA
        }
        
        gradient_val <- TSE_reps_cp_grad(x = pars_old_ordered, fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length)
        
        # add gradient for additional parnames
        # for (parname_add in par_names[which(par_names %in% c("k1'/k2", "k1'k2'/k2"))]){
        #   if (parname_add == "k1'/k2"){
        #     gradient_val <- c(gradient_val, gradient_val[par_names_old == "k1'"] - gradient_val[par_names_old == "k2"])
        #   }else if(parname_add == "k1'k2'/k2"){
        #     if(no_np){
        #       gradient_val <- c(gradient_val, gradient_val[par_names_old == "k1'k2'"] - gradient_val[par_names_old == "k2"])
        #     }else{
        #       gradient_val <- c(gradient_val, gradient_val[par_names_old == "k1'"] + gradient_val[par_names_old == "k2'"] - gradient_val[par_names_old == "k2"])
        #     }
        #   }
        # }
        # # remove 
        # gradient_val <- gradient_val[-which(par_names_old %in% c("k1'k2'", "k1'", "k2'"))]
        # 
        
        if (!no_np){
          trans_mat <- matrix(0, nrow=length(par_names), ncol=length(par_names_old))
          trans_mat[1:4,1:4] <- diag(1,nrow = 4) # span, sigma_t, sigma_b, rep1_ca0 unchanged
          trans_mat[, par_names_old == 'k2'] <- as.numeric(par_names == 'k2')
          trans_mat[, par_names_old == 'kdeg'] <- as.numeric(par_names == 'kdeg')
          trans_mat[, par_names_old == "k1'"] <- as.numeric(par_names == "k1'/k2") + as.numeric(par_names == "k2")
          trans_mat[, par_names_old == "k2'"] <- as.numeric(par_names == "k1'k2'/k2")-as.numeric(par_names == "k1'/k2")
        }else{
          trans_mat <- matrix(0, nrow=length(par_names), ncol=length(par_names_old))
          trans_mat[1:4,1:4] <- diag(1,nrow = 4) # span, sigma_t, sigma_b, rep1_ca0 unchanged
          trans_mat[, par_names_old == 'k2'] <- as.numeric(par_names == 'k2')
          trans_mat[, par_names_old == 'kdeg'] <- as.numeric(par_names == 'kdeg')
          trans_mat[, par_names_old == "k1'k2'"] <- as.numeric(par_names == "k1'k2'/k2") + as.numeric(par_names == "k2")
        }
        
        gradient_val <- trans_mat %*% gradient_val
        
        
        hessian_val <- TSE_reps_cp_hessian(x = pars_old_ordered,fn_fit=fn_all_ode_cp, data = data_rpkm_all, time_data=time_data, convert_param=T, log_in=T, u=input_fun_cp, replicate=replicate, results = "sum", library_size = lib_mat, length = length)
        
        hessian_val <- trans_mat %*% hessian_val %*% t(trans_mat)
        
        # add hessian for additional parnames
        # for (parname_add in par_names[which(par_names %in% c("k1'/k2", "k1'k2'/k2"))]){
        #   if (parname_add == "k1'/k2"){
        #     vec <- rep(0,ncol(hessian_val))
        #     vec[which(par_names_old == "k1'")] <- 1
        #     vec[which(par_names_old == "k2")] <- -1
        #   }else if(parname_add == "k1'k2'/k2"){
        #     vec <- rep(0,ncol(hessian_val))
        #     if(no_np){
        #       vec[which(par_names_old == "k1'k2'")] <- 1
        #       vec[which(par_names_old == "k2")] <- -1
        #     }else{
        #       vec[which(par_names_old == "k1'")] <- 1
        #       vec[which(par_names_old == "k2'")] <- 1
        #       vec[which(par_names_old == "k2")] <- -1
        #     }
        #   }
        #   hessian_val_add <- hessian_val %*% vec
        #   hessian_val <- cbind(hessian_val, hessian_val_add)
        #   hessian_val <- rbind(hessian_val, c(hessian_val_add, t(hessian_val_add) %*% vec))
        # }
        # 
        # # hessian remove
        # hessian_val <- hessian_val[,-which(par_names_old %in% c("k1'k2'", "k1'", "k2'"))][-which(par_names_old %in% c("k1'k2'", "k1'", "k2'")),]
        # 
        # hessian_val <- matrix(0, ncol=length(pars_new_ordered), nrow=length(pars_new_ordered))
        
        if ( !is.null(fixed) ){
          # 2* because not included in the chisq
          res <- objlist(
            value = scale * val,
            gradient = scale * gradient_val[-fixed_index],
            hessian = scale * hessian_val[-fixed_index,][,-fixed_index]
          )
        }else{
          # print(val)
          res <- objlist(
            value =  scale * val, 
            gradient = scale * gradient_val,
            hessian = scale * hessian_val
          )
        }
        
        #attr(res, controls$attr.name) <- res$value
        attr(res, "env") <- NULL
        return(res)
      }
      
      class(obj_vgh_transformed) <- c("objfn", "fn")
      attr(obj_vgh_transformed, "conditions") <- NULL
      attr(obj_vgh_transformed, "parameters") <- names(par_best_old_transformed)
      
      # fit_best <- trust(objfun = obj_vgh, parinit = par_best_old, rinit = 0.1, rmax= 10)
      # test.approx <- profile(obj = obj_vgh, pars = par_best_old, whichPar = c(5), limits = c(-3,3), method = "integrate", cores = 4, verbose = T)
      
      if(file.exists(tmp_file)){
        load(tmp_file)
      }else{
        dat_g <- c()
      }
      
      
      param_to_do <- which(names(par_best_old) %in% c("k1'","k1'k2'", "k2", "k2'", "kdeg"))
      i_to_do <- sample(param_to_do, length(param_to_do),replace=F)
      for(i in i_to_do){
        print(paste0('Doing Profile for: ', names(par_best_old)[i]))
        
        if(names(par_best_old)[i] %in% names(dat_g)){
          next
        }else{
          cpt1 <- proc.time()
          # test.exact <- profile(obj = obj_vgh, pars = par_best_old, whichPar = i, limits = c(-3,3), method = "optimize", cores = 1, verbose = T, scale = 2, no_np = no_np, no_ca0 = no_ca0)
          # test.exact <- profile_dMod_with_save(obj = obj_vgh, pars = par_best_old, whichPar = i, limits = c(-3,3), method = "optimize", cores = 1, verbose = T, scale = 2, no_np = no_np, no_ca0 = no_ca0, file_name = "/tmp/dMod_test_save.Rdata")
          test.approx <- profile_dMod_with_save(obj = obj_vgh_main, pars = par_best_old, whichPar = i, limits = c(-10,10), method = "integrate", cores = 1, verbose = T, scale = 2, no_np = no_np, no_ca0 = no_ca0, file_name = paste0(out_dir_path,'/temp_par_',names(par_best_old)[i], '_', replicate, '_', gene_name, '.Rdata',sep=""), stepControl=list(limit=100000), algoControl = list(reoptimize = T))
          
          # plotProfile(test.exact)
          # plotPaths(test.exact, sort = TRUE)
          cf <- confint(test.approx, val.column = "value")
          cpt2 <- proc.time()
          
          dat_g <- c(dat_g, list(list(Profile = test.approx, CI=cf, computing_time = cpt2 - cpt1)))
          names(dat_g)[length(dat_g)] <- names(par_best_old)[i]
          save(dat_g, file=tmp_file)
          
          if(file.exists(paste0(out_dir_path,'/temp_par_',names(par_best_old)[i], '_', replicate, '_', gene_name, '.Rdata',sep=""))){
            file.remove(paste0(out_dir_path,'/temp_par_',names(par_best_old)[i], '_', replicate, '_', gene_name, '.Rdata',sep=""))
          }
          # quit('no') # start the next param in a clean 24h session
        }  
      }
      
      # additional parameters
      param_to_do <- which(names(par_best_old_transformed) %in% c( "k1'/k2", "k1'k2'/k2"))
      i_to_do <- sample(param_to_do, length(param_to_do),replace=F)
      for(i in i_to_do){
        print(paste0('Doing Profile for: ', names(par_best_old_transformed)[i]))
        
        if(names(par_best_old_transformed)[i] %in% names(dat_g)){
          next
        }else{
          cpt1 <- proc.time()
          # test.exact <- profile(obj = obj_vgh, pars = par_best_old, whichPar = i, limits = c(-3,3), method = "optimize", cores = 1, verbose = T, scale = 2, no_np = no_np, no_ca0 = no_ca0)
          # test.exact <- profile_dMod_with_save(obj = obj_vgh, pars = par_best_old, whichPar = i, limits = c(-3,3), method = "optimize", cores = 1, verbose = T, scale = 2, no_np = no_np, no_ca0 = no_ca0, file_name = "/tmp/dMod_test_save.Rdata")
          test.approx <- profile_dMod_with_save(obj = obj_vgh_transformed, pars = par_best_old_transformed, whichPar = i, limits = c(-3,3), method = "integrate", cores = 1, verbose = T, scale = 2, no_np = no_np, no_ca0 = no_ca0, file_name = paste0(out_dir_path,'/temp_par_',gsub('/','over',names(par_best_old_transformed)[i]), '_', replicate, '_', gene_name, '.Rdata',sep=""), stepControl=list(limit=100000))
          
          # plotProfile(test.exact)
          # plotPaths(test.exact, sort = TRUE)
          cf <- confint(test.approx, val.column = "value")
          cpt2 <- proc.time()
          
          dat_g <- c(dat_g, list(list(Profile = test.approx, CI=cf, computing_time = cpt2 - cpt1)))
          names(dat_g)[length(dat_g)] <- names(par_best_old_transformed)[i]
          save(dat_g, file=tmp_file)
          
          if(file.exists(paste0(out_dir_path,'/temp_par_',gsub('/','over',names(par_best_old_transformed)[i]), '_', replicate, '_', gene_name, '.Rdata',sep=""))){
            file.remove(paste0(out_dir_path,'/temp_par_',gsub('/','over',names(par_best_old_transformed)[i]), '_', replicate, '_', gene_name, '.Rdata',sep=""))
          }
          # quit('no') # start the next param in a clean 24h session
        }  
      }

      sink()
      save(dat_g, file = paste0(out_dir_path,'/profile_dMod_',replicate, '_', gene_name, '_Finished.Rdata',sep=""))
      file.remove(tmp_file)
    }
  # }

  # if all finished merge
  dat_profile <- list()
  all_finished <- TRUE
  for (replicate in c('rep1', 'rep2', 'all', 'rep3')){
    if(!file.exists(paste0(out_dir_path,'/profile_dMod_', replicate, '_', gene_name, '_Finished.Rdata',sep=""))){
      all_finished <- FALSE
      break
    }
  }
  if(all_finished){
    dat_sa <- c()
    for (replicate in c('rep1', 'rep2', 'all', 'rep3')){
      load(file = paste0(out_dir_path,'/profile_dMod_', replicate, '_', gene_name, '_Finished.Rdata',sep=""))
      dat_sa <- c(dat_sa, list(dat_g))
      names(dat_sa)[length(dat_sa)] <- replicate
    } 
    save(dat_sa, file = paste0(out_dir_path,'/profile_dMod_merged_', gene_name, '_Finished.Rdata',sep=""))
  }
}

quit('no')
