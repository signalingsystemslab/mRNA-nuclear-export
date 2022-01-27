## Get input arguments
print('Starting the Run')
args = commandArgs(TRUE)

print('Given Arguments:')
print(args)

i <- as.numeric(args[1])
iter <- as.numeric(args[2])
condition <- as.character(args[3])
folder <- as.character(args[4])
batch <- as.character(args[5])

print(paste0('i: ', i))
print(paste0('iter: ', iter))
print(paste0('condition: ', condition))
print(paste0('folder: ', folder))
print(paste0('batch: ', batch))

## Load Packages
print('Loading packages')
for (pkg in c('compiler', 'deSolve')){
  if (!(pkg %in% installed.packages())){
    install.packages(pkg, repos = 'https://cran.cnr.berkeley.edu/')
  }
}

library(compiler)
library(deSolve)
source('functions_weighting_for_smoothing_include_negbinom_error_in_model_v13.R') # change library to matrix to be faster and use mu instead of prob for the neg binomial (less issues for small numbers)

################
## Load Rdata ##
################
print('Creating needed folders')

scratch_dir <- Sys.getenv("SCRATCH")
print(scratch_dir)
out_dir <- paste0(scratch_dir,'/',folder)
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir, '/Optimisation')
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir, '/BFGS')
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir,'/', condition)
dir.create(file.path(out_dir), showWarnings = FALSE)
out_dir <- paste0(out_dir,'/Run_log')
dir.create(file.path(out_dir), showWarnings = FALSE)

out_dir_path <- paste0(scratch_dir,'/',folder,'/Optimisation/BFGS/', condition)

####################################
## Get data for all compartiments ##
####################################
print('Loading data')

## Exons last 5kb ----
if (condition == 'naive'){
  load('caRNA_Naive_exons_rpkm_5kb.Rdata')
  load('npRNA_Naive_exons_rpkm_5kb.Rdata')
  load('cytoRNA_Naive_exons_rpkm_5kb.Rdata')
  g <- rownames(caRNA_Naive_exons_rpkm_5kb)[i]
}else if (condition == 'lpa'){
  load('caRNA_LPA_exons_rpkm_5kb.Rdata')
  load('npRNA_LPA_exons_rpkm_5kb.Rdata')
  load('cytoRNA_LPA_exons_rpkm_5kb.Rdata')
  g <- rownames(caRNA_LPA_exons_rpkm_5kb)[i]
}

load('gene_infos_with_length.Rdata')
gene_name <- gene_infos$external_gene_name[gene_infos$ensembl_gene_id == substr(g,1,18)]

load('lib.size.Rdata')
lib <- lib_exons[,c("lib.size", "norm.factors")]
row.names(lib) <- lib_exons$samples


print(g)
print(gene_name)

if (condition == "naive"){
  data_rpkm_all <- list( caRNA   = as.matrix(create_dat_g_rpkm_cp(caRNA_Naive_exons_rpkm_5kb,   g, batch=c("b1", "b2", "b3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         npRNA   = as.matrix(create_dat_g_rpkm_cp(npRNA_Naive_exons_rpkm_5kb,   g, batch=c("b1", "b2", "b3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         cytoRNA = as.matrix(create_dat_g_rpkm_cp(cytoRNA_Naive_exons_rpkm_5kb, g, batch=c("b1", "b2", "b3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
  length <- gene_infos$Length_5kb_naive[gene_infos$ensembl_gene_id == substr(g,1,18)]
}else{
  data_rpkm_all <- list( caRNA   = as.matrix(create_dat_g_rpkm_cp(caRNA_LPA_exons_rpkm_5kb,   g, batch=c("b1", "b2", "b3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         npRNA   = as.matrix(create_dat_g_rpkm_cp(npRNA_LPA_exons_rpkm_5kb,   g, batch=c("b1", "b2", "b3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))), 
                         cytoRNA = as.matrix(create_dat_g_rpkm_cp(cytoRNA_LPA_exons_rpkm_5kb, g, batch=c("b1", "b2", "b3"), time = c("000","010","015","020","025","030","035","040","050","060","075","090","120"))))
  length <- gene_infos$Length_5kb_lpa[gene_infos$ensembl_gene_id == substr(g,1,18)]
}

time_data <- sort(unique(unlist(lapply(data_rpkm_all, function(x){x[,"t"]}))))

# Do the optimisation
print('Starting the optimisation')

# if not finished continue
if (!file.exists(paste0(out_dir_path,'/param_bfgs_', iter, 'ri_', gene_name, '_Finished.Rdata',sep=""))){
    if (!file.exists(paste0(out_dir_path,'/param_bfgs_', iter, 'ri_',batch, '_', gene_name, '_Finished.Rdata',sep=""))){
        sink(file = paste0(out_dir_path, '/Run_log/', gene_name, '_', batch, '.txt'), append = TRUE, split = FALSE)
           
        tmp_file <- paste0(out_dir_path,'/param_bfgs_', iter, 'ri_', batch, '_', gene_name, '.Rdata',sep="")
            
        lib_mat <- as.matrix(lib[grep(condition, row.names(lib), ignore.case = T),])
        row.names(lib_mat) <- gsub(pattern = condition, replacement = '', x = row.names(lib_mat), ignore.case = T)
        row.names(lib_mat) <- gsub(pattern = '..', replacement = '.', x = row.names(lib_mat), fixed = T)
      
        tmp_bfgs <- optim_all_bs_cp(n_init = iter, obj_fn = TSE_bs_cp, data_rpkm = data_rpkm_all, time_data=time_data, batch = batch, save_file=tmp_file, length = length, libsize = lib_mat)

        cpt_total <- sum(tmp_bfgs$time)

        bfgs <- list(all = tmp_bfgs, computing_time = list(cpt_total))
        names(bfgs)[1] <- batch
            
        save(bfgs, file = paste0(out_dir_path,'/param_bfgs_', iter, 'ri_', batch, '_', gene_name, '_Finished.Rdata',sep=""))
        file.remove(tmp_file)
        sink()
    }
  
  all_finished <- TRUE
  # if all finished merge
  bfgs_merged <- list()
  cpt <- c()
  
  for (batch in c('all', 'b1', 'b2', 'b3')){
    if(!file.exists(paste0(out_dir_path,'/param_bfgs_', iter, 'ri_',batch, '_', gene_name, '_Finished.Rdata',sep=""))){
        all_finished <- FALSE
    }
  }
  
  if(all_finished){
    for (batch in c('all', 'b1', 'b2', 'b3')){
        load(file = paste0(out_dir_path,'/param_bfgs_', iter, 'ri_',batch, '_', gene_name, '_Finished.Rdata',sep=""))
        bfgs_merged <- c(bfgs_merged, list(bfgs[[1]]))
        names(bfgs_merged)[length(bfgs_merged)] <- batch
        cpt <- c(cpt, as.vector(bfgs[[2]]))
        names(cpt)[length(cpt)] <- batch
    } 
    bfgs_merged <- c(bfgs_merged, computing_time = list(cpt))
    save(bfgs_merged, file = paste0(out_dir_path,'/param_bfgs_', iter, 'ri_merged_', gene_name, '_Finished.Rdata',sep=""))
  }
}

quit('no')
