
## Convert old to new parametrisation
parameters_convert <- function(par, log_in, log_out=log_in, direction=-1, parameter="row", no_np = FALSE){
  if(is.numeric(par) & is.null(dim(par))){
    par <- matrix(par, ncol=1)
  }
  
  if(parameter=="col"){
    par <- t(par)
  }
  if (log_in) {
    tmp <-  par
    par <- 10^tmp
  }
  
  par_old <- par
  par_new <- par
  if(no_np){
    if(direction == -1){ # convert from new to old
      par_old[1,] <- par_new[1,] # k1'* k2' = k_ca2np * (mass_ca / mass_np) * k_np2cyt * mass_np/mass_cyt
      par_old[2,] <- (par_new[2,] + sqrt(par_new[2,]^2-4* par_new[3,]))/2 # k2 or kdeg_cyt
      par_old[3,] <- (par_new[2,] - sqrt(par_new[2,]^2-4* par_new[3,]))/2 # k2 or kdeg_cyt
      res <- par_old
    }else if(direction == 1) { # convert from old to new
      par_new[3,] <- par_old[2,] * par_old[3,] # par_new[3] = knp2cyt * kdeg_cyto
      par_new[2,] <- par_old[2,] + par_old[3,] # par_new[2] = knp2cyt + kdeg_cyto 
      par_new[1,] <- par_old[1,]  # par_new[1] = k_ca2np * (mass_ca / mass_np) * k_np2cyt * mass_np/mass_cyt
      res <- par_new
    }
  }else{
    if(direction == -1){ # convert from new to old
      par_old[1,] <- par_new[1,] # k_ca2np * mass_ca / mass_np
      par_old[2,] <- par_old[1,] / par_new[2,] # (k_np2cyt + k_np_deg)
      par_old[3,] <- par_old[2,] * par_new[3,] # k_np2cyt * mass_np/mass_cyt
      par_old[4,] <- par_old[3,] / par_new[4,] # kdeg_cyt
      res <- par_old
    }else if(direction == 1) { # convert from old to new
      par_new[4,] <- par_old[3,] / par_old[4,] # par_new[11] = cyt_0/np_0 = k_np2cyto * mass_np / mass_cyt / k_deg_cyto 
      par_new[3,] <- par_old[3,] / par_old[2,] # par_new[10] = k_np2cyto / (k_np2cyto + k_np_deg) * (mass_np / mass_cyt)
      par_new[2,] <- par_old[1,] / par_old[2,] # par_new[9] = np_0/ca_0 = k_ca2np * mass_ca/mass_np / (k_np2cyto + k_np_deg)
      par_new[1,] <- par_old[1,]  # par_new[8] = k_ca2np * (mass_ca / mass_np)
      res <- par_new
    }
  }
  if(log_out){
    tmp <-  res
    res <- log10(tmp)
  }
  
  if(parameter=="col"){
    res <- t(res)
  }
  
  return(res)
}
parameters_convert_cp <- cmpfun(parameters_convert)

## Input function
input_fun <- function(t, input, start = 0, end = 120){
  t_index <- floor((t-start)/(1/60))+1
  t_index[ t > 360 ] <- ((360-start)/(1/60) + 1 + floor((log(t[t > 360]-360,60)+1)/0.05) + 1)
  
  # t_index[t<start | t>end] <- NA
  t_index[ t < start ] <- NA
  t_index[ t > end ] <- NA
  
  if (length(t) == 1){
    it <- t(input[t_index,])
    itp1 <- t(input[t_index+1,])
  }else{
    it <- input[t_index,]
    itp1 <- input[t_index+1,]
  }
  alpha <- (t-it[,1])/(itp1[,1]-it[,1])
  
  # return(sapply(1:length(t_index),function(i){if(t[i] <= start) { input[1,-1] }else if(t[i] <= end){colSums(input[t_index[i] +0:1,] * c(1-alpha[i],alpha[i]))[-1]}else{rep(NA, ncol(input)-1)}}))
  res <- matrix(it[,-1]* (1-alpha) + itp1[,-1] * alpha, ncol=length(t), byrow=T)
  res [,t <= start] <- input[1,-1]
  return(res)
}
input_fun_cp <- cmpfun(input_fun)

## Steady State
fn_all_ss <-function(par, caRNA_0, no_np=FALSE){
   # npRNA_0 <- par[1] * caRNA_0 / par[2] # npRNA(0) = k_ca2np * caRNA(0) / (k_np2cyto + k_np_deg) * mass_ca / mass_np
  # cytoRNA_0 <- par[3] * npRNA_0 / par[4] # cytoRNA(0) = k_np2cyto * npRNA(0) / k_cyto_deg * mass_np / mass_cyto
  
  if(no_np){
    res<- c(cytoRNA = par[1] / (par[2] * par[3]) * caRNA_0, dcytoRNA = 0) 
  }else{
    res<- c(npRNA = 1  , cytoRNA = par[3] / par[4]  ) * par[1] / par[2] * caRNA_0
  }
  return(res)
}
fn_all_ss_cp <- cmpfun(fn_all_ss)

## ODE systems
fn_all_derivative <- function(t, time_data, y, par, u, extrapolated_input, no_np=FALSE){ # faster
  if (no_np){
    dcytoRNA <- y[2]
    d2cytoRNA <- par[1] * u(t, extrapolated_input) - (par[2]+par[3])*y[2] - par[2]*par[3]*y[1]
    return(list(c(dcytoRNA, d2cytoRNA)))
  }else{
    dnpRNA <- par[1] * u(t, extrapolated_input)  - par[2] * y[1] # dnpRNA/dt = mass_ca / mass_np * k_ca2np * caRNA - (k_np2cyto + k_np_deg) * npRNA
    dcytoRNA <- par[3] * y[1] - par[4] * y[2] # dcytoRNA/dt = mass_np/mass_cyto * k_np2cyto * npRNA - k_cyto_deg * cytoRNA
    return(list(c(dnpRNA, dcytoRNA)))
  }
  # return(list(c(par[1] * u(t - time_data[1], extrapolated_input)  - par[2] * y[1],  # dnpRNA/dt = mass_ca / mass_np * k_ca2np * caRNA - (k_np2cyto + k_np_deg) * npRNA
  #               par[3] * y[1] - par[4] * y[2]))) # dcytoRNA/dt = mass_np/mass_cyto * k_np2cyto * npRNA - k_cyto_deg * cytoRNA
} 
fn_all_derivative_cp <- cmpfun(fn_all_derivative)

fn_all_ode <- function(t, par, convert_param = F, log_in = F, extrapolated_input, u, no_np=FALSE) {
  
  if (convert_param){
    par <- parameters_convert_cp(par, log_in = log_in, log_out = F, direction = -1, no_np = no_np)
  }
  options(warn=-1)
  on.exit(options(warn = 0))
  capture.output(
    res_ode <- ode(y = fn_all_ss_cp(par = par, u(0, extrapolated_input), no_np=no_np) , 
                   times = t, 
                   func = fn_all_derivative_cp, 
                   par = par, 
                   tcrit=120,  
                   u = u, 
                   extrapolated_input = extrapolated_input, time_data=t, no_np=no_np)
  )
  if(no_np){
   res_ode <- res_ode[,-which(colnames(res_ode)=="dcytoRNA")]
  }
  return(res_ode)
}
fn_all_ode_cp <- cmpfun(fn_all_ode)

# get derivative of cyto when no_np
# fn_all_ode_2 <- function(t, par, convert_param = F, log_in = F, extrapolated_input, u, no_np=FALSE) {
#  
#   if (convert_param){
#     par <- parameters_convert_cp(par, log_in = log_in, log_out = F, direction = -1, no_np = no_np)
#   }
#   options(warn=-1)
#   on.exit(options(warn = 0))
#   capture.output(
#     res_ode <- ode(y = fn_all_ss_cp(par = par, u(0, extrapolated_input), no_np = no_np),
#                    times = t, 
#                    func = fn_all_derivative_cp, 
#                    par = par, 
#                    tcrit=120,  
#                    u = u, 
#                    extrapolated_input = extrapolated_input, time_data=t, no_np=no_np)
#   )
#   return(res_ode)
# }
# fn_all_ode_2_cp <- cmpfun(fn_all_ode_2)


## Cost functions

TSE_bs <- function(par, fn_fit, data, time_data, batch, results="sum", library_size, length, ...){
  # par:
  # fn_fit:
  # data: data of all replicates
  # time_data: 
  # batch: 
  # results: 
  # library_size: library size for each sample, need to be subsetted for that condition 
  # length:
  # ...: 
  
  # check if the batch has np data or not
  no_np <- FALSE
  if(batch == 'b2'){
    no_np <- TRUE
  }
  
  # extract the smoothing parameter and variance from the parameter list
  span <- par[1]
  sigmat <- par[2]
  sigmab <- par[3]
  
  par <- par[-c(1:3)]
  
  # keep al data and extract the data corresponding to the wanted batch
  data_rpkm_all <- data
  if(batch != 'all'){
    b <- as.numeric(gsub('b','',batch)) 
    data <- lapply(data_rpkm_all, function(l){l[,c(1,1+b)]})
  }
  
  # create data.frame of caRNA for the weight for the smooting
  dat_g_exp <- data.frame(Time=rep(data_rpkm_all$caRNA[,1], ncol(data_rpkm_all$caRNA)-1), RPKM = as.vector(data_rpkm_all$caRNA[,-1]), batch = rep(c("b1", "b2", "b3"), each=nrow(data_rpkm_all$caRNA)))
 
  tryCatch({
    # times to extrapolate, every second
    extra_t_end <- time_data[length(time_data)]
    extra_t_seq <- seq(time_data[1],extra_t_end + 1, by=1/60)
    
    # extract data to extrapolate
    data_extrapolation <- data$caRNA
    
    # replace missing t=0 by parameters 
    if(any(is.na(data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])]))){
      data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])] <- par[1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,]))]
      missing_ca <- par[c(1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
      par <- par[-c(1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
    }
    
    # remove time point if all na
    to_remove <- apply(data_extrapolation, 1, function(x){all(is.na(x[-1]))})
    data_extrapolation <- data_extrapolation[!to_remove,]
    
    # extract library size #takes time because dataframe -> changed to mat
    # libsize <- as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''), as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))
    dat_g_exp_pattern <- paste0('chromatin', '.',sapply(nchar(as.character(dat_g_exp$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp$Time),'.',dat_g_exp$batch)
    # dat_g_exp_pattern_2 <- sapply(1:nrow(dat_g_exp),function(i){paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''), as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i])})
    # lib_indices <- sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = dat_g_exp_pattern[i], x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}})
    lib_indices <- match(dat_g_exp_pattern, row.names(library_size))
    libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[lib_indices]))
    
    # Calculate prob of p using a beta distribution: k ~ binom(p,N), p ~ beta(k+1, N-k+1) in log f is t*f(t)
    proba <- sapply(1:nrow(dat_g_exp), function(i){ 2^dat_g_exp$RPKM[i] * length / 10^3 / 10^6 * dbeta( x = 2^dat_g_exp$RPKM[i] * length / 10^3 / 10^6, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1)})

    # prob_rpkm_ub <- log2(sapply(1:nrow(dat_g_exp), function(i){10^6*10^3/length * qbeta( p = 0.25, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, lower.tail=F)}))
    # prob_rpkm_lb <- log2(sapply(1:nrow(dat_g_exp), function(i){10^6*10^3/length * qbeta( p = 0.25, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, lower.tail=T)}))
    
    res <- vector(mode = "list", length=ncol(data$caRNA)-1)
    # for each batch
    for (b in 2:(ncol(data$caRNA))){
      # extrapolate caRNA for specific batch 
      data_extra <- data.frame(x=as.vector(data_extrapolation[,1]), y=as.vector(data_extrapolation[,b]))
      data_extra <- data_extra[!is.na(data_extra$y),]

      # calculate weight as 1/IQR(p) saturating at 10
      if(ncol(data$caRNA) > 2){
        weight <- proba[dat_g_exp$batch == paste0('b',b-1)]
        # weight <- sapply(1:length(prob_rpkm_ub),function(i){1/max(0.1,prob_rpkm_ub[i]-prob_rpkm_lb[i])})[dat_g_exp$batch == paste0('b',b-1)] # don't go to low values
      }else{
        weight <- proba[dat_g_exp$batch == batch]
        # weight <- sapply(1:length(prob_rpkm_ub),function(i){1/max(0.1,prob_rpkm_ub[i]-prob_rpkm_lb[i])})[dat_g_exp$batch == batch] # don't go to low values
      }
      weight <- weight/max(weight, na.rm=T)
      # add max weigth to the missing t=0 values replaced by parameter (i.e. sure values)
      if(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))){
        weight[which(is.na(data$caRNA[data$caRNA[,1]==0,-1]))] <- 1
      }
      weight <- weight[!is.na(weight)]
      
      # order data to extrapolate and weight
      weight <- weight[order(data_extra$x)]
      data_extra <- data_extra[order(data_extra$x),]
      
      # extrapolate
      extrapolated_input <-  cbind(extra_t_seq, 2^predict(smooth.spline(x=data_extra, w=weight, spar = span), extra_t_seq)$y)
      
      # simulate data for every minutes to get derivative 
      fit <- fn_fit(t = seq(time_data[1], time_data[length(time_data)], by=1), par, extrapolated_input = extrapolated_input, no_np = no_np, ...)
      if (any(is.nan(fit) | is.na(fit) |  fit < 0 | nrow(fit)<length(seq(time_data[1], time_data[length(time_data)], by=1)))) {
        # print(p)
        # print(fit)
        return(10^16)
      }
      fit[,-1] <- log2(fit[,-1])

      # calculate the derivative (not using the model, because if quassi steady state derivative from model will be null)
      fit_deriv_ca <- data.frame(Time = extrapolated_input[,1], caRNA=c((extrapolated_input[2,2]-extrapolated_input[1,2])/(extrapolated_input[2,1]-extrapolated_input[1,1]), 
                                 0.5*((extrapolated_input[,2][-c(nrow(extrapolated_input)-1:0)]-extrapolated_input[,2][-c(1,nrow(extrapolated_input))])/(extrapolated_input[,1][-c(nrow(extrapolated_input)-1:0)]-extrapolated_input[,1][-c(1,nrow(extrapolated_input))])+(extrapolated_input[,2][-c(1:2)]-extrapolated_input[,2][-c(1,nrow(extrapolated_input))])/(extrapolated_input[,1][-c(1:2)]-extrapolated_input[,1][-c(1,nrow(extrapolated_input))])),
                                 (extrapolated_input[nrow(extrapolated_input),2]-extrapolated_input[nrow(extrapolated_input)-1,2])/(extrapolated_input[nrow(extrapolated_input),1]-extrapolated_input[nrow(extrapolated_input)-1,1])))
      
      fit_deriv <- data.frame(Time = fit_deriv_ca$Time, caRNA = fit_deriv_ca$caRNA/extrapolated_input[,2])[(seq(min(time_data),max(time_data)) - time_data[1])*60+1,]
       
      for (col in colnames(fit)){
        if (col == "npRNA"){
          fit_deriv_col <- log(2)*data.frame(npRNA=c((fit[2,col]-fit[1,col])/(fit[2,1]-fit[1,1]), 
                                             0.5*((fit[,col][-c(nrow(fit)-1:0)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(nrow(fit)-1:0)]-fit[,1][-c(1,nrow(fit))])+
                                             (fit[,col][-c(1:2)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(1:2)]-fit[,1][-c(1,nrow(fit))])),
                                    (fit[nrow(fit),col]-fit[nrow(fit)-1,col])/(fit[nrow(fit),1]-fit[nrow(fit)-1,1])))
          fit_deriv <- cbind(fit_deriv, npRNA = fit_deriv_col)
        }
        if (col == "cytoRNA"){
          fit_deriv_col <- log(2)*data.frame(cytoRNA=c((fit[2,col]-fit[1,col])/(fit[2,1]-fit[1,1]), 
                                             0.5*((fit[,col][-c(nrow(fit)-1:0)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(nrow(fit)-1:0)]-fit[,1][-c(1,nrow(fit))])+
                                                  (fit[,col][-c(1:2)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(1:2)]-fit[,1][-c(1,nrow(fit))])),
                                           (fit[nrow(fit),col]-fit[nrow(fit)-1,col])/(fit[nrow(fit),1]-fit[nrow(fit)-1,1])))
          fit_deriv <- cbind(fit_deriv, cytoRNA = fit_deriv_col)
        }
      }
      
      # fit[,-1][fit[,-1] <= -20] <- -20
      
      # extract time matching time data
      fit <- fit[match(time_data, fit[,1]),]  
      fit <- cbind(fit, caRNA=log2(extrapolated_input[(time_data - time_data[1])*60+1,2]))
      
      fit_deriv <- fit_deriv[match(time_data, fit_deriv[,1]),]  
 
      
      # compare data to simulations
      for (cpt in names(data)){
        if(!all(is.na(data[[cpt]][,b]))){
          if(cpt == "caRNA"){
            cpt_name <- "chromatin"
          }else if(cpt == "npRNA"){
            cpt_name <- "Nuc"
          }else if(cpt == "cytoRNA"){
            cpt_name <- "cytoplasmic"
          }
          # print(cbind(fit[match(data[[cpt]][,"t"], time_data), cpt], data[[cpt]]))
          
          # takes time because data.frame
          # libsize <- as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0(cpt_name,'.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$batch == paste0('b',b-1)]
          # libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0(cpt_name, '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$batch == paste0('b',b-1)]
          dat_g_exp_pattern <- paste0(cpt_name, '.',sapply(nchar(as.character(dat_g_exp$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp$Time),'.', dat_g_exp$batch)
          # lib_indices <- sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = dat_g_exp_pattern[i], x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}})
          lib_indices <- match(dat_g_exp_pattern, row.names(library_size))
          libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[lib_indices]))[dat_g_exp$batch == paste0('b',b-1)]
          
          res[[b-1]] <- c(res[[b-1]], TSE_neg_binom_cp(fit[match(data[[cpt]][,"t"], time_data), cpt], fit_deriv[match(data[[cpt]][,"t"], time_data), cpt],  data[[cpt]][,b], libsize, length, sigmab, sigmat))
        }else{
          res[[b-1]] <- c(res[[b-1]], NA)
        }
      }
    }
    
    if (results != "detailed"){
      res <- sum(unlist(lapply(res, sum, na.rm=T)))
      # add prior of sigmat and sigmab ~ regularisation (2 times higher because can only be positive)
      # v1
      # res <- res + log(2*10*sqrt(2*pi)) + sigmat^2/(2*10^2) + log(2*0.5*sqrt(2*pi)) + sigmab^2/(2*0.5^2) + log(sqrt(2*pi)*0.1) + (span - 0.45)^2/(2*0.2^2)
      
      # v2
      res <- res - log(2*1/(5*sqrt(2*pi))) + sigmat^2/(2*5^2) - log(2/(0.1*sqrt(2*pi))) + sigmab^2/(2*0.1^2) - log(1/(sqrt(2*pi)*0.05)) + (span - 0.45)^2/(2*0.05^2)
      
      # add regul on b1_ca0 based the other ca_0
      if(any(is.na(data$caRNA[data$caRNA[,1]==0,]))){
        ca_other_0 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==0,-1] 
        ca_other_10 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==10,-1] 
        
        delta <- mean(ca_other_0 - ca_other_10, na.rm=T)
        
        # libsize <- mean(as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$Time == 0], na.rm=T)
    
        # proba <- pbeta( q = 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])])* length / 10^3 / 10^6, shape1 = 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]) * length/10^3 * libsize/10^6 + 1, shape2 = libsize - 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]) * length/10^3 * libsize/10^6 + 1)

        # fit_deriv <- 0
        
        # sigma <- sqrt(fit_deriv^2*sigmat^2 + sigmab^2)
        # m <- 1/(exp(sigma^2) - 1) 
        # theta <- exp(log(2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])])) + 1/2*sigma^2)*1/m
  
        # prob <- sapply(1:length(missing_ca), function(i){dnbinom(x = round(2^missing_ca[i] * libsize/10^6 * length/10^3,0), size = m[i], prob = 1/(1+theta[i]), log=TRUE)}) 

        # res <- res - prob
        
        # res <- res - sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(sqrt(2*pi)*0.5) + (missing_ca - (delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]))^2/(2*0.5^2)
        
        # Mistake but doesn't change optim because fixed valued
        res <- res - sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(1/(sqrt(2*pi)*0.5)) + (missing_ca - (delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]))^2/(2*0.5^2)

        
        # ca_next <- data$caRNA[apply(data$caRNA, 2, function(col){which.min(is.na(col))})[which(is.na(data$caRNA[data$caRNA[,1]==0,]))]+nrow(data$caRNA)*(which(is.na(data$caRNA[data$caRNA[,1]==0,]))-1)] 
        # res <- res + sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(sqrt(2*pi)*0.1) + (missing_ca - ca_next)^2/(2*0.1^2)
      }
    }
    return(res)
  }, error=function(e){return(10^16)})
}
TSE_bs_cp <- cmpfun(TSE_bs)


TSE_bs_theo <- function(par, data, time_data, batch, results="sum", library_size, length, ...){
  # par:
  # fn_fit:
  # data: data of all replicates
  # time_data: 
  # batch: 
  # results: 
  # library_size: library size for each sample, need to be subsetted for that condition 
  # length:
  # ...: 
  
  # check if the batch has np data or not
  no_np <- FALSE
  if(batch == 'b2'){
    no_np <- TRUE
  }
  
  # extract the smoothing parameter and variance from the parameter list
  span <- par[1]
  sigmat <- par[2]
  sigmab <- par[3]
  
  par <- par[-c(1:3)]
  
  # keep al data and extract the data corresponding to the wanted batch
  data_rpkm_all <- data
  if(batch != 'all'){
    b <- as.numeric(gsub('b','',batch)) 
    data <- lapply(data_rpkm_all, function(l){l[,c(1,1+b)]})
  }
  
  # create data.frame of caRNA for the weight for the smooting
  dat_g_exp <- data.frame(Time=rep(data_rpkm_all$caRNA[,1], ncol(data_rpkm_all$caRNA)-1), RPKM = as.vector(data_rpkm_all$caRNA[,-1]), batch = rep(c("b1", "b2", "b3"), each=nrow(data_rpkm_all$caRNA)))
  
  tryCatch({
    # times to extrapolate, every second
    extra_t_end <- time_data[length(time_data)]
    extra_t_seq <- seq(time_data[1],extra_t_end + 1, by=1/60)
    
    # extract data to extrapolate
    data_extrapolation <- data$caRNA
    
    # replace missing t=0 by parameters 
    if(any(is.na(data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])]))){
      data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])] <- par[1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,]))]
      missing_ca <- par[c(1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
      par <- par[-c(1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
    }
    
    # remove time point if all na
    to_remove <- apply(data_extrapolation, 1, function(x){all(is.na(x[-1]))})
    data_extrapolation <- data_extrapolation[!to_remove,]
    
    # extract library size 
    # libsize <- as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''), as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))
    dat_g_exp_pattern <- paste0('chromatin', '.',sapply(nchar(as.character(dat_g_exp$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp$Time),'.',dat_g_exp$batch)
    # dat_g_exp_pattern_2 <- sapply(1:nrow(dat_g_exp),function(i){paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''), as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i])})
    # lib_indices <- sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = dat_g_exp_pattern[i], x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}})
    lib_indices <- match(dat_g_exp_pattern, row.names(library_size))
    libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[lib_indices]))
    
    
    # Calculate prob of p using a beta distribution: k ~ binom(p,N), p ~ beta(k+1, N-k+1) in log f is t*f(t)
    proba <- sapply(1:nrow(dat_g_exp), function(i){ 2^dat_g_exp$RPKM[i] * length / 10^3 / 10^6 * dbeta( x = 2^dat_g_exp$RPKM[i] * length / 10^3 / 10^6, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1)})
    
    # prob_rpkm_ub <- log2(sapply(1:nrow(dat_g_exp), function(i){10^6*10^3/length * qbeta( p = 0.25, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, lower.tail=F)}))
    # prob_rpkm_lb <- log2(sapply(1:nrow(dat_g_exp), function(i){10^6*10^3/length * qbeta( p = 0.25, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, lower.tail=T)}))
    
    res <- vector(mode = "list", length=ncol(data$caRNA)-1)
    # for each batch
    for (b in 2:(ncol(data$caRNA))){
      # extrapolate caRNA for specific batch 
      data_extra <- data.frame(x=as.vector(data_extrapolation[,1]), y=as.vector(data_extrapolation[,b]))
      data_extra <- data_extra[!is.na(data_extra$y),]
      
      # calculate weight as 1/IQR(p) saturating at 10
      if(ncol(data$caRNA) > 2){
        weight <- proba[dat_g_exp$batch == paste0('b',b-1)]
        # weight <- sapply(1:length(prob_rpkm_ub),function(i){1/max(0.1,prob_rpkm_ub[i]-prob_rpkm_lb[i])})[dat_g_exp$batch == paste0('b',b-1)] # don't go to low values
      }else{
        weight <- proba[dat_g_exp$batch == batch]
        # weight <- sapply(1:length(prob_rpkm_ub),function(i){1/max(0.1,prob_rpkm_ub[i]-prob_rpkm_lb[i])})[dat_g_exp$batch == batch] # don't go to low values
      }
      weight <- weight/max(weight, na.rm=T)
      # add max weigth to the missing t=0 values replaced by parameter (i.e. sure values)
      if(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))){
        weight[which(is.na(data$caRNA[data$caRNA[,1]==0,-1]))] <- 1
      }
      weight <- weight[!is.na(weight)]
      
      # order data to extrapolate and weight
      weight <- weight[order(data_extra$x)]
      data_extra <- data_extra[order(data_extra$x),]
      
      # extrapolate
      extrapolated_input <-  cbind(extra_t_seq, 2^predict(smooth.spline(x=data_extra, w=weight, spar = span), extra_t_seq)$y)
      
      # simulate data for every minutes to get derivative 
      fit  <- cbind(t=time_data, caRNA = data$caRNA[,b], npRNA=data$npRNA[,b], cytoRNA = data$cytoRNA[,b])
      
      # compare data to simulations
      for (cpt in names(data)){
        if(!all(is.na(data[[cpt]][,b]))){
          if(cpt == "caRNA"){
            cpt_name <- "chromatin"
          }else if(cpt == "npRNA"){
            cpt_name <- "nuc"
          }else if(cpt == "cytoRNA"){
            cpt_name <- "cytoplasmic"
          }
          # print(cbind(fit[match(data[[cpt]][,"t"], time_data), cpt], data[[cpt]]))
          
          # libsize <- as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0(cpt_name,'.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$batch == paste0('b',b-1)]
          
          dat_g_exp_pattern <- paste0(cpt_name, '.',sapply(nchar(as.character(dat_g_exp$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp$Time),'.', dat_g_exp$batch)
          # lib_indices <- sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = dat_g_exp_pattern[i], x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}})
          lib_indices <- match(dat_g_exp_pattern, row.names(library_size))
          libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[lib_indices]))[dat_g_exp$batch == paste0('b',b-1)]
          
          
          res[[b-1]] <- c(res[[b-1]], TSE_poisson_cp(fit[match(data[[cpt]][,"t"], time_data), cpt],  data[[cpt]][,b], libsize, length))
        }else{
          res[[b-1]] <- c(res[[b-1]], NA)
        }
      }
    }
    
    if (results != "detailed"){
      res <- sum(unlist(lapply(res, sum, na.rm=T)))
      # add prior of sigmat and sigmab ~ regularisation (2 times higher because can only be positive)
      # v1
      # res <- res + log(2*10*sqrt(2*pi)) + sigmat^2/(2*10^2) + log(2*0.5*sqrt(2*pi)) + sigmab^2/(2*0.5^2) + log(sqrt(2*pi)*0.1) + (span - 0.45)^2/(2*0.2^2)
      
      # v2
      res <- res - log(2*1/(5*sqrt(2*pi))) + sigmat^2/(2*5^2) - log(2/(0.1*sqrt(2*pi))) + sigmab^2/(2*0.1^2) - log(1/(sqrt(2*pi)*0.05)) + (span - 0.45)^2/(2*0.05^2)
      
      # add regul on b1_ca0 based the other ca_0
      if(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))){
        ca_other_0 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==0,-1] 
        ca_other_10 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==10,-1] 
        
        delta <- mean(ca_other_0 - ca_other_10, na.rm=T)
        
        # libsize <- mean(as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$Time == 0], na.rm=T)
        
        # proba <- pbeta( q = 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])])* length / 10^3 / 10^6, shape1 = 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]) * length/10^3 * libsize/10^6 + 1, shape2 = libsize - 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]) * length/10^3 * libsize/10^6 + 1)
        
        # fit_deriv <- 0
        
        # sigma <- sqrt(fit_deriv^2*sigmat^2 + sigmab^2)
        # m <- 1/(exp(sigma^2) - 1) 
        # theta <- exp(log(2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])])) + 1/2*sigma^2)*1/m
        
        # prob <- sapply(1:length(missing_ca), function(i){dnbinom(x = round(2^missing_ca[i] * libsize/10^6 * length/10^3,0), size = m[i], prob = 1/(1+theta[i]), log=TRUE)}) 
        
        # res <- res - prob
        res <- res - sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(1/(sqrt(2*pi)*0.5)) + (missing_ca - (delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]))^2/(2*0.5^2)
        
        # ca_next <- data$caRNA[apply(data$caRNA, 2, function(col){which.min(is.na(col))})[which(is.na(data$caRNA[data$caRNA[,1]==0,]))]+nrow(data$caRNA)*(which(is.na(data$caRNA[data$caRNA[,1]==0,]))-1)] 
        # res <- res + sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(sqrt(2*pi)*0.1) + (missing_ca - ca_next)^2/(2*0.1^2)
      }
    }
    return(res)
  }, error=function(e){return(10^16)})
}
TSE_bs_theo_cp <- cmpfun(TSE_bs_theo)

TSE_bs_theo_bt <- function(par, fn_fit, data, time_data, batch, results="sum", library_size, length, ...){
  # par:
  # fn_fit:
  # data: data of all replicates
  # time_data: 
  # batch: 
  # results: 
  # library_size: library size for each sample, need to be subsetted for that condition 
  # length:
  # ...: 
  
  # check if the batch has np data or not
  no_np <- FALSE
  if(batch == 'b2'){
    no_np <- TRUE
  }
  
  # extract the smoothing parameter and variance from the parameter list
  span <- par[1]
  sigmat <- par[2]
  sigmab <- par[3]
  
  par <- par[-c(1:3)]
  
  # keep al data and extract the data corresponding to the wanted batch
  data_rpkm_all <- data
  if(batch != 'all'){
    b <- as.numeric(gsub('b','',batch)) 
    data <- lapply(data_rpkm_all, function(l){l[,c(1,1+b)]})
  }
  
  # create data.frame of caRNA for the weight for the smooting
  dat_g_exp <- data.frame(Time=rep(data_rpkm_all$caRNA[,1], ncol(data_rpkm_all$caRNA)-1), RPKM = as.vector(data_rpkm_all$caRNA[,-1]), batch = rep(c("b1", "b2", "b3"), each=nrow(data_rpkm_all$caRNA)))
  
  tryCatch({
    # times to extrapolate, every second
    extra_t_end <- time_data[length(time_data)]
    extra_t_seq <- seq(time_data[1],extra_t_end + 1, by=1/60)
    
    # extract data to extrapolate
    data_extrapolation <- data$caRNA
    
    # replace missing t=0 by parameters 
    if(any(is.na(data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])]))){
      data_extrapolation[data_extrapolation[,1]==0, is.na(data_extrapolation[data_extrapolation[,1]==0,])] <- par[1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,]))]
      missing_ca <- par[c(1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
      par <- par[-c(1:sum(is.na(data_extrapolation[data_extrapolation[,1]==0,])))]
    }
    
    # remove time point if all na
    to_remove <- apply(data_extrapolation, 1, function(x){all(is.na(x[-1]))})
    data_extrapolation <- data_extrapolation[!to_remove,]
    
    # extract library size 
    dat_g_exp_pattern <- paste0('chromatin', '.',sapply(nchar(as.character(dat_g_exp$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp$Time),'.',dat_g_exp$batch)
    # dat_g_exp_pattern_2 <- sapply(1:nrow(dat_g_exp),function(i){paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''), as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i])})
    # lib_indices <- sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = dat_g_exp_pattern[i], x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}})
    lib_indices <- match(dat_g_exp_pattern, row.names(library_size))
    libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[lib_indices]))
    
    # libsize <- as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''), as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))
    
    # Calculate prob of p using a beta distribution: k ~ binom(p,N), p ~ beta(k+1, N-k+1) in log f is t*f(t)
    proba <- sapply(1:nrow(dat_g_exp), function(i){ 2^dat_g_exp$RPKM[i] * length / 10^3 / 10^6 * dbeta( x = 2^dat_g_exp$RPKM[i] * length / 10^3 / 10^6, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1)})
    
    # prob_rpkm_ub <- log2(sapply(1:nrow(dat_g_exp), function(i){10^6*10^3/length * qbeta( p = 0.25, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, lower.tail=F)}))
    # prob_rpkm_lb <- log2(sapply(1:nrow(dat_g_exp), function(i){10^6*10^3/length * qbeta( p = 0.25, shape1 = 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, shape2 = libsize[i] - 2^dat_g_exp$RPKM[i] * length/10^3 * libsize[i]/10^6 + 1, lower.tail=T)}))
    
    res <- vector(mode = "list", length=ncol(data$caRNA)-1)
    # for each batch
    for (b in 2:(ncol(data$caRNA))){
      # extrapolate caRNA for specific batch 
      data_extra <- data.frame(x=as.vector(data_extrapolation[,1]), y=as.vector(data_extrapolation[,b]))
      data_extra <- data_extra[!is.na(data_extra$y),]
      
      # calculate weight as 1/IQR(p) saturating at 10
      if(ncol(data$caRNA) > 2){
        weight <- proba[dat_g_exp$batch == paste0('b',b-1)]
        # weight <- sapply(1:length(prob_rpkm_ub),function(i){1/max(0.1,prob_rpkm_ub[i]-prob_rpkm_lb[i])})[dat_g_exp$batch == paste0('b',b-1)] # don't go to low values
      }else{
        weight <- proba[dat_g_exp$batch == batch]
        # weight <- sapply(1:length(prob_rpkm_ub),function(i){1/max(0.1,prob_rpkm_ub[i]-prob_rpkm_lb[i])})[dat_g_exp$batch == batch] # don't go to low values
      }
      weight <- weight/max(weight, na.rm=T)
      # add max weigth to the missing t=0 values replaced by parameter (i.e. sure values)
      if(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))){
        weight[which(is.na(data$caRNA[data$caRNA[,1]==0,-1]))] <- 1
      }
      weight <- weight[!is.na(weight)]
      
      # order data to extrapolate and weight
      weight <- weight[order(data_extra$x)]
      data_extra <- data_extra[order(data_extra$x),]
      
      # extrapolate
      extrapolated_input <-  cbind(extra_t_seq, 2^predict(smooth.spline(x=data_extra, w=weight, spar = span), extra_t_seq)$y)
      
      # simulate data for every minutes to get derivative 
      fit <- fn_fit(t = seq(time_data[1], time_data[length(time_data)], by=1), par, extrapolated_input = extrapolated_input, no_np = no_np, ...)
      if (any(is.nan(fit) | is.na(fit) |  fit < 0 | nrow(fit)<length(seq(time_data[1], time_data[length(time_data)], by=1)))) {
        # print(p)
        # print(fit)
        return(10^16)
      }
      fit[,-1] <- log2(fit[,-1])
      
      # calculate the derivative (not using the model, because if quassi steady state derivative from model will be null)
      fit_deriv_ca <- data.frame(Time = extrapolated_input[,1], caRNA=c((extrapolated_input[2,2]-extrapolated_input[1,2])/(extrapolated_input[2,1]-extrapolated_input[1,1]), 
                                                                        0.5*((extrapolated_input[,2][-c(nrow(extrapolated_input)-1:0)]-extrapolated_input[,2][-c(1,nrow(extrapolated_input))])/(extrapolated_input[,1][-c(nrow(extrapolated_input)-1:0)]-extrapolated_input[,1][-c(1,nrow(extrapolated_input))])+(extrapolated_input[,2][-c(1:2)]-extrapolated_input[,2][-c(1,nrow(extrapolated_input))])/(extrapolated_input[,1][-c(1:2)]-extrapolated_input[,1][-c(1,nrow(extrapolated_input))])),
                                                                        (extrapolated_input[nrow(extrapolated_input),2]-extrapolated_input[nrow(extrapolated_input)-1,2])/(extrapolated_input[nrow(extrapolated_input),1]-extrapolated_input[nrow(extrapolated_input)-1,1])))
      
      fit_deriv <- data.frame(Time = fit_deriv_ca$Time, caRNA = fit_deriv_ca$caRNA/extrapolated_input[,2])[(seq(min(time_data),max(time_data)) - time_data[1])*60+1,]
      
      for (col in colnames(fit)){
        if (col == "npRNA"){
          fit_deriv_col <- log(2)*data.frame(npRNA=c((fit[2,col]-fit[1,col])/(fit[2,1]-fit[1,1]), 
                                                     0.5*((fit[,col][-c(nrow(fit)-1:0)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(nrow(fit)-1:0)]-fit[,1][-c(1,nrow(fit))])+
                                                            (fit[,col][-c(1:2)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(1:2)]-fit[,1][-c(1,nrow(fit))])),
                                                     (fit[nrow(fit),col]-fit[nrow(fit)-1,col])/(fit[nrow(fit),1]-fit[nrow(fit)-1,1])))
          fit_deriv <- cbind(fit_deriv, npRNA = fit_deriv_col)
        }
        if (col == "cytoRNA"){
          fit_deriv_col <- log(2)*data.frame(cytoRNA=c((fit[2,col]-fit[1,col])/(fit[2,1]-fit[1,1]), 
                                                       0.5*((fit[,col][-c(nrow(fit)-1:0)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(nrow(fit)-1:0)]-fit[,1][-c(1,nrow(fit))])+
                                                              (fit[,col][-c(1:2)]-fit[,col][-c(1,nrow(fit))])/(fit[,1][-c(1:2)]-fit[,1][-c(1,nrow(fit))])),
                                                       (fit[nrow(fit),col]-fit[nrow(fit)-1,col])/(fit[nrow(fit),1]-fit[nrow(fit)-1,1])))
          fit_deriv <- cbind(fit_deriv, cytoRNA = fit_deriv_col)
        }
      }
      
      # fit[,-1][fit[,-1] <= -20] <- -20
      
      # extract time matching time data
      fit <- fit[match(time_data, fit[,1]),]  
      fit <- cbind(fit, caRNA=log2(extrapolated_input[(time_data - time_data[1])*60+1,2]))
      
      fit_deriv <- fit_deriv[match(time_data, fit_deriv[,1]),]  
      
      # reset the points to the correct ones
      fit  <- cbind(t=time_data, caRNA = data$caRNA[,b], npRNA=data$npRNA[,b], cytoRNA = data$cytoRNA[,b])
      
      
      # compare data to simulations
      for (cpt in names(data)){
        if(!all(is.na(data[[cpt]][,b]))){
          if(cpt == "caRNA"){
            cpt_name <- "chromatin"
          }else if(cpt == "npRNA"){
            cpt_name <- "nuc"
          }else if(cpt == "cytoRNA"){
            cpt_name <- "cytoplasmic"
          }
          # print(cbind(fit[match(data[[cpt]][,"t"], time_data), cpt], data[[cpt]]))
          dat_g_exp_pattern <- paste0(cpt_name, '.',sapply(nchar(as.character(dat_g_exp$Time)),function(x){paste0(rep(0,3-x),collapse='')}), as.character(dat_g_exp$Time),'.', dat_g_exp$batch)
          # lib_indices <- sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = dat_g_exp_pattern[i], x = row.names(library_size), ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}})
          lib_indices <- match(dat_g_exp_pattern, row.names(library_size))
          libsize <- as.integer(floor((library_size[,"lib.size"]*library_size[,"norm.factors"])[lib_indices]))[dat_g_exp$batch == paste0('b',b-1)]
          
          # libsize <- as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0(cpt_name,'.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$batch == paste0('b',b-1)]
          
          res[[b-1]] <- c(res[[b-1]], TSE_neg_binom_cp(fit[match(data[[cpt]][,"t"], time_data), cpt], fit_deriv[match(data[[cpt]][,"t"], time_data), cpt],  data[[cpt]][,b], libsize, length, sigmab, sigmat))
        }else{
          res[[b-1]] <- c(res[[b-1]], NA)
        }
      }
    }
    
    if (results != "detailed"){
      res <- sum(unlist(lapply(res, sum, na.rm=T)))
      # add prior of sigmat and sigmab ~ regularisation (2 times higher because can only be positive)
      # v1
      # res <- res + log(2*10*sqrt(2*pi)) + sigmat^2/(2*10^2) + log(2*0.5*sqrt(2*pi)) + sigmab^2/(2*0.5^2) + log(sqrt(2*pi)*0.1) + (span - 0.45)^2/(2*0.2^2)
      
      # v2
      res <- res - log(2*1/(5*sqrt(2*pi))) + sigmat^2/(2*5^2) - log(2/(0.1*sqrt(2*pi))) + sigmab^2/(2*0.1^2) - log(1/(sqrt(2*pi)*0.05)) + (span - 0.45)^2/(2*0.05^2)
      
      # add regul on b1_ca0 based the other ca_0
      if(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))){
        ca_other_0 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==0,-1] 
        ca_other_10 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==10,-1] 
        
        delta <- mean(ca_other_0 - ca_other_10, na.rm=T)
        
        # libsize <- mean(as.integer(floor((library_size$lib.size*library_size$norm.factors)[unlist(sapply(1:nrow(dat_g_exp),function(i){tmp <- grep(pattern = paste0('chromatin.*', '.',paste0(rep(0,3-nchar(as.character(dat_g_exp$Time[i]))), collapse = ''),  as.character(dat_g_exp$Time[i]),'.',dat_g_exp$batch[i]), x = library_size$samples, ignore.case = T); if (length(tmp)==0){return(NA)}else{return(tmp)}}))]))[dat_g_exp$Time == 0], na.rm=T)
        
        # proba <- pbeta( q = 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])])* length / 10^3 / 10^6, shape1 = 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]) * length/10^3 * libsize/10^6 + 1, shape2 = libsize - 2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]) * length/10^3 * libsize/10^6 + 1)
        
        # fit_deriv <- 0
        
        # sigma <- sqrt(fit_deriv^2*sigmat^2 + sigmab^2)
        # m <- 1/(exp(sigma^2) - 1) 
        # theta <- exp(log(2^(delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])])) + 1/2*sigma^2)*1/m
        
        # prob <- sapply(1:length(missing_ca), function(i){dnbinom(x = round(2^missing_ca[i] * libsize/10^6 * length/10^3,0), size = m[i], prob = 1/(1+theta[i]), log=TRUE)}) 
        
        # res <- res - prob
        
        res <- res - sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(1/sqrt(2*pi)*0.5) + (missing_ca - (delta + data$caRNA[data$caRNA[,1]==10, is.na(data$caRNA[data$caRNA[,1]==0,])]))^2/(2*0.5^2)
        
        # ca_next <- data$caRNA[apply(data$caRNA, 2, function(col){which.min(is.na(col))})[which(is.na(data$caRNA[data$caRNA[,1]==0,]))]+nrow(data$caRNA)*(which(is.na(data$caRNA[data$caRNA[,1]==0,]))-1)] 
        # res <- res + sum(any(is.na(data$caRNA[data$caRNA[,1]==0, is.na(data$caRNA[data$caRNA[,1]==0,])]))) * log(sqrt(2*pi)*0.1) + (missing_ca - ca_next)^2/(2*0.1^2)
      }
    }
    return(res)
  }, error=function(e){return(10^16)})
}
TSE_bs_theo_bt_cp <- cmpfun(TSE_bs_theo_bt)



# Total sum of square
TSE_raw <- function(fit, data){
  tmp <- (data[, -1] - fit) ^ 2
  MSE <- sum(tmp , na.rm = T) 
  return(MSE)
}
TSE_raw_cp <- cmpfun(TSE_raw)

# Error using negative binomial model 
TSE_neg_binom <- function(fit, fit_deriv, data, libsize, length, sigmab, sigmat){
  sigma <- sqrt(fit_deriv^2*sigmat^2 + sigmab^2)
  m <- 1/(exp(sigma^2) - 1) 
  m[exp(sigma^2) == 1] <- (1/sigma^2)[exp(sigma^2) == 1] # approximation if sigma numerically to close to zero
  
  # theta <- exp(log(2^fit) + log(libsize/10^6*length/10^3) + 1/2*sigma^2)*1/m
  mu <- exp(log(2^fit) + log(libsize/10^6*length/10^3) + 1/2*sigma^2)
  
  # prob <- sapply(1:length(data), function(i){dnbinom(x = round(2^data[i] * libsize[i]/10^6 * length/10^3,0), size = m[i], prob = ifelse(1/(1+theta[i]) == 1, 1-theta[i] , 1/(1+theta[i])), log=TRUE)}) 
  prob <- sapply(1:length(data), function(i){dnbinom(x = round(2^data[i] * libsize[i]/10^6 * length/10^3,0), size = m[i], mu = mu[i], log=TRUE)}) 
  
  # prob[m == Inf & theta == 0 & prob == 0] <- NaN
  
  # convert to rpkm prob
  # prob <- log( log(2) * (2^data[i] * libsize[i]/10^6 * length/10^3) * exp(prob))
  prob <- log(log(2)) + log(2^data * libsize/10^6 * length/10^3) + prob

  
  prob[is.nan(prob)] <- -Inf
  ll <- sum(-prob, na.rm=T)
  # print(prob)
  # print(ll)
  return(ll)
}
TSE_neg_binom_cp <- cmpfun(TSE_neg_binom)


TSE_poisson <- function(fit, data, libsize, length){
  
  prob <- sapply(1:length(data), function(i){dpois(x = round(2^data[i] * libsize[i]/10^6 * length/10^3,0), lambda = round(2^data[i] * libsize[i]/10^6 * length/10^3,0), log=TRUE)}) 
  prob <- log(log(2)) + log(2^data * libsize/10^6 * length/10^3) + prob
  
  #  prob[m == Inf & theta == 0 & prob == 0] <- NaN
  
  prob[is.nan(prob)] <- -Inf
  ll <- sum(-prob, na.rm=T)
  # print(prob)
  # print(ll)
  return(ll)
}
TSE_poisson_cp <- cmpfun(TSE_poisson)



TSE_bs_sa_fixed <- function(par, fixed, fixed_value, fixed_name, fn_fit, data_all, data_rpkm, time_data, batch, results="sum", ...){
  if (fixed_name != "k1'/k2" & fixed_name != "k1'k2'/k2"){
    par <- c(par[1:(3+sum(is.na(data_rpkm$caRNA[1,])) + fixed-1)], fixed_value, par[(3+sum(is.na(data_rpkm$caRNA[1,]))+fixed):length(par)])
  }else if (fixed_name == "k1'/k2" ){
    par <- c(par[1:(3+sum(is.na(data_rpkm$caRNA[1,])))], fixed_value + par[3+sum(is.na(data_rpkm$caRNA[1,]))+1], par[(3+sum(is.na(data_rpkm$caRNA[1,]))+1):length(par)])
  }else if( fixed_name == "k1'k2'/k2" ) {
    par <- c(par[1:(3+sum(is.na(data_rpkm$caRNA[1,])))], fixed_value + par[3+sum(is.na(data_rpkm$caRNA[1,]))+1] - par[3+sum(is.na(data_rpkm$caRNA[1,]))+2] , par[(3+sum(is.na(data_rpkm$caRNA[1,]))+1):length(par)])
  }
  
  if (batch != 'b2'){
    par_new <- c(par[1:(3+sum(is.na(data_rpkm$caRNA[1,])))],parameters_convert_cp(par[3+sum(is.na(data_rpkm$caRNA[1,]))+1:4],log_in = T, log_out = T, direction = 1, no_np = F))
    res <- TSE_bs_cp(par = par_new, fn_fit = fn_fit, data = data_all, time_data = time_data, batch=batch, results = results, ...)
    return(res)
  }else{
    par_new <- c(par[1:(3+sum(is.na(data_rpkm$caRNA[1,])))],parameters_convert_cp(par[3+sum(is.na(data_rpkm$caRNA[1,]))+1:3],log_in = T, log_out = T, direction = 1, no_np = T))
    res <- TSE_bs_cp(par = par_new, fn_fit = fn_fit, data = data_all, time_data = time_data, batch=batch, results = results, ...)
    return(res)
  }
}
TSE_bs_sa_fixed_cp <- cmpfun(TSE_bs_sa_fixed)


## optimisation wrappers
optim_all_bs <- function(n_init, obj_fn, data_rpkm, batch, save_file, length, libsize, ...){
   data_rpkm_all <- data_rpkm
    
  if(batch != 'all'){
    b <- as.numeric(gsub('b','',batch)) 
    data_rpkm <- lapply(data_rpkm_all, function(l){l[,c(1,1+b)]})
  }
    
  # resume if took too long and not finished
  new_table <- TRUE
  if (file.exists(save_file)){
    if(file.size(save_file) != 0){
      new_table <- FALSE
      
      list_object <- c(ls(),"list_object")
      load(save_file)
      list_object_new <- ls()
      
      if(list_object_new[!(list_object_new %in% list_object)] == "bfgs"){
        best_all <- bfgs[[batch]]
      }else if (list_object_new[!(list_object_new %in% list_object)] == "best_all") {
      }
      n_init <- length(best_all$value)
    }
  }
    
  if (new_table){
    if(batch != 'b2'){
      best_all <- list(value=as.numeric(rep(NA, n_init)), 
                       start_value = as.numeric(rep(NA, n_init)), 
                       seed = as.integer(rep(NA, n_init)), 
                       start_param = matrix(nrow = n_init, ncol=7+sum(is.na(data_rpkm$caRNA[1,]))), 
                       parameter = matrix(nrow = n_init, ncol=7+sum(is.na(data_rpkm$caRNA[1,]))), 
                       convergence=as.numeric(rep(NA, n_init)),
                       time = as.numeric(rep(NA, n_init)))
    }else{
      best_all <- list(value=as.numeric(rep(NA, n_init)), 
                       start_value = as.numeric(rep(NA, n_init)), 
                       seed = as.integer(rep(NA, n_init)), 
                       start_param = matrix(nrow = n_init, ncol=6+sum(is.na(data_rpkm$caRNA[1,]))), 
                       parameter = matrix(nrow = n_init, ncol=6+sum(is.na(data_rpkm$caRNA[1,]))), 
                       convergence=as.numeric(rep(NA, n_init)),
                       time = as.numeric(rep(NA, n_init)))
    }
  }
    

  while (sum(is.na(best_all$value)) != 0 ) { 
    
    print(which(is.na(best_all$value))[1])
    
    seed <- floor(runif(1) * 10^8)
    set.seed(seed)
    
    if( batch != 'b2'){
      p <- vector(length=7+sum(is.na(data_rpkm$caRNA[1,])))
    }else{
      p <- vector(length=6+sum(is.na(data_rpkm$caRNA[1,])))
    }
      
    p[1] <- rnorm(1,0.45,0.05) # spar
    p[2] <- abs(rnorm(1,0,5)) # sigmat
    p[3] <- abs(rnorm(1,0,0.1)) # sigmab
    
    if(any(is.na(data_rpkm$caRNA[1,]))){i
      ca_other_0 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==0,-1]
      ca_other_10 <- data_rpkm_all$caRNA[data_rpkm_all$caRNA[,1]==10,-1]

      delta <- mean(ca_other_0 - ca_other_10, na.rm=T)
      p[3+1:sum(is.na(data_rpkm$caRNA[1,]))] <- sapply(which(is.na(data_rpkm$caRNA[1,])), function(idx){data_rpkm$caRNA[,idx][which(!is.na(data_rpkm$caRNA[,idx]))[1]]}) + delta + rnorm(mean = 0, sd = 0.5, n =sum(is.na(data_rpkm$caRNA[1,])))
    }
    
    if( batch != 'b2'){
      p[3+sum(is.na(data_rpkm$caRNA[1,])) + 1:4]  <- runif(4, -5, 5)
    }else{
      p[3+sum(is.na(data_rpkm$caRNA[1,])) + 1:3]  <- runif(3, -5, 5)
    }

    init_val <- obj_fn(par = p, fn_fit = fn_all_ode_cp, data = data_rpkm_all, convert_param = T, log_in = T, u=input_fun_cp, batch = batch, length = length, library_size=libsize, ... )
    
    if (!is.null(init_val) & init_val != 10^16){
      cpt_1 <- proc.time()
        
      fit_all_tot <- NULL
      tryCatch({
        fit_all_tot <- optim(par = p, fn = obj_fn, fn_fit = fn_all_ode_cp,  data = data_rpkm_all, convert_param = T, log_in = T,  method = "BFGS", control = list(maxit=1000, trace=5, REPORT=10), u=input_fun_cp, batch = batch, length=length, library_size=libsize, ...)  # print(fit_np_tot$value)
      }, error = function(e){print('Error for initial parameter:'); print(p); fit_all_tot <- NULL})
      
      cpt_2 <- proc.time()
      
      if(is.null(fit_all_tot)){
      }else if (fit_all_tot$value == 10^16){
      }else{
        best_all$start_value[which(is.na(best_all$value))[1]] <- init_val
        best_all$seed[which(is.na(best_all$value))[1]] <- seed
        best_all$start_param[which(is.na(best_all$value))[1],] <- p
        best_all$parameter[which(is.na(best_all$value))[1],] <- fit_all_tot$par
        best_all$convergence[which(is.na(best_all$value))[1]] <- fit_all_tot$convergence
        best_all$time[which(is.na(best_all$value))[1]] <- (cpt_2-cpt_1)[3]
        best_all$value[which(is.na(best_all$value))[1]] <- fit_all_tot$value

        save(best_all, file=save_file)
      }    
    }else{
    }
  }
  return(best_all)
}
optim_all_bs_cp <- cmpfun(optim_all_bs)

# second optim
second_optim_all_bs <- function(n_init, obj_fn, data_rpkm, batch, first_best, ...){
  # if(batch == "b1"){
  #   best_all <- list(value=vector(length=n_init), parameter = matrix(nrow = n_init, ncol=6), convergence=vector(length=n_init))
  # }else{
  best_all <- list(value=vector(length=n_init), parameter = matrix(nrow = n_init, ncol=5), convergence=vector(length=n_init))
  # }
  # pp <- signal:::pchip(x=data_rpkm_all$caRNA[,1], y = rowMeans(data_rpkm_all$caRNA[,-1], na.rm=T))
  
  # optim span of loess such that sum of mse by batch is min
  
  # span <- 0.65
  # span_optim <- optim(span, MSE_span, data=data_rpkm_all, method="Nelder-Mead")
  
  # loess_caRNA <- loess(caRNA ~ Time,data = data.frame(Time=as.vector(matrix(rep(data_rpkm_all$caRNA[,1], 3),byrow = F)), caRNA=as.vector(data_rpkm_all$caRNA[,-1]), Batch=as.vector(matrix(c("b1", "b2", "b3"),ncol=col(data_rpkm_all$caRNA[,-1]),byrow = T))), span=1)
  # lines(predict(loess_caRNA, 0:120))
  
  # tmp<-sapply(seq(0.1,1,by=0.01), function(s){MSE_span(s, data_rpkm_all)})
  
  i <- 0
  while (i < n_init ) { 
    i <- i + 1
    
    previous <- first_best[[batch]]$parameter[order(first_best[[batch]]$value, decreasing = F), ][i,]
    # print(previous)
    best <- list(value=first_best[[batch]]$value[order(first_best[[batch]]$value, decreasing = F)][i], parameter = previous, convergence=first_best[[batch]]$convergence[order(first_best[[batch]]$value, decreasing = F)][i])
    
    j <- 0
    
    while (j < 10 ) { 
      j <- j + 1
      
      # if(batch == "b1"){
      #   p <- vector(length=6)
      #   p <- c(runif(1,0.2,1), runif(1,-5,5),runif(4, -3, 3))
      # } else {
      p <- vector(length=5)
      p <- c(runif(1,0.2,1),log10(runif(4, 0.5, 2)) + previous[-1])
      # }
      
      fit_all_tot <- tryCatch({
        optim(par = p, fn = obj_fn, fn_fit = fn_all_ode,  data = data_rpkm, convert_param = T, log_in = T,  method = "BFGS", control = list(maxit=1000, trace=5, REPORT=1000), u=interpolate_caRNA, batch = batch, ...)  # print(fit_np_tot$value)
      }, error = function(e){e; print('Error for initial parameter:'); print(p); print(e); return(NULL)})
      
      if(is.null(fit_all_tot) | fit_all_tot$value == 10^16){
        j <- j-1
      }else if (fit_all_tot$value < best$value){
        # print(fit_all_tot)
        best$value <- fit_all_tot$value
        best$parameter <- fit_all_tot$par
        best$convergence <- fit_all_tot$convergence
      }    
    }
    
    # fit_all_tot <- tryCatch({
    #   optim(par = previous, fn = obj_fn, fn_fit = fn_all_ode,  data = data_rpkm, convert_param = T, log_in = T,  method = "BFGS", control = list(maxit=1000, trace=5, REPORT=1000), u=interpolate_caRNA, batch = batch, ...)  # print(fit_np_tot$value)
    # }, error = function(e){e; print('Error for initial parameter:'); print(p); print(e); return(NULL)})
    
    best_all$value[i] <- best$value
    best_all$parameter[i,] <- best$par
    best_all$convergence[i] <- best$convergence
  }
  return(best_all)
}

## extract gene data
create_dat_g_rpkm <- function(rpkm, gene, batch = sort(unique(unlist(lapply(strsplit(colnames(rpkm), ".", fixed = T), function(x){x[4]})))), time = sort(unique(unlist(lapply(strsplit(colnames(rpkm), ".", fixed = T), function(x){x[3]}))))){
  dat_g_rpkm <- data.frame()
  for (t in time) {
    tmp <- c()     
    for (b in batch) {
      if (length(grep(t, grep(b, colnames(rpkm), value = T))) == 0) {
        tmp <- c(tmp, NA)
      }else{
        tmp <- c(tmp, rpkm[row.names(rpkm) == gene, grep(t, grep(b, colnames(rpkm), value = T), value = T)])
      }
    }
    dat_g_rpkm <- rbind(dat_g_rpkm, data.frame(t = as.numeric(t), t(tmp)))
  }
  names(dat_g_rpkm) <- c('t', batch)
  return(dat_g_rpkm)
}
create_dat_g_rpkm_cp <- cmpfun(create_dat_g_rpkm)
