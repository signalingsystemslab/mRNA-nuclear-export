profile_dMod_with_save <- function (obj, pars, whichPar, alpha = 0.05, limits = c(lower = -Inf, 
                                                       upper = Inf), method = c("integrate", "optimize"), stepControl = NULL, 
         algoControl = NULL, optControl = NULL, verbose = FALSE, cores = 1, file_name=NULL,
         ...) 
{
  dotArgs <- list(...)
  sanePars <- dMod:::sanitizePars(pars, dotArgs$fixed)
  pars <- sanePars$pars
  fixed <- sanePars$fixed
  method <- match.arg(method)
  if (method == "integrate") {
    sControl <- list(stepsize = 1e-04, min = 1e-04, max = Inf, 
                     atol = 0.01, rtol = 0.01, limit = 500, stop = "value")
    aControl <- list(gamma = 1, W = "hessian", reoptimize = FALSE, 
                     correction = 1, reg = .Machine$double.eps)
    oControl <- list(rinit = 0.1, rmax = 10, iterlim = 10, 
                     fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  }
  if (method == "optimize") {
    sControl <- list(stepsize = 0.01, min = 1e-04, max = Inf, 
                     atol = 0.1, rtol = 0.1, limit = 100, stop = "value")
    aControl <- list(gamma = 0, W = "identity", reoptimize = TRUE, 
                     correction = 1, reg = 0)
    oControl <- list(rinit = 0.1, rmax = 10, iterlim = 100, 
                     fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  }
  cores <- min(length(whichPar), cores)
  cores <- dMod:::sanitizeCores(cores)
  if (!is.null(stepControl)) 
    sControl[match(names(stepControl), names(sControl))] <- stepControl
  if (!is.null(algoControl)) 
    aControl[match(names(algoControl), names(aControl))] <- algoControl
  if (!is.null(optControl)) 
    oControl[match(names(optControl), names(oControl))] <- optControl
  if (cores > 1) {
    if (Sys.info()[["sysname"]] == "Windows") {
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl = cluster)
      parallel::clusterCall(cl = cluster, function(x) .libPaths(x), 
                            .libPaths())
      varlist <- ls()
      varlist <- c("obj", "whichPar", "alpha", "limits", 
                   "method", "verbose", "cores", "pars", "fixed", 
                   "dotArgs", "sControl", "aControl", "oControl")
      parallel::clusterExport(cluster, envir = environment(), 
                              varlist = varlist)
    }
    else {
      doParallel::registerDoParallel(cores = cores)
    }
    "%mydo%" <- foreach::"%dopar%"
  }
  else {
    "%mydo%" <- foreach::"%do%"
  }
  if (is.character(whichPar)) 
    whichPar <- which(names(pars) %in% whichPar)
  
  loaded_packages <- .packages()
  out <- foreach::foreach(whichIndex = whichPar, .packages = loaded_packages, 
                          .inorder = TRUE, .options.multicore = list(preschedule = FALSE)) %mydo% 
    {
      if(!is.null(file_name)){
        save_iteration <- TRUE
      }
      if(file.exists(file_name)){
        load(file_name, envir = environment())
        environment(init) <- environment(continue) <- environment()
        environment(obj.prof) <- environment(pseudoinverse) <- environment(constraint) <- environment(lagrange) <- environment()
        environment(doIteration) <- environment(doAdaption) <- environment()
        
        if ( direction == 1){
          continue()
          
          direction <- -1 
          
          init()
          continue()
          
        }else if (direction == -1){
          continue()
        }
        
        out <- as.data.frame(out)
        out$whichPar <- whichPar.name
        parframe(out, parameters = names(pars), metanames = c("value", 
                                                              "constraint", "stepsize", "gamma", "whichPar"), 
                 obj.attributes = names(out.attributes))
        
      }else{
        loadDLL(obj)
        whichPar.name <- names(pars)[whichIndex]
        obj.opt <- obj
        obj.prof <- function(p, ...) {
          out <- obj(p, ...)
          Id <- diag(1/.Machine$double.eps, length(out$gradient))
          Id[whichIndex, whichIndex] <- 1
          colnames(Id) <- rownames(Id) <- names(out$gradient)
          W <- match.arg(aControl$W[1], c("hessian", "identity"))
          out$hessian <- switch(W, hessian = out$hessian, 
                                identity = Id)
          return(out)
        }
        pseudoinverse <- function(m, tol) {
          msvd <- svd(m)
          index <- which(abs(msvd$d) > max(dim(m)) * max(msvd$d) * 
                           tol)
          if (length(index) == 0) {
            out <- array(0, dim(m)[2:1])
          }
          else {
            out <- msvd$u[, index] %*% (1/msvd$d[index] * 
                                          t(msvd$v)[index, ])
          }
          attr(out, "valid") <- 1:length(msvd$d) %in% index
          return(out)
        }
        constraint <- function(p) {
          value <- p[whichIndex] - pars[whichIndex]
          gradient <- rep(0, length(p))
          gradient[whichIndex] <- 1
          return(list(value = value, gradient = gradient))
        }
        lagrange <- function(y) {
          p <- y
          lambda <- 0
          out <- do.call(obj.prof, c(list(p = p), dotArgs))
          g.original <- constraint(p)
          g <- direction * g.original$value
          gdot <- direction * g.original$gradient
          ldot <- out$gradient
          lddot <- out$hessian
          M <- rbind(cbind(lddot, gdot), matrix(c(gdot, 
                                                  0), nrow = 1))
          v <- c(-rep(gamma, length(p)) * (ldot + lambda * 
                                             gdot), 1)
          v0 <- c(-rep(0, length(p)) * (ldot + lambda * 
                                          gdot), 1)
          W <- pseudoinverse(M, tol = aControl$reg)
          valid <- attr(W, "valid")
          if (any(!valid)) {
            dy <- try(as.vector(W %*% v)[1:length(p)], 
                      silent = FALSE)
            dy0 <- try(as.vector(W %*% v0)[1:length(p)], 
                       silent = FALSE)
            dy[!valid[1:length(p)]] <- dy0[!valid[1:length(p)]] <- 0
            dy[whichIndex] <- dy0[whichIndex] <- direction
            warning(paste0("Iteration ", i, ": Some singular values of the Hessian are below the threshold. Optimization will be performed."))
          }
          else {
            dy <- try(as.vector(W %*% v)[1:length(p)], 
                      silent = FALSE)
            dy0 <- try(as.vector(W %*% v0)[1:length(p)], 
                       silent = FALSE)
          }
          if (!inherits(dy, "try-error")) {
            names(dy) <- names(y)
            correction <- sqrt(sum((dy - dy0)^2))/sqrt(sum(dy^2))
          }
          else {
            dy <- NA
            correction <- 0
            warning(paste0("Iteration ", i, ": Impossible to invert Hessian. Trying to optimize instead."))
          }
          out.attributes <- attributes(out)[sapply(attributes(out), 
                                                   is.numeric)]
          out.attributes.names <- names(out.attributes)
          return(c(list(dy = dy, value = out$value, gradient = out$gradient, 
                        correction = correction, valid = valid, attributes = out.attributes.names), 
                   out.attributes))
        }
        doIteration <- function() {
          optimize <- aControl$reoptimize
          if (is.na(dy[1])) {
            optimize <- TRUE
            y.try <- y
            y.try[whichIndex] <- y[whichIndex] + direction * 
              stepsize
            rinit <- oControl$rinit
          }
          else {
            dy.norm <- sqrt(sum(dy^2))
            rinit <- min(c(oControl$rinit, 3 * dy.norm))
            y.try <- y + dy
            if (any(!lagrange.out$valid)) 
              optimize <- TRUE
          }
          if (optimize) {
            parinit.opt <- y.try[-whichIndex]
            fixed.opt <- c(fixed, y.try[whichIndex])
            arglist <- c(list(objfun = obj.opt, parinit = parinit.opt, 
                              fixed = fixed.opt, rinit = rinit), oControl[names(oControl) != 
                                                                            "rinit"], dotArgs[names(dotArgs) != "fixed"])
            myfit <- try(do.call(trust, arglist), silent = FALSE)
            if (!inherits(myfit, "try-error")) {
              y.try[names(myfit$argument)] <- as.vector(myfit$argument)
            }
            else {
              warning("Optimization not successful. Profile may be erroneous.")
            }
          }
          return(y.try)
        }
        doAdaption <- function() {
          lagrange.out.try <- lagrange(y.try)
          valid <- TRUE
          dobj.pred <- sum(lagrange.out$gradient * (y.try - 
                                                      y))
          dobj.fact <- lagrange.out.try$value - lagrange.out$value
          correction <- lagrange.out.try$correction
          if (correction > aControl$correction) 
            gamma <- gamma/2
          if (correction < 0.5 * aControl$correction) 
            gamma <- min(c(aControl$gamma, gamma * 2))
          if (abs(dobj.fact - dobj.pred) > sControl$atol & 
              stepsize > sControl$min) {
            stepsize <- max(c(stepsize/1.5, sControl$min))
            valid <- FALSE
          }
          if (abs(dobj.fact - dobj.pred) < 0.3 * sControl$atol | 
              abs((dobj.fact - dobj.pred)/dobj.fact) < 0.3 * 
              sControl$rtol) {
            stepsize <- min(c(stepsize * 2, sControl$max))
          }
          if (verbose) {
            diff.thres <- diff.steps <- diff.limit <- 0
            if (threshold < Inf) 
              diff.thres <- 1 - max(c(0, min(c(1, (threshold - 
                                                     lagrange.out.try$value)/delta))))
            if (sControl$limit < Inf) 
              diff.steps <- i/sControl$limit
            diff.limit <- switch(as.character(sign(constraint.out$value)), 
                                 `1` = 1 - (limits[2] - constraint.out$value)/limits[2], 
                                 `-1` = diff.limit <- 1 - (limits[1] - constraint.out$value)/limits[1], 
                                 `0` = 0)
            percentage <- max(c(diff.thres, diff.steps, 
                                diff.limit), na.rm = TRUE) * 100
            progressBar(percentage)
            myvalue <- format(substr(lagrange.out$value, 
                                     0, 8), width = 8)
            myconst <- format(substr(constraint.out$value, 
                                     0, 8), width = 8)
            mygamma <- format(substr(gamma, 0, 8), width = 8)
            myvalid <- all(lagrange.out$valid)
            cat("\tvalue:", myvalue, "constraint:", myconst, 
                "gamma:", mygamma, "valid:", myvalid)
          }
          return(list(lagrange = lagrange.out.try, stepsize = stepsize, 
                      gamma = gamma, valid = valid))
        }
        init <- function(){
          e <- parent.env(environment())
          if (verbose & direction == 1) {
            cat("Compute right profile\n")
          }else if (verbose & direction == -1) {
            cat("\nCompute left profile\n")
          }
          
          e$i <- 0
          e$y <- ini
          lagrange_ini <- lagrange(ini)
          constraint_ini <- constraint(pars)
          e$threshold <- lagrange_ini[[sControl$stop]] + delta
          if(direction == 1){
            e$out.attributes <- unlist(lagrange_ini[lagrange_ini$attributes])
            e$out <- c(value = lagrange_ini$value, constraint = as.vector(constraint_ini$value), 
                     stepsize = stepsize, gamma = gamma, whichPar = whichIndex, 
                     out.attributes, ini)
          }
          e$lagrange.out <- lagrange_ini
          e$constraint.out <- constraint_ini
        }
        continue <- function(){
          e <- parent.env(environment())
          while (i < sControl$limit) {
            e$sufficient <- FALSE
            e$retry <- 0
            
            while (!sufficient & retry < 5) {
              # print(stepsize)
              e$dy <- stepsize * lagrange.out$dy
              e$y.try <- try(doIteration(), silent = TRUE)
              e$out.try <- try(doAdaption(), silent = TRUE)
              if (inherits(y.try, "try-error") | inherits(out.try, 
                                                          "try-error")) {
                e$sufficient <- FALSE
                e$stepsize <- stepsize/1.5
                e$retry <- retry + 1
              }else {
                e$sufficient <- out.try$valid
                e$stepsize <- out.try$stepsize
             }
            }
            if (inherits(y.try, "try-error") | inherits(out.try, 
                                                        "try-error")) 
              break
            
            e$y <- y.try
            e$lagrange.out <- out.try$lagrange
            e$constraint.out <- constraint(y.try)
            e$stepsize <- out.try$stepsize
            e$gamma <- out.try$gamma
            e$out.attributes <- unlist(lagrange.out[lagrange.out$attributes])
            
            if (direction == 1){
              e$out <- rbind(out, c(value = lagrange.out$value, 
                                  constraint = as.vector(constraint.out$value), 
                                  stepsize = stepsize, gamma = gamma, whichPar = whichIndex, 
                                  out.attributes, y))
            }else if (direction == -1){
              e$out <- rbind(c(value = lagrange.out$value, constraint = as.vector(constraint.out$value), 
                             stepsize = stepsize, gamma = gamma, whichPar = whichIndex, 
                             out.attributes, y), out)
            }
            # print(c(value = lagrange.out$value, 
            #         constraint = as.vector(constraint.out$value), 
            #         stepsize = stepsize, gamma = gamma, whichPar = whichIndex, 
            #         out.attributes, y))
            e$value <- lagrange.out[[sControl$stop]]
            
            if(direction == 1){
              if (value > threshold | constraint.out$value > 
                  limits[2]) 
                break
            }else if (direction == -1){
              if (value > threshold | constraint.out$value < 
                  limits[1]) 
                break
            }
            e$i <- i + 1
            if(save_iteration){
              save(list = ls(all=T,envir = e), envir = e, file = file_name)
            }
          }
        }

        gamma <- aControl$gamma
        stepsize <- sControl$stepsize
        ini <- pars
        delta <- qchisq(1 - alpha, 1)
        
        for( direction in c(1,-1)){
          init()
          
          continue()
          
        }
        out <- as.data.frame(out)
        out$whichPar <- whichPar.name
        parframe(out, parameters = names(pars), metanames = c("value", 
                                                            "constraint", "stepsize", "gamma", "whichPar"), 
                obj.attributes = names(out.attributes))
      }
    }
  if (Sys.info()[["sysname"]] == "Windows" & cores > 1) {
    parallel::stopCluster(cluster)
    doParallel::stopImplicitCluster()
  }
  do.call(rbind, out)
}