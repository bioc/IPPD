setGeneric("fitModelParameters", function(mz, intensities, model = c("Gaussian", "EMG"), 
fitting = c("most_intense", "model"), formula.alpha = formula(~1), formula.sigma = formula(~1), 
formula.mu = formula(~1), control = list(window = 6, threshold = NULL))
           standardGeneric("fitModelParameters"))

setMethod("fitModelParameters", signature(mz = "numeric", intensities = "numeric"), function(mz,
                                                            intensities,
             model = c("Gaussian", "EMG"),
                                                            fitting = c("most_intense", "model"),
             formula.alpha = formula(~1),
             formula.sigma = formula(~1),
                                                            formula.mu = formula(~1),
                                                            control = list(window = 6, threshold = NULL)){
  ################################################################
  x <- mz
  if(any(is.na(x))) stop("'mz' contains missing values \n")
  y <- intensities
  if(any(is.na(y))) stop("'intensities' contains missing values \n")
  n <- length(x)
  if(length(y) != n) stop("Length of 'mz' and length of 'intensities differ \n")
  if(any(y < 0))
    stop("'y' must be nonnegative \n")
model <-  match.arg(model)
if(!is.element(model, c("Gaussian", "EMG")))
  stop("'model' must be one of 'Gaussian' or 'EMG' \n")
  ################################################################   
  window <- control$window
  if(is.null(window)) window <- 6
  if(window <= 0 | !(as.integer(window) == window))
    stop("Control parameter 'window' has to be positive \n")
  threshold <- control$threshold
  if(is.null(threshold)){
    warning("'control$threshold' not specifed; set to 'max(intensities) - 1e-05' \n") 
    threshold <- max(y) - 1e-05
  }
   if(threshold < 0)
     stop("Control parameter 'threshold' has to be nonnegative \n")
  #################################################################
  fitting <- match.arg(fitting)
  if(!is.element(fitting, c("most_intense", "model")))
    stop("'fitting' must be one of 'most_intense' or 'model' \n")
  #################################################################
  if(fitting == "most_intense"){
    if(any(formula.alpha != formula(~1), formula.sigma != formula(~1), formula.mu != formula(~1)))
      warning("'fitting = 'most intense'', but non-default values for one of the formulae used. In the case a model should be fitted, set 'fitting = 'model'' \n")
    detection <-  simplepeakdetect(cbind(x,y), window = window, threshold = threshold)
    if(nrow(detection) < 2*window)
      stop("No peak of the chosen width ('window') found. Try to reduce 'window' \n")
    else{
      if(model == "Gaussian"){
        fitt <- try(fit.gauss(detection[,1], detection[,2]), silent = TRUE)
        if(inherits(fitt, "try-error"))
        stop("Fitting failed. \n")
        sigma <- fitt$sigma
        sigmafunction <- function(mz) {}
        body(sigmafunction) <- eval(substitute(expression(rep(sigmavar, length(mz))), list(sigmavar = sigma)))

        bestpeak <- list(mz = detection[,1], intensities  = detection[,2],  sigma = fitt$sigma, mu = fitt$mu)
        peakfitresults <- matrix(nrow = 1, ncol = 4, data = c(nrow(detection), fitt$rss, fitt$sigma, fitt$mu), byrow = TRUE)
        colnames(peakfitresults) <- c("datapoints", "rss", "sigma", "mz")
       }
      if(model == "EMG"){
        fitt <- try(fit.EMG(detection[,1], detection[,2], gridsearch = TRUE), silent = TRUE)
        if(inherits(fitt, "try-error"))
          stop("Fitting failed. \n")
        alpha <- fitt$alpha
        sigma <- fitt$sigma
        mu <- fitt$mu
        alphafunction <- function(mz){}
        sigmafunction <- function(mz){}
        mufunction <- function(mz){}
        body(alphafunction) <- eval(substitute(expression(rep(alphavar, length(mz))), list(alphavar = alpha)))
        body(sigmafunction) <- eval(substitute(expression(rep(sigmavar, length(mz))), list(sigmavar = sigma)))
        body(mufunction) <- eval(substitute(expression(rep(muvar, length(mz))),
                                               list(muvar = mu)))
        
      

      bestpeak <- list(mz = detection[,1], intensities  = detection[,2],  alpha = fitt$alpha, sigma = fitt$sigma, mu = fitt$mu)
      peakfitresults <- matrix(nrow = 1, ncol = 6, data = c(nrow(detection), fitt$rss, fitt$alpha, fitt$sigma, fitt$mu, mean(detection[,1])), byrow = TRUE)
      colnames(peakfitresults) <- c("datapoints", "rss", "alpha", "sigma", "mu", "mz")
    }
  }
  }
  if(fitting == "model"){ 
    require(MASS)
    detection <-  peakdetect(cbind(x,y), window = window, threshold = threshold)
    detection <- detection
    if(model == "Gaussian"){
       varnames.sigma <- all.vars(formula.sigma)
       if(length(unique(varnames.sigma)) > 1)
        stop("'formula.sigma' is invalid: only one variable is allowed \n'")
       if(length(unique(varnames.sigma)) > 0 && varnames.sigma != "mz")
          stop("'formula.sigma' is invalid: one or several
variables not equal to 'mz' are present \n")
       charsigma <- strsplit(as.character(formula.sigma), split = "")
        if(any(is.element(c(".", ":", "*"), charsigma)))
       stop("Invalid characters in one of the formulae: interaction terms/multiplications are not allowed \n")
       intercept.sigma <- attr(terms(formula.sigma), "intercept")

       peakfitresults <- matrix(nrow = 0, ncol = 4)
       colnames(peakfitresults) <- c("datapoints", "rss", "sigma", "mz")
    
    for(i in seq(along = detection)){
        x.i <- detection[[i]][,1]
        if(length(x.i) < 2*window) next
        y.i <- detection[[i]][,2]
         fitt <- try(fit.gauss(x.i, y.i), silent = TRUE)
        if(inherits(fitt, "try-error")){
          warning("Fitting failed. \n")
          next
        }
        newcol <- c(length(x.i), fitt$rss, fitt$sigma, fitt$mu)
        peakfitresults <- rbind(peakfitresults, newcol)
  }
  if(nrow(peakfitresults) == 0)
     stop("No peak of the chosen width ('window') found. Try to reduce 'window' \n")

       mz <- peakfitresults[,"mz"]
      sigmavec <- peakfitresults[,"sigma"]
       formula.sigmanew <- as.formula(paste("sigmavec", as.character(formula.sigma[2]), sep = "~")) 
     l1fitsigma <- try(rlm(formula.sigmanew, k = 1e-6), silent = TRUE)
     if(inherits(l1fitsigma, "try-error"))
     stop("Error in linear model estimation for parameter 'sigma' \n")
     coefsigma <- coef(l1fitsigma)
       
     sigmafunction <- formulacoef2function(formula.sigma, coef = coefsigma, intercept = intercept.sigma)

       MSE <- peakfitresults[,"rss"]/peakfitresults[,"datapoints"]
       bestpeakind <- which.min(MSE)
       bestpeak <- list(mz = detection[[bestpeakind]][,1], intensities  = detection[[bestpeakind]][,2],  sigma = peakfitresults[bestpeakind,"sigma"],
                        mu = peakfitresults[bestpeakind,"mz"])
     }
    if(model == "EMG"){
     varnames.alpha <- all.vars(formula.alpha)
     varnames.sigma <- all.vars(formula.sigma)
     varnames.mu <- all.vars(formula.mu)
     varnames <- c(varnames.alpha, varnames.sigma, varnames.mu)
     if(length(unique(varnames)) > 1)
        stop('One or several of the formulae are invalid: only one variable is allowed \n')
        if(length(unique(varnames)) > 0 && varnames != "mz")
          stop("One or several of the formulae are invalid: one or several
variables not equal to 'mz' are present \n")
     charalpha <- strsplit(as.character(formula.alpha), split = "")
     charsigma <- strsplit(as.character(formula.sigma), split = "")
     charmu <- strsplit(as.character(formula.mu), split = "")
     if(any(is.element(c(".", ":", "*"), c(charalpha, charsigma, charmu))))
       stop("Invalid characters in one of the formulae: interaction terms/multiplications are not allowed \n")
     intercept.alpha <- attr(terms(formula.alpha), "intercept")
     intercept.sigma <- attr(terms(formula.sigma), "intercept")
     intercept.mu <- attr(terms(formula.mu), "intercept")   
     
    peakfitresults <- matrix(nrow = 0, ncol = 6)
colnames(peakfitresults) <- c("datapoints", "rss", "alpha", "sigma", "mu", "mz")

     grid.alpha.basis <- grid.alpha <- 10^((seq(from = -5, to = 5, length = 100)))
     grid.sigma.basis <- grid.sigma <- 10^((seq(from = -5, to = 5, length = 100)))
     grid.mu <- seq(from = -1, to = 1, length = 100)
     
    for(i in seq(along = detection)){
        x.i <- detection[[i]][,1]
        if(length(x.i) < 2*window) next
        y.i <- detection[[i]][,2]
         fitt <- try(fit.EMG(x.i, y.i, gridsearch = TRUE, grid.alpha = grid.alpha, grid.sigma = grid.sigma, grid.mu = grid.mu), silent = TRUE)
        if(inherits(fitt, "try-error")){
          warning("Fitting failed. \n")
          next
        }
        newcol <- c(length(x.i), fitt$rss, fitt$alpha, fitt$sigma, fitt$mu, mean(x.i))
        peakfitresults <- rbind(peakfitresults, newcol)

        dist.alpha <- abs( fitt$alpha - grid.alpha.basis )
        dist.sigma <- abs( fitt$sigma - grid.sigma.basis )
        o.alpha <- order(dist.alpha)[1:20]
        o.sigma <- order(dist.sigma)[1:20]
        grid.alpha <- sort(grid.alpha.basis[o.alpha])
        grid.sigma <- sort(grid.sigma.basis[o.sigma])  
  }
  if(nrow(peakfitresults) == 0)
     stop("No peak of the chosen width ('window') found. Try to reduce 'window' \n")

     mz <- peakfitresults[,"mz"]

     alphavec <- peakfitresults[,"alpha"]
     formula.alphanew <- as.formula(paste("alphavec", as.character(formula.alpha[2]), sep = "~")) 
     l1fitalpha <- try(rlm(formula.alphanew, k = 1e-6), silent = TRUE)
     if(inherits(l1fitalpha, "try-error"))
       stop("Error in linear model estimation for parameter 'alpha' \n")
     coefalpha <- coef(l1fitalpha)
     alphafunction <- formulacoef2function(formula.alpha, coef = coefalpha, intercept = intercept.alpha)

     sigmavec <- peakfitresults[,"sigma"]
     formula.sigmanew <- as.formula(paste("sigmavec", as.character(formula.sigma[2]), sep = "~")) 
     l1fitsigma <- try(rlm(formula.sigmanew, k = 1e-6), silent = TRUE)
     if(inherits(l1fitsigma, "try-error"))
       stop("Error in linear model estimation for parameter 'sigma' \n")
     coefsigma <- coef(l1fitsigma)
     sigmafunction <- formulacoef2function(formula.sigma, coef = coefsigma, intercept = intercept.sigma)

     muvec <- peakfitresults[,"mu"]
     formula.munew <- as.formula(paste("muvec", as.character(formula.mu[2]), sep = "~")) 
     l1fitmu <- try(rlm(formula.munew, k = 1e-6), silent = TRUE)
     if(inherits(l1fitmu, "try-error"))
       stop("Error in linear model estimation for parameter 'mu' \n")
     coefmu <- coef(l1fitmu)
     mufunction <- formulacoef2function(formula.mu, coef = coefmu, intercept = intercept.mu)

     MSE <- peakfitresults[,"rss"]/peakfitresults[,"datapoints"]
     bestpeakind <- which.min(MSE)
     bestpeak <- list(mz = detection[[bestpeakind]][,1], intensities  = detection[[bestpeakind]][,2],  alpha = peakfitresults[bestpeakind,"alpha"], sigma = peakfitresults[bestpeakind,"sigma"], mu = peakfitresults[bestpeakind,"mu"])
   }
  }

    if(model == "Gaussian") alphafunction <- mufunction <- function(mz){}  
    new("modelfit", model = model, fitting = fitting, alphafunction = alphafunction,
        sigmafunction = sigmafunction, mufunction = mufunction, peakfitresults = peakfitresults,
        bestpeak = bestpeak)

   
 })
