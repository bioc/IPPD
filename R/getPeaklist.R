setGeneric("getPeaklist", function(mz, intensities, model = c("Gaussian", "EMG"), 
                                   model.parameters = list(alpha = function(){},
                                   sigma = function(){},
                                   mu = function(){}),
                                   loss = c("L2", "L1"), 
                                   binning = FALSE,
                                   postprocessing = TRUE,
                                   trace = TRUE,
                                   returnbasis = TRUE,
                                   control.basis = list(charges = c(1,2,3,4), eps = 1e-05),
                                   control.localnoise = list(quantile = 0.5, factor.place = 1.5, 
                                     factor.post = 0, window = NULL, subtract = FALSE),
                                   control.postprocessing = list(mzfilter = FALSE, prune = FALSE,
                                     factor.prune = NULL,   ppm = NULL, goodnessoffit = FALSE),         
                                   control.binning = list(tol = 0.01))
           standardGeneric("getPeaklist"))

setMethod("getPeaklist", signature(mz = "numeric", intensities = "numeric"),
          function(mz, intensities,
                   model = c("Gaussian", "EMG"),
                   model.parameters = list(alpha = function(){},
                     sigma = function(){},
                     mu = function(){}),
                   loss = c("L2", "L1"),
                   binning = FALSE,
                   postprocessing = TRUE,
                   trace = TRUE,
                   returnbasis = TRUE,
                   control.basis = list(charges = c(1,2,3,4), eps = 1e-05),
                   control.localnoise = list(quantile = 0.5, factor.place = 1.5, #factor.cutoff = 4,
                     factor.post = 0, window = NULL, subtract = FALSE),
                   control.postprocessing = list(mzfilter = FALSE, prune = FALSE,
                     factor.prune = NULL,   ppm = NULL, goodnessoffit = FALSE),#,  ###refit = TRUE),                 
                   control.binning = list(tol = 0.01)){
#############################################################################
x <- mz
if(any(is.na(x))) stop("'mz' contains missing values \n")
y <- intensities
if(any(is.na(y))) stop("'intensities' contains missing values \n")
n <- length(x)
if(length(y) != n) stop("Length of 'mz' and length of 'intensities' differ \n")
if(any(y < 0))
  stop("'y' must be nonnegative \n")
model <-  match.arg(model)
if(!is.element(model, c("Gaussian", "EMG")))
  stop("'model' must be one of 'Gaussian' or 'EMG' \n")
############################################################################

charges <- control.basis$charges
if(is.null(control.basis$charges))
   stop("'control.basis$charges' has to be specified \n")
 eps <- control.basis$eps
if(is.null(control.basis$eps))
   eps <- 1e-05
if(model == "Gaussian"){
    if(class(model.parameters) == "list"){
            sigma <- model.parameters$sigma
            if(!is.function(sigma)) stop("'model.parameters$sigma' must be a function \n")
            alpha <- function(mz){}
            mu <- function(mz){}
          }  
          else{
            if(class(model.parameters) == "modelfit"){
            sigma <- slot(model.parameters, "sigmafunction")
            alpha <- slot(model.parameters, "alphafunction")
            mu <- slot(model.parameters, "mufunction")
          }
            else
              stop("'model.parameters' specified in an invalid way \n")
          }
}

if(model == "EMG"){
 if(class(model.parameters) == "list"){
  alpha <- model.parameters$alpha
  sigma <- model.parameters$sigma
  mu <- model.parameters$mu
  if(any(!is.function(alpha), !is.function(sigma), !is.function(mu)))
     stop("'model.parameter$alpha',  'model.parameter$sigma',' modelparameter$mu' must be functions \n")
}
else{
  if(class(model.parameters) == "modelfit"){
  alpha <- slot(model.parameters, "alphafunction")
  sigma <- slot(model.parameters, "sigmafunction")
  mu <- slot(model.parameters, "mufunction")
  }
  else
   stop("'model.parameters' specified in an invalid way \n")
}
}  
#############################################################################
if(trace) cat("Computing local noise level ... \n")
         #if(is.null(control.localnoise$window)){
           chargemin <- min(charges)
           ### from 26/11/09: always computed
           if(model == "Gaussian"){
             temp <- calculatebasis.gaussian(x, positions = x[floor(n/2) + 1], sigma = sigma, charges = chargemin, eps = eps,
                                        uppernonzero = length(x))
             supp <- as.logical(temp$Phi > 0)
             windowtemp <- sum(supp)
             windowmztemp <- diff(x[range(which(supp))]) 
           }
           if(model == "EMG"){
             temp <- calculatebasis.emg(x, positions = x[floor(n/2) + 1], alpha = alpha, sigma = sigma, mu = mu, charges = chargemin, eps = eps,
                                        uppernonzero = length(x))
             supp <- as.logical(temp$Phi > 0)
             windowtemp <- sum(supp)
             windowmztemp <- diff(x[range(which(supp))])
           }
         #}
         #else{
           ## if window is pre-specified by the user
           if(!is.null(control.localnoise$window)){
             window.mz <- control.localnoise$window
             window <- ceiling(window.mz / windowmztemp * windowtemp)
           }
           else window <- windowtemp 
           # window.mz <- control.localnoise$window
           #spacing.median <- median(diff(x))
           #spacing.mz <- x[floor(n/2) + 1] - x[floor(n/2)]
           #if(spacing.mz > 1.5 * spacing.median) spacing.mz <- spacing.median
           #window <- ceiling(window.mz / spacing.mz)
         #}
         #}

         if(window <= 0 | window > n)
            stop("'window' must be between 0 and length('mz') \n")

  if((window %% 2) == 0) window <- window + 1

###############################################################################

if(is.null(control.localnoise$factor.place))
  factor.place <- 1.5
else
  factor.place <- control.localnoise$factor.place

if(is.null(control.localnoise$quantile))
  quant <- 0.5
else
  quant <- control.localnoise$quantile


                                        #if(is.null(control.localnoise$quantiledist))
  quantiledist <- 0.1
#else
#  quantiledist <- control.localnoise$quantiledist

quantseqq <- seq(from = 0, to = 1, by = quantiledist)
lquantseqq <- length(quantseqq)
quantindex<- min(findInterval(quant, quantseqq), lquantseqq)


#if(is.null(control.localnoise$factor.cutoff))
#  factor.cutoff <- 4
#else
#  factor.cutoff <- control.localnoise$factor.cutoff


locnoise  <- localnoise(y, window = window, quantiledist = quantiledist)
locnoisequant <- locnoise[,quantindex]
biggereps <- locnoisequant > eps
if(sum(biggereps) == 0){
#  medianouter <- 0
  stop("Local noise level is zero everywhere \n")
}
else medianouter <- median(locnoisequant[biggereps]) 
#medianinner <- median(locnoisequant[locnoisequant > medianouter * 0.25])
locnoise.cutoff <- pmax(locnoisequant, medianouter * 0.25) ### !lower bound !
#whichzerolocnoisequant <- which(locnoisequant < eps)
                                        #nozeros <- length(whichzerolocnoisequant)
#qplus <- quantindex
#while(nozeros > 0){
#  qplus <- qplus + 1
#  if(qplus > ncol(locnoise)) stop("Local noise level does not exceed zero everywhere. \n")
#  locnoisequant[whichzerolocnoisequant] <- locnoise[whichzerolocnoisequant,qplus]
#  whichzerolocnoisequant <- which(locnoisequant < eps)
  #nozeros <- length(whichzerolocnoisequant)
#}  
locnoise.basis <- locnoise.cutoff * factor.place

#locnoise.cutoffpost <- locnoise.cutoff * factor.cutoff

if(is.null(control.localnoise$factor.post))
  factor.post <- 0
else
  factor.post <- control.localnoise$factor.post
#browser()
##############################################################################
if(is.null(control.localnoise$subtract))
  subtract <- FALSE
else
  subtract <- control.localnoise$subtract
##############################################################################
if(!binning){
  if(trace) cat("Computing basis functions ... \n")
positions <- x[y > locnoise.basis & locnoise.basis > eps]
npositions <- length(positions)
if(npositions == 0)
             stop("No potential peak locations found \n")
        if(npositions <= floor(1/1000 * n))
          warning("Very few potential peak locations (<= 0.001 * sampling points) \n")
########################################################################################
deltamean <- mean(diff(x))
##############################################################################
        if(model == "Gaussian"){
           meansigma <- mean(sigma(x))
           testgr <- c(-(100:0) * deltamean, (1:100) * deltamean)
           proposalpeak <- sum(gaussfun(testgr, mu = 0, sigma = meansigma) > eps)
           
           uppernonzero <- npositions * (10 * length(charges)) * proposalpeak 

           
           bas <- try(calculatebasis.gaussian(x, positions = positions, sigma = sigma, charges = charges, eps = eps, uppernonzero = uppernonzero), silent = TRUE)
        
}

if(model == "EMG"){
              sigmax <- sigma(x)
              alphax <- alpha(x)
              mux <- mu(x)
              meansigma <- sort(sigmax)[ceiling(n/2)]
              meanalpha <- sort(alphax)[ceiling(n/2)]
              ###
              meansigma.x <- x[which(sigmax == meansigma)[1]]
              meanmu.sigma <- mu(meansigma.x)
              meanalpha.sigma <- alpha(meansigma.x)
              ###
              meanalpha.x <- x[which(alphax == meanalpha)[1]]
              meanmu.alpha <- mu(meanalpha.x)
              meansigma.alpha <- sigma(meanalpha.x)                
              ####
              testgr <- c(-(100:0) * deltamean, (1:100) * deltamean)
              ###
              proposalpeak.alpha <- sum(EMG(testgr, mu = meanmu.alpha, sigma = meansigma.alpha, alpha = meanalpha) > eps)
              proposalpeak.sigma <- sum(EMG(testgr, mu = meanmu.sigma, sigma = meansigma, alpha = meanalpha.sigma) > eps)
              ###
              proposalpeak <- mean(proposalpeak.alpha, proposalpeak.sigma)
         
           uppernonzero <- npositions * (10 * length(charges)) * proposalpeak
              bas <- try(calculatebasis.emg(x, positions = positions, alpha = alpha, sigma = sigma, mu = mu, charges = charges, eps = eps, uppernonzero = uppernonzero), silent = TRUE)
            }


 if(inherits(bas, "try-error"))
          stop(paste("Error when computing basis function matrix:", as.character(bas), sep = " "))

 basis <- bas$Phi
 book <- bas$book
 rm(bas); gc()

  G <- forceSymmetric(t(basis) %*% basis)

  
  loss <- match.arg(loss)
  if(!is.element(loss, c("L2", "L1")))
    stop("'loss' must be either 'L2' or 'L1' \n")

if(loss == "L2"){  
if(subtract) scale.y <- max(y - locnoise.cutoff)
  else scale.y <- max(y)

if(subtract) C <- drop(t(basis) %*% (y - locnoise.cutoff)/(scale.y))
  else C <- drop(t(basis) %*% y/scale.y) 


if(subtract) nnlssol <- try(nnlslogbarrier((y - locnoise.cutoff)/scale.y, betastart = rep(0.01, ncol(G)), trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-06), silent = TRUE)

else nnlssol <- try(nnlslogbarrier(y/scale.y, betastart = rep(0.01, ncol(G)), trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-06), silent = TRUE)

  
        if(inherits(nnlssol, "try-error"))
          stop("Error in non-negative least squares estimation:", as.character(nnlssol), sep = " ")

  beta <- (nnlssol$beta * scale.y)
}
else{
  if(subtract) scale.y <- max(y - locnoise.cutoff)
  else scale.y <- max(y)

  if(subtract) nnladsol <- nnladlogbarrier((y - locnoise.cutoff)/scale.y, betastart = rep(0.01, ncol(basis)), trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-06)

  else nnladsol <- try(nnladlogbarrier(y/scale.y, betastart = rep(0.01, ncol(basis)), trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-06), silent = TRUE)

  
        if(inherits(nnladsol, "try-error"))
          stop("Error in non-negative least absolute deviation estimation:", as.character(nnladsol), sep = " ")

  beta <- (nnladsol$beta * scale.y)


}
}

else{  
if(trace) cat("Generating binning... \n")
returnbasis <- FALSE
        tol <- control.binning$tol
        if(is.null(tol)) tol <- 0.01
        minx <- min(x)
maxx <- max(x)
knots <- seq(from = minx, to = maxx, by =(1+tol))
        if(!is.element(maxx, knots)) knots <- c(knots, maxx)
        intervals <- cbind(knots[-length(knots)], knots[-1])
        binsC <- .C("generatebinning", x = x, y = y, locnoise = locnoise.basis,
		     intervals = intervals, binstart = double(length(knots)),
                     binend = double(length(knots)),  
                     nintervals = as.integer(nrow(intervals)),
                     n = as.integer(length(x)),
                     flag = as.integer(0), nbins = as.integer(0))
        bins <- cbind(c(minx, binsC$binstart[2:(binsC$nbins+ 1)]), binsC$binend[1:(binsC$nbins+ 1)])
        lbins <- nrow(bins)
        book <- NULL
        beta <- NULL
        if(model == "Gaussian"){
           if(class(model.parameters) == "list"){
            sigma <- model.parameters$sigma
            if(!is.function(sigma)) stop("'model.parameters$sigma' must be a function \n")
            alpha <- function(mz){}
            mu <- function(mz){}
          }  
          else{
            if(class(model.parameters) == "modelfit"){
            sigma <- slot(model.parameters, "sigmafunction")
            alpha <- slot(model.parameters, "alphafunction")
            mu <- slot(model.parameters, "mufunction")
          }
            else
              stop("'model.parameters' specified in an invalid way \n")
          } 
        }  
        if(model == "EMG"){
           if(class(model.parameters) == "list"){
           alpha <- model.parameters$alpha
           sigma <- model.parameters$sigma
           mu <- model.parameters$mu
           if(any(!is.function(alpha), !is.function(sigma), !is.function(mu)))
             stop("'model.parameter$alpha',  'model.parameter$sigma','modelparameter$mu' must be functions \n")
         }
           else{
             if(class(model.parameters) == "modelfit"){
               alpha <- slot(model.parameters, "alphafunction")
               sigma <- slot(model.parameters, "sigmafunction")
               mu <- slot(model.parameters, "mufunction")
             }
             else
               stop("'model.parameters' specified in an invalid way \n")
           }
        }
        for(i in 1:lbins){
          if(trace) cat("bin:", i, "\n")
          f1 <- x >= bins[i,1] &  x < bins[i,2]
          if(sum(f1) == 0){
            next
          }
  
          positions <- x[f1][y[f1] > locnoise.basis[f1] & locnoise.basis[f1] > 0]
          npositions <- length(positions)
          
          if(trace) cat("number of positions:", npositions, "\n")

          if(npositions == 0) next

          deltamean <- mean(diff(x[f1]))
          
          if(model == "Gaussian"){
            meansigma <- mean(sigma(x[f1]))
            testgr <- c(-(100:0) * deltamean, (1:100) * deltamean)
            proposalpeak <- sum(gaussfun(testgr, mu = 0, sigma = meansigma) > eps)
           
            uppernonzero <- npositions * (10 * length(charges)) * proposalpeak 
            bas <-  try(calculatebasis.gaussian(x[f1], positions = positions, sigma = sigma, charges = charges, eps = eps,
                                 uppernonzero = uppernonzero), silent = TRUE)
          }
          if(model == "EMG"){
             sigmax <- sigma(x[f1])
              alphax <- alpha(x[f1])
              mux <- mu(x[f1])
              meansigma <- sort(sigmax)[ceiling(sum(f1)/2)]
              meanalpha <- sort(alphax)[ceiling(sum(f1)/2)]
              ###
              meansigma.x <- x[f1][which(sigmax == meansigma)[1]]
              meanmu.sigma <- mu(meansigma.x)
              meanalpha.sigma <- alpha(meansigma.x)
              ###
              meanalpha.x <- x[f1][which(alphax == meanalpha)[1]]
              meanmu.alpha <- mu(meanalpha.x)
              meansigma.alpha <- sigma(meanalpha.x)                
              ####
              testgr <- c(-(100:0) * deltamean, (1:100) * deltamean)
              ###
              proposalpeak.alpha <- sum(EMG(testgr, mu = meanmu.alpha, sigma = meansigma.alpha, alpha = meanalpha) > eps)
              proposalpeak.sigma <- sum(EMG(testgr, mu = meanmu.sigma, sigma = meansigma, alpha = meanalpha.sigma) > eps)
              ###
              proposalpeak <- mean(proposalpeak.alpha, proposalpeak.sigma)
         
             uppernonzero <- npositions * (10 * length(charges)) * proposalpeak
             bas <-  try(calculatebasis.emg(x[f1], positions = positions, alpha = alpha, sigma = sigma, mu = mu, charges = charges, eps = eps,
                                            uppernonzero = uppernonzero), silent = TRUE)

          }   
          if(inherits(bas, "try-error"))  stop(paste("Error when computing basis function matrix:", as.character(bas), sep = " "))
  
          basis <- bas$Phi
          bookbin <- bas$book
          G <- crossprod(basis)

          if(loss == "L2"){
          
          if(subtract) C <- drop(t(basis) %*% (y[f1] - locnoise.cutoff[f1])/max((y[f1] - locnoise.cutoff[f1])))
            else C <- drop(t(basis) %*% y[f1]/max(y[f1]))

          if(subtract) nnlssol <- try(nnlslogbarrier(response = (y[f1] - locnoise.cutoff[f1])/max((y[f1] - locnoise.cutoff[f1])), betastart = rep(0.01, ncol(G)), trace = trace,
                                        alpha =  0.01, gammastart = 10,
                                        gammamax = 10^15, gammamult = 20, eps = 1e-6)) 
            else nnlssol <- try(nnlslogbarrier(response = y[f1]/max(y[f1]), betastart = rep(0.01, ncol(G)), trace = trace,
                                        alpha =  0.01, gammastart = 10,
                                        gammamax = 10^15, gammamult = 20, eps = 1e-6))

          if(inherits(nnlssol, "try-error")){
           warning("Error in non-negative least squares estimation:", as.character(nnlssol), sep = " ")
           next
          }
          if(subtract) betabin <- nnlssol$beta * max(y[f1] - locnoise.cutoff[f1])
            else betabin <- nnlssol$beta * max(y[f1])
        }

        else{
           if(subtract) nnladsol <- try(nnladlogbarrier(response = (y[f1] - locnoise.cutoff[f1])/max((y[f1] - locnoise.cutoff[f1])), betastart = rep(0.01, ncol(G)), trace = trace,
                                                        alpha =  0.01, gammastart = 10,
                                        gammamax = 10^15, gammamult = 20, eps = 1e-6)) 
            else nnladsol <- try(nnladlogbarrier(response = y[f1]/max(y[f1]), betastart = rep(0.01, ncol(G)), trace = trace,
                                        alpha =  0.01, gammastart = 10,
                                        gammamax = 10^15, gammamult = 20, eps = 1e-6))

          if(inherits(nnladsol, "try-error")){
           warning("Error in non-negative least absolute deviation estimation:", as.character(nnlssol), sep = " ")
           next
          }
          if(subtract) betabin <- nnladsol$beta * max(y[f1] - locnoise.cutoff[f1])
            else betabin <- nnladsol$beta * max(y[f1])
        }  
          book <- rbind(book, bookbin)
          beta <- c(beta, betabin)
        }
}

 if(trace) cat("Generating peaklist ... \n")
        locnoiseindices <- findInterval(book[,2], x)
        cutoff <- locnoise.cutoff[locnoiseindices]
        peaklist <-   cbind(book, beta)
        orr <- order(peaklist[,1])
        peaklisto <- peaklist[orr, , drop = FALSE]
        colnames(peaklisto) <- c("loc_init", "loc_most_intense", "charge", "amplitude")
        if(postprocessing){
         ##################################################################################
            presel <- peaklisto[,"amplitude"] > factor.post * cutoff
            if(sum(presel) == 0){
              warning("Post-processing failed; there are no peaks entering postprocessing \n")
              peaklistprocessed <- matrix(0)
              gof <- numeric(0)
            }
            else{
              if(is.null(control.postprocessing$mzfilter))
                mzfilter <- FALSE
              else mzfilter <- control.postprocessing$mzfilter

              if(is.null(control.postprocessing$prune))
                prune <- FALSE
              else prune <- control.postprocessing$prune

              if(is.null(control.postprocessing$prune))
                prune <- FALSE
              else prune <- control.postprocessing$prune

              if(prune){
                factor.prune <- control.postprocessing$factor.prune
                if(is.null(factor.prune))
                  factor.prune <- 0.05
                if(factor.prune >= 1 | factor.prune <= 0)
                  stop("'factor.prune' is not between 0 and 1 \n")
             }
              ppm <- control.postprocessing$ppm     
      
              if(is.null(ppm)){
                deltainit <- (x[2] - x[1])
                ppm <- deltainit / mean(c(x[1:2])) * 10^6
              }

              goodnessoffit <- control.postprocessing$goodnessoffit    
      
              if(is.null(goodnessoffit))
                goodnessoffit <- FALSE
              else
                goodnessoffit <- control.postprocessing$goodnessoffit    
         ##################################################################################   
            if(model == "Gaussian"){ 
              peaklistprocessed <- try(postprocess.gaussian(ppm = ppm, sigma = sigma, peaklist = peaklisto[presel, , drop = FALSE], trace = trace), silent = TRUE)
            }
              if(model == "EMG"){
                delta.x <- max(diff(x)) ### to be corrected in the future
                peaklistprocessed <- try(postprocess.emg(ppm = ppm, delta.x = delta.x, alpha = alpha, sigma = sigma, mu = mu, peaklist = peaklisto[presel, , drop = FALSE], trace = trace, npoints = 100), silent = TRUE)
         }
              
        if(inherits(peaklistprocessed, "try-error")){
          warning("Post-processing failed; try to reduce 'ppm' \n")
          peaklistprocessed <- matrix(0)
          gof <- numeric(0) 
        }
            
        else{
          if(trace) cat("Postprocessing successful... \n")     
      
          if(mzfilter){
            if(trace) cat("Applying m/z filter... \n")
            filtered <- .C("mzfilter", positions = as.double(peaklistprocessed[,1]), npositions = as.integer(nrow(peaklistprocessed)),
                            charge = as.integer(peaklistprocessed[,3]),  filteredout = integer(nrow(peaklistprocessed)))
              peaklistprocessed <- peaklistprocessed[!as.logical(filtered$filteredout), ,drop = FALSE]
          }
           
          if(prune){
            if (trace) cat("Pruning... \n")
            if(is.null(factor.prune)) factor.prune <- 0.05
            indexprune <- findInterval(peaklistprocessed[,2], x)
            cutoff.prune <- locnoise[,"1"][indexprune]
            pruned <- peaklistprocessed[,4] < (cutoff.prune * factor.prune)
            peaklistprocessed <- peaklistprocessed[!pruned,,drop = FALSE]
          }
          #################################################################
          if(goodnessoffit){
          if(trace) cat("Evaluating goodness-of-fit... \n")
          if(!binning){
          basispattern <- basis
          rm(basis);gc()
          uppernonzero <- min(proposalpeak * n * mean(diff(mz))/min(diff(mz)), 10^7)
          if(model == "Gaussian"){
            basis <- try(calculatedenoisingbasis.gaussian(x, sigma, eps = 1e-05, uppernonzero = uppernonzero), silent = TRUE)
          }  
          if(model == "EMG"){
            basis <- try(calculatedenoisingbasis.emg(x, alpha, sigma, mu, eps = 1e-05, uppernonzero = uppernonzero), silent = TRUE)
          }
          #################################################################
          if(inherits(basis, "try-error"))
            stop(paste("Error when computing basis function matrix for goodness-of-fit:", as.character(basis), sep = " "))
          #################################################################
          G <- forceSymmetric(t(basis) %*% basis)
          C <- as.numeric(t(basis) %*% (y/scale.y))
          if(loss == "L2"){
          nnlsgof <- try(nnlslogbarrier(response = y / scale.y, ### Error in non-negative least absolute deviation estimation
                                    betastart = rep(0.01, n),
                                    trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-6), silent = TRUE)
          if(inherits(nnlsgof, "try-error")) 
            stop(paste("Error in non-negative least squares estimation for goodness-of-fit:", as.character(nnlsgof), sep = " ")) 
          residualsgof <- (y - drop(basis %*% nnlsgof$beta) * scale.y)^2
          mov.residualsgof <- as.numeric(filter(residualsgof, filter = rep(1/window, window)))                  
          whichna <- which(is.na(mov.residualsgof))
          compwhichna <- setdiff(1:n, whichna)
          mincompwhichna <- min(compwhichna)
          maxcompwhichna <- max(compwhichna)
          mov.residualsgof[whichna[whichna < mincompwhichna]] <- mov.residualsgof[mincompwhichna]
          mov.residualsgof[whichna[whichna > maxcompwhichna]] <- mov.residualsgof[maxcompwhichna]   
          mov.y <- as.numeric(filter(y^2, filter = rep(1/window, window)))
          mov.y[whichna[whichna < mincompwhichna]] <- mov.y[mincompwhichna]
          mov.y[whichna[whichna > maxcompwhichna]] <- mov.y[maxcompwhichna]  
          gof <- pmax(1 - mov.residualsgof / mov.y, 0.5)
 #         browser()
        }
          if(loss == "L1"){
            nnladgof <- nnladlogbarrier(y / scale.y, betastart = rep(0.01, n), trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-06)
             if(inherits(nnladgof, "try-error")) 
               stop(paste("Error in non-negative least absolute deviation estimation for goodness-of-fit:", as.character(nnladgof), sep = " ")) 
            residualsgof <- abs(y - drop(basis %*% nnladgof$beta) * scale.y)
            mov.residualsgof <- as.numeric(filter(residualsgof, filter = rep(1/window, window)))
            whichna <- which(is.na(mov.residualsgof))
            compwhichna <- setdiff(1:n, whichna)
            mincompwhichna <- min(compwhichna)
            maxcompwhichna <- max(compwhichna)  
            mov.residualsgof[whichna[whichna < mincompwhichna]] <- mov.residualsgof[mincompwhichna]
            mov.residualsgof[whichna[whichna > maxcompwhichna]] <- mov.residualsgof[maxcompwhichna]   
            mov.y <- as.numeric(filter(abs(y), filter = rep(1/window, window)))
            mov.y[whichna[whichna < mincompwhichna]] <- mov.y[mincompwhichna]
            mov.y[whichna[whichna > maxcompwhichna]] <- mov.y[maxcompwhichna]  
            gof <- pmax(1 - mov.residualsgof / mov.y, 0.5)
          }
        }
        else{### compute goodness-of-fit separately for each bin
           gof <- rep(0, length(x))
           for(i in 1:lbins){
             if(trace) cat("bin:", i, "\n")
             f1 <- x >= bins[i,1] &  x < bins[i,2]
             if(sum(f1) == 0){
               next
             }
             if(trace) cat("number of positions:", sum(f1), "\n")
             deltamean <- mean(diff(x[f1]))

             basispattern <- basis
             
             if(model == "Gaussian"){
               meansigma <- mean(sigma(x[f1]))
               testgr <- c(-(100:0) * deltamean, (1:100) * deltamean)
               proposalpeak <- sum(gaussfun(testgr, mu = 0, sigma = meansigma) > eps)
               uppernonzero <- sum(f1) * proposalpeak 
               basis <-  try(calculatedenoisingbasis.gaussian(x[f1], sigma = sigma, eps = eps,
                                 uppernonzero = uppernonzero), silent = TRUE) 
             }
             if(model == "EMG"){
               sigmax <- sigma(x[f1])
               alphax <- alpha(x[f1])
               mux <- mu(x[f1])
               meansigma <- sort(sigmax)[ceiling(sum(f1)/2)]
               meanalpha <- sort(alphax)[ceiling(sum(f1)/2)]
               ###
               meansigma.x <- x[f1][which(sigmax == meansigma)[1]]
               meanmu.sigma <- mu(meansigma.x)
               meanalpha.sigma <- alpha(meansigma.x)
               ###
               meanalpha.x <- x[f1][which(alphax == meanalpha)[1]]
               meanmu.alpha <- mu(meanalpha.x)
               meansigma.alpha <- sigma(meanalpha.x)                
               ####
               testgr <- c(-(100:0) * deltamean, (1:100) * deltamean)
               ###
               proposalpeak.alpha <- sum(EMG(testgr, mu = meanmu.alpha, sigma = meansigma.alpha, alpha = meanalpha) > eps)
               proposalpeak.sigma <- sum(EMG(testgr, mu = meanmu.sigma, sigma = meansigma, alpha = meanalpha.sigma) > eps)
               ###
               proposalpeak <- mean(proposalpeak.alpha, proposalpeak.sigma)
               uppernonzero <- sum(f1) * proposalpeak
               basis <-  try(calculatedenoisingbasis.emg(x[f1], alpha = alpha, sigma = sigma, mu = mu, eps = eps,
                                                       uppernonzero = uppernonzero), silent = TRUE)
             }   
          if(inherits(basis, "try-error"))  stop(paste("Error when computing basis function matrix for goodness-of-fit:", as.character(basis), sep = " "))
             G <- forceSymmetric(t(basis) %*% basis)
             scale.y <- max(y[f1])
             C <- as.numeric(t(basis) %*% (y[f1]/scale.y))
             if(loss == "L2"){
          nnlsgof <- try(nnlslogbarrier(response = y[f1] / scale.y, ### Error in non-negative least absolute deviation estimation
                                        betastart = rep(0.01, sum(f1)),
                                        trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-6), silent = TRUE)
          if(inherits(nnlsgof, "try-error")) 
            stop(paste("Error in non-negative least squares estimation for goodness-of-fit:", as.character(nnlsgof), sep = " ")) 
          residualsgof <- (y[f1] - drop(basis %*% nnlsgof$beta) * scale.y)^2
          mov.residualsgof <- as.numeric(filter(residualsgof, filter = rep(1/min(window, sum(f1)), min(sum(f1), window))))                  
          whichna <- which(is.na(mov.residualsgof))
          compwhichna <- setdiff(1:sum(f1), whichna)
          mincompwhichna <- min(compwhichna)
          maxcompwhichna <- max(compwhichna)
          mov.residualsgof[whichna[whichna < mincompwhichna]] <- mov.residualsgof[mincompwhichna]
          mov.residualsgof[whichna[whichna > maxcompwhichna]] <- mov.residualsgof[maxcompwhichna]   
          mov.y <- as.numeric(filter(y^2, filter = rep(1/min(window, sum(f1)), min(sum(f1), window))))
          mov.y[whichna[whichna < mincompwhichna]] <- mov.y[mincompwhichna]
          mov.y[whichna[whichna > maxcompwhichna]] <- mov.y[maxcompwhichna]  
          gof[f1] <- pmax(1 - mov.residualsgof / mov.y, 0.5)
        }
          if(loss == "L1"){
            nnladgof <- nnladlogbarrier(y[f1] / scale.y, betastart = rep(0.01, sum(f1)), trace = trace, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-06)
             if(inherits(nnladgof, "try-error")) 
               stop(paste("Error in non-negative least absolute deviation estimation for goodness-of-fit:", as.character(nnladgof), sep = " ")) 
            residualsgof <- abs(y - drop(basis %*% nnladgof$beta) * scale.y)
            mov.residualsgof <- as.numeric(filter(residualsgof, filter = rep(1/min(window, sum(f1)), min(window, sum(f1)))))
            whichna <- which(is.na(mov.residualsgof))
            compwhichna <- setdiff(1:sum(f1), whichna)
            mincompwhichna <- min(compwhichna)
            maxcompwhichna <- max(compwhichna)  
            mov.residualsgof[whichna[whichna < mincompwhichna]] <- mov.residualsgof[mincompwhichna]
            mov.residualsgof[whichna[whichna > maxcompwhichna]] <- mov.residualsgof[maxcompwhichna]   
            mov.y <- as.numeric(filter(abs(y), filter = rep(1/min(window, sum(f1)), min(window, sum(f1)))))
            mov.y[whichna[whichna < mincompwhichna]] <- mov.y[mincompwhichna]
            mov.y[whichna[whichna > maxcompwhichna]] <- mov.y[maxcompwhichna]  
            gof[f1] <- pmax(1 - mov.residualsgof / mov.y, 0.5)
          }   
           }
         }
          basis <- basispattern
        }
        else gof <- numeric(0)  
          #################################################################
          locnoiseindices   <- findInterval(peaklistprocessed[,2], x)
          peaklistprocessed <- cbind(peaklistprocessed, locnoise.cutoff[locnoiseindices])
          peaklistprocessed <- cbind(peaklistprocessed, peaklistprocessed[,4] / peaklistprocessed[,5])
          colnames(peaklistprocessed) <- c("loc_init", "loc_most_intense", "charge", "amplitude", "localnoise", "ratio") 
          #################################################################
          if(goodnessoffit){
          peaklistprocessed <- cbind(peaklistprocessed, goodness_of_fit = gof[locnoiseindices],
                                     ratio_adj = peaklistprocessed[,6] * gof[locnoiseindices])                               
         }
        }
            }
          }
          
else{
     peaklistprocessed <- matrix(0)
     gof <- numeric(0)
   }         
            
if(returnbasis)
 new("peaklist",  peaklist = peaklisto, peaklistprocessed = peaklistprocessed,
     model = model, loss = loss,
     alpha = alpha, sigma = sigma,
     mu = mu, charges = charges, basis = basis,          
                                                    book = book,
                                                    beta = beta,
                                                    locnoise = locnoise,
                                                    noiselevel = locnoise.cutoff,
                                                    goodnessoffit = gof,
                                                    data = list(x = x, y = y))
else
  new("peaklist",  peaklist = peaklisto, peaklistprocessed = peaklistprocessed,
     model = model, loss = loss,
     alpha = alpha, sigma = sigma,
     mu = mu, charges = charges, basis = sparseMatrix(i = 1, j = 1, x = 0),          
                                                    book = book,
                                                    beta = beta,
                                                    locnoise = locnoise,
                                                    noiselevel = locnoise.cutoff,
                                                    goodnessoffit = gof,
                                                    data = list(x = x, y = y))
  
     
})
