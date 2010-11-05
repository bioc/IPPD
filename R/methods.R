################################################################################
setMethod("show", signature="peaklist", function(object){
     processed <- any(dim(object@peaklistprocessed) != dim(matrix(0)))
     text1 <- ifelse(processed,
                     "An object of class 'peaklist'(with postprocessing) \n",
                     "An object of class 'peaklist' (without postprocessing) \n")
     text1a <- paste(paste("Peak model used: ", object@model, sep=""), "\n")
     textnew <- paste(paste("Loss function used: ", object@loss, sep=""), "\n")
     nopeaks <- ifelse(processed, nrow(object@peaklistprocessed),
                                       nrow(object@peaklist))
     text3 <- paste(paste("number of peaks: ", nopeaks, sep = ""), "\n")
     text4 <- paste("charge states used:", paste(object@charges, collapse =","),
                    sep = " ")
     cat(text1); cat(textnew); cat(text1a); cat(text3); cat(text4); cat("\n")
          
   })
###############################################################################

setGeneric("threshold", function(object, threshold,  ratio = c("ratio", "ratio_adj"),
                                 refit = FALSE, trace = TRUE, eps = 1e-05, ...)
           standardGeneric("threshold"))  ### method has still to be documented

setMethod("threshold", signature="peaklist", function(object, threshold,
                         ratio = c("ratio", "ratio_adj"), refit = FALSE,
                         trace = TRUE, eps = 1e-05, ...){
     processed <- any(dim(object@peaklistprocessed) != dim(matrix(0)))
     noiselevel <- object@noiselevel
     x <- object@data$x
     y <- object@data$y
     model <- object@model
     alpha <- object@alpha
     sigma <- object@sigma
     mu <- object@mu
     ratio <- match.arg(ratio)
     if(!is.element(ratio, c("ratio", "ratio_adj")))
       stop("'ratio' must be one of 'ratio' or 'ratio_adj' \n")
     if(length(object@goodnessoffit) == 0 & ratio == "ratio_adj"){
       ratio <- "ratio"
       warning("Peaklist lacks goodness-of-fit statistic; 'ratio_adj' changed to 'ratio'. \n")
     }
     if(!processed){
       locnoise <- noiselevel[findInterval(object@peaklist[,"loc_most_intense"], x)]
       ratio <- object@peaklist[,"amplitude"] / locnoise
       sel <- object@peaklist[,"amplitude"] >= threshold
       }
     else{
        if(ratio == "ratio") sel <- object@peaklistprocessed[,"ratio"] >= threshold
        if(ratio == "ratio_adj") sel <- object@peaklistprocessed[,"ratio_adj"] >= threshold
      }
     if(sum(sel) == 0){
         warning("There are no peaks exceeding the chosen threshold \n")
         return(numeric())             
       }
     if(!processed){
       peaklisttrunc <- cbind(peaklist[sel,], locnoise[sel], ratio[sel])
       colnames(peaklisttrunc) <- c(colnames(peaklisttrunc), "localnoise", "ratio")
     }                          
     else  peaklisttrunc <- object@peaklistprocessed[sel, , drop = FALSE]
     if(refit){
       if(trace) cat("Re-fitting... \n")
         positionspro <- peaklisttrunc[,"loc_most_intense"]
         chargespro <- peaklisttrunc[,"charge"]
         uniqcharge <- unique(chargespro)
         perm <- c() 
         basis <- NULL
          for(l in seq(along = uniqcharge)){
            ll <- which(chargespro == uniqcharge[l])
            if(length(ll) == 0) next
            perm <- c(perm, ll)
            positionsprol <- positionspro[ll]
            if(model == "Gaussian"){
              if(is.null(basis))
                basis <- calculatebasis.gaussian(x, positions = positionsprol,
                                                 sigma = sigma, charges = uniqcharge[l], eps = eps,
                                                 uppernonzero = length(x) * length(positionsprol))$Phi
              else basis <- cbind2(basis, calculatebasis.gaussian(x, positions = positionsprol,
                                                                  sigma = sigma, charges = uniqcharge[l], eps = eps,
                                                                  uppernonzero = length(x) * length(positionsprol))$Phi) 
            }
            if(model == "EMG"){
              if(is.null(basis)) basis <- calculatebasis.emg(x, positions = positionsprol, alpha = alpha, sigma = sigma, mu = mu,
                                                             charges = uniqcharge[l], eps = eps,
                                                             uppernonzero = length(x) * length(positionsprol))$Phi
              else basis <- cbind2(basis, calculatebasis.emg(x, positions = positionsprol, alpha = alpha, sigma = sigma, mu = mu,
                                                             charges = uniqcharge[l], eps = eps,
                                                             uppernonzero = length(x) * length(positionsprol))$Phi) 
            }
          }
          if(is.null(basis)) stop("Error in re-fitting: no basis function exceeds cutoff \n") 
            G <- crossprod(basis)
            C <- drop(t(basis) %*% y/max(y))
            nlssol <- nnlslogbarrier(response = y/max(y), betastart = rep(mean(noiselevel), ncol(G)), trace = trace,
                                        alpha =  0.01, gammastart = 10,
                                        gammamax = 10^15, gammamult = 20, eps = 1e-6)
            peaklisttrunc[perm,"amplitude"] <- (nlssol$beta * max(y))
            peaklisttrunc[,"ratio"] <- peaklisttrunc[,"amplitude"] / peaklisttrunc[,"localnoise"]
            if(processed & length(object@goodnessoffit) > 0) peaklisttrunc[,"ratio_adj"] <- peaklisttrunc[,"ratio"] * peaklisttrunc[,"goodness_of_fit"]
            if(!processed) removeafterrefit <- peaklisttrunc[,"ratio"] < threshold
            if(processed){
              if(ratio == "ratio")
                removeafterrefit <- peaklisttrunc[,"ratio"] < threshold
              if(ratio == "ratio_adj")
                removeafterrefit <- peaklisttrunc[,"ratio_adj"] < threshold
              }
            peaklisttrunc <- peaklisttrunc[!removeafterrefit,,drop = FALSE]
            rownames(peaklisttrunc) <- NULL
     }
     return(peaklisttrunc)
   }
)


###############################################################################

setMethod("show", signature = "modelfit", function(object){
  fitting <- object@fitting
  if(fitting == "most_intense"){

    text1 <- paste(
          paste("Peak model '", paste(object@model, "'", sep=""), sep=""),
          "fitted using the most intense peak \n")

    cat(text1)

  }
  else{

     text1 <- paste(
          paste("Peak model '", paste(object@model, "'", sep=""), sep=""),
          "fitted as fuction of m/z \n")
     text2 <- paste(paste("number of peaks used:", nrow(object@peakfitresults), sep = " "), "\n", sep="")
     cat(text1); cat(text2)
  }
})

################################################################################


setGeneric("visualize", function(object, mz, intensities, ...) standardGeneric("visualize"))

setMethod("visualize", signature("peaklist", "numeric", "numeric"),
          function(object, mz, intensities, lower, upper,
                   truth = FALSE,
                   signal = TRUE,
                   fitted = TRUE,
                   postprocessed = TRUE,
                   fittedfunction = FALSE,
                   fittedfunction.cut = FALSE, ### to be modified 
                   quantile = NULL,
                   booktrue = NULL,
                   cutoff.eps = NULL,
                   cutoff.functions = 10,
                   ...){
            
model <- object@model
basisavailable <- is.character(all.equal(object@basis@x, 0))

if(model == "Gaussian")
  sigma <- object@sigma
if(model == "EMG"){
  alpha <- object@alpha
  sigma <- object@sigma
  mu <- object@mu
}

charges <- object@charges
ncharges <- length(charges)

x <- mz
ff <-  x >= lower & x <= upper
if(sum(ff) == 0)
  stop("'lower'/'upper' does not match 'mz'. \n")

if(fitted | fittedfunction | fittedfunction.cut){
  if(!basisavailable){
      section <- object@book[,2] >= lower & object@book[,2] <= upper
      pos <- object@book[section,2]
      opos <- order(object@book[section,3], pos)
      pos <- pos[opos]
      if(model == "Gaussian"){
        basis <- calculatebasis.gaussian(x[ff], positions = pos, sigma = sigma, charges = charges, eps = 1e-05, uppernonzero = length(pos) * ncharges * sum(ff))
      }
      if(model == "EMG"){
        basis <- calculatebasis.emg(x[ff], positions = pos, alpha = alpha, sigma = sigma, mu = mu, charges = charges, eps = 1e-05, uppernonzero = length(pos) * ncharges * sum(ff))
      }
      #sectionbeta
      if(nrow(basis$book) != length(object@beta[section]))
         stop("Matrix of basisfunctions needed for visualization cannot be computd \n")
      bookfitted <- cbind(basis$book, object@beta[section][opos])
      #browser()
  }
  else  bookfitted <- cbind(object@book, object@beta)
}


maxplots <- truth + signal + fitted + postprocessed
if(maxplots == 0) stop("Nothing has been selected for plotting \n")

layout(matrix(nrow = maxplots, ncol = 1, data = 1:maxplots))

if(!is.null(quantile)){
  quant <- as.character(quantile)
  if(length(x) != nrow(object@locnoise)) stop("Invalid 'mz' specified. \n")
 if(any(!is.element(quant, colnames(object@locnoise))))
 stop("Invalid 'quantile' chosen. Check the column names of 'slot(object, locnoise)' \n")
 localnoiselevel <- object@locnoise[findInterval(lower, x):findInterval(upper, x), quant]
 if(!is.matrix(localnoiselevel)) localnoiselevel <- matrix(localnoiselevel, nrow = length(localnoiselevel), ncol = length(quant))
 auxiliarysequence <- seq(from = min(x[ff]), max(x[ff]), length = nrow(localnoiselevel))
 
}

if(truth){
if(is.null(booktrue)) stop("If 'truth = TRUE', 'booktrue' has to be specified \n")
if(!is.matrix(booktrue))
  stop("Wrong format for 'booktrue' \n")
else
  if(ncol(booktrue) != 4)
    stop("Wrong format for 'booktrue' \n")
otrue <- order(booktrue[,2])
booktrue <- booktrue[otrue,,drop=FALSE]
booktrueind <- which(booktrue[,2] >= lower & booktrue[,2] <= upper)
howmanyfunctions <- length(booktrueind)
positionstrue <- booktrue[booktrueind,2]
chargetrue <- booktrue[booktrueind,3]
ampltrue <- booktrue[booktrueind,4]
if(model == "Gaussian")
  basistrue <- calculatebasis.gaussian(x[ff], positions = positionstrue, sigma = sigma, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff))
if(model == "EMG")
  basistrue <- calculatebasis.emg(x[ff], positions = positionstrue, alpha = alpha, sigma = sigma, mu = mu, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff))
for(i in seq(along = positionstrue)){

  colindex <- apply(basistrue$book[,2:3], 1, function(z) all(z == c(positionstrue[i], chargetrue[i])))
    
  if(i == 1){
      plot(x[ff], ampltrue[i] * basistrue$Phi[,colindex], col = i, main = "truth", type = "l", ylim = c(0, max(intensities[ff])), ...)
 }
 else lines(x[ff], ampltrue[i] * basistrue$Phi[,colindex], col = i, ...)
  }
if(!is.null(quantile)) matlines(auxiliarysequence, localnoiselevel, lty = "dotted", col = "black") 
}

####
if(signal){
  plot(x[ff], intensities[ff], main = "signal", ylim = c(0, max(intensities[ff])))
  if(fittedfunction | fittedfunction.cut){

    if(!basisavailable)
      stop("'fittedfunction' or 'fittedfunction.cut' = TRUE, but object of class peaklist contains
an empty basis matrix \n. Re-run with 'getPeaklist' with option 'returnbasis  = TRUE'. Note that this option is not supported if 'binning = TRUE. \n")
    
  if(fittedfunction & !fittedfunction.cut){
    fittedall <- drop(object@basis %*% bookfitted[,4])
    lines(x[ff], fittedall[ff], type = "l", ...)
  }

  if(!fittedfunction & fittedfunction.cut){
     if(!is.character(all.equal(object@peaklistprocessed, matrix(0))))
     bookpostprocessed <- object@peaklist
  else
     bookpostprocessed <- object@peaklistprocessed
     bookproind <- which(bookpostprocessed[,2] >= lower & bookpostprocessed[,2] <= upper)
     howmanyfunctions <- length(bookproind)
     amplpro <- bookpostprocessed[bookproind,4]
     chargepro <- bookpostprocessed[bookproind,3]
     positionspro <- bookpostprocessed[bookproind,2]
     if(length(positionspro) > 0){
     datapro <- cbind(positionspro, chargepro)
     if(model == "Gaussian")
       bascut <- calculatebasis.gaussian(x[ff], positions = positionspro, sigma = sigma, charges = charges, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff))
     if(model == "EMG")
       bascut <- calculatebasis.emg(x[ff], positions = positionspro, alpha = alpha, sigma = sigma, mu = mu, charges = charges, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff))
 
     basiscut <- bascut$Phi
     basiscutbook <- bascut$book
     sel <- c()
      for(i in 1:nrow(basiscutbook)){
       basiscutbook.i <- basiscutbook[i,2:3]
        for(j in 1:nrow(datapro)){
          if(all(basiscutbook.i == datapro[j,]))
            sel <- c(sel, i)
      }      
     }
     if(length(sel) > 0){
     permute <- match(basiscutbook[sel,2], datapro[,1])
     fittedcut <- drop(basiscut[,sel,drop = FALSE] %*% amplpro[permute])
     lines(x[ff], fittedcut, type = "l", ...)
   }
  }
 }    

   if(fittedfunction & fittedfunction.cut){
     fittedall <- drop(object@basis %*% bookfitted[,4])
     lines(x[ff], fittedall[ff], type = "l", main = "fitted function", ...)
     if(!is.character(all.equal(object@peaklistprocessed, matrix(0))))
       bookpostprocessed <- object@peaklist
     else
       bookpostprocessed <- object@peaklistprocessed
     bookproind <- which(bookpostprocessed[,2] >= lower & bookpostprocessed[,2] <= upper)
     howmanyfunctions <- length(bookproind)
     amplpro <- bookpostprocessed[bookproind,4]
     chargepro <- bookpostprocessed[bookproind,3]
     positionspro <- bookpostprocessed[bookproind,2]
     if(length(positionspro) > 0){
       datapro <- cbind(positionspro, chargepro)
     if(model == "Gaussian")
     bascut <- calculatebasis.gaussian(x[ff], positions = positionspro, sigma = sigma, charges = charges, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff))
     if(model == "EMG")
     bascut <- calculatebasis.emg(x[ff], positions = positionspro, alpha = alpha, sigma = sigma, mu = mu, charges = charges, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff))

     basiscut <- bascut$Phi
     basiscutbook <- bascut$book
     sel <- c()
      for(i in 1:nrow(basiscutbook)){
       basiscutbook.i <- basiscutbook[i,2:3]
        for(j in 1:nrow(datapro)){
          if(all(basiscutbook.i == datapro[j,]))
            sel <- c(sel, i)
      }      
     }
     if(length(sel) > 0){
     permute <- match(basiscutbook[sel,2], datapro[,1])
     fittedcut <- drop(basiscut[,sel, drop = FALSE] %*% amplpro[permute])
     lines(x[ff], fittedcut, type = "l", ...)
   }
   }
   }
    if(!is.null(quantile)) matlines(auxiliarysequence, localnoiselevel, lty = "dotted", col = "black") 
  }
}

  if(fitted){
    bookfittedind <- which(bookfitted[,2] >= lower & bookfitted[,2] <= upper)
    howmanyfitted <- length(bookfittedind)
    amplfitted <- bookfitted[bookfittedind,4]
    oo <- order(amplfitted, decreasing = TRUE)
    amplfitted <- amplfitted[oo]
    positionsfitted <- bookfitted[bookfittedind,2][oo]
    chargefitted <- bookfitted[bookfittedind,3][oo]

if(is.null(cutoff.eps)) looplimit <- min(cutoff.functions, howmanyfitted)
else looplimit <- sum(amplfitted > cutoff.eps)
if(looplimit > 0){
for(i in 1:looplimit){
 cand <- which(bookfitted[bookfittedind,2] == positionsfitted[i])
 colindex <- bookfittedind[cand[bookfitted[bookfittedind[cand],3] == chargefitted[i]]] 
 if(i == 1){
   if(basisavailable) plot(x[ff], amplfitted[i] * object@basis[ff,colindex], ylim = c(0, max(intensities[ff])), col = i, main = "fit", type = "l")
   else plot(x[ff], amplfitted[i] * basis$Phi[,colindex], ylim = c(0, max(intensities[ff])), col = i, main = "fit", type = "l")
 }
 else{
   if(basisavailable) lines(x[ff], amplfitted[i] * object@basis[ff,colindex], col = i, ...)
   else lines(x[ff], amplfitted[i] * basis$Phi[,colindex], col = i, ...)
 }
}
}
 else plot(x[ff], intensities[ff], main = "fitted", type = "n")   
    if(!is.null(quantile)) matlines(auxiliarysequence, localnoiselevel, lty = "dotted", col = "black")    
}

if(postprocessed){

  if(!is.character(all.equal(object@peaklistprocessed, matrix(0))))
     bookpostprocessed <- object@peaklist
  else
     bookpostprocessed <- object@peaklistprocessed
  
  opro <- order(bookpostprocessed[,2])
  bookpostprocessed <- bookpostprocessed[opro,,drop=FALSE]

  bookproind <- which(bookpostprocessed[,2] >= lower & bookpostprocessed[,2] <= upper)
  if(length(bookproind) > 0){
  howmanyfunctions <- length(bookproind)
  amplpro <- bookpostprocessed[bookproind,4]
  ooo <- order(amplpro, decreasing = TRUE)
  amplpro <- amplpro[ooo]
  positionspro <- bookpostprocessed[bookproind,2][ooo]
  chargepro <- bookpostprocessed[bookproind,3][ooo]
  if(model == "Gaussian")
  basispro <- calculatebasis.gaussian(x[ff], positions = positionspro, sigma = sigma, eps = 1e-05, uppernonzero = howmanyfunctions * max(charges) * sum(ff), charges = charges)
  if(model == "EMG")
    basispro <- calculatebasis.emg(x[ff], positions = positionspro, alpha = alpha, sigma = sigma, mu = mu, eps = 1e-05, uppernonzero = howmanyfunctions * ncharges * sum(ff), charges = charges)  
  if(is.null(cutoff.eps)) looplimit <- min(cutoff.functions, howmanyfunctions)
  else looplimit <- sum(amplpro > cutoff.eps)
  
  for(i in 1:looplimit){
   colindex <- apply(basispro$book[,2:3], 1, function(z) all(z == c(positionspro[i], chargepro[i])))
    if(i == 1){
      plot(x[ff], amplpro[i] * basispro$Phi[,colindex], col = i, main = "postprocessed", type = "l", ylim = c(0, max(intensities[ff])), ...)
 }
 else lines(x[ff], amplpro[i] * basispro$Phi[,colindex], col = i, ...)
}
}
else plot(x[ff], intensities[ff], type = "n", main = "postprocessed")  
  if(!is.null(quantile)) matlines(auxiliarysequence, localnoiselevel, lty = "dotted", col = "black") 
}
layout(matrix(1))
})

###############################################################################

setMethod("visualize", signature("modelfit", "missing", "missing"),
          function(object, type = c("peak", "model"),
                           parameters = c("alpha", "sigma", "mu"),
                           modelfit = FALSE, ...){
            
          type <- match.arg(type)  
            
          if(object@fitting == "most_intense"){
            if(type == "model")
              stop("Only the most intense peak has been fitted. In this case, only 'type = peak' is a valid option \n")
          }

          if(type == "peak"){
            bestpeak <- object@bestpeak

            if(object@model == "Gaussian"){
              x <- bestpeak$mz
              y <- bestpeak$intensities
              scale.y <- max(y)
              sigma <- bestpeak$sigma
              mu <- bestpeak$mu
              xgrid <- seq(from = min(x), to = max(x), length = 1000)
              plot(x, y, ...)
              lines(xgrid, gaussfun(xgrid, mu = mu, sigma = sigma) * scale.y)
          }
            if(object@model == "EMG"){
              x <- bestpeak$mz
              y <- bestpeak$intensities
              alpha <- bestpeak$alpha
              sigma <- bestpeak$sigma
              mu <- bestpeak$mu
              mutilde <- mu + mean(x)
              scale.mode <- determinemode(alpha, sigma, mu)$val
              scale.y <- max(y)
              scale.EMG <- scale.y/scale.mode    
              xgrid <- seq(from = min(x), to = max(x), length = 1000)
              plot(x, y, ...)
              lines(xgrid, EMG(xgrid, alpha = alpha, sigma = sigma, mu = mutilde) * scale.EMG)
            }
          }

          else{
            peakfitresults <- object@peakfitresults
            if(object@model == "Gaussian"){
              plot(peakfitresults[,"mz"], peakfitresults[,"sigma"], xlab = "m/z", ylab = expression(sigma(m/z)), ...)
              if(modelfit){
                foo <- object@sigmafunction
                xgrid <- seq(from = min(peakfitresults[,"mz"]), to = max(peakfitresults[,"mz"]))
                lines(xgrid, foo(xgrid))
              }
            }
            if(object@model == "EMG"){
              pars <- sort(unique(parameters[parameters %in% c("alpha", "sigma", "mu")]))
              lpars <- length(pars)
              opar <- par(mfrow = c(1,lpars))
              for(i in seq(along = pars)){
                yann <- paste(paste(pars[i], "m/z", sep="("), ")", sep = "")
                plot(peakfitresults[,"mz"], peakfitresults[,pars[i]], xlab = "m/z", ylab = parse(text = yann), ...)
                if(modelfit){
                  xgrid <- seq(from = min(peakfitresults[,"mz"]), to = max(peakfitresults[,"mz"]))
                  foo <- switch(pars[i], alpha = object@alphafunction,
                                         sigma = object@sigmafunction,
                                         mu = object@mufunction)
                  lines(xgrid, foo(xgrid))
                }
              }
              on.exit(par(opar))
            }
          }
          })
                   
          
          
          
          
          
          
          
          
          
          
          
