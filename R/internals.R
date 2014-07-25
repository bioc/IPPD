### 1. Computation of localnoise level: localnoise()

localnoise <- function(y, window, quantiledist = 0.1) {
    if ((window %% 2) == 0) {
        stop("window width should be odd \n")
    }

    quantileseq <- seq(from = 0, to = 1, by = quantiledist)
    qindex <- pmax(floor(window * quantileseq), 1)
    nquantiles <- length(qindex)
    halfwindow <- max(floor(window/2), 1)
  
    n <- length(y)
   
    q <- matrix(nrow = n, ncol = nquantiles, data = 0)
   
    indexsort <- rank(y[1:window], ties.method = "first")
    sorty <- sort(y[1:window])
  
    q[1:(halfwindow + 1),] <- rep(sorty[qindex], each= halfwindow + 1)
  
  
    Cres <- .C("localquantile",
               y = as.double(y),
               sorty = as.double(sorty),
               indexsort = as.integer(indexsort - 1),
               window = as.integer(window),
               halfwindow = as.integer(halfwindow),
               n = as.integer(n),
               flag = as.integer(0),
               q = as.double(q),
               qindex = as.integer(qindex - 1),
               nquantiles = as.integer(nquantiles))

    q <- matrix(Cres$q, nrow = n, ncol = nquantiles)              

    q[((n- halfwindow)+1):n,] <- rep(q[n - halfwindow,], each = halfwindow)

    colnames(q) <- as.character(round(quantileseq, 2))

    return(q)
}

### 2. calculation of peakheights according to the averagine model

getpeakheights <- function(mass, averagine.table) {
    data(averagine.table)
    tabb <- as.matrix(tableaveragine)
    endpoints <- tabb[,1]
    amplitudes <- tabb[,-1]

    ### currently: at most 31 peaks
    peakheights <- matrix(0, nrow = length(mass), ncol = 31)
    dismass <- as.numeric(cut(mass, breaks = c(-Inf, endpoints, Inf)))
    peakheights <- .C("interpolatepeakheights",
                      peakheights = peakheights,
                      dismass = as.integer(dismass),
                      mass = as.double(mass),
                      endpoints = as.double(endpoints),
                      amplitudes = amplitudes,
                      maxnpeaks = as.integer(31),
                      numbermass = as.integer(length(mass)),
                      maxdismass = as.integer(length(endpoints)))$peakheights

    peakheights
}

### 3. Computation of Phi-matrix for model = Gaussian

calculatebasis.gaussian <- function(x, positions = x, sigma, charges = c(1,2,3,4),
                                    eps = 1e-05, uppernonzero = 10^7, averagine.table) {
    n <- length(x)

    positions <- unique(positions)
    positions <- sort(positions)

    npositions <- length(positions)

    kappa <- 1.008
    book <- NULL

    for (l in seq(along = charges)) {
        chargel <- charges[l]
        amplitudes <- getpeakheights(positions * chargel, averagine.table)

        ### ell-infinity normalization 
        amplitudes <- t(apply(amplitudes, 1, function(z) {maxz = max(z); c(z/maxz, maxz)}))
        ###

        npeaks <- ncol(amplitudes) - 1
        ainf <- amplitudes[, (npeaks + 1), drop = FALSE] 
        centers <- positions + kappa * t(apply(amplitudes, 1, 
                                               function(z) {
                                                   origin <- which.max(z)
                                                   
                                                   if ((origin > 1) & (origin < npeaks)) { 
                                                       ret <- c((-1)*((origin-1):1), 0, 1:(npeaks-origin))
                                                   }
                                                   
                                                   if (origin == 1) {
                                                       ret <- (0:(npeaks - 1))
                                                   }

                                                   if ((origin == npeaks)) {
                                                       ret <- (-1) * ((npeaks-1):0)
                                                   }
                                                   
                                                   ret})) / chargel
        initials <- apply(centers, 1, min)
        ### filtering (optional)
        #if(filter){
        #
        #  filtered <- .C("mzfilter", positions = as.double(initials), npositions = as.integer(length(initials)),
        #                 charge = as.integer(chargel),  filteredout = integer(length(initials)))
        #   centers[as.logical(filtered$filteredout),] <- 0
        # }

        #if(!is.null(epslocal) & !is.null(y)){
        #  if(length(epslocal) != length(x))
        #    stop("Length of 'epslocal' has to be equal to the length of 'x' \n")
        #  start <- apply(amplitudes, 1, which.max) - 1
        #  centersC <- .C("patterntruncate", centers = centers, amplitudes = amplitudes, x = as.double(x), y = as.double(y),
        #                localeps = as.double(epslocal), start = as.integer(start), n = as.integer(n), npostions = as.integer(npositions),
        #                npeaks = as.integer(npeaks), flag = as.integer(0))
        #  centers <- matrix(nrow = npositions, ncol = npeaks, data = centersC$centers)
        #}
        #else{
        centers[centers > max(x)] <- 0
        centers[centers < min(x)] <- 0
#}
        block <- matrix(0, nrow = uppernonzero, ncol = 3)

        sigmal <- sigma(positions)
   
        rangecenters <- apply(centers, 1, function(z) {
            if (all(z == 0))
                return(c(0,0))
            else
                return(range(z[z > 0]))})

        fence <- sqrt(-log(eps) * sigmal)
    
        resultC <- .C("gaussbasis", 
                       block = block,
                       x = as.double(x),
                       amplitudes = amplitudes,
                       centers = centers,
                       npositions = as.integer(npositions),
                       n = as.integer(n),
                       maxnpeaks = as.integer(npeaks),
                       sigma = as.double(sigmal),
                       nonzero = as.integer(-1),
                       eps = as.double(eps),
                       uppernonzero = as.integer(uppernonzero),
                       lower = rangecenters[1,] - fence,
                       upper = rangecenters[2,] + fence,
                       colmax = double(npositions))

        if ((resultC$nonzero +1) == 0)
            next

        if (!exists("Phi")) {
            temp <- sparseMatrix(i = as.integer(resultC$block[1:(resultC$nonzero + 1), 1]),
                                 j = as.integer(resultC$block[1:(resultC$nonzero + 1), 2]),
                                 x = resultC$block[1:(resultC$nonzero + 1), 3],
                                 dims = c(n, npositions))

            remove <- resultC$colmax == 0
            Phi <- temp[,!remove,drop = FALSE] %*% Diagonal(n = sum(!remove), 1/resultC$colmax[!remove])
        } else {
            temp <- sparseMatrix(i = as.integer(resultC$block[1:(resultC$nonzero + 1), 1]),
                                 j = as.integer(resultC$block[1:(resultC$nonzero + 1), 2]),
                                 x = resultC$block[1:(resultC$nonzero + 1), 3],
                                 dims = c(n, npositions))
       
            remove <- resultC$colmax == 0
            temp <- temp[,!remove,drop = FALSE] %*% Diagonal(n = sum(!remove), 1/resultC$colmax[!remove])
            Phi <- cbind2(Phi, temp) 
        }
        ###
        bookl <- cbind(initials, positions, chargel, ainf)[!remove,,drop = FALSE]
        book <- rbind(book, bookl)
    }

    colnames(book) <- c("initial", "most_intense", "charge", "a_inf")
    return (list(Phi = Phi, book = book))
}


### 4. Computation of Phi-matrix for model = EMG

calculatebasis.emg <- function(x, positions = x, alpha, sigma, mu, charges = c(1,2,3,4),
                               eps = 1e-05, uppernonzero = 10^7, averagine.table) {
    n <- length(x)
    delta <- max(diff(x))
    positions <- unique(positions)
    positions <- sort(positions)
    npositions <- length(positions)

    alpha <- alpha(positions)
    sigma <- sigma(positions)
    mu <- mu(positions)
  
    moderes <- apply(cbind(alpha, sigma, mu), 1, function(z) {
                                                    getmode <- determinemode(z[1], z[2], z[3])
                                                    c(-getmode$mode + z[3], 1/getmode$val)})

    mu <- moderes[1,] 
    scale <- moderes[2,] 
    kappa <-  1.008 
    book <- NULL
  
    for (l in seq(along = charges)) {
        chargel <- charges[l]
        amplitudes <- getpeakheights(positions * chargel, averagine.table)
        amplitudes <- t(apply(amplitudes, 1, function(z) {maxz = max(z); c(z/maxz, maxz)}))
        ###
        npeaks <- ncol(amplitudes) - 1
        ainf <- amplitudes[, (npeaks + 1), drop = FALSE] 
        centers <- positions + kappa * 
                    t(apply(amplitudes, 1, function(z) {
                        origin <- which.max(z)
                        if ((origin > 1) & (origin < npeaks)) 
                            ret <- c((-1)*((origin-1):1), 0, 1:(npeaks-origin))
                        if (origin == 1)
                            ret <- (0:(npeaks - 1))
                        if ((origin == npeaks))
                            ret <- (-1) * ((npeaks-1):0)
                        ret}))/chargel
        initials <- apply(centers, 1, min) 
    #if(filter){
    #  filtered <- .C("mzfilter", positions = as.double(initials), npositions = as.integer(length(initials)),
    #                 charge = as.integer(chargel),  filteredout = integer(length(initials)))
    #   centers[as.logical(filtered$filteredout),] <- 0
    # }
    #if(!is.null(epslocal) & !is.null(y)){
    #  if(length(epslocal) != length(x))
    #    stop("Length of 'epslocal' has to be equal to the length of 'x' \n")

    #  start <- apply(amplitudes, 1, which.max) - 1
    #  centersC <- .C("patterntruncate", centers = centers, amplitudes = amplitudes, x = as.double(x), y = as.double(y),
    #                localeps = as.double(epslocal), start = as.integer(start), n = as.integer(n), npostions = as.integer(npositions),
    #                npeaks = as.integer(npeaks), flag = as.integer(0))
    #  centers <- matrix(nrow = npositions, ncol = npeaks, data = centersC$centers)
    #}
    #else{
     
        centers[centers > max(x)] <- 0
        centers[centers < min(x)] <- 0
   #}
    
        rangecenters <- apply(centers, 1, function(z) {
                if (all(z == 0))
                    return (c(0,0))
                else 
                    return (range(z[z > 0]))})

        fence <- sigma^2 / (2 * alpha) + mu - alpha * log(alpha * eps/(sqrt(2*pi) * sigma))

        block <- matrix(0, nrow = uppernonzero, ncol = 3)
        resultC <- .C("emgbasis",
                      block = block,
                      x = as.double(x),
                      amplitudes = amplitudes,
                      centers = centers,
                      npositions = as.integer(npositions),
                      n = as.integer(n),
                      maxnpeaks = as.integer(npeaks),
                      alpha = as.double(alpha),
                      sigma = as.double(sigma),
                      mu = as.double(mu),
                      nonzero = as.integer(-1),
                      eps = as.double(eps),
                      uppernonzero = as.integer(uppernonzero),
                      lower = rangecenters[1,] - fence,
                      upper = rangecenters[2,] + fence, colmax = double(npositions)) 

        if ((resultC$nonzero +1) == 0)
            next 
    
        if (!exists("Phi")) {
            temp <- sparseMatrix(i = as.integer(resultC$block[1:(resultC$nonzero + 1), 1]),
                                 j = as.integer(resultC$block[1:(resultC$nonzero + 1), 2]),
                                 x = resultC$block[1:(resultC$nonzero + 1), 3],
                                 dims = c(n, npositions))
            remove <- resultC$colmax == 0
            Phi <- temp[,!remove,drop = FALSE] %*% Diagonal(n = sum(!remove), 1/resultC$colmax[!remove])
        } else{
            temp <- sparseMatrix(i = as.integer(resultC$block[1:(resultC$nonzero + 1), 1]),
                                 j = as.integer(resultC$block[1:(resultC$nonzero + 1), 2]),
                                 x = resultC$block[1:(resultC$nonzero + 1), 3],
                                 dims = c(n, npositions))
       
            remove <- resultC$colmax == 0
            temp <- temp[,!remove,drop = FALSE] %*% Diagonal(n = sum(!remove), 1/resultC$colmax[!remove])
            Phi <- cbind2(Phi, temp) 
        }

        bookl <- cbind(initials, positions, chargel, ainf)[!remove,,drop = FALSE]
    
        book <- rbind(book, bookl)  
    }

    colnames(book) <- c("initial", "most_intense", "charge", "a_inf")
    return (list(Phi = Phi, book = book)) 
}

### 5. Nonnegative least squares

nnlslogbarrier <- function(response, betastart, trace = FALSE, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-6){ 
  betaold <- betastart
  gamma <- gammastart
  rss <- Inf
 while(gamma < gammamax){
  converged <- FALSE 
  if(trace) cat("gamma: ", gamma, "\n")
  while(!converged){
  if(any(betaold < 0)) stop("'beta' has to be non-negative \n")
   gradf <- (drop(get("G", parent.frame(1)) %*% betaold) - get("C", parent.frame(1)))
   grad <-   gradf + (-1/betaold)/gamma
   diagadd <- .C("addiagonal", colptr = as.integer(get("G", parent.frame(1))@p), rowind = as.integer(get("G", parent.frame(1))@i), val = as.double(get("G", parent.frame(1))@x), add = as.double( (1/betaold^2 )/ gamma), p = as.integer(get("G", parent.frame(1))@Dim[1]))
   Hess <- new("dsCMatrix", uplo = "U", i = diagadd$rowind, p = diagadd$colptr, x = diagadd$val, Dim = get("G", parent.frame(1))@Dim, factors = list())
  
   R <- try(Cholesky(Hess), silent = TRUE)
   if(inherits(R, "try-error")){
      warning("Error in nonnegative least squares estimation: Hessian not positive definite.  \n")
     gamma <- gammamax
     break  
   }
   dir <- solve(R, -grad)@x
   stepsize <- try(linesearch(response, alpha), silent = TRUE)
   if(inherits(stepsize, "try-error")){
     warning("Error in nonnegative least squares estimation: stepsize selection failed \n")
     gamma <- gammamax
     break    
   }
   if(trace) cat("stepsize: ", stepsize, "\n")
   beta <- betaold + dir * stepsize
   if(trace) cat("convergence inner loop:", abs(sum(grad*drop(dir))), "\n")
   if((stepsize == 1) | abs(sum(grad*drop(dir))) < eps){
        converged <- TRUE
        gamma <- gamma * gammamult
      }
   betaold <- beta
  }
 }
  if(is.function(beta)) beta <- betaold 
 return(list(beta = beta, gradnorm = sqrt(sum(grad * grad)), rss = 2 * rss))
}

### 6. linesearch [for nonnegative least squares]

linesearch <- function(response, alpha){
  stepsize <- 1
  sigma <- 0.95
  sigmavec <- sigma^(0:500)
  whichsmaller0 <- get("dir", parent.frame(1)) < 0
  if(sum(whichsmaller0) > 0)
    maxstepsize <- min(- (get("betaold", parent.frame(1)) / get("dir", parent.frame(1)))[whichsmaller0])
  else maxstepsize <- 1
  if(maxstepsize < 1)
      maxstepsize <- max(sigmavec[sigmavec < maxstepsize])
  stepsize <- min(c(1, maxstepsize)) 
  finished <- FALSE
  objectiveold <- 0.5 * sum((response - get("basis", parent.frame(2)) %*% get("betaold", parent.frame(1)))^2) - (1/get("gamma", parent.frame(1))) * sum(log(get("betaold", parent.frame(1))))
  while(!finished){
    if(get("trace", parent.frame(1))) cat("linesearch...", stepsize, "\n")
    beta <-  get("betaold", parent.frame(1)) + stepsize * get("dir", parent.frame(1))
    objectiverss <-  0.5 * sum((response - get("basis", parent.frame(2)) %*% beta)^2)
    objectivelogbarrier <- -(1/get("gamma", parent.frame(1))) * sum(log(beta))
    objective <- objectiverss + objectivelogbarrier
    if(is.nan(objective)) objective <- Inf
    rhs <- (objectiveold + alpha * stepsize * sum(get("grad", parent.frame(1)) * get("dir", parent.frame(1))))
    if(  objective >= rhs){  
     stepsize <- sigma * stepsize
     if(stepsize < (0.5 * .Machine$double.eps)){
       stop("Canceled stepsize selection. \n")
     }
     }  
    else finished <- TRUE
  }
  if(get("trace", parent.frame(1))){ cat("objectiverss:", objectiverss, "\t") ; cat("objectivelogbarrier:", objectivelogbarrier, "\n")}
  assign("rss", objectiverss, parent.frame(1))
  stepsize
}

### 7. postprocessing Gaussian

postprocess.gaussian <- function(ppm, sigma, peaklist, eps = 1e-10, maxruns = 10, trace = FALSE){


orr <- order(peaklist[,2])

peaklisto <- peaklist[orr,,drop = FALSE]


sigmas <- sigma(peaklisto[,2])
          
  
peaklisto <- cbind(peaklisto, sigmas)

peaklistrem <- NULL
peaklistnew <- NULL


l <- 1
while(l <= nrow(peaklisto)){
  indexplus = 0
  if((l + indexplus + 1) > nrow(peaklisto)){
    peaklistrem <- rbind(peaklistrem, peaklisto[l,])
    break
  }
  else dif <- (peaklisto[l+ indexplus + 1,2] - peaklisto[l + indexplus,2])
  
  while(dif < eps){
    indexplus <- indexplus + 1
    if((l + indexplus + 1) > nrow(peaklisto)) break
    else dif <- (peaklisto[l+ indexplus + 1,2] - peaklisto[l + indexplus,2]) 
  }
  if(indexplus == 0){
    peaklistrem <- rbind(peaklistrem, peaklisto[l,])
    l <- l + 1
  }
  else{
     whichmaxampl <- which.max(peaklisto[l:(l+indexplus),5])
     indtoadd <- (l:(l + indexplus))[whichmaxampl]
     peaklistrem <- rbind(peaklistrem, peaklisto[indtoadd,])
     l <- l + indexplus + 1
 }   
}

outer <- 1
converged <- FALSE
while(!converged & outer <= maxruns){ 
if(trace) cat("postprocessing: round ",  outer, "\n")
focus <- 1:nrow(peaklistrem)
while(length(focus) >= 1){
  loc <- focus[1]
  locmz <- peaklistrem[loc,2]
  delta <- ppm * locmz / 10^6
  indeltawindow <- focus[(peaklistrem[focus,2] - locmz) <= delta]
  s_indeltawindow <- length(indeltawindow)
  if(s_indeltawindow == 1){ 
    peaklistnew <- rbind(peaklistnew, peaklistrem[loc,])
    focus <- setdiff(focus, loc)
  }
  else{
     chargeloc <- peaklistrem[loc,3]
     samecharge  <- which(peaklistrem[indeltawindow,3] ==  chargeloc)
     if(length(samecharge) == 0){
       peaklistnew <- rbind(peaklistnew, peaklistrem[loc,])
       focus <- setdiff(focus, loc)
     }
     else{
       collect <- indeltawindow[samecharge]
       initial <- peaklistrem[collect, 1]
       mus <- peaklistrem[collect,2]
       betas <- peaklistrem[collect,5]
       ainf <- peaklistrem[collect, 4]
       sigmaloc <- peaklistrem[loc,6]
       img <- intermediategaussian(mus, betas,  sigma = sigmaloc)
       shift <- img$mu - mus[which.max(betas)]
       initialnew <- initial[which.max(betas)] + shift 
       peaklistnew <- rbind(peaklistnew, c(initialnew, img$mu, chargeloc, mean(ainf), img$beta, sigmaloc))    
       focus <- setdiff(focus, collect)
     }
  }
  
}
  if(nrow(peaklistnew) == nrow(peaklistrem)) converged <- TRUE
  peaklistrem <- peaklistnew
  peaklistnew <- NULL
  outer <- outer + 1

}

return(peaklistrem[,-6, drop = FALSE])

}

### 8. intermediategaussian [for postprocess.Gaussian]

intermediategaussian <- function(mus, betas, muinit = NULL, sigma, eps = 1e-6){
  N <- length(mus) 
  if(length(betas) != N) stop("Length of 'betas' does not match length of 'mus' \n")
  if(all(betas == 0)) stop("All betas are equal to zero \n")
  if(is.null(muinit)){
    wsum <- sum(betas)
    muinit <- sum(betas*mus)/wsum
  }
  mu0 <- muinit
  beta0 <- mean(betas)
  ws <- gaussfun(mus, mu = mu0, sigma = 2 * sigma) 
  iter <- 0
  maxiter <- 100
  converged <- FALSE
  while(!converged & iter <= maxiter){
    terms <- betas * ws
    beta0new <-  sum(terms)
    mu0new <- sum(terms * mus) / beta0new
    wsnew <- gaussfun(mus, mu = mu0new,  sigma = 2 * sigma)   
    if(all(c(abs(mu0new - mu0), abs(beta0new - beta0), abs(wsnew - ws)) < eps) | iter == maxiter){
      converged <- TRUE
      eq1 <- abs(beta0new - sum(betas * wsnew))
      eq2 <- abs(sum(betas * mus * wsnew) - mu0new * beta0new)
      #if(!exists("eq1")) browser()
      #bug <- try(print(eq1), silent = TRUE)
      #if(inherits(bug, "try-error")) browser()
      if(eq1 > eps) converged <- FALSE
      if(eq2 > eps) converged <- FALSE
    }
    if(iter == (maxiter-1)) warning("Maximum iteration count reached \n")
    mu0 <- mu0new
    beta0 <- beta0new
    ws <- wsnew
    iter <- iter + 1
    #if(!exists("eq1")) browser()
    #bug <- try(print(eq1), silent = TRUE)
    #if(inherits(bug, "try-error")) browser()
  }
  return(list(mu = mu0, beta = beta0, eps1 = eq1, eps2 = eq2, iter = iter))
}



### 9. gaussfun [for intermediategaussian]

gaussfun <- function(z, mu, sigma){
   factor <- sqrt(2*pi) * sqrt(sigma/2)
   factor * dnorm(z, mean = mu, sd = sqrt(sigma/2))
}


### 10. postprocessing EMG

postprocess.emg <- function(ppm, delta.x, alpha, sigma, mu, peaklist, eps = 1e-10, maxruns = 10, npoints = 100, trace = FALSE){
  

orr <- order(peaklist[,2])

peaklisto <- peaklist[orr,,drop = FALSE]


peaklistrem <- NULL
peaklistnew <- NULL

l <- 1
while(l <= nrow(peaklisto)){
  indexplus = 0
  if((l + indexplus + 1) > nrow(peaklisto)){
    peaklistrem <- rbind(peaklistrem, peaklisto[l,])
    break
  } 
  else dif <- (peaklisto[l+ indexplus + 1,2] - peaklisto[l + indexplus,2])
  
  while(dif < eps){
    indexplus <- indexplus + 1
    if((l + indexplus + 1) > nrow(peaklisto)) break
    else dif <- (peaklisto[l+ indexplus + 1,2] - peaklisto[l + indexplus,2]) 
  }
  if(indexplus == 0){
    peaklistrem <- rbind(peaklistrem, peaklisto[l,])
    l <- l + 1
  }
  else{
    
     whichmaxampl <- which.max(peaklisto[l:(l+indexplus),4])
     indtoadd <- (l:(l + indexplus))[whichmaxampl]
     peaklistrem <- rbind(peaklistrem, peaklisto[indtoadd,])
     l <- l + indexplus + 1
 }   
}

outer <- 1
converged <- FALSE

while(!converged & outer <= maxruns){
 
if(trace) cat("postprocessing: round ",  outer, "\n")
focus <- 1:nrow(peaklistrem)
while(length(focus) >= 1){
  loc <- focus[1]
  locmz <- peaklistrem[loc,2]
  delta <- locmz * ppm / 10^6
  indeltawindow <- focus[(peaklistrem[focus,2] - locmz) <= delta]
  s_indeltawindow <- length(indeltawindow)
  if(s_indeltawindow == 1){ 
    peaklistnew <- rbind(peaklistnew, peaklistrem[loc,])
    focus <- setdiff(focus, loc)
  }
  else{
     chargeloc <- peaklistrem[loc,3]
     samecharge  <- which(peaklistrem[indeltawindow,3] ==  chargeloc)
     if(length(samecharge) <= 1){
       peaklistnew <- rbind(peaklistnew, peaklistrem[loc,])
       focus <- setdiff(focus, loc)
     }
     else{
       collect <- indeltawindow[samecharge]
       initial <- peaklistrem[collect, 1]
       betas <- peaklistrem[collect,5]
       alphaloc <- alpha(peaklistrem[loc,2])
       sigmaloc <- sigma(peaklistrem[loc,2])
       muloc <- mu(peaklistrem[loc,2])
       modelocres <- determinemode(alphaloc, sigmaloc, muloc)
       shiftloc <- -modelocres$mode
       scaleloc <- 1/modelocres$val
       ainfloc <- mean(peaklistrem[collect,4])
       mus <- peaklistrem[collect,2] + (muloc + shiftloc)
      
       lowmu <- min(mus)
       uppmu <- max(mus)
       fence <- sigmaloc^2 / (2 * alphaloc) + (muloc + shiftloc) - alphaloc * log(alphaloc * sqrt(eps)/(sqrt(2*pi) * sigmaloc))
       low <-  lowmu - fence            
       upp <- max(mus) + fence

       ### find better bounds using uniroot
       if(alphaloc > 1.5 * sigmaloc)
       low.interval.left  <- low
       low.interval.right <- min(peaklistrem[collect,2])  
       low.corr <- try(uniroot(function(z) EMG(z, alpha = alphaloc, sigma = sigmaloc, mu = lowmu) - sqrt(eps),
                           lower = low.interval.left, upper = low.interval.right)$root, silent = TRUE)
       if(!inherits(low.corr, "try-error")) low <- low.corr
       ###
       npointsfac <- max(1, floor(delta / delta.x))
       
       imemg <- intermediateemg(mus, betas * scaleloc,  alpha = alphaloc, sigma = sigmaloc, npoints = npoints * npointsfac, low = low, upp = upp)
       shiftmu <- imemg$mu - mus[which.max(betas)]
       initialnew <- initial[which.max(betas)] + shiftmu
       newrow <- c(initialnew, imemg$mu - muloc - shiftloc, chargeloc, ainfloc, imemg$beta / scaleloc)
       acceptnewrow = TRUE
       if(any(newrow[5] + sqrt(eps) < betas)){
         acceptnewrow <- FALSE
         #browser()
       }
         
       #if(newrow[2] < min(mus) | nrow[2] > max(mus))
       #  acceptnewrow <- FALSE
       if(acceptnewrow) peaklistnew <- rbind(peaklistnew, newrow)
       else peaklistnew <- rbind(peaklistnew, peaklistrem[collect,])
       focus <- setdiff(focus, collect)
     }
  }
  
}
  if(nrow(peaklistnew) == nrow(peaklistrem)) converged <- TRUE
  peaklistrem <- peaklistnew
  peaklistnew <- NULL
  outer <- outer + 1

}

return(peaklistrem)

}

### 11. erfc

erfc <- function(z) pnorm(z*sqrt(2), lower.tail = FALSE)

### 12. P1
P1 <- function(x, alpha, sigma, mu) exp(sigma^2/(2*alpha^2) + (mu - x)/alpha)
### 13. P2
P2 <- function(alpha, sigma) sqrt(2 * pi) * sigma/alpha
### 14. P3
P3 <- function(x, alpha, sigma, mu){
 z <- (sigma/alpha + (mu - x)/sigma)/sqrt(2)
 erfc(z)
}
### 15. P3prime
P3prime <- function(x, alpha, sigma, mu){
 z <- sigma/alpha + (mu - x)/sigma
 dnorm(z)
}
### 16. EMG
EMG <- function(x, alpha = 1, sigma = 1, mu = 0) P1(x, alpha, sigma, mu) * P2(alpha, sigma) * P3(x, alpha, sigma, mu)
### 17. determinemode
determinemode <- function(alpha, sigma, mu, lower = -0.5, upper = 0.5){

  opt <- optimize(f = EMG, alpha = alpha, sigma = sigma, mu = mu, lower = lower,
              upper = upper, maximum = TRUE, tol = .Machine$double.eps^0.5)

  list(mode = opt$maximum, val = opt$objective)
}
### 18. dalpha
dalpha <- function(x, alpha, sigma, mu){
 P1(x, alpha, sigma, mu) * P2(alpha, sigma) * (P3(x, alpha, sigma, mu)*((x - mu)/alpha^2 - 1/alpha - sigma^2/(alpha^3)) + P3prime(x, alpha, sigma, mu) * sigma/alpha^2)
}
### 19. dsigma
dsigma <- function(x, alpha, sigma, mu){
 P1(x, alpha, sigma, mu) * P2(alpha, sigma) * (P3(x, alpha, sigma, mu)*(1/sigma + sigma/(alpha^2)) - P3prime(x, alpha, sigma, mu) * (1/alpha + (x - mu)/sigma^2))
}


### 20.dmu 
dmu <- function(x, alpha, sigma, mu){
 P1(x, alpha, sigma, mu) * P2(alpha, sigma) * (P3(x, alpha, sigma, mu)/alpha - P3prime(x, alpha, sigma, mu)/sigma) 
}

### 21. objective
objective <- function(theta, x, y){
 alpha <- theta[1]
 sigma <- theta[2]
 mu <- theta[3]
 y <- y/max(y)
 x <- x - mean(x)
 yhat <- EMG(x, alpha, sigma, mu)
 yhat <- yhat/max(yhat)
 sum((y - yhat)^2)
}

### 22.grad
grad <- function(theta, x, y){
 alpha <- theta[1]
 sigma <- theta[2]
 mu <- theta[3]
 y <- y/max(y)
 x <- x - mean(x)
 yhat <- EMG(x, alpha, sigma, mu)
 yhat <- yhat/max(yhat)
 res <- (y - yhat)
 gradalpha <- dalpha(x, alpha, sigma, mu) 
 gradsigma <- dsigma(x, alpha, sigma, mu)
 gradmu <- dmu(x, alpha, sigma, mu)
 -2 * c(sum(res * gradalpha), sum(res * gradsigma), sum(res * gradmu))  
}

### 23. gridsearch
gridsearch <- function(x, y, grid.length = c(100, 100, 100), grid.alpha = NULL, grid.sigma = NULL, grid.mu = NULL){
  
   if(is.null(grid.alpha)) grid.alpha <- 10^((seq(from = -5, to = 5, length = floor(grid.length[1]))))
   if(is.null(grid.sigma)) grid.sigma <- 10^((seq(from = -5, to = 5, length = floor(grid.length[2]))))
   if(is.null(grid.mu)) grid.mu <- seq(from = -1, to = 1, length = floor(grid.length[3]))
   gridsearchres <- .C("gridsearchemg", x = x, y = y, alpha = grid.alpha, sigma = grid.sigma, mu = grid.mu, n = as.integer(length(x)),
       alphagridlength = as.integer(length(grid.alpha)), sigmagridlength = as.integer(length(grid.sigma)), mugridlength = as.integer(length(grid.mu)),
       alphawinner = as.double(1), sigmawinner = as.double(1), muwinner = as.double(0))  

  start <- c(gridsearchres$alphawinner, gridsearchres$sigmawinner, gridsearchres$muwinner)

  return(start)
}

### 24. fit.EMG
fit.EMG <- function(x, y, alpha.start = 0.03, sigma.start = 0.03, mu.start = -0.03, gridsearch = FALSE, ...){
  if(gridsearch) start <- gridsearch(x = x - mean(x), y = y/max(y), ...)
  else start <- c(alpha.start, sigma.start, mu.start)
  optcall <- optim(par = start, fn = objective, gr = grad, method = "BFGS", x=x, y=y)
  if(optcall$convergence != 0) warning("Convergence failed \n")
  outp <- list(alpha = optcall$par[1], sigma = optcall$par[2], mu = optcall$par[3])
  yhat <- EMG(x - mean(x), alpha = outp$alpha, sigma = outp$sigma, mu = outp$mu)
  outp$yhat <- yhat/max(yhat)
  outp$rss <- sum((yhat/max(yhat) - y/max(y))^2)
  outp$converged <- optcall$convergence == 0
  outp$meanx <- mean(x)
  return(outp)
}

### 25. objective.combine
objective.combine <- function(theta, alpha, sigma, x, y){
 mu <- theta[1]
 beta <- theta[2]
 yhat <- beta * EMG(x, alpha, sigma, mu)
 sum((y - yhat)^2)
}

### 26. grad.combine
grad.combine <- function(theta, alpha, sigma, x, y){
 mu <- theta[1]
 beta <- theta[2]
 yhat <- beta * EMG(x, alpha, sigma, mu)
 res <- (y - yhat)
 gradmu <- beta * dmu(x, alpha, sigma, mu)
 gradbeta <- EMG(x, alpha, sigma, mu)
 -2 * c(sum(res * gradmu), sum(res * gradbeta))  
}

### 27. fit.combine
fit.combine <- function(x, y, alpha = 0.03, sigma = 0.03, mu.start = -0.03, beta.start = 1, lower, upper){
  start <- c(mu.start, beta.start)
  optcall <- optim(par = start, fn = objective.combine, gr = grad.combine, method = "L-BFGS-B", x = x, y = y, alpha = alpha, sigma = sigma,
                   lower = lower, upper = upper)
  if(optcall$convergence != 0) warning("Convergence failed \n")
  outp <- list(mu = optcall$par[1], beta = optcall$par[2])
  outp$converged <- optcall$convergence == 0
  return(outp)
}

### 28. intermediatemg [for postprocess.emg]

intermediateemg <- function(mus, betas, muinit = NULL, betainit = NULL, alpha, sigma, npoints = 100, low, upp, ...){
  N <- length(mus) 
  if(length(betas) != N) stop("Length of 'betas' does not match length of 'mus' \n")
  if(all(betas == 0)) stop("All betas are equal to zero \n")
  if(is.null(muinit)) muinit <- mean(mus)
  if(is.null(betainit)) betainit <- max(betas)

  xx <- seq(from = low, to = upp, length = npoints)
  yy <- drop( sapply(mus, function(z) EMG(xx, alpha = alpha, sigma = sigma, mu = z )) %*% betas)

  outp <- fit.combine(xx, yy, alpha = alpha, sigma = sigma, mu.start = muinit, beta.start = betainit,
                      lower = c(min(mus), max(betas)), upper = c(max(mus), sum(betas)))
  if(!outp$converged) warning("Nonlinear least squares did not converge \n")
  return(list(mu = outp$mu, beta = outp$beta, converged = outp$converged))
}

### 29. peakdetection

peakdetect <- function(spectrum, window, threshold){
  l <- length(spectrum[,1])
  delta <- diff(spectrum[,2])/diff(spectrum[,1]) 
  Cres <- .C("peakdetect", spectrum = spectrum[,2], delta = delta,
             l = as.integer(l), window = as.integer(window),
             threshold = threshold, peakind = integer(l),
             counter = as.integer(0))
  if(Cres$counter == 0) stop("No peaks found \n")
  else{
    peakind <- Cres$peakind[1:Cres$counter]+1
    peaklist <- vector(mode = "list", length = length(peakind))
    for(j in seq(along = peakind)){
     bestpeak <- peakind[j] ####+ window + 1
     xpeak <- spectrum[bestpeak,1]  
     xenvelopinit <- xpeak - 1
     xenvelopend <- xpeak + 1
     indinit <- which(spectrum[,1] >= xenvelopinit & spectrum[,1] < xpeak)
     ###indinit <- c(indinit[1]-1, indinit)
     indend <-  which(spectrum[,1] > xpeak  & spectrum[,1] <= xenvelopend)
     ###inend <- c(indend, indend[length(indend)] + 1)
     diffinit <- rev(delta[indinit]) ### rev(diff(spectrum[indinit,2])/diff(spectrum[indinit,1]))
     diffend <-  delta[indend - 1] ###diff(spectrum[indend,2])/diff(spectrum[indend,1])
     xxind <- bestpeak###, bestpeak+1)
      ### go left
      ileft <- 1
        while(diffinit[ileft] > 0){
          xxind <- c(bestpeak - ileft, xxind)
          ileft <- ileft +1
          if(ileft > length(diffinit)) break
         }     
      ### go right
     iright <- 1 
      while(diffend[iright] < 0){
          xxind <- c(xxind, bestpeak + iright)
          iright <- iright +1
          if(iright > length(diffend)) break
         }
       peaklist[[j]] <- cbind(x = spectrum[xxind,1], y = spectrum[xxind,2])
    }
  }
   return(peaklist)
}


### 30. fit

fit.gauss <- function(x,y,...){
  
  data <- data.frame(x = x - mean(x), yy = y/max(y))
  nls.fit <- nls(yy ~ exp(-(x-mu)^2/sigma), 
       data = data, start = list(sigma = 0.001, mu = 0))
  rss <- sum(residuals(nls.fit)^2)
  return(list(maxy = max(y), mu = coef(nls.fit)["mu"] + mean(x),
              sigma = coef(nls.fit)["sigma"], rss = rss))
}

### 31. formulacoef2function

formulacoef2function <- function(formula, coef, intercept){
 char <- as.character(formula[2])
 splitt <- unlist(strsplit(char, split = ""))
 critical <- c("+", "-", "1")
 matchcritical <- which(splitt %in% critical)
 criticalpoints <- sort(c(0, matchcritical, length(splitt)+1))
 if(intercept) {offset <- coef[1]; coef <- coef[-1] }
 k <- 1
 for(i in 1:(length(criticalpoints)- 1)){
   lower <- criticalpoints[i]+1
   upper <- criticalpoints[i+1]
   if(lower == upper) next
   term <- paste(splitt[lower:(upper-1)], collapse="")
   containsmz <- length(grep("mz", term)) > 0
   if(!containsmz) next
   termcoef <- paste(coef[k], term, sep = "*")
   if(k == 1) functionbodystring <- termcoef
   else functionbodystring <- paste(functionbodystring, termcoef, sep = "+")
   k <- k + 1
 }
 if(exists("functionbodystring")){
 if(intercept)
  functionbodystring <- paste(functionbodystring, offset, sep="+")
  foo <- function(mz){}
  body(foo) <- parse(text = functionbodystring) 
}
 else{
    if(!intercept) stop("Error occured in conversion from formula to function \n")
    foo <- function(mz){}
    body(foo) <-  eval(substitute(expression(rep(estimate, length(mz))), list(estimate = offset))) 
 }

  return(foo)
}


### 32. simplepeakdetect

simplepeakdetect <- function(spectrum, window, threshold){
  l <- length(spectrum[,1])
  delta <- diff(spectrum[,2])/diff(spectrum[,1])
  Cres <- .C("peakdetect", spectrum = spectrum[,2], delta = delta,
             l = as.integer(l), window = as.integer(window),
             threshold = threshold, peakind = integer(l),
             counter = as.integer(0))
  if(Cres$counter == 0) stop("No peaks found \n")
  else{
     peakind <- Cres$peakind[1:Cres$counter]+1 
     bestpeak <- peakind[which.max(spectrum[peakind,2])] 
     xpeak <- spectrum[bestpeak,1]  
     xenvelopinit <- xpeak - 1 
     xenvelopend <- xpeak + 1
     indinit <- which(spectrum[,1] >= xenvelopinit & spectrum[,1] < xpeak)
     indend <-  which(spectrum[,1] > xpeak  & spectrum[,1] <= xenvelopend)
     diffinit <- rev(delta[indinit]) 
     diffend <-  delta[indend - 1] 
     xxind <- bestpeak
      ileft <- 1
        while(diffinit[ileft] > 0){
          xxind <- c(bestpeak - ileft, xxind)
          ileft <- ileft +1
          if(ileft > length(diffinit)) break
         }     
     iright <- 1 
      while(diffend[iright] < 0){
          xxind <- c(xxind, bestpeak + iright)
          iright <- iright + 1
          if(iright > length(diffend)) break
         }
       return(cbind(x = spectrum[xxind,1], y = spectrum[xxind,2]))
    }
}

################################################################################

nnladlogbarrier <- function(response, betastart, trace = FALSE, alpha = 0.01, gammastart = 10, gammamax = 10^15, gammamult = 20, eps = 1e-6){
  betaold <- betastart
  fitted <- drop( get("basis", parent.frame(1)) %*% betaold)
  rold <- abs(response - fitted) + betastart[1]
  gamma <- gammastart
  obj <- Inf  
  while(gamma < gammamax){
  converged <- FALSE 
  if(trace) cat("gamma: ", gamma, "\n")
  while(!converged){
  if(any(betaold < 0)) stop("'beta' has to be non-negative \n")
  if(any(rold < 0)) stop("'rold' has to be non-negative \n")
  fitted <- drop( get("basis", parent.frame(1)) %*% betaold)
  res <- fitted - response
  Xiplus  <-  1/(res  + rold)
  Ximinus <-  1/(-res  + rold)
  grad.r <- 1 - Xiplus/gamma - Ximinus/gamma  - 1/(gamma * rold)  
  grad.beta <- -drop((Xiplus - Ximinus) %*% get("basis", parent.frame(1)))/gamma - 1/(gamma * betaold)
  Ainv <- 1/(Xiplus^2 + Ximinus^2 + rold^(-2))
  D.x <- pmax( (Xiplus^2 + Ximinus^2 - (Ainv * (Xiplus^2 - Ximinus^2)^2)), .Machine$double.eps)
  R <- .C("multiply", colptr = as.integer(get("G", parent.frame(1))@p),
                       rowind = as.integer(get("G", parent.frame(1))@i),
                       val = as.double(get("G", parent.frame(1))@x),
                       colptrphi = as.integer(get("basis", parent.frame(1))@p),
                       rowindphi = as.integer(get("basis", parent.frame(1))@i),
                       valphi = as.double(get("basis", parent.frame(1))@x),
                       vec = as.double(D.x / gamma),
                       p = as.integer(ncol(get("G", parent.frame(1)))))
  R <- .C("addiagonal", colptr = as.integer(R$colptr), rowind = as.integer(R$rowind), val = as.double(R$val), add = as.double( (1/betaold^2 )/ gamma), p = as.integer(ncol(get("basis", parent.frame(1)))))
  R <- new("dsCMatrix", p = R$colptr, i = R$rowind, x = R$val, uplo = "U", Dim = get("G", parent.frame(1))@Dim, factors = list()) 
  R <- try(Cholesky(R), silent = TRUE)
  if(inherits(R, "try-error")){
     warning("Error in nonnegative least absolute deviation estimation: Hessian not positive definite.  \n")
     gamma <- gammamax
     break    
   } 
  rhs.beta <- (-grad.beta  + drop (((Xiplus^2 - Ximinus^2) * Ainv * grad.r) %*% get("basis", parent.frame(1))))
  dir.beta <- solve(R, rhs.beta)@x
  dir.r <- (gamma * (grad.r) + (Xiplus^2 - Ximinus^2) * drop(get("basis", parent.frame(1)) %*% dir.beta)) * (-Ainv)
  stepsize <- try(linesearchlad(response, alpha), silent = TRUE)
   if(inherits(stepsize, "try-error")){
     warning("Error in nonnegative least absolute deviation estimation: stepsize selection failed \n")
     gamma <- gammamax
     break    
   }
   if(trace) cat("stepsize: ", stepsize, "\n")
   beta <- betaold + dir.beta * stepsize
   r <- rold + dir.r * stepsize
   conv <- abs(sum(grad.beta * drop(dir.beta)) + sum(grad.r * drop(dir.r)))
   if(trace) cat("convergence inner loop:", conv, "\n")
   if((stepsize == 1) | conv < eps){
        converged <- TRUE
        gamma <- gamma * gammamult
      }
   betaold <- drop(beta)
   rold <- drop(r)
}
 }
  if(!exists("r")) r <- rold
  if(is.function(beta)) beta <- betaold
 return(list(beta = drop(beta), r = drop(r), gradnorm = sqrt(sum(grad.beta * grad.beta) + sum(grad.r * grad.r)), obj = obj))
}

### linesearch [for nonnegative least absolute deviation]

linesearchlad <- function(response, alpha){
  stepsize <- 1
  sigma <- 0.95
  sigmavec <- sigma^(0:500)
  whichsmaller0.beta <- get("dir.beta", parent.frame(1)) < 0
  if(sum(whichsmaller0.beta) > 0)
    maxstepsize.beta <- min(- (get("betaold", parent.frame(1)) / get("dir.beta", parent.frame(1)))[whichsmaller0.beta])
  else maxstepsize.beta <- 1
  if(maxstepsize.beta < 1)
      maxstepsize.beta <- max(sigmavec[sigmavec < maxstepsize.beta])
  
  whichsmaller0.r <- get("dir.r", parent.frame(1)) < 0
  if(sum(whichsmaller0.r) > 0)
    maxstepsize.r <- min(- (get("rold", parent.frame(1)) / get("dir.r", parent.frame(1)))[whichsmaller0.r])
  else maxstepsize.r <- 1
  if(maxstepsize.r < 1)
      maxstepsize.r <- max(sigmavec[sigmavec < maxstepsize.r])

  stepsize <- min(c(1, maxstepsize.beta, maxstepsize.r)) 
  finished <- FALSE
  lossold <- sum(get("rold", parent.frame(1)))
  barrierold <- -(1/get("gamma", parent.frame(1))) * (sum(log(get("betaold", parent.frame(1)))) +
                   - sum(log(get("Xiplus", parent.frame(1)))) - sum(log(get("Ximinus", parent.frame(1)))) +
                                                      sum(log(get("rold", parent.frame(1)))))
  objectiveold <- lossold + barrierold
  while(!finished){
    if(get("trace", parent.frame(1))) cat("linesearch...", stepsize, "\n")
    beta <-  get("betaold", parent.frame(1)) + stepsize * get("dir.beta", parent.frame(1))
    r <- get("rold", parent.frame(1)) + stepsize * get("dir.r", parent.frame(1))
    
    loss <-  sum(r)

    fitted <- drop( get("basis", parent.frame(2)) %*% beta) 
    barrier <- -(1/get("gamma", parent.frame(1))) * (sum(log(fitted - response + r)) + sum(log(response - fitted + r)) + sum(log(beta)) + sum(log(r)))
    objective <- loss + barrier
    if(is.nan(objective)) objective <- Inf
    rhs <- objectiveold + alpha * stepsize * (sum(get("grad.beta", parent.frame(1)) * get("dir.beta", parent.frame(1))) + sum(get("grad.r", parent.frame(1)) * get("dir.r", parent.frame(1)))) 
    if(  objective >= rhs){  
     stepsize <- sigma * stepsize
     if(stepsize < (0.5 * .Machine$double.eps)){
       stop("Canceled stepsize selection. \n")
     }
     }  
    else finished <- TRUE
  }
  if(get("trace", parent.frame(1))){ cat("loss:", loss, "\t") ; cat("logbarrier:", barrier, "\n")}
  assign("obj", loss, parent.frame(1))
  stepsize
}

###


calculatedenoisingbasis.gaussian <- function(x, sigma, eps = 1e-05, uppernonzero = 10^6){ ### find better estimate for uppernonzero.

  n <- length(x) 
  require(Matrix)
  #
  sigma <- sigma(x)
  #dyn.load("getblocksparseold.so")
  block <- matrix(0, nrow = uppernonzero, ncol = 3)
  resultC <- .C("getblocksparse", block = block, x = as.double(x),
                               n = as.integer(n),
                               sigma = as.double(sigma),
                 eps = as.double(eps), uppernonzero = as.integer(uppernonzero),
                 nonzero = as.integer(-1))
  return(forceSymmetric(
                     sparseMatrix(i = as.integer(resultC$block[1:(resultC$nonzero + 1), 1]),
                                  j = as.integer(resultC$block[1:(resultC$nonzero + 1), 2]),
                                  x = resultC$block[1:(resultC$nonzero + 1), 3],
                                  dims = c(n, n))))}


calculatedenoisingbasis.emg <- function(x, alpha, sigma, mu, eps = 1e-05, uppernonzero = 10^6){ ### find better estimate for uppernonzero.
  n <- length(x) 
  require(Matrix)
  # Determine modes and so on
  alpha <- alpha(x)
  sigma <- sigma(x)
  mu <- mu(x)
  moderes <- apply(cbind(alpha, sigma, mu), 1, function(z) {getmode <- determinemode(z[1], z[2], z[3]); c(-getmode$mode + z[3], 1/getmode$val)})
  # to be passed to the C function
  mu <- moderes[1,] 
  scale <- moderes[2,] 
  #dyn.load("getblocksparseold.so")
  block <- matrix(0, nrow = uppernonzero, ncol = 3)
  resultC <- .C("getblocksparseemg", block = block, x = as.double(x),
                               n = as.integer(n),
                               alpha = as.double(alpha),
                               sigma = as.double(sigma),
                               mu = as.double(mu),
                               scale = as.double(scale),
                 eps = as.double(eps),
                 uppernonzero = as.integer(uppernonzero),
                 nonzero = as.integer(-1))
                 return(sparseMatrix(i = as.integer(resultC$block[1:(resultC$nonzero + 1), 1]),
                              j = as.integer(resultC$block[1:(resultC$nonzero + 1), 2]),
                              x = resultC$block[1:(resultC$nonzero + 1), 3],
                              dims = c(n ,n)))}


####

base64decode <- function(z, what, size = NA, signed = TRUE, endian = .Platform$endian) 
{

     require(bitops)
    if (!is.character(z))                                               
        stop("base64decode: Input argument 'z' is suppose to be a string")
    if (length(z) == 1)                                                   
        z = strsplit(z, NULL)[[1]]
    if (length(z)%%4 != 0)
        warning("In base64decode: Length of base64 data (z) not a multiple of 4.")
    alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
    alpha = strsplit(alpha, NULL)[[1]]
    y = match(z, alpha, nomatch = -1) - 1
    if (any(y == -1))
        stop("base64decode: Input string is not in Base64 format")
    if (any(y == 64))
        y = y[y != 64]
    neByte = length(y)
    nBlock = ceiling(neByte/4)
    ndByte = 3 * nBlock
    if (neByte < 4 * nBlock)
        y[(neByte + 1):(4 * nBlock)] = 0
    dim(y) = c(4, nBlock)
    x = matrix(as.integer(0), 3, nBlock)
    x[1, ] = bitOr(bitShiftL(y[1, ], 2), bitShiftR(y[2, ], 4))
    x[2, ] = bitOr(bitShiftL(y[2, ], 4), bitShiftR(y[3, ], 2))
    x[3, ] = bitOr(bitShiftL(y[3, ], 6), y[4, ])
    x = bitAnd(x, 255)
    if (neByte%%4 == 2)
        x = x[1:(ndByte - 2)]
    if (neByte%%4 == 3)
        x = x[1:(ndByte - 1)]
    r = as.raw(x)
    TypeList = c("logical", "integer", "double", "complex", "character",
        "raw", "numeric", "int")
    if (!is.character(what) || length(what) != 1 || !(what %in%
        TypeList))
        what <- typeof(what)
    if (what == "raw")
        return(r)
    if (is.na(size))
        size = switch(match(what, TypeList), 4, 4, 8, 16, 2,
            1, 8, 4)
    n = length(r)
    if (n%%size)
        stop("raw2bin: number of elements in 'r' is not multiple of 'size'")
    x = readBin(r, what, n = n%/%size, size = size, signed = signed,
        endian = endian)
    if (what == "character")
        x = paste(x, collapse = "")
    return(x)
}
