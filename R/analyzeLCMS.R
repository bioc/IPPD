analyzeLCMS <- function(data, arglist.fitModelParameters = list(),
                              arglist.getPeaklist = list(),
                        arglist.threshold = list(),
                              arglist.sweepline = list(),
                              trace = TRUE) {
  
  ###------------------------------------------------------------------------------------------####
### PART 0: parse input data.
  ###------------------------------------------------------------------------------------------####
  if(trace){
    cat("Parsing input data...\n")
  }
  if(is.character(data)){
    # get file suffix
    strparts <- unlist(strsplit(data, ".", fixed = TRUE))
    if(length(strparts) != 2)
      stop("Argument 'data' has to be of the form 'filename.suffix' \n")
    if(strparts[2] == "mzXML"){
      data <- read.mzXML(data)
      format <- "mzXML"      
    }
    else{# invoke read.table() function
      data <- read.table(data)
      if(ncol(data) != 3){
        stop("Argument 'data' does not have a proper format; it
              has to be either an mzXML file or
              a three-column matrix \n")
      }
      data <- try(apply(data, 2, as.numeric), silent = TRUE)
      if(inherits(data, "try-error")){
         stop("Argument 'data' does not have a proper format; it
              has to be either an mzXML file or
              a three-column matrix \n")
      }
      format <- "txt"
   }   
  }
  
  else{
    if(class(data) == "mzXML"){
      format <- "mzXML"
    }
    else{
     if(is.matrix(data) || is.data.frame(data) ){
       if(ncol(data) != 3){
        stop("Argument 'data' does not have a proper format; it
              has to be either an mzXML file or
              a three-column matrix \n")  
       }
     }
     else{
       stop("Object 'data' has be either an object of class 'mzXML'
             or three-column matrix \n")
     }
     data <- try(apply(data, 2, as.numeric), silent = TRUE)
     if(inherits(data, "try-error")){
         stop("Argument 'data' does not have a proper format; it
              has to be either an mzXML file or
              a three-column matrix \n")
      }
      format <- "txt"
   }   
  }

  if(format == "mzXML"){
 
    scans <- data$scan
    msone <- unlist( lapply(scans, function(z) z$msLevel) )
    if( sum(msone == 1) == 0)
      stop("No valid spectra (msLevel == 1) found \n")
    else
      scans <- scans[msone == 1]

    nscans <- length(scans)

    getRT <- function(scans){
      rtlist <- lapply(scans, function(x)
                   as.numeric(sub("([^0-9]*)([0-9|.]+)([^0-9]*)", "\\2", x$scanAttr)))
      unlist(rtlist)
    }

    rt <- getRT(scans)
    perm <- order(rt) # order scans by increasing retention times
  }
  else{
    rt <- unique(data[,1]) ## retention time always in first column.
    sortrt <- sort(rt)
    nscans <- length(rt)
  }

  ###------------------------------------------------------------------------------------------####
  ### PART 1: Fit model parameters (peak shape)
  ###------------------------------------------------------------------------------------------####

  if(trace){
    cat("Fit model parameters...\n")
  } 
  
    fitmodelparameters <- is.null(arglist.getPeaklist$model.parameters)
    if(fitmodelparameters){

      l_fitpar <- vector(mode = "list", length = nscans)

      # set parameters (if not pre-specified)
      if(is.null(arglist.fitModelParameters$control$window))
        arglist.fitModelParameters$control$window = 5

      if(is.null(arglist.fitModelParameters$model))
        arglist.fitModelParameters$model = "Gaussian"

      if(is.null(arglist.fitModelParameters$formula.alpha))
        arglist.fitModelParameters$formula.alpha = formula(~1)

      if(is.null(arglist.fitModelParameters$formula.sigma))
        arglist.fitModelParameters$formula.sigma = formula(~1)

      if(is.null(arglist.fitModelParameters$formula.mu))
        arglist.fitModelParameters$formula.mu = formula(~1)
      #

    # fit peak-shape parameters
    for(i in 1:nscans){

      if(format == "mzXML"){
        ix <- perm[i]
        yi <- scans[[ix]]$peaks
        xi <- scans[[ix]]$mass
      }
      else{
        rti <- data[,1] == sortrt[i]
        xi <- data[rti, 2]
        yi <- data[rti, 3]
      }  

      arglist.fitModelParameters$mz <- xi
      arglist.fitModelParameters$intensities <- yi


      ### set threshold scan-wise
      if(is.null(arglist.fitModelParameters$control$threshold))
        arglist.fitModelParameters$control$threshold = median(yi)
      ###
      l_fitpar[[i]] <- try(do.call(fitModelParameters,
                                   args = arglist.fitModelParameters),
          silent = TRUE)
    }
    # check for errors
    error_fitpar <- unlist(lapply(l_fitpar, function(z) inherits(z, "try-error")))
    if(all(error_fitpar))
      stop("Model parameters could not be estimated. Try to (re)-set arglist.fitModelParameters \n")

      which_noerror_fitpar <- which(!error_fitpar)
   }
    

   ###------------------------------------------------------------------------------------------####
   ### PART 2: process scan by scan.
   ###------------------------------------------------------------------------------------------####

   if(trace){
    cat("Process scans...\n")
   }
  
   peaklists <- vector(mode = "list", length = nscans)

   # set parameters (if not pre-specified)
   if(is.null(arglist.getPeaklist$model)) 
     arglist.getPeaklist$model <- "Gaussian"

    if(is.null(arglist.getPeaklist$loss)) 
     arglist.getPeaklist$loss <- "L2"

    if(is.null(arglist.getPeaklist$binning)) 
      arglist.getPeaklist$binning <- FALSE

    if(is.null(arglist.getPeaklist$postprocessing)) 
      arglist.getPeaklist$postprocessing <- TRUE

    if(is.null(arglist.getPeaklist$trace)) 
      arglist.getPeaklist$trace <- FALSE 

    if(is.null(arglist.getPeaklist$returnbasis)) 
      arglist.getPeaklist$returnbasis <- FALSE 
   #

    for(i in 1:nscans){
      if(trace){
        cat("Processing scan ", paste(i, sep = ""), "...\n")
      }
      
      if(fitmodelparameters){ 
        if(error_fitpar[i]){ # find nearest neighbour
          nei <- which_noerror_fitpar[ which.min(abs(which_noerror_fitpar  - i)) ]
          arglist.getPeaklist$model.parameters <- l_fitpar[[nei]]
        }
        else
          arglist.getPeaklist$model.parameters <- l_fitpar[[i]]         
      }

      if(format == "mzXML"){
        ix <- perm[i]
        arglist.getPeaklist$mz <- scans[[ix]]$mass 
        arglist.getPeaklist$intensities <- scans[[ix]]$peaks
      }
      else{
         rti <- data[,1] == sortrt[i]
         arglist.getPeaklist$mz <- data[rti, 2]
         arglist.getPeaklist$intensities <- data[rti, 3]
      }  
      
      
      # get peak list
      peaklists[[i]] <- try(do.call("getPeaklist", args = arglist.getPeaklist),
                            silent = TRUE)
      # 
    }

   ###------------------------------------------------------------------------------------------####
   ### PART 3: threshold results.
   ###------------------------------------------------------------------------------------------####

   # set parameters

   if(is.null(arglist.threshold$threshold))
     arglist.threshold$threshold <- 3
   if(is.null(arglist.threshold$ratio))
      arglist.threshold$ratio <- "ratio"
   if(is.null(arglist.threshold$eps))
     arglist.threshold$eps <- 1e-05
   if(is.null(arglist.threshold$refit))
     arglist.threshold$refit <- FALSE
  
   peakliststhresholded <- lapply(peaklists, function(z){ if(inherits(z, "try-error")){ return(NULL)}
                                                             else{
                                                              threshold(z, threshold = arglist.threshold$threshold,
                                  ratio = arglist.threshold$ratio,
                                  refit = arglist.threshold$refit,
                                  eps = arglist.threshold$eps)}}) 
   
   ###------------------------------------------------------------------------------------------####
   ### PART 4: Agglomerate results using sweep line.
   ###------------------------------------------------------------------------------------------####

  if(trace){
    cat("Agglomerate results of scans using sweepline...\n")
   }


  if(is.null(arglist.sweepline$tol))
    arglist.sweepline$tol <-100
  if(is.null(arglist.sweepline$gap))
    arglist.sweepline$gap <- 2
  if(is.null(arglist.sweepline$minboxlength))
    arglist.sweepline$minboxlength <- 5
  
  boxes <- sweepline(peakliststhresholded, rt, tol = arglist.sweepline$tol, gap = arglist.sweepline$gap,
                     minboxlength = arglist.sweepline$minboxlength)

  return(list(peaklists = peakliststhresholded, boxes = boxes))

  ### END of function ###
  
}
