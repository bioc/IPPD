### pl is already passed as an input (i.e. the finished peaklists)

sweepline <- function(peaklists, rt, tol = 100, gap = 2, minboxlength = 5){
  
  rt_int <-  as.numeric(as.factor(rt)) 
  
  rt_int_uq <- unique(rt_int)
  
  nscans <- length(rt_int_uq)
  
  flag <- 1
  init <- 1
  while(flag & init <= nscans){
    pl_init <- peaklists[[init]] 
   if(length(pl_init) > 0){
    flag  <- 0
   }
     init <- init + 1  
  }

  clmnames <- c("loc", "charge", "quant", "rt_begin", "rt_end", "npeaks", "gapcount", "processed", "closed")
  
  if(flag){
    warning("All peaklists are empty \n")
    boxes <- matrix(nrow = 0, ncol = length(clmnames))
    colnames(boxes) <- clmnames
  }

  init <- init - 1
  boxes  <- cbind(pl_init[,"loc_init"],  pl_init[,"charge"], pl_init[,"quant"], rt_int_uq[init], rt_int_uq[init], 1, 0, 0, 0)
  colnames(boxes) <- clmnames

  init <- init + 1
  ### auxiliary matching functions

  ppm <- function(m1, m2)  abs((m1 - m2)/(( 0.5 *  m1 + 0.5 * m2)) * 10^6) 
  

  ###

  boxes_closed <- NULL
  
  for(ii in init:nscans){

    pl_i <- peaklists[[ii]]
    if(length(pl_i) > 0) 
       N <- nrow(pl_i)
    else
       N <- 0 
    rtcur <- rt_int_uq[ii] 

    closed <- boxes[,"closed"]
    
    # remove boxes
    boxes_closed <- rbind(boxes_closed, boxes[closed == 1,])
    #
    boxes <- boxes[!closed,, drop = FALSE]
    boxes[,"processed"] <- 0 ### initialize all boxes as 'not processed'

    if(N > 0){
    for(jj in 1:N){

      peak_jj <- pl_i[jj,]

      loc_jj <- peak_jj["loc_init"]
      ch_jj <- peak_jj["charge"]

      match <- FALSE
      cand_jj <- which((boxes[,"charge"] == ch_jj) & (boxes[,"rt_end"] < rtcur))
      #cand_jj <- which((boxes[,"charge"] == ch_jj))

      if( length(cand_jj) > 0){
       match_ix <- cand_jj[which.min(abs(boxes[cand_jj,"loc"] - loc_jj))]

       ppm_jj <- ppm(boxes[match_ix, "loc"] * ch_jj, loc_jj * ch_jj)

       match <- ppm_jj < tol

      }

      if(match){
         npeaks <- boxes[match_ix, "npeaks"]
         wold <-  sum(boxes[match_ix, "quant"])
         wnew <-  peak_jj["quant"]
         sumw <- wold + wnew
         newloc <- boxes[match_ix, "loc"] * (wold/sumw)  + loc_jj * (wnew/sumw)
         newgapcount <- rtcur - boxes[match_ix, "rt_end"] - 1
         boxes[match_ix, "npeaks"] <-  npeaks + 1 
         boxes[match_ix, "loc"] <- newloc
         boxes[match_ix, "gapcount"] <- boxes[match_ix, "gapcount"] + newgapcount
         boxes[match_ix, "rt_end"] <- rtcur 
         boxes[match_ix, "processed"] <- 1
         boxes[match_ix, "quant"] <- boxes[match_ix, "quant"] + peak_jj["quant"]
      }
      else # open a new box
          boxes <- rbind(boxes,
                         c(peak_jj["loc_init"],  peak_jj["charge"], peak_jj["quant"], rtcur, rtcur, 1, 0, 1, 0))
    }
  }

    boxes[,"closed"] <- !boxes[,"processed"] & (rtcur - boxes[,"rt_end"]) > gap 

     
  }

  boxes <- rbind(boxes_closed, boxes)

  surv <-  boxes[,"npeaks"] >= minboxlength
  boxes <- boxes[surv,, drop = FALSE]
  ### convert rt back
  rtsort <- sort(unique(rt))
  boxes[,"rt_begin"] <-  rtsort[boxes[,"rt_begin"]]
  boxes[,"rt_end"] <- rtsort[boxes[,"rt_end"]]
  ###  
  
  oo <- order(boxes[,"loc"])
  
  

  return(boxes[oo,c("loc", "charge", "quant", "rt_begin", "rt_end", "npeaks", "gapcount"), drop = FALSE])
  

}
