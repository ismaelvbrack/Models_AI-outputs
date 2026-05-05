
simul_counts <- function(
    nsites = 15, # number of sites
    
    N = round(runif(15, 50,300)), # colony abundance
    
    pa = runif(15, 0.6,0.9), # algorithm detection probability
    ph = 0.7, # human detection probability
    
    precision = runif(15, 0.4,0.8),
    
    prop2obs = 1/10, # proportion of the area with verification
    
    mu.scores = c(qlogis(0.5), qlogis(.75)), # mean score for FPs and TPs

    sd.scr = c(0.5, 0.5) # sd for FP and TP scores
){
  
  stopifnot(all.equal(nsites, length(N), length(pa), length(precision)) )
  
  stopifnot(mu.scores[2] > mu.scores[1] )
  
  # Simulate Counts ----------------------------------------------------------------
  
  N1 = round(prop2obs * N) # true abundance in the verification subset
  N2 = N - N1 #  true abundance in the remaining site
  
  #* Verified counts
  Y <- array(dim=c(nsites,4),dimnames=list(1:nsites,c("11","10","01","00")))
  
  for(s in 1:nsites){
    obs.human <- rbinom(N1[s],1,ph) # human detections for each individual in N1
    tp.algo <- rbinom(N1[s],1,pa)  # algorithm detections for each individual in N1
    
    temp.c <- table(factor(paste0(obs.human,tp.algo),  # combine detections into counts
                           levels=colnames(Y))
    )
    
    Y[s,] <- temp.c
    
  } #s
  
  # Algorithm true positives
  TP1 <- rowSums(Y[,c("11","01")])
  TP2 <- rbinom(nsites,N2,pa)
  
  # Algorithm false positives
  FP1 <- rpois(nsites, round(TP1*(1-precision) / precision))
  FP2 <- rpois(nsites, round(TP2*(1-precision) / precision))
  
  #cbind(N, N1, N2, Y, TP1, TP2, FP1, FP2)
  
  counts <- cbind(Y[,-4],FP1, unverified= TP2 + FP2)
  
  # Simulate scores data -----------------------------------------------------------------------------------------------------
  scr.data <- list()
  
  oo=0
  for(i in 1:nsites){
    oo=oo+1
    # True/False detection
    tab <- data.frame(
      class=c(rep(1,TP1[i] + TP2[i]) ,rep(0,FP1[i] + FP2[i])),
      obs.class=c(rep(c(1,NA), c(TP1[i],TP2[i])),
                  rep(c(0,NA), c(FP1[i],FP2[i]))),
      score=NA,
      site=i
    )
    # Scores for each detected object
    tab[which(tab$class==1),"score"] <- rnorm(TP1[i] + TP2[i], mu.scores[2], sd.scr[2]) # TP scores
    
    tab[which(tab$class==0),"score"] <- rnorm(FP1[i] + FP2[i], mu.scores[1], sd.scr[1]) # FP scores
    
    scr.data[[oo]] <- tab
    
  }#i
  
  scr.data <- do.call(rbind, scr.data)
  
  return(list(
    true = list(
      N = cbind(N, N1, N2), 
      det = cbind(pa, ph, precision),
      scores= rbind(mu.scores, sd.scr)
    ),
    counts = counts,
    scr.data = scr.data
  ))
  
} # function
