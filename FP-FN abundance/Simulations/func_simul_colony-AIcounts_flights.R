

simul_counts <- function(
    nsites = 20, # number of sites
    nflights = 8, #number of flights per site
    
    N = t(sapply(round(runif(20, 40,200)), rep, each=8)), # colony abundance
    
    pa = runif(20, 0.6,0.6), # algorithm detection probability
    ph = 0.7, # human detection probability
    
    mean.prec = runif(20, 0.2,0.8), # mean precision for each site
    sd.prec = 0.4, # std. dev of precision within flights of the same site
    
    prop2obs = 2/10, # proportion of the area with verification
    
    mu.scores = c(qlogis(0.5), qlogis(.75)), # mean score for FPs and TPs
    
    sd.scr = c(0.5, 0.5) # sd for FP and TP scores
){
  
  stopifnot(all.equal(nsites, nrow(N), length(pa), length(mean.prec)) )
  
  stopifnot(all.equal(nflights, ncol(N)))
  
  stopifnot(mu.scores[2] > mu.scores[1] )
  
  # Simulate Counts ----------------------------------------------------------------
  
  N1 = TP1 = TP2 = FP1 = FP2 = 
    matrix(NA, nsites, nflights)
  #* Verified flight
  for(j in 1:nflights){
    N1[,j] = rbinom(nsites,N[,j],prop2obs)
  }
  N2 = N - N1 #  true abundance in the remaining site
  
  #* Verified counts
  Y <- array(dim=c(nsites,4,nflights),dimnames=list(1:nsites,c("11","10","01","00"),1:nflights) )
  
  for(s in 1:nsites){
    for(j in 1:nflights){
      obs.human <- rbinom(N1[s,j],1,ph) # human detections for each individual in N1
      tp.algo <- rbinom(N1[s,j],1,pa)  # algorithm detections for each individual in N1
      
      temp.c <- table(factor(paste0(obs.human,tp.algo),  # combine detections into counts
                             levels=colnames(Y))
      )
      
      Y[s,,j] <- temp.c
    } #j

  } #s
  
  #* Algorithm precision
  mu.prec <- qlogis(mean.prec)
  
  lprec <- t(sapply(mu.prec, function(x) rnorm(nflights,x,sd.prec)))
  
  precision <- plogis(lprec)
  
  counts <- array(dim=c(nsites,5,nflights),
                  dimnames=list(1:nsites,c("11","10","01","fp","unverified"),1:nflights) )
  for(j in 1:nflights){
    # Algorithm true positives
    TP1[,j] <- rowSums(Y[,c("11","01"),j])
    TP2[,j] <- rbinom(nsites,N2[,j],pa)
    
    # Algorithm false positives
    FP1[,j] <- rpois(nsites, round(TP1[,j] * (1-precision[,j]) / precision[,j]))
    FP2[,j] <- rpois(nsites, round(TP2[,j] * (1-precision[,j]) / precision[,j]))
    
    counts[,,j] <- cbind(Y[,-4,j],fp=FP1[,j], unverified= TP2[,j] + FP2[,j])
  }
  
  #cbind(N[,1], N1[,1], N2[,1], Y[,,1], TP1[,1], TP2[,1], FP1[,1], FP2[,1])
  
  # Simulate scores data -----------------------------------------------------------------------------------------------------
  scr.data <- list()
  
  oo=0
  for(i in 1:nsites){
    for(j in 1:nflights){
      oo=oo+1
      # True/False detection
      tab <- data.frame(
        class=c(rep(1,TP1[i,j] + TP2[i,j]) ,rep(0,FP1[i,j] + FP2[i,j])),
        obs.class=c(rep(c(1,NA), c(TP1[i,j],TP2[i,j])),
                    rep(c(0,NA), c(FP1[i,j],FP2[i,j]))),
        score=NA,
        site=i,
        flight=j
      )
      # Scores for each detected object
      tab[which(tab$class==1),"score"] <- rnorm(TP1[i,j] + TP2[i,j], mu.scores[2], sd.scr[2]) # TP scores
      
      tab[which(tab$class==0),"score"] <- rnorm(FP1[i,j] + FP2[i,j], mu.scores[1], sd.scr[1]) # FP scores
      
      scr.data[[oo]] <- tab
    } #j
  }#i
  
  scr.data <- do.call(rbind, scr.data)
  
  return(list(
    true = list(
      abund = list(N = N, N1 = N1, N2 = N2), 
      det = list(pa = pa, ph = ph, mean.prec = mean.prec, sd.prec = sd.prec, precision = precision),
      scores= rbind(mu.scores, sd.scr)
    ),
    counts = counts,
    scr.data = scr.data
  ))
  
} # function
