


model{
  
  # Priors ------------------------------------------------------------------ 
  
  pa ~ dunif(0,1) # algorithm detection prob.
  ph ~ dunif(0,1) # human detection prob.
  
  # Algorithm precision (TP prop.) - random effects
  mu.prec ~ dnorm(0, 0.01)
  sd.prec ~ dunif(0,5)
  tau.prec <- 1 / sd.prec^2
  
  # Scores
  mu.scr[1] ~ dnorm(0, 0.01) # mean score values for class=0
  mu.dif ~ dunif(0,5)     # impose restriction: mu2 > mu1
  mu.scr[2] <- mu.scr[1] + mu.dif  # mean score value for class=1
  
  sd.scr[1] ~ dunif(0, 10) # sd for scores
  sd.scr[2] ~ dunif(0, 10) # sd for scores
  tau.scr[1] <- 1 / sd.scr[1]^2
  tau.scr[2] <- 1 / sd.scr[2]^2
  
  # Likelihood ---------------------------------------------------------------
  
  for(s in 1:nsites){
    # Define multinomial cell probabilities
    pi[s,1] <- pa*ph       # 11
    pi[s,2] <- (1-pa)*ph   # 10 
    pi[s,3] <- pa*(1-ph)   # 01
    
    # Conditional multinomial probs.
    for(k in 1:3){
      pic[s,k] <- pi[s,k] / sum(pi[s,])
    }
    
    pcap[s] <-  1 - (1-pa)*(1-ph) # overall detection prob.
    
    Y[s,] ~ dmulti(pic[s,], ndet1[s]) # observarion of encounter histories
    
    TP2[s,1] ~ dbin(precision[s], ndet2[s,1])
    
    #* Derived abundance estimates
    N1[s] <- ndet1[s] / pcap[s]
    N2[s] <- TP2[s,1] / pa
    N[s,1] <- N1[s] + N2[s]
    
    # Algorithm precision random effects
    alpha[s] ~ dnorm(mu.prec,tau.prec)
    logit(precision[s]) <- alpha[s]
    
    for(j in 2:nflight[s]){
      # True positives in the unverified counts
      TP2[s,j] ~ dbin(precision[s], ndet2[s,j])
      
      #* Derived abundance estimates
      N[s,j] <- TP2[s,j] / pa
      
    } #j
    
  } #s 
  
  #* Scores
  for(i in 1:nobj){
    class[i] ~ dbern(precision[site[i]]) # true latent class (TP of FP)
    
    score[i] ~ dnorm(mu.scr[(class[i]+1)], tau.scr[(class[i]+1)]) # object-level scores
  } #i
  
}



