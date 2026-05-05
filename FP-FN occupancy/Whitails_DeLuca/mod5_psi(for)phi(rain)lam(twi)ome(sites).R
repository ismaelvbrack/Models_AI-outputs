

model{
  
  ### psi(forest) phi(.) lam(twi) ome(.)
  
  # Priors
  
  # Psi
  alpha0 ~ dnorm(0,0.01)
  alpha1 ~ dnorm(0,0.01)
  
  # availability prob.
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  
  # Lambda
  delta0 ~ dnorm(0,0.01)
  delta1 ~ dnorm(0,0.01)
  
  # expected number of detected FP objects per image
  mu.ome ~ dnorm(0,0.01)
  sd.ome ~ dunif(0,10)
  tau.ome <- 1 / sd.ome^2
  
  # Scores
  mu[1] ~ dnorm(0, 0.01) # mean score values for class=0
  mu.dif ~ dunif(0,5)     # impose restriction: mu2 > mu1
  mu[2] <- mu[1] + mu.dif  # mean score value for class=1
  
  sd ~ dunif(0, 10) # sd for scores
  tau <- 1 / sd^2
  
  # Likelihood
  for(s in 1:nsites){
    
    logit(psi[s]) <- alpha0 + alpha1*forest[s]
    z[s] ~ dbern(psi[s]) # latent occupancy state
    
    log(lam[s]) <- delta0 + delta1*twi[s]
    
    # random effects for precision
    gamma[s] ~ dnorm(mu.ome, tau.ome)
    log(ome[s]) <- gamma[s]
    
    for(j in 1:nvisits){
      
      u[s,j] ~ dbern(z[s]*phi[s,j]) # species available in the occasion/day
      
      logit(phi[s,j]) <- beta0 + beta1*rain[s,j]
      
      for(i in 1:nimgs[s,j]){
        
        ndets[s,j,i] ~ dpois(lam[s]*u[s,j] + ome[s]) # number of objects detected per image
        
        tp.prop[s,j,i] <- lam[s]*u[s,j] / (lam[s]*u[s,j] + ome[s]) # proportion of true positives
      } #i
    } #j
  } #s
  
  # Scores given class
  for(k in 1:nobj){
    class[k] ~ dbern(tp.prop[site[k],visit[k],image[k]]) # true latent class (TP of FP)
    
    score[k] ~ dnorm(mu[(class[k]+1)], tau) # object-level scores
  } #k
  
}