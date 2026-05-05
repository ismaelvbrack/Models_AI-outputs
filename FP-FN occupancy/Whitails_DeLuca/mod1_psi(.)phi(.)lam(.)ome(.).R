

model{
  
  # Priors
  psi ~ dunif(0,1) # occupancy prob.
  phi ~ dunif(0,1) # availability prob.
  
  lam ~ dunif(0, 10) # expected number of detected TP objects per image
  ome ~ dunif(0, 10) # expected number of detected FP objects per image
  
  # Scores
  mu[1] ~ dnorm(0, 0.01) # mean score values for class=0
  mu.dif ~ dunif(0,5)     # impose restriction: mu2 > mu1
  mu[2] <- mu[1] + mu.dif  # mean score value for class=1
  
  sd ~ dunif(0, 10) # sd for scores
  tau <- 1 / sd^2
  
  # Likelihood
  for(s in 1:nsites){
    z[s] ~ dbern(psi) # latent occupancy state
    
    for(j in 1:nvisits){
      u[s,j] ~ dbern(z[s]*phi) # species available in the occasion/day
      
      for(i in 1:nimgs[s,j]){
        ndets[s,j,i] ~ dpois(lam*u[s,j] + ome) # number of objects detected per image
        tp.prop[s,j,i] <- lam*u[s,j] / (lam*u[s,j] + ome) # proportion of true positives
      } #i
    } #j
  } #s
  
  # Scores given class
  for(k in 1:nobj){
    class[k] ~ dbern(tp.prop[site[k],visit[k],image[k]]) # true latent class (TP of FP)
    
    score[k] ~ dnorm(mu[(class[k]+1)], tau) # object-level scores
  } #k
  
}