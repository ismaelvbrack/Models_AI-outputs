
model{
  
  # Priors
  # Psi
  alpha0 ~ dnorm(0,0.01)
  alpha1 ~ dnorm(0,0.01)
  
  # detection prob.
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  
  # Likelihood
  for(s in 1:nsites){
    z[s] ~ dbern(psi[s]) # latent occupancy state
    
    # occupancy prob.
    logit(psi[s]) <- alpha0 + alpha1 * forest[s]
    
    for(j in 1:ndays){
      # detection prob.
      logit(p[s,j]) <- beta0 + beta1 * twi[s]
      
      Y[s,j] ~ dbern(z[s]*p[s,j])
      
    } #j
  } #s
  
}