

model{
  
  # Priors
  psi ~ dunif(0,1) # occupancy prob.
  p ~ dunif(0,1) # detection prob.
  
  # Likelihood
  for(s in 1:nsites){
    z[s] ~ dbern(psi) # latent occupancy state
    
    for(j in 1:ndays){
      Y[s,j] ~ dbern(z[s]*p) # species presence in the image
      
    } #j
  } #s

}