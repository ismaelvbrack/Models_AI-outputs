
library(jagsUI)

# Simulate data
#source(file.path("Occupancy AI objects","Simulations","simul_occuAI_3-level_normal scores.R"))

source(file.path("Simulations","func_simul_occuFPFN.R"))
tab <- simul_occuFPFN(
  nsites = 30, 
  nvisits = 10,
  mean.nimgs = 8,
  prop.verif = 0.1,
  
  psi = 0.2,
  phi = 0.2,
  lam = 1,
  ome = 0.2,
  
  mu.scores = qlogis(c(0.5,0.8)),
  sd.scr = rep(0.5,2)
)

table(tab$obj.data$class.obs)
# Including verified data -------------------------------------------------

# Bundle data for JAGS
dat <- list(
  # Sampling
  nsites = nrow(tab$nimgs), 
  nvisits = ncol(tab$nimgs),
  nimgs = tab$nimgs,
  
  ndets=tab$ndets,
  
  # Object-level information
  class=tab$obj.data$class.obs,
  score=tab$obj.data$score,
  site=tab$obj.data$site,
  visit=tab$obj.data$visit,
  image=tab$obj.data$image,
  nobj=nrow(tab$obj.data)
)

# Initial values
inits <- function() list(
  z=rep(1,dat$nsites),
  u=matrix(1,dat$nsites,dat$nvisits),
  class=ifelse(is.na(tab$obj.data$class.obs),1,NA)
)

# Parameters monitored
params <- c("psi","phi","lam","ome","mu")

# MCMC settings
ni <- 20000; nt <- 2; nb <- 10000; nc <- 3; na <- 500

out2 <- jags(dat, inits, params,
             model.file=file.path("Simulations","occu_model_3level_normal scores.R"),
             n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
             parallel=T)

beepr::beep(2)

plot(out2)

MCMCvis::MCMCtrace(out3,
                   gvals=c(psi,phi,lam,ome,fp.mu,tp.mu,0),
                   Rhat=T,n.eff=T,pdf=F)

few.samps <- sample(1:out3$mcmc.info$n.samples, 10000)
pairs(cbind(psi=out3$sims.list$psi[few.samps],
            phi=out3$sims.list$phi[few.samps],
            lambda=out3$sims.list$lam[few.samps],
            omega=out3$sims.list$ome[few.samps],
            mu=out3$sims.list$mu[few.samps,]
            ),panel=panel.smooth,
      col=ggplot2::alpha("gray60",0.7))

# see posteriors with true values
par(mfrow=c(2,2))
plot(density(out3$sims.list$psi),lwd=2,col="cyan3",main="Occupancy prob.")
lines(density(out2$sims.list$psi),lwd=2,col="blue")
lines(density(out3$sims.list$psi),lwd=2,col="darkblue")
abline(v=psi,col="red",lwd=2)

plot(density(out3$sims.list$phi),lwd=2,col=4,main="Availability prob.")
abline(v=phi,col=2)
plot(density(out3$sims.list$lam),lwd=2,col=4,main="TP rate")
abline(v=lam,col=2)
plot(density(out3$sims.list$ome),lwd=2,col=4,main="FP rate")
abline(v=ome,col=2)


# Without verified data ---------------------------------------------------

# # Bundle data
# dat <- list(
#   score=scr.data$score,
#   site=scr.data$site,
#   image=scr.data$image,
#   nsites=nsites,nimgs=nimgs,
#   nobj=nrow(scr.data),
#   ndets=ndets
# )
# 
# # Initial values
# inits <- function() list(
#   z=rep(1,nsites),
#   u=matrix(1,nsites,nimgs),
#   psi=runif(1),
#   lam=runif1),
#   ome=runif(1),
#   class=rep(1,nrow(scr.data))
# )
# 
# # Parameters monitored
# params <- c("psi","phi","theta","lam","ome","mu")
# 
# # MCMC settings
# ni <- 40000; nt <- 1; nb <- 20000; nc <- 3; na <- 5000
# 
# out2 <- jags(dat, inits, params,
#              model.file=file.path("occu_model_normal-scores.R"),
#              n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
#              parallel=T)
# 
# out2
