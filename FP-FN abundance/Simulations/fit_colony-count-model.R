
library(jagsUI)
library(MCMCvis)

source(file.path("Simulations","simul_colony-count-model_scores.R"))

# Bundle data
dat <- list(
  nsites=nsites,
  
  class=scr.data$obs.class,
  score=scr.data$score,
  site=scr.data$site,
  nobj=nrow(scr.data),
  
  Y=counts[,1:3],
  ndet1=rowSums(counts[,1:3]),
  #TP1=TP1,
  ndet2=counts[,"unverified"]
)

# Initial values
inits <- function() list(
  class=ifelse(is.na(scr.data$obs.class),1,NA)
)

# Parameters monitored
params <- c("pa","ph","mu.prec","mu.scr","precision","N","N1","N2")

# MCMC settings
ni <- 20000; nt <- 1; nb <- 5000; nc <- 3; na <- 1000

out1 <- jags(dat, inits, params,
             model.file=file.path("model1_colony-count.R"),
             n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
             parallel=T)



# See results -----------------------------------------------------------------------------------------------------

MCMCtrace(out1,
          params=c("pa","ph","mu.prec","mu.scr"),
          gvals=c(pa,ph,qlogis(mean(precision)),fp.mu,tp.mu),
          Rhat=T,n.eff=T,pdf=F)

rango <- range(c(N,out1$mean$N))


par(mfrow=c(1,3))

plot(N,out1$mean$N,pch=19,
     xlab="True N", ylab="Estimated N",
     xlim=rango,ylim=rango)
abline(a=0,b=1,col="red",lwd=2)
abline(lm(out1$mean$N~N),lwd=2,lty=2)

plot(N1,out1$mean$N1,pch=19,
     xlab="True N (Verified area)", ylab="Estimated N (Verified area)")
abline(a=0,b=1,col="red",lwd=2)

plot(N2,out1$mean$N2,pch=19,
     xlab="True N (Unverified area)", ylab="Estimated N (Unverified area)")
abline(a=0,b=1,col="red",lwd=2)



