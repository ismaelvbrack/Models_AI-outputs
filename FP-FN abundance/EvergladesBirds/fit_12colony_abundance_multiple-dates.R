
library(jagsUI)
library(ggplot2)

# Import data -------------------------------------------------------------------------------------------

counts <- read.table(file.path("EvergladesBirds","counts_12colonies_all-dates.txt"),h=T,check.names=F)

counts$date <- as.Date(counts$date)

scr.data <- read.csv(file.path("EvergladesBirds","AI-detections_12colonies_all-dates.csv"))

table(scr.data$det.type, scr.data$site)

ggplot(data=counts,
       aes(x=date, y=site, col=ifelse(TP>0 | FP>0,"HQ",NA))) +
  geom_point(size=3) +
  theme_bw(base_size=14) +
  theme(legend.position="none")

#* Separate count data
ver.sites <- counts[which(!is.na(counts$TP)),] # sites with verification data

unver.sites <- counts[which(is.na(counts$TP)),] # sites wiithout verification


#* Separate scores data
aqui <- which(paste0(scr.data$site,scr.data$date) %in%
        paste0(ver.sites$site,ver.sites$date))

scr.versites <- scr.data[aqui,] # sites with verification data
 
scr.unvers <- scr.data[-aqui,] # sites wiithout verification 

unver.sites$flight.id <- factor(paste0(unver.sites$site,"_", format(unver.sites$date,"%b%d")),
                                levels=paste0(unver.sites$site,"_", format(unver.sites$date,"%b%d")))

scr.unvers$flight.id <- unver.sites[match(paste0(scr.unvers$site, scr.unvers$date),
                                          paste0(unver.sites$site, unver.sites$date)), "flight.id"]

any(unver.sites$unverified != table(scr.unvers$flight.id)) # FALSE!

# Only flight to predict ------------------------------------------------------------------------------------------

resu1 <- readRDS(file.path("EvergladesBirds","resu_model2.rds"))

#plot(resu1, par="precision")

#* Transform precision probability estimates into beta distribution parameters
p.mu = resu1$mean$precision
p.var = (resu1$sd$precision)^2

beta.pars <- cbind(
  p.mu * ((p.mu*(1-p.mu)/p.var) - 1),
  (1-p.mu)*((p.mu*(1-p.mu)/p.var) - 1)
)
# See if they match
par(mfrow=c(3,4))
for(i in 1:nrow(beta.pars)){
  plot(density(resu1$sims.list$precision[,i]))
  lines(density(rbeta(length(resu1$sims.list$precision[,i]), beta.pars[i,1], beta.pars[i,2])),col=2)
}


# Exclude Hidden colony
unver.sites <- unver.sites[-which(unver.sites$site=="Hidden"),]

scr.unvers <- scr.unvers[-which(scr.unvers$site=="Hidden"),]

# Bundle data
dat <- list(
  # estimated pars
  beta.pars=beta.pars,
  mu.scr=resu1$mean$mu.scr,
  sd.scr=resu1$mean$sd.scr,
  pa=resu1$mean$pa,
  
  # unverified flights
  ndets=unver.sites$unverified,
  score=scr.unvers$score,
  nobj=nrow(scr.unvers),
  nflight=nrow(unver.sites),
  
  nsites=length(ver.sites$site),
  site=as.numeric(factor(unver.sites$site, levels=ver.sites$site)),
  scr.site=as.numeric(factor(scr.unvers$site, levels=ver.sites$site))
  
)

# Initial values
inits <- function() list(
  class=rep(1, dat$nobj)
)

# Parameters monitored
params <- c("N","TP","precision")

# MCMC settings
ni <- 8000; nt <- 1; nb <- 4000; nc <- 3; na <- 500

# RUN!
out2 <- jags(dat, inits, params,
             model.file=file.path("model2_predict unver flights.R"),
             n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
             parallel=T)

beepr::beep(2)

saveRDS(out2,file.path("EvergladesBirds","resu3_model2_predictions.rds"))

