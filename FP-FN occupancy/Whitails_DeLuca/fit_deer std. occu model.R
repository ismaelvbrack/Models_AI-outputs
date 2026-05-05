
library(jagsUI)

root <- file.path("Whitails_DeLuca")
# Import data -----------------------------------------------------------------------------------------------------
site.imgs <- read.csv(file.path(root,"data","sites_images_annotations.csv"))

siteCovs <- read.table(file.path(root,"data","siteCovars.txt"),h=T)

temps <- read.table(file.path(root,"data","Temperatures60cm.txt"),h=T)
rain <- read.table(file.path(root,"data","Rainfall.txt"),h=T)

rain <- ifelse(rain>0,1,0)

nsites = nrow(siteCovs)
ndays = ncol(temps)

# Organize data ---------------------------------------------------------------------------------------------------

Y <- table(ifelse(site.imgs$common_name=="White-tailed Deer",1,0),
           factor(site.imgs$site, levels=siteCovs$placename),
           factor(site.imgs$day.id, levels=1:ndays)
)[2,,]

Y <- ifelse(Y>0,1,0)

dat <- list(
  Y=Y,
  forest=as.numeric(scale(siteCovs$forest)),
  twi=as.numeric(scale(siteCovs$twi)),
  rain=rain,
  
  ndays=ndays,
  nsites=nsites
)

inits <- function() list(
  z=rep(1,nsites)
)

# Parameters monitored
params <- c("alpha0","alpha1",
              "beta0","beta1","beta2",
              "psi","p")

# MCMC settings
ni <- 50000; nt <- 1; nb <- 20000; nc <- 3; na <- 1000

out1 <- jags(dat, inits, params,
             #model.file=file.path(root,"std.occu_psi(.)p(.).R"),
             model.file=file.path(root,"std.occu_psi(for)p(twi+rain).R"),
             n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
             parallel=T)

saveRDS(out1,file.path(root,"resu_std.occu5.rds"))

out1 <- readRDS(file.path(root,"resu_std.occu5.rds"))


  