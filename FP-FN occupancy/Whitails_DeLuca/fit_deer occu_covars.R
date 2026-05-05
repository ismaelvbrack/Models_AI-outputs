

library(jagsUI)
library(ggplot2)

root <- file.path("Whitails_DeLuca")
# Import data -----------------------------------------------------------------------------------------------------

deer <- read.table(file.path(root,"data","deer_detections_scores.txt"),h=T)
nimgs <- read.table(file.path(root,"data","nimgs_site-visit.txt"),h=T)
ndets <- readRDS(file.path(root,"data","3Darray_ndets_deer.rds"))
siteCovs <- read.table(file.path(root,"data","siteCovars.txt"),h=T)

temps <- read.table(file.path(root,"data","Temperatures60cm.txt"),h=T)
rain <- read.table(file.path(root,"data","Rainfall.txt"),h=T)

rain <- ifelse(rain>0,1,0)

# Oraganize data for JAGS -----------------------------------------------------------------------------------------

# Bundle data
dat <- list(
  ndets=ndets,
  
  class=deer$obs.class,
  score=deer$confidence,
  
  site=deer$site.id,
  visit=deer$day.id,
  image=deer$image.id,
  
  forest=as.numeric(scale(siteCovs$forest)),
  twi=as.numeric(scale(siteCovs$twi)),
  open=as.numeric(scale(siteCovs$open)),
  temp=temps,
  rain=rain,
  
  nimgs=nimgs,
  
  nsites=nrow(siteCovs),
  nvisits=ncol(nimgs),
  nobj=nrow(deer)
)

# Initial values

inits <- function() list(
  z=rep(1,nrow(siteCovs)),
  u=matrix(1,nrow(siteCovs),ncol(nimgs)),
  class=ifelse(is.na(deer$obs.class),1,NA)
)

# Parameters monitored
params <- c("alpha0","alpha1",
            "beta0","beta1",
            "mu.ome","sd.ome",
            'delta0','delta1',
            "psi","phi","lam","ome",
            "mu","sd","class")

# MCMC settings
ni <- 100000; nt <- 1; nb <- 40000; nc <- 3; na <- 1000


# Run model! ------------------------------------------------------------------------------------------------------

out7 <- jags(dat, inits, params,
             #model.file=file.path(root,"mod1_psi(.)phi(.)lam(.)ome(.).R"),
             model.file=file.path(root,"mod7_psi(for)phi(.)lam(open)ome(sites).R"),
             n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
             parallel=T)

plot(out7)

saveRDS(out7,file.path(root,"resu_mod7.rds"))


# -----------------------------------------------------------------------------------------------------------------
out2 <- readRDS(file.path(root,"resu_mod2.rds"))
out4 <- readRDS(file.path(root,"resu_mod4.rds"))
out5 <- readRDS(file.path(root,"resu_mod5.rds"))

out4$summary[1:8,]

source("C:/Users/i.verrastrobrack/OneDrive/CODING/predictJAGS_func/func_predictJAGS.R")

newfor <- data.frame(
  forest=seq(min(scale(siteCovs$forest)), max(scale(siteCovs$forest)), ,100)
)

pred.psifor <- predictJAGS(out4, c("alpha0","alpha1"),newdata=newfor,link="logit")

pred.psifor$forest <- pred.psifor$forest * attr(scale(siteCovs$forest),"scaled:scale") +
  attr(scale(siteCovs$forest),"scaled:center")

ggplot(data=pred.psifor,
       aes(x=forest,y=mean,ymin=`2.5%`,ymax=`97.5%`)) +
  geom_line(col="darkgreen", linewidth=1.2) +
  geom_ribbon(fill="darkgreen", alpha=0.3) + 
  labs(x="% of Upland Forest (500m buffer)",
       y="Occupancy probability") +
  theme_classic(base_size = 14)

