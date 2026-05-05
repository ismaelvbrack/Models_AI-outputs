
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###*** Simulate scenarios with different number of sites, different abundance, 
###***    and varying precision between sites
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(jagsUI)
library(foreach)
library(doParallel)

source(file.path("Simulations","func_simul_occuFPFN.R"))

# Define scenarios for simulation ---------------------------------------------------------------------------------

n.sites = c(30,60) # number of colonies
n.visits = 10 # number of occasions
prop.verif = c(0.1,0.3) # proportion of the site with verified samples
mean.nimgs = 8 # mean number of photos per occasion per day

psis = c(0.2, 0.6)
phis = c(0.2, 0.6)
lam = 1
ome = 0.4

mean.FPs = c(qlogis(0.5)) # mean score values (normal dist.) for FALSE positives 
mean.TPs = qlogis(0.8) # mean score values (normal dist.) for TRUE positives 
sd.scores = 0.5 # std. deviation for the scores normal distribution

scenarios <- expand.grid(
  nsites = n.sites,
  nvisits = n.visits,
  mean.nimgs = mean.nimgs,
  prop.verif = prop.verif,
  psi = psis, phi = phis,
  lam = lam, ome = ome,
  mean.FP=mean.FPs, mean.TP=mean.TPs,
  sd.scr = sd.scores
)

# Simulation settings ---------------------------------------------------------------------------------------------

n.simul <- 200

# Path to model code
model.pth <- file.path("Simulations","occu_model_3level_normal scores.R")

# Parameters monitored
params <- c("psi","phi","lam","ome","mu")

# MCMC settings
ni <- 18000; nt <- 1; nb <- 10000; nc <- 3; na <- 100

# Object to contains final results
resu <- list()

##** Running in parallel
registerDoParallel(makeCluster(15))

# RUUUUUUN! -------------------------------------------------------------------------------------------------------

for(scen in 1:nrow(scenarios)){
  
  cat("\nScenario:",scen,"/",nrow(scenarios))
  print(Sys.time())
  
  resu.inn <- foreach(w=1:n.simul, .errorhandling = 'pass') %dopar% {
    
    require(jagsUI)
    source(file.path("Simulations","func_simul_occuFPFN.R"))
    
    # Simulate data ---------------------------------------------------------------------------------------------------

    tab <- simul_occuFPFN(
      nsites = scenarios[scen,"nsites"], 
      nvisits = scenarios[scen,"nvisits"],
      mean.nimgs = scenarios[scen,"mean.nimgs"],
      prop.verif = scenarios[scen,"prop.verif"],
      
      psi = scenarios[scen,"psi"],
      phi = scenarios[scen,"phi"],
      lam = scenarios[scen,"lam"],
      ome = scenarios[scen,"ome"],
      
      mu.scores = c(scenarios[scen,"mean.FP"], scenarios[scen,"mean.TP"]), 
      sd.scr = rep(scenarios[scen,"sd.scr"],2)
    )
    
    # Bundle data for JAGS
    dat <- list(
      # Sampling
      nsites = scenarios[scen,"nsites"], 
      nvisits = scenarios[scen,"nvisits"],
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
    
    # Run model -------------------------------------------------------------------------------------------------------
    
    out1 <- jags(dat, inits, params,
                     model.file=model.pth,
                     n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
                     parallel=F,verbose=F)

    
    outs <- list(scen = scenarios[scen,], data = tab, result = out1$summary)
    
    return(outs)
    
  } # iters
  
  resu <- c(resu, resu.inn)
  
} # scenario

saveRDS(resu, file.path("Simulations", "results_occu-scenarios.rds") )
