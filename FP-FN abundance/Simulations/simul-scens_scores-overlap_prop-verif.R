

##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###*** Simulate scenarios with different overlaps between the scores values
###***  and different proportions of verified data
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(jagsUI)
library(foreach)
library(doParallel)

source(file.path("Simulations","func_simul_colony-AIcounts.R"))

# Define scenarios for simulation ---------------------------------------------------------------------------------

abund.rng = c(10,100) # abundance range between colonies

n.sites = 20 # number of colonies
prop2obs = c(0.1, 0.3, 0.5) # proportion of the site with verified samples
ph = 0.8 # human detection probability
pa = 0.7 # algorithm detection probability (Recall)
precision = 0.5 # algorithm proportion of true positives (Precision)
mean.FPs = c(qlogis(0.4), qlogis(0.6), qlogis(0.8)) # mean score values (normal dist.) for FALSE positives 
mean.TPs = qlogis(0.85) # mean score values (normal dist.) for TRUE positives 
sd.scores = 0.5 # std. deviation for the scores normal distribution

scenarios <- expand.grid(
  nsites = n.sites,
  prop2obs = prop2obs,
  ph=ph, pa=pa, precision=precision,
  mean.FP=mean.FPs, mean.TP=mean.TPs,
  sd.scr = sd.scores
)


# Simulation settings ---------------------------------------------------------------------------------------------

n.simul <- 200

# Path to model code
model.pth <- file.path("model1_colony-count.R")

# Parameters monitored
params <- c("pa","ph","mu.scr","precision","N","N1","N2")

# MCMC settings
ni <- 10000; nt <- 1; nb <- 4000; nc <- 3; na <- 100

# Object to contains final results
resu <- list()

##** Running in parallel
registerDoParallel(makeCluster(15))

# RUUUUUUN! -------------------------------------------------------------------------------------------------------

for(scen in 1:nrow(scenarios)){
  
  cat("\nScenario:",scen,"/",nrow(scenarios))
  print(Sys.time())
  
  resu.inn <- foreach(w=1:n.simul, .errorhandling = 'remove') %dopar% {
    
    require(jagsUI)
    source(file.path("Simulations","func_simul_colony-AIcounts.R"))
    
    ns <- scenarios[scen,"nsites"]
    
    # Simulate data ---------------------------------------------------------------------------------------------------
    
    abund <- round(runif(ns, abund.rng[1],abund.rng[2]))
    
    tab <- simul_counts(
      nsites = ns, 
      N = abund, 
      pa = rep(scenarios[scen,"pa"], ns), 
      ph = 0.7, 
      precision = rep(scenarios[scen,"precision"], ns),
      prop2obs = scenarios[scen,"prop2obs"],
      mu.scores = c(scenarios[scen,"mean.FP"], scenarios[scen,"mean.TP"]), 
      sd.scr = rep(scenarios[scen,"sd.scr"], 2) 
    )
    
    # Bundle data for JAGS
    dat <- list(
      nsites= ns,
      
      class= tab$scr.data$obs.class,
      score= tab$scr.data$score,
      site= tab$scr.data$site,
      nobj= nrow(tab$scr.data),
      
      Y= tab$counts[,1:3],
      ndet1= rowSums(tab$counts[,1:3]),
      
      ndet2=tab$counts[,"unverified"]
    )
    
    # Initial values
    inits <- function() list(
      class=ifelse(is.na(tab$scr.data$obs.class),1,NA)
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

saveRDS(resu, file.path("Simulations", "results_scores-prop.rds") )


