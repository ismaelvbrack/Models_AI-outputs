

##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###*** Simulate scenarios varying precision between sites and flights
###***    and compare estimates of models that account or not for flight variations 
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(jagsUI)
library(foreach)
library(doParallel)

source(file.path("Simulations","func_simul_colony-AIcounts_flights.R"))

# Define scenarios for simulation ---------------------------------------------------------------------------------
abund.rng = c(20,100) # abundance range between colonies

n.sites = c(10,40) # number of colonies
n.flights = 8
prop2obs = c(0.1) # proportion of the site with verified samples
ph = 0.8 # human detection probability
pa = 0.7 # algorithm detection probability (Recall)

prec.rng = c(0.2, 0.8) # algorithm proportion of true positives (Precision)
prec.sd = c(0.2,0.5) # std. dev. of precision within sites

mean.FPs = c(qlogis(0.5)) # mean score values (normal dist.) for FALSE positives 
mean.TPs = qlogis(0.8) # mean score values (normal dist.) for TRUE positives 
sd.scores = 0.5 # std. deviation for the scores normal distribution

scenarios <- expand.grid(
  nsites = n.sites,
  nflights = n.flights,
  prop2obs = prop2obs,
  ph=ph, pa=pa, sd.prec=prec.sd,
  mean.FP=mean.FPs, mean.TP=mean.TPs,
  sd.scr = sd.scores
)

# Simulation settings ---------------------------------------------------------------------------------------------

n.simul <- 200

# Path to model code
model.pth <- file.path("model2_combined flights.R")

# Parameters monitored
params <- c("ph","pa","mu.prec","sd.prec",
            "mu.scr","sd.scr", 
            "precision",
            "N","N1","N2")

# MCMC settings
ni <- 10000; nt <- 1; nb <- 4000; nc <- 3; na <- 100

# Object to contains final results
resu <- list()

##** Running in parallel
registerDoParallel(makeCluster(8))


# RUUUUUUN! -------------------------------------------------------------------------------------------------------

#for(scen in 1:nrow(scenarios)){
  
  scen = 4
  cat("\nScenario:",scen,"/",nrow(scenarios))
  print(Sys.time())
  
  resu.inn <- foreach(w=1:n.simul, .errorhandling = 'remove') %dopar% {
    
    require(jagsUI)
    source(file.path("Simulations","func_simul_colony-AIcounts_flights.R"))
    
    ns <- scenarios[scen,"nsites"]
    
    # Simulate data ---------------------------------------------------------------------------------------------------
    
    abund <- t(sapply(round(runif(ns, abund.rng[1],abund.rng[2])), rep, each=8))
    
    tab <- simul_counts(
      nsites = ns, 
      N = abund, 
      pa = rep(scenarios[scen,"pa"], ns), 
      ph = 0.7, 
      mean.prec = runif(ns, prec.rng[1],prec.rng[2]),
      sd.prec = scenarios[scen,"sd.prec"], 
      prop2obs = scenarios[scen,"prop2obs"],
      mu.scores = c(scenarios[scen,"mean.FP"], scenarios[scen,"mean.TP"]), 
      sd.scr = rep(scenarios[scen,"sd.scr"], 2) 
    )
    
    # IMPORTANT: only the first flight is verified!
    # observed class for flights>1 must be NA
    tab$scr.data[which(tab$scr.data$flight > 1),"obs.class"] <- NA
    
    # combine all detections for flights>1
    ndets2 <- cbind(tab$counts[,"unverified",1],
                    apply(tab$counts[,,2:n.flights], c(1,3), sum))
    
    # Bundle data for JAGS
    dat <- list(
      nsites= ns,
      nflight=rep(n.flights,ns),
      
      class= tab$scr.data$obs.class,
      score= tab$scr.data$score,
      site= tab$scr.data$site,
      nobj= nrow(tab$scr.data),
      
      Y= tab$counts[,1:3,1],
      ndet1= rowSums(tab$counts[,1:3,1]),
      
      # combine counts for flights >1
      ndet2=ndets2
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
  
#} # scenario

saveRDS(resu, file.path("Simulations", "results4_very-prec-prop.rds") )







