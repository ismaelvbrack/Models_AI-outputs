
library(jagsUI)
library(ggplot2)

# Import data -------------------------------------------------------------------------------------------

counts <- read.table(file.path("EvergladesBirds","counts_12colonies_all-dates.txt"),h=T,check.names=F)

counts$date <- as.Date(counts$date)

scr.data <- read.csv(file.path("EvergladesBirds","AI-detections_12colonies_all-dates.csv"))

# Exclude colony without verification data
counts <- counts[-which(counts$site == "Hidden"),]
scr.data <- scr.data[-which(scr.data$site == "Hidden"),]


ggplot(data=counts, aes(x=log(unverified), y=TP / (TP+FP)) ) +
         geom_point(size=2) +
  #geom_smooth(method="gam") +
  theme_bw(base_size=14) +
  labs(x="Unverified counts", y="Poportion of TPs in verified counts") +
  scale_y_continuous(breaks=seq(0,1,.2))

ggplot(scr.data, aes(x=site, y=score)) +
  geom_boxplot(fill="gray80") +
  theme_bw(base_size=14) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  theme(axis.text.x=element_text(angle=90))

# Organize data to JAGS ----------------------------------------------------------------------------------

##* Organize objects
sites <- unique(counts$site)

ndets2 <- matrix(NA, nrow=length(sites), ncol=max(table(counts$site)))
ndets1 <- as.numeric()
Y <- matrix(NA, nrow=length(sites), ncol=3)
scr.data$site.id <- NA

for(s in 1:length(sites)){
  aqui <- which(counts$site == sites[s])
  
  Y[s,] <- as.numeric(counts[aqui[1],c("11","10","01")])
  ndets1[s] <- counts[aqui[1], "TP"]
  ndets2[s,1:length(aqui)] <- counts[aqui, "unverified"]
  
  scr.data$site.id[which(scr.data$site == sites[s])] <- s
}

any(ndets2[,1] != counts[1:11,"unverified"]) #FALSE!

# Bundle data
dat <- list(
  nsites=length(sites),
  
  class=ifelse(scr.data$det.type=="TP",1,0),
  score=qlogis(scr.data$score),
  site=scr.data$site.id,
  nobj=nrow(scr.data),
  
  Y=Y,
  ndet1=ndets1,
  ndet2=ndets2,
  
  nflight=as.numeric(table(counts$site))
)

# Initial values
inits <- function() list(
  class=ifelse(is.na(dat$class),1,NA)
)

# Parameters monitored
params <- c("ph","mu.prec","sd.prec",
            "mu.pa","sd.pa",
            "mu.scr","sd.scr", # normal dist
            #"mu.scr","phi", # beta dist
            "pa","precision",
            "TP2",
            "N","N1","N2")

# MCMC settings
ni <- 10000; nt <- 1; nb <- 4000; nc <- 3; na <- 500

# RUN!
out2 <- jags(dat, inits, params,
             model.file=file.path("model2_combined flights.R"),
             n.chains=nc,n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na,
             parallel=T)

beepr::beep(2)

saveRDS(out2,file.path("EvergladesBirds","resu_model2_combined flights.rds"))








