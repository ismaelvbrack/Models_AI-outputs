
library(ggplot2)
library(gridExtra)


resu <- readRDS(file.path("Simulations", "results_occu-scenarios.rds"))

do.call(rbind, lapply(resu[which(sapply(resu, length)==3)], "[[", 1))

unique(do.call(rbind, lapply(resu, "[[", 1)))

resu <- resu[which(sapply(resu, length)==3)]

# Mean error in abundance estimates -------------------------------------------------------------------------------

occu.bias <- list()

for(i in 1:length(resu)){
  mean.est <- resu[[i]]$result[grep("psi", rownames(resu[[i]]$result), fixed=T),"mean"]
  error <- (resu[[i]]$result[grep("psi", rownames(resu[[i]]$result), fixed=T),"mean"] - 
              resu[[i]]$data$true$params["psi"])
  
  sd <- resu[[i]]$result[grep("psi", rownames(resu[[i]]$result), fixed=T),"sd"] 
  
  occu.bias[[i]] <- cbind(resu[[i]]$scen, mean.est, error, sd)
}

occu.bias <- do.call(rbind, occu.bias)


# Figures ---------------------------------------------------------------------------------------------------------

tab2 <- occu.bias

names(tab2)[names(tab2) %in% c("psi","phi")] <- c("Occupancy","Availability")

# f1 <- ggplot(data=tab2, aes(x=as.factor(nsites), y=error, fill=as.factor(prop.verif))) +
#   geom_boxplot(outlier.shape=1) +
#   facet_grid(Availability~Occupancy,labeller="label_both") +
#   geom_hline(yintercept=0,color="red2", linetype="dashed",linewidth=1) +
#   scale_fill_manual(values=c("#20B2AA", "#FFB90F")) +
#   labs(x="Number of sampled sites", y="Absolute error in occupancy estimate",
#        fill="Verified proportion",psi="Occupancy",phi="Availability") +
#   theme_bw(base_size=14) +
#   theme(legend.position="top")
f1 <- ggplot(data=tab2, aes(x=as.factor(nsites), y=mean.est, fill=as.factor(prop.verif))) +
  geom_boxplot(outlier.shape=1) +
  geom_hline(aes(yintercept=Occupancy),color="red2", linetype="dashed",linewidth=1) +
  facet_grid(Availability~Occupancy,labeller="label_both") +
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
  scale_fill_manual(values=c("#20B2AA", "#FFB90F")) +
  labs(x="Number of sampled sites", y="Occupancy estimates",
       fill="Verified proportion",psi="Occupancy",phi="Availability") +
  theme_bw(base_size=14) +
  theme(legend.position="top")

ggsave(file.path("Fig simul-results_occu simuls.png"),f1,
       width=20,height=16,unit="cm",dpi=400)

#* Coefficient of variation in N 

ggplot(data=occu.bias, aes(x=as.factor(nsites), y=sd, fill=as.factor(prop.verif))) +
  geom_boxplot(outlier.shape=1) +
  facet_grid(phi~psi,labeller="label_both") +
  scale_fill_manual(values=c("#20B2AA", "#FFB90F")) +
  labs(x="Number of sampled sites", y="Absolute error in occupancy estimate",
       fill="Verified proportion",psi="Occupancy",phi="Availability") +
  theme_bw(base_size=14) +
  theme(legend.position="top")

# Root mean squared error -----------------------------------------------------------------------------------------
results <- unique(do.call(rbind, lapply(resu, "[[", 1)))

rmse <- as.data.frame.table(tapply(occu.bias$error,
                                   list(psi=occu.bias$psi, phi=occu.bias$phi,
                                        nsites=occu.bias$nsites, prop.verif=occu.bias$prop.verif),
                                   function(x) sqrt(mean(x^2)) ))

results <- merge(results, rmse, by=c("psi","phi","nsites","prop.verif"))

write.table(results, "clipboard", sep = "\t", row.names = FALSE)

