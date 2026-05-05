
library(ggplot2)
library(gridExtra)


resu <- readRDS(file.path("Simulations", "results_sites-abund-prop.rds"))

do.call(rbind, lapply(resu[which(sapply(resu, length)==3)], "[[", 1))

resu <- resu[which(sapply(resu, length)==3)]

# Mean error in abundance estimates -------------------------------------------------------------------------------

N.bias <- list()

for(i in 1:length(resu)){
  mean.est <- 
  error <- (resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"mean"] - 
              resu[[i]]$data$true$N[,"N"])
  rel.error <- error / resu[[i]]$data$true$N[,"N"]
  
  cv <- resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"sd"] / 
    resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"mean"]
  
  N.bias[[i]] <- cbind(resu[[i]]$scen, error, rel.error, cv)
}

N.bias <- do.call(rbind, N.bias)


f1 <- 
ggplot(data=N.bias, aes(x=abund, y=rel.error, fill=as.factor(prop2obs))) +
  geom_boxplot(outlier.shape=1) +
  geom_hline(yintercept=0,color="red3", linetype="dashed",linewidth=1) +
  scale_fill_manual(values=c("#20B2AA", "#FFB90F", "gray80")) +
  coord_cartesian(ylim=c(-1,1.5)) +
  facet_wrap(.~nsites,labeller="label_both") +
  labs(x="Abundance", y="Relative error in abundance estimates",
       fill="Verified proportion") +
  theme_bw(base_size=14) +
  theme(legend.position="top")

#* Coefficient of variation in N 

f2 <- ggplot(data=N.bias, aes(x=abund, y=cv, fill=as.factor(prop2obs))) +
  geom_boxplot(outlier.shape=1) +
  labs(x="Abundance", y="Coefficient of variation in abundance estimates",
       fill="Verified proportion") +
  coord_cartesian(ylim=c(0,2)) +
  facet_wrap(.~nsites,labeller="label_both") +
  theme_bw(base_size=14) +
  #geom_hline(yintercept=0.2, col=2, linetype="dashed") +
  theme(legend.position="none")

grid.arrange(f1,f2,ncol=1)


ggsave(file.path("Fig simul-results_abund prop-ver.png"),f1,
       width=24,height=14,unit="cm",dpi=400)

# Root mean squared error -----------------------------------------------------------------------------------------
results <- unique(do.call(rbind, lapply(resu, "[[", 1)))

rmse <- as.data.frame.table(tapply(N.bias$error,
                                   list(abund=N.bias$abund, nsites=N.bias$nsites, prop2obs=N.bias$prop2obs),
                     function(x) sqrt(mean(x^2)) ))

results <- merge(results, rmse, by=c("abund","nsites","prop2obs"))

write.table(results, "clipboard", sep = "\t", row.names = FALSE)
