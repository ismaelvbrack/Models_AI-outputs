
library(ggplot2)
library(gridExtra)


resu <- c(
  readRDS(file.path("Simulations", "results1_very-prec-prop.rds")),
  readRDS(file.path("Simulations", "results2_very-prec-prop.rds")),
  readRDS(file.path("Simulations", "results3_very-prec-prop.rds")),
  readRDS(file.path("Simulations", "results4_very-prec-prop.rds"))
)

unique(do.call(rbind, lapply(resu, "[[", 1)))

# Mean error and CVin abundance estimates -------------------------------------------------------------------------------

N.bias <- list()

for(i in 1:length(resu)){
  error <- (resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"mean"] - 
              as.vector(resu[[i]]$data$true$abund$N)
            )  / 
    as.vector(resu[[i]]$data$true$abund$N)
  
  cv <- resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"sd"] / 
    resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"mean"]
  
  N.bias[[i]] <- cbind(resu[[i]]$scen, error, cv)
}

N.bias <- do.call(rbind, N.bias)

N.bias$prec.var <- factor(N.bias$sd.prec, labels=c("low","high"))


# Figures ---------------------------------------------------------------------------------------------------------

#* Relative mean error
f1 <- ggplot(data=N.bias, aes(x=prec.var, y=error, fill=as.factor(nsites))) +
  geom_violin(position = position_dodge(0.8)) +
  geom_boxplot(position = position_dodge(0.8),width=0.2,outlier.size=0.7) +
  labs(x="Variation in precision between flights", y="Relative error in abundance estimates",
       fill="Number of sites") +
  geom_hline(yintercept=0,color="red2", linetype="dashed",linewidth=1) +
  theme_bw(base_size=14) +
  theme(legend.position="top")

#* Coefficient of variation in N 

f2 <- ggplot(data=N.bias, aes(x=as.factor(nsites), y=cv,fill=prec.var)) +
  geom_boxplot() +
  labs(fill="Variation in precision between flights", y="Coefficient of variation in abundance estimates",
      x="Number of sites") +
  theme_bw(base_size=14) +
  theme(legend.position="top")

# Root mean squared error -----------------------------------------------------------------------------------------
results <- unique(do.call(rbind, lapply(resu, "[[", 1)))

rmse <- as.data.frame.table(tapply(N.bias$error,
                                   list(nsites=N.bias$nsites, sd.prec=N.bias$sd.prec),
                                   function(x) sqrt(mean(x^2)) ))

results <- merge(results, rmse, by=c("nsites","sd.prec"))

write.table(results, "clipboard", sep = "\t", row.names = FALSE)
