
library(ggplot2)
library(gridExtra)


resu <- readRDS(file.path("Simulations", "results_scores-prop.rds"))

do.call(rbind, lapply(resu, "[[", 1))



# Mean error in abundance estimates -------------------------------------------------------------------------------

N.bias <- list()

for(i in 1:length(resu)){
  error <- (resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"mean"] - 
               resu[[i]]$data$true$N[,"N"])  / 
                                          resu[[i]]$data$true$N[,"N"]
  
  cv <- resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"sd"] / 
              resu[[i]]$result[grep("N[", rownames(resu[[i]]$result), fixed=T),"mean"]
  
  N.bias[[i]] <- cbind(resu[[i]]$scen, error, cv)
}

N.bias <- do.call(rbind, N.bias)

N.bias$over <- factor(N.bias$mean.FP, labels=c("low","medium","high"))

f1 <- ggplot(data=N.bias, aes(x=over, y=error, fill=as.factor(prop2obs))) +
  geom_boxplot() +
  labs(x="Overlap in scores distributions", y="Relative error in abundance estimates",
       fill="Verified proportion") +
  theme_bw(base_size=14) +
  theme(legend.position="top")


#* Coefficient of variation in N 

f2 <- ggplot(data=N.bias, aes(x=as.factor(prop2obs), y=cv, fill=over)) +
  geom_boxplot() +
  labs(x="Verified proportion", y="Coefficient of variation in abundance estimates",
       fill="Overlap in scores distributions") +
  theme_bw(base_size=14) +
  theme(legend.position="top")


grid.arrange(f1, f2, ncol=2)



# Root mean squared error -----------------------------------------------------------------------------------------
results <- unique(do.call(rbind, lapply(resu, "[[", 1)))

results$over <- factor(results$mean.FP, labels=c("low","medium","high"))

rmse <- as.data.frame.table(tapply(N.bias$error, list(over=N.bias$over, prop2obs=N.bias$prop2obs),
                                   function(x) sqrt(mean(x^2)) ))

results <- merge(results, rmse, by=c("over","prop2obs"))

write.table(results, "clipboard", sep = "\t", row.names = FALSE)



# Scores distributions --------------------------------------------------------------------------------------------

hist(resu[[1]]$data$scr.data$score)

scr.seq <- seq(-5, 5, 0.05)

FP.seq <- unique(do.call(rbind, lapply(resu, "[[", 1))$mean.FP)

TP.val <- unique(do.call(rbind, lapply(resu, "[[", 1))$mean.TP)

sd.val <- unique(do.call(rbind, lapply(resu, "[[", 1))$sd.scr)


tab <- data.frame(
  dens = c(  dnorm(scr.seq, TP.val, sd.val),
      dnorm(scr.seq, FP.seq[1], sd.val),
      dnorm(scr.seq, FP.seq[2], sd.val),
      dnorm(scr.seq, FP.seq[3], sd.val)),
  score = rep(scr.seq, 4),
  type = factor(rep(c("TP","FP.low","FP.mid","FP.hig"),
             each=length(scr.seq)), levels=c("TP","FP.low","FP.mid","FP.hig"))
)


ggplot(tab, aes(x=plogis(score), y=dens, col=type)) +
  geom_line(linewidth=1.4, alpha=0.4) +
  scale_color_manual(values=c("#00CD00","#FF7256", "#FF0000", "#8B2500")) +
  labs(x="Score values", y="Density",
       color="Type") +
  theme_classic(base_size=14) +
  theme(legend.position="right")

ggplot(tab, aes(x=plogis(score),y=dens,ymin=0, ymax=dens, group=type, col=type, fill=type)) +
  geom_line(linewidth=1, alpha=0) +
  geom_ribbon(alpha=0.2) +
  scale_color_manual(values=c("#00CD00","#FF7256", "#FF0000", "#8B2500")) +
  scale_fill_manual(values=c("#00CD00","#FF7256", "#FF0000", "#8B2500")) +
  labs(x="Score values", y="Density") +
  theme_classic(base_size=14) +
  theme(legend.position="right")

