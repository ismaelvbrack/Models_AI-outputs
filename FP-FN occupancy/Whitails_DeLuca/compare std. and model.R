

library(ggplot2)
library(jagsUI)

root <- file.path("Whitails_DeLuca")

# Comparison with covariates --------------------------------------------------------------------------------------
out1 <- readRDS(file.path(root,"resu_mod6.rds"))
out2 <- readRDS(file.path(root,"resu_std.occu2.rds"))

##** All parameters
pars <- as.data.frame(rbind(
  out1$summary[c(1:7,494:496),],
  out2$summary[1:4,]
))

pars$par <- factor(c(
  # FP-FN model
  "Occupancy Int.",
  "Forest Slope",
  "Availability Int.",
  "FP RanEff Mean",
  "FP RanEff SD",
  "TP rate Int.",
  "TWI Slope",
  "FP Scores Mean",
  "TP Scores Mean",
  "Scores SD",
  # Standard model
  "Occupancy Int.",
  "Forest Slope",
  "Detection Int.",
  "TWI Slope"
),
levels=rev(c(
  "Occupancy Int.",
  "Forest Slope",
  "Detection Int.",
  "TWI Slope",
  "Availability Int.",
  "TP rate Int.",
  "FP Normal Mean",
  "FP RanEff Mean",
  "FP RanEff SD",
  "FP Scores Mean",
  "TP Scores Mean",
  "Scores SD"
)))

pars$model <- rep(c("FP-FN Occupancy model","Standard Occupancy model"),
                  c(10,4))

fig.pars <- ggplot(data=pars, aes(x=mean,xmin=`2.5%`,xmax=`97.5%`,y=par,color=model)) +
  geom_point(size=3,position=position_dodge(0.7)) +
  geom_errorbarh(linewidth=.8,height=0.3,position=position_dodge(0.7)) +
  scale_color_manual(values=c("gold3","darkblue")) +
  #scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  labs(x="Estimate", y="Parameter",color="") + #title="Whitetails occurrence model"
  theme_bw(base_size=14) +
  theme(legend.position="top",
        legend.text=element_text(size=14)) 

ggsave(file.path("Fig Param ests Whitails DeLuca.png"),fig.pars,
       width=22,height=12,unit="cm",dpi=400)

###*** Parameter values
plogis(as.matrix(pars[c("alpha0","beta0"), c("mean","2.5%","97.5%")]))

exp(as.matrix(pars[c("delta0","mu.ome"), c("mean","2.5%","97.5%")]))

fps <- out1$summary[grep("ome[",rownames(out1$summary),fixed=T),
                              c("mean","2.5%","97.5%")]
range(fps[,"mean"])

tps <- out1$summary[grep("lam[",rownames(out1$summary),fixed=T),
                    c("mean","2.5%","97.5%")]
range(tps[,"mean"])

# Comparison with the constant occupancy model -----------------------------------
out1 <- readRDS(file.path(root,"resu_mod1.rds"))
out2 <- readRDS(file.path(root,"resu_std.occu.rds"))

psi.samps <- data.frame(
  samps=c(out1$sims.list$psi,out2$sims.list$psi),
  model=rep(c("AI-scores (FP-FN)","Manual Std. (FN)"),
            c(out1$mcmc.info$n.samples,out2$mcmc.info$n.samples))
)

psi.resu <- as.data.frame(rbind(
  out1$summary["psi",],
  out2$summary["psi",]
))
psi.resu$model <- c("AI-scores (FP-FN)","Manual Std. (FN)")


ggplot(data=psi.samps, aes(x=samps,fill=model,color=model)) +
  geom_density(alpha=0.3) +
  scale_color_manual(values=c("gold3","darkblue")) +
  scale_fill_manual(values=c("gold3","darkblue")) +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  geom_point(data=psi.resu, aes(x=mean,y=-c(0.4,0.2),color=model), size=2) +
  geom_errorbarh(data=psi.resu,
                 aes(xmin=`2.5%`,xmax=`97.5%`,y=-c(0.4,0.2),color=model),
                 linewidth=1,inherit.aes=F) +
  labs(x="Occupancy probability", y="Density",title="White tail deer - De Luca Forest") +
  theme_classic(base_size=14) +
  theme(legend.position="bottom") 

##** All parameters
pars <- as.data.frame(rbind(
  out1$summary[1:4,],
  out2$summary[-3,]
))

pars$par <- factor(c(rownames(out1$summary[1:4,]),rownames(out2$summary[-3,])),
                   levels=rev(c("psi","phi","p","lam","ome")),
                   labels=rev(c("Occurence prob.",
                            "Avaliabiltiy prob.","Overall detection prob.",
                            "True-positive rate",
                            "False-positive rate"))
              )

pars$model <- rep(c("AI-scores (FP-FN)","Manual Std. (FN)"),
                  c(4,2))

fig.pars <- ggplot(data=pars, aes(x=mean,xmin=`2.5%`,xmax=`97.5%`,y=par,color=model)) +
  geom_point(size=3,position=position_dodge(0.4)) +
  geom_errorbarh(linewidth=1,height=0.2,position=position_dodge(0.4)) +
  scale_color_manual(values=c("black","gray60")) + #c("gold3","darkblue")
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
  labs(x="Estimate", y="Parameter") + #title="Whitetails occurrence model"
  theme_bw(base_size=16) +
  theme(legend.position="top",plot.title=element_text(face="bold",size=20),
        legend.text=element_text(size=18)) 

ggsave("Parameter estimates.png",fig.pars,width=22,height=15,unit="cm",dpi=400)




pars["phi","mean"] * pars["lam","mean"]

pars["p","mean"]


x = seq(0,1,length.out=1000)
y = 25 ^ -x

par(mar=c(3,5,2,2))
plot(y~x,cex.lab=2,cex.axis=1.6,type="l",lwd=4,xlab="",ylab="Performance")

legend("topright","Algorithm counts",fill="gray80",bty="n",cex=3)
legend("topright","Abundance estimates",pch=19,bty="n",cex=2.5)



