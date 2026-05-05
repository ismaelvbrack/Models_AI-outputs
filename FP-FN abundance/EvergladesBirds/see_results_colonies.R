
library(jagsUI)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggbreak)

# Import data --------------------------------------------------------------------------------

counts <- read.table(file.path("EvergladesBirds","counts_12colonies_all-dates.txt"),h=T,check.names=F)

counts$date <- as.Date(counts$date)

counts <- counts[-which(counts$site=="Hidden"),]

scr.data <- read.csv(file.path("EvergladesBirds","AI-detections_12colonies_all-dates.csv"))

scr.data <- scr.data[-which(scr.data$site=="Hidden"),]

table(scr.data$det.type, scr.data$site)

# Flights in colonies through time
ggplot(data=counts,
       aes(x=date, y=site, col=ifelse(TP>0 | FP>0,"HQ",NA))) +
  geom_point(size=3) +
  theme_bw(base_size=14) +
  theme(legend.position="none")

#resu1 <- readRDS(file.path("EvergladesBirds","resu_model2.rds"))

#resu2.1 <- readRDS(file.path("EvergladesBirds","resu_model2_predictions.rds"))
#resu2.2 <- readRDS(file.path("EvergladesBirds","resu2_model2_predictions.rds"))
resu2 <- readRDS(file.path("EvergladesBirds","resu_model2_combined flights.rds"))

# plot(resu1$mean$N ~ resu2$mean$N[,1])
# abline(a=0,b=1)

#plot(resu1$mean$precision ~ resu2$mean$precision)
#abline(a=0, b=1)

#resu1$summary[grep("mu.scr",rownames(resu1$summary)), c("mean","2.5%","97.5%")]
resu2$summary[grep("mu.scr",rownames(resu2$summary)), c("mean","2.5%","97.5%")]

# Abundance per colony --------------------------------------------------------------------------------------------
# counts1 <- counts[which(!is.na(counts$TP)),] # sites with verification data
# 
# aqui <- which(paste0(scr.data$site,scr.data$date) %in%
#                 paste0(counts1$site,counts1$date))
# 
# scr.data1 <- scr.data[aqui,] # sites with verification data
# 
# table(scr.data1[which(scr.data1$score >= 0.5 &
#                         is.na(scr.data1$det.type)),"site"])
# 
# 
# tab1 <- data.frame(
#   colony=counts1$site,
#   count=counts1$unverified + counts1$TP,
#   count50=as.numeric(table(scr.data1[which(scr.data1$score >= 0.5 &
#                                   is.na(scr.data1$det.type)),"site"])),
#   mean=resu1$mean$N,
#   lcl=resu1$q2.5$N,
#   ucl=resu1$q97.5$N
# )
# 
# ggplot(data=tab1, 
#        aes(x=colony,y=mean,ymin=lcl,ymax=ucl)) +
#   geom_col(aes(x=colony,y=count),fill="gray80",width=0.3) +
#   #geom_col(aes(x=colony,y=count50),fill="gray20",width=0.3) +
#   geom_point(size=2) +
#   geom_errorbar(linewidth=1,width=0.2) +
#   labs(y="Population size",x="Colony",
#        title="Spring / 2023") +
#   theme_classic(base_size=14) +
#   theme(axis.text.x=element_text(angle=90))
# 
# ggplot(data=tab1, 
#        aes(x=count,y=mean,ymin=lcl,ymax=ucl)) +
#   geom_point(size=2) +
#   geom_errorbar(linewidth=1,width=0.2) +
#   labs(y="Estimated abundance",x="Observed counts",
#        title="Spring / 2023") +
#   theme_classic(base_size=14) +
#   theme(axis.text.x=element_text(angle=90))
#   

# Abundance per colony combined approach ------------------------------------------------------------------------

sites <- unique(counts$site)
counts2 <- counts

counts2$meanN <- counts2$lclN <- counts2$uclN <- NA

for(s in 1:length(sites)){
  aqui <- which(counts2$site ==  sites[s])
  
  counts2$meanN[aqui] <- na.omit(resu2$mean$N[s,])
  counts2$lclN[aqui] <- na.omit(resu2$q2.5$N[s,])
  counts2$uclN[aqui] <- na.omit(resu2$q97.5$N[s,])
}

tabs <- data.frame(
  site=counts$site,
  date=as.Date(counts$date),
  verif=ifelse(is.na(counts$TP), "only-unverif", "with-verif"),
  mean=counts2$meanN,
  lcl=counts2$lclN,
  ucl=counts2$uclN,
  count=counts$unverified + ifelse(is.na(counts$TP),0,counts$TP)
)

abund2 <- ggplot(data=tabs,
       aes(x=date, y=mean, ymin=lcl, ymax=ucl, col=verif)) +
  geom_point(size=1.5) +
  geom_errorbar(width=1,linewidth=0.8) +
  scale_color_manual(values=c("gray10", "red3")) +
  coord_cartesian(ylim=c(0,NA)) +
  labs(y="Colony abundance",x="") +
  facet_wrap(.~site,scales="free_y",nrow=11) +
  theme_bw(base_size=14) +
  theme(legend.position="none", axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=14)) 

ggsave(file.path("abundance_11colonies_multiple-dates.png"),abund2,
       width=20,height=30,unit="cm",dpi=400)


#* Estimates ~ Counts
ggplot(data=tabs,
       aes(x=count, y=mean, ymin=lcl, ymax=ucl, col=verif)) +
  geom_point(size=2) +
  geom_errorbar(width=50,linewidth=0.8) +
  scale_color_manual(values=c("gray10", "red3")) +
  coord_cartesian(xlim=c(0,3000),ylim=c(0,3000)) +
  geom_abline(intercept=0,slope=1,col="gray40",linewidth=1,linetype="dashed") +
  labs(y="Estimated abundance",x="Observed counts") +
  theme_bw(base_size=14) +
  theme(legend.position="none", axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14)) 

tabs[which(tabs$verif=='with-verif'), "count"] / tabs[which(tabs$verif=='with-verif'),"mean"]

# Figure Abundance vs Counts & Precision ----------------------------------------------------------------

fig.abund <- ggplot(data=tabs[which(tabs$verif=='with-verif'),], 
                    aes(x=site,y=mean,ymin=lcl,ymax=ucl)) +
  geom_col(aes(x=site,y=count),fill="gray80",width=0.3) +
  geom_point(size=2) +
  geom_errorbar(linewidth=1,width=0.2) +
  labs(y="Colony abundance",x="Bird colony",title="b)") +
  
  scale_y_break(c(1200,6000),space=0.5) +
  
  scale_y_continuous(breaks=seq(0,7000,500), sec.axis = dup_axis(breaks = NULL)) +
  #title="Everglades birds abundance model") +
  theme_classic(base_size=14) +
  theme(axis.text.x=element_text(angle=90),
        axis.title.y=element_text(size=14),
        plot.title=element_text(face="bold",size=16))

fig.b <- ggplot(data=tabs[which(tabs$site == "6thBridge" |
                                  tabs$site == "CypressCity"),],
                 aes(x=date, y=mean, ymin=lcl, ymax=ucl, col=verif)) +
  geom_point(size=2) +
  geom_errorbar(width=1,linewidth=0.8) +
  scale_color_manual(values=c("gray10", "red2")) +
  coord_cartesian(ylim=c(0,NA)) +
  labs(y="Colony abundance",x="",,title="c)") +
  facet_wrap(.~site,scales="free_y",nrow=11) +
  theme_bw(base_size=14) +
  theme(legend.position="none", axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=14),
        plot.title=element_text(face="bold",size=16))
# 
# grid.arrange(arrangeGrob(fig.a, fig.b, fig.c,
#                          ncol = 1,
#                          heights = c(6, 2, 2)) ,
#              left=textGrob("Colony abundance",rot=90,
#                            gp = gpar(fontface = "bold", fontsize = 12))
# )


# ggsave(file.path("EvergladesBirds","abundance-counts_11colonies.svg"),fig.abund,
#        width=22,height=15,unit="cm",dpi=400)
# ggsave(file.path("EvergladesBirds","abundance_6thBridge-CypressCity.svg"),fig.b,
#        width=22,height=10,unit="cm",dpi=400)

# Abundance per colony fixed precision ----------------------------------------------------------------------------
# 
# tabs <- data.frame(
#   site=counts$site,
#   date=as.Date(counts$date),
#   verif=ifelse(is.na(counts$TP), "only-unverif", "with-verif"),
#   mean=c(resu1$mean$N, resu2.1$mean$N),
#   lcl=c(resu1$q2.5$N, resu2.1$q2.5$N),
#   ucl=c(resu1$q97.5$N, resu2.1$q97.5$N),
#   count=counts$unverified + ifelse(is.na(counts$TP),0,counts$TP)
# )
# 
# 
# abund1 <- ggplot(data=tabs,
#        aes(x=date, y=mean, ymin=lcl, ymax=ucl, col=verif)) +
#   geom_point(size=2) +
#   geom_errorbar(width=1,linewidth=0.8) +
#   scale_color_manual(values=c("gray10", "red3")) +
#   coord_cartesian(ylim=c(0,NA)) +
#   labs(y="Colony abundance",x="") +
#   facet_wrap(.~site,scales="free_y",nrow=11) +
#   theme_bw(base_size=14) +
#   theme(legend.position="none", axis.text.y=element_text(size=10),
#         axis.text.x=element_text(size=14)) 


# Compare precision -----------------------------------------------------------------------------------------------

# prec.tab <- data.frame(
#   site=rep(unique(counts$site), 2),
#   mean=c(resu1$mean$precision, resu2$mean$precision),
#   lcl=c(resu1$q2.5$precision,  resu2$q2.5$precision),
#   ucl=c(resu1$q97.5$precision, resu2$q97.5$precision),
#   model=rep(c("fixed prec","combined"),each=11)
# )
# 
# ggplot(data=prec.tab,
#        aes(x=model, y=mean, ymin=lcl, ymax=ucl)) +
#   geom_point(size=3) +
#   geom_errorbar(width=0.1,linewidth=0.8) +
#   coord_cartesian(ylim=c(0,1)) +
#   labs(y="Precision",x="Model") +
#   scale_x_discrete(limits = c("fixed prec","combined")) +
#   facet_wrap(.~site,nrow=11) +
#   theme_bw(base_size=14) +
#   theme(legend.position="none", axis.text.y=element_text(size=10),
#         axis.text.x=element_text(size=14)) 
#   
# 
# plot(prec.tab[which(prec.tab$model=="combined"),"mean"] ~ log(tabs$mean[1:11]),
#      xlab="ln(Abundance)",ylab=("Precision"), pch=19)


# Composed figure: Precision and Abundance ------------------------------------------------------------------

prec.tab <- data.frame(
  site=unique(counts$site),
  mean=resu2$mean$precision,
  lcl=resu2$q2.5$precision,
  ucl=resu2$q97.5$precision,
  model=rep("combined",11)
)

prec.order <- order(prec.tab[which(prec.tab$model=="combined"),"mean"])
fig.abund <- fig.abund + scale_x_discrete(limits=prec.tab$site[prec.order])

fig.prec <- ggplot(data=prec.tab[which(prec.tab$model=="combined"),],
                   aes(x=site,y=mean,ymin=lcl,ymax=ucl)) +
  geom_point(size=2) +
  geom_errorbar(width=0.3,linewidth=1) +
  labs(x="",y="Estimated algorithm precision", title="a)") +
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
  theme_classic(base_size=14) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        plot.title=element_text(face="bold",size=16)) +
scale_x_discrete(limits=prec.tab$site[prec.order])

# tabs22 <- merge(tabs[which(tabs$verif=='with-verif'),],
#                 prec.tab[which(prec.tab$model=="combined"),], by="site")

# # Difference between counts and abundance vs. precision
# fig.prec <- ggplot(data=tabs22, aes(y=(count-mean.x)/count,
#                                     ymin=ifelse((count-lcl.x)/count>0,(count-lcl.x)/count,0),
#                                     ymax=ifelse((count-ucl.x)/count>0,(count-ucl.x)/count,0),
#                         x=mean.y,xmin=lcl.y,xmax=ucl.y)) +
#   geom_errorbar(width=0.02,linewidth=1,color="gray60") +
#   geom_errorbarh(height=0.02,linewidth=1,color="gray60") +
#   geom_point(size=2) +
#   #geom_text(label=tabs22$site) +
#   labs(y="(Counts - Abundance) / Counts",x="Algorithm precision",title="b)") +
#   scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
#   scale_x_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
#   theme_classic(base_size=14) +
#   theme(plot.title=element_text(face="bold",size=16))
# 
# # Abundance vs. Prcision
# fig.abundprec <- ggplot(data=tabs22, aes(x=log(mean.x),xmin=ifelse(log(lcl.x)>0,log(lcl.x),0),xmax=log(ucl.x),
#                                          y=mean.y,ymin=lcl.y,ymax=ucl.y)) +
#   geom_errorbar(width=0.1,linewidth=1,color="gray60") +
#   geom_errorbarh(height=0.02,linewidth=1,color="gray60") +
#   geom_point(size=2) +
#   #geom_text(label=tabs22$site) +
#   labs(x="ln(Abundance)",y="Algorithm precision",title="c)") +
#   scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
#   scale_x_continuous(breaks=seq(0,8,2)) +
#   theme_classic(base_size=14) +
#   theme(plot.title=element_text(face="bold",size=16))
# 
# 
# grid.arrange(fig.abund,fig.prec,fig.abundprec)
# ggsave(file.path("EvergladesBirds","Figb_counts-abund_precision.svg"),fig.prec,
#        width=20,height=12,unit="cm",dpi=400)
# ggsave(file.path("EvergladesBirds","Figa_abundance.svg"),fig.abund,
#        width=20,height=15,unit="cm",dpi=400)
# ggsave(file.path("EvergladesBirds","Figc_prec~abund.svg"),fig.abundprec,
#        width=20,height=12,unit="cm",dpi=400)

grid.arrange(fig.prec,fig.abund,fig.b)

ggsave(file.path("EvergladesBirds","Fig5a_precision.svg"),fig.prec,
       width=20,height=12,unit="cm",dpi=400)
ggsave(file.path("EvergladesBirds","Fig5b_abundance.svg"),fig.abund,
       width=20,height=15,unit="cm",dpi=400)
ggsave(file.path("EvergladesBirds","Fig5c_time-series.svg"),fig.b,
       width=20,height=12,unit="cm",dpi=400)





