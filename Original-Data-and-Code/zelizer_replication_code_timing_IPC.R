require(plyr)
require(stringr)
require(ggplot2)
require(data.table)
require(dplyr)

setwd("C:/")
load("zelizer_replication_data_dates_IPC.RData")

table(df$cosp , !is.na(df$days)) ### Note two cosp withdrawals; six missing dates.

nsims <- 1000

df$treatment <- interaction(df$t,df$t.spill)
levels(df$treatment)

nDays <- seq(0,98,7)
days.ri <- data.frame(matrix(NA,nrow=length(nDays) , ncol = 3))
days.ri.SE <- data.frame(matrix(NA,nrow=length(nDays) , ncol = 3))

for(i in nDays){
  Ydays <- as.numeric(!is.na(df$days) & df$days <= i & df$days >=0 & df$cosp==1)
  days.ri[match(i,nDays),1:3] <- summary(lm(Ydays ~ treatment + bill + leg ,
                                            data = df,
                                            weights=ipw))$coef[c('treatment1.0','treatment0.1','treatment1.1'),'Estimate']
  temp <- data.frame(matrix(NA,nrow=nsims,ncol=3))
  for(j in 1:nsims){
    treatment.sim <- interaction(t.for.ri2[,j],t.spill.for.ri2[,j])
    levels(treatment.sim)
    temp[j,1:3] <- summary(lm(Ydays ~ treatment.sim + bill + leg ,
                              data = df,
                              weights=ipw.for.ri2[,j]))$coef[c('treatment.sim1.0','treatment.sim0.1','treatment.sim1.1'),'Estimate']
  }
  temp <- apply(temp , 2 , sd)
  days.ri.SE[match(i,nDays),1:3] <- temp
  cat("\r", i, "of", nDays) ### How far loop has progressed
}

days.ri
days.ri.SE


pdf("redefine_DV_direct_treatment.pdf",height=3,width=6,onefile=FALSE)
tempnum <- 1
data.frame(ATE = days.ri[,tempnum],
           CIlow = days.ri[,tempnum] - 1.65*days.ri.SE[,tempnum],
           CIhigh = days.ri[,tempnum] + 1.65*days.ri.SE[,tempnum],
           X = nDays/7) %>%
  ggplot(., aes(x = factor(X) )) +
  geom_errorbar(aes(ymin = CIlow*100 , ymax = CIhigh*100 , width = 0)) +
  geom_point(aes(y = ATE*100)) +
  scale_y_continuous(minor_breaks=waiver(), labels = waiver(), breaks=c(-2,0,2,4,6,8) , limits = c(-2,9)) +
  labs(x = "Week of Session", y = "Est. ATE (perc. pts.)") +
  geom_smooth(aes(y = ATE*100) , col = "black" , method = "loess" , se = FALSE) +
  geom_vline(xintercept=0 , color = "#535353" , lwd = 1) +
  geom_hline(yintercept=0 , color = "#535353" , lwd = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(color = "#535353" , size = 12, face="bold"),
        axis.title.y = element_text(color = "#535353" , size = 12, face="bold"),
        axis.text.y = element_text(size = 12, colour = "grey50"),
        axis.text.x = element_text(size = 8, colour = "grey50"),
        panel.border = element_blank(),
        legend.key = element_rect(colour = "white",fill = "white"),
        legend.text = element_text(size=10,color="grey50"),
        legend.title = element_text(size = 11, face = "bold", colour = "#535353"),
        strip.text.x = element_text(size = 11, face = "bold", colour = "#535353"),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.position="right",
        plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
dev.off()

pdf("redefine_DV_indirect_treatment.pdf",height=3,width=6,onefile=FALSE)
tempnum <- 2
data.frame(ATE = days.ri[,tempnum],
           CIlow = days.ri[,tempnum] - 1.65*days.ri.SE[,tempnum],
           CIhigh = days.ri[,tempnum] + 1.65*days.ri.SE[,tempnum],
           X = nDays/7) %>%
  ggplot(., aes(x = factor(X) )) +
  geom_errorbar(aes(ymin = CIlow*100 , ymax = CIhigh*100 , width = 0)) +
  geom_point(aes(y = ATE*100)) +
  scale_y_continuous(minor_breaks=waiver(), labels = waiver(), breaks=c(-2,0,2,4,6,8) , limits = c(-2,9)) +
  labs(x = "Week of Session", y = "Est. ATE (perc. pts.)") +
  geom_smooth(aes(y = ATE*100) , col = "black" , method = "loess" , se = FALSE) +
  geom_vline(xintercept=0 , color = "#535353" , lwd = 1) +
  geom_hline(yintercept=0 , color = "#535353" , lwd = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(color = "#535353" , size = 12, face="bold"),
        axis.title.y = element_text(color = "#535353" , size = 12, face="bold"),
        axis.text.y = element_text(size = 12, colour = "grey50"),
        axis.text.x = element_text(size = 8, colour = "grey50"),
        panel.border = element_blank(),
        legend.key = element_rect(colour = "white",fill = "white"),
        legend.text = element_text(size=10,color="grey50"),
        legend.title = element_text(size = 11, face = "bold", colour = "#535353"),
        strip.text.x = element_text(size = 11, face = "bold", colour = "#535353"),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.position="right",
        plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
dev.off()

pdf("redefine_DV_combined_treatment.pdf",height=3,width=6,onefile=FALSE)
tempnum <- 3
data.frame(ATE = days.ri[,tempnum],
           CIlow = days.ri[,tempnum] - 1.65*days.ri.SE[,tempnum],
           CIhigh = days.ri[,tempnum] + 1.65*days.ri.SE[,tempnum],
           X = nDays/7) %>%
  ggplot(., aes(x = factor(X) )) +
  geom_errorbar(aes(ymin = CIlow*100 , ymax = CIhigh*100 , width = 0)) +
  geom_point(aes(y = ATE*100)) +
  labs(x = "Week of Session", y = "Est. ATE (perc. pts.)") +
  geom_smooth(aes(y = ATE*100) , col = "black" , method = "loess" , se = FALSE) +
  geom_vline(xintercept=0 , color = "#535353" , lwd = 1) +
  geom_hline(yintercept=0 , color = "#535353" , lwd = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(color = "#535353" , size = 12, face="bold"),
        axis.title.y = element_text(color = "#535353" , size = 12, face="bold"),
        axis.text.y = element_text(size = 12, colour = "grey50"),
        axis.text.x = element_text(size = 8, colour = "grey50"),
        panel.border = element_blank(),
        legend.key = element_rect(colour = "white",fill = "white"),
        legend.text = element_text(size=10,color="grey50"),
        legend.title = element_text(size = 11, face = "bold", colour = "#535353"),
        strip.text.x = element_text(size = 11, face = "bold", colour = "#535353"),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.position="right",
        plot.margin = unit(c(0.75,0.75,0.75,0.75), "cm"))
dev.off()