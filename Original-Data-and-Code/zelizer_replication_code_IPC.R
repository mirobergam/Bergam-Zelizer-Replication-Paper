require(plyr)
require(stringr)
require(ggplot2)
require(boot)
require(grid)
require(RColorBrewer)
require(reshape2)
require(grid)
require(gridExtra)
require(data.table)
require(blockTools)
require(dplyr)
require(Hmisc)

setwd("C:/")

#############################################################################
#######################          Import Data       ##########################
#############################################################################

load("zelizer_replication_data_IPC.RData")

df$t <- relevel(as.factor(df$t) , ref = 2) ### Relevel factors so control is base category; Direct treatment: -1 is advocate; 0 is control; 1 is staffer
df$t.spill <- relevel(as.factor(df$t.spill) , ref = 2) ### Indirect treatment:  -1 is advocate; 0 is control; 1 is staffer; 2 is double staffer
for(i in 1:ncol(t.for.ri)){
  t.for.ri[,i] <- relevel(as.factor(t.for.ri[,i]) , ref = 2)
  t.spill.for.ri[,i] <- relevel(as.factor(t.spill.for.ri[,i]) , ref = 2)
  cat("\r", i, "of", ncol(t.for.ri)) ### How far loop has progressed
}
df$treatment <- interaction(df$t,df$t.spill)
levels(t.for.ri[,1])
levels(t.spill.for.ri[,1])
levels(df$treatment)

#############################################################################
#######################          Summary Stats     ##########################
#############################################################################
### Table 3:
subset(df , office.size > 1) %>%
  ddply(. , .(treatment) , function(x) out = data.frame(
    Weighted.Mean = round(weighted.mean(x$cosp , x$ipw),3) ,
    N = length(x$cosp)))


### Examine unweighted results by study and office size.
subset(df , office.size > 1) %>%
  ddply(. , .(study , office.size , treatment) , function(x) out =
      data.frame(
          Weighted.Mean = round(mean(x$cosp),3) ,
          N = length(x$cosp)))


#############################################################################
##########           Estimated Treatment Effects (Table 4)         ##########
#############################################################################

# Preferred specification
estimates <- data.frame(itt = matrix(NA,nrow=3,ncol=1))
row.names(estimates) <- c("t","t.spill","t+t.spill")

estimates[1:3,] <- summary(lm(cosp ~ treatment + bill + leg ,
                              data = subset(df , office.size > 1),
                              weights=ipw))$coef[
                              c('treatment1.0','treatment0.1','treatment1.1'),
                              'Estimate']
estimates
estimates[2,1]/estimates[1,1] ### Contagion rate
nrow(subset(df , office.size > 1))

### Alternative specifications
# No spillover - create weights as if no spillover
df$ipw2 <- 1 / 0.75 ### For control obs in Study 1
df$ipw2[df$study==1 & df$t!=0] <- 1 / 0.25 ### For treated obs in Study 1
df$ipw2[df$study==2 & df$t==0] <- 1 / (1/3 + 2/3*0.75) ### For control obs in Study 2
df$ipw2[df$study==2 & df$t!=0] <- 1 / (2/3 * 0.25) ### For treated obs in Study 2

summary(lm(cosp ~ t + bill + leg,
           data = subset(df , office.size > 1),
           weights = ipw2))$coef[c('t1'),'Estimate']
# No legislator fixed effects
summary(lm(cosp ~ treatment + bill ,
           data = subset(df , office.size > 1),
           weights=ipw))$coef[
           c('treatment1.0','treatment0.1','treatment1.1'),'Estimate']
# No bill effects
summary(lm(cosp ~ treatment + leg ,
           data = subset(df , office.size > 1),
           weights=ipw))$coef[
           c('treatment1.0','treatment0.1','treatment1.1'),'Estimate']
# No fixed effects
summary(lm(cosp ~ treatment ,
           data = subset(df , office.size > 1),
           weights=ipw))$coef[
           c('treatment1.0','treatment0.1','treatment1.1'),'Estimate']

#######################               RI           ##########################

ri.t <- data.frame(matrix(NA, nrow = 3 , ncol = ncol(t.for.ri)))
row.names(ri.t) <- row.names(estimates)
ci.t <- ri.t
tempdf <- df

tempdf$Y0 <- df$cosp
tempdf$Y0[tempdf$treatment=="1.0"] <- tempdf$cosp[tempdf$treatment=="1.0"] -
                                            estimates[1,1]
tempdf$Y0[tempdf$treatment=="0.1"] <- tempdf$cosp[tempdf$treatment=="0.1"] -
                                            estimates[2,1]
tempdf$Y0[tempdf$treatment=="1.1"] <- tempdf$cosp[tempdf$treatment=="1.1"] -
                                            estimates[3,1]

for(i in 1:ncol(ri.t)){
  tempdf$t <- t.for.ri[,i]
  tempdf$t.spill <- t.spill.for.ri[,i]
  tempdf$ipw <- ipw.for.ri[,i]
  tempdf$treatment <- interaction(tempdf$t,tempdf$t.spill)

  ri.t[1:3,i] <- summary(lm(cosp ~ treatment + bill + leg ,
                            data = subset(tempdf , office.size > 1),
                            weights=ipw))$coef[
                            c('treatment1.0','treatment0.1','treatment1.1'),
                            'Estimate']

  tempdf$cosp1 <- tempdf$Y0
  tempdf$cosp1[tempdf$treatment=="1.0"] <- tempdf$Y0[tempdf$treatment=="1.0"] +
                                                estimates[1,1]
  tempdf$cosp1[tempdf$treatment=="0.1"] <- tempdf$Y0[tempdf$treatment=="0.1"] +
                                                estimates[2,1]
  tempdf$cosp1[tempdf$treatment=="1.1"] <- tempdf$Y0[tempdf$treatment=="1.1"] +
                                                estimates[3,1]

  ci.t[1:3,i] <- summary(lm(cosp1 ~ treatment + bill + leg ,
                            data = subset(tempdf , office.size > 1),
                            weights=ipw))$coef[
                            c('treatment1.0','treatment0.1','treatment1.1'),
                            'Estimate']
  cat("\r", i, "of", ncol(ri.t)) ### How far loop has progressed
}

#### Table 4:
apply(cbind.data.frame(ri.t[1:3,], estimates) , 1 , function(x)
  out = unlist(data.frame(
    Est = round(x[ncol(ri.t)+1],3),
#    mean.sim.ITT = round(mean(x[1:ncol(ri.t)]),4),
    SE = round(sd(x[1:ncol(ri.t)]),3)
#    p.one.sided = round(mean(x[1:ncol(ri.t)] > x[ncol(ri.t)+1]),4)
  )))

#### Test of interaction effect vs. null of additive effect:
mean(apply(matrix(1:10000 , ncol = 1) , 1 , function(x)
  out = unlist(data.frame(
    p.one.sided = as.numeric(estimates[3,1] < ci.t[1,x] + ci.t[2,x])
  ))))

#############################################################################
###############     Heterogeneous effects by bill progress  #################
#############################################################################

estimates <- data.frame(itt = matrix(NA,nrow=3,ncol=2))
row.names(estimates) <- c("t","t.spill","t+t.spill")
colnames(estimates) <- c("no_vote","vote")

estimates[1:3,1] <- unlist(summary(lm(cosp ~ treatment + bill + leg,
                                      data = subset(df , office.size > 1 & passed==0),
                                      weights=ipw))$coef[c('treatment1.0',
                                                           'treatment0.1',
                                                           'treatment1.1'),
                                                           'Estimate'])


estimates[1:3,2] <- unlist(summary(lm(cosp ~ treatment + bill + leg,
                                      data = subset(df , office.size > 1 & passed==1),
                                      weights=ipw))$coef[c('treatment1.0',
                                                           'treatment0.1',
                                                           'treatment1.1'),
                                                           'Estimate'])

estimates

#######################               RI           ##########################

ri.t3 <- array(NA, dim = c(3 , 2, ncol(t.spill.for.ri)),
               dimnames = list(
                 c("t","t.spill","t+t.spill"),
                 c("no_vote","vote"),
                 c(1:10000)))

for(i in 1:ncol(t.spill.for.ri)){
  tempdf <- df[,c("bill","cosp","office.size","leg","passed")]
  tempdf$t <- t.for.ri[,i]
  tempdf$t.spill <- t.spill.for.ri[,i]
  tempdf$ipw <- ipw.for.ri[,i]

  tempdf$treatment <- interaction(tempdf$t,tempdf$t.spill)
  ri.t3[1:3,1,i] <- unlist(summary(lm(cosp ~ treatment + bill + leg,
                                      data = subset(tempdf , office.size > 1 & passed==0),
                                      weights=ipw))$coef[c('treatment1.0',
                                                           'treatment0.1',
                                                           'treatment1.1'),
                                                           'Estimate'])

  ri.t3[1:3,2,i] <- unlist(summary(lm(cosp ~ treatment + bill + leg,
                                      data = subset(tempdf , office.size > 1 & passed==1),
                                      weights=ipw))$coef[c('treatment1.0',
                                                           'treatment0.1',
                                                           'treatment1.1'),
                                                           'Estimate'])

  cat("\r", i, "of", ncol(t.spill.for.ri)) ### How far loop has progressed
}

estimates
### Table 5:
for(j in 1:2){
  print(colnames(estimates)[j])
  for(i in 1:3){
    print(row.names(estimates)[i])
    print(data.frame(
#      Mean = mean(ri.t3[i,j,]) ,
      Est = estimates[i,j] ,
      SE = sd(ri.t3[i,j,])
#      pone = mean(ri.t3[i,j,] > estimates[i,j])
    ))
  }}

mean(ri.t3[2,2,] - ri.t3[2,1,] > 0.0513)

#############################################################################
#######################        Alternative Models  ##########################
#############################################################################

### This function creates simulation-based inverse probability weights.
ipw.function <- function(DATASET , T , T.SPILL , RI.T , RI.T.SPILL , VARNAME){
	temp <- apply(matrix(1:ncol(RI.T.SPILL) , nrow=1) , 2 , function(z) {
		temp1 <- as.numeric(RI.T.SPILL[,z]==DATASET[,T.SPILL])
		temp2 <- as.numeric(RI.T[,z]==DATASET[,T])
		return(as.numeric(temp1+temp2==2))
	})

	DATASET$V1 <- 1/apply(temp , 1 , mean)
	setnames(DATASET,"V1",VARNAME)
	return(DATASET)
}

load("replication_data_alternative_spillover_simulations.RData")
### Add IPWs for observed assignments
df <- ipw.function(df , "t" , "t.spill.desk" , t.for.ri , ri.t.spill.desk , "ipw.desk")
df <- ipw.function(df , "t" , "t.spill.dist" , t.for.ri , ri.t.spill.dist , "ipw.dist")
df <- ipw.function(df , "t" , "t.spill.ideo" , t.for.ri , ri.t.spill.ideo , "ipw.ideo")

table(df$t.spill,useNA="always")
table(df$t.spill.desk,useNA="always")
table(df$t.spill.dist,useNA="always")
table(df$t.spill.ideo,useNA="always")

df$treatment.desk <- relevel(interaction(df$t,df$t.spill.desk) , ref = 4)
df$treatment.dist <- relevel(interaction(df$t,df$t.spill.dist) , ref = 4)
df$treatment.ideo <- relevel(interaction(df$t,df$t.spill.ideo) , ref = 4)

estimates <- data.frame(itt = matrix(NA,nrow=3,ncol=4))
row.names(estimates) <- c("t","t.spill","t+t.spill")
colnames(estimates) <- c("office","desk","distance","ideology")

### Office Spillover
estimates[1:3,1] <- unlist( summary(lm(cosp ~ treatment + bill + leg,
                      data = subset(df , office.size > 1 &
                                      !is.na(t.spill.desk) &
                                      !is.na(t.spill.dist) &
                                      !is.na(t.spill.ideo)),
                      weights=ipw))$coef[c('treatment1.0',
                                           'treatment0.1',
                                           'treatment1.1'),
                                           'Estimate'] )

### Desk Spillover
estimates[1:3,2] <- unlist( summary(lm(cosp ~ treatment.desk + bill + leg,
                      data = subset(df , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw.desk))$coef[c('treatment.desk1.0',
                                                'treatment.desk0.1',
                                                'treatment.desk1.1'),
                                                'Estimate'] )

 ### Distance Spillover
estimates[1:3,3] <- unlist( summary(lm(cosp ~ treatment.dist + bill + leg,
                      data = subset(df , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw.dist))$coef[c('treatment.dist1.0',
                                                'treatment.dist0.1',
                                                'treatment.dist1.1'),
                                                'Estimate'] )

### Ideology Spillover
estimates[1:3,4] <- unlist( summary(lm(cosp ~ treatment.ideo + bill + leg,
                      data = subset(df , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw.ideo))$coef[c('treatment.ideo1.0',
                                                'treatment.ideo0.1',
                                                'treatment.ideo1.1'),
                                                'Estimate'] )

estimates

#######################               RI           ##########################

ri.t2 <- array(NA, dim = c(3 , 4, ncol(ri.t.spill.desk)),
               dimnames = list(
                 c("t","t.spill","t+t.spill"),
                 c("office","desk","distance","ideology"),
                 c(1:1000)))

for(i in 1:ncol(ri.t2)){
  tempdf <- df[,c("bill","cosp","office.size","leg")]
  tempdf$t <- t.for.ri[,i]
  tempdf$t.spill <- t.spill.for.ri[,i]
  tempdf$ipw <- ipw.for.ri[,i]

  tempdf$t.spill.desk <- ri.t.spill.desk[,i]
  tempdf$t.spill.dist <- ri.t.spill.dist[,i]
  tempdf$t.spill.ideo <- ri.t.spill.ideo[,i]


  tempdf <- ipw.function(tempdf , "t" , "t.spill.desk" , t.for.ri , ri.t.spill.desk , "ipw.desk")
  tempdf <- ipw.function(tempdf , "t" , "t.spill.dist" , t.for.ri , ri.t.spill.dist , "ipw.dist")
  tempdf <- ipw.function(tempdf , "t" , "t.spill.ideo" , t.for.ri , ri.t.spill.ideo , "ipw.ideo")
  tempdf$treatment <- interaction(tempdf$t,tempdf$t.spill)
  tempdf$treatment.desk <- interaction(tempdf$t,tempdf$t.spill.desk)
  tempdf$treatment.dist <- interaction(tempdf$t,tempdf$t.spill.dist)
  tempdf$treatment.ideo <- interaction(tempdf$t,tempdf$t.spill.ideo)
  print(levels(tempdf$treatment.ideo)[1])
  print(levels(tempdf$treatment.desk)[1])
  print(levels(tempdf$treatment.dist)[1])

  ### Office Spillover
  ri.t2[1:3,1,i] <- unlist( summary(lm(cosp ~ treatment + bill + leg,
                      data = subset(tempdf , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw))$coef[c('treatment1.0',
                                           'treatment0.1',
                                           'treatment1.1'),
                                           'Estimate'] )

  ### Desk Spillover
  ri.t2[1:3,2,i] <- unlist( summary(lm(cosp ~ treatment.desk + bill + leg,
                      data = subset(tempdf , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw.desk))$coef[c('treatment.desk1.0',
                                                'treatment.desk0.1',
                                                'treatment.desk1.1'),'Estimate'] )

  ### Distance Spillover
  ri.t2[1:3,3,i] <- unlist( summary(lm(cosp ~ treatment.dist + bill + leg,
                      data = subset(tempdf , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw.dist))$coef[c('treatment.dist1.0',
                                                'treatment.dist0.1',
                                                'treatment.dist1.1'),
                                                'Estimate'] )

  ### Ideology Spillover
  ri.t2[1:3,4,i] <- unlist( summary(lm(cosp ~ treatment.ideo + bill + leg,
                      data = subset(tempdf , office.size > 1 &
                                     !is.na(t.spill.desk) &
                                     !is.na(t.spill.dist) &
                                     !is.na(t.spill.ideo)),
                      weights=ipw.ideo))$coef[c('treatment.ideo1.0',
                                                'treatment.ideo0.1',
                                                'treatment.ideo1.1'),
                                                'Estimate'] )

  cat("\r", i, "of", ncol(ri.t.spill.desk)) ### How far loop has progressed
}

estimates
### Table 6:
for(j in 1:4){
  print(colnames(estimates)[j])
  for(i in 1:3){
  print(row.names(estimates)[i])
    print(data.frame(
#  	Mean = mean(ri.t2[i,j,]) ,
  	Est = estimates[i,j] ,
  	SE = sd(ri.t2[i,j,])
#	  pone = mean(ri.t2[i,j,] > estimates[i,j])
    ))
  }}

tempdf2 <- subset(tempdf , office.size > 1 &
                    !is.na(t.spill.desk) & !is.na(t.spill.dist) & !is.na(t.spill.ideo))
### Note lack of restriction in alternative models leads to fewer spillover assignments
ddply(tempdf2 , .(treatment) , nrow)
ddply(tempdf2 , .(treatment.desk) , nrow)
ddply(tempdf2 , .(treatment.dist) , nrow)
ddply(tempdf2 , .(treatment.ideo) , nrow)

### Big weights due to the possibility of interacting these low probability treatments.
range(tempdf2$ipw)
range(tempdf2$ipw.desk)
range(tempdf2$ipw.dist)
range(tempdf2$ipw.ideo)