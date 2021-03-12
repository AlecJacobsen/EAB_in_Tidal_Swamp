#### Ash Research Part 2 ####
rm(list = ls())
setwd("~/Documents")
library(clusrank)
library(reshape2)
library(plotrix)
library(vegan)
library(tidyr)
library(tidyverse)
library(lme4)
library(jtools)
library(sjPlot)
library(ggplot2)
library(dplyr)

all.params <- read.csv(file="Jug.Bay.Ash.Reasearch.Organized.Datas.only-updated.csv")
all.trees <- read.csv(file="Jug.Bay.Ash.Reasearch.Organized.Datas.only-all.trees.csv")

all.trees$Canopy[is.na(all.trees$Canopy)] <- 100
all.trees$BA <- ((all.trees$dbh..cm./2)^2)*pi 
all.trees$smith.index <- ifelse(all.trees$Canopy >= 99, 1,
                                ifelse(all.trees$Canopy >= 71 & all.trees$Canopy <= 98, 2, 
                                       ifelse(all.trees$Canopy >= 51 & all.trees$Canopy <= 70, 3,
                                              ifelse(all.trees$Canopy >= 10 & all.trees$Canopy <= 50, 4,
                                                     ifelse(all.trees$Canopy >= 0 & all.trees$Canopy <= 9, 5, NA)))))

all.trees$simple <- ifelse(all.trees$spp == "green" & all.trees$Canopy >= 90, "Healthy",
                           ifelse(all.trees$spp == "green" & all.trees$Canopy >= 1 & all.trees$Canopy <= 89, "Infested",
                                  ifelse(all.trees$spp == "green" & all.trees$Canopy == 0, "Dead", "Non-ash")))

all.ash <- subset(all.trees, spp == "green")



################ Tree Params ################
#### Density of Ash ####
clusWilcox.test(all.params$ash.density, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -1.8631, p-value = 0.06245
mean(aggregate(ash.density ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 325.9259
std.error(aggregate(ash.density ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 103.7037
mean(aggregate(ash.density ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 696.2963
std.error(aggregate(ash.density ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 147.9258


# no clustering
shapiro.test(all.params$ash.density.forest)
shapiro.test(all.params$ash.density.marsh)
# both normal 
var.test(all.params$ash.density.forest,all.params$ash.density.marsh)
# variances not equal 
wilcox.test(all.params$ash.density.forest,all.params$ash.density.marsh)
boxplot(all.params$ash.density.forest,all.params$ash.density.marsh)

mean(all.params$ash.density.forest)
std.error(all.params$ash.density.forest)
mean(all.params$ash.density.marsh, na.rm = TRUE)
std.error(all.params$ash.density.marsh, na.rm = TRUE)
# p-value = 0.002286 - significant difference in distributions


#### Stand Density ####
clusWilcox.test(all.params$stand.density, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -1.732, p-value = 0.08327
mean(aggregate(stand.density ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 582.716
std.error(aggregate(stand.density ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 172.8395
mean(aggregate(stand.density ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 943.7037
std.error(aggregate(stand.density ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 162.5179



# no clustering 
shapiro.test(all.params$stand.density.forest)
# not normal
shapiro.test(all.params$stand.density.marsh)
# normal
var.test(all.params$stand.density.forest, all.params$stand.density.marsh)
# Equal Variance
wilcox.test(all.params$stand.density.forest, all.params$stand.density.marsh)
boxplot(all.params$stand.density.forest, all.params$stand.density.marsh)
# p-value = 0.00735 - Marsh greater

mean(all.params$stand.density.forest)
std.error(all.params$stand.density.forest)
mean(all.params$stand.density.marsh, na.rm = TRUE)
std.error(all.params$stand.density.marsh, na.rm = TRUE)
#### Relative Density of Ash ####
clusWilcox.test(all.params$relative.density, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -2.1925, p-value = 0.02834
mean(aggregate(relative.density ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.5767047
std.error(aggregate(relative.density ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.04159341
mean(aggregate(relative.density ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.7278109
std.error(aggregate(relative.density ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 0.03968046


#### Ash BA ####
clusWilcox.test(all.params$ash.BA, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -0.57599, p-value = 0.5646
mean(aggregate(ash.BA ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 221068.6
std.error(aggregate(ash.BA ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 66007.64
mean(aggregate(ash.BA ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 262645.7
std.error(aggregate(ash.BA ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 50378.52


# not clust
shapiro.test(all.params$ash.BA.forest)
shapiro.test(all.params$ash.BA.marsh)
# Both normal
var.test(all.params$ash.BA.forest,all.params$ash.BA.marsh)
# Variances equal
t.test(all.params$ash.BA.forest,all.params$ash.BA.marsh, var.equal = TRUE)
# p-value = 0.1913 - no difference 

mean(all.params$ash.BA.forest)
std.error(all.params$ash.BA.forest)
mean(all.params$ash.BA.marsh, na.rm = TRUE)
std.error(all.params$ash.BA.marsh, na.rm = TRUE)
#### Stand BA ####
clusWilcox.test(all.params$stand.BA, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = 1.3188, p-value = 0.1872
mean(aggregate(stand.BA ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 382904.9
std.error(aggregate(stand.BA ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 52243.64
mean(aggregate(stand.BA ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 305664.2
std.error(aggregate(stand.BA ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 48969.39



# not clust
shapiro.test(all.params$stand.BA.forest)
shapiro.test(all.params$stand.BA.marsh)
# both normal
var.test(all.params$stand.BA.forest,all.params$stand.BA.marsh)
# equal variances
t.test(all.params$stand.BA.forest,all.params$stand.BA.marsh, var.equal = TRUE)
# p-value = 0.2824 - no difference 

mean(all.params$stand.BA.forest)
std.error(all.params$stand.BA.forest)
mean(all.params$stand.BA.marsh, na.rm = TRUE)
std.error(all.params$stand.BA.marsh, na.rm = TRUE)

#### Relative BA of Ash ####
clusWilcox.test(all.params$relative.BA, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -2.2958, p-value = 0.02169
mean(aggregate(relative.BA ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.5732761
std.error(aggregate(relative.BA ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.08643668
mean(aggregate(relative.BA ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.8456794
std.error(aggregate(relative.BA ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 0.03993284

#### DBH ####
clusWilcox.test(all.ash$dbh..cm., cluster = all.ash$Stand, group = as.numeric(all.ash$habitat), method = "ds")
# Z = 2.0963, p-value = 0.03606
mean(aggregate(dbh..cm. ~ Stand, all.ash[all.ash$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 26.95215
std.error(aggregate(dbh..cm. ~ Stand, all.ash[all.ash$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.2957726
mean(aggregate(dbh..cm. ~ Stand, all.ash[all.ash$habitat == "TS",], FUN = mean)[,2]) #TS mean - 19.40055
std.error(aggregate(dbh..cm. ~ Stand, all.ash[all.ash$habitat == "TS",], FUN = mean)[,2]) # SE - 1.46273


#### relative BA and density of sick and dead ash ####
class.BA <- aggregate(BA ~ simple + Plot + Stand + habitat, data = all.trees, FUN = "sum")
for (i in 1:length(class.BA$BA)) {
  class.BA$relative[i] <- class.BA$BA[i]/sum(class.BA$BA[class.BA$Plot == class.BA$Plot[i]])
}
class.BA.exp <- merge(expand.grid(simple = c("Dead","Infested","Healthy","Non-ash"), Plot = levels(as.factor(class.BA$Plot))), class.BA, by=c('Plot','simple'), all = T)
class.BA.exp.2 <- class.BA.exp %>% fill(Stand, habitat, .direction = "downup")
class.BA.exp.2[is.na(class.BA.exp.2)] <- 0

# relative BA healthy 
clusWilcox.test(class.BA.exp.2$relative[class.BA.exp.2$simple == "Healthy"], cluster = class.BA.exp.2$Stand[class.BA.exp.2$simple == "Healthy"], group = as.numeric(class.BA.exp.2$habitat[class.BA.exp.2$simple == "Healthy"]), method = "ds")
# Z = -0.59776, p-value = 0.55
mean(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Healthy" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.1118103
std.error(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Healthy" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.07680469
mean(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Healthy" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.08086032
std.error(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Healthy" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) # SE - 0.03481725

# relative BA Sick/dead
for (i in 1:(length(class.BA.exp.2$relative)/4)) {
  class.BA.exp.2$sick.dead[((4*i)-3)] <- class.BA.exp.2$relative[class.BA.exp.2$Plot == class.BA.exp.2$Plot[((4*i)-3)] & class.BA.exp.2$simple == "Infested"] + class.BA.exp.2$relative[class.BA.exp.2$Plot == class.BA.exp.2$Plot[((4*i)-3)] & class.BA.exp.2$simple == "Dead"]
}
clusWilcox.test(class.BA.exp.2$sick.dead[class.BA.exp.2$simple == "Dead"], cluster = class.BA.exp.2$Stand[class.BA.exp.2$simple == "Dead"], group = as.numeric(class.BA.exp.2$habitat[class.BA.exp.2$simple == "Dead"]), method = "ds")
# Z = -1.8147, p-value = 0.06956
mean(aggregate(sick.dead ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.5654052
std.error(aggregate(sick.dead ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.03474932
mean(aggregate(sick.dead ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.7306246
std.error(aggregate(sick.dead ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) # SE - 0.05878097

# sick
mean(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Infested" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.2396373
std.error(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Infested" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.09045691
mean(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Infested" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.3766763
std.error(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Infested" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) # SE - 0.09342973

# dead
mean(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.2911889
std.error(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.1351639
mean(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.3765655
std.error(aggregate(relative ~ Stand, class.BA.exp.2[class.BA.exp.2$simple == "Dead" & class.BA.exp.2$habitat == "TS",], FUN = mean)[,2]) # SE - 0.07994917



# relative Density healthy
class.dense <- all.trees %>% count(habitat,Stand,Plot,simple)
class.dense <- merge(expand.grid(simple = c("Dead","Infested","Healthy","Non-ash"), Plot = levels(as.factor(class.dense$Plot))), class.dense, by=c('Plot','simple'), all = T)
class.dense <- class.dense %>% fill(Stand, habitat, .direction = "downup")
class.dense[is.na(class.dense)] <- 0
for (i in 1:length(class.dense$n)) {
  class.dense$relative[i] <- class.dense$n[i]/sum(class.dense$n[class.dense$Plot == class.dense$Plot[i]]) 
}
for (i in 1:(length(class.dense$relative)/4)) {
  class.dense$sick.dead[((4*i)-3)] <- class.dense$relative[class.dense$Plot == class.dense$Plot[((4*i)-3)] & class.dense$simple == "Infested"] + class.dense$relative[class.dense$Plot == class.dense$Plot[((4*i)-3)] & class.dense$simple == "Dead"]
}

clusWilcox.test(class.dense$relative[class.dense$simple == "Healthy"], cluster = class.dense$Stand[class.dense$simple == "Healthy"], group = as.numeric(class.dense$habitat[class.dense$simple == "Healthy"]), method = "ds")
# Z = -1.0797, p-value = 0.2803
mean(aggregate(relative ~ Stand, class.dense[class.dense$simple == "Healthy" & class.dense$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.1109034
std.error(aggregate(relative ~ Stand, class.dense[class.dense$simple == "Healthy" & class.dense$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.04571666
mean(aggregate(relative ~ Stand, class.dense[class.dense$simple == "Healthy" & class.dense$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.1639425
std.error(aggregate(relative ~ Stand, class.dense[class.dense$simple == "Healthy" & class.dense$habitat == "TS",], FUN = mean)[,2]) # SE - 0.04983024

#relative Density sick/dead

clusWilcox.test(class.dense$sick.dead[class.dense$simple == "Dead"], cluster = class.dense$Stand[class.dense$simple == "Dead"], group = as.numeric(class.dense$habitat[class.dense$simple == "Dead"]), method = "ds")
# Z = 0.12016, p-value = 0.9044
mean(aggregate(sick.dead ~ Stand, class.dense[class.dense$simple == "Dead" & class.dense$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.5451896
std.error(aggregate(sick.dead ~ Stand, class.dense[class.dense$simple == "Dead" & class.dense$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.1209173
mean(aggregate(sick.dead ~ Stand, class.dense[class.dense$simple == "Dead" & class.dense$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.5404586
std.error(aggregate(sick.dead ~ Stand, class.dense[class.dense$simple == "Dead" & class.dense$habitat == "TS",], FUN = mean)[,2]) # SE - 0.0.0549509



#### average dbh of healthy vs sick/dead in TS ####
all.ash$sickornot <- as.factor(ifelse(all.ash$simple == "Infested" | all.ash$simple == "Dead", "Sick/Dead", "Healthy"))

mean(aggregate(dbh..cm. ~ Stand, all.ash[all.ash$sickornot == "Sick/Dead" & all.ash$habitat == "TS",], FUN = mean)[,2]) # Sick/Dead - 21.92953
mode(aggregate(dbh..cm. ~ Stand, all.ash[all.ash$sickornot == "Sick/Dead" & all.ash$habitat == "TS",], FUN = mean)[,2])

model <- glmer(sickornot ~ dbh..cm.*habitat + (1|Stand/Plot), data = all.ash, family = binomial)
summary(model)
effect_plot(model,pred = habitat, plot.points = T, interval = T)

par(mfcol = c(2,2))
plot(model)
par(mfcol = c(1,1))

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: binomial  ( logit )
#Formula: sickornot ~ dbh..cm. * habitat + (1 | Stand/Plot)
#Data: all.ash
#
#AIC      BIC   logLik deviance df.resid 
#238.6    259.7   -113.3    226.6      243 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-8.4189  0.1709  0.3495  0.4885  1.2711 
#
#Random effects:
#  Groups     Name        Variance  Std.Dev. 
#Plot:Stand (Intercept) 3.681e-01 6.067e-01
#Stand      (Intercept) 8.072e-09 8.984e-05
#Number of obs: 249, groups:  Plot:Stand, 21; Stand, 8
#
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)         1.92285    0.81214   2.368 0.017902 *  
#  dbh..cm.           -0.01158    0.02485  -0.466 0.641264    
#habitatTS          -2.53542    0.92881  -2.730 0.006338 ** 
#  dbh..cm.:habitatTS  0.12877    0.03573   3.604 0.000314 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Correlation of Fixed Effects:
#  (Intr) dbh... hbttTS
#dbh..cm.    -0.866              
#habitatTS   -0.873  0.757       
#dbh..cm.:TS  0.613 -0.695 -0.808
#convergence code: 0
#boundary (singular) fit: see ?isSingular


all.ash$fit <- fitted(model, all.ash, type="response")

ggplot()+geom_line(data = all.ash[all.ash$habitat == "TS",], aes(x = dbh..cm., y = fit), col = "green")+
  geom_line(data = all.ash[all.ash$habitat == "BHF",], aes(x = dbh..cm., y = fit), col = "red")


#### slowed ash mortality ####
all.ash <- merge(all.ash, all.params[,c(2,4)], by.x ="Plot", by.y = "plot", all = T)

infected.ash <- class.dense[class.dense$simple == "Infested" | class.dense$simple == "Dead",-c(6:7)]
for (i in 1:length(infected.ash$n)) {
  infected.ash$relative[i] <- infected.ash$n[i]/sum(infected.ash$n[infected.ash$Plot == infected.ash$Plot[i]]) 
}
infected.ash <- infected.ash[infected.ash$simple == "Infested",]
infected.ash$arcsine <- asin(sqrt(infected.ash$relative))
infected.ash$density <- all.params$ash.density

density.mortality <- lmer(arcsine ~ density*habitat + (1|Stand), data = infected.ash)
summary(density.mortality)

par(mfcol = c(2,2))
plot(density.mortality)
par(mfcol = c(1,1))

#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
#Formula: arcsine ~ density * habitat + (1 | Stand)
#   Data: infected.ash
#
#REML criterion at convergence: 55.9
#
#Scaled residuals: 
#     Min       1Q   Median       3Q      Max 
#-1.55451 -0.66397  0.02979  0.56939  1.46486 
#
#Random effects:
# Groups   Name        Variance Std.Dev.
# Stand    (Intercept) 0.07727  0.2780  
# Residual             0.18334  0.4282  
#Number of obs: 21, groups:  Stand, 8
#
#Fixed effects:
#                    Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept)        0.7312511  0.3193733 13.4341007   2.290   0.0388 *
#density            0.0002693  0.0006263 15.6676624   0.430   0.6730  
#habitatTS          0.2907961  0.4197650 14.3573840   0.693   0.4995  
#density:habitatTS -0.0006360  0.0007008 15.5041062  -0.908   0.3780  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Correlation of Fixed Effects:
#            (Intr) densty hbttTS
#density     -0.752              
#habitatTS   -0.761  0.572       
#dnsty:hbtTS  0.672 -0.894 -0.728


#### percentage of habitat for each species ####
spp.counts <- count(all.trees, spp, Plot, Stand, habitat) %>% 
  group_by(Plot,) %>% 
  complete(spp, fill = list(n = 0)) %>% 
  fill(Stand, habitat, .direction = "downup") %>%
  mutate( percentage = n/sum(n))


# ash
mean(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "green" & spp.counts$habitat == "TS",], FUN = mean)[,2]) 
std.error(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "green" & spp.counts$habitat == "TS",], FUN = mean)[,2]) 


mean(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "green" & spp.counts$habitat == "BHF",], FUN = mean)[,2]) 
std.error(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "green" & spp.counts$habitat == "BHF",], FUN = mean)[,2]) 

# red maple
mean(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "red maple" & spp.counts$habitat == "TS",], FUN = mean)[,2]) 
std.error(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "red maple" & spp.counts$habitat == "TS",], FUN = mean)[,2]) 

mean(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "red maple" & spp.counts$habitat == "BHF",], FUN = mean)[,2]) 
std.error(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "red maple" & spp.counts$habitat == "BHF",], FUN = mean)[,2]) 

# sweetgum
mean(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "sweetgum" & spp.counts$habitat == "TS",], FUN = mean)[,2]) 
std.error(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "sweetgum" & spp.counts$habitat == "TS",], FUN = mean)[,2]) 

# persimmon
mean(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "persimmon" & spp.counts$habitat == "BHF",], FUN = mean)[,2]) 
std.error(aggregate(percentage ~ Stand, spp.counts[spp.counts$spp == "persimmon" & spp.counts$habitat == "BHF",], FUN = mean)[,2]) 


################ Understory ################
#### Total Groundcover ####
clusWilcox.test(all.params$total.groundcover, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = 0.82627, p-value = 0.4087
mean(aggregate(total.groundcover ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 70.10648
std.error(aggregate(total.groundcover ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 4.019313
mean(aggregate(total.groundcover ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 57.48
std.error(aggregate(total.groundcover ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 11.7111


# no clust
shapiro.test(all.params$total.groundcover.forest)
# normal
shapiro.test(all.params$total.groundcover.marsh)
# normal
var.test(all.params$total.groundcover.forest,all.params$total.groundcover.marsh)
# Variances Equal
t.test(all.params$total.groundcover.forest,all.params$total.groundcover.marsh)
# p-value = 0.9776 no difference 

#### Invasive Groundcover ####
clusWilcox.test(all.params$invas.groudcover, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = 2.0131, p-value = 0.0441
mean(aggregate(invas.groudcover ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 84.84722
std.error(aggregate(invas.groudcover ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 17.48353
mean(aggregate(invas.groudcover ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 23.68333
std.error(aggregate(invas.groudcover ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 13.68897



# no clust
shapiro.test(all.params$invas.groudcover.forest)
shapiro.test(all.params$invas.groundcover.marsh)
# not normal
var.test(all.params$invas.groudcover.forest,all.params$invas.groundcover.marsh)
# Variances Equal
wilcox.test(all.params$invas.groudcover.forest,all.params$invas.groundcover.marsh)
boxplot(all.params$invas.groudcover.forest,all.params$invas.groundcover.marsh)
# p-value = 0.004193 Forest Greater 
#### ANOSIM for understory comm ####
# run nmds for herb datas first
anosim(total.abun, grouping = site.scores$community, permutations = 999, distance = "bray")
## R 0.8125 p = 0.001 
################ Regen ################
#### Ash Seeding Density ####
clusWilcox.test(all.params$ash.seed.dens, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -2.1964, p-value = 0.02807
mean(aggregate(ash.seed.dens ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 52.09877
std.error(aggregate(ash.seed.dens ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 23.96589
mean(aggregate(ash.seed.dens ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 700.7407
std.error(aggregate(ash.seed.dens ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 284.9648


# no clust
shapiro.test(all.params$ash.seed.dens.forest)
# not normal
shapiro.test(all.params$ash.seed.dens.marsh)
# not normal
var.test(all.params$ash.seed.dens.forest,all.params$ash.seed.dens.marsh)
# not normal
wilcox.test(all.params$ash.seed.dens.forest,all.params$ash.seed.dens.marsh)
boxplot(all.params$ash.seed.dens.forest,all.params$ash.seed.dens.marsh)
# p-value = 0.0008339 marsh greater
#### Ash Seedlings Relative Density ####
clusWilcox.test(all.params$ash.seed.relative.dens, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -2.2461, p-value = 0.0247
mean(aggregate(ash.seed.relative.dens ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.255262
std.error(aggregate(ash.seed.relative.dens ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.1154856
mean(aggregate(ash.seed.relative.dens ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.8260192
std.error(aggregate(ash.seed.relative.dens ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 0.0908896


# no cust
shapiro.test(all.params$porportion.ash.seed.forest)
#not normal
shapiro.test(all.params$porportion.ash.seed.marsh)
# not normal
var.test(all.params$porportion.ash.seed.forest,all.params$porportion.ash.seed.marsh)
# normal
wilcox.test(all.params$porportion.ash.seed.forest,all.params$porportion.ash.seed.marsh)
boxplot(all.params$porportion.ash.seed.forest,all.params$porportion.ash.seed.marsh)
# p-value = 0.002199 higher proportion of ash seedlings in marsh
#### Basal Shoot Density ####
clusWilcox.test(all.params$basal.shoot.dens, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -2.0588, p-value = 0.03951
mean(aggregate(basal.shoot.dens ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 322.2222
std.error(aggregate(basal.shoot.dens ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 54.22269
mean(aggregate(basal.shoot.dens ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 1248.889
std.error(aggregate(basal.shoot.dens ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 489.1335


# no clust
shapiro.test(all.params$basal.shoot.dens.forest)
# not normal
shapiro.test(all.params$basal.shoot.dens.marsh)
# not normal
var.test(all.params$basal.shoot.dens.forest,all.params$basal.shoot.dens.marsh)
# not equal
wilcox.test(all.params$basal.shoot.dens.forest,all.params$basal.shoot.dens.marsh)
boxplot(all.params$basal.shoot.dens.forest,all.params$basal.shoot.dens.marsh)
# p-value = 0.01186 Marsh Greater

#### Proportion of Ash With Shoots ####
clusWilcox.test(all.params$proportion.trees.w.shoots, cluster = all.params$stand, group = as.numeric(all.params$habitat),  method = "ds")
# Z = -1.7091, p-value = 0.08743
mean(aggregate(proportion.trees.w.shoots ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) #BHF mean - 0.1931397
std.error(aggregate(proportion.trees.w.shoots ~ stand, all.params[all.params$habitat == "BHF",], FUN = mean)[,2]) # SE - 0.03932024
mean(aggregate(proportion.trees.w.shoots ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) #TS mean - 0.3870632
std.error(aggregate(proportion.trees.w.shoots ~ stand, all.params[all.params$habitat == "TS",], FUN = mean)[,2]) # SE - 0.1032269


#### ANOSIM for regen ####
BHFanosim <- anosim(BHF, grouping = BHFcommunity, permutations = 999, distance = "bray")
BHFanosim
## R = 0.3235, p = 0.001 ##
TSanosim <- anosim(TS, grouping = TScommunity, permutations = 999, distance = "bray")
TSanosim
## R= 0.1874, p = 0.015 ##
#### Basal Shoots per Ash - canopy health ####

# per canopy class
basal.shoots <- aggregate(Shoots ~ simple + Plot + Stand + habitat, data = all.ash, FUN = mean)
basal.shoots.cast <- dcast(basal.shoots, ... ~ simple)

all.ash$habitat <- relevel(all.ash$habitat, ref = "BHF")

basal.shoot.model <- glmer(Shoots ~ simple*habitat + (1|Stand/Plot), data = all.ash, family = poisson)
summary(basal.shoot.model)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: poisson  ( log )
#Formula: Shoots ~ simple * habitat + (1 | Stand/Plot)
#Data: all.ash
#
#AIC      BIC   logLik deviance df.resid 
#969.9    998.1   -477.0    953.9      241 
#
#Scaled residuals: 
#  Min     1Q Median     3Q    Max 
#-2.608 -1.037 -0.507  0.000  9.194 
#
#Random effects:
#  Groups     Name        Variance  Std.Dev. 
#Plot:Stand (Intercept) 1.480e+00 1.217e+00
#Stand      (Intercept) 7.126e-13 8.442e-07
#Number of obs: 249, groups:  Plot:Stand, 21; Stand, 8
#
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)  
#(Intercept)               -0.5223     0.5136  -1.017   0.3092  
#simpleHealthy            -20.5499    61.6377  -0.333   0.7388  
#simpleInfested            -0.5536     0.2968  -1.865   0.0622 .
#habitatTS                  1.1345     0.6338   1.790   0.0735 .
#simpleHealthy:habitatTS   18.2812    61.6377   0.297   0.7668  
#simpleInfested:habitatTS  -0.2833     0.3206  -0.884   0.3769  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Correlation of Fixed Effects:
#  (Intr) smplHl smplIn hbttTS smH:TS
#simpleHlthy  0.000                            
#simplInfstd -0.226  0.000                     
#habitatTS   -0.795  0.000  0.181              
#smplHlth:TS  0.000 -1.000  0.000  0.000       
#smplInfs:TS  0.215  0.000 -0.927 -0.199  0.000
#convergence code: 0
#boundary (singular) fit: see ?isSingular

clusWilcox.test(basal.shoots.cast$Dead, cluster = basal.shoots.cast$Stand, group = as.numeric(basal.shoots.cast$habitat),  method = "ds")
# Z = -1.6145, p-value = 0.1064
clusWilcox.test(basal.shoots.cast$Infested, cluster = basal.shoots.cast$Stand, group = as.numeric(basal.shoots.cast$habitat),  method = "ds")
# Z = -1.4975, p-value = 0.1343
clusWilcox.test(basal.shoots.cast$Healthy, cluster = basal.shoots.cast$Stand, group = as.numeric(basal.shoots.cast$habitat),  method = "ds")
# Z = -1.4061, p-value = 0.1597

# per tree
clusWilcox.test(all.ash$Shoots, cluster = all.ash$Stand, group = as.numeric(all.ash$habitat), method = "ds")
# Z = -0.93261, p-value = 0.351

# for trees with shoots
trees.w.shoots <- all.ash[all.ash$Shoots > 0,]
clusWilcox.test(trees.w.shoots$Shoots, cluster = trees.w.shoots$Stand, group = as.numeric(trees.w.shoots$habitat), method = "ds")
# Z = 1.4818, p-value = 0.1384

################ Figures ################
#### canopyloss x BA & dense plot ####
library(ggplot2)
library(dplyr)
library(ggpubr)

#simple index - healthy/infested/dead/non-ash
# stand level - BA
simple.BA.2 <-aggregate(BA ~ simple + Stand + habitat, data=all.trees, FUN = sum)
for (i in 1:length(simple.BA.2$BA)) {
  simple.BA.2$relative[i] <- simple.BA.2$BA[i]/sum(simple.BA.2$BA[simple.BA.2$Stand == simple.BA.2$Stand[i]])
}
simple.BA.2$simple <- as.factor(simple.BA.2$simple)
simple.BA.2$simple <- relevel(simple.BA.2$simple, ref = "Non-ash")
simple.BA.2$simple <- relevel(simple.BA.2$simple, ref = "Healthy")
simple.BA.2$simple <- relevel(simple.BA.2$simple, ref = "Infested")
simple.BA.2$simple <- relevel(simple.BA.2$simple, ref = "Dead")

simple.ba.plot <- ggplot(data = simple.BA.2, aes(x=Stand,y=relative,fill=simple))+
  geom_bar(position = "stack", stat = "identity")+facet_grid(~habitat, scales = "free", space = "free_x")+
  theme_bw()+scale_fill_manual("Canopy Class",values = c("grey90","grey80","grey30","black"))+
  ylab("Relative Basal Area")+ggtitle("") + theme(text = element_text(size = 20))
simple.ba.plot

# stand level - density
simple.dens <- all.trees %>% count(Stand,simple,habitat)
for (i in 1:length(simple.dens$n)) {
  simple.dens$relative[i] <- simple.dens$n[i]/sum(simple.dens$n[simple.dens$Stand == simple.dens$Stand[i]])
}
simple.dens$simple <- as.factor(simple.dens$simple)
simple.dens$simple <- relevel(simple.dens$simple, ref = "Non-ash")
simple.dens$simple <- relevel(simple.dens$simple, ref = "Healthy")
simple.dens$simple <- relevel(simple.dens$simple, ref = "Infested")
simple.dens$simple <- relevel(simple.dens$simple, ref = "Dead")

simple.dens.plot <- ggplot(data = simple.dens, aes(x=Stand,y=relative,fill=simple))+
  geom_bar(position = "stack", stat = "identity")+facet_grid(~habitat, scales = "free", space = "free_x")+
  theme_bw()+scale_fill_manual("Canopy Class",values = c("grey90","grey80","grey30","black"))+
  ylab("Relative Density") + theme(text = element_text(size = 20))
simple.dens.plot

#both 
both.plot <- ggarrange(simple.ba.plot,simple.dens.plot, common.legend = T, legend = "right", align = "hv")

both.plot


# not using these
# combo BA and density
combo <- merge(simple.BA.2,simple.dens, by=c('simple','Stand'))
combo <- combo[,-c(3,5)]
colnames(combo) <- c("simple","Stand","rel.BA","rel.dens")
combo.melt <- melt(data = combo, id.vars = c("simple","Stand"), measure.vars = c("rel.BA","rel.dens"))
ggplot(data = combo.melt, aes(x=variable, y=value, fill=simple))+geom_bar(position = "stack", stat = "identity")+facet_grid(~ Stand)

# Addative BA with Smith Index
add.BA <- aggregate(BA ~ smith.index + Stand, data = all.trees, FUN = sum)
add.BA.ash <- aggregate(BA ~ smith.index + Stand, data = all.ash, FUN = sum)
add.Dens <- all.ash %>% count(Stand,smith.index)

# all trees BA
ggplot(data = add.BA, aes(y=BA,x=Stand,fill=smith.index))+geom_bar(position = "stack", stat = "identity")
# just ash BA
ggplot(data = add.BA.ash, aes(y=BA,x=Stand,fill=smith.index))+geom_bar(position = "stack", stat = "identity")
# all trees density
ggplot(data = add.Dens, aes(y=n,x=Stand,fill=smith.index))+geom_bar(position = "stack", stat = "identity")

#Relative BA
ave.can <- aggregate(binned.canopy ~ Stand, all.ash, FUN = mean)
ave.can$proportion <- (ave.can$binned.canopy/5)
BA <- all.params[,c(1:3,9)]
stand.BA <- aggregate(relative.BA ~ stand, BA, FUN = mean)
stand.BA$otherSPP <- 1-stand.BA$relative.BA
stand.BA$ash.remain <- stand.BA$relative.BA*ave.can$proportion
stand.BA$ash.lost <- stand.BA$relative.BA-stand.BA$ash.remain

barplot(t(as.matrix(stand.BA[,c(3:5)])))
ggplot(data = stand.BA)



#### Histograms of Canopy Class vs size ####
healthy.ts <- subset(all.ash, simple == "Healthy" & habitat == "TS")
hist(healthy.ts$dbh..cm.)
healthy.bhf <- subset(all.ash, simple != "Healthy" & habitat == "BHF")
hist(healthy.bhf$dbh..cm.)

ggplot()+geom_histogram(data = healthy.ts, aes(x = dbh..cm.), binwidth = 5, center = 2.5, fill = "black", color= "white" ) +
  geom_histogram(data = healthy.bhf, aes(x = dbh..cm.), binwidth = 5, center = 2.5, fill = "grey50", color = "white", alpha = 0.6 ) +
  theme_bw() 

all.ash$sickornot <- ifelse(all.ash$simple == "Healthy", "Healthy", "Infested/Dead")

ggplot(data = all.ash, aes(x=dbh..cm.)) + geom_histogram(binwidth = 5, center = 2.5, color = "white") + 
  facet_grid(habitat ~ sickornot) + theme_bw() + ylab("Count") + xlab("DBH") + theme(text = element_text(size = 20))
  

#### Regen Plot ####

regenplot <- all.params[,c(1:3,13,15)]
regenplot <- melt(all.params, id.vars = c('habitat','plot','stand'), measure.vars = c('ash.seed.dens','basal.shoot.dens'))

regen <- aggregate(value ~ stand + variable + habitat, data = regenplot, FUN = mean)


ggplot(regen, aes(y = value, x = variable, fill=habitat)) + 
    geom_boxplot(width=0.8) + scale_fill_grey(start = 1, end = 0.4) + xlab("")+
    theme_bw(base_size = 20) + ylab(bquote(~Stems~ha^-1)) + scale_x_discrete(breaks=NULL)+ 
    theme(legend.text = element_text(size = 15),legend.position = c(0.15,0.82), legend.title = element_blank(),legend.spacing.y = unit(0, "mm"), legend.background = element_rect(color = "black", fill = "white", size = 0.35))+
    annotate(geom = "text", x = c(1,2), y = -300 , label = c("Seedlings", "Basal Shoots"), size = 6)+
    coord_cartesian(ylim = c(0,3000), expand = T, clip = "off") + theme(panel.border = element_rect(size = 0.35))
  
  
  
  theme(legend.text = element_text(size = 10),legend.position = c(0.19,0.90), legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),panel.border = element_rect(colour = "black", fill=NA),aspect.ratio = 1, axis.text = element_text(colour = 1, size = 8), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))
#ggplot(regenplot, aes(y = value, x = variable, fill=variable)) + 
 # geom_boxplot(width=0.8) + theme_bw(base_size = 25) + theme(axis.text.x = element_text(angle = 340, hjust = 0, size = 18)) + scale_fill_grey(start = 1, end = 0.4)+theme(legend.position = "none") + labs(x = element_blank(), y = "Stems per Hectare")
  



#### nmds for clustering ####
rm(list = ls())
library(vegan)
library(ggrepel)
library(ggpubr)
library(ggplot2)

setwd("~/Documents/Old Documents")
ash_plot_speciescomp <- read.csv("ash_plot_speciescomp.csv")

total.abun <- ash_plot_speciescomp[-c(7,14),c(1,5,44:80)]
colnames(total.abun) <- colnames(ash_plot_speciescomp[,c(1,5,7:43)])
total.abun$plot.type <- ifelse(total.abun$plot.type == "forest", "BHF", "TS")
plots <- total.abun[,1]
type <- total.abun[,2]
total.abun <- total.abun[,-c(1:2)]

total.abun.mds <- metaMDS(total.abun)
stressplot(total.abun.mds)

ordiplot(total.abun.mds)
orditorp(total.abun.mds, display = "species", col = "red")
orditorp(total.abun.mds, display = "sites", col = "black", cex = 1)

ordiplot(total.abun.mds)
ordihull(total.abun.mds, groups = type, draw = "polygon", label = T)

# ggplot
spec.scores <- as.data.frame(scores(total.abun.mds,"species"))
spec.scores$species <- rownames(spec.scores)
site.scores <- as.data.frame(scores(total.abun.mds,"sites"))
site.scores$community <- type
grp.a <- site.scores[site.scores$community == "BHF", ][chull(site.scores[site.scores$community == "BHF", c("NMDS1", "NMDS2")]), ]  
grp.b <- site.scores[site.scores$community == "TS", ][chull(site.scores[site.scores$community == "TS", c("NMDS1", "NMDS2")]), ]  
hulls <- rbind(grp.a,grp.b)


ggplot()+geom_polygon(data = hulls,aes(x=NMDS1,y=NMDS2,fill=community,group=community), alpha = 0.8)+
  scale_fill_manual(values = c("grey80","grey5"))+
  geom_label_repel(data=spec.scores, aes(x=NMDS1,y=NMDS2,label=species), label.padding = 0.2, point.padding = 0, label.size = 0.17)+
  theme_bw() + theme(plot.title = element_text(size = 20)) + labs(x="",y="") + 
  ggtitle("Ground Cover Species") + scale_colour_grey(start = 0, end = 0.6, name = "") + 
  theme(legend.text = element_text(size = 13),legend.position = c(0.11,0.92), legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),panel.border = element_rect(colour = "black", fill=NA),aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))



#### Making ordinations for seedlings and overstory species ####
rm(list = ls())
library(vegan)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(tidyr)
setwd("~/Documents")
TreesAndSeeds <- read.csv(file = "JugBay_TreeNMDS2.csv",header = TRUE)
TreesAndSeeds <- TreesAndSeeds[c(1:44),-19]
TreesAndSeeds[is.na(TreesAndSeeds)]=0
row.names(TreesAndSeeds) <- TreesAndSeeds[,2]
BHF <- TreesAndSeeds[c(1:22),-c(1:2)]
TS <- TreesAndSeeds[c(23:44),-c(1:2)]
## BHF Ordination ##
BHFcommunity <- c("Seedlings","Seedlings","Seedlings","Seedlings","Seedlings","Seedlings","Seedlings","Seedlings","Seedlings","Seedlings","Trees","Trees","Trees","Trees","Trees","Trees","Trees","Trees","Trees","Trees","Trees","Trees")
#BHFcolors <- c(rep("red",10),rep("blue",12))

BHFnmds <-metaMDS(BHF, distance = "bray", k=2)
# stressplot(BHFnmds)
# ordiplot(BHFnmds,type = "n")
# orditorp(BHFnmds,display="species",col="black",air=1)
# BHFelipse1 <- ordiellipse(BHFnmds, groups = BHFcommunity, draw = "polygon", kind = "sd", col = "red")
# BHFelipse2 <- ordiellipse(BHFnmds, groups = BHFcommunity, draw = "polygon", label = TRUE, kind = "ehull")

## doing it in ggplot ##
BHFscores <- as.data.frame(scores(BHFnmds,"species"))
BHFscores$species <- rownames(BHFscores)
BHFsite.scores <- as.data.frame(scores(BHFnmds,"sites"))
BHFsite.scores$community <- BHFcommunity
grp.a <- BHFsite.scores[BHFsite.scores$community == "Seedlings", ][chull(BHFsite.scores[BHFsite.scores$community == "Seedlings", c("NMDS1", "NMDS2")]), ]  
grp.b <- BHFsite.scores[BHFsite.scores$community == "Trees", ][chull(BHFsite.scores[BHFsite.scores$community == "Trees", c("NMDS1", "NMDS2")]), ]  
BHFhulls <- rbind(grp.a,grp.b)


BHFplot <- ggplot()+geom_polygon(data = BHFhulls,aes(x=NMDS1,y=NMDS2,fill=community,group=community), alpha = 0.8)+
            scale_fill_manual(values = c("grey5","grey80"))+
            geom_label_repel(data=BHFscores, aes(x=NMDS1,y=NMDS2,label=species), label.padding = 0.2, point.padding = 0, label.size = 0.25)+
            theme_bw() + 
            theme(plot.title = element_text(size = 20)) + 
            labs(x="",y="") + ggtitle("Bottomland Hardwood Forest") + scale_colour_grey(start = 0, end = 0.6, name = "") + 
            theme(legend.text = element_text(size = 11),legend.position = c(0.19,0.90), legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),panel.border = element_rect(colour = "black", fill=NA),aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))

## TS Ordination ##
TScommunity <- c(rep("Seedlings",11),rep("Trees",11))
#TScolors <- c(rep("red",10),rep("blue",12))

TSnmds <-metaMDS(TS, distance = "bray", k=2)
#stressplot(TSnmds)
#ordiplot(TSnmds,type = "n")
#orditorp(TSnmds,display="species",col="red",air=0.01)
#ordiellipse(TSnmds, groups = TScommunity, draw = "polygon", label = TRUE)

## doing it in ggplot ##
TSscores <- as.data.frame(scores(TSnmds,"species"))
TSscores$species <- rownames(TSscores)
TSsite.scores <- as.data.frame(scores(TSnmds,"sites"))
TSsite.scores$community <- TScommunity
grp.c <- TSsite.scores[TSsite.scores$community == "Seedlings", ][chull(TSsite.scores[TSsite.scores$community == "Seedlings", c("NMDS1", "NMDS2")]), ]  
grp.d <- TSsite.scores[TSsite.scores$community == "Trees", ][chull(TSsite.scores[TSsite.scores$community == "Trees", c("NMDS1", "NMDS2")]), ]  
TShulls <- rbind(grp.c,grp.d)


TSplot <- ggplot()+geom_polygon(data = TShulls,aes(x=NMDS1,y=NMDS2,fill=community,group=community), alpha = 0.8) +
            scale_fill_manual(values = c("grey5","grey80"))+
            geom_label_repel(data=TSscores, aes(x=NMDS1,y=NMDS2,label=species), label.padding = 0.2, point.padding = 0, label.size = 0.25)+
            theme_bw()+
            theme(plot.title = element_text(size = 20)) + 
            labs(x="",y="") + ggtitle("Tidal Swamp") + 
            theme(legend.text = element_text(size = 11),legend.position = c(0.19,0.90), legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),panel.border = element_rect(colour = "black", fill=NA),aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))
TSplot
## plotting them together ##
ggarrange(TSplot,BHFplot, align = "hv")

#### common seedlings density stats ####
seedlings <- separate(TreesAndSeeds, Plot, into = c('plot','type'), sep = '\ ')
seedlings <- separate(seedlings, Community, into = c('habitat','type1'), sep = '\ ')
seedlings <- seedlings[seedlings$type1 == "Seedlings", -c(2,4)]
seedlings <- merge(all.params[,1:3], seedlings, by = c('habitat','plot'))

# Blackhaw - TS
mean(aggregate(VIPR ~ stand, seedlings[seedlings$habitat == "TS",], FUN = mean)[,2])/0.0225
std.error(aggregate(VIPR ~ stand, seedlings[seedlings$habitat == "TS",], FUN = mean)[,2])/0.0225

# Swamp Tupelow - TS
mean(aggregate(NYBI ~ stand, seedlings[seedlings$habitat == "TS",], FUN = mean)[,2])/0.0225
std.error(aggregate(NYBI ~ stand, seedlings[seedlings$habitat == "TS",], FUN = mean)[,2])/0.0225

# american beech - BHF
mean(aggregate(FAGR ~ stand, seedlings[seedlings$habitat == "BHF",], FUN = mean)[,2])/0.0225
std.error(aggregate(FAGR ~ stand, seedlings[seedlings$habitat == "BHF",], FUN = mean)[,2])/0.0225

#sweetgum - BHF
mean(aggregate(LIST ~ stand, seedlings[seedlings$habitat == "BHF",], FUN = mean)[,2])/0.0225
std.error(aggregate(LIST ~ stand, seedlings[seedlings$habitat == "BHF",], FUN = mean)[,2])/0.0225

# persimmon - BHF
mean(aggregate(DIVI ~ stand, seedlings[seedlings$habitat == "BHF",], FUN = mean)[,2])/0.0225
std.error(aggregate(DIVI ~ stand, seedlings[seedlings$habitat == "BHF",], FUN = mean)[,2])/0.0225


#######################################
#####
#####
