#Microcolony analysis

#loading libraries
library(ggplot2)
library(ggsurvfit) #for survival pltos
library(dplyr)
library(ggpubr)
library(plyr)
library(tidyverse) 
library(car)       # Companion to applied regression 
library(tseries)   # For timeseries analysis
library(lmtest)    #For hetoroskedasticity analysis
library(mosaic)    #For means and sds
library(janitor)
library(lme4)
library(DescTools) #for Duncan analysisicrocolony data 
library(survival)
library(survminer)
library(knitr)

#Pollen consumption per Microcolony
#loading file
pollen<-read.csv("Trial2_pollen.csv")
View(pollen)

#transforming to correct fromat 
pollen$Treatment<-as.factor(pollen$Treatment)
pollen$Temperature<-as.factor(pollen$Temperature)
pollen$Humidity<-as.factor(pollen$Humidity)
pollen$Microcolony<-as.factor(pollen$Microcolony)
pollen$Colony_source<-as.factor(pollen$Colony_source)
pollen$Feeding_event<-as.factor(pollen$Feeding_event)


#checking linearity of data
#linearity for number of workers
model3<-lm(Pollen_consumed_g~Temperature*Humidity,data=pollen) 
plot(model3,1) #trend line looks linear

#checking variance of data
bartlett.test(Pollen_consumed_g~interaction(Temperature,Humidity),data=pollen) #p-value = 0.7278

#getting averages of pollen consumption per microcolony
pollen2=ddply(pollen, .(Treatment,Temperature,Humidity,Microcolony,Colony_source), numcolwise(mean),na.rm=T)
View(pollen2)

#doing a linear mixed model 
pollen_lmer<-lmer(Pollen_consumed_g~Temperature*Humidity + (1|Colony_source),
                  data=pollen2)
summary(pollen_lmer) 

#REML criterion at convergence: -11.5
#Estimate Std. Error t value
#(Intercept)                                 0.91914    0.06613  13.899
#pollen2$TemperatureYes                      0.22489    0.08002   2.811
#pollen2$HumidityYes                         0.16181    0.08257   1.960
#pollen2$TemperatureYes:pollen2$HumidityYes -0.41288    0.11498  -3.591


Anova(pollen_lmer, type = 3, test.statistic="F") 
#                                  F Df  Df.res    Pr(>F)    
#(Intercept)          193.1954  1  8.7302 2.965e-07 ***
#Temperature            7.8989  1 29.0004  0.008772 ** 
#Humidity               3.8268  1 29.0938  0.060100 .  
#Temperature:Humidity  12.8703  1 29.0488  0.001208 ** 
---

  hist(resid(pollen_lmer)) #slightly right skewed
car::qqPlot(resid(pollen_lmer))

#plotting pollen consumption on boxplot 
plw1 = ggplot(pollen2,aes(x=Treatment,y=Pollen_consumed_g, fill=Treatment))+
  geom_boxplot()+
  labs(y="Pollen Consumption per Microcolony (g)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave &
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08),size=2.25) +
  scale_fill_brewer(palette = "Set2")+ng4;plw1

#lipids per microcolony (mg)
#transforming my data into appropriate formats
lipids$Microcolonyl_Id<-as.factor(lipids$Microcolonyl_Id)
lipids$Colony_source<-as.factor(lipids$Colony_source)
lipids$Treatment<-as.factor(lipids$Treatment)
lipids$Temperature<-as.factor(lipids$Temperature)
lipids$Humidity<-as.factor(lipids$Humidity)
lipids$Worker_id<-as.factor(lipids$Worker_id)



#In the data frame, I have individual lipids content, but my experimnetal unit was microcolony
#below is a new data frame with that information based on microcolony ID 

l=ddply(lipids, .(Treatment,Microcolonyl_Id,Temperature, Humidity,Colony_source), numcolwise(mean),na.rm=T)
View(l)
l$Thx_Wet_weight<-as.numeric

#---------------------------------------- will first look into abdominal lipids
#Noramlity - looks good
#Normality of three categories of weight - wet, initial dry, final dry

#wet weight 
hist(l$Abd_Wet_weight)
shapiro.test(l$Abd_Wet_weight)

#inital dry mass 
hist(l$Abd_Initial_dry_weight) 
shapiro.test(l$Abd_Initial_dry_weight)

#final dry mass
hist(l$Abd_48h_dry_weight)
shapiro.test(l$Abd_48h_dry_weight)

#---------------------------------------- Variance - looks good

#wet weigh 
bartlett.test(l$Abd_Wet_weight~l$Treatment,
              data=l)


#inital dry weight 
bartlett.test(l$Abd_Initial_dry_weight~l$Treatment,
              data=l)

#final dry weight 
bartlett.test(l$Abd_48h_dry_weight~l$Treatment,
              data=l)

#---------------------------------------- Linearity -looks good

#wet weigh 
model1<-lm(Abd_Wet_weight~Treatment,data=l) 
plot(model1,1) 

#initial dry weight 
model2<-lm(Abd_Initial_dry_weight~Treatment,data=l) 
plot(model2,1) 


#final dry weight 
model3<-lm(Abd_48h_dry_weight~Treatment, data=l)
plot(model3,1)

#--------------------------------------- will do a lmer, but first neet to remove outliers
#will remove outliers!
summary(l$Diff_Abd_dry_weight)
IQR(l$Diff_Abd_dry_weight)

#getting min and max lipid values for thresholds
lmin = 0.001685-(1.5*0.001075) 
lmax =  0.002760+(1.5*0.001075) 

# removing outlier
lo<-l[!(l$Diff_Abd_dry_weight<lmin | l$Diff_Abd_dry_weight>lmax),]
View(lo)

#checking assumptions to do lmer
modello<-lm(Diff_Abd_dry_weight~Treatment, data=lo)
plot(modello,1) #linear

lipidso_lmer<-lmer(lo$Diff_Abd_dry_weight*1000~lo$Temperature*lo$Humidity+lo$Abd_Wet_weight+
                     (1|lo$Colony_source))
summary(lipidso_lmer)

#Random effects:
#Groups           Name        Variance Std.Dev.
#lo$Colony_source (Intercept) 0.0000   0.0000  
#Residual                     0.3615   0.6012  
#Number of obs: 30, groups:  lo$Colony_source, 3

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept)                        0.3020     0.7652   0.395
#lo$TemperatureYes                  0.4425     0.3096   1.429
#lo$HumidityYes                     0.2762     0.3044   0.907
#lo$Abd_Wet_weight                 21.3885     9.1484   2.338
#lo$TemperatureYes:lo$HumidityYes  -0.8816     0.4435  -1.988

Anova(lipidso_lmer, test.statistic="F")
#Response: lo$Diff_Abd_dry_weight
#F Df Df.res  Pr(>F)  
#lo$Temperature             0.0263  1 24.329 0.87251  
#lo$Humidity                0.3750  1 23.447 0.54616  
#lo$Abd_Wet_weight          5.3712  1 23.355 0.02957 *
#lo$Temperature:lo$Humidity 3.8619  1 23.480 0.06133 .

#plotting residuals
plot(resid(lipidso_lmer)~fitted(lipidso_lmer)) 
hist(resid(lipidso_lmer))
car::qqPlot(resid(lipidso_lmer))


#plotting
plwl=ggplot(lo,aes(x=Treatment,y=1000*Diff_Abd_dry_weight,fill= Treatment))+
  geom_boxplot()+
  labs(y="Abdominal Lipid per Microcolony (mg)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave &
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08), size=2.25) +
  scale_fill_brewer(palette = "Set2")+ng4;plwl

ggsave("20240802 average abdominal lipids.jpg", plot = plwl,  height=25,width=28,units = c("cm"))

#getting medians
abd_lipids_meds <- ddply(lo, .(Treatment), summarise, med = median(Diff_Abd_dry_weight*1000))
#Treatment    med
#1        CN 2.1050
#2        HH 2.4580
#3        HW 2.3850
#4      HWHH 1.9475

#getting number of samples
ddply(lo, .(Treatment), summarise, 
      n = length(Diff_Abd_dry_weight) )
#Treatment n
#1        CN 9
#2        HH 7
#3        HW 7
#4      HWHH 7

#making table with SE and means
abd_lipids <- ddply(lo, .(Treatment), summarise, 
                    M = mean(Diff_Abd_dry_weight*1000), SE = sd(Diff_Abd_dry_weight*1000) / sqrt((length(Diff_Abd_dry_weight*1000))), 
                    SD = sd(Diff_Abd_dry_weight*1000))
abd_lipids
#Treatment        M        SE        SD
#1        CN 2.028611 0.2130577 0.6391732
#2        HH 2.373286 0.3237735 0.8566243
#3        HW 2.322357 0.1990962 0.5267589
#4      HWHH 1.879405 0.1997639 0.5285256

#---------------------------------------  will now look at thoracic
#none three categories of weight - wet, initial dry, final dry were not normal, linear, or homoesdatic
#will remove outliers to see if it works

#removing outliers
summary(l$Diff_Thx_dry_weight)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0001700 0.0004925 0.0007875 0.0010200 0.0012225 0.0060025 

IQR(l$Diff_Thx_dry_weight*1000) #0.00073

#getting min and max lipid values for thresholds
lamin = 0.0004925-(1.5*0.00073) 
lamax =  0.0012225+(1.5*0.00073) 

loa<-l[!(l$Diff_Thx_dry_weight<lamin | l$Diff_Thx_dry_weight>lamax),]
View(loa)

#variance - homo p-value = 0.4401
bartlett.test(loa$Diff_Thx_dry_weight~loa$Treatment,
              data=loa)

#linearity
modeloa<-lm(Diff_Thx_dry_weight~Treatment, data=loa)
plot(modeloa,1) #linear

#doing a lmm
lipidsoa_lmer<-lmer(loa$Diff_Thx_dry_weight*1000~loa$Temperature*loa$Humidity+
                      loa$Thx_Wet_weight+(1|loa$Colony_source))
summary(lipidsoa_lmer)

#Estimate Std. Error t value
#(Intercept)                         -0.1153     0.4613  -0.250
#loa$TemperatureYes                   0.2338     0.1711   1.366
#loa$HumidityYes                      0.1507     0.1648   0.915
#loa$Thx_Wet_weight                  16.4312     8.3342   1.972
#loa$TemperatureYes:loa$HumidityYes  -0.5569     0.2442  -2.280

car::qqPlot(resid(lipidsoa_lmer)) #look good
plot(resid(lipidsoa_lmer)~fitted(lipidsoa_lmer))  #look good
hist(resid(lipidsoa_lmer)) #look good

Anova(lipidsoa_lmer, type = 3, test.statistic="F") 
#                                 F Df Df.res  Pr(>F)  
#(Intercept)                  0.0614  1 25.075 0.80638  
#loa$Temperature              1.7865  1 24.969 0.19340  
#loa$Humidity                 0.8298  1 24.199 0.37132  
#loa$Thx_Wet_weight           3.8172  1 24.388 0.06229 .
#loa$Temperature:loa$Humidity 5.0931  1 24.441 0.03323 

#------------------------
thx_lipids_meds <- ddply(loa, .(Treatment), summarise, med = median(Diff_Thx_dry_weight*1000))
thx_lipids_meds
#Treatment      med
#1        CN 0.777500
#2        HH 1.030000
#3        HW 1.166667
#4      HWHH 0.502500

#getting number of samples
ddply(loa, .(Treatment), summarise, 
      n = length(Diff_Thx_dry_weight) )
#Treatment n
#1        CN 9
#2        HH 8
#3        HW 7
#4      HWHH 7

#making table with SE and means
thx_lipids <- ddply(loa, .(Treatment), summarise, 
                    M = mean(Diff_Thx_dry_weight*1000), SE = sd(Diff_Thx_dry_weight*1000) / sqrt((length(Diff_Thx_dry_weight*1000))), 
                    SD = sd(Diff_Thx_dry_weight*1000))
thx_lipids
#Treatment         M         SE        SD
#1        CN 0.7667593 0.08848829 0.2654649
#2        HH 0.9432500 0.14205721 0.4017985
#3        HW 0.9682143 0.17278286 0.4571405
#4      HWHH 0.5770238 0.10541786 0.2789094


plwloa=ggplot(loa,aes(x=Treatment,y=1000*Diff_Thx_dry_weight,fill= Treatment))+
  geom_boxplot()+
  labs(y="Mean Thoracic Lipids per Microcolony (mg)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave &
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08), size=2.25) +
  scale_fill_brewer(palette = "Set2")+ng4;plwloa


ggsave("20240802 average thoracic lipids.jpg", plot = plwloa,  height=25,width=28,units = c("cm"))

#-----------------------------------------------Drone Numbers per Microcolonies

drone<-read.csv("drone weight.csv")
drone$Treatment<-as.factor(drone$Treatment)
drone$Temperature<-as.factor(drone$Temperature)
drone$Humidity<-as.factor(drone$Humidity)
drone$Colony_source<-as.factor(drone$Colony_source)
drone$Microcolony_id<-as.factor(drone$Microcolony_id)

#removing unecessary columns
d1<-subset(drone, select = -c(color,Notes))
summary(d1)

#making a table with drone counts 
drone_count<-dplyr::count(d1, Treatment,Temperature,Humidity, Colony_source, Microcolony_id) 
drone_count$Treatment<-as.factor(drone_count$Treatment)
drone_count$Colony_source<-as.factor(drone_count$Colony_source)
drone_count$Microcolony_id<-as.factor(drone_count$Microcolony_id)
View(drone_count)

#making table with SE and means
dmeans <- ddply(drone_count, .(Treatment), summarise, 
                M = mean(n), SE = sd(n) / sqrt((length(n))), 
                SD = sd(n))

View(dmeans)

#summarizing
drone_count=(ddply(drone, .(Treatment, Temperature, Humidity,Microcolony_id,Colony_source), summarise, 
                   n = length(Colony_source) ))
view(drone_count)
plot(drone_count$n~drone_count$Treatment)

#checking for linearity
model8<-lm(n~Treatment,data=drone_count) 
plot(model8,1)#pretty linear

#checking for normality
hist(drone_count$n) #normal

#checking for variance
bartlett.test(drone_count$n~drone_count$Treatment,
              data=drone_count) # p-value = 0.6352, no variance

#do a poisson family with glmer!

######### mixed effects model with colony source as random effect
droneglmer<-glmer(drone_count$n~drone_count$Temperature*drone_count$Humidity + 
                    (1|drone_count$Colony_source),
                  family = poisson)
summary(droneglmer)
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                         2.49412    0.09578  26.039   <2e-16 ***
#drone_count$TemperatureYes                          0.26448    0.12734   2.077   0.0378 *  
#drone_count$HumidityYes                            -0.05178    0.14158  -0.366   0.7146    
#drone_count$TemperatureYes:drone_count$HumidityYes -0.35545    0.19746  -1.800   0.0718 .  

plot(resid(droneglmer)~fitted(droneglmer)) 
hist(resid(droneglmer))
car::qqPlot(resid(droneglmer))

Anova(droneglmer, type = 3)
#Chisq Df Pr(>Chisq)    
#(Intercept)                                  678.0509  1    < 2e-16 ***
#drone_count$Temperature                        4.3134  1    0.03781 *  
#drone_count$Humidity                           0.1337  1    0.71458    
#drone_count$Temperature:drone_count$Humidity   3.2404  1    0.07184 . 


#box plot with drones emerged to look into difference
plw3=ggplot(drone_count,aes(x=Treatment,y=n,fill= Treatment))+
  geom_boxplot()+
  labs(y="Average Drones Emerged per Microcolony",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave &
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08), size=2.25) +
  scale_fill_brewer(palette = "Set2")+ng4;plw3

ggsave("20240802 drone emergence.jpg", plot = plw3,   height=25,width=28,units = c("cm"))


#-------------------------------------------------------------Microcolony and  Drones Weight
#getting average drone weights per microcolony
dweights=ddply(d1, .(Treatment,Temperature, Humidity, Microcolony_id,Colony_source), numcolwise(mean),na.rm=T)
View(dweights)

dweights_<-merge(dweights,dev1)
View(dweights_)

#checking normality, linearity and variance
hist(dweights$Weight_emergence) #normal but doing a shapiro to confirm
shapiro.test(dweights$Weight_emergence) #normal
bartlett.test(dweights$Weight_emergence~dweights$Treatment)#p-value = 0.7932
fligner.test(Weight_emergence ~ Treatment, data = dweights) #  p-value = 0.8054

model7<-lm(Weight_emergence~Treatment,data=dweights) 
plot(model7,1) # linear!
(1|dev1$total_drones)

#taking into to account all the number of drones in the microcolony, as workers would have to divide provisions to all drones
dweight_lmer<-lmer(dweights_$Weight_emergence*1000~dweights_$Temperature*dweights_$Humidity+
                     dweights_$total_drones+(1|dweights_$Colony_source) )

summary(dweight_lmer) 
#Estimate Std. Error t value
#(Intercept)                                    Estimate Std. Error t value
#(Intercept)                                    126.5979    50.6234   2.501
#dweights_$TemperatureYes                        -1.1486     6.7172  -0.171
#dweights_$HumidityYes                           -7.6086     7.6548  -0.994
#dweights_$total_drones                           0.6772     1.6435   0.412
#dweights_$TemperatureYes:dweights_$HumidityYes  14.0243    11.4692   1.223


Anova(dweight_lmer,type=3, test.statistic="F")
#                                         F Df Df.res  Pr(>F)  
#(Intercept)                              6.0444  1 27.959 0.02041 *
#dweights_$Temperature                    0.0285  1 27.631 0.86711  
#dweights_$Humidity                       0.9839  1 27.111 0.33002  
#dweights_$total_drones                   0.1640  1 27.862 0.68855  
#dweights_$Temperature:dweights_$Humidity 1.4834  1 27.222 0.23370  

#looking at residuals 
plot(resid(dweight_lmer)~fitted(dweight_lmer)) #same problem as before
hist(resid(dweight_lmer)) #does NOT look good
car::qqPlot(resid(dweight_lmer)) #look oo

#getting number of samples
ddply(dweights_, .(Treatment), summarise, 
      n = length(Weight_emergence) )
#Treatment n
#1        CN 9
#2        HH 8
#3        HW 9
#4      HWHH 8

#making table with SE and means
dweights_mse <- ddply(dweights_, .(Treatment), summarise, 
                      M = mean(Weight_emergence *1000), SE = sd(Weight_emergence*1000) / sqrt((length(Weight_emergence*1000))), 
                      SD = sd(Weight_emergence*1000))
#Treatment        M       SE       SD
#1        CN 147.3647 4.881139 14.64342
#2        HH 138.3825 4.219388 11.93423
#3        HW 145.8790 5.000393 15.00118
#4      HWHH 153.4183 4.921588 13.92035


#plotting drone weight - in which heat wave was statistically different from control and high humdity
plw5=ggplot(dweights, aes(x=Treatment,y=Weight_emergence *1000,fill= Treatment))+
  geom_boxplot()+
  labs(y="Average Drone Weight at Emergence (mg)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave & High Humidity"))+ #assigning labels for x-axis
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave
    & High Humidity"))+
  geom_jitter()+
  scale_fill_brewer(palette = "Set2")+ng4;plw5

ggsave("20240208 drone weight at emergence.jpg", plot = plw5,  height=25,width=28,units = c("cm"))

#-------------------------------------------------------drone development time
dev<-read.csv(" Developmental_time.csv") #developmental time

#transforming into proper format
dev$Microcolony_id<-as.factor(dev$Microcolony_id)
dev$Colony_source<-as.factor(dev$Colony_source)
dev$Treatment<-as.factor(dev$Treatment)

#separating drone per micro = drone_micro
drone_micro<-drone %>% group_by(Microcolony_id, .drop = FALSE) %>% count() 
print(drone_micro,n=34)
View(drone_micro)

#merging to see total drones per microcolony, as well as the developmental time for one
dev1<-merge(dev,drone_micro)
View(dev1)

##renaming the column n
colnames(dev1)[7] <- "total_drones"

#looking into normallity
hist(dev1$Development_time) #slightly right skewed
shapiro.test(dev1$Development_time) #p value is 0.0823

model20<-lm(Development_time~Temperature*Humidity,data=dev1)  #most simple model
plot(model20,1) # linear

bartlett.test(dev1$Development_time~dev1$Treatment,data=dev1) #p-value = 0.9405

#will use non transformed data since only normality is slightly off but this isn't that importat

#including colony source
develop_lmer<-lmer(dev1$Development_time~dev1$Temperature*dev1$Humidity+ 
                     (1|dev1$Colony_source))
summary(develop_lmer) 
# Estimate Std. Error t value
#(Intercept)                           26.0000     0.5641  46.087
#dev1$TemperatureYes                    0.7778     0.7978   0.975
#dev1$HumidityYes                       1.1250     0.8224   1.368
#dev1$TemperatureYes:dev1$HumidityYes   1.3472     1.1630   1.158

Anova(develop_lmer,type=3,test.statistic="F")
#F Df Df.res Pr(>F)    
#(Intercept)                    2124.0407  1 18.107 <2e-16 ***
#dev1$Temperature                  0.9314  1 28.533 0.3426    
#dev1$Humidity                     1.8554  1 28.225 0.1839    
#dev1$Temperature:dev1$Humidity    1.3417  1 28.013 0.2565   

hist(resid(develop_lmer)) #look good
car::qqPlot(resid(develop_lmer)) #looks good


#----------------------------- given state of residuals in the plot, will try a more complex model
#and choose one doing model selection

#If number of drones affects pollen provision, maybe it will affect development time?
develop_drones_lmer<-lmer(dev1$Development_time~dev1$Temperature*dev1$Humidity+
                            dev1$total_drones+
                            (1|dev1$Colony_source))
summary(develop_drones_lmer) 
#                                       Estimate Std. Error t value
#(Intercept)                           -0.4421     3.6756  -0.120
#dev1$TemperatureYes                    1.1610     0.4878   2.380
#dev1$HumidityYes                       2.8854     0.5561   5.189
#dev1$total_drones                      0.8622     0.1193   7.225
#dev1$TemperatureYes:dev1$HumidityYes  -1.8383     0.8331  -2.207


Anova(develop_drones_lmer,type=3,test.statistic="F")
#                                   F Df Df.res    Pr(>F)    
#(Intercept)                    275.7147  1 27.988 5.053e-16 ***
#dev1$Temperature                 5.6966  1 27.660   0.02409 *  
#dev1$Humidity                   26.0926  1 27.118 2.251e-05 ***
#dev1$total_drones               46.8038  1 27.898 2.003e-07 ***
#dev1$Temperature:dev1$Humidity   4.9425  1 27.232   0.03469 *  
---               

#look good 
plot(resid(develop_drones_lmer)~fitted(develop_drones_lmer))
hist(resid(develop_drones_lmer))
car::qqPlot(resid(develop_drones_lmer))

#plotting developmental time - in which heat wave was satitically different from control and high humdity
plw13=ggplot(dev1,aes(x=Treatment,y=Development_time,fill= Treatment))+
  geom_boxplot()+
  labs(y="First Drone Developmental Time (days)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave &
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08), size=2.25) +
  scale_fill_brewer(palette = "Set2")+ng4;plw13

ggsave("20240802 dev time.jpg", plot = plw13,  height=25,width=28,units = c("cm"))

dev_meds <- ddply(dev1, .(Treatment), summarise, med = median(Development_time))
dev_meds
#Treatment  med
#1        CN 26.0
#2        HH 27.0
#3        HW 26.0
#4      HWHH 29.5
ddply(dev1, .(Treatment), summarise, 
      n = length(Development_time) )
#Treatment n
#1        CN 9
#2        HH 8
#3        HW 9
#4      HWHH 8

ddply(dev1, .(Treatment), summarise, 
      M = mean(Development_time), SE = sd(Development_time) / sqrt((length(Development_time))), 
      SD = sd(Development_time))
#Treatment        M        SE       SD
#1        CN 26.00000 0.5773503 1.732051
#2        HH 27.12500 0.5805632 1.642081
#3        HW 26.77778 0.4937886 1.481366
#4      HWHH 29.25000 0.6748016 1.908627

#model selection

#getting AIC to compare models
dev_models <- list(develop_drones_lmer, develop_lmer)

#specify model names
dev_mod.names <- c('develop_drones_lmer', 'develop_lmer')

#calculate AIC of each model
aictab(cand.set = dev_models, modnames = dev_mod.names)

#Model selection based on AICc:

#K   AICc Delta_AICc AICcWt Cum.Wt Res.LL
#develop_drones_lmer 7 115.11       0.00      1      1 -48.40
#develop_lmer        6 140.37      25.26      0      1 -62.63


anova(develop_lmer,develop_drones_lmer) 
#best model is the that includes number of drones!

r.squaredGLMM(develop_lmer)
#         R2m       R2c
#[1,] 0.3376558 0.3376558

r.squaredGLMM(develop_drones_lmer)
#         R2m       R2c
#[1,] 0.7475985 0.7475985

#-------------------------------------------------------growth rate
micro_gr<-read.csv("Trial2_gr.csv")
micro_gr$Treatment<-as.factor(micro_gr$Treatment)
micro_gr$Microcolony_id<-as.factor(micro_gr$Microcolony_id)
micro_gr$Colony_source<-as.factor(micro_gr$Colony_source)
micro_gr$Temperature<-as.factor(micro_gr$Temperature)
micro_gr$Humidity<-as.factor(micro_gr$Humidity)

#creating a new column for growth rate with diff weight/time passed
gr2=micro_gr$Diff_weight/micro_gr$Time_elapsed_days
mgr2<-cbind(micro_gr,gr2)
View(mgr2)
mgr2 <- na.omit(mgr2)

#looking at assumptions
hist(mgr2$gr2) #does not look good

#looking at linearity of data
model4<-lm(gr2~Treatment,data=mgr2) 
plot(model4,1) #trend line does NOT looks linear

#seeing whether sqrt transformation workers for both linearity and normality
mgr2$gr_sqrt<-sqrt(mgr2$gr2)
hist(mgr2$gr_sqrt)
shapiro.test(mgr2$gr_sqrt) #transformation works!  p-value = 0.5061

model5<-lm(gr_sqrt~Temperature*Humidity,data=mgr2) 
plot(model5,1) #trend line looks linear

#checking variance of data
bartlett.test(gr_sqrt~interaction(Temperature,Humidity),
              data=mgr2)  #no significant variance, p-value = 0.4997

#given the state of the data and that no transformations work, I am removing outliers
mgr2$Humidity<-as.factor(mgr2$Humidity)

mm23<-lmer(mgr2$gr_sqrt~mgr2$Temperature*mgr2$Humidity+(1|factor(mgr2$Colony_source)) )
summary(mm23)

plot(resid(mm23)~fitted(mm23)) #look bad!
hist(resid(mm23))#look kinda bad
shapiro.test(resid(mm23))# p-value = 0.1014
car::qqPlot(resid(mm23)) #looks kinda good


#--removing outliers
summary(mgr2$gr2)
IQR(mgr2$gr2)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.2886  0.1629  0.2400  0.2373  0.3420  0.4655 

#getting min and max lipid values for thresholds
mgrmin = 0.1629-(1.5*0.1791071) 
mgrmax =  0.3420+(1.5*0.1791071) 

# removing outlier
mgro<-mgr2[!(mgr2$gr2<mgrmin | mgr2$gr2>mgrmax),]
View(mgro)

#merging with drone counts to account for number of drones produced in this model
mgro_<-merge(mgro,drone_count)
View(mgro_)

#looking at assumptions
hist(mgro_$gr2)
shapiro.test(mgro_$gr2) # data is normal p-value = 0.429

#looking at linearity of data
modelmgro<-lm(gr2~Temperature*Humidity,data=mgro_) 
plot(modelmgro,1) #not really linear - shoudl I chooose an alternative model?

#checking variance of data
bartlett.test(mgro_$gr2~interaction(mgro_$Temperature,mgro_$Humidity),
              data=mgro_)  #no significant variance, p-value = 0.3417


#testing models
mgr_n_lmer<- lmer(mgro_$gr2~mgro_$Temperature*mgro_$Humidity+
                    (1|mgro_$Colony_source)+(1|mgro_$n))

plot(resid(mgr_lmer)~fitted(mgr_lmer)) #even including total number of drones it's not working
#i think its not working because data is NOT linear

#simple model
mgr_lmer<- lmer(gr2~Temperature*Humidity+
                  (1|Colony_source),data=mgro)
summary(mgr_lmer)
#                           Estimate Std. Error t value
#(Intercept)                 0.20698    0.02076   9.972
#TemperatureYes              0.17233    0.02935   5.871
#HumidityYes                 0.06444    0.03026   2.130
#TemperatureYes:HumidityYes -0.31396    0.04359  -7.202

Anova(mgr_lmer, type = 3, test.statistic="F") 
# F Df Df.res    Pr(>F)    
#(Intercept)                      96.578  1 16.614 2.472e-08 ***
#mgro_$Temperature                32.333  1 26.215 5.391e-06 ***
#mgro_$Humidity                    4.366  1 26.215   0.04651 *  
#mgro_$Temperature:mgro_$Humidity 49.197  1 26.341 1.765e-07 ***


#plotting microcolony growth rates
plwmgro=ggplot(mgro,aes(x=Treatment,y=gr2,fill=Treatment))+
  geom_boxplot()+
  labs(y="Microcolony Growth Rate (day^-1)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08),size=2.25) +
  geom_hline(yintercept=0, linetype="dashed")+
  scale_fill_brewer(palette = "Set2")+ng4;plwmgro


ggsave("20240802 microgr without outliers.jpg", plot = plwmgro,   height=25,width=28,units = c("cm"))

mgr_meds <- ddply(mgro, .(Treatment), summarise, med = median(gr2))
#Treatment       med
#1        CN 0.2057143
#2        HH 0.2928571
#3        HW 0.3965517
#4      HWHH 0.1314286


ddply(mgro, .(Treatment), summarise, 
      n = length(gr2) )
#Treatment n
#1        CN 9
#2        HH 8
#3        HW 9
#4      HWHH 7

gr_means <- ddply(mgro, .(Treatment), summarise, 
                  M = mean(gr2), SE = sd(gr2) / ((length(gr2))), 
                  SD = sd(gr2))
#Treatment         M          SE         SD
#1        CN 0.2069841 0.004981981 0.04483783
#2        HH 0.2714286 0.010567123 0.08453698
#3        HW 0.3793103 0.006962752 0.06266476
#4      HWHH 0.1297959 0.007123812 0.04986668





