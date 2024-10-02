##Microcolony analysis

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
library(ggplot2)
library(dplyr)
library(tidyverse) 
library(emmeans)
library(MuMIn)

#Queen survival probability

#loading
qsurv<-read.csv("queen_survival.csv")

#excluding colonies after death
qt3<-read.csv("queen_survival - try3.csv")

qt3survfit_cox<-coxph(Surv(Days_since_start, Queen_status) ~ Temperature*Humidity, data=qt3)
print(qt3survfit_cox)
#                             coef exp(coef) se(coef)      z       p
#TemperatureYes              -2.0697    0.1262   1.1536 -1.794 0.07280
#HumidityYes                  1.5477    4.7008   1.0343  1.496 0.13454
#TemperatureYes:HumidityYes   5.1702  175.9548   1.7375  2.976 0.00292

#the above are the results for the analysis, and unfortunately not able to plot. 
#after some reading, i can trust the plot below for my analysis, but have to use  treatment  because "interaction terms are not valid for this function"
qt3survfit <- survfit(Surv(Days_since_start, Queen_status) ~ Treatment, data=qt3)
summary(qt3survfit)
autoplot(qt3survfit)

#plotting with ggplot not working cox - error is below
#Error in ggsurvplot(qt3survfit_cox, conf.int = TRUE, risk.table = TRUE,  :  object 'ggsurv' not found
ggsurvplot(qt3survfit,conf.int = TRUE,risk.table = TRUE,
           risk.table.col="strata",
           linetype = "strata",
           surv.median.line="hv",
           ggtheme=theme_bw(),
           legend.labs =c("Control", "High Humidity","Heat Wave", "Heat Wave & High Humidity"),
           xlab = "Time in days",
           cumevents.title = "Percentage of alive workers",data=qt3)

#colony foraging

#loading traffic data file
traffic_w<-read.csv("colony_traffic_copy.csv")
traffic_w$Date = as.Date(traffic_w$Date,format='%m/%d/%Y')
traffic_w$julian = yday(traffic_w$Date)
traffic_w

#removing colonies that died early
traffic_w<-traffic_w[!(traffic_w$Colony_id==14 | traffic_w$Colony_id==13),]
na.omit(traffic_w)
View(traffic_w)

#making a table of number of flights per colony overall
tw_total=(ddply(traffic_w, .(Treatment,Colony_id,Temperature,Humidity,Room), summarise, 
                n = length(Colony_id) ))
tw_total<-na.omit(tw_total)

#getting the average number of workers throughout experiment
traffic_workers=ddply(traffic_w, .(Colony_id,Temperature,Temperature,Humidity), 
                      numcolwise(mean),na.rm=T)
View(traffic_workers)

#making table with the total numebr of flights with the average number of workers and 
#average weight across experiment
t<-merge(tw_total,traffic_workers)
View(t)

#putting things into proper format
t$Treatment<-as.factor(t$Treatment)
t$Window<-as.factor(t$Window)
t$Room<-as.factor(t$Room)
t$Period<-as.factor(t$Period)
t$Weather<-as.factor(t$Weather)
t$Colony_id<-as.factor(t$Colony_id)
t$Pollen<-as.factor(t$Pollen)
t$Direction<-as.integer(t$Direction)
t$Temperature<-as.factor(t$Temperature)
t$Humidity<-as.factor(t$Humidity)

#looking at assumptions
hist(t$n) #not normal
hist(log(t$n)) #looks better
shapiro.test(log(t$n)) # p-value 0.3362

#testing variance, homo!
bartlett.test(log(n)~interaction(Temperature,Humidity),data=t) 
#p-value = 0.8015

#checking linearity of simplest model
modelt<-lmer(log(n)~Temperature*Humidity+(1|Room), data=t) 
car::qqPlot(resid(modelt)) #looks ok

#----------------------------------------------------------------------
#from raw data, isolating trips to when all colonies were still alive 
#(earliest colony death was at age 15 from HWHH) 
tw_period_ <- traffic_w[-which(traffic_w$Days_Since_Start > 15),]
View(tw_period_)

#summarizing number of trips per colony during that time
tw_period_a=(ddply(tw_period_, .(Temperature,Humidity, Treatment,Colony_id, Room), summarise, 
                   n = length(Colony_id) ))
View(tw_period_a)

#getting average number of workers and colony weight at first 15 days
traffic_workers_15=ddply(traffic_w, .(Colony_id,Temperature,Temperature,Humidity), 
                         numcolwise(mean),na.rm=T)
traffic_workers_15

#merging it so we have an idea with the number of workers
tp_15<-merge(tw_period_a,traffic_workers_15)
tp_15<-na.omit(tp_15)
View(tp_15)

#---------------------- normalizing data sets - all using 15 days datasets 
#this, I will divide the number of trips (n)/num of workers to see how many trips per "worker"
tp_15<-tp_15%>%mutate(pt=n/Num_workers)
hist(tp_15$pt) #rightly skewed
shapiro.test(tp_15$pt) #p-value = 0.08433

#analyzing with lmer bc isn't that strict with nearlity and residuals look better 
pt_lmer<-lmer(pt~Temperature*Humidity+(1|Room),data=tp_15)
car::qqPlot(resid(pt_lmer)) #looks good
plot(resid(pt_lmer)~fitted(pt_lmer)) #there might be clumping
hist(resid(pt_lmer)) #looks weird  - 

summary(pt_lmer)
#                           Estimate Std. Error t value
#(Intercept)                  3.0463     0.9005   3.383
#TemperatureYes               0.3868     0.7473   0.518
#HumidityYes                  1.0721     0.6535   1.641
#TemperatureYes:HumidityYes  -1.8615     0.8850  -2.103

Anova(pt_lmer,type=3, test.statistic = "F")
#                           F Df Df.res  Pr(>F)  
#(Intercept)          11.1592  1 3.0510 0.04328 *
#Temperature           0.2371  1 6.5976 0.64205  
#Humidity              2.6272  1 6.1086 0.15529  
#Temperature:Humidity  4.4241  1 6.0019 0.08009 .

pt=ggplot(tp_15,aes(x=Treatment,y=pt,fill= Treatment))+
  geom_boxplot()+
  labs(y="Colony Traffic per Worker in 15 days",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave &
High Humidity"))+
  geom_jitter(position = position_jitter(height = .03, width = .08)) +
  scale_fill_brewer(palette = "Set2")+ng4;pt
ggsave("20240306 per worker fifteen days traffic.jpg", plot = pt,   height=25,width=28,units = c("cm"))


#getting number of samples
ddply(tp_15, .(Treatment), summarise, 
      n = length(pt) )
#Treatment n
#1        CN 3
#2        HH 3
#3        HW 3
#4      HWHH 3

#making table with SE and means
ddply(tp_15, .(Treatment), summarise, 
      M = mean(pt), SE = sd(pt) / sqrt((length(pt))), 
      SD = sd(pt))
#Treatment        M        SE        SD
#1        CN 2.827827 0.8707575 1.5081962
#2        HH 3.330698 0.4850908 0.8402019
#3        HW 4.220904 0.2177515 0.3771566
#4      HWHH 2.862308 1.1094184 1.9215690


#---------------------------------------------------------colony growth rate
colony_gr<-read.csv("Colony_Growth_Rate.csv")

#putting data into correct format
colony_gr$Room<-as.factor(colony_gr$Room)
colony_gr$Treatment<-as.factor(colony_gr$Treatment)
colony_gr$Colony_id<-as.factor(colony_gr$Colony_id)

#creating a new column for growth rate with diff weight/time passed
gr=colony_gr$Diff_weight/colony_gr$Time_elapsed_days
colony_gr2<-cbind(colony_gr,gr)
View(colony_gr2)

#looking at nested anova assumptions
shapiro.test(colony_gr2$gr) #data is normal
bartlett.test(gr~interaction(Temperature,Humidity), data=colony_gr2)  #no variance, p-value = 0.4982

#looking at linearity
model14<-lm(gr~Temperature*Humidity,data=colony_gr2)
plot(model14,1) #looks linear

#colony 
cgr_lmer<-lmer(colony_gr2$gr~colony_gr2$Temperature*colony_gr2$Humidity+(1|colony_gr2$Room))
summary(cgr_lmer)  #this summary and the anova are giving different results
#Estimate Std. Error t value
#(Intercept)                                       -0.4798     1.0157  -0.472
#colony_gr2$TemperatureYes                          1.0681     0.8007   1.334
#colony_gr2$HumidityYes                             0.1412     0.6724   0.210
#colony_gr2$TemperatureYes:colony_gr2$HumidityYes  -1.6993     0.9166  -1.854

Anova(cgr_lmer, type=3,test.statistic="F")
#F Df Df.res Pr(>F)
#(Intercept)                                0.2189  1 2.9433 0.6723
#colony_gr2$Temperature                     1.6038  1 7.6024 0.2428
#colony_gr2$Humidity                        0.0428  1 7.1670 0.8419
#colony_gr2$Temperature:colony_gr2$Humidity 3.4349  1 7.0052 0.1062

#i think it looks good 
shapiro.test(resid(cgr_glmer)) #good enough
car::qqPlot(resid(cgr_glmer)) #look good



#plotting growth rates
plw7=ggplot(colony_gr2,aes(x=Treatment,y=gr,fill=Treatment))+geom_boxplot()+
  labs(y="Colony Growth Rate (g)",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave 
    & High Humidity"))+ 
  geom_jitter()+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_fill_brewer(palette = "Set2")+ng4;plw7

ggsave("20240211 colony growth rate.jpg", plot = plw7,   height=25,width=28,units = c("cm"))

ddply(colony_gr2, .(Treatment), summarise, 
      M = mean(gr), SE = sd(gr) / sqrt((length(gr))), 
      SD = sd(gr))

#Treatment          M        SE        SD
#1        CN -1.0962963 0.7285007 1.2618003
#2        HH -1.3472222 0.3173766 0.6347531
#3        HW  1.5185185 0.8522382 1.4761199
#4      HWHH -0.3532922 1.0960703 1.8984494


#------------------------------------colony workers
performance<-read.csv("performance.csv")
performance$Date = as.Date(performance$Date,format='%m/%d/%Y')
performance$julian = yday(performance$Date)

#removing unnecessary "average" column which was used to indicate that some of the 
#data, such as colony weight, number of workers and brood were averaged. 
p1<-subset(performance, select = -Averaged.Column)
summary(p1)

#Transforming data to appropriate format
p1$Room<-as.factor(p1$Room)
p1$Treatment<-as.factor(p1$Treatment)
View(p1)
p1$Weather<-as.factor(p1$Weather)
p1$Colony_id<-as.factor(p1$Colony_id)
p1$Humidity<-as.factor(p1$Humidity)
p1$Temperature<-as.factor(p1$Temperature)
p1$Num_workers<-as.numeric(p1$Num_workers)
p1$Weight_g<-as.numeric(p1$Weight_g)
p1$Num_emerged_drones<-as.numeric(p1$Num_emerged_drones)
p1$Num_emerged_gynes<-as.numeric(p1$Num_emerged_gynes)
p1$Num_gyne_cells<-as.numeric(p1$Num_gyne_cells)

View(p1)

#Will remove colonies that were substituted early because of my malpractice or unknown reasons
p2<-p1[!(p1$Colony_id==14 | p1$Colony_id==13),]

#p2 therefore is the complete table!

#getting an average of workers throughout experiment
p3=ddply(p2, .(Treatment,Colony_id,Room, Temperature, Humidity), numcolwise(mean),na.rm=T)
View(p3)

#will transform data to try to make it linear
hist(p3$Num_workers)
shapiro.test(p3$Num_workers)# non normalp-value = 0.05222
p3$Num_workers_log<-log(p3$Num_workers)
shapiro.test(p3$Num_workers_log)  #p-value = 0.4609

pworkers_glmer<-glmer(p3$Num_workers_log~p3$Temperature*p3$Humidity+(1|factor(p3$Room)))

summary(pworkers_glmer) 
#Estimate Std. Error t value
#(Intercept)                        51.113     12.852   3.977
#p3$TemperatureYes                  40.804     12.576   3.245
#p3$HumidityYes                      8.121     11.133   0.729
#p3$TemperatureYes:p3$HumidityYes  -45.881     15.117  -3.035


Anova(pworkers_glmer, test.statistic="F")
#                                F Df Df.res  Pr(>F)   
#p3$Temperature             2.4012  1 7.3841 0.16296  
#p3$Humidity                3.0532  1 6.3331 0.12857  
#p3$Temperature:p3$Humidity 9.2112  1 6.0047 0.02293 *
---
  
#box plot with drones emerged to look into difference
plw4=ggplot(p3,aes(x=Treatment,y=Num_workers,fill= Treatment))+
  geom_boxplot()+
  labs(y="Average Number of Workers per Colony",x="Treatment")+
  scale_x_discrete(labels=c("CN" = "Control", "HW" = "Heat Wave","HH"="High Humidity",
                            "HWHH"="Heat Wave
    & High Humidity"))+
  geom_jitter()+
  annotate(geom='text',x=-Inf, y=-Inf,hjust=-.1,vjust=-.5,size=8, label= 'p<0.001')+
  scale_fill_brewer(palette = "Set2")+ng4;plw4

ggsave("20231210 workers count.jpg", plot = plw4,   height=25,width=28,units = c("cm"))



#--------------------------gyne and drone analysis
p2_<-subset(p2, select = -c(Temperature.1,Realtive.Humidity,Weather,
                            External_Temp_C,ExternaL_Hum))
View(p2_)

p3=ddply(p2_, .(Treatment,Colony_id,Room, Temperature, Humidity), numcolwise(mean),na.rm=T)

#getting mean gynes
p_no0_g<-filter(p2_,Num_emerged_gynes>0)
View(p_no0_g)

mean_gynes<-ddply(p_no0_g, .(Treatment,Colony_id,Room, Temperature, Humidity), numcolwise(mean),na.rm=T)
mean_gynes 
#Colony_id Num_emerged_gynes
#17           2.666667
#25           2.333333
#22           9.000000
#26           13.333333

#----------------getting mean gyne cells
p_no0_gc<-filter(p2_,Num_gyne_cells>0)
View(p_no0_gc)

mean_gynes_cells<-ddply(p_no0_gc, .(Treatment,Colony_id,Room, Temperature, Humidity), numcolwise(mean),na.rm=T)
mean_gynes_cells

#Colony_id  Num_gyne_cells
#22           11.16667  
#26           18.73810


#---------------------------getting mean drones since I'm already here
p_no0_d<-filter(p2_,Num_emerged_drones>0)
View(p_no0_d)

mean_colony_drones<-ddply(p_no0_d, .(Treatment,Colony_id,Room, Temperature, Humidity), numcolwise(mean),na.rm=T) 
mean_colony_drones
#Treatment Colony_id Num_emerged_drones
#1        CN        17    1.000000
#2        CN        25    1.000000
#3        CN        27    1.000000
#4        HH        24    2.000000
#5        HW        22    1.766667
#6        HW        26   3.277778

#will all of the information above, I made a new table with mean gyne cell, gyne, and drone produciton in each colony!

dg<-read.csv("drone_gyne.csv")
View(dg)

#---------------------------------------------------- first putting new dataset into correct format
dg$Room<-as.factor(dg$Room)
dg$Treatment<-as.factor(dg$Treatment)
dg$Colony_id<-as.factor(dg$Colony_id)
dg$Humidity<-as.factor(dg$Humidity)
dg$Temperature<-as.factor(dg$Temperature)


#will do a binomial glm - so first need to change 'yes' and 'no' to '0' and '1'
dg$Gyne_cell_status<-ifelse(dg$Gyne_cell_status=="Yes",1,0)
dg$Gyne_status<-ifelse(dg$Gyne_status=="Yes",1,0)
dg$Drone_stauts<-ifelse(dg$Drone_stauts=="Yes",1,0)

dg$Gyne_cell_status<-as.integer(dg$Gyne_cell_status)
dg$Gyne_status<-as.integer(dg$Gyne_status)
dg$Drone_stauts<-as.integer(dg$Drone_stauts) #notice the mistype!

#-------------------------------------------------checking assumptions
#normality
hist(dg$Mean_gyne_cells) 
hist(dg$Mean_gynes)
hist(dg$Mean_drones)

#none of them are normal, normal with count data - will try transforming

#linearity
modelgynec<-lm(Mean_gyne_cells~Temperature*Humidity, data=dg)
summary(modelgynec)
plot(modelgynec,1) #linear

modeldronec<-lm(Mean_drones~Temperature*Humidity, data=dg)
plot(modeldronec,1) #linear


#variance
bartlett.test(Mean_drones~interaction(Temperature,Humidity),
              data=dg) #not homo  p-value < 2.2e-16


bartlett.test(Mean_gyne_cells~interaction(Temperature,Humidity),
              data=dg) #no homo p-value < 2.2e-16

#data is therfore "heteroskadatic", not linnear

#since these are not working, I need to account  account for;

#since colonies will only produce drones later in colony life, i will take into account age
#loading file with age of colony death
queen<-read.csv("queen.csv")

#making a table with the drone information + age
dg_age<-merge(queen,dg)

#--------------------------  drone prodution
dg_age_lmer <- lmer(Mean_drones ~ Temperature*Humidity + Days_since_start+ (1|Room),
                    data = dg_age)
summary(dg_age_lmer)
#Estimate Std. Error t value
#(Intercept)                -1.29493    1.22885  -1.054
#TemperatureYes              0.59059    0.69521   0.850
#HumidityYes                 0.59827    0.83097   0.720
#Days_since_start            0.03408    0.01673   2.037
#TemperatureYes:HumidityYes -0.84826    1.01138  -0.839


car::qqPlot(resid(dg_age_lmer)) # look good
r.squaredGLMM(dg_age_lmer)
#R2m       R2c
#[1,] 0.4824545 0.4824545

Anova(dg_age_lmer, type = 3, test.statistic="F") 

#                     F Df Df.res  Pr(>F)  
#(Intercept)          0.8859  1 6.6900 0.37931  
#Temperature          0.4286  1 6.8049 0.53419  
#Humidity             0.3869  1 6.1452 0.55631  
#Days_since_start     3.9353  1 5.4025 0.09982 .
#Temperature:Humidity 0.7012  1 5.1483 0.43951  

#---------------------  gyne cell production 
gweight_lmer <- lmer(Mean_gyne_cells ~ Temperature * Humidity + Weight_g+
                       (1|Room),
                     data = dg_weights)
summary(gweight_lmer)
#Estimate Std. Error t value
#(Intercept)                -77.27161   30.07540  -2.569
#TemperatureYes               2.58789    3.62601   0.714
#HumidityYes                 -1.10436    2.73543  -0.404
#Weight_g                     0.12801    0.04904   2.610
#TemperatureYes:HumidityYes  -1.39213    5.00052  -0.278

Anova(gweight_lmer, type = 3, test.statistic="F") 
#F Df Df.res  Pr(>F)  
#(Intercept)          5.1462  1 6.2852 0.06182 .
#Temperature          0.4486  1 5.6419 0.52942  
#Humidity             0.1373  1 5.6517 0.72447  
#Weight_g             5.4813  1 6.1118 0.05697 .
#Temperature:Humidity 0.0701  1 5.5370 0.80070  

car::qqPlot(resid(gweight_lmer)) # look good
r.squaredGLMM(gweight_lmer)

#R2m       R2c
#[1,] 0.6603749 0.7253806


#----------------------------- colony growth rate






