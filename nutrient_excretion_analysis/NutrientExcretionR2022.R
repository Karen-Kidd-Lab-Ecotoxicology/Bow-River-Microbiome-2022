#Title: Nutrient Excretion Analyses 2022
#Author: Emilie Diesbourg
#DAte: 2023-03-20

#rm(list=ls()) #clears R environment

#### Loading packages and data #### 
#{.tabset .tabset-dropdown}
#setwd("Nutrient_excretion_ED")
library(dplyr)
library(FSA)
library(ggplot2)
library(tidyverse)
library(car)
library(lme4)
library(lmerTest)
library(multcomp)
library(AICcmodavg)
library(ggbreak)
library(emmeans)
library(ggsignif)

Nut_data<-read.csv("C:\\Users\\Emilie\\Documents\\McMaster University 2022-2024\\Bow River Thesis\\BowRiverRProject\\Nutrient_excretion_ED\\NutrientExcretionR2022.csv")
str(Nut_data)
Nut_data$Site <- factor(Nut_data$Site, levels = c("Cochrane", "Sunalta", "Graves Bridge", "Policeman Flats", "BBR2", "PCR1", "PCR3"))
Nut_data$Type<-factor(Nut_data$Type, levels = c("Bow River", "ACWA Stream"))
Nut_data$Length<-as.numeric(Nut_data$Length)
Nut_data$Mass<-as.numeric(Nut_data$Mass)
Nut_data$M.L<-as.numeric(Nut_data$M.L)
Nut_data$P_Mass<-as.numeric(Nut_data$P_Mass)
Nut_data$NH3_Mass<-as.numeric(Nut_data$NH3_Mass)
Nut_data$P_Norm<-as.numeric(Nut_data$P_Norm)
Nut_data$NH3_Norm<-as.numeric(Nut_data$NH3_Norm)
Nut_data$Excr_rate_P<-as.numeric(Nut_data$Excr_rate_P)
Nut_data$Excr_rate_NH3<-as.numeric(Nut_data$Excr_rate_NH3)
#Nut_data$Water_Temp<-as.factor(Nut_data$Water_Temp)
Nut_data$Dissolved_oxygen<-as.numeric(Nut_data$Dissolved_oxygen)
Nut_data$TP<-as.numeric(Nut_data$TP)

str(Nut_data)

Nut_samples<-filter(Nut_data, Sample_Type=="Sample")


July<-filter(Nut_samples, Month=="July")
September<-filter(Nut_samples, Month=="September")

Cochrane<-filter(Nut_samples, Site=="Cochrane")
Cochjuly<-filter(July, Site=="Cochrane")
Cochsept<-filter(September, Site=="Cochrane")

Sunalta<-filter(Nut_samples, Site=="Sunalta")
Sunjuly<-filter(July, Site=="Sunalta")
Sunsept<-filter(September, Site=="Sunalta")

pmf<-filter(Nut_samples, Site=="Policeman Flats")
pmfjuly<-filter(July, Site=="Policeman Flats")
pmfsept<-filter(September, Site=="Policeman Flats")

#### Phosphorus Data ####

#Total phosphorus concentration (- blank phosphorus concentration) divided by the total mass 
#of organisms in each tube at each site, divided by the time the organisms were in each tube. 
#(Rate of excretion) 

ggplot(Nut_samples, aes(x=Site, y=Excr_rate_P, colour=Month))+ 
  geom_boxplot()+
  theme_classic()+
  ylab("Phosphorus (ug/L/min)")+
  labs(colour="Sampling Month")+
  scale_y_log10()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  facet_grid(~Type, scales = "free")
#Looks like not much difference across sampling months but a difference at Policeman Flats compared to other sites

ggplot(Nut_samples, aes(x=Mass, y=org_rate_P, colour=Site))+ 
  geom_point()+
  theme_classic()+
  ylab("Phosphorus (ug/L/min)")+
  labs(colour="Sampling Month")+
  scale_y_log10()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  geom_smooth(method="glm", se=F)+
  facet_grid(~Type, scales="free")

#Testing whether excretion rate is different across sampling time points
shapiro.test(Cochjuly$Excr_rate_P) #Normal
shapiro.test(Cochsept$Excr_rate_P) #Not normal

#Have to use non-parametric test to test differences between excretion rates
wilcox.test(Excr_rate_P ~ Date, data=Cochrane) #p=0.7394. Not sig. different.

shapiro.test(Sunjuly$Excr_rate_P) #Normal
shapiro.test(Sunsept$Excr_rate_P) #Not normal

#Have to use non-parametric test
wilcox.test(Excr_rate_P ~ Date, data=Sunalta) #p=0.8392. Not sig. different.

shapiro.test(pmfjuly$Excr_rate_P) #Normal
shapiro.test(pmfsept$Excr_rate_P) #Normal

#Can use a t-test. 
wilcox.test(Excr_rate_P ~ Month, data=pmf) #p=0.9654. Not sig. different.
t.test(Excr_rate_P~Month, data=pmf) #p=0.7646 Not sig. different. 

hist(Nut_samples$Excr_rate_P) #Data are heavily right skewed
hist(log(Nut_samples$Excr_rate_P)) #Looks much better when log transformed

logPExcr<-log(Nut_samples$Excr_rate_P)
cbind(logPExcr, Nut_samples)

#Looking for outliers

ggplot(Nut_samples, aes(x=Mass, y=Phosphorus))+
  geom_point()+
  geom_smooth(method="glm", se=T)+
  geom_point(data=Nut_samples %>% filter(Phosphorus>600), pch=21,size=10, colour="purple")

#2 major outliers: COCH_NUT_H2O_4, SUN_NUT_H2O_8


#### Mixed effects model phosphorus####

hist(Nut_samples$org_rate_P) #Right skewed
hist(log(Nut_samples$org_rate_P)) #Looks normal when log transformed

logPrate<-log(Nut_samples$org_rate_P)
Nut_samples<-cbind(logPrate, Nut_samples)


summary(m1<-lmer(data=Nut_samples, logPrate~Mass * Site * Water_Temp + (1|Month)))
summary(glht(m1, linfct=mcp(Site="Tukey")))
plot(m1)

summary(m2<-lmer(data=Nut_samples, logPrate~Mass + Site + Water_Temp + (1|Month)))
summary(glht(m2, linfct=mcp(Site="Tukey")))
plot(m2)

summary(m3<-lmer(data=Nut_samples, logPrate~Mass + Site + (1|Month))) 
#This one is the best using AIC and p-value method
summary(glht(m3, linfct=mcp(Site="Tukey")))
plot(m3)

summary(m4<-lmer(data=Nut_samples, logPrate~ Site + (1|Month)))
summary(glht(m4, linfct=mcp(Site="Tukey")))
plot(m4)

mixedmodels=list(m1, m2, m3, m4) 
aictab(cand.set = mixedmodels) 
#Model 3 without water temperature and the interaction term between Mass and Site is the best predictor of nutrient excretion

#Using multiple linear regression instead with Month as a fixed term.

summary(lm1<-aov(data=Nut_samples, logPrate~Mass*Site*Month*Water_Temp)) #No significant interaction terms so I will drop them
summary(lm2<-aov(data=Nut_samples, logPrate~Mass+Site+Month+Water_Temp)) #Drop Water_Temp because it's the next most insignificant
summary(lm3<-aov(data=Nut_samples, logPrate~Mass+Site+Month)) #Drop Month. 
summary(lm4<-aov(data=Nut_samples, logPrate~Site+Mass)) 
summary(lm5<-aov(data=Nut_samples, logPrate~Site))

mixedmodels=list(lm1, lm2, lm3, lm4, lm5) 
aictab(cand.set = mixedmodels) 
#Model 4 was the best with Mass and Site being significant predictors of PO4-P excretion
plot(lm4) #Meets the assumptions of linear model. One outlier but not terrible
summary(glht(lm4, linfct=mcp(Site="Tukey")))

#What makes more sense? Linear mixed model with Month as a random factor or a general linear model with Month included in the model?
#The results are the same nonetheless. 

#Linear model graphs
#Phosphorus excretion rate by site and month
ggplot(Nut_samples, aes(x=Mass, y=org_rate_P, colour=Site))+
  geom_point()+
  geom_smooth(method="glm", se=F)+
  theme_bw()+
  scale_y_log10()+
  facet_grid(.~Month)+
  theme(legend.position="bottom")+
  ylab(expression(paste("Excretion Rate" ~PO[4]~"-P ("~mu~"g/L/minute)")))
  
#ggsave(filename =  "Dry mass predicts PO4-P excretion rate july vs september.png", height=6, width=7)

#PO4-P excretion rate by site and stream type
ggplot(Nut_samples, aes(x=Mass, y=org_rate_P, colour=Site))+
  geom_point()+
  geom_smooth(method="lm", se=F)+
  theme_bw()+
  scale_y_log10()+
  facet_grid(.~Type)+
  theme(legend.position="bottom")+
  ylab(expression(paste("Excretion Rate" ~PO[4]~"-P ("~mu~"g/L/minute)")))

#ggsave(filename =  "Dry mass predicts PO4-P excretion rate river vs stream.png", height=6, width=7)

#### NH3 Data ####
#Ammonia

ggplot(Nut_samples, aes(x=Site, y=org_rate_NH3, colour=Month))+ 
  geom_boxplot()+
  theme_classic()+
  ylab(expression(paste("Excretion Rate" ~NH[3]~"-N ("~mu~"g/L/minute)")))+
  labs(colour="Sampling Month")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
#Looks like there are differences in nitrogen excretion across sampling months 

#Log transforming y axis, will remove all negative values (16 of them). 
ggplot(Nut_samples, aes(x=Site, y=org_rate_NH3, colour=Month))+ 
  geom_boxplot()+
  theme_classic()+
  ylab("Nitrogen Excretion Rate (ug/L/mg caddis/min)")+
  labs(colour="Sampling Month")+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  scale_y_log10()


hist(Nut_samples$org_rate_NH3)
hist(log(Nut_samples$org_rate_NH3))
logNH3Excr<-log(Nut_samples$org_rate_NH3)
Nut_samples<-cbind(logNH3Excr, Nut_samples)


summary(lm1<-lm(data=Nut_samples, Excr_rate_P~Excr_rate_NH3)) 
cor(Nut_samples$Excr_rate_P, Nut_samples$Excr_rate_NH3) #0.018. Not a high correlation between these variables
#Phosphorus excretion does not predict ammonia excretion
plot(lm1) #One very large outlier

#Testing to see if nitrogen excretion rate is sig different across sampling points
shapiro.test(Cochjuly$logNH3Excr) #Normal
shapiro.test(Cochsept$logNH3Excr) #Not normal

#Have to use non-parametric test
wilcox.test(logNH3Excr ~ Month, data=Cochrane) #p=0.015. Different across sampling points so they cannot be combined.

shapiro.test(Sunjuly$logNH3Excr) #Normal
shapiro.test(Sunsept$logNH3Excr) #Normal

#Can use T-test
t.test(logNH3Excr~Month, data=Sunalta) #p=0.1767. Not sig. different


shapiro.test(pmfjuly$logNH3Excr) #Not normal
shapiro.test(pmfsept$logNH3Excr) #Not normal

#Have to use non-parametric test
wilcox.test(logNH3Excr ~ Month, data=pmf) #p=0.2667. Sig. different.


#### Mixed effects model NH3 ####

summary(NH3_m1<-lmer(data=Nut_samples, logNH3Excr~Site * Mass * Water_Temp + (1|Month)))
summary(glht(NH3_m1, linfct=mcp(Site="Tukey")))
plot(NH3_m1)

summary(NH3_m2<-lmer(data=Nut_samples, logNH3Excr~Site + Mass + Water_Temp+(1|Month)))
summary(glht(NH3_m2, linfct=mcp(Site="Tukey")))
plot(NH3_m2)

summary(NH3_m3<-lmer(data=Nut_samples, logNH3Excr~Water_Temp+Site+(1|Month)))
summary(glht(NH3_m3, linfct=mcp(Site="Tukey")))
plot(NH3_m3)

summary(NH3_m4<-lmer(data=Nut_samples, logNH3Excr~Water_Temp+(1|Month)))
summary(glht(NH3_m4, linfct=mcp(Site="Tukey")))
plot(NH3_m4)

summary(NH3_m5<-lmer(data=Nut_samples, logNH3Excr~Site+(1|Month)))
summary(glht(NH3_m5, linfct=mcp(Site="Tukey")))
plot(NH3_m5)
#Best model based on AIC

NH3mixedmodels=list(NH3_m1, NH3_m2, NH3_m3, NH3_m4, NH3_m5) 
aictab(cand.set = NH3mixedmodels)

#Linear regression

summary(NH3lm1<-lm(data=Nut_samples, logNH3Excr~Site * Mass * Water_Temp * Month)) #No sig interactions
summary(NH3lm2<-lm(data=Nut_samples, logNH3Excr~Site + Mass + Water_Temp + Month)) 
#Best model using AIC
summary(NH3lm3<-lm(data=Nut_samples, logNH3Excr~Site + Water_Temp + Month))
summary(NH3lm4<-lm(data=Nut_samples, logNH3Excr~Site + Water_Temp))

linearmodels<-list(NH3lm1, NH3lm2, NH3lm3, NH3lm4)
aictab(cand.set = linearmodels)

summary(glht(NH3lm2, linfct=mcp(Site="Tukey")))
plot(NH3lm2) #Meets assumptions of linear model

#Nitrogen excretion rate at each site, individual sampling month. 
ggplot(Nut_samples, aes(x=Mass, y=org_rate_NH3, colour=Site))+ 
  geom_point()+
  geom_smooth(method="glm", se=F)+
  scale_y_log10()+
  facet_grid(.~Month)+
  theme_bw()+
  ylab(expression(paste("Excretion Rate" ~NH[3]~"-N ("~mu~"g/L/minute)")))+
  labs(colour="Sampling Month")+
  scale_y_log10()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))

#Combined by sampling month
ggplot(Nut_samples, aes(x=Mass, y=org_rate_NH3, colour=Site))+ 
  geom_point()+
  geom_smooth(method="glm", se=F)+
  scale_y_log10()+
  facet_grid(.~Type)+
  theme_bw()+
  ylab(expression(paste("Excretion Rate" ~NH[3]~"-N ("~mu~"g/L/minute)")))+
  labs(colour="Sampling Month")+
  scale_y_log10()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))


#Looking for outliers

ggplot(Nut_samples, aes(x=Mass, y=Ammonia))+
  geom_point()+
  geom_smooth(method="glm")+
 geom_point(data=Nut_samples %>% filter(Ammonia==444), pch=21,size=10, colour="purple")+
  geom_point(data=Nut_samples %>% filter(Ammonia>500), pch=21,size=10, colour="purple")

#2 major outliers are COCH_NUT_H2O_4 and COCH_NUT_H2O_9

#### 2 panel plot with P and NH3 excretion rate, normalized by mass ####

#PO4-P excretion rate normalized by mass of organism (both time points combined)

excr_P_plot<-ggplot(Nut_samples, aes(x=Site, y=Excr_rate_P))+ 
  geom_boxplot()+
  theme_classic()+
  ylab(expression(paste("Excretion Rate" ~PO[4]~"-P ("~mu~"g/L/mg caddis/minute)")))+
  labs(colour="Sampling Month")+
  scale_y_log10()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(size =16, angle = 40, vjust = 1, hjust=1),
        axis.text.y = element_text(size =16))+
  facet_grid(~Type, scales='free')
#geom_signif(comparisons = list(c("Cochrane", "Sunalta", "Graves Bridge", "Policeman Flats", "BRR2", "PCR1", "PCR3")), 
#map_signif_level=TRUE)

excr_P_plot2<-excr_P_plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#NH3-N excretion rate normalized by mass of organism (both time points combined)
excr_NH3_plot<-ggplot(Nut_samples, aes(x=Site, y=Excr_rate_NH3))+ 
  geom_boxplot()+
  scale_y_log10()+
  theme_classic()+
  ylab(expression(paste("Excretion Rate" ~NH[3]~"-N ("~mu~"g/L/mg caddis/minute)")))+
  labs(colour="Sampling Month")+
  scale_y_log10()+
  theme(legend.position="bottom")+
  theme(axis.text.x = element_text(size=16, angle = 40, vjust = 1, hjust=1),
        axis.text.y = element_text(size=16))+
  facet_grid(~Type, scales="free")

excr_NH3_plot

excr_NH3_plot2<-excr_NH3_plot + theme(axis.title.y = element_blank(), axis.title.x = element_blank())

excr_rate_nutrient<-ggarrange(excr_P_plot2, excr_NH3_plot2, ncol=2, nrow=1, align = "v")

annotate_figure(excr_rate_nutrient, bottom = text_grob("Site", size = 18))

