#Stable Isotope Analysis July 2022
#Emilie Diesbourg
#2023-03-14

#### Loading packages/data ####
library(rjags)
library(MixSIAR)
library(dplyr)
library(plyr)
library(ggplot2)
library(R2WinBUGS)
library(splancs)
library(devtools)
library(multcomp)
library(vegan)
library(tidyverse)

isotope<-read.csv("C:\\Users\\Emilie\\Documents\\McMaster University 2022-2024\\Bow River Thesis\\BowRiverRProject\\Stable_isotope_ED\\StableIsotopeAvg.csv")
isotope$Site<-factor(isotope$Site, levels=c("Cochrane", "Sunalta", "Cushing Bridge", "Graves Bridge", "Policeman Flats"))
str(isotope)

#### Preliminary look at data ####

#Checking for indication of carbonates

plant<-isotope[isotope$Taxa %in% c("Biofilm", "Shrub"), ]

biofilm<-filter(isotope, Taxa=="Biofilm")

ggplot(biofilm, aes(y=X.C, x=d13C))+ 
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(.~Site)+
  xlab(expression(paste(delta^13, "C (\u2030)", sep="")))+
  ylab("%C")+
  theme_classic()

biofilmhcl<-subset(biofilm, Well>"D6")
ggplot(biofilmhcl, aes(y=X.C, x=d13C))+ 
  geom_point()+
  geom_smooth(method="lm")+
  xlab(expression(paste(delta^13, "C (\u2030)", sep="")))+
  ylab("%C")+
  theme_classic()

biofilmuntreated<-subset(biofilm, Well < "D6")
ggplot(biofilmuntreated, aes(y=X.C, x=d13C))+ 
  geom_point()+
  geom_smooth(method="lm")+
  xlab(expression(paste(delta^13, "C (\u2030)", sep="")))+
  ylab("%C")+
  theme_classic()

#Checking for indication of effect of lipid
inverts<-filter(isotope, Tissue.Type=="Whole Invertebrate")

ggplot(inverts, aes(y=C.N, x=d13C))+ #Only looks slightly depleted at Policeman Flats so probably not an issue for this data
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(.~Site)+
  xlab(expression(paste(delta^13, "C (\u2030)", sep="")))+
  ylab("Carbon : Nitrogen")+
  theme_classic()

ggplot(inverts, aes(x=Taxa, y=C.N))+
  geom_boxplot()+
  geom_hline(aes(yintercept=3.5), colour="red")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))+
  ylab("Carbon : Nitrogen")


#Apparently the C:N ration shouldn't be above 3.5 but almost all samples are above that threshold...not sure what that means
#Brian says that could indicate a lipid effect

#Check for nitrogen enrichment across sites (all inverts included)
ggplot(inverts, aes(x=Site, y=d15N))+
  geom_boxplot()+
  theme_classic()
#Looks like nitrogen is enriched at Graves Bridge which is interesting because that site was 
#right below the wastewater outfall from Bonnybrook. This could indicate a shift in food web structure due to 
#municipal wastewater exposure. However, Policeman flats dips down again. Potentially from dissipation of 
#wastewater?

#Statistically test this
hist(inverts$d15N) #looks pretty normal
shapiro.test(inverts$d15N) #p=0.1138, normal distribution. We can use ANOVA. 

summary(a1<-aov(data=inverts, d15N~Site)) #p<0.001. Differences btw sites. 
summary(glht(a1, linfct=mcp(Site="Tukey")))
#Differences between:
#Graves Bridge and Cochrane
#Graves Bridge and Sunalta
#Graves Bridge and Cushing Bridge
#Graves Bridge and Policeman Flats
#Almost between Policeman Flats and Cushing Bridge but not quite (p=0.055)

plot(a1) #residuals look good. 

#Nitrogen enrichment in larval chloroperlidae
larvae<-filter(inverts, Stage=="Larvae")
larchloro<-filter(larvae, Taxa=="Chloroperlidae")

ggplot(larchloro, aes(x=Site, y=d15N))+
  geom_boxplot()+
  theme_classic()
#Graves Bridge is much higher in d15N than the other two sites (Sunalta + Cushing Bridge). Indicates an 
#influence of municipal wasteater. 

larhep<-filter(larvae, Taxa=="Heptageniidae")
ggplot(larhep, aes(x=Site, y=d15N))+
  geom_boxplot()+
  theme_classic()
#Graves Bridge is much higher in d15N but Policeman Flats is the lowest.

araneidae<-filter(inverts, Taxa=="Araneidae")
ggplot(araneidae, aes(x=Site, y=d15N))+
  geom_boxplot()+
  theme_classic()
#Much higher in d15N at Graves Bridge. 

tetragnathidae<-filter(inverts, Taxa=="Tetragnathidae")
ggplot(tetragnathidae, aes(x=Site, y=d15N))+
  geom_boxplot()+
  theme_classic()
#Much higher in d15N at Graves Bridge.

#### Biplots of individual data ####

#d13C vs d15N at all sites
ggplot(isotope, aes(x=d13C, y=d15N, color=Taxa, shape=Tissue.Type))+
  geom_point()+
  facet_grid(.~Site)+
  theme_bw()

#Sites separated out
Cochrane<-filter(isotope, Site=="Cochrane")
ggplot(Cochrane, aes(x=d13C, y=d15N, color=Taxa))+
  geom_point()

Sunalta<-filter(isotope, Site=="Sunalta")
ggplot(Sunalta, aes(x=d13C, y=d15N, color=Taxa, shape=Tissue.Type))+
  geom_point()


Cushbr<-filter(isotope, Site=="Cushing Bridge")
ggplot(Cushbr, aes(x=d13C, y=d15N, color=Taxa))+
  geom_point()

Grvbr<-filter(isotope, Site=="Graves Bridge")
ggplot(Grvbr, aes(x=d13C, y=d15N, color=Taxa))+
  geom_point()

Polflt<-filter(isotope, Site=="Policeman Flats")
ggplot(Polflt, aes(x=d13C, y=d15N, color=Taxa))+
  geom_point()

#### Averaged data + SDs biplots ####

#Sumamrizing a dataframe with mean and sd C and N values
data.sum<-ddply(isotope, c("Site", "Taxa", "FFG"), summarise,
                d13Cmn=mean(d13C), #mean
                d13Csd=sd(d13C), #standard deviation
                d13Cse=sd(d13C)/sqrt(length(d13C)), #standard error
                d15Nmn=mean(d15N), #mean
                d15Nsd=sd(d15N), #standard deviation
                d15Nse=sd(d15N)/sqrt(length(d15N))) #standard error

data.sum1<-ddply(isotope, c("Taxa", "FFG"), summarise,
                d13Cmn=mean(d13C), #mean
                d13Csd=sd(d13C), #standard deviation
                d13Cse=sd(d13C)/sqrt(length(d13C)), #standard error
                d15Nmn=mean(d15N), #mean
                d15Nsd=sd(d15N), #standard deviation
                d15Nse=sd(d15N)/sqrt(length(d15N))) #standard error

#Making errorbars
Ylims<-aes(ymax=d15Nmn +d15Nsd, ymin=d15Nmn-d15Nsd)
Xlims<-aes(xmax=d13Cmn+d13Csd, xmin=d13Cmn-d13Csd)

#All sites together

ggplot(data.sum1, aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(paste(delta^15, "N (\u2030)", sep="")))+
  xlab(expression(paste(delta^13, "C (\u2030)", sep="")))+
  theme_classic()

#All sites (faceted)

ggplot(data.sum, aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(delta^15~N))+
  xlab(expression(delta^13~C))+
  facet_grid(.~Site)+
  theme_bw()

ggplot(subset(data.sum, Site=="Cochrane"), aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(delta^15~N))+
  xlab(expression(delta^13~C))+
  theme_classic()

ggplot(subset(data.sum, Site=="Sunalta"), aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(delta^15~N))+
  xlab(expression(delta^13~C))+
  theme_classic()

ggplot(subset(data.sum, Site=="Cushing Bridge"), aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(delta^15~N))+
  xlab(expression(delta^13~C))+
  theme_classic()

ggplot(subset(data.sum, Site=="Graves Bridge"), aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(delta^15~N))+
  xlab(expression(delta^13~C))+
  theme_classic()

ggplot(subset(data.sum, Site=="Policeman Flats"), aes(x=d13Cmn, y=d15Nmn, colour=Taxa, shape=FFG))+
  geom_point(size=3)+
  geom_errorbar(Ylims, width=0.2)+
  geom_errorbarh(Xlims, height=0.2)+
  ylab(expression(delta^15~N))+
  xlab(expression(delta^13~C))+
  theme_classic()
          

#### MixSIAR models ####

#Model with only larvae consumers including stonefly larvae (Perlidae and Chloroperlidae)
mixlar<-load_mix_data(filename="Stable_isotope_ED\\larvae_consumer.csv", 
                   iso_names=c("d13C","d15N"), 
                   factors=c("Site","Taxa"), 
                   fac_random=c(TRUE,TRUE), 
                   fac_nested = c(FALSE, FALSE),
                   cont_effects=NULL)

sourcelar<-load_source_data(filename="Stable_isotope_ED\\larvae_sources.csv",
                         source_factors = "Site",
                         conc_dep = FALSE, 
                         data_type = "means",
                         mixlar)

discrlar<-load_discr_data(filename="Stable_isotope_ED\\isotope_tef.csv", mixlar)

#isospace plot
plot_data(filename="isospace_plot", plot_save_pdf = F, plot_save_png = F, mixlar, sourcelar, discrlar)

#default uninformative/generalist prior (alpha=1)
plot_prior(alpha.prior = c(1), source)

#Informative prior
plot_prior(alpha.prior = c(0.5, 1), source)

#How do you know which weighting to use? Seems like the inverts consume a more terrestrial source so should I weight it higher?

#JAGS model file
model_filename<-"MixSIAR_model1.txt"
resid_err<-TRUE
process_err<-TRUE
write_JAGS_model(model_filename, resid_err, process_err, mixlar, sourcelar)

#Test model to see how long we need to run the real model
#jags.1<-run_model(run="test", mixlar, sourcelar, discrlar, model_filename, alpha.prior = 1, resid_err, process_err)

output_options<-list(summary_save=TRUE,
                     summary_name="summary_statistics",
                     sup_post=TRUE,
                     plot_post_save_pdf=TRUE,
                     plot_post_name="posterior_density",
                     sup_pairs=TRUE,
                     plot_pairs_save_pdf=TRUE,
                     plot_pairs_name="pairs_plot",
                     sup_xy=TRUE,
                     plot_xy_save_pdf=TRUE,
                     plot_xy_name="xy_plot",
                     gelman=TRUE, 
                     heidel=TRUE,
                     geweke=TRUE,
                     diag_save=TRUE,
                     diag_name="diagnostics",
                     indiv_effect=FALSE,
                     plot_post_save_png=FALSE, 
                     plot_pairs_save_png=FALSE,
                     plot_xy_save_png=FALSE)

#output_JAGS(jags.1, mixlar, sourcelar, output_options)
#Important to look at diagnostics and summary statistics (gives mean resource use (%) and confidence intervals)

#Running real model with 'normal' chain length and burn-in.
run<-list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
jags<-run_model(run, mixlar, sourcelar, discrlar, model_filename, alpha.prior = 1, resid_err, process_err)
output_JAGS(jags, mixlar, sourcelar, output_options)

#Model with larvae consumers but without stonefly larvae (because they are predators)
mixlarnostonefly<-load_mix_data(filename="Stable_isotope_ED\\larvae_consumer_no_stonefly.csv", 
                      iso_names=c("d13C","d15N"), 
                      factors=c("Site","Taxa"), 
                      fac_random=c(TRUE,TRUE), 
                      fac_nested = c(FALSE, FALSE),
                      cont_effects=NULL)
run<-list(chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE)
jags_lar_nostonefly<-run_model(run, mixlarnostonefly, sourcelar, discrlar, model_filename, alpha.prior = 1, resid_err, process_err)
output_JAGS(jags_lar_nostonefly, mixlarnostonefly, sourcelar, output_options1)

#### Adonis calculation between invertebrate taxa ####


