#title: "R code for regional variation analyses for Chimpanzee buttress drumming shows rhythmicity and subspecies variation"
#author: Vesta Eleuteri
#date: 31/08/2024

#INSTALL NECESSARY PACKAGES AND LOAD FUNCTIONS----

##Packages collinearity
install.packages("corrplot")
library(corrplot)
install.packages("car")

library(car)

##Packages repDFA and GLMM
install.packages("devtools") 
library(devtools)
devtools::install_github("gobbios/cfp", force=TRUE) 
install_github("gobbios/cfp@v.0.1.0") 
library(cfp)

install.packages("MASS")
library(MASS)

install.packages("lmerTest") 
library(lmerTest)
library(lme4)
library(MuMIn)

library(dplyr)

##Plot theme package
library(scico)

##Functions and packages needed for pDFA, repDFA and GLMM (developed by Mundry and Neumann)
source("pdfa_functions.r")
source("repDFA_nested.r")
source("diagnostic_fcns.r")
source("glmm_stability.r")
source("boot_glmm.r")



#DATA CLEANING----
xdata=read.table("drumming_wide.csv", sep=",", dec=".", header=T, fill=T, stringsAsFactors=TRUE, row.names=NULL)
N_drums_ID=table(xdata$Com_indiv_code) 
N_drums_Com=table(xdata$Community)
N_drum_beats_Com=table(xdata$Com_indiv_code, xdata$N_beats)
N_drum_beats_Sub=table(xdata$Subspecies, xdata$N_beats)
str(xdata)
N_drums_ID_Com=table(xdata$Com_indiv_code, xdata$Community)
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results")
write.table(N_drums_ID_Com, file="Drums_ID_Com_table.csv", row.names=T, col.names=T, sep="\t", dec=".")


#SUBSPECIES, POPULATIONS, ALL COMMUNITIES PDFAS----
##Sort Subspecies, Populations, and Communities Data for pDFAs----
range(xdata$Total_bout_duration)
range(xdata$N_beats)
range(xdata$Bout_cv)
range(xdata$Bout_npvi)
range(xdata$Bout_entropy)

##Add constants to variables to logtransform all non-normal variables with 0 below
c=0.0001
xdata$Total_bout_duration_c=xdata$Total_bout_duration + c
xdata$N_beats_c=xdata$N_beats + c
xdata$Bout_cv_c=xdata$Bout_cv + c
xdata$Bout_npvi_c=xdata$Bout_npvi + c
xdata$Bout_entropy_c=xdata$Bout_entropy + c

range(xdata$Total_bout_duration_c)
range(xdata$N_beats_c)
range(xdata$Bout_cv_c)
range(xdata$Bout_npvi_c)
range(xdata$Bout_entropy_c)

##Distributions and transformations of skewed variables----
hist(xdata$N_beats_c) 
hist(log(xdata$N_beats_c)) 
xdata$log.N_beats_c=log(xdata$N_beats_c) 

hist(xdata$Total_bout_duration_c)
hist(log(xdata$Total_bout_duration_c)) 
xdata$log.Total_bout_duration_c=log(xdata$Total_bout_duration_c) 

hist(xdata$Bout_cv_c) 
hist(log(xdata$Bout_cv_c)) 

hist(xdata$Bout_npvi_c) 
hist(log(xdata$Bout_npvi_c)) 

hist(xdata$Bout_entropy_c) 
hist(log(xdata$Bout_entropy_c)) 

##Communities pDFA----
set.seed(131)
pdfa.res.communities=pDFA.nested(test.fac="Pop_com_code", contr.fac="Com_indiv_code",
                                  variables=c("log.Total_bout_duration_c", "log.N_beats_c", 
                                               "Bout_cv_c", "Bout_npvi_c", "Bout_entropy_c"), restrict.by=NULL,
                                  n.contr.fac.levels.to.sel=NULL,
                                  n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000,
                                  pdfa.data=xdata) 
str(xdata)

pdfa.res.communities$result 
pdfa.res.communities
All_comms_allvabs_pDFA_res=pdfa.res.communities$result
All_comms_allvabs_pDFA_res=round(All_comms_allvabs_pDFA_res, 3)

setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/pDFA")
write.table(All_comms_allvabs_pDFA_res, file="Com_pDFA_res.csv", row.names=T, col.names=T, sep="\t", dec=".")

##Follow up tests -  Communities pDFA----
###All Communities repDFA----
###Extract the features that distinguish communities in pDFA (keep those allowing for discrimination)
set.seed(5)
ImpVars_all_comms_allvabs=repDFA_nested(xdata, testfactor="Pop_com_code", 
                                  balancefactor = c("Pop_com_code", "Com_indiv_code"), 
                                  varnames = c("log.Total_bout_duration_c", "log.N_beats_c", 
                                 "Bout_cv_c", "Bout_npvi_c", "Bout_entropy_c"),
                                  npercomb = 3, nrand = 1000)

ImpVars_all_comms_allvabs
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/repDFA")
All_comms_allvabs_repDFA_res_df1=table(ImpVars_all_comms_allvabs$df1_best) 
All_comms_allvabs_repDFA_res_df1
write.table(All_comms_allvabs_repDFA_res_df1, file="Com_repDFA_df1_res.csv", row.names=T, col.names=T, sep="\t", dec=".")

All_comms_allvabs_repDFA_res_df2=table(ImpVars_all_comms_allvabs$df2_best) 
All_comms_allvabs_repDFA_res_df2
write.table(All_comms_allvabs_repDFA_res_df2, file="Com_repDFA_df2_res.csv", row.names=T, col.names=T, sep="\t", dec=".")

###Communities GLMMs for most contributing repDFA variables: nPVI, N beats, Total bout duration, Entropy----

###Check data
plot(table(table(xdata$Community))) #8-73 observations x community

xx=table(xdata$Individual, xdata$Community)
table(apply(X=xx>0, MARGIN=2, FUN=sum)) #2 to 8 subjects x community

###Check responses
hist(xdata$Bout_npvi_c) #ok a bit skewed
hist(xdata$Bout_entropy_c) #skewed
hist(xdata$log.Total_bout_duration_c) #ok normal

xdata$Community=relevel(xdata$Community, ref="Sonso") 
levels(xdata$Community)

  #Model: nPVI----
Model_npvi=lmer(Bout_npvi_c ~ Community + 
                     (1|Individual), data=xdata, REML=F)

  ##Full-null model comparison
null_Model_npvi=lmer(Bout_npvi_c ~ (1|Individual), data =xdata, REML=F)
Chisq_Model_npvi=as.data.frame(anova(null_Model_npvi, Model_npvi, test="Chisq")) 
Chisq_Model_npvi=round(Chisq_Model_npvi, 3)
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/GLMMs/repDFA variables/npvi")
write.table(Chisq_Model_npvi, file="Com_Chisq_Model_npvi.csv", row.names=T, col.names=T, sep="\t", dec=".")

  #Model: Entropy----
Model_entropy=lmer(Bout_entropy_c ~ Community + 
                  (1|Individual), data =xdata, REML=F)

  ##Full-null model comparison
null_Model_entropy=lmer(Bout_entropy_c ~ (1|Individual), data =xdata, REML=F)
Chisq_Model_entropy=as.data.frame(anova(null_Model_entropy, Model_entropy, test="Chisq")) 
Chisq_Model_entropy=round(Chisq_Model_entropy, 3) #not reported as not loading very high
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/GLMMs/repDFA variables/entropy")
write.table(Chisq_Model_entropy, file="Com_Chisq_Model_entropy.csv", row.names=T, col.names=T, sep="\t", dec=".")

  #Model: Drumming bout duration----
Model_drumduration=lmer(log.Total_bout_duration_c ~ Community + 
                     (1|Individual), data =xdata, REML=F)

  ##Full-null model comparison
null_Model_drumduration=lmer(log.Total_bout_duration_c ~ (1|Individual), data =xdata, REML=F)
Chisq_Model_drumduration=as.data.frame(anova(null_Model_drumduration, Model_drumduration, test="Chisq")) 
Chisq_Model_drumduration=round(Chisq_Model_drumduration, 3)
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/GLMMs/repDFA variables/drum duration")
write.table(Chisq_Model_drumduration, file="Com_Chisq_Model_drum_duration.csv", row.names=T, col.names=T, sep="\t", dec=".")

  #Model: N beats----
Model_numbeats=lmer(log.N_beats_c ~ Community + 
                          (1|Individual), data =xdata, REML=F)

  ##Full-null model comparison
null_Model_numbeats=lmer(log.N_beats_c ~ (1|Individual), data =xdata, REML=F)
Chisq_Model_numbeats=as.data.frame(anova(null_Model_numbeats, Model_numbeats, test="Chisq")) 
Chisq_Model_numbeats=round(Chisq_Model_numbeats, 3)
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/GLMMs/repDFA variables/n beats")
write.table(Chisq_Model_numbeats, file="Com_Chisq_Model_Chisq_Model_numbeats.csv", row.names=T, col.names=T, sep="\t", dec=".")


##Populations pDFA----
set.seed(128)
pdfa.res.populations=pDFA.nested(test.fac="Population", contr.fac="Com_indiv_code",
                                 variables=c("log.Total_bout_duration_c", "log.N_beats_c", 
                                             "Bout_cv_c", "Bout_npvi_c", "Bout_entropy_c"), restrict.by=NULL,
                                 n.contr.fac.levels.to.sel=NULL,
                                 n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000,
                                 pdfa.data=xdata) 


pdfa.res.populations$result 

Pop_allvabs_pDFA_res=pdfa.res.populations$result
Pop_allvabs_pDFA_res=round(Pop_allvabs_pDFA_res, 3)
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/pDFA")
write.table(Pop_allvabs_pDFA_res, file="Pop_pDFA_res.csv", row.names=T, col.names=T, sep="\t", dec=".")


##Subspecies pDFA----
set.seed(126)
pdfa.res.subspecies=pDFA.nested(test.fac="Subspecies", contr.fac="Com_indiv_code",
                                variables=c("log.Total_bout_duration_c", "log.N_beats_c", 
                                            "Bout_cv_c", "Bout_npvi_c", "Bout_entropy_c"), restrict.by=NULL,
                                n.contr.fac.levels.to.sel=NULL,
                                n.to.sel.per.contr.fac.level=NULL, n.sel=100, n.perm=1000,
                                pdfa.data=xdata) 

pdfa.res.subspecies$result 
pdfa.res.subspecies

Sub_allvabs_pDFA_res=pdfa.res.subspecies$result
Sub_allvabs_pDFA_res=round(Sub_allvabs_pDFA_res, 3)

setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/pDFA")
write.table(Sub_allvabs_pDFA_res, file="Sub_pDFA_res.csv", row.names=T, col.names=T, sep="\t", dec=".")

##Follow up tests - Subspecies----
###Subspecies repDFA----
###Extract the features that distinguish subspecies in pDFA (keep those allowing for discrimination)
set.seed(5)
ImpVars_sub_allvabs=repDFA_nested(xdata, testfactor="Subspecies", 
                                  balancefactor = c("Subspecies", "Com_indiv_code"), 
                                  varnames = c("log.Total_bout_duration_c", "log.N_beats_c", 
                                               "Bout_cv_c", "Bout_npvi_c", "Bout_entropy_c"),
                                  npercomb = 3, nrand = 1000)

ImpVars_sub_allvabs
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/repDFA")
Sub_allvabs_repDFA_res_df1=table(ImpVars_sub_allvabs$df1_best) #Bout_npvi_c=901 (90% of 100%) 
Sub_allvabs_repDFA_res_df1
write.table(Sub_allvabs_repDFA_res_df1, file="Sub_repDFA_df1_res.csv", row.names=T, col.names=T, sep="\t", dec=".") 

###Subspecies GLMMs for most contributing repDFA variables----
  #Model: nPVI----
Model_npvi_sub=lmer(Bout_npvi_c ~ Subspecies + 
                      (1|Individual), data=xdata, REML=F)

  ##Full-null model comparison
null_Model_npvi_sub=lmer(Bout_npvi_c ~ (1|Individual), data =xdata, REML=F)
Chisq_Model_npvi_sub=as.data.frame(anova(null_Model_npvi_sub, Model_npvi_sub, test="Chisq"))
Chisq_Model_npvi_sub=round(Chisq_Model_npvi_sub, 3)
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/GLMMs/repDFA variables/npvi")
write.table(Chisq_Model_npvi_sub, file="Sub_Chisq_Model_npvi.csv", row.names=T, col.names=T, sep="\t", dec=".")

  #Model: Number of beats----
Model_nbeats_sub=lmer(log.N_beats_c ~ Subspecies + 
                        (1|Individual), data=xdata, REML=F)

  ##Full-null model comparison
null_Model_nbeats_sub=lmer(log.N_beats_c ~ (1|Individual), data =xdata, REML=F)
Chisq_Model_nbeats_sub=as.data.frame(anova(null_Model_nbeats_sub, Model_nbeats_sub, test="Chisq")) #sig
Chisq_Model_nbeats_sub=round(Chisq_Model_nbeats_sub, 3) 
setwd("~/Desktop/Drumming/Acoustic data analyses 2023/Analyses/results/GLMMs/repDFA variables/n beats")
write.table(Chisq_Model_npvi_sub, file="Sub_Chisq_Model_npvi.csv", row.names=T, col.names=T, sep="\t", dec=".")


##PANT-HOOT DRUMMING COMBINATION GLMMs----
##Sort data
str(xdata)
table(xdata$Start_Ph, xdata$Start_Ph_N)
xdata_ph=subset(xdata, Start_Ph %in% c("Build-Up", "Climax", "Let-Down", "roar", "Roar")) 
xdata_ph=droplevels(xdata_ph) 
table(xdata_ph$Start_Ph, xdata_ph$Start_Ph_N)

levels(xdata_ph$Start_Ph)[levels(xdata_ph$Start_Ph) %in% c("Build-Up", "roar", "Roar")]="Before Climax Start" 
levels(xdata_ph$Start_Ph)[levels(xdata_ph$Start_Ph) %in% c("Climax", "Let-Down")]="After Climax Start" 
str(xdata_ph)
N_drums_Com_ph=table(xdata_ph$Community) #ok all with 3 drums 

##Check frequencies of response variable
table(xdata_ph$Start_Ph) 
table(xdata_ph$Start_Ph, xdata_ph$Community) #South has 0 drums After Climax Start so remove as no sufficient variation
table(xdata_ph$Start_Ph, xdata_ph$Population) 
xdata_ph=subset(xdata_ph, Community != c("S")) 
xdata_ph=droplevels(xdata_ph) 

##Dummy code response variable
xdata_ph$Start_Ph=as.numeric(xdata_ph$Start_Ph=="Before Climax Start") 
str(xdata_ph)

##Sorting Community and Population baselines 
xdata_ph$Community=relevel(xdata_ph$Community, ref="Sonso") 
xdata_ph$Population=relevel(xdata_ph$Population, ref="Budongo") 

  #Model: Community----
Model_start_ph_com=glmer(Start_Ph ~ Community + 
                       (1|Individual), data=xdata_ph, family=binomial, control=contr) 

  ##Full-null model comparison
null_Model_start_ph_com=glmer(Start_Ph ~
                            (1|Individual), data=xdata_ph, family=binomial, control=contr)
Chisq_Model_start_ph_com=as.data.frame(anova(null_Model_start_ph_com, Model_start_ph_com, test="Chisq")) 
Chisq_Model_start_ph_com=round(Chisq_Model_start_ph_com, 3)
write.table(Chisq_Model_start_ph_com, file="Com_Chisq_Model_Start_ph.csv", row.names=T, col.names=T, sep="\t", dec=".")

  #Model: Population----
Model_start_ph_pop=glmer(Start_Ph ~ Population +
                       (1|Individual), data=xdata_ph, family=binomial, control=contr) 

  ##Full-null model comparison
null_Model_start_ph_pop=glmer(Start_Ph ~
                            (1|Individual), data=xdata_ph, family=binomial, control=contr)
Chisq_Model_start_ph_pop=as.data.frame(anova(null_Model_start_ph_pop, Model_start_ph_pop, test="Chisq")) 
Chisq_Model_start_ph_pop=round(Chisq_Model_start_ph_pop, 3)
write.table(Chisq_Model_start_ph_pop, file="Pop_Chisq_Model_Start_ph.csv", row.names=T, col.names=T, sep="\t", dec=".")

  #Model: Subspecies----
Model_start_ph_sub=glmer(Start_Ph ~ Subspecies +
                           (1|Individual), data=xdata_ph, family=binomial, control=contr) 

  ##Full-null model comparison
null_Model_start_ph_sub=glmer(Start_Ph ~
                                (1|Individual), data=xdata_ph, family=binomial, control=contr)
Chisq_Model_start_ph_sub=as.data.frame(anova(null_Model_start_ph_sub, Model_start_ph_sub, test="Chisq")) 
Chisq_Model_start_ph_sub=round(Chisq_Model_start_ph_sub, 3) 
write.table(Chisq_Model_start_ph_sub, file="Sub_Chisq_Model_Start_ph.csv", row.names=T, col.names=T, sep="\t", dec=".")
