######## PROJECT: Laboratory experiment: Germination and seedling performance after freezing 
#### PURPOSE:Examine germination success and seedling survival and growth in response to freezing treatments .
#### AUTHOR: Inam Jameel
# AUTHOR: Inam Jameel
#### DATE LAST MODIFIED: 2 Nov 23

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4) #for running linear mixed models
require(ggplot2) #for plotting 
require(visreg) # for plotting
require(car) # to run ANOVA on model output
require(plyr) # for data wrangling
require(dplyr) # for data wrangling
require(effects) # for plotting
require(emmeans) #for plotting
require(glmmTMB) # for running survival model, have to load twice
require(gamlss) # for running phenology model
require(broom.mixed) #for making tables
require(multcomp) #for pairwise comparisons
require(vioplot) #for violin plots


# set working directory
setwd("/Users/inam/Library/CloudStorage/OneDrive-UniversityofGeorgia/FrztolOnedrive/Freeze tolerance/data") #this is where you specify the folder where you have the data on your computer


#*******************************************************************************
#### 1. Field temperature data and chamber schedule for laboratory experiment #####
#*******************************************************************************

#to plot the air temperature data in the window of snow melt for 2021 for garden located at 2553m

tempdata <-read.csv("Estess20202021_air.csv", header=T)

#filter for temperature during window of snow melt
temp_window <- tempdata %>% filter(order %in% (8477:9548) )

#Figure 2a  line graph of temperature during window of snow melt 
ggplot(data = temp_window) + theme_bw() +
  geom_line(aes(x = order, y = temp)) + scale_x_continuous("Date (2021)",breaks=seq(8477,9548,48),labels=c("Mar 28","29","30","31","Apr 1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))+scale_y_continuous("Temperature (C°)")


#to plot freezing chamber temperature schedule for laboratory experiment

chamber_temp <-read.csv("Laboratory_treatment_schedule.csv", header=T)

chamber_temp$Tmin <- as.factor(chamber_temp$Tmin)


#Figure 2b line graph of treatment levels for laboratory experiment 
ggplot(data = chamber_temp) +
  geom_line(aes(x = hour, y = Temp, group = Tmin, color = Tmin, linetype = Tmin),linewidth=0.9) + theme_bw()+scale_linetype_manual(values=c("solid", "dashed","dotted","dotdash")) + scale_color_manual(name = "Tmin", values = c("#882255","darkorange","#6699cc","#117733"))+ scale_x_continuous("Hour of schedule ",breaks=seq(0,14,1))+scale_y_continuous("Temperature (C°)") +   theme(legend.position="top")

#*******************************************************************************
#### 2. Probability of seed germination ###
#*******************************************************************************

Fseed <-read.csv("LaboratoryExp_Seed.csv", header=T)

head(Fseed,21)

sapply(Fseed,class)

Fseed$Elev_km<-Fseed$elevation/1000 #converts elevation m to km
Fseed$Selev<-scale(Fseed $elevation ,center=TRUE, scale=TRUE)  #scales elevation

Fseed$genotype<-as.factor(Fseed$genotype)
Fseed$tray<-as.factor(Fseed$tray)
Fseed$cell<-as.factor(Fseed$cell)
Fseed$plantID<-as.factor(Fseed$plantID)
Fseed$Tmin<-as.factor(Fseed$Tmin)

#logistic regression for germination
FGermination_model<-glmmTMB(Germinated~Elev_km*Tmin +(1|genotype)+(1| tray) +(1|cell:tray),data=Fseed,family=binomial(link="logit"))

Anova(FGermination_model,type="III")

#save(FGermination_model,file="FGermination_model.rda")

##Testing the random effect of genotype
#FGermination_model_nogeno<-glmmTMB(Germinated~Elev_km*Tmin +(1| Tray) +(1|cell:Tray),data=freeze,family=binomial(link="logit"))
#anova(FGermination_model, FGermination_model_nogeno)

##Testing the random effect of block (tray)
#FGermination_model_noblock<-glmmTMB(Germinated~Elev_km*Tmin + (1|Genotype) +(1|cell:Tray),data=freeze,family=binomial(link="logit"))
#anova(FGermination_model, FGermination_model_noblock)

##Testing the random effect of cell:Tray
#FGermination_model_noCT<-glmmTMB(Germinated~Elev_km*Tmin + (1|Genotype) +(1|Tray) ,data=freeze,family=binomial(link="logit"))
#anova(FGermination_model, FGermination_model_noCT)


#Figure 4a: probability of germination for laboratory experiment

plot(predictorEffects(FGermination_model, ~ Elev_km), type="response",partial.residuals=FALSE, confint=list(style="auto"), xlab="Population source elevation (Km)", ylab="Probability of seed germination",line=list(multiline=TRUE, lty=1:4,col=c("#882255","darkorange","#6699cc","#117733")))


# Obtain slopes for each treatment level. These slopes are exponentiated to calculate odds ratios

FSeedOdds<- emtrends(FGermination_model, specs = c("Tmin"), var = "Elev_km")
FSeedOdds<- as.data.frame(summary(FSeedOdds))[c('Tmin', 'Elev_km.trend', 'SE')]
FSeedOdds %>% mutate(
          OddsRatio = exp(Elev_km.trend),
          Lower95 = OddsRatio * exp(-1.96*SE),
          Upper95 = OddsRatio * exp(1.96*SE))

#*******************************************************************************
#### 3. Seedling analysis: Survival ###
#*******************************************************************************


seedling <-read.csv("LaboratoryExp_Seedling.csv", header=T)
sapply(seedling,class)

seedling$Elev_km<-seedling$elevation/1000
seedling$Selev<-scale(seedling $elevation ,center=TRUE, scale=TRUE)

seedling$initsize<-scale(seedling $init_size ,center=TRUE, scale=TRUE)
seedling$genotype<-as.factor(seedling$genotype) 
seedling$tray<-as.factor(seedling$tray) #Block 
seedling$leaf<-as.factor(seedling$leaf)#seedling class, two or six leaves

seedling $Tmin<-factor(seedling $Tmin)


#Logistic regression for Suvival till end of experiment

Fsurvival_model<-glmmTMB(Final_status ~initsize+Tmin*Selev +(1|tray)+(1|genotype),data=seedling, family=binomial)

Anova(Fsurvival_model,type="III")
summary(Fsurvival_model)

#tukey test
summary(glht(Fsurvival_model, linfct = mcp(Tmin = "Tukey")))


##Testing the random effect of genotype
#Fsurvival_nogeno<-glmmTMB(Final_status ~initsize+Selev*Tmin +(1|Tray),data=seedling, family=binomial)
#anova(Fsurvival_model, Fsurvival_nogeno)

##Testing the random effect of block (tray)
#Fsurvival_noblock<-glmmTMB(Final_status ~initsize+Selev*Tmin +(1|Genotype),data=seedling, family=binomial)
#anova(Fsurvival_model, Fsurvival_noblock)


#Figure 4b: Seedling survival

vioplot(Final_status ~ Tmin, data= seedling, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,ylim=c(0,1),colMed = "black",colMed2 = c("#117733","#6699cc","darkorange","#882255"), col=c("#117733","#6699cc","darkorange","#882255"), ylab="Probability of seedling survival", xlab="Minimum temperature (°C)")



#*******************************************************************************
#### 5. Seedling analysis: Relative growth rate ###
#*******************************************************************************

seedling <-read.csv("LaboratoryExp_Seedling.csv", header=T)
sapply(seedling,class)

seedling$Elev_km<-seedling$elevation/1000
seedling$Selev<-scale(seedling $elevation ,center=TRUE, scale=TRUE)

seedling$initsize<-scale(seedling $init_size ,center=TRUE, scale=TRUE)
seedling$genotype<-as.factor(seedling$genotype) 
seedling$tray<-as.factor(seedling$tray) #Block 
seedling$leaf<-as.factor(seedling$leaf)#seedling class, two or six leaves


RGR <-lmer (RGR_total~ Tmin*Selev+leaf 
            +(1|tray)
            +(1|genotype)
            , data = seedling,
            control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))

Anova(RGR,type="III")

summary(glht(RGR, linfct = mcp(Tmin = "Tukey")))

#Testing the random effect of genotype

#FRGRnogeno <-lmer (RGR_total~ leaf+Tmin*Selev+(1|Tray), data = seedling,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))


#anova(RGR, FRGRnogeno)

##Testing the random effect of block (tray)
#FRGRnoblock <-lmer (RGR_total~ leaf+Tmin*Selev+(1|Genotype), data = seedling,control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))


#anova(RGR, FRGRnoblock)

#Figure 4C

seedlingRGR<-na.omit(seedling[,c('RGR_total','Tmin','treatment')])

seedlingRGR $Tmin<-factor(seedlingRGR $Tmin, levels = c("4", "-5.7","-10.7"))


vioplot(RGR_total ~ Tmin, data= seedlingRGR, plotCentre = "point",  pchMed = 23,  horizontal= FALSE,colMed = "black",colMed2 = c("#117733","#6699cc","darkorange"), col=c("#117733","#6699cc","darkorange"), ylab="Relative growth rate", xlab="Minimum temperature (°C)")+stripchart(RGR_total ~ Tmin, data= seedlingRGR,  method = "jitter", col = alpha("black", 0.2), pch=16 ,vertical = TRUE, add = TRUE)




