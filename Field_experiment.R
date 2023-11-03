######## PROJECT: Field experiment: Germination and seedling performance under climate change
#### PURPOSE:Examine germination success, seedling survival and germination phenology in response to snowmelt timing in the field.
#### AUTHOR: Jill Anderson
# AUTHOR:Inam Jameel
#### DATE LAST MODIFIED: 2 Nov 23

# remove objects and clear workspace
rm(list = ls(all=TRUE))


#require packages
require(lme4)
require(car)
require(lmerTest)
require(MuMIn)
require(MASS)
require(ggplot2)
require(optimx)
require(blme)
require(dplyr)
require(visreg)
require(effects)
require(performance)
require(gamlss)
require(glmmTMB)
require(glmmTMB) #have to load twice
require(emmeans) #for extracting slopes
require(broom.mixed) #for extracting slopes for phenology model

# set working directory
#setwd("~/Documents/freezing tolerance")
setwd("/Users/inam/OneDrive - University of Georgia/FrztolOnedrive/Freeze tolerance/data")


#*******************************************************************************
#### 1. Feed in the datasets #####
#*******************************************************************************

#Read in climatic data
climate <-read.csv("climate_all_gardens.csv", header=T)
climate $year<-as.factor(climate $Year)
climate $Garden_Block<-as.factor(climate $Garden_Block)
climate $Garden<-as.factor(climate $Garden)
sapply(climate,class)

##Standardize climatic variables to a mean of 0 and SD of 1
climate$snow_melt <- (climate $SnowMeltDOY - mean(climate $SnowMeltDOY, na.rm = TRUE)) / sd(climate $SnowMeltDOY,na.rm = TRUE)
climate$FDD <- (climate $FDD_early - mean(climate $FDD_early, na.rm = TRUE)) / sd(climate $FDD_early,na.rm = TRUE)
climate$GDD <- (climate $GDD_nosnow - mean(climate $GDD_nosnow, na.rm = TRUE)) / sd(climate $GDD_nosnow,na.rm = TRUE)
climate$CWD <- (climate $CWD_snow_min - mean(climate $CWD_snow_min, na.rm = TRUE)) / sd(climate $CWD_snow_min,na.rm = TRUE)


##Read in field data. has climate data for each year, block, garden
seed <-read.csv("FieldExp.csv", header=T,stringsAsFactors=TRUE)

sapply(seed,class)
seed$Garden_Block<-as.factor(seed$Garden_Block)
seed$Garden<-as.factor(seed$Garden)
seed$Treatment<-as.factor(seed$treatment)
seed$Selev <- (seed $Elevation_Population - mean(seed $Elevation_Population, na.rm = TRUE)) / sd(seed $Elevation_Population,na.rm = TRUE)
seed$elev_km <- (seed $Elevation_Population) / 1000
seed$Germ_days<-seed $Date_of_germination - seed $SnowMeltDOY
seed$germ_phen<-log(seed $Date_of_germination)


##Concatenate garden block and year

seed $Garden_Block_year <-interaction(seed $Garden_Block, seed $year,sep = "_")


#*******************************************************************************
#### 2. Analyze freezing degree days across gardens and generate a figure #####
#*******************************************************************************
##This analysis uses all gardens from 2014-2022. To use a Gamma distribution,  we have to increase the value of the zeros

climate$FDD<-climate$FDD_early+0.001

freezing<-glmer(FDD ~ Garden+ SnowMeltDOY +(1|year),family=Gamma(link="log"), data= climate)
Anova(freezing,type="III")

#To test random effect
freezing_noyear<-glm(FDD ~ Garden+ snow_melt, family=Gamma,data= climate)
anova(freezing, freezing_noyear)

#Figure 1
snowmelt_FDD =ggplot(climate, aes(x= SnowMeltDOY,y= FDD,shape=Garden, color=Garden))+geom_point(aes(shape= Garden),size=4)+scale_shape_manual(values = c(16,2,1,17,5)) +scale_x_continuous("Day of snowmelt (ordinal day of year)")+ scale_y_continuous("Accumulated Freezing Degree Units (Growing Degree Days) 0-60 days post snowmelt")
snowmelt_FDD +theme_bw()+ theme_bw()+theme(text = element_text(size=15),axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), panel.border = element_blank(),panel.grid.major =element_blank(), panel.grid.minor=element_blank(),legend.position = "bottom")+ scale_color_manual(values=c("#D55E00", "#009E73", "#CC79A7", "#E69F00","#56B4E9"))+geom_smooth(method="glm", method.args=list(family="poisson"),se=FALSE,  size=1, aes(group=1),colour="black")



###To calculate the average snowmelt times and projected snowmelt times in inset table in figure 3

control<-subset(climate, treat =='control')

c14_19<-subset(control, !(Year %in% c("2020","2021","2022"))) #subset to focus on years with field data

current_snow <- aggregate(c14_19 $SnowMeltDOY, list(c14_19 $Garden), FUN=mean) 

current_snow <- current_snow %>% 
  
  rename(Garden = Group.1,   SnowMeltDOY= x)

###If projected snowmelts are 32 days earlier than these values by 2070

current_snow$ Projected_snowmelt <-current_snow$SnowMeltDOY-32

current_snow


#*******************************************************************************
#### 3. Germination phenology as a function of snowmelt timing #####
#*******************************************************************************
#have to remove NA's for using gamlss
pheno_data<-na.omit(seed[,c('Garden_Block_year','year','Garden','elev_km','Genotype','SnowMeltDOY','Date_of_germination')])
pheno_data$Genotype<-as.factor(pheno_data$Genotype)
pheno_data$year<-as.factor(pheno_data$year)

phenology<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden+ re(random = ~1|Genotype)+ re(random = ~1| year), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))


##Testing random effects 
#phenology_test1<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden, family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))


#phenology_test2<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden+ re(random = ~1|Genotype), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))

#phenology_test3<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden+ re(random = ~1|year), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))

##This one does not converge
#phenology_test4<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden+ re(random = ~1| Garden_Block_year), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))

##This model also fails to converge
#phenology_test5<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden+ re(random = ~1|Garden_Block_year)+ re(random = ~1|Genotype)+ re(random = ~1| Garden_Block_year), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))


#phenology_test6<-gamlss(Date_of_germination ~ SnowMeltDOY*Garden+elev_km*Garden+ re(random = ~1|year)+ re(random = ~1|Genotype), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))

#lapply(list(phenology_test1, phenology_test2, phenology_test3, phenology_test6),AIC)

#figure 3a
visreg(phenology,"SnowMeltDOY", by="Garden", overlay=TRUE,  xlab="Day of snowmelt (ordinal day)", ylab="Germination phenology (ordinal day of first germination)", partial=TRUE, type="conditional",line=list(lty=1:5,col=c("#D55E00","#009E73","#E69F00")), ylim=c(78,162), points=list(col=c("#D55E00","#009E73","#E69F00")),fill=list(col=grey(c(0.8), alpha=0.4)))


#Extract slopes
##Estess: Garden at 2553m
summary(phenology)
coef(phenology)
tidy(phenology,conf.int=TRUE,exponentiate=FALSE,effects="fixed")


##Changing baseline to extract slopes
pheno_data$PM<-factor(pheno_data $Garden, levels = c("PeanutMine","Estess","Gothic"))
pheno_data$G<-factor(pheno_data $Garden, levels = c("Gothic","PeanutMine","Estess"))

phenologyPM<-gamlss(Date_of_germination ~ SnowMeltDOY* PM +elev_km* PM + re(random = ~1|year)+ re(random = ~1|Genotype), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))
coef(phenologyPM)
summary(phenologyPM)
tidy(phenologyPM,conf.int=TRUE,exponentiate=FALSE,effects="fixed")


phenologyG<-gamlss(Date_of_germination ~ SnowMeltDOY* G +elev_km* G + re(random = ~1|year)+ re(random = ~1|Genotype), family=LO(mu.link = "identity", sigma.link = "log"), data=pheno_data,control = gamlss.control(n.cyc = 500))
coef(phenologyG)
summary(phenologyG)
tidy(phenologyG,conf.int=TRUE,exponentiate=FALSE,effects="fixed")



#*******************************************************************************
#### 4. Seed germination ###
#*******************************************************************************
Germination_model<-glmmTMB(Germination~SnowMeltDOY*Garden+elev_km*Garden+I(elev_km^2) +(1|Genotype)+(1|Garden_Block_year)+(1|year),data=seed,family=binomial(link="logit"))
Anova(Germination_model,type="III")


#figure 3b
visreg(Germination_model,"SnowMeltDOY", by="Garden", overlay=TRUE,  scale = "response", xlab="Day of snowmelt (ordinal day)", ylab="Probability of germination", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#D55E00","#009E73","#E69F00")), points=list(col=c("#D55E00","#009E73","#E69F00")),fill=list(col=grey(c(0.8), alpha=0.4)))

#save(Germination_model,file="Germination_model.rda")

##Testing the random effect of genotype
#Germination_model_nogeno<-glmmTMB(Germination~SnowMeltDOY*Garden+elev_km*Garden+I(elev_km^2) +(1|Garden_Block_year)+(1|year),data=seed,family=binomial(link="logit"))
#anova(Germination_model, Germination_model_nogeno)

##Testing the random effect of block
#Germination_model_noblock<-glmmTMB(Germination~SnowMeltDOY*Garden+elev_km*Garden+I(elev_km^2) +(1|Genotype)+(1|year),data=seed,family=binomial(link="logit"))
#anova(Germination_model, Germination_model_noblock)

##Testing the random effect of year
#Germination_model_noyear<-glmmTMB(Germination~SnowMeltDOY*Garden+elev_km*Garden+I(elev_km^2) +(1|Genotype)+(1|Garden_Block_year),data=seed,family=binomial(link="logit"))
#anova(Germination_model, Germination_model_noyear)



# Obtain slopes for each garden. These slopes are exponentiated to calculate odds ratios

germinationOdds<- emtrends(Germination_model, specs = c("Garden"), var = "SnowMeltDOY")
germinationOdds<- as.data.frame(summary(germinationOdds))[c('Garden', 'SnowMeltDOY.trend', 'SE')]
germinationOdds$Garden_elevation <- c(2553, 2890, 2710) #add garden elevations, estess, gothic, peanut mine
germinationOdds <- germinationOdds[order(germinationOdds$Garden_elevation, decreasing = FALSE), ]%>% mutate(
  OddsRatio = exp(SnowMeltDOY.trend),
  Lower95 = OddsRatio * exp(-1.96*SE),
  Upper95 = OddsRatio * exp(1.96*SE))

germinationOdds


#*******************************************************************************
#### 5. Seedling survival ###
#*******************************************************************************
seedling<-subset(seed,Germination=="1")

Established<-glmmTMB(Established_first_year~SnowMeltDOY*Garden+elev_km*Garden +(1|Genotype)+(1|Garden_Block_year)+(1|year),data=seedling,family=binomial(link="logit"))
Anova(Established,type="III")
#save(Established,file="Established.rda")

##Testing the random effect of genotype
#Established_first_year_model_nogeno<-glmmTMB(Established_first_year~SnowMeltDOY*Garden+elev_km*Garden +(1|Garden_Block_year)+(1|year),data= seedling,family=binomial(link="logit"))
#anova(Established, Established_first_year_model_nogeno)

##Testing the random effect of block
#Established_first_year_model_noblock<-glmmTMB(Established_first_year~SnowMeltDOY*Garden+elev_km*Garden +(1|Genotype)+(1|year),data= seedling,family=binomial(link="logit"))
#anova(Established, Established_first_year_model_noblock)

##Testing the random effect of year
#Established_first_year_model_noyear<-glmmTMB(Established_first_year~SnowMeltDOY*Garden+elev_km*Garden +(1|Genotype)+(1|Garden_Block_year),data= seedling,family=binomial(link="logit"))
#anova(Established, Established_first_year_model_noyear)

#figure 3c
visreg(Established,"SnowMeltDOY", by="Garden", overlay=TRUE,   scale = "response", xlab="Day of snowmelt (ordinal day)", ylab="Probability of seedling survival", partial=FALSE,type="conditional",line=list(lty=1:3,col=c("#D55E00","#009E73","#E69F00")), points=list(col=c("#D55E00","#009E73","#E69F00")),fill=list(col=grey(c(0.8), alpha=0.4)))


# Obtain slopes for each garden. These slopes are exponentiated to calculate odds ratios

establishmentOdds<- emtrends(Established, specs = c("Garden"), var = "SnowMeltDOY")
establishmentOdds<- as.data.frame(summary(establishmentOdds))[c('Garden', 'SnowMeltDOY.trend', 'SE')]
establishmentOdds$Garden_elevation <- c(2553, 2890, 2710) #add garden elevations, estess, gothic, peanut mine

establishmentOdds<- establishmentOdds[order(establishmentOdds$Garden_elevation, decreasing = FALSE), ]%>% mutate(
  OddsRatio = exp(SnowMeltDOY.trend),
  Lower95 = OddsRatio * exp(-1.96*SE),
  Upper95 = OddsRatio * exp(1.96*SE))

establishmentOdds

#*******************************************************************************
#### 6. Seedling growth ###
#*******************************************************************************
seedling<-subset(seed,Germination=="1")

growth<-lmer(final_true_leaves ~ Date_of_germination +SnowMeltDOY*Garden+elev_km*Garden +(1|Genotype)+(1|Garden_Block_year)+(1|year),data=seedling)
Anova(growth,type="III")

#figure 3d
visreg(growth,"SnowMeltDOY", by="Garden", overlay=TRUE,   scale = "response", xlab="Day of snowmelt (ordinal day)", ylab="Final size (leaf number)", partial=TRUE,type="conditional",line=list(lty=1:3,col=c("#D55E00","#009E73","#E69F00")), points=list(col=c("#D55E00","#009E73","#E69F00")),fill=list(col=grey(c(0.8), alpha=0.4)))

##Testing the random effect of genotype
#growth_first_year_model_nogeno<-lmer(final_true_leaves ~Date_of_germination +SnowMeltDOY*Garden+elev_km*Garden +(1|Garden_Block_year)+(1|year),data= seedling)
#anova(growth, growth_first_year_model_nogeno)

##Testing the random effect of block
#growth_first_year_model_noblock<-lmer(final_true_leaves ~Date_of_germination +SnowMeltDOY*Garden+elev_km*Garden +(1|Genotype)+(1|year),data= seedling)
#anova(growth, growth_first_year_model_noblock)

##Testing the random effect of year
#growth_first_year_model_noyear<-lmer(final_true_leaves ~Date_of_germination +SnowMeltDOY*Garden+elev_km*Garden +(1|Genotype)+(1|Garden_Block_year),data= seedling)
#anova(growth, growth_first_year_model_noyear)

# Obtain slopes for each garden. 

growthslopes<- emtrends(growth, specs = c("Garden"), var = "SnowMeltDOY")
growthslopes<- as.data.frame(summary(growthslopes))[c('Garden', 'SnowMeltDOY.trend', 'SE')]
growthslopes$Garden_elevation <- c(2553, 2890, 2710) #add garden elevations, estess, gothic, peanut mine
growthslopes<- growthslopes[order(growthslopes$Garden_elevation, decreasing = FALSE), ]%>% mutate(
  Lower95 = SnowMeltDOY.trend - (1.96*SE),
  Upper95 = SnowMeltDOY.trend + (1.96*SE))

growthslopes





