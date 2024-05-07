################################################################################
#
#   Cod Otoliths Document Accelerating Climate Impacts in the Baltic Sea
#
#   Yvette Heimbrand, Karin Limburg, Karin Hüssy, Tomas Næraa & Michele Casini
#
#   Currently under revision at Scientific reports 2024-04-29
#
#  *****************************************************************************
#
#  Figure 1      Hypoxia exposure (Mn:Mg) over time ~ anoxic volume and salinity
#                Repeated measurement mixed model
#                 
#  *****************************************************************************
#
#  Figure 2       Hypoxia exposure (Mn:Mg) ~ salinity (Sr:Ca) for different age classes
#                 Pairwise comparisons using Wilcoxon rank sum test
#
#  *****************************************************************************
#
#  Figure 3       Metabolic status (Mg:Ca) over time ~ BY15 DO sat % 30-90 m + BY15 Temp 30-90 m Q4
#                 Repeated measurement mixed model - Age class 1-8 
#
#  *****************************************************************************
#
#  Figure 4       Predicted mean length at age
#
#  *****************************************************************************  
#
#  Figure 5      Environmental parameters (Fig 5 A)
#                The R code and data for the maps in Fig 5 D-F are available upon request. 
#
#
#  ***************************************************************************** 
#
#              R code by Yvette Heimbrand SLU Aqua 2024-04-02
#
#
########################################################################################################
#
#  Figure 1       Hypoxia exposure (Mn:Mg) over time versus anoxic volume and salinity
#                 Repeated measurement mixed model
#
#########################################################################################################
#version[['version.string']] # 4.1.2
#
rm(list = ls()) # Remove all variables from the current environment
options("install.lock"=FALSE) # Run in case packages will not install in library

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RcppRoll)
library(zoo)
library(varhandle)
library(plyr)

#Read data files
dfr<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes
env<- read.csv("Data\\Data_environmental_factors.csv", sep = ";") 

# Subset Neolithic Stone Age data
Neo <- subset(dfr, Decades =="Neolithic")
Neo_old <- subset(Neo, Age >"0")
summary(Neo_old$Avg_MnMg06)
sd(Neo_old$Avg_MnMg06)

Cal_Yearn <- c(1925:2020)
df <- data.frame(Cal_Yearn)

N1<-df %>%
  mutate(mean = 0.16271) #Mean value for Neolithic time period
head(N1)

#Standard deviation (sd) min
df2 <- data.frame(Cal_Yearn)
sd_min<-(0.16271-0.08785458)
head(sd_min) #0.07485542

#Standard deviation (sd) max
df3 <- data.frame(Cal_Yearn)
sd_max<-(0.16271+0.08785458)
head(sd_max) #0.2505646 

#Create data frame for mean and sd
N1<-df %>%
  mutate(mean = 0.16271) #Mean value for Neolithic time period
head(N1)
N2<-N1 %>%
  mutate(sd_min = 0.07485542) #Mean value for Neolithic time period
head(N2)
NEO<-N2 %>%
  mutate(sd_max = 0.2505646) #Mean value for Neolithic time period
head(NEO)
NEOb<-NEO %>%
  mutate(name = "rib") #Mean value for Neolithic time period
head(NEOb)

# Exclude Neolithic samples from data frame
dg <- subset(dfr, Decades !="Neolithic")

# Adjust data frame for modern samples
data <- subset(dg, Age !="NA") # Exclude NA
data <- subset(data, Age >"0") # Exclude young-of-year=Age 0
data$Cal_Years<-as.factor(data$Cal_Year) 
data$Cal_Yearn<-as.numeric(data$Cal_Year) 
data$Yearhyp <-as.factor(data$Year_hyp)
data$Agef <-as.factor(data$Age)

#List data
str(data, list.len=Inf)

# Moving average 3 for the anoxic volume
env$roll_anox <- roll_mean(env$Hypoxic.volume.0, n = 3, align = "center", fill = NA)
env$Cal_Years<-as.factor(env$Cal_Year)
env$Cal_Yearn<-as.numeric(env$Cal_Year)

# Moving average 3 for Mn:Mg for age class 3
data_age_3 <- subset(data, Age =="3")
ab<-aggregate(Avg_MnMg06 ~ Cal_Yearn, data = data_age_3, mean, na.rm=TRUE)
list(ab)

Cal_Yearn <- c(1929:2019)
DF <- data.frame(Cal_Yearn)
DF2 <- merge(DF, ab, by=c('Cal_Yearn'),all=T)
list(DF2)
DF3<-transform(DF2, mean_demand = rollmeanr(Avg_MnMg06, 3, fill = NA)) 
list(DF3)

data_N<-merge(data, NEOb, by='Cal_Yearn')
str(data_N)

data$Cal_Yearf<-as.factor(data$Cal_Yearn)
tre <- subset(data, Age =="3")
tre_mean<- aggregate(Avg_MnMg06 ~ Cal_Yearf, data = tre, mean, na.rm=TRUE)

tre_mean$Cal_Yearn<-unfactor(tre_mean$Cal_Yearf)
head(tre_mean)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
 
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
tgc <- summarySE(tre, measurevar="Avg_MnMg06", groupvars=c("Cal_Yearn"))
tgc

# Adjust secondary y-axis in plot
ylim.prim <- c(0,2)   
ylim.sec <- c(-1100,2800) 
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

# Plot to figure 1
p<-ggplot(data_N,aes(x=Cal_Yearn, y=Avg_MnMg06)) + 
  geom_hline(yintercept = 0.16271, color = "red",linetype = "dashed", size=1.1)+
  geom_ribbon(aes(ymin =sd_min, ymax = sd_max), fill = "red",alpha=0.1)+
  geom_jitter(aes(colour=Agef), alpha=0.8, width=0.1,size=2) + 
  scale_color_viridis_d(option = "mako",name = "Age",direction=-1)+
  geom_line(aes(y = mean_demand),data=DF3, color="blue",size=1.3,linetype="solid")+
  geom_point(data=tre_mean,aes(x=Cal_Yearn, y=Avg_MnMg06)) + 
  geom_errorbar(data=tgc, aes(ymin=Avg_MnMg06-se, ymax=Avg_MnMg06+se), width=.1) +
  geom_line(aes(y=a+roll_anox*b, group=1),data=env, color="black",size=1.6, stat = "summary", fun="mean",linetype="solid")+
  coord_cartesian(ylim=c(0,1.7))+
  theme_classic()+
  scale_size_manual(values = c(5,1))+ 
  scale_y_continuous("Annual mean Mn:Mg",sec.axis = sec_axis(~(.-a)/b,name=expression(Anoxic~volume~km^3)))+ 
  theme(axis.text.x = element_text(angle = 90))+ 
  scale_x_continuous(breaks=seq(1925,2020,by=5))+
  xlab("Year")


p + theme(axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,10,0,0),face="plain"),
          axis.title.y.right = element_text(family = "sans", size = 15, margin=margin(0,10,0,0),face="bold"),
          panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title.x = element_text(family = "sans", size = 15, margin=margin(10,0,0,0),face="plain"), 
          axis.text = element_text(family = "sans", size = 13,face="plain"),
          axis.text.y.left = element_text(family = "sans", size = 13, margin=margin(10,0,0,0),face="plain"),
          axis.text.y.right = element_text(family = "sans", size = 13, margin=margin(0,10,0,0),face="plain"),
          plot.title = element_text(hjust = 0.5,family = "sans", size = 15, margin=margin(0,0,10,0),face="bold"),
          panel.background = element_rect(fill = NA))+theme(legend.text = element_text(size = 14))


# Save plot
ggsave(path = "figs", filename = "PLOT_1_Anoxic_MnMg_with_se.png", dpi = 600, width = 12, height = 7)

######################################################################################################
#
#  Repeated measures mixed model for figure 1
#
######################################################################################################
#
# Empty data frame
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library
# Load libraries
library(car)
library(ggplot2)
library(lme4)
library(report)
library(effects)
library(sjPlot)
library(r2glmm)
library(report)

#Read file
dd<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes
env<- read.csv("Data\\Data_environmental_factors.csv", sep = ";") 

# merge two data frames by Year
Tot <- merge(dd,env,by="Cal_Year", all.x = TRUE)
list(Tot)

str(Tot)
dg <- subset(Tot, Decades !="Neolithic")
df <- subset(dg, Cal_Year >"1959") 
str(df) 

# Check distribution of Mn:Mg
summary(df$Avg_MnMg06)
hist(df$Avg_MnMg06)
qqPlot(df$Avg_MnMg06)

# Log transform Mn:Mg
df$MnMg_log<-log(df$Avg_MnMg06)
hist(df$MnMg_log) # Better
#ggsave(path = "figs", filename = "Log_MnMg_histogram.png", dpi = 600, width = 6, height = 5)

qqPlot(df$MnMg_log)
#ggsave(path = "figs", filename = "Log_MnMg_qqplot.png", dpi = 600, width = 6, height = 5)

df$Age_num<-as.numeric(df$Ages)
df$Age_n<-(df$Age_num-1)
summary(df$Age_n)

x_subset <- df[df$Age_n >=1 ,]
summary(x_subset$Age_n)
str(x_subset)
x_subset$Anox_vol<-(x_subset$Hypoxic.volume.0.y/1000) # Dividing with 1000 for km3
str(x_subset, list.len=Inf)
#
### Repeated measurements mixed model for log Mn:Mg with salinity (stn BY15 30-90 m), anoxic volume and age as fixed factors and calendar year and fish ID nested within age class as the random effects.
#
# Best model from comparing models with this being the model with lowest AIC and significant results
g3.mixed<-lmer(MnMg_log~ BY15.Q4.salinitet.30_90m.y + Anox_vol +Age_n+(Age_n|Fish_ID)+(1|Cal_Year),data=x_subset,na.action=na.exclude)

Anova(g3.mixed)

#Check Multicollinearity
vif(g3.mixed) #All vif < 5 = No multicolinearity

# Check effects
sjPlot::plot_model(g3.mixed)
ggsave(path = "figs", filename = "Log_MnMg_model_sjplot.png", dpi = 600, width = 6, height = 5)

plot(allEffects(g3.mixed))
ggsave(path = "figs", filename = "Log_MnMg_model_alleffects.png", dpi = 600, width = 6, height = 5)

r2 <- r2beta(model=g3.mixed,partial=TRUE,method='sgv')
print(r2)
plot(x=r2)
ggsave(path = "figs", filename = "Log_MnMg_model R2.png", dpi = 600, width = 6, height = 5)

# Report model statistics
report(g3.mixed) # Report function for model statistics


########################################################################################################
########################################################################################################
#
#  Figure 2       Hypoxia exposure (Mn:Mg) versus salinity (Sr:Ca) for different age classes
#                 Pairwise comparisons using Wilcoxon rank sum test
#*******************************************************************************************************
#*
# Empty data frame
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library

# Load libraries
library(ggplot2)

# Find Neolithic baselines
data<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes

str(data)
summary(data$Age)
# All age classes
newdata27 <- subset(data, ICES_SD=="27")
str(newdata27)

# Age class 0
newdata270 <- subset(newdata27, Age<"1")
summary(newdata270$Age)

Cal_Yearn <- c(1925:2020)
df <- data.frame(Cal_Yearn)

# Sr:Ca
summary(newdata270$Avg_SrCa)  # Mean Sr:Ca = 3.992
sd(newdata270$Avg_SrCa) ## sd Sr:Ca = 0.5474722

#Sd min Sr:Ca
sd_min_Sr_0<-(3.992-0.5474722)
head(sd_min_Sr_0) # 3.444528

#Sd max Sr:Ca
sd_max_Sr_0<-(3.992+0.5474722)
head(sd_max_Sr_0) # 4.539472 

# Mn:Mg
summary(newdata270$Avg_MnMg06) # Mean Mn:Mg = 0.24568
sd(newdata270$Avg_MnMg06) # 0.102091

#Sd min Sr:Ca
sd_min_Mn_Mg_0<-(0.24568-0.102091)
head(sd_min_Mn_Mg_0) # 0.143589

#Sd max Sr:Ca
sd_max_Mn_Mg_0<-(0.24568+0.102091)
head(sd_max_Mn_Mg_0) # 0.347771 

#Create data frame for mean and sd
N1<-df %>%
  mutate(mean_Sr = 3.992) #Mean value for Neolithic time period
head(N1)
N2<-N1 %>%
  mutate(sd_min_Sr_0 = 3.444528) #Mean value for Neolithic time period
head(N2)
N3<-N2 %>%
  mutate(sd_max_Sr_0 = 4.539472) #Mean value for Neolithic time period
head(N3)
N4<-N3 %>%
  mutate(mean_Mn_Mg = 0.24568) #Mean value for Neolithic time period
head(N4)
N5<-N4 %>%
  mutate(sd_min_Mn_Mg_0 = 0.143589) #Mean value for Neolithic time period
head(N5)
N6<-N5 %>%
  mutate(sd_max_Mn_Mg_0 = 0.347771) #Mean value for Neolithic time period
head(N6)

# Subset data for young of year = Age class 0 for modern samples
newdata0 <- subset(data, Age<"1")
newdata0w <- subset(newdata0, ICES_SD!="27") #Exclude Neolithic samples
str(newdata0w)

newdata0w$Cal_Yearn<-as.numeric(newdata0w$Cal_Year)
data_df<-merge(newdata0w, N6, by='Cal_Yearn')

#Plot Sr:Ca for Age class 0 over time in SD 25 and 28
ggplot(newdata0w, aes(x=as.numeric(as.character(Cal_Year)), y= Avg_SrCa)) + 
  geom_hline(yintercept = 3.992, color = "red",linetype = "dashed", size=1.2)+
  geom_ribbon(aes(ymin =sd_min_Sr_0, ymax = sd_max_Sr_0), fill = "red",alpha=0.1)+
  geom_point(colour = "lightskyblue2",size=2) +
  geom_smooth(method="loess",colour = "red", size = 2,fill = "orange1")+
  facet_wrap(.~ICES_SD)+ylab("Mean Sr:Ca")+xlab("Year")+coord_cartesian(ylim=c(1,10))+ 
  theme_classic(base_size = 23)+ theme(axis.text.x = element_text(angle = 90))+
  theme(text = element_text(size = 30)) +theme(legend.position="none")

#Save plot
ggsave(path = "figs", filename = "SrCa_YOY_plot_2.png", dpi = 600, width = 8, height = 6)

#Plot Mn:Mg for Age class 0 over time in SD 25 and 28
ggplot(newdata0w, aes(x=as.numeric(as.character(Cal_Year)), y= Avg_MnMg06)) + 
  geom_hline(yintercept = 0.24568, color = "red",linetype = "dashed", size=1.2)+
  geom_ribbon(aes(ymin =sd_min_Mn_Mg_0, ymax = sd_max_Mn_Mg_0), fill = "red",alpha=0.1)+
  geom_point(colour = "lightskyblue2",size=2) +
  geom_smooth(method="loess",colour = "red", size = 1.3,fill = "orange1")+
  facet_wrap(.~ICES_SD)+ylab("Mean Mn:Mg")+xlab("Year")+coord_cartesian(ylim=c(0, 1.5))+ 
  theme_classic(base_size = 23)+ theme(axis.text.x = element_text(angle = 90))+
  theme(text = element_text(size = 30)) +theme(legend.position="none")

#Save plot
ggsave(path = "figs", filename = "MnMg_YOY_plot_2.png", dpi = 600, width = 8, height = 6)
####

# Subset data for young of year = Age class 0.
newdata0 <- subset(data, Age<"1")

### SrCa_SD25
newdata_25_0 <- subset(newdata0, ICES_SD=="25")
cod_Sr_0_25<-newdata_25_0[,c("Avg_SrCa","Decades")]
head(cod_Sr_0_25)
cod_Sr_0_25$Decade_f<-as.factor(cod_Sr_0_25$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_SrCa ~ Decades, data = cod_Sr_0_25) 

## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_Sr_0_25$Avg_SrCa, cod_Sr_0_25$Decades, p.adjust.method="bonferroni")

###
# MnMg_SD25
###

cod_MnMg_0_25<-newdata_25_0[,c("Avg_MnMg06","Decades")]
head(cod_MnMg_0_25)
cod_MnMg_0_25$Decade_f<-as.factor(cod_MnMg_0_25$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_MnMg06 ~ Decades, data = cod_MnMg_0_25) 


## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_MnMg_0_25$Avg_MnMg06, cod_MnMg_0_25$Decades, p.adjust.method="bonferroni")

###
# Sr_SD28
###
newdata_28_0 <- subset(newdata0, ICES_SD=="28")
cod_Sr_0_28<-newdata_28_0[,c("Avg_SrCa","Decades")]
head(cod_Sr_0_28)
cod_Sr_0_28$Decade_f<-as.factor(cod_Sr_0_28$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_SrCa ~ Decades, data = cod_Sr_0_28) 


## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_Sr_0_28$Avg_SrCa, cod_Sr_0_28$Decades, p.adjust.method="bonferroni")

###
# MnMg_SD28
###
newdata_28_0 <- subset(newdata0, ICES_SD=="28")
cod_MnMg_0_28<-newdata_28_0[,c("Avg_MnMg06","Decades")]
head(cod_MnMg_0_28)
cod_MnMg_0_28$Decade_f<-as.factor(cod_MnMg_0_28$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_MnMg06 ~ Decades, data = cod_MnMg_0_28) 


## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_MnMg_0_28$Avg_MnMg06, cod_MnMg_0_28$Decades, p.adjust.method="bonferroni")

################################################################################
# Age class > 0
################################################################################
# Empty data frame
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library

# Load libraries
library(ggplot2)

data<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes

# Age class > 0 Neolithic
dff<- subset(data, ICES_SD=="27")
newdata_Neo_old <- subset(dff, Age>"0")
summary(newdata_Neo_old$Age) 

Cal_Yearn <- c(1925:2020)
df <- data.frame(Cal_Yearn)

# Sr:Ca
summary(newdata_Neo_old$Avg_SrCa)  # Mean Sr:Ca = # 5.775
sd(newdata_Neo_old$Avg_SrCa) ## sd Sr:Ca = 1.179022

#Sd min Sr:Ca
sd_min_Sr_old<-(5.775-1.179022)
head(sd_min_Sr_old) # 4.595978

#Sd max Sr:Ca
sd_max_Sr_old<-(5.775+1.179022)
head(sd_max_Sr_old) # 6.954022 

# Mn:Mg
summary(newdata_Neo_old$Avg_MnMg06) # Mean Mn:Mg = 0.16271
sd(newdata_Neo_old$Avg_MnMg06) # 0.08785458

#Sd min Mn:Mg
sd_min_Mn_Mg_old<-(0.16271-0.08785458)
head(sd_min_Mn_Mg_old) # 0.07485542

#Sd max Sr:Ca
sd_max_Mn_Mg_old<-(0.16271+0.08785458)
head(sd_max_Mn_Mg_old) # 0.2505646 


#Create data frame for mean and sd
N1<-df %>%
  mutate(mean_Sr_old = 5.775) #Mean value for Neolithic time period
head(N1)
N2<-N1 %>%
  mutate(sd_min_Sr_old = 4.595978) #Mean value for Neolithic time period
head(N2)
N3<-N2 %>%
  mutate(sd_max_Sr_old = 6.954022) #Mean value for Neolithic time period
head(N3)
N4<-N3 %>%
  mutate(mean_Mn_Mg_old = 0.16271) #Mean value for Neolithic time period
head(N4)
N5<-N4 %>%
  mutate(sd_min_Mn_Mg_old = 0.07485542) #Mean value for Neolithic time period
head(N5)
N6<-N5 %>%
  mutate(sd_max_Mn_Mg_old = 0.2505646) #Mean value for Neolithic time period
head(N6)

newdata1 <- subset(data, Age>"0")
newdata1w <- subset(newdata1, ICES_SD!="27")
summary(data$Age)
summary(newdata1w$Age)
str(newdata1w)

newdata1w$CalY_n<-as.numeric(newdata1w$Cal_Year)
newdata1w$Age_f<-as.factor(newdata1w$Age)
str(newdata1w)

newdata1w$Cal_Yearn<-as.numeric(newdata1w$Cal_Year)
data_df<-merge(newdata1w, N6, by='Cal_Yearn')

ggplot(newdata1w, aes(x=CalY_n, y= Avg_MnMg06)) + 
  geom_jitter(aes(colour=Age_f), alpha=0.8,width=0.1,size=2) + 
  scale_color_viridis_d(option = "mako",name = "Age",direction=-1)+
  geom_hline(yintercept = 0.16271, color = "red",linetype = "dashed", size=1.2)+
  geom_ribbon(aes(ymin =sd_min_Mn_Mg_old, ymax = sd_max_Mn_Mg_old), fill = "red",alpha=0.1)+
  geom_smooth(method="loess",colour = "red", size = 1.3,fill = "orange1")+
  facet_wrap(.~ICES_SD)+ylab("Mean Mn:Mg")+xlab("Year")+coord_cartesian(ylim=c(0, 1.5))+ 
  theme_classic(base_size = 23)+ theme(axis.text.x = element_text(angle = 90))+
  theme(text = element_text(size = 30)) +theme(legend.position="none")

ggsave(path = "figs", filename = "MnMg_old_plot_2.png", dpi = 600, width = 8, height = 6)

ggplot(newdata1w, aes(x=as.numeric(as.character(Cal_Year)), y= Avg_SrCa)) + 
  geom_hline(yintercept = 5.775, color = "red",linetype = "dashed", size=1.2)+
  geom_ribbon(aes(ymin =sd_min_Sr_old, ymax = sd_max_Sr_old), fill = "red",alpha=0.1)+
  geom_jitter(aes(colour=Age_f), alpha=0.8,width=0.1,size=2) + 
  scale_color_viridis_d(option = "mako",name = "Age",direction=-1)+
  geom_smooth(method="loess",colour = "red", size = 1.3,fill = "orange1")+facet_wrap(.~ICES_SD)+
  ylab("Mean Sr:Ca")+xlab("Year")+coord_cartesian(ylim=c(1,10))+ theme_classic(base_size = 23)+ 
  theme(axis.text.x = element_text(angle = 90))+
  theme(text = element_text(size = 30))+theme(legend.position="none")

ggsave(path = "figs", filename = "SrCa_old_plot_2.png", dpi = 600, width = 8, height = 6)

################################################################################
### Kruskal-Wallis & Wilcoxon rank sum test 
################################################################################
#*
# Empty data frame
rm(list = ls()) 
options("install.lock"=FALSE) # Run in case packages will not install in library

# Read data
data<- read.csv("Data_otolithchemistry.csv", sep = ";") # All age classes

newdata1 <- subset(data, Age>"0")

### SrCa_SD25
newdata_25_1 <- subset(newdata1, ICES_SD=="25")
cod_Sr_1_25<-newdata_25_1[,c("Avg_SrCa","Decades")]
head(cod_Sr_1_25)
cod_Sr_1_25$Decade_f<-as.factor(cod_Sr_1_25$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_SrCa ~ Decades, data = cod_Sr_1_25) 


## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_Sr_1_25$Avg_SrCa, cod_Sr_1_25$Decades, p.adjust.method="bonferroni")

###
# MnMg_SD25
###

cod_MnMg_1_25<-newdata_25_1[,c("Avg_MnMg06","Decades")]
head(cod_MnMg_1_25)
cod_MnMg_1_25$Decade_f<-as.factor(cod_MnMg_1_25$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_MnMg06 ~ Decades, data = cod_MnMg_1_25) 

## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_MnMg_1_25$Avg_MnMg06, cod_MnMg_1_25$Decades, p.adjust.method="bonferroni")

###
# Sr_SD28
###
newdata_28_1 <- subset(newdata1, ICES_SD=="28")
cod_Sr_1_28<-newdata_28_1[,c("Avg_SrCa","Decades")]
head(cod_Sr_1_28)
cod_Sr_1_28$Decade_f<-as.factor(cod_Sr_1_28$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_SrCa ~ Decades, data = cod_Sr_1_28) 


## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_Sr_1_28$Avg_SrCa, cod_Sr_1_28$Decades, p.adjust.method="bonferroni")

###
# MnMg_SD28
###
newdata_28_1 <- subset(newdata1, ICES_SD=="28")
cod_MnMg_1_28<-newdata_28_1[,c("Avg_MnMg06","Decades")]
head(cod_MnMg_1_28)
cod_MnMg_1_28$Decade_f<-as.factor(cod_MnMg_1_28$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Avg_MnMg06 ~ Decades, data = cod_MnMg_1_28) 


## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod_MnMg_1_28$Avg_MnMg06, cod_MnMg_1_28$Decades, p.adjust.method="bonferroni")
#
########################################################################################################
#        SALINITY CLASS versus HYPOXIA EXPOSURE
########################################################################################################
# Salinity class       Sr:Ca           Salinity
#       1               <3                <6
#       2               3 - 4.06         6-12
#       3               4.06 - 4,6        12-16
#       4              > 4.6              > 16
# 
#
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library

# Load libraries
library(ggplot2)
#
# Equation:  Srppm<-SrCa*380000/1000
#
# Equation:  salinity<-exp((log(Srppm)-6.241)/0.442)
#
# Salinity group 2
#  Srppm<-SrCa*380000/1000
Srppm<-3*380000/1000
Srppm #1140

salinity<-exp((log(1140)-6.241)/0.442)
salinity
#6.079608 PSU

# Salinity group 3
Srppm<-3.75*380000/1000
Srppm
#1425
salinity<-exp((log(1425)-6.241)/0.442)
salinity
#10.07231 PSU

# Salinity group 4
Srppm<-4.6*380000/1000
Srppm
# 1748
salinity<-exp((log( 1748)-6.241)/0.442)
salinity
#15.99071 PSU

# Read data file
data<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes

str(data)

# Creating baselines per salinity group for Neolithic samples
Neodata<-subset(data, Decade=="Neolithic")
summary(Neodata$New_Salinity_class)
Neodata$Salinity_classf<-as.factor(Neodata$New_Salinity_class)
levels(Neodata$Salinity_classf) <- c("Salinity < 6 PSU","Salinity 6 - 12 PSU","Salinity 12 - 16 PSU","Salinity > 16 PSU")

# Mean
grand.means <- aggregate(Avg_MnMg06 ~ Salinity_classf, data = Neodata, FUN = mean)
grand.means

# Sd
Neo_sd <- aggregate(Avg_MnMg06 ~ Salinity_classf, data = Neodata, FUN = sd)
Neo_sd

# Merging mean and sd per salinity group for Neolithic samples
data_merge<-merge(grand.means, Neo_sd, by='Salinity_classf')
data_merge
data_merge$sd_min<-(data_merge$Avg_MnMg06.x - data_merge$Avg_MnMg06.y)
data_merge$sd_max<-(data_merge$Avg_MnMg06.x + data_merge$Avg_MnMg06.y)
data_merge


# Create new data frame without Neolithic samples
newdata<-subset(data, Decade!="Neolithic")
newdata$Salinity_classf<-as.factor(newdata$New_Salinity_class)
str(newdata)
# Rename all levels for Salinity_classf
levels(newdata$Salinity_classf) <- c("Salinity < 6 PSU","Salinity 6 - 12 PSU","Salinity 12 - 16 PSU","Salinity > 16 PSU")
newdata <- newdata[!is.na(newdata$Salinity_classf), ] 
list(newdata)

data_all<-merge(newdata, data_merge, by='Salinity_classf')
str(data_all)

#*****************************
# Plot 2 Salinity classes    
#*****************************
gg<-ggplot(data_all, aes(y = Avg_MnMg06, x = Decades,fill=Salinity_classf)) + 
  geom_boxplot(outlier.shape=NA)+
  geom_ribbon(aes(ymin =sd_min, ymax = sd_max), fill = "red",alpha=0.1,group = "supp")+
  facet_grid(.~ Salinity_classf)+ theme(axis.text.x = element_text(angle = 90))+
    ylab("Mean Mn:Mg ")+ggtitle("")+ theme_classic(base_size = 15) +
  geom_hline(data = grand.means, aes(yintercept = Avg_MnMg06),color=("red"),size=1.1,
             linetype = 2,
             group = "supp")

gg+ scale_fill_manual(values=c("lightskyblue2","mediumturquoise","palegreen3","dodgerblue"),
                      name="Salinity groups",
                      breaks=c("Salinity < 6 PSU","Salinity 6 - 12 PSU","Salinity 12 - 16 PSU","Salinity > 16 PSU"),
                      labels=c("1  Mean Sr:Ca < 3 ", "2  Mean Sr:Ca 3 - 4", "3  Mean Sr:Ca 4 - 4.6","4  Mean Sr:Ca > 4.6"))+
  theme(axis.text.x = element_text(angle = 90))+ylab("Annual mean Mn:Mg") + xlab("Decade")+
  theme(plot.title = element_text(hjust = 0.5))+ylim(0,1.7)+
  theme(legend.title=element_text(size=12),legend.text = element_text(size=15),legend.position="bottom", legend.box = "horizontal")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+ theme(text = element_text(size = 20))

ggsave(path = "figs", filename = "Salinity_groups_plot_2.png", dpi = 600, width = 10, height = 6)

######################################################################################
#
# Statistic test to compare salinity groups and decades
#
######################################################################################
# Empty dataframe
rm(list=ls())
options("install.lock"=FALSE) # Run in case packages will not install in library

# Read data file
cod<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes
str(cod)# Check data

# Turn salinity group into a factor
cod$Salinity_classf<-as.factor(cod$New_Salinity_class)

str(cod)

# Check distribution of Avg_MnMg06
summary(cod$Avg_MnMg06)
hist(cod$Avg_MnMg06)
qqPlot(cod$Avg_MnMg06)

cod$Mn_Mg_log<-log(cod$Avg_MnMg06) #Log transform

# Check distribution of Avg_MnMg06
hist(cod$Mn_Mg_log)
qqPlot(cod$Mn_Mg_log)

ab<-aggregate(Avg_MnMg06 ~ Salinity_classf+Decade, data = cod, mean, all.x = TRUE)
list(ab)

# Salinity group 1
codx <- subset(cod, Decades!="Neolithic")
codxy <- subset(codx, Decades!="1950s")
cod1 <- subset(codxy, Salinity_classf=="1")

cod1<-cod1[,c("Mn_Mg_log","Decades")]
head(cod1)
str(cod1)
cod1$Decade_f<-as.factor(cod1$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Mn_Mg_log ~ Decades, data = cod1) 

## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod1$Mn_Mg_log, cod1$Decades, p.adjust.method="bonferroni")


# Salinity group 2
cod2 <- subset(cod, Salinity_classf=="2")
cod2$Decade_f<-as.factor(cod2$Decades)

# Kruskal-Wallis rank sum test 
kruskal.test(Mn_Mg_log ~ Decades, data = cod2) 

## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod2$Mn_Mg_log, cod2$Decades, p.adjust.method="bonferroni")

# Salinity group 3
cod3 <- subset(cod, Salinity_classf=="3")
# Kruskal-Wallis rank sum test 
kruskal.test(Mn_Mg_log ~ Decades, data = cod3) 

## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod3$Mn_Mg_log, cod3$Decades, p.adjust.method="bonferroni")

# Salinity group 4
cod4 <- subset(cod, Salinity_classf=="4")
# Kruskal-Wallis rank sum test 
kruskal.test(Mn_Mg_log ~ Decades, data = cod4) 

## Pairwise comparisons using Wilcoxon rank sum test 
pairwise.wilcox.test(cod4$Mn_Mg_log, cod4$Decades, p.adjust.method="bonferroni")

# All these pairwise comparisons are inserted into Supplementary information in the manuscript

############################################################################################################
#
#  Figure 3   Metabolic status (Mg:Ca) over time versus BY15 DO sat % 30-90 m + BY15 Temp 30-90 m Q4
#             Repeated measurements model - Age class 1-8 
############################################################################################################
##
# Empty data frame
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library
# Load libraries
library(car)
library(ggplot2)
library(lme4)
library(report)
library(effects)
library(sjPlot)
library(r2glmm)

#Read in file
dd<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes
env<- read.csv("Data\\Data_environmental_factors.csv", sep = ";") 

# merge two data frames by Year
Tot <- merge(dd,env,by="Cal_Year", all.x = TRUE)
list(Tot)

str(Tot)
dg <- subset(Tot, Decades !="Neolithic")
df <- subset(dg, Cal_Year >"1959") 
str(df) 

# Check mean Mg:Ca per age class
#Modern cod
Mg_age_modern<- aggregate(Mg.ESF06 ~ Age, data = df, mean, na.rm=TRUE)
Mg_age_modern

# Neolithic cod
Mg_age_Neo <- subset(Tot, Decades =="Neolithic")
Mg_age_Neolithic<- aggregate(Mg.ESF06 ~ Age, data = Mg_age_Neo, mean, na.rm=TRUE)
Mg_age_Neolithic


# Check distribution of Mg:Ca
summary(df$Mg.ESF06)
hist(df$Mg.ESF06)
qqPlot(df$Mg.ESF06)

# Log transform Mg:Ca
df$Mg_log<-log(df$Mg.ESF06)
hist(df$Mg_log) # Better
ggsave(path = "figs", filename = "Log_MgCa_histogram.png", dpi = 600, width = 6, height = 5)
qqPlot(df$Mg_log)
ggsave(path = "figs", filename = "Log_MgCa_qqplot.png", dpi = 600, width = 6, height = 5)


df$Age_num<-as.numeric(df$Ages)
df$Age_n<-(df$Age_num-1)
summary(df$Age_n)
x_subset <- df[df$Age_n >=1 ,]
summary(x_subset$Age_n)
str(x_subset)

#####################################################################################
#
###  Repeated measurements mixed model for log Mg:Ca temperature and oxygen saturation
#
#####################################################################################
#
str(x_subset, list.len=Inf)
summary(x_subset$T.BY15.Q4..30.90.m.y)
# Best model from comparing models with this being the model with lowst AIC and significant results
g3.mixed<-lmer(Mg_log~T.BY15.Q4..30.90.m.y+O.sat.30.90m.BY15 +Age_n+(Age_n|Fish_ID)+(1|Cal_Year),data=x_subset)
Anova(g3.mixed)
summary(g3.mixed)

# Check multicollinearity
vif(g3.mixed) #All vif < 5 = No multicollinearity

sjPlot::plot_model(g3.mixed)
ggsave(path = "figs", filename = "Log_MgCa_model_sjplot.png", dpi = 600, width = 6, height = 5)

# Effects
plot(allEffects(g3.mixed))
ggsave(path = "figs", filename = "Log_MgCa_model_alleffects.png", dpi = 600, width = 6, height = 5)

r2 <- r2beta(model=g3.mixed,partial=TRUE,method='sgv')
print(r2)
plot(x=r2)
ggsave(path = "figs", filename = "Log_MgCa_model R2.png", dpi = 600, width = 6, height = 5)

# Report model statistics
report(g3.mixed)


#### Plot 3 (Not logged Mg)
dd<-read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes
env<- read.csv("Data\\Data_environmental_factors.csv", sep = ";") 

# merge two data frames by Year
Tot <- merge(dd,env,by="Cal_Year", all.x = TRUE)
list(Tot)

str(Tot)
df <- subset(Tot, Decades !="Neolithic")

df$Age_num<-as.numeric(df$Ages)
df$Age_n<-(df$Age_num-1)
summary(df$Age_n)

x_subset <- df[df$Age_n >=1,]
summary(x_subset$Age_n)
# Remove NA
x_subset <- x_subset[!is.na(x_subset$Age_n), ] 
summary(x_subset$Age_n)

x_subset$Cal_Yearf<-as.factor(x_subset$Cal_Year)
x_subset$Agef<-as.factor(x_subset$Age_n)
summary(x_subset$Age_n)

str(x_subset, list.len=Inf)

x_subset$Cal_Yearn<-unfactor(x_subset$Cal_Yearf)
x_subset$agen<-unfactor(x_subset$Agef)

trend <- subset(x_subset, Cal_Yearn >"1959")
str(trend)

# Neolithic cod
Mg_age_Neo <- subset(Tot, Decades =="Neolithic")
Mg_age_Neo$Age_num<-as.numeric(Mg_age_Neo$Ages)
Neo_subset <- Mg_age_Neo[Mg_age_Neo$Age_n >=1 ,]
summary(Neo_subset$Age_n)

# Mean
summary(Neo_subset$Mg.ESF06) # 0.09409

# Sd
sd(Neo_subset$Mg.ESF06) # 0.02524345

Cal_Yearn <- c(1925:2020)
df <- data.frame(Cal_Yearn)

#Sd min
df2 <- data.frame(Cal_Yearn)
sd_min<-(0.09409-0.02524345)
head(sd_min) # 0.06884655

#Sd max
df3 <- data.frame(Cal_Yearn)
sd_max<-(0.09409+0.02524345)
head(sd_max) # 0.1193334


N1<-df %>%
  mutate(mean = 0.09409) #Mean value for Neolithic time period
head(N1)
N2<-N1 %>%
  mutate(sd_min = 0.06884655) #Mean value for Neolithic time period
head(N2)
N3<-N2 %>%
  mutate(sd_max = 0.1193334) #Mean value for Neolithic time period
head(N3)

data_merge<-merge(x_subset, N3, by='Cal_Yearn')
str(data_merge)

#### PLOT 3
### Repeated measurements mixed model for log Mg:Ca temperature and oxygen saturation
str(x_subset, list.len=Inf)

x_subset$Cal_Yearn<-unfactor(x_subset$Cal_Yearf)
x_subset$agen<-unfactor(x_subset$Agef)

trend <- subset(x_subset, Cal_Yearn >"1959")
str(trend)

# Mean ox BY15 80 m per year
ox<- aggregate(O.sat.30.90m.BY15 ~ Cal_Yearf, data = trend, mean, na.rm=TRUE)
head(ox)

# Moving average 3
ox$rollavgtre <- roll_mean(ox$O.sat.30.90m.BY15, n = 3, align = "center",fill=NA)
head(ox)

ox$Cal_Yearn<-unfactor(ox$Cal_Yearf)

# ox points
trend2 <- subset(x_subset, Cal_Yearn >"1959")
str(trend2)

tre <- subset(x_subset, Age =="3")
# Mean ox BY15 80 m per year
tre_mean<- aggregate(Mg.ESF06 ~ Cal_Yearf, data = tre, mean, na.rm=TRUE)

tre_mean$Cal_Yearn<-unfactor(tre_mean$Cal_Yearf)
head(tre_mean)
# Moving average 3

tre_mean$roll<- roll_mean(tre_mean$Mg.ESF06, n = 3, align = "center",fill=NA)
head(tre_mean)

tre_mean$Cal_Yearn<-unfactor(tre_mean$Cal_Yearf)

# Mean ox BY15 80 m per year
oxy<- aggregate(O.sat.30.90m.BY15 ~ Cal_Yearf, data = trend2, mean, na.rm=TRUE)
head(oxy)

oxy$Cal_Yearn<-unfactor(oxy$Cal_Yearf)
head(oxy)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
 
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
tgc <- summarySE(tre, measurevar="Mg.ESF06", groupvars=c("Cal_Yearn"))
tgc

data_merge<-merge(x_subset, N3, by='Cal_Yearn')
str(data_merge)

ylim.prim <- c(1,2)   # in this example, Age >0 Mg:Ca
ylim.sec <- c(18,420) 
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

p<-ggplot(data_merge,aes(x=Cal_Yearn, y=Mg.ESF06)) + 
  geom_hline(yintercept = 0.0962, color = "red",linetype = "dashed", size=1.1)+
  geom_ribbon(aes(ymin =sd_min, ymax = sd_max), fill = "red",alpha=0.1)+
  geom_jitter(aes(colour=Agef), alpha=0.8, width=0.1,size=2) + 
  scale_color_viridis_d(option = "mako",name = "Age",direction=-1)+
  geom_line(aes(y = roll),data=tre_mean, color="black",size=1.1,linetype="solid")+
  coord_cartesian(ylim=c(0.02,0.18))+
  geom_point(data=tre_mean,aes(x=Cal_Yearn, y=Mg.ESF06)) + 
  geom_errorbar(data=tgc, aes(ymin=Mg.ESF06-se, ymax=Mg.ESF06+se), width=.1) +
  geom_line(data=ox,aes(x=Cal_Yearn, y=a+rollavgtre*b, group=1), color="dodgerblue",size=1.6, stat = "summary", fun="mean")+
  scale_size_manual(values = c(5,1))+ theme_classic()+
  scale_y_continuous("Annual mean Mg:Ca",sec.axis = sec_axis(~(.-a)/b,name=expression("DO saturation (%)")))+ 
  theme(axis.text.x = element_text(angle = 90))+ 
  scale_x_continuous(limits = c(1925, 2020), breaks = seq(1925, 2020, by = 5))+xlab("Year")
  



p + theme(axis.title.y = element_text(family = "sans", size = 15, margin=margin(0,10,0,0),face="plain"),
          axis.title.y.right = element_text(family = "sans", size = 15, margin=margin(0,10,0,0),face="bold"),
          panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title.x = element_text(family = "sans", size = 15, margin=margin(10,0,0,0),face="plain"), 
          axis.text = element_text(family = "sans", size = 13,face="plain"),
          axis.text.y.left = element_text(family = "sans", size = 13, margin=margin(10,0,0,0),face="plain"),
          axis.text.y.right = element_text(family = "sans", size = 13, margin=margin(0,10,0,0),face="plain"),
          plot.title = element_text(hjust = 0.5,family = "sans", size = 15, margin=margin(0,0,10,0),face="bold"),
          panel.background = element_rect(fill = NA))+theme(legend.text = element_text(size = 14))

ggsave(path = "figs", filename = "Plot_3.png", dpi = 600, width = 10, height = 8)
#####################################################################################################
#  *****************************************************************************
#   Figure 4        Predicted mean length at age
#
###*****************************************************************************
### Figure 4 Decades only
####****************************************************************************
# Empty data frame
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library

# Load libraries
library(ggplot2)
library(lme4)
library(nlme)
library(writexl)
library(readxl)
library (pastecs)

#Read in file
data<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes

#*******************************************************************************
#1930
data1930<-subset(data, Decades=="1930s")
str(data1930)

fitlm1930 <- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1930)
data1930$predlm1930 = predict(fitlm1930)
anova(fitlm1930)
summary(fitlm1930)

plot(fitlm1930)
qqnorm(fitlm1930, ~ranef(.))

data1930A <- subset(data1930, Ages=="1") #Subset Age
#CI_t(data1930A$Ages, ci = 0.95)

#Summarize means etc per age
data1930A <- subset(data1930, Ages=="12") #Subset Age

summary(data1930A$predlm1930)


stat.desc(data1930A$predlm1930)
#CI(data1930A$predlm1930, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1930), size = 2,data=data1930,se = FALSE)

#*******************************************************************************
#1940
data1940<-subset(data, Decades=="1940s")
str(data1940)

fitlm1940 <- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1940)
data1940$predlm1940 = predict(fitlm1940)
plot(fitlm1940)
qqnorm(fitlm1940, ~ranef(.))

data1940A <- subset(data1940, Ages=="1") #Subset Age
#CI_t(data1940A$Ages, ci = 0.95)

#Summarize means etc per age
data1940A <- subset(data1940, Ages=="12") #Subset Age

summary(data1940A$predlm1940)


stat.desc(data1940A$predlm1940)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1940), size = 2,data=data1940,se = FALSE)

#*******************************************************************************
#1950
data1950<-subset(data, Decades=="1950s")
str(data1950)

fitlm1950<- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1950)
data1950$predlm1950 = predict(fitlm1950)
summary(fitlm1950)
plot(fitlm1950)
qqnorm(fitlm1950, ~ranef(.))
data1950A <- subset(data1950, Ages=="1") #Subset Age
#CI_t(data1950A$Ages, ci = 0.95)

#Summarize means etc per age
data1950A <- subset(data1950, Ages=="10") #Subset Age

summary(data1950A$predlm1950)
stat.desc(data1950A$predlm1950)
#CI(data1950A$predlm1950, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1950), size = 2,data=data1950,se = FALSE)

#*******************************************************************************
#1960
data1960<-subset(data, Decades=="1960s")
str(data1960)

fitlm1960<- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1960)
data1960$predlm1960 = predict(fitlm1960)

plot(fitlm1960)
qqnorm(fitlm1960, ~ranef(.))

data1960A <- subset(data1960, Ages=="1") #Subset Age
#CI_t(data1960A$Ages, ci = 0.95)
#
#Summarize means etc per age
data1960A <- subset(data1960, Ages=="9") #Subset Age

summary(data1960A$predlm1960)
stat.desc(data1960A$predlm1960)
#CI(data1960A$predlm1960, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1960), size = 2,data=data1960,se = FALSE)

#*******************************************************************************
#1970
data1970<-subset(data, Decades=="1970s")
str(data1970)

fitlm1970<- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1970)
data1970$predlm1970 = predict(fitlm1970)

plot(fitlm1970)
qqnorm(fitlm1970, ~ranef(.))
data1970A <- subset(data1970, Ages=="1") #Subset Age
#CI_t(data1970A$Ages, ci = 0.95)

#Summarize means etc per age
data1970A <- subset(data1970, Ages=="9") #Subset Age

summary(data1970A$predlm1970)
stat.desc(data1970A$predlm1970)
#CI(data1970A$predlm1970, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1970), size = 2,data=data1970,se = FALSE)

#*******************************************************************************

#1980
data1980<-subset(data, Decades=="1980s")
str(data1980)

fitlm1980<- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1980)
data1980$predlm1980 = predict(fitlm1980)

plot(fitlm1980)
qqnorm(fitlm1980, ~ranef(.))

#Summarize means etc per age
data1980A <- subset(data1980, Ages=="1") #Subset Age
#CI_t(data1980A$Ages, ci = 0.95)

data1980A <- subset(data1980, Ages=="14") #Subset Age
summary(data1980A$predlm1980)
stat.desc(data1980A$predlm1980)
#CI(data1980A$predlm1980, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1980), size = 2,data=data1980,se = FALSE)

#*******************************************************************************
#
#1990
data1990<-subset(data, Decades=="1990s")
str(data1990)

fitlm1990 <- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data1990)
data1990$predlm1990 = predict(fitlm1990)
plot(fitlm1990)
qqnorm(fitlm1990, ~ranef(.))

data1990A <- subset(data1990, Ages=="1") #Subset Age
#CI_t(data1990A$Ages, ci = 0.95)

#Summarize means etc per age
data1990A <- subset(data1990, Ages=="8") #Subset Age
summary(data1990A$predlm1990)
stat.desc(data1990A$predlm1990)
#CI(data1990A$predlm1990, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm1990), size = 2,data=data1990,se = FALSE)

#*******************************************************************************

#2000s
data2000<-subset(data, Decades=="2000s")
str(data2000)

fitlm2000 <- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data2000)
data2000$predlm2000 = predict(fitlm2000)
plot(fitlm2000)
qqnorm(fitlm2000, ~ranef(.))
data2000A <- subset(data2000, Ages=="1") #Subset Age
#CI_t(data2000A$Ages, ci = 0.95)

#Summarize means etc per age
data2000A <- subset(data2000, Ages=="10") #Subset Age
summary(data2000A$predlm2000)
stat.desc(data2000A$predlm2000)
#CI(data2000A$predlm2000, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm2000), size = 2,data=data2000,se = FALSE)

#*******************************************************************************
#
#2010s
data2010<-subset(data, Decades=="2010s")
str(data2010)

fitlm2010 <- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = data2010)
data2010$predlm2010 = predict(fitlm2010)
plot(fitlm2010)
qqnorm(fitlm2010, ~ranef(.))

data2010A <- subset(data2010, Ages=="1") #Subset Age
#CI_t(data2010A$Ages, ci = 0.95)

#Summarize means etc per age
data2010A <- subset(data2010, Ages=="10") #Subset Age
summary(data2010A$predlm2010)
stat.desc(data2010A$predlm2010)
#CI(data2010A$predlm2010, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlm2010), size = 2,data=data2010,se = FALSE)

#*******************************************************************************
#Neo
dataNeo<-subset(data, Decades=="Neolithic")
str(dataNeo)

fitlmNeo <- lme(New.L_age ~ Ages,random=~1|Fish_ID, data = dataNeo)
dataNeo$predlmNeo = predict(fitlmNeo)

plot(fitlmNeo)
qqnorm(fitlmNeo, ~ranef(.))

dataNeoA <- subset(dataNeo, Ages=="1") #Subset Age
#CI_t(dataNeoA$Ages, ci = 0.95)

dataNeo$group<-paste(dataNeo$predlmNeo,dataNeo$Ages)

dataNeo$groupf<-as.factor(dataNeo$group)
head(dataNeo)
table(dataNeo$predlmNeo,dataNeo$Ages)

#Summarize means etc per age
dataNeo1 <- subset(dataNeo, Ages=="2") #Subset Age
summary(dataNeo1$predlmNeo)
stat.desc(dataNeo1$predlmNeo)
#CI(dataNeo1$predlmNeo, ci=0.95)

ggplot(data,aes(x = Ages, y = New.L_age, colour=Decades))+ylab("Total length (mm)") + 
  xlab("Age (years)")+scale_x_continuous(breaks = seq(0, 13, 1))+
  geom_smooth(aes(y = predlmNeo), size = 2,data=dataNeo,se = FALSE)

######################################################################################################
### *****************************************************************************  
###   Oxygen saturation data
###
####################################################################################################
###
env<- read.csv("Data\\Data_environmental_factors.csv", sep = ";") 
str(env)
env$Decades<-as.factor(env$Decade_ox)
str(env)
ox_mean<-aggregate(O.sat.30.90m.BY15 ~ Decades, data = env, mean, na.rm=TRUE)
ox_mean
###########################################################################################
###########  Decades mean oxygen Neolithic
###########################################################################################

dataNeo$Decades<-as.factor(dataNeo$Decades)
df_Neo<-dataNeo[c("predlmNeo", "Ages","Decades")]
str(df_Neo)

data_merge_Neo <- merge(df_Neo, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_Neo)

write_xlsx(data_merge_Neo, "Data\\Neo_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 1930s
###########################################################################################
data1930$Decades<-as.factor(data1930$Decades)
df_1930<-data1930[c("predlm1930", "Ages","Decades")]
str(df_1930)

data_merge_1930 <- merge(df_1930, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1930)

write_xlsx(data_merge_1930, "Data\\1930_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 1940s
###########################################################################################
data1940$Decades<-as.factor(data1940$Decades)
df_1940<-data1940[c("predlm1940", "Ages","Decades")]
str(df_1940)

data_merge_1940 <- merge(df_1940, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1940)

write_xlsx(data_merge_1940, "Data\\1940_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 1950s
###########################################################################################
data1950$Decades<-as.factor(data1950$Decades)
df_1950<-data1950[c("predlm1950", "Ages","Decades")]
str(df_1950)

data_merge_1950 <- merge(df_1950, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1950)

write_xlsx(data_merge_1950, "Data\\1950_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 1960s
###########################################################################################
data1960$Decades<-as.factor(data1960$Decades)
df_1960<-data1960[c("predlm1960", "Ages","Decades")]
str(df_1960)

data_merge_1960 <- merge(df_1960, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1960)

write_xlsx(data_merge_1960, "Data\\1960_ox.xlsx")

###########################################################################################
###########  Decades mean oxygem 1970s
###########################################################################################

data1970$Decades<-as.factor(data1970$Decades)
df_1970<-data1970[c("predlm1970", "Ages","Decades")]
str(df_1970)

data_merge_1970 <- merge(df_1970, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1970)

write_xlsx(data_merge_1970, "Data\\1970_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 1980s
###########################################################################################
data1980$Decades<-as.factor(data1980$Decades)
df_1980<-data1980[c("predlm1980", "Ages","Decades")]
str(df_1980)

data_merge_1980 <- merge(df_1980, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1980)

write_xlsx(data_merge_1980, "Data\\1980_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 1990s
###########################################################################################
data1990$Decades<-as.factor(data1990$Decades)
df_1990<-data1990[c("predlm1990", "Ages","Decades")]
str(df_1990)

data_merge_1990 <- merge(df_1990, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_1990)

write_xlsx(data_merge_1990, "Data\\1990_ox.xlsx")

###########################################################################################
###########  Decades mean oxygen 2000s
###########################################################################################
data2000$Decades<-as.factor(data2000$Decades)
df_2000<-data2000[c("predlm2000", "Ages","Decades")]
str(df_2000)

data_merge_2000 <- merge(df_2000, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_2000)

write_xlsx(data_merge_2000, "Data\\2000_ox.xlsx")

##############################################################################################
###########  Decades mean oxygen 2010s
##############################################################################################

data2010$Decades<-as.factor(data2010$Decades)
df_2010<-data2010[c("predlm2010", "Ages","Decades")]
str(df_2010)

data_merge_2010 <- merge(df_2010, ox_mean, by = c("Decades"),all.x = TRUE)
head(data_merge_2010)

write_xlsx(data_merge_2010, "Data\\2010_ox.xlsx")

##############################################################################################
###    Decades mean oxygen  merged to file "ox_growth.xlsx"
##############################################################################################

Data_ox_growth <- read_excel("Data\\ox_growth.xlsx")
str(Data_ox_growth)
Data_ox_growth$Decades<-as.factor(Data_ox_growth$Decades)
Data_ox_growth$group<-as.factor(Data_ox_growth$group)
Data_ox_growth$Ages<-as.factor(Data_ox_growth$Ages)

str(Data_ox_growth)
summary(Data_ox_growth$O.sat.30.90m.BY15)

Data_ox_growth_sub <- subset(Data_ox_growth, Decades !="1940s")
Data_ox_growth_1940 <- subset(Data_ox_growth, Decades =="1940s")
Data_ox_growth_sub2 <- subset(Data_ox_growth_sub, Decades !="Neolithic")
Data_ox_growth_Neo <- subset(Data_ox_growth, Decades =="Neolithic")

p<-ggplot(Data_ox_growth_sub2, aes(x = Ages, y = value, colour =O.sat.30.90m.BY15, group = group))+
  geom_smooth(size=2,se = FALSE)+
  scale_colour_gradient2(low = "red", mid = "yellow", high = "dodgerblue", midpoint =70)+
  guides(col = guide_colourbar(title = "Dissolved oxygen (%)"))+
  geom_smooth(aes(y = value, group = group), Data_ox_growth_1940,colour="black", size = 2,se = FALSE,linetype="dotted")+
  geom_smooth(aes(y = value, group = group), Data_ox_growth_Neo, colour="black", size = 2,se = FALSE,linetype="dashed")+
  ylab("Total length (mm)") + 
  xlab("Age (year)")

print(p)+theme_classic(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = c(.1, .9),
        legend.justification = c("left", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))

ggsave(path = "figs", filename = "Decades_oxygen_Plot_4.png", dpi = 600, width = 10, height = 8)
###########################################################################################################
##############################################################################################

#
#  Plot 5        Environmental parameters 
#
#
#  ***************************************************************************** 

# Empty data frame
rm(list = ls())
options("install.lock"=FALSE) # Run in case packages will not install in library
# Load libraries
library(ggplot2)
library(RcppRoll)


#Read in file
dd<- read.csv("Data\\Data_otolithchemistry.csv", sep = ";") # All age classes
env<- read.csv("Data\\Data_environmental_factors.csv", sep = ";") 

# merge two data frames by Year
Tot <- merge(dd,env,by="Cal_Year", all.x = TRUE)
str(Tot, list.len=ncol(Tot))

df <- subset(Tot, Cal_Year >"1959")

# Moving average 3
df$roll_ox <- roll_mean(df$O.sat.30.90m.BY15, n = 3, align = "center", fill = NA)
df$roll_anox <- roll_mean(df$Hypoxic.volume.0.y, n = 3, align = "center", fill = NA)
df$roll_hypox <- roll_mean(df$Hypoxic.volume.2.y, n = 3, align = "center", fill = NA)


ylim.prim <- c(0,1)   # secondary y-axis. In this example, hypoxic and anoxic volume
ylim.sec <- c(-30,2) 
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

#Plot oxygen saturation (BY15 30-90 m depth), hypoxic volume and anoxic volume over time
p<-ggplot(df,aes(x=Cal_Year, y=O.sat.30.90m.BY15)) + 
  geom_line(aes(x=Cal_Year, y=roll_ox, group=1), color="dodgerblue",size=1.5, stat = "summary", fun="mean")+
  geom_line(aes(x=Cal_Year, y=a+roll_anox*b, group=1), color="black",size=1.5, stat = "summary", fun="mean")+
  geom_line(aes(x=Cal_Year, y=a+roll_hypox*b, group=1), color="gray",size=1.5, stat = "summary", fun="mean")+
  scale_size_manual(values = c(5,1))+ 
  scale_y_continuous("Oxygen saturation (%)",sec.axis = sec_axis(~(.-a)/b,name=expression(Hypoxic~"/"~anoxic~volume~km^3)))+
  theme(axis.text.x = element_text(angle = 90))+ scale_x_discrete(breaks=seq(1960,2020,by=5))+xlab("Year")

p + theme(axis.title.y = element_text(family = "sans", size = 14, margin=margin(0,10,0,0),face="plain"),
          axis.title.y.right = element_text(family = "sans", size = 14, margin=margin(0,10,0,0),face="bold"),
          panel.grid.major = element_line(linetype = "blank"), 
          panel.grid.minor = element_line(linetype = "blank"), 
          axis.title.x = element_text(family = "sans", size = 14, margin=margin(10,0,0,0),face="plain"), 
          axis.text = element_text(family = "sans", size = 12,face="plain"),
          axis.text.y.left = element_text(family = "sans", size = 12, margin=margin(10,0,0,0),face="plain"),
          axis.text.y.right = element_text(family = "sans", size = 12, margin=margin(0,10,0,0),face="plain"),
          plot.title = element_text(hjust = 0.5,family = "sans", size = 14, margin=margin(0,0,10,0),face="bold"),
          panel.background = element_rect(fill = NA))+
  theme(legend.position="bottom")+
  theme(legend.key.size = unit(4, 'mm'),legend.key.width = unit(2.5, 'cm'))

ggsave(path = "figs", filename = "Environment_panel_1_Plot_5.png", dpi = 600, width = 10, height = 8)


##################

