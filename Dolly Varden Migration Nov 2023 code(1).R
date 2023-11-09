
library(ggplot2)
library(grid)
require(mgcv)
require(lme4)
require(lmerTest)
require(tidyverse)
require(tidyr)
require(plyr)

#define the look of plots for use later
MSTheme <-  theme_bw(16)+  
  theme(
    legend.title=element_blank(),
    legend.position="right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key=element_blank())

#set working directory
setwd("C:/Data/Documents/Graduate Committees/Fuller")
#read in data file
run.timing.lm <- read.csv("DailyFishCountsLinearModel_final.csv")

#select the days of year (DOY) that are consistently available 
sel<-run.timing.lm[run.timing.lm$DOY>176 & run.timing.lm$DOY<217,]

#remove rows with NAs in selected columns
#Just considering the counts of Dolly Varden, Chum salmon, sockeye salmon
sel2<-na.omit(sel[,c(1:3,5,7,9)])

#ensure that variables are stored as the correct data type
sel2$Dolly<-as.numeric(sel2$Dolly)
sel2$Sockeye<-as.numeric(sel2$Sockeye)
sel2$Chum<-as.numeric(sel2$Chum)
sel2$year<-as.factor(sel2$year)
sel2$count_date<-as.Date(sel2$count_date,format="%Y-%m-%d")


#Set of three possible models
# random effect is "(Sockeye|year)" the parentheses around it denotes this
#is a random effect and including "Sockeye" before the line tells the
#model that you want to include a random slope effect for the Sockeye salmon coefficient

M.null<-lmer(Dolly~1 +(1|year),data=sel2,REML = T)
M1<-lmer(Dolly~scale(Sockeye) +(scale(Sockeye)|year),data=sel2,REML=T)
M2<-lmer(Dolly~scale(Chum) +(scale(Chum)|year),data=sel2,REML=T)

#read outputs for models
summary(M.null)
summary(M1)
summary(M2)

#Test Akaike's information criteria for model fit
AIC(M.null, M1, M2) #this output provided the degrees of freedom (df)
                    #dfs appear in Table 1
require(AICcmodavg)
#AICc now used because of small sample sizes
AICc(M.null)
AICc(M1)
AICc(M2)
#M1 is the best model so we need to calculated changes in AIC values compared to that model
#These values are reported as the Delta AICc in Table 1
AICc(M2)-AICc(M1)
AICc(M.null)-AICc(M1)

#refit models with maximum likelihood (not REML)
#these models are used following model comparisons above
M1<-lmer(Dolly~scale(Sockeye) +(scale(Sockeye)|year),data=sel2,REML=F)
M2<-lmer(Dolly~scale(Chum) +(scale(Chum)|year),data=sel2,REML=F)
#view model summary for the information on t-values and
#coefficients that appear in the results section of this paper
summary(M1)
summary (M2)


### function to calculate the r2 values of mixed effect models
#Provided in Table 1
require(piecewiseSEM)
rsquared.glmm=function(modlist) {
  do.call(rbind,lapply(modlist,function(i) {
    if(inherits(i,"merMod") | class(i)=="merLmerTest") {
      VarF=var(as.vector(fixef(i) %*% t(i@pp$X)))
      VarRand=colSums(do.call(rbind,lapply(VarCorr(i),function(j) j[1])))
      VarResid=attr(VarCorr(i),"sc")^2
      Rm=VarF/(VarF+VarRand+VarResid)
      Rc=(VarF+VarRand)/(VarF+VarRand+VarResid)
      Rsquared.mat=data.frame(Class=class(i),Marginal=Rm,Conditional=Rc,
                              AIC=AIC(update(i,REML=F))) }
    else { print("Function requires models of class lm, lme, mer, or
merMod")
    } } ) ) }

#list of models provided in the order shown in Table 1
models<-list(M1,M2,M.null)
#output does not include model names, but output order matches previous line of code
rsquared.glmm(models)


#calculate the r values for annual specific correlations 
#shown in Figure 2
annual.cor<-as.data.frame(sel2 %>%
                            group_by(year) %>%
                            summarize(sockeye.cor=cor(Sockeye, Dolly),
                                      chum.cor=cor(Chum, Dolly),
                                      salmon.cor=cor(Sockeye,Chum)))
annual.cor

#calculate the p values for each correlation to see which are significant 
#determines the significance stars in Figure 2
#Dolly Varden and sockeye salmon
require(forestmangr)
detach(package:plyr)
require(tidyverse)
require(broom)
#Dolly Varden and sockeye salmon models for each year
by_year <- group_by(sel2, year)
sockeye.cor<-as.data.frame(do(by_year, 
    glance( 
    lm(scale(Dolly) ~ scale(Sockeye), data = .))))
sockeye.cor

#Dolly Varden and chum salmon models for each year
group_by(sel2, year)
chum.cor<-as.data.frame(do(by_year, 
                              glance( 
                                lm(scale(Dolly) ~ scale(Chum), data = .))))
chum.cor




#define colors to use on plots
colors <- c("Dolly Varden" = "black", "Sockeye Salmon" = "#F8766D", "Chum Salmon" = "#619CFF")

#plot of fish counts using scaled values
#Figure 2
ggplot(sel2)+
  geom_line(aes(x=DOY,y=scale(Dolly),col="Dolly Varden"),size=1.25)+
  geom_line(aes(x=DOY,y=scale(Sockeye),col="Sockeye Salmon"))+
  geom_line(aes(x=DOY,y=scale(Chum),col="Chum Salmon"))+
  facet_wrap(~year)+
  scale_color_manual(values=colors)+
  scale_y_continuous("Scaled Fish Counts")+
  labs(color="Legend")+
  MSTheme+
  theme(legend.position = "bottom")
ggsave("DailyFishCounts.pdf",width=8, height = 9)


#plot of fish counts using absolute values
#Figure not included in manuscript
ggplot(sel2)+
  geom_line(aes(x=DOY,y=Dolly,col="Dolly Varden"),size=1.25)+
  geom_line(aes(x=DOY,y=Sockeye,col="Sockeye Salmon"))+
  geom_line(aes(x=DOY,y=Chum,col="Chum Salmon"))+
    facet_wrap(~year)+
  scale_color_manual(values=colors)+
  scale_y_continuous("Fish Counts")+
  labs(color="Legend")+
  MSTheme+
  theme(legend.position = "bottom")


####check residuals
#extract predictions and residuals from both fitted models
sel2$fit.M1<-predict(M1) #add model fits to the dataframe
sel2$resid.M1<-resid(M1, type = "pearson")
sel2$fit.M2<-predict(M2) #add model fits to the dataframe
sel2$resid.M2<-resid(M2, type = "pearson")

#residual plots by year
#Sockeye salmon model
ggplot(sel2, aes(x=fit.M1,y=resid.M1, group=year,col=year))+
  facet_wrap(~year,nrow=3)+
  geom_point()+
  geom_hline(yintercept =0, col="red")
#chum salmon model
ggplot(sel2, aes(x=fit.M2,y=resid.M2, group=year,col=year))+
  facet_wrap(~year,nrow=3)+
  geom_point()+
  geom_hline(yintercept =0, col="red")




