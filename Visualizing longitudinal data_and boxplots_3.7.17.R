library(foreign)#import from SAS
library(stats)
library(epicalc)
library(Hmisc)
library(survey)
library (plyr) #use to revalue variables
library(reshape)#use to rename variables
library (ggplot2)#use to plot
library (psych)
library (ICC)
library(mgcv)#use for GAMMs
library(nlme)#use for LMMs

#########################################
#       Data Import and Reshape 
#########################################
#Changes all variable names in wide dataset to lower case
names(aim2_wide)<-tolower(names(aim2_wide))
names(aim2_wide) #check that names weren't truncated

#Check import of data
table(aim2_wide$preterm)#check to make sure numbers of cases and controls correct (116 cases, 323 controls)
aim2_wide[1:5,1:8]

#Create dataset with only thyroid hormones and key covariates
thyroid_wide<-aim2_wide[,c("id","ft4_1","ft4_2","ft4_3","ft4_4","tsh_1","tsh_2","tsh_3","tsh_4", 
                   "t3_1","t3_2","t3_3","t3_4","t4_1","t4_2","t4_3","t4_4","gat1","gat2","gat3","gat4",
                   "age_cat2","race_cat2","bmi_cat2","edu_cat2", "insur2","preterm")]

#Rename gestational age for reshaping
which(colnames(thyroid_wide)==c("gat1"))#determine the column numbers of gat variables
names(thyroid_wide[,c(18:21)])#check to make sure these are gat1-gat4
names(thyroid_wide)[18:21] <- c("gat_1","gat_2","gat_3","gat_4") #rename columns
names(thyroid_wide[,c(18:21)])#check to make sure column names changed

#Reshape from wide to long (with each hormone and gestational age as its own column)
thyroid_long1<-reshape(thyroid_wide, 
               varying=c("ft4_1","ft4_2","ft4_3","ft4_4","tsh_1","tsh_2","tsh_3","tsh_4",
                         "t3_1","t3_2","t3_3","t3_4","t4_1","t4_2","t4_3","t4_4",
                         "gat_1","gat_2","gat_3","gat_4"), #specific variables that are time-varying
               direction="long", #specific direction of long data set
               idvar="id", #keep as identifier
               sep="_", #time_varying variables separated by underscore so will drop
               timevar=c("visit"))#will had visit as time variable variable (1:4)

#Order dataset by id and visit (ascending order)
thyroid_long.order<-thyroid_long1[order(thyroid_long1$id,thyroid_long1$visit),]

#Check to make sure dataset is in correct order
thyroid_long.order[thyroid_long.order$id==1,]#compare new long dataset with wide to see if #'s match
thyroid_wide[thyroid_wide$id==1,]

#########################################
#       Create Spaghetti Plot
#########################################

#Spaghetti plot used to visualize variablility in hormones between subjects 
#Also used to visualize variance within a subject over time

#Create complete dataset for spaghetti
thyroid.plot<-na.omit(thyroid_long1[,c("tsh","gat","preterm","id","visit")])#need to use first long dataset
attach(thyroid.plot)

summary(thyroid.plot$tsh)

#Relabel preterm variable for facet
thyroid.plot$preterm<-factor(thyroid.plot$preterm, levels=c(0,1), labels=c("Term","Preterm"))
levels(thyroid.plot$preterm)#check to make sure relabel worked

##ggplot for spaghetti plot
plot.1<-ggplot(data=thyroid.plot, aes(gat, tsh, group=id))+ #introduce main plot
  geom_point(size=0.2)+ #add points for each measurement
  geom_line(color=id, size=0.2)+ #create lines for each person
  facet_grid(. ~ preterm)+ #facet by preterm birth status
  scale_x_continuous(name="Gestational Age (weeks)") + 
  scale_y_continuous(limits=c(0,6), name="Concentration (ng/dL)") +
  theme_bw()+
  ggtitle("Spaghetti Plot Showing Individual Change In TSH Across Pregnancy")+
  theme(plot.title=element_text(face="bold",size=14, family="sans", hjust = 0.5), 
        axis.text.x = element_text(angle = 0, hjust = 1,vjust=1, size=10,family="sans"), 
        axis.text.y=element_text(vjust=1, size=10,family="sans"),
        text=element_text(size=14,family="sans"), 
        strip.text.x=element_text(size=14, face="bold",family="sans"),
        legend.title=element_blank())
plot.1

#########################################
# Visualize general trend over time
#########################################

#Same code as spaghetti plot but drop geom_line and add a stat_smooth term and look at general trend
plot.2<-ggplot(data = thyroid.plot, aes(x = gat, y = tsh, group = id))+
  geom_point(aes(color=factor(visit)))+
  stat_smooth(aes(group=1))+#need to aesthetic to fit over all day (instead of by group= id)
  facet_grid(. ~ preterm)+
  scale_color_brewer(palette="Dark2", direction=-1)+
  labs(title = "TSH Levels in Pregnancy", x = "Gestational Age (Weeks)", 
       y = "Concentration (uIU/mL)", color = "Study Visit") +
  scale_y_log10(limits=c(0.1,10))+ #add limits to zoom into trends (get rid of limits to see outliers)
  theme_bw()+
  theme(plot.title=element_text(face="bold",size=14, family="sans", hjust = 0.5), 
        axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size=10,family="sans"), 
        axis.text.y=element_text(vjust=0.5, size=10,family="sans"),
        text=element_text(size=14,family="sans"), 
        strip.text.x=element_text(size=14, face="bold",family="sans"))
plot.2

########################################################
# Visualize hormone distributions by Age Group: Boxplots
########################################################

#Reshape long dataset again so that concentration in one column and hormone label in another
thyroid_long2<-reshape(thyroid_long1,
                       varying=c("ft4","tsh","t3","t4"),
                       v.name="conc",
                       direction="long",
                       idvar=c("id", "visit"), #keep as an identifier bc already in long format
                       timevar="hormone",
                       times=c("ft4","tsh","t3","t4")) 
thyroid_long2.order<-thyroid_long2[order(thyroid_long2$id,thyroid_long2$visit),]

#Check to make sure dataset is correctly formed
thyroid_long2.order[thyroid_long2.order$id==1 | thyroid_long2.order$id==4,]
thyroid_long1[thyroid_long1$id==1 | thyroid_long1$id==4,]

###calculate quantiles and mean for each hormone by each age group
quantiles2 <- ddply(thyroid_long2, .(age_cat2, hormone), summarize,
                   lower=quantile(conc, 0.25, na.rm=T),
                   middle=quantile(conc, 0.50, na.rm=T),
                   upper=quantile(conc, 0.75, na.rm=T),
                   ymin=quantile(conc, 0.05, na.rm=T), 
                   ymax=quantile(conc, 0.95, na.rm=T),
                   mean=mean(conc, na.rm=T))

##Rename hormones for facet labels
quantiles2$hormone<-as.factor(quantiles2$hormone)
levels(quantiles2$hormone)

quantiles2$hormone<- factor(quantiles2$hormone, levels=c("ft4","t3","t4","tsh"), 
             labels=c("FT4","T3","T4","TSH"))
levels(quantiles2$hormone)

##revalue age for x-axis labels
quantiles2$age_cat2<-as.numeric(quantiles2$age_cat2)
quantiles2$age_cat2<- factor(quantiles2$age_cat2, levels=c("1","2","3","4"), 
                            labels=c("18-24[ref]","25-29","30-34","35+"))
levels(quantiles2$age_cat2)

#Change order of Hormones for facets
library("plyr")
neworder <- c("TSH","FT4","T3","T4")
hormone.order<- arrange(transform(quantiles2,hormone=factor(hormone,levels=neworder)),hormone)
hormone.order$hormone2<-factor(hormone.order$hormone, levels=c("TSH","FT4","T3","T4"),
                               labels = c("TSH (uIU/mL)", "FT4 (ng/dL)", "T3 (ng/dL)", "T4 (ug/dL)"))

#Create data frame to add astricks for significance levels (LMMs with random intercept for ID and slope for GA to test differences)
sig <- data.frame(hormone2 = c(rep("T3 (ng/dL)",3), rep("T4 (ug/dL)",2)), p = paste("*"), 
                  x=c(2, 3, 4, 3, 4), y=c(2.28, 2.4, 2.32, 13.5, 13.3))
##Boxplot
hormone_plot.age <- ggplot(hormone.order, aes(x=factor(age_cat2))) +
         geom_boxplot(aes(lower=lower, middle=middle, upper=upper,
                          ymin=ymin, ymax=ymax, width = 0.75, fill=hormone),
                      color="black", stat="identity", outlier.shape = NA) + 
  scale_fill_brewer(palette="Blues", guide=FALSE)+
  theme_bw()+
  facet_wrap(~hormone2, scales="free")+
  ylab("Concentration")+
  xlab("Maternal Age at Enrollment (years)")+
  ggtitle("Thyroid Hormone Concentrations by Age Group")+
  geom_text(data=sig, aes(x, y, label=p), size=8, colour="black",family="sans")+
  geom_point(data=hormone.order, aes(y=mean), 
             color="red", shape=20, size=3)+
  theme(plot.title=element_text(face="bold",size=14, family="sans", hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=1, size=10,family="sans"), 
        axis.text.y=element_text(vjust=0.5,, size=10,family="sans"),
        text=element_text(size=14,family="sans"), 
        strip.text.x=element_text(size=14, face="bold",family="sans"))
homorone_plot.age
  
  
  
  
  