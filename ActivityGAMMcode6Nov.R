#functions
se<-function(x) { sqrt((var(na.omit(x))/length(na.omit(x))))}   #calculates std error
length.noNA<-function(x){length(x[which(is.na(x)!=T)])}

#packages
library(lubridate)
library(tidyverse)
library(readxl)
library(dplyr)
library (mgcv)
library (nlme)
library(ggplot2)
library(MuMIn)

##read in the data
#raw data
a<-read.csv("/Users/gracevaziri/OneDrive/Grace Thesis/Analyses/GAMMs/WA Temp/All birds_best antenna_with temp.csv")
mapping_file=read_tsv("16S_mapping_20180518.txt")                                                          #has condition data
field_data = read_excel("~/OneDrive/Box Sync/WA field season/Data entry/Daily data entry_ALL_DATA.xlsx",   #load field data
                        sheet = "recap")

##create a column with the date and time that the receiver started collecting data called "data_date_time
##create another column with the date and time that the transmitter was actually placed on the bird called "on_date_time"
a$date_time<-paste(a$Date,a$WA.time)
a$on_date_time<-paste(a$Transmitter.date,a$Transmitter_on)

##convert with lubridate 
#lubridate compatible date/time of observations
a$date_time<-mdy_hms(a$date_time, tz = "US/Pacific")
a$date_time<-round_date(a$date_time,unit = "minute")

#lubridate compatible date/time of start of data collection
a$on_date_time<-mdy_hms(a$on_date_time, tz = "US/Pacific")
a$on_date_time<-round_date(a$on_date_time, unit = "minute")

#cut this down so that I only have what I need
b <- a %>%
  dplyr::select(bird_id,comb_trt, date_time, on_date_time, Power)

#create a matrix with average power per minute per bird
powmin<-tapply(b$Power,list(b$bird_id,b$date_time),mean,na.rm=T)
dim(powmin)

# taking row names and repeating 29282 times
bird_id<-rep(dimnames(powmin)[[1]], 29282)  
length(bird_id)

# taking col names repeating each one the same # of time as there were rows in powmin
time<-rep((dimnames(powmin)[[2]]),each=nrow(powmin))        

#make powmin a vector
powmin<-as.vector(powmin)

#make a datafram with powmin vector, bird_id, and time
df2<-data.frame(bird_id,time,powmin)

#make a rounded time column that goes to the half hour
df2$time<-as.POSIXct(df2$time)
df2$rounded_time<-round_date(df2$time, unit = "30 mins")
df2$rounded_time<-as.factor(df2$rounded_time)
df2$rounded_time<-as.character(df2$rounded_time)
df2$rounded_time<-as.POSIXct(df2$rounded_time)

#Match the start time with the bird
df2$inj_time<-b$on_date_time[match(df2$bird_id,b$bird_id)]


#make a rounded injection time column that goes to the half hour
df2$inj_rounded_time<-round_date(df2$inj_time, unit = "30 mins")
df2$inj_rounded_time<-as.factor(df2$inj_rounded_time)
df2$inj_rounded_time<-as.character(df2$inj_rounded_time)
df2$inj_rounded_time<-as.POSIXct(df2$inj_rounded_time)

#get rid of rows with NA in powmin
df2<-na.omit(df2)
length(df2$bird_id)

#make some new columns to denote whether an activity reading is preceded by a reading that occurred exactly one minute before
#find the interval of time in minutes between when a bird was treated and when a measurement was recorded
str(df2)
df2$inj_time<-as.factor(df2$inj_time)
df2$inj_time_POSIXct<-as.character(df2$inj_time)
df2$inj_time_POSIXct<-as.POSIXct.default(df2$inj_time_POSIXct)
df2$time_span_min<-interval(df2$inj_time_POSIXct,df2$time)


#find the number of minutes in the interval
df2$min_since_inj<-df2$time_span_min/ dminutes(x=1)


#find the interval of time in hours between when a bird was treated and when a measurement was recorded
df2$time_span_hours<-interval(df2$inj_rounded_time,df2$rounded_time)


#find the number of hours in the interval
df2$hours_since_inj<-df2$time_span_hours/ dhours(x=1)
range(df2$hours_since_inj)


#df2$condition <- mapping_file$Condition[match(df2$bird_id,mapping_file$phinchID)]

#add the number of days it's been since a bird was treated with anthelmintic or water
df2$days_since_treat <- field_data$Days_since_treat[match(df2$bird_id,field_data$bird_id)]
dplyr::filter(df2,hours_since_inj< 0)

#show columns I want
df2<-df2 %>%
  dplyr::select(bird_id, time, inj_time, powmin, min_since_inj, hours_since_inj,days_since_treat)

#show the number of minutes that have elapsed since the previous reading, by bird
df2$elapsed_time<- ave(df2$min_since_inj, df2$bird_id,FUN = function(x) c(0, diff(x)))

#make a column detnoting whether the number of minutes elapsed beween one minute and the next was one or not one
df2$elapsed_min_1<-ifelse(df2$elapsed_time==1, 1,NA)
head(df2)
df2[100:200,]
unique(df2$elapsed_min_1)

#make a column that shows the difference in power from one minute to the next
df2$power_diff<- ave(df2$powmin, df2$bird_id, FUN = function(x) c(0, diff(x)))
head(df2)

#make a column that evaluates whether the previous observation occured one minute ago, and if so,
#whether the bird was active between a minute and the previous minute
df2$active<-ifelse(is.na(df2$elapsed_min_1 == T), NA,
                  ifelse(df2$power_diff >= 4 |df2$power_diff <= -4, 1,0))


#arrange by bird to make it easier to check
df2<-arrange(df2,bird_id)
head(df2)

#look at who's causing the weirdness(having so many datapoints for weird hours)
df2[df2$hours_since_inj >50,]

#remove the rows of data for more than 48hrs
df2<- df2 %>%
  filter(!hours_since_inj > 48)
head(df2)

#look at the count of minutes recorded for each hour each bird
obs_count<-tapply(df2$active, list(df2$bird_id,df2$hours_since_inj), length.noNA )


#look at the sum of minutes active for each hour each bird
obs_active_sum<-tapply(df2$active, list(df2$bird_id,df2$hours_since_inj), sum, na.rm=T)


# taking row names and repeating 96 times
bird_id<-rep(dimnames(obs_active_sum)[[1]],96)  
length(bird_id)

# taking col names repeating each one the same # of time as there were rows in powmin
hour<-rep((dimnames(obs_active_sum)[[2]]),each=nrow(obs_active_sum))        
length(hour)

#make obs_active_sum a vector
obs_active_sum<-as.vector(obs_active_sum)
head(obs_active_sum)
length(obs_active_sum)

#make obs_count a vector
obs_count<-as.vector(obs_count)
length(obs_count)

#make a dataframe with the bird, the hour, the number of obs and the sum of activity
df3<-data.frame(bird_id, hour, obs_count, obs_active_sum)

#convert hours from factor into numeric
df3$hour<-as.character(df3$hour)
df3$hour<-as.numeric(df3$hour)


#reorder dataframe by bird_id
df3<-arrange(df3,bird_id)
head(df3,100)

#get rid of rows with NA in obs_count
df3<-na.omit(df3)
length(df3$bird_id)

#only keep rows with more than 10 observations for each half hour
df3<-df3 %>%
  filter(!obs_count < 11)
head(df3)
length(df3$bird_id)

#create a column for the proportion of an hour spent active
df3$propact<- df3$obs_active_sum/df3$obs_count
head(df3)

#Match the start time with the bird
head(b)
df3$inj_time<-b$on_date_time[match(df3$bird_id,b$bird_id)]
head(df3)

#add in treatment information for each bird
df3$trt<-b$comb_trt[match(df3$bird_id,b$bird_id)]
head(df3)

#add in condition and days since AH trt for each bird
df3$condition <- mapping_file$Condition[match(df3$bird_id,mapping_file$phinchID)] 

df3$days_since_AH_trt <- field_data$Days_since_treat[match(df3$bird_id, field_data$bird_id)]
unique(df3$condition)


#just a note "hour" is the number of hours since a bird was injected
class(df3$hour)
df3$hour<-as.character(df3$hour)
head(df3$hour)
df3$hour<-as.numeric(df3$hour)
head(df3)


#now lets find hour of the day and which day it is
df3$hourassec<-df3$hour*60*60
head(df3$hourassec,100)

#make a column for observation time
df3$obstime<-df3$inj_time + df3$hourassec
head(df3,100)
str(df3)

#check on a couple birds and make sure everything looks good still
filter(df3,bird_id==38611)
filter(df3,bird_id==38604)

#make a column for day and hour of experimentation
df3$injectionday<-as.numeric(day(df3$inj_time))

df3$obsday<-as.numeric(day(df3$obstime))

df3$dayofexp<-df3$obsday-df3$injectionday

hours<-hour(df3$obstime)
df3$hourofday<-hours

df3$dayhour<-paste(df3$dayofexp, hours, sep = "_")


#make seperate fen and lps columns
df3$fen<-substr(df3$trt, 1,1)
df3$lps<-substr(df3$trt, 3,3)
head(df3)

#create a variable called hr_of_inj
df3$hr_of_inj<-as.numeric(hour(df3$inj_time))
head(df3)

#create a variable called dayhour 
df3$dayhour<-paste(df3$dayofexp, df3$hourofday, sep = "_")

#make sure that the hour variable is numeric
df3$hour<-as.numeric(df3$hour)

#bin the time since AH trt applied into weeks
df3$days_since_AH_trt = as.character(df3$days_since_AH_trt)
df3$days_since_AH_trt = as.numeric(df3$days_since_AH_trt)
df3$weeks_post_AH_trt <- floor(df3$days_since_AH_trt/7)

#make dummy variables for different smoothers

#alt3 is all fen NO, Y_N, and Y_Y
df3<-df3 %>%
  mutate(alt3 = ifelse(fen =="N", "No_fen",
                       ifelse(fen == "Y" & lps == "Y", "Y_Y",
                              ifelse(fen == "Y" & lps == "N", "Y_N", NA))))
unique(df3$alt3)


#alt4 is all fen YES, N_Y, and N_N
df3<-df3 %>%
  mutate(alt4 = ifelse(fen =="Y", "Fen_yes",
                       ifelse(fen == "N" & lps == "Y", "N_Y",
                              ifelse(fen == "N" & lps == "N", "N_N", NA))))
unique(df3$alt4)



#alt5 is all LPS NO, N_Y and Y_Y
df3<-df3 %>%
  mutate(alt5 = ifelse(lps =="N", "No_lps",
                       ifelse(lps == "Y" & fen == "Y", "Y_Y",
                              ifelse(lps == "Y" & fen == "N", "N_Y", NA))))
unique(df3$alt5)


#alt 6 is all LPS YES, Y_N, and N_N
df3<-df3 %>%
  mutate(alt6 = ifelse(lps =="Y", "Yes_lps",
                       ifelse(lps == "N" & fen == "Y", "Y_N",
                              ifelse(lps == "N" & fen == "N", "N_N", NA))))
unique(df3$alt6)


#look at some snapshots in the data to make sure nothing looks off
filter(df3,dayofexp==0 & hourofday ==2)
filter(df3,bird_id==38646)
filter(df3, bird_id==38631)
filter(df3,bird_id==38625)


#take a rough look at the data, by bird
ggplot(data=df3, aes(x=hour, y= propact)) +         
  geom_point()+                                                 
  geom_line()+                                                     
  facet_wrap(~bird_id)  

# getting transmitter on time from original database
a$on_time_decimal<-hour(a$on_date_time) + (minute(a$on_date_time)/60)
bybird<-group_by(a,bird_id)
birdstart<-summarize(bybird,count = n (), on_time = mean (on_time_decimal, na.rm = T))
df3$on_time_decimal<-birdstart$on_time[match(df3$bird_id, birdstart$bird_id)]
head(df3)

#Make sure condition is a numeric variable
df3$condition = as.numeric(df3$condition)

#write a csv with df3 for future 
# write.csv(df3,file = "df3_activity_by_half_hour.csv", row.names = F)
# getwd()

#filter df3 to exclude hours 0 and 0.5 because birds are freaking out during these times
df3<-df3 %>%
  filter(hour > 0.5)

#only examine the first 8 hours
df_8hr<-df3 %>%
  filter(hour < 8.5)

str(df_8hr)
#gamms 

# diff smoothers for each group (Y_Y,Y_N, N_Y, N_N) and Diff intercepts for group*Site2:
act8.gamm1<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                 random=list(bird_id=~1),data=df_8hr)
act8.gamm1.Exp<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                    random=list(bird_id=~1),data=df_8hr,corr=corExp(form=~hour|bird_id))
act8.gamm1.Exp.nug<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                        random=list(bird_id=~1),data=df_8hr,corr=corExp(form=~hour|bird_id, nugget = T))
act8.gamm1.Gaus<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                     random=list(bird_id=~1),data=df_8hr,corr=corGaus(form=~hour|bird_id))
act8.gamm1.Gaus.nug<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                         random=list(bird_id=~1),data=df_8hr,corr=corGaus(form=~hour|bird_id, nugget = T))
act8.gamm1.Lin<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                    random=list(bird_id=~1),data=df_8hr,corr=corLin(form=~hour|bird_id))
act8.gamm1.Lin.nug<-gamm(propact~fen*lps+weeks_post_AH_trt+condition+s(hour,by = as.factor(trt)),
                        random=list(bird_id=~1),data=df_8hr,corr=corLin(form=~hour|bird_id, nugget = T))
act8.gamm1.Ratio<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                      random=list(bird_id=~1),data=df_8hr,corr=corRatio(form=~hour|bird_id))
act8.gamm1.Ratio.nug<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                          random=list(bird_id=~1),data=df_8hr,corr=corRatio(form=~hour|bird_id, nugget = T))
act8.gamm1.Spher<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                      random=list(bird_id=~1),data=df_8hr,corr=corSpher(form=~hour|bird_id))
act8.gamm1.Spher.nug<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(trt)),
                          random=list(bird_id=~1),data=df_8hr,corr=corSpher(form=~hour|bird_id, nugget = T))


summary(act8.gamm1$lme)


AICc(act8.gamm1,act8.gamm1.Exp,act8.gamm1.Exp.nug,
     act8.gamm1.Gaus,
     act8.gamm1.Lin,act8.gamm1.Lin.nug,
     act8.gamm1.Spher,act8.gamm1.Spher.nug)


# df      AICc
# act8.gamm1           16 -333.9646
# act8.gamm1.Exp       17 -332.5005
# act8.gamm1.Exp.nug   18 -330.9735
# act8.gamm1.Gaus      17 -332.4083
# act8.gamm1.Lin       17 -332.4122
# act8.gamm1.Lin.nug   18 -330.2347
# act8.gamm1.Spher     17 -332.4122
# act8.gamm1.Spher.nug 18 -330.2347



#looks like no structure is the best 
plot(Variogram (act8.gamm1$lme,form=~hour,maxDist=10,resType="n",robust=F))
plot(Variogram (act8.gamm1.Exp.nug$lme,form=~hour,maxDist=10,resType="n",robust=F))
plot(Variogram(act8.gamm1.Ratio$lme,form = ~hour,maxDist = 10,resType = "n",robust = F))
plot(Variogram(act8.gamm1.Ratio.nug$lme,form = ~hour,maxDist = 10,resType = "n",robust = F))


summary(act8.gamm1$lme)



par(mfrow=c(2,2))
plot (act8.gamm1$gam, select =c(1),main="N_N")
plot (act8.gamm1$gam, select =c(2),main="N_Y")
plot (act8.gamm1$gam, select =c(3),main="Y_N")
plot (act8.gamm1$gam, select =c(4),main="Y_Y")


# diff smoothers for each Fen neg group, only one for all Fen pos

act8.gamm2<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(alt4)),
                        random=list(bird_id=~1),data=df_8hr)

par(mfrow=c(2,2))
plot (act8.gamm2$gam, select =c(1),main="Fen_yes")
plot (act8.gamm2$gam, select =c(2),main="N_N")
plot (act8.gamm2$gam, select =c(3),main="N_Y")


summary(act8.gamm2$lme)

#without weeks post AH trt
act8.gamm2.1<-gamm(propact~fen*lps+condition+s(hour,by = as.factor(alt4)),
                 random=list(bird_id=~1),data=df_8hr)


# diff smoothers for each LPS pos group, only one for all LPS neg:
act8.gamm3<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(alt5)),
                    random=list(bird_id=~1),data=df_8hr)
par(mfrow=c(2,2))
plot (act8.gamm3$gam, select =c(2),main="No_lps")
plot (act8.gamm3$gam, select =c(3),main="Y_Y")
plot (act8.gamm3$gam, select =c(1),main="N_Y")

summary(act8.gamm3$lme)

# diff smoothers for each Fen pos group, only one for Fen neg:
act8.gamm4<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(alt3)),
                    random=list(bird_id=~1),data=df_8hr)
par(mfrow=c(2,2))
plot (act8.gamm4$gam, select =c(1),main="No_fen")
plot (act8.gamm4$gam, select =c(2),main="Y_N")
plot (act8.gamm4$gam, select =c(3),main="Y_Y")

summary(act8.gamm4$lme)

# diff smoothers for each LPS neg group, only one for LPS pos:
act8.gamm5<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(alt6)),
                        random=list(bird_id=~1),data=df_8hr)

par(mfrow=c(2,2))
plot (act8.gamm5$gam, select =c(1),main="N_N")
plot (act8.gamm5$gam, select =c(2),main="Y_N")
plot (act8.gamm5$gam, select =c(3),main="Yes_lps")

summary(act8.gamm5$lme)
# Fixed effects: y ~ X - 1 
# Value  Std.Error  DF   t-value p-value
# X(Intercept)                        0.6645415 0.03562310 397 18.654794  0.0000
# XfenY                               0.0171464 0.03621635  27  0.473445  0.6397
# XlpsY                              -0.1101642 0.03315720  27 -3.322482  0.0026
# Xcondition                         -0.0304263 0.01347571  27 -2.257864  0.0322
# Xweeks_post_AH_trt                 -0.0242532 0.02224840  27 -1.090110  0.2853
# XfenY:lpsY                          0.0115463 0.04609928  27  0.250466  0.8041
# Xs(hour):as.factor(alt6)N_NFx1      0.0488609 0.01601345 397  3.051241  0.0024
# Xs(hour):as.factor(alt6)Y_NFx1     -0.0084355 0.01696505 397 -0.497228  0.6193
# Xs(hour):as.factor(alt6)Yes_lpsFx1 -0.0743736 0.05920863 397 -1.256128  0.2098
gam.check(act8.gamm5$gam)

#without weeks post AH trt
act8.gamm5.1<-gamm(propact~fen*lps+condition+s(hour,by = as.factor(alt6)),
                 random=list(bird_id=~1),data=df_8hr)
summary(act8.gamm5.1$lme)
# X(Intercept)                        0.6364851 0.02506113 397 25.397305  0.0000
# XfenY                               0.0171101 0.03678835  28  0.465096  0.6455
# XlpsY                              -0.1177709 0.03291789  28 -3.577716  0.0013
# Xcondition                         -0.0287900 0.01360438  28 -2.116230  0.0433
# XfenY:lpsY                          0.0115084 0.04682461  28  0.245777  0.8076
# Xs(hour):as.factor(alt6)N_NFx1      0.0484977 0.01599752 397  3.031577  0.0026
# Xs(hour):as.factor(alt6)Y_NFx1     -0.0082537 0.01694819 397 -0.486995  0.6265
# Xs(hour):as.factor(alt6)Yes_lpsFx1 -0.0746606 0.05938739 397 -1.257179  0.2094

#2 smoothers, one for fen + one for fen-
act8.gamm6<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(fen)),
                    random=list(bird_id=~1),data=df_8hr)
#Error in lme.formula(y ~ X - 1, random = rand, data = strip.offset(mf),  : 
#                       nlminb problem, convergence error code = 1
#                     message = singular convergence (7)


#2 smoothers, one for lps+ one for lps-
act8.gamm7<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(lps)),
                    random=list(bird_id=~1),data=df_8hr)

par(mfrow=c(2,2))
plot (act8.gamm7$gam, select =c(1),main="lps+")
plot (act8.gamm7$gam, select =c(2),main="lps-")


summary(act8.gamm7$lme)


#without weeks post AH trt
act8.gamm7.1<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour,by = as.factor(lps)),
                 random=list(bird_id=~1),data=df_8hr)

summary(act8.gamm7.1$lme)


#1 smoother for all birds
act8.gamm8<-gamm(propact~fen*lps+condition+weeks_post_AH_trt+s(hour),
                    random=list(bird_id=~1),data=df_8hr)

par(mfrow=c(2,2))
plot (act8.gamm8$gam, select =c(1),main="all")

#without weeks post AH trt
act8.gamm8.1<-gamm(propact~fen*lps+condition+s(hour),
                 random=list(bird_id=~1),data=df_8hr)


summary(act8.gamm8$lme)
# Fixed effects: y ~ X - 1 
# Value  Std.Error  DF   t-value p-value
# X(Intercept)        0.6636150 0.03532014 399 18.788570  0.0000
# XfenY               0.0158871 0.03591103  27  0.442402  0.6617
# XlpsY              -0.1121946 0.03289489  27 -3.410699  0.0021
# Xcondition         -0.0304748 0.01335915  27 -2.281191  0.0306
# Xweeks_post_AH_trt -0.0222847 0.02206332  27 -1.010036  0.3214
# XfenY:lpsY          0.0131015 0.04572153  27  0.286549  0.7766
# Xs(hour)Fx1         0.0126592 0.00766158 399  1.652298  0.0993

AICc(act8.gamm1,act8.gamm2,act8.gamm2.1,act8.gamm3,act8.gamm4,
     act8.gamm5,act8.gamm5.1,act8.gamm6,act8.gamm7,act8.gamm7.1,act8.gamm8,act8.gamm8.1)

AICc(act8.gamm2.1,act8.gamm5.1,act8.gamm7.1,act8.gamm8.1)

#========================for CHANGE=========================================
####Graphing change in activity###
#find average prop act for each treatment group by hours since injection
head(df3)
act_mean<-tapply(df3$propact, list(df3$trt,df3$hour), mean, na.rm=T)

#find ses by treatment group by hours since injection
act_ses<-tapply(df3$propact, list(df3$trt, df3$hour), se)
act_ses<-as.vector(act_ses)
length(act_ses)


#take row names (lps Y or N) from tch_avs and repeat the same number of times as there are columns in tch_avs
dim(tch_avs)
trt<-rep(dimnames(act_mean)[[1]],93)
length(trt)

#take column names from tch_avs (time interval in half hours since injection) and repeat the same number of times are there are rows (2)
hours_post_inj<-rep(dimnames(act_mean)[[2]], each = nrow(act_mean))
length(hours_post_inj)
hours_post_inj<-as.character(hours_post_inj)
hours_post_inj<-as.numeric(hours_post_inj)

#make act_mean into a vector
act_mean<-as.vector(act_mean)

#make a new dataframe
act_by_trt<-data.frame(hours_post_inj, act_mean,act_ses, trt)

#separate out different treatments from "trt" variable
act_by_trt$fen<-substr(act_by_trt$trt, 1,1)
act_by_trt$lps<-substr(act_by_trt$trt, 3,3)

#make variables that describe what treatment a bird got for graphing
act_by_trt$fenlevel<-factor(act_by_trt$fen, labels = c("Control (water)", "Fenbendazole"))
act_by_trt$lpslevel<-factor(act_by_trt$lps, labels = c("Control (no LPS)", "LPS"))

#add in variables that describe the upper and lower mean+/-SE boundaries
act_by_trt$meanANDse<-act_by_trt$act_mean + act_by_trt$act_ses
act_by_trt$meanMINUSse<-act_by_trt$act_mean - act_by_trt$act_ses

#only include first 8 hours of data
act_by_trt<-act_by_trt %>%
  filter(hours_post_inj < 8.5)

#ugly plot, just to check
act<-ggplot(act_by_trt, aes(x=hours_post_inj, y = act_mean, colour = lps))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=act_mean-act_ses, ymax=act_mean+act_ses), width=.2)+
  facet_wrap(~fen)

#make a file for making a pretty figure later
write.csv(act_by_trt, file = "Act_data_for_figure3.csv", row.names = F)



#ugly plot just to get a preview
ggplot(act_by_trt, aes(x=hours_post_inj, y = act_mean, colour = lps))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=act_mean-act_ses, ymax=act_mean+act_ses), width=.2)+
  facet_wrap(~fen)
