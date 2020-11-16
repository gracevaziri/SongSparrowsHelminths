#functions
se<-function(x) { sqrt((var(na.omit(x))/length(na.omit(x))))}   #calculates std error

#packages
library(lubridate)
library(tidyverse)
library(dplyr)
library(nlme)
library(ggplot2)
library(mgcv)
library(MuMIn)

##read in the data
#raw data
a<-read.csv("/Users/gracevaziri/OneDrive/Grace Thesis/Analyses/GAMMs/WA Temp/All birds_best antenna_with temp.csv")
#bring in the information on bird condition
mapping_file=read_tsv("16S_mapping_20180518.txt")  


##create a column with the date and time that the receiver started collecting data called "data_date_time
##create another column with the date and time that the transmitter was actually placed on the bird called "start_date_time"
a$date_time<-paste(a$Date,a$WA.time)
a$on_date_time<-paste(a$Transmitter.date,a$Transmitter_on)


##convert with lubridate 
#lubridate compatible date/time of observations
a$date_time<-mdy_hms(a$date_time, tz = "US/Pacific")
a$date_time<-round_date(a$date_time,unit = "30 mins")
#lubridate compatible date/time of start of data collection
a$on_date_time<-mdy_hms(a$on_date_time, tz = "US/Pacific")
a$on_date_time<-round_date(a$on_date_time, unit = "30 mins")

#add data on condition          
a$condition <- mapping_file$Condition[match(a$bird_id,mapping_file$phinchID)]
a2<-as_tibble(a)

#cut this down so that I only have what I need
b <- a2 %>%
  dplyr::select(bird_id,comb_trt,temp, date_time, on_date_time,condition)


##get rid of any rows for which temperature was greater than 45 or less than 36
dplyr::filter(b, b$temp <=45 & b$temp >=35)
c<-filter(b, b$temp <=45 & b$temp >=35)

#find average temp  per bird per time interval
av_temp<-tapply(c$temp,list(c$bird_id, c$date_time), mean, na.rm=T)
# ncol(av_temp)
# dim(av_temp)
# length(av_temp)

#take row names (bird_id) from bird_av_temp and repeat the same number of times as there are columns in bird_av_temp
bird_id<-rep(dimnames(av_temp)[[1]], 996)
length(bird_id)

#take column names from bird_av_temp (day intervals) and repeat the same number of times are there are rows (39)
time<-rep(dimnames(av_temp)[[2]], each = nrow(av_temp))
av_temp<-as.vector(av_temp)

#make a new dataframe with the average temp per day per bird, bird ID, and day
df_ch<-data.frame(bird_id,av_temp,time)

#Match the start time with the bird
df_ch$inj_time<-c$on_date_time[match(df_ch$bird_id,c$bird_id)]

#add in treatment information for each bird
df_ch$trt<-c$comb_trt[match(df_ch$bird_id,c$bird_id)]

#add in the birds' conditions
df_ch$condition <- c$condition[match(df_ch$bird_id, c$bird_id)]

#find the interval of time between when a bird was treated and when a measurement was recorded
df_ch$time_span<-interval(df_ch$inj_time,df_ch$time)

#find the number of hours in the interval
df_ch$hours_since_inj<-df_ch$time_span/ dhours(x=1)

#find the time of day an observation was made
t.lub<-ymd_hms(df_ch$time)
h.lub<-hour(t.lub) + minute(t.lub)/60
df_ch$obs_time<-h.lub
df_ch$obs_time_div_100<-as.numeric(df_ch$obs_time/100)


#cleanup df_ch so that we just have the columns we're interested in 
df_ch<-df_ch %>%
  dplyr::select(bird_id, av_temp, trt, hours_since_inj, obs_time, obs_time_div_100, condition)
head(df_ch)

#for purposes of making an exploratory plot, get rid of rows with NA in av_temp_ch
df_ch_no_NA<-na.omit(df_ch)

#plot birds individually to see if anyone is looking weird
ggplot(data=df_ch_no_NA, aes(x=hours_since_inj, y= av_temp)) +         ##throw out 38608,
  geom_point()+                                                 ##only use first 24 hours of data, also, throw out 38635,
  geom_line()+                                                  ## NOTE: figure out whether we can use 38637 (lots of missing data in beginning)    
  facet_wrap(~bird_id)                                          ##also, make a rule to get rid of weird spikes

#arrange the dataframe df_ch by each bird
arrange(df_ch, bird_id)

#get rid of any row of data for which the value in 'hours_since_inj' variable is 0.5 (because bird is still recovering from being handled at this time)
df_ch<-df_ch %>%
  filter(hours_since_inj > 0.5)

#make a column with average av_temp observations for hours 1.0, 1.5, and 2.0 for all birds that have all three, 
#or any combination of the three, else fill the column with the first temp reading for the bird

#first, make 'hours_since_inj' a factor
df_ch$hours_since_inj<-as.factor(df_ch$hours_since_inj)

#then make a column with the first observation for each group

#make a dataframe that combines the av_temp_ch for the 1.0, 1.5, and 2.0 time intervals for each bird
hour1<-df_ch %>%
  filter(hours_since_inj == "1")

hour1.5<-df_ch %>%
  filter(hours_since_inj=="1.5")

hour2<-df_ch %>%
  filter(hours_since_inj=="2")

#make a vector called birds that is a list of all the bird names
birds<-unique(df_ch$bird_id)

#make the vector into a dataframe
start.times<-data.frame(birds)
start.times

#match the appropriate av temp ch reading for each hour for each bird in the start.times dataframe
start.times$hour1<-hour1$av_temp[match(start.times$birds,hour1$bird_id)]
start.times$hour1.5<-hour1.5$av_temp[match(start.times$birds,hour1.5$bird_id)]
start.times$hour2<-hour2$av_temp[match(start.times$birds,hour2$bird_id)]

#get the row mean for each of temp changes, and now identify which birds are missing all three reading for the first three time intervals
start.times$avstarttemp<- rowMeans(start.times[,2:4], na.rm = T)             #birds 38604,38608,38637, and 38653 are missing all three

#manually add in the first start temp for the birds that are missing the avstarttemp
#create a dataframe with bird id and the first temp and the hours since injection that the first temp was measured
#make the bird id variable
bird_ID<-c('38604','38608','38637',  '38653')
#find the first temp for each of these birds
filter(df_ch, bird_id=="38604")            #first temp ch measurement was  0.37843285 at 2.5 hrs post injection
filter(df_ch, bird_id=="38608")            #first temp ch measurement was 14.53022 at 8.5 hrs post injection<-throw out this bird, has very little data
filter(df_ch, bird_id=="38637")            #first temp ch measurement was 0.63144171 at 4 hrs post injection                   
filter(df_ch, bird_id=="38653")            #first temp ch measurement was 0.198404351 at 2.5 hrs post injection    

firsttemp<-c(42.19910,37.77574,42.06932,42.71940 )
timeoffirsttemp<-c(2.5,8.5,4,2.5)

#now make the dataframe
oddballstarttemps<-data.frame(bird_ID, firsttemp, timeoffirsttemp)
head(oddballstarttemps)

#match the firsttemp from oddballs to start.time dataframe for the missing birds
start.times$firsttemp<-oddballstarttemps$firsttemp[match(start.times$birds,oddballstarttemps$bird_ID)]

#combining the avstarttemp and firsttemp columns into a new column called first with all the first temp change recordings for birds
start.times$first<-rowMeans(start.times[,c(5:6)], na.rm = T)

#Match the start temp time with the bird
df_ch$start_temp<-start.times$first[match(df_ch$bird_id,start.times$birds)]

#check that it worked
filter(df_ch, bird_id=="38605")
filter(df_ch, bird_id=="38604")   

#tidy up the dataframe
df_ch<-df_ch %>%
  dplyr::select(bird_id, av_temp, trt, hours_since_inj, obs_time, obs_time_div_100, start_temp, condition)

#remove all the data for hours since injection < 2.5 from dataframe 
#first change hours since inj back to numeric
df_ch$hours_since_inj<-as.character(df_ch$hours_since_inj)
df_ch$hours_since_inj<-as.numeric(df_ch$hours_since_inj)

df3.1<-df_ch %>%
  filter(hours_since_inj >2)

range(df3.1$hours_since_inj)  #check, and it looks good

#create a column for the hours since first temp ch measure
#go back to the start.times dataframe
#first temp ch recording happened at hour 2, unless it was an oddball bird
start.times<-start.times %>%
  mutate(time1 = ifelse(avstarttemp =="NaN", NA, 2))

#now add in hours for the oddballs
oddballstarttemps
start.times$time_other<-oddballstarttemps$timeoffirsttemp[match(start.times$birds,oddballstarttemps$bird_ID)]

#now merge columns time1 and time_other
start.times$first_time<-rowMeans(start.times[,c(8:9)], na.rm=T)

#match the first_time column to the df3 dataframe
df3.1$first_time<-start.times$first_time[match(df3.1$bird_id, start.times$birds)]

#get rid of the data before the start time for odball birds
oddballstarttemps

df3.1<-df3.1 %>%
  filter(!(bird_id == "38604" & hours_since_inj < 3))
filter(df3.1, bird_id == "38604")

df3.1<-df3.1 %>%
  filter(!(bird_id == "38608" & hours_since_inj < 9))
filter(df3.1, bird_id == "38608")

df3.1<-df3.1 %>%
  filter(!(bird_id == "38637" & hours_since_inj < 4.5))
filter(df3.1, bird_id == "38637")


df3.1<-df3.1 %>%
  filter(!(bird_id == "38653" & hours_since_inj < 3))
filter(df3.1, bird_id == "38653")

df3.1<-arrange(df3.1, bird_id)
df3.1<-na.omit(df3.1)

#get rid of everything past 48 hours
df4.1<-df3.1 %>%
  filter(hours_since_inj <49)

##bin the time since AH trt applied into weeks
df3.1$days_since_AH_trt <- field_data$Days_since_treat[match(df3.1$bird_id, field_data$bird_id)]
df3.1$days_since_AH_trt = as.character(df3.1$days_since_AH_trt)
df3.1$days_since_AH_trt = as.numeric(df3.1$days_since_AH_trt)

df3.1$weeks_post_AH_trt <- floor(df3.1$days_since_AH_trt/7)

#reorder columns to make this easier to look at
df4.1<-df4.1 %>%
  dplyr::select(bird_id, trt, obs_time, obs_time_div_100, hours_since_inj, first_time,start_temp, av_temp, condition)


#find the number of hours that have elapsed since the time used to calculate the average start temp, and the time since a bird was injected
df4.1$hours_since_first_measure<-df4.1$hours_since_inj-df4.1$first_time
df4.1$diff_hr<-ave(df4.1$hours_since_first_measure, df4.1$bird_id, FUN=function(x) c(0, diff(x)))


#find the difference in temperature between a current temp and the start temp
df4.1$tempdiff<-df4.1$av_temp-df4.1$start_temp

#lets check it out
ggplot(data=df4.1, aes(x=hours_since_first_measure, y= tempdiff)) +         ##still need to throw out 38608,
  geom_point()+                                                 ##only use first 24 hours of data, also, get rid of 38635,
  geom_line()+                                                 ##also get rid of any rows for which tempdifff is greater than 4 or less that -4
  facet_wrap(~bird_id)                                         

#now lets get rid of birds with malfunctioning transmitters and crazy spikes
df4.1<-filter(df4.1,bird_id!= "38608" &
              bird_id!= "38635")

df4.1<-filter(df4.1,tempdiff <= 4 &
              tempdiff >=-4)

range(df4.1$tempdiff)

#get rid of stuff past 24 hours
df4.1<-filter(df4.1,hours_since_inj <= 24)



#lets check it again
ggplot(data=df4.1, aes(x=hours_since_first_measure, y= tempdiff)) +         ##this is looking better
  geom_point()+                                                 
  geom_line()+                                                 
  facet_wrap(~bird_id)                                         

#make some smoothers for gamms
#make a fen column
df4.1$fen<-substr(df4.1$trt,1,1)

#make an lps column
df4.1$lps<-substr(df4.1$trt, 3,3)

#make a variable for 3 smoothers, treat all birds that didn't get fen as a group (No_fen), 
#and treat birds that got fen + lps (Y_Y) and fen+h20 as separate groups (Y_N)
df5.1<-df4.1 %>%
  mutate(alt3 = ifelse(fen =="N", "No_fen",
                       ifelse(fen == "Y" & lps == "Y", "Y_Y",
                              ifelse(fen == "Y" & lps == "N", "Y_N", NA))))

#make a variable for 3 smoothers, treat all birds that DID get fen as a group (Fen_yes), 
#and treat birds that got h20 + lps (N_Y) and h20+h20 as separate groups (N_N)
df5.1<-df5.1 %>%
  mutate(alt4 = ifelse(fen =="Y", "Fen_yes",
                       ifelse(fen == "N" & lps == "Y", "N_Y",
                              ifelse(fen == "N" & lps == "N", "N_N", NA))))

#make a variable for 3 smoothers, treat all birds that didn't get lps (so the fen+h20 and the h20+h20) as a group (No_lps), 
#and treat birds that got fen + lps (Y_Y) and h20+lps as separate groups (N_Y)
df5.1<-df5.1 %>%
  mutate(alt5 = ifelse(lps =="N", "No_lps",
                       ifelse(lps == "Y" & fen == "Y", "Y_Y",
                              ifelse(lps == "Y" & fen == "N", "N_Y", NA))))


#make a variable for 3 smoothers, treat all birds that DID get lps (so the fen+lps and the h20+lps) as a group (Yes_lps), 
#and treat birds that got fen + h20 (Y_N) and h20+h20 as separate groups (N_N)
df5.1<-df5.1 %>%
  mutate(alt6 = ifelse(lps =="Y", "Yes_lps",
                       ifelse(lps == "N" & fen == "Y", "Y_N",
                              ifelse(lps == "N" & fen == "N", "N_N", NA))))

# getting transmitter on time from original data base
a$on_time_decimal<-hour(a$on_date_time)+(minute(a$on_date_time)/60)
bybird<-group_by(a,bird_id)
birdstart <- summarize(bybird, count = n(), on_time = mean(on_time_decimal, na.rm = T))
df5.1$on_time_decimal<-birdstart$on_time[match(df5.1$bird_id,birdstart$bird_id)]

df5.1$condition<-as.character(df5.1$condition)
df5.1$condition<-as.numeric(df5.1$condition)


df5.1$weeks_post_AH_trt = df3.1$weeks_post_AH_trt[match(df5.1$bird_id,df3.1$bird_id)]
str(df5.1)

#do the gamms
fix.temp.ch.gamm1<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                    random=list(bird_id=~1),data=df5.1)
fix.temp.ch.gamm1.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                        random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))
fix.temp.ch.gamm1.Exp.nug<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                            random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id, nugget = T))
fix.temp.ch.gamm1.Gaus<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                         random=list(bird_id=~1),data=df5.1,corr=corGaus(form=~hours_since_inj|bird_id))
fix.temp.ch.gamm1.Gaus.nug<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                             random=list(bird_id=~1),data=df5.1,corr=corGaus(form=~hours_since_inj|bird_id, nugget = T))
fix.temp.ch.gamm1.Lin<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                        random=list(bird_id=~1),data=df5.1,corr=corLin(form=~hours_since_inj|bird_id))
fix.temp.ch.gamm1.Lin.nug<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                            random=list(bird_id=~1),data=df5.1,corr=corLin(form=~hours_since_inj|bird_id, nugget = T))
fix.temp.ch.gamm1.Ratio<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                          random=list(bird_id=~1),data=df5.1,corr=corRatio(form=~hours_since_inj|bird_id))
fix.temp.ch.gamm1.Ratio.nug<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                              random=list(bird_id=~1),data=df5.1,corr=corRatio(form=~hours_since_inj|bird_id, nugget = T))
fix.temp.ch.gamm1.Spher<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                          random=list(bird_id=~1),data=df5.1,corr=corSpher(form=~hours_since_inj|bird_id))
fix.temp.ch.gamm1.Spher.nug<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(trt)),
                              random=list(bird_id=~1),data=df5.1,corr=corSpher(form=~hours_since_inj|bird_id, nugget = T))


summary(fix.temp.ch.gamm1$lme)


anova(fix.temp.ch.gamm1$lme,fix.temp.ch.gamm1.Exp$lme,fix.temp.ch.gamm1.Exp.nug$lme,
      fix.temp.ch.gamm1.Gaus$lme,fix.temp.ch.gamm1.Gaus.nug$lme,
      fix.temp.ch.gamm1.Lin$lme,fix.temp.ch.gamm1.Lin.nug$lme,
      fix.temp.ch.gamm1.Ratio$lme,fix.temp.ch.gamm1.Ratio.nug$lme,
      fix.temp.ch.gamm1.Spher$lme,fix.temp.ch.gamm1.Spher.nug$lme)

require(MuMIn)
AICc(fix.temp.ch.gamm1,fix.temp.ch.gamm1.Exp,fix.temp.ch.gamm1.Exp.nug,
     fix.temp.ch.gamm1.Gaus,fix.temp.ch.gamm1.Gaus.nug,
     fix.temp.ch.gamm1.Lin, fix.temp.ch.gamm1.Lin.nug,
     fix.temp.ch.gamm1.Ratio,fix.temp.ch.gamm1.Ratio.nug,
     fix.temp.ch.gamm1.Spher,fix.temp.ch.gamm1.Spher.nug)

#Exp is best correlation structure
plot(Variogram (fix.temp.ch.gamm1$lme,form=~hours_since_inj,maxDist=10,resType="n",robust=F))
plot(Variogram (fix.temp.ch.gamm1.Exp$lme,form=~hours_since_inj,maxDist=10,resType="n",robust=F))

summary(fix.temp.ch.gamm1.Exp$lme)

par(mfrow=c(2,2))
plot (fix.temp.ch.gamm1.Exp$gam, select =c(1),main="N_N")
plot (fix.temp.ch.gamm1.Exp$gam, select =c(2),main="N_Y")
plot (fix.temp.ch.gamm1.Exp$gam, select =c(3),main="Y_N")
plot (fix.temp.ch.gamm1.Exp$gam, select =c(4),main="Y_Y")

# diff smoothers for each Fen neg group, only one for all Fen pos
fix.temp.ch.gamm2.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(alt4)),
                        random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))

par(mfrow=c(2,2))
plot (fix.temp.ch.gamm2.Exp$gam, select =c(1),main="Fen_yes")
plot (fix.temp.ch.gamm2.Exp$gam, select =c(2),main="N_N")
plot (fix.temp.ch.gamm2.Exp$gam, select =c(3),main="N_Y")


summary(fix.temp.ch.gamm2.Exp$lme)


# diff smoothers for each LPS pos group, only one for all LPS neg:
fix.temp.ch.gamm3.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(alt5)),
                        random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))
par(mfrow=c(2,2))
plot (fix.temp.ch.gamm3.Exp$gam, select =c(2),main="No_lps")
plot (fix.temp.ch.gamm3.Exp$gam, select =c(3),main="Y_Y")
plot (fix.temp.ch.gamm3.Exp$gam, select =c(1),main="N_Y")

summary(fix.temp.ch.gamm3.Exp$lme)

# diff smoothers for each Fen pos group, only one for Fen neg:
fix.temp.ch.gamm4.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(alt3)),
                        random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))
par(mfrow=c(2,2))
plot (fix.temp.ch.gamm4.Exp$gam, select =c(1),main="No_fen")
plot (fix.temp.ch.gamm4.Exp$gam, select =c(2),main="Y_N")
plot (fix.temp.ch.gamm4.Exp$gam, select =c(3),main="Y_Y")

summary(fix.temp.ch.gamm4.Exp$lme)

# diff smoothers for each LPS neg group, only one for LPS pos:
fix.temp.ch.gamm5.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(alt6)),
                        random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))

par(mfrow=c(2,2))
plot (fix.temp.ch.gamm5.Exp$gam, select =c(1),main="N_N")
plot (fix.temp.ch.gamm5.Exp$gam, select =c(2),main="N_Y")
plot (fix.temp.ch.gamm5.Exp$gam, select =c(3),main="Yes_lps")

summary(fix.temp.ch.gamm5.Exp$lme)

#wihout weeks post AH trt
fix.temp.ch.gamm5.1.Exp<-gamm(tempdiff~fen*lps+condition+s(hours_since_inj,by = as.factor(alt6)),
                            random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))

#diff smoothers by fen pos or neg:
fix.temp.ch.gamm6.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(fen)),
                            random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))

par(mfrow=c(1,2))
plot (fix.temp.ch.gamm6.Exp$gam, select =c(1),main="Fen yes")
plot (fix.temp.ch.gamm6.Exp$gam, select =c(2),main="Fen no")

summary(fix.temp.ch.gamm6.Exp$lme)

#diff smoothers by LPS:
fix.temp.ch.gamm7.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj,by = as.factor(lps)),
                        random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))
par(mfrow=c(1,2))
plot (fix.temp.ch.gamm7.Exp$gam, select =c(1),main="LPS no")
plot (fix.temp.ch.gamm7.Exp$gam, select =c(2),main="LPS yes")

summary(fix.temp.ch.gamm7.Exp$lme)

#wihout weeks post AH trt
fix.temp.ch.gamm7.1.Exp<-gamm(tempdiff~fen*lps+condition+s(hours_since_inj,by = as.factor(lps)),
                            random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))

#only 1 smoother:
fix.temp.ch.gamm8.Exp<-gamm(tempdiff~fen*lps+condition+weeks_post_AH_trt+s(hours_since_inj),
                         random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))
par(mfrow=c(1,2))
plot (fix.temp.ch.gamm8.Exp$gam, select =c(1),main="All birds")
summary(fix.temp.ch.gamm8.Exp$lme)


#only 1 smoother:#wihout weeks post AH trt
fix.temp.ch.gamm8.1.Exp<-gamm(tempdiff~fen*lps+condition+s(hours_since_inj),
                            random=list(bird_id=~1),data=df5.1,corr=corExp(form=~hours_since_inj|bird_id))
par(mfrow=c(1,2))
plot (fix.temp.ch.gamm8.1.Exp$gam, select =c(1),main="All birds")
summary(fix.temp.ch.gamm8.1.Exp$lme)



AICc(fix.temp.ch.gamm1.Exp,fix.temp.ch.gamm2.Exp,fix.temp.ch.gamm3.Exp,
     fix.temp.ch.gamm4.Exp,fix.temp.ch.gamm5.1.Exp,
     fix.temp.ch.gamm6.Exp,fix.temp.ch.gamm7.1.Exp,
     fix.temp.ch.gamm8.Exp,
     fix.temp.ch.gamm8.1.Exp)

#this says the best model uses only one smoother for all birds, and doesn't include weeks since trt
#                         df     AICc
# fix.temp.ch.gamm1.Exp   17 1244.729
# fix.temp.ch.gamm2.Exp   15 1233.277
# fix.temp.ch.gamm3.Exp   15 1244.211
# fix.temp.ch.gamm4.Exp   15 1244.385
# fix.temp.ch.gamm5.1.Exp 14 1227.524
# fix.temp.ch.gamm6.Exp   13 1232.528
# fix.temp.ch.gamm7.1.Exp 12 1227.237
# fix.temp.ch.gamm8.Exp   11 1215.393
# fix.temp.ch.gamm8.1.Exp 10 1213.640
                  
summary(fix.temp.ch.gamm8.1.Exp$lme)
# Fixed effects: y ~ X - 1 
# Value Std.Error   DF    t-value p-value
# X(Intercept)           -0.4758804 0.2131810 1364 -2.2322833  0.0258
# XfenY                  -0.4363120 0.3128487   28 -1.3946423  0.1741
# XlpsY                  -0.1689220 0.2777047   28 -0.6082792  0.5479
# Xcondition              0.2159157 0.1148104   28  1.8806284  0.0705
# XfenY:lpsY              1.1684797 0.3944065   28  2.9626281  0.0062
# Xs(hours_since_inj)Fx1  0.8025187 0.2982244 1364  2.6909899  0.0072


#========================getting data for temp graph=========================================
tch<-df5.1
head(tch)
names(tch)

#make a column that separates the fenY lpsY treatment birds from the rest
tch<-tch %>%
  mutate(Y_Y = ifelse(lps =="N", "N",
                      ifelse(lps == "Y" & fen == "Y", "Y_Y",
                             ifelse(fen == "N", "N", NA))))
unique(tch$Y_Y)


#find average temperatures for each treatment group by hours since injection
head(tch)
tch_avs<-tapply(tch$tempdiff, list(tch$trt, tch$hours_since_inj), mean, na.rm=T)
head(tch_avs)
dim(tch_avs)

#find ses by treatment group by hours since injection
tch_ses<-tapply(tch$tempdiff, list(tch$trt, tch$hours_since_inj), se)
tch_ses<-as.vector(tch_ses)

#take row names (lps Y or N) from tch_avs and repeat the same number of times as there are columns in tch_avs
trt<-rep(dimnames(tch_avs)[[1]], 44)
length(trt)

#take column names from tch_avs (time interval in half hours since injection) and repeat the same number of times are there are rows (2)
hours_since_inj<-rep(dimnames(tch_avs)[[2]], each = nrow(tch_avs))
length(hours_since_inj)
hours_since_inj<-as.character(hours_since_inj)
hours_since_inj<-as.numeric(hours_since_inj)
#make tch_avs into a vector
tch_avs<-as.vector(tch_avs)

#make a new dataframe
av_by_trt<-data.frame(hours_since_inj, tch_avs,tch_ses, trt)
str(av_by_trt)
head(av_by_trt)

av_by_trt$fen<-substr(av_by_trt$trt, 1,1)
av_by_trt$lps<-substr(av_by_trt$trt, 3,3)

av_by_trt$meanANDse<-av_by_trt$tch_avs + av_by_trt$tch_ses

av_by_trt$meanMINUSse<-av_by_trt$tch_avs - av_by_trt$tch_ses

head(av_by_trt)

write.csv(av_by_trt, file = "Temp_data_for_figure2.csv", row.names = F)

