library(ggplot2)
library(grid)
library(dplyr)
require(RColorBrewer)

#read in data
temp<- read.csv("Temp_data_for_figure2.csv")
act<-read.csv("Act_data_for_figure3.csv")

#look at names and make sure data all look right
head(temp)
head(act)

#change name of hours column in act data to match the temp data
colnames(act)[1]<-"hours_since_inj"

#combine act and temp
tempact<-right_join(act, temp, by = c("hours_since_inj", "trt"))
dim(temp) #176 8
dim(act) #60 10
dim(tempact) #176 16


#get rid of unnecessary columns
tempact2<-tempact %>%
  select(hours_since_inj, act_mean, act_ses, meanANDse.x, meanMINUSse.x, tch_avs, tch_ses, meanANDse.y, meanMINUSse.y, trt, fenlevel, lpslevel)

#break the trt column into seperate treatment columns
tempact2$fenlevel = substr(tempact2$trt,1,1)
tempact2$lpslevel = substr(tempact2$trt,3,3)

#rename columns to make it easier to read
colnames(tempact2)[4]<-"upper_se_act"
colnames(tempact2)[5]<-"lower_se_act"
colnames(tempact2)[6]<-"av_temp_change"
colnames(tempact2)[7]<-"se_temp_change"
colnames(tempact2)[8]<-"upper_se_temp"
colnames(tempact2)[9]<-"lower_se_temp"
colnames(tempact2)[10]<-"trt"
colnames(tempact2)[11]<-"fen"
colnames(tempact2)[12]<-"lps"

#make the fen and lps columns factors, rename the levels, and put them in the order you want them to display in the graph
tempact2$lps<-as.factor(tempact2$lps)
levels(tempact2$lps)
tempact2$lps<-factor(tempact2$lps, labels = c("Control (no LPS)","LPS"))
tempact2$lps<-ordered(tempact2$lps,levels = c("Control (no LPS)", "LPS"))

tempact2$fen<-as.factor(tempact2$fen)
levels(tempact2$fen)
tempact2$fen<-factor(tempact2$fen, labels = c("Control (water)","Anthelminthic"))
tempact2$fen<-ordered(tempact2$fen, levels =c("Control (water)", "Anthelminthic"))

#=====================================================================================================

#separate graphs
tempplot<-tempact2 %>%
  select(hours_since_inj, av_temp_change, upper_se_temp, lower_se_temp, lps, fen) %>%
  na.omit() %>%
  ggplot() +
  facet_wrap(~fen, scales = "free")+
  geom_line(aes(x = hours_since_inj, y = av_temp_change, group= lps, linetype = lps),size=1) +
  geom_ribbon(aes(x = hours_since_inj, ymin = lower_se_temp, ymax= upper_se_temp, fill = lps),alpha=0.65)+
  ylim(-2.7,1.5)+
  scale_linetype_manual(values = c(1,3))+
  scale_fill_manual(values = c("gray","slategray")) +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.25,0.9),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.line = element_line(),
        axis.title = element_text(size = 14),
         panel.background = element_blank())+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=1.5)+
  labs(x = "Time since injection (hours)", y = "Change in skin temperature (°C)")+
  theme(plot.title = element_text(hjust = 0.5))

tempplot


# tiff("Figure2.tiff", units="in", width=6.5, height=7, res=300)
# tempplot
# dev.off()

actplot<-tempact2 %>%
  select(hours_since_inj, act_mean, upper_se_act, lower_se_act, fen, lps) %>%
  filter(hours_since_inj<8.5)%>%
  na.omit() %>%
  ggplot() +
  facet_wrap(~fen, scales = "free")+
  ylim(0.23,0.83)+
  geom_line(aes(x = hours_since_inj, y = act_mean, group= lps, linetype = lps), size=1) +
  geom_ribbon(aes(x = hours_since_inj, ymin = lower_se_act, ymax= upper_se_act, fill = lps), alpha=0.65)+
  scale_linetype_manual(values = c(1,3))+
  scale_fill_manual(values = c("gray","slategray")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.25,0.9),
        legend.text = element_text(size = 12),
        axis.line = element_line(),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank())+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=0.8)+
  labs(x="Time since injection (hours)", y = "Proportion of time spent active")
actplot 

# tiff("Figure3.tiff", units="in", width=6.5, height=7, res=300)
# actplot 
# dev.off()


