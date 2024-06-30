#make precip and NDVI figures, divide seasons by NDVI

library(lubridate, warn.conflicts = FALSE)
library("dplyr")
library(ggplot2)
library(cowplot)

setwd("C:/Users/Sara//Dropbox/diet_specialization/R_code") #laptop

NDVI<-read.csv(file = "whitewater_ndvi_evi_colnames_19891001_20231218.csv") 
precip<-read.csv(file = "whitewater_precip_colnames_19801001_202311XX.csv")

str(NDVI)

NDVI$date<-as.Date(NDVI$date, format= "%m/%d/%Y")

NDVI$month<-month(NDVI$date)
NDVI$monthF <- factor(format(NDVI$date, "%b"), month.abb, ordered = TRUE)
NDVI$monthF2<-factor(NDVI$monthF , 
	levels= c( "Dec","Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
			"Aug", "Sep", "Oct", "Nov") )
NDVI$year<-year(NDVI$date)

NDVI$JulDate<- as.integer(format(NDVI$date, "%j")) #get day of year

mn_ndvi<-mean(NDVI$ndvi) #0.1850742
sd_ndvi<-sd(NDVI$ndvi) # 0.1131309

tabMonth<-NDVI%>% group_by(monthF) %>% 
  	summarize(mean = mean(ndvi))
dfM<-data.frame(tabMonth)
dfM$Zscore <-( dfM$mean - mn_ndvi ) /sd_ndvi
#use dfM to determine green and dry seasons, by monthly averages

#add scaled z-score value to NDVI dataframe
NDVI$Z<-( NDVI$ndvi - mn_ndvi ) /sd_ndvi


library(viridis)
ggplot()+
	geom_line(data= NDVI, aes(x=JulDate, y=Z, color=as.character(year)),
					 size=1, alpha=0.3)+
	scale_color_viridis(discrete = TRUE, begin =.5, end= .8) +
	theme_classic() +
	theme(legend.position = "none")
head(NDVI)



#For plotting NDVI 
tabM<-NDVI%>% group_by(monthF) %>% 
  	summarize(mean = mean(Z), sd=sd(Z), n = n())
tabM$SE<- tabM$sd/sqrt(tabM$n)
tabM$ord<-c(seq(2,12,1),1)
tabM2<-data.frame(tabM)
veg<-ggplot(tabM2, aes(x=ord, y=mean))+
	geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), fill = "olivedrab2", alpha=0.4) +
	geom_line(color="darkgreen", size=1)+
	labs(y= "Scaled NDVI", x = "") +
	scale_x_continuous(breaks=seq(1, 12, 1), limits=c(1, 12))+
	scale_y_continuous(position = "right")+
	theme_bw()+
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		axis.text.x=element_blank())

#############Precipitation
precip$date<-as.Date(precip$date, format= "%m/%d/%Y")
precip$month<-month(precip$date)
precip$year<-year(precip$date)
precip$monthF <- factor(format(precip$date, "%b"), month.abb, ordered = TRUE)
precip$monthF2<-factor(precip$monthF , 
	levels= c( "Dec","Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
			"Aug", "Sep", "Oct", "Nov") )
precip2<-na.omit(precip)#one NA in dec caused problems

grouped <- group_by(precip2, year, monthF2)
tabP2<-summarise(grouped, sum = sum(pcp_mm_chirps))

#summary of monthly average summed precip
tabP3<-tabP2%>%  group_by(monthF2) %>% 
  	summarize(mean = mean(sum))

#########Plot total rain by month
rain<-ggplot(tabP2, aes(x = monthF2, y= sum)) + 
    	stat_summary(fun = mean, geom = "bar", color= "navyblue", fill= "steelblue3") + 
    	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = .2)+ #95%CI
	labs(y= "Total precipitation (mm)", x = "Month") +
		theme_bw()+
		theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

pdf("ndvi_rain_draft.pdf" , width=8, height=6)
cowplot::plot_grid( rain , 
			veg,
                  ncol = 1,
                  labels = "auto")
dev.off()

#Combine to multipanel fig in Illustrator

#use NVDI z scores for each date to assign 
#resource values to specific sampling dates
write.csv(NVDI, "WW_scaled_NVDI.csv")

