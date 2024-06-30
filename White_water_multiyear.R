#Multi-year sampling of White Water population

library("phyloseq")
library("ggplot2")
library("vegan")#for adonis
library(dplyr)

library("betareg")
library(multcomp) #for posthoc comparison
library("RInSp") #PSI; for diet overlaps

setwd("C:/location/here") #laptop

psRF<-readRDS(file = "ps_rarF_Bigtrnl_29Nov23.rds") 

#subset to only White Water
psWW<-subset_samples(psRF,Site == "White Water")
psWW<-prune_taxa(taxa_sums(psWW) > 0, psWW)#get rid of taxa not in this dataset

########################################################
#Is there a difference between species? (Sup mat)
#ancestry is continuous
#psWWLB<-subset_samples(psWW,Phylo_Code != "")#remove some unknown ancestry
psWWLB<-subset_samples(psWW,ancestry != "") #162 animals

orduB2 = ordinate(psWWLB, "PCoA", "bray")
evals <- orduB2$values$Eigenvalues
p3<-plot_ordination(psWWLB, orduB2, color="ancestry")
p3+theme_classic()+
	ggtitle("Diet by ancestry (Bray-Curtis Dissimilarity)")+
	scale_color_gradient(low="blue", high="red")+
  	coord_fixed(sqrt(evals[2] / evals[1]))

# Calculate bray curtis distance matrix
sp_bray <- phyloseq::distance(psWWLB, method = "bray")

# make a data frame from the sample_data
sampledf2 <- data.frame(sample_data(psWWLB))

# Adonis test: Permutational Multivariate Analysis of Variance Using Distance Matrices
#adonis2(sp_bray ~ sample_data(psWWLB)$Phylo_Code, data = sampledf2) #p=0.451
adonis2(sp_bray ~ sample_data(psWWLB)$ancestry, data = sampledf2) #p=0.276


sp_jac <- phyloseq::distance(psWWLB, method = "jaccard")
#adonis2(sp_jac ~ sample_data(psWWLB)$Phylo_Code, data = sampledf2) #p=0.446
adonis2(sp_jac ~ sample_data(psWWLB)$ancestry, data = sampledf2) #p=0.34

#####no diff between species, proceed as one pop
#######################################################
#look at diet patterns by season
#order months as factor for plots
Months<- c("January", "February", "March", 
		"April", "May", "June", "July", "August", "September", 
		"October","November", "December")
sample_data(psWW)$Month<- factor(sample_data(psWW)$Month, levels = Months)
	
#Use Months Jan/Feb/March (rainy),May/June/July (early dry),  August/Sep/October/November (late dry)
#this leaves months (april, dec) between each "season" where diets may be transitioning

#feeding patterns differ between seasons
#Diet comp (BC, J) and diversity differs between season (D, R, PR), justifying splitting
#include august in late-dry season
#add extra column to meta_data for season
df<-sample_data(psWW)
df$Season<- ifelse(df$Month == "July" | df$Month == "June" | df$Month == "May" , "ED","")
df$Season<- ifelse(df$Month == "January" | df$Month == "February" | df$Month == "March" , "W",df$Season)
df$Season<- ifelse(df$Month == "October" | df$Month == "November" | df$Month == "September"| df$Month == "August", "LD",df$Season)

sample_data(psWW)<-df  #replace sample data

data.frame(data.frame(df) %>% count(Month)) #no april or dec

sample_data(psWW)$Season<- factor(sample_data(psWW)$Season, 
					levels = c("W", "ED", "LD"))

#add column for plants in each sample
rich<-estimate_richness(psWW, measures = c("Observed"))
sample_data(psWW)$Rich<-rich$Observed

###################
orduB = ordinate(psWW, "PCoA", "bray")
evals <- orduB$values$Eigenvalues
p2<-plot_ordination(psWW, orduB, color="Season")
#BC plot for all three seasons

pdf("dietPCOA2.pdf", width=4, height=4)
p2+theme_classic()+
	geom_point(mapping = aes(size = Rich), alpha=0.5) +
	scale_color_manual(values=c("#076612", "#D1B305", "#914223"))+
	ggtitle("Diet by seasons (BC)")+
  	coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

orduJ = ordinate(psWWs, "PCoA", "jaccard")
p1<-plot_ordination(psWWs, orduJ, color="Season")

#Do the three seasons have different centroids?
# Calculate bray curtis (bray) distance matrix or jaccard (jaccard)
se_bray <- phyloseq::distance(psWWs, method = "jaccard")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psWWs))

# Adonis test: Permutational Multivariate Analysis of Variance Using Distance Matrices
adonis2(se_bray ~ sample_data(psWWs)$Season, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(se_bray, sampledf$Season)
permutest(beta) #sig- groups have diff dispersions
################
#Do early dry and late dry differ?- yes, more subtle than Rainy season dif, but sig difference (including aug as late dry)
psPRD<-subset_samples(psWW, Season != "W") #130 samples

orduB = ordinate(psPRD, "PCoA", "bray")
evals <- orduB$values$Eigenvalues
p2<-plot_ordination(psPRD, orduB, color="Season")

p2+theme_classic()+
	scale_color_manual(values=c( "#D1B305", "#914223"))+
	ggtitle("Diet by seasons (BC)")+
  	coord_fixed(sqrt(evals[2] / evals[1]))

#Do the two  seasons have different centroids?
# Calculate bray curtis (bray) distance matrix or jaccard (jaccard)
se_bray <- phyloseq::distance(psPRD, method = "bray")
se_bray <- phyloseq::distance(psPRD, method = "jaccard")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psPRD))

# Adonis test: Permutational Multivariate Analysis of Variance Using Distance Matrices
adonis2(se_bray ~ sample_data(psPRD)$Season, data = sampledf)#sig p =  0.001 (also sig for Jac)

# Homogeneity of dispersion test
beta <- betadisper(se_bray, sampledf$Season)
permutest(beta) #NS- p = 0.097 groups do not have diff dispersions
#########################################################
#get creosote consumed by each animal, compare relative amounts between seasons
psWWr= transform_sample_counts(psWW, function(x) x/sum(x))
psWWrZ<- subset_taxa(psWWr, Family=="Zygophyllaceae")
rZ<-data.frame(otu_table(psWWrZ))
colnames(rZ)[1]<-"rZyg"
sample_data(psWW)$rZyg<-rZ$rZyg
Zdf<-data.frame(sample_data(psWW))

sumZ<-Zdf%>% group_by(Season) %>% 
  	summarize(mean = mean(rZyg), sd=sd(rZyg), max=max(rZyg), min=min(rZyg), n=n())

model <- aov(rZyg~Season, data=Zdf)
summary(model)
TukeyHSD(model, conf.level=.95)

#how much Krameriaceae consumed each season?
psWWrK<- subset_taxa(psWWr, Family=="Krameriaceae")
rK<-data.frame(otu_table(psWWrK))
colnames(rK)[1]<-"rKram"

Kdf<-data.frame(sample_data(psWW),rKram= rK$rKram )

sumK<-Kdf%>% group_by(Season) %>% 
  	summarize(mean = mean(rKram), sd=sd(rKram), max=max(rKram), min=min(rKram), n=n())


################
#make Euler plot (sup mat)
merg1 = merge_samples(psWW, "Season") #Just want family lists for seasons

OT<-data.frame(otu_table(merg1 ))
TOT<-data.frame(t(OT)) #transpose for euler package
TOT2<-TOT %>% mutate_all(as.logical) #uses dplyr

library(eulerr)
fit2 <- euler(TOT2)
pdf("wwDieteuler.pdf", width=4, height=4)
plot(fit2, quantities = TRUE, fill= c("#076612", "#D1B305", "#914223")) 
dev.off()

#######################################
#look at animal diets (not included in manuscript/sup mat)
#make barplots of diet components relative abundance each season
#color by which seasons consumed "#076612", "#D1B305", "#914223")) 

merg1r= transform_sample_counts(merg1 , function(x) x/sum(x))
MOT<-data.frame(otu_table(merg1r ))
MTOT<-data.frame(t(MOT)) #transpose for euler package
tt<-data.frame(tax_table(merg1r))
MTOT$Fam<-tt$Family

MTOT <- MTOT[order(-MTOT$R),]

sub<-MTOT[1:15,] #use this to force plots to same number of bars, bar size
sub$Fam<-factor(sub$Fam, levels= sub$Fam)
Rbar<-ggplot(sub, aes(x=Fam, y=R)) + 
 	geom_bar(stat = "identity", fill="#076612")+
	ylim(0, 0.45)+
	theme_classic()+
	theme(axis.title.x=element_blank(), axis.text.x=element_blank())
	#theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

MTOT <- MTOT[order(-MTOT$PR),]
sub<-MTOT[1:15,] #use this to force plots to same number of bars, bar size
sub$Fam<-factor(sub$Fam, levels= sub$Fam)
PRbar<-ggplot(sub, aes(x=Fam, y=PR)) + 
 	geom_bar(stat = "identity", fill="#D1B305")+
	ylim(0, 0.45)+
	theme_classic()+
	theme(axis.title.x=element_blank(), axis.text.x=element_blank())

MTOT <- MTOT[order(-MTOT$D),]
sub<-MTOT[1:15,] #use this to force plots to same number of bars, bar size
sub$Fam<-factor(sub$Fam, levels= sub$Fam)
Dbar<-ggplot(sub, aes(x=Fam, y=D)) + 
 	geom_bar(stat = "identity", fill="#914223")+
	ylim(0, 0.45)+
	theme_classic()+
	theme(axis.title.x=element_blank(), axis.text.x=element_blank())

library(cowplot)
pdf("season_bar_diets.pdf", width=3, height=6)
plot_grid(Rbar, Dbar, PRbar, ncol = 1, align = 'v')
dev.off()
#######################
#for each diet sample, get scaled NVDI closest to sampling date
#library(survival)
#install.packages("birk")
library(birk)
library(lubridate)
WW_GI<-read.csv( "WW_scaled_NVDI.csv")
WW_GI$date<-as.Date(WW_GI$date)
WW_GI$JulDate<- as.integer(format(WW_GI$date, "%j")) #get day of year
WW_GI$year<-year(WW_GI$date)

TdfWW <- data.frame(sample_data(psWW))
TdfWW$Date<-as.Date(TdfWW$Date, format= "%m/%d/%y")

#function to give NVDI value closest to each sampling date
func2 <- function(y) {
		x= WW_GI$date
		q= WW_GI$Z
		ind<-which.closest(x,y)
		q[ind]
		}
TdfWW$GI<-sapply(TdfWW$Date, func2) #Z=GI

sample_data(psWW)<-TdfWW  #replace sample data with NVDI included version

####plot
#use TdfWW to add points of when samples were taken, color by season?
subTdf<-TdfWW[, c("Date","Month","Year", "Season","GI")]

dfpoints<-subTdf%>% group_by(Date) %>% 
  		summarize(max = max(GI), n=n())
dfpoints<-data.frame(Date=dfpoints$Date, Z=dfpoints$max, n=dfpoints$n)
dfpoints$JulDate<- as.integer(format(dfpoints$Date, "%j")) #get day of year
dfpoints$Month<-month(dfpoints$Date)
dfpoints$year<-year(dfpoints$Date)
dfpoints$Season<- ifelse(dfpoints$Month <4  , "W","")
dfpoints$Season<- ifelse(dfpoints$Month >4 & dfpoints$Month < 9 , "ED",dfpoints$Season)
dfpoints$Season<- ifelse(dfpoints$Month >  9, "LD",dfpoints$Season)
dfpoints$Season<-factor(dfpoints$Season, levels= c("W", "ED", "LD"))

ML<-substr(month.name, start = 1, stop = 1) #month abbreviations

#### add stronger NVDI lines for years of sampling 2017-2021
GI_S<-subset(WW_GI, WW_GI$year > 2016 & WW_GI$year < 2022 )

library(ggnewscale)
#install.packages("ggnewscale")#allows two color scales in one plot
library(viridis)

###plot sampling trips by NDVI data, with 30 year NDVI data under
library(ggrepel)

pdf("sampling_NDVI.pdf", width=4, height=4)
ggplot()+
	geom_line(data= WW_GI, aes(x=JulDate, y=Z, color=as.character(year)),
					 size=1, alpha=0.3)+
	scale_color_viridis(discrete = TRUE, begin =.5, end= .8) +
	geom_hline(yintercept=0, linetype="dashed")+
	guides(color="none") +
	
	new_scale_colour()+
	geom_point(data=dfpoints, aes(x=JulDate, y=Z, size=n, fill =Season), shape=21)+
	scale_fill_manual(values=c("#076612",  "#D1B305","#914223" ))+
	scale_size(range = c(2,5.5)) +
	
	scale_x_continuous(breaks=seq(15,365,30), labels= ML)+
	theme_classic() +
	xlab(NULL)+
	ylab("Scaled NDVI")+
	guides(fill="none") +
	theme(legend.position = "top")
dev.off()

##################
#alternative labeled by year version for sup mat
ggplot()+
	scale_color_viridis(discrete = TRUE, begin =.3, end= .8) +
	geom_line(data= GI_S, aes(x=JulDate, y=Z, color=as.character(year)),
					 size=0.5, alpha=1, lty =1)+
	geom_point(data=dfpoints, aes(x=JulDate, y=Z, size=n, fill =Season), shape=21)+
	facet_grid(rows = vars(year))+
	scale_fill_manual(values=c("#076612",  "#D1B305","#914223" ))+
	scale_size(range = c(3,5.5)) +
	
	scale_x_continuous(breaks=seq(15,365,30), labels= ML)+
	geom_hline(yintercept=0, linetype="dashed")+
	theme_bw() +
	xlab(NULL)+
	ylab("Scaled NDVI")+
	guides( color = "none") +
	theme(legend.position = "right")

#####################
#save psWW updated version
saveRDS(psWW, file = "psWW_27Dec23.rds")

#pattern of creo consumption by ndvi?
row.names(Zdf)==row.names(TdfWW) #same order
ZGI<-data.frame(GI=TdfWW$GI, rZyg=Zdf$rZyg)
ZGI3<-ZGI[which(ZGI$GI<2.5),]

p<-ggplot(ZGI, aes(x=GI, y=rZyg ))+
	geom_point()+
	geom_line(aes(y = predict(zyg_logit, ZGI)), linewidth=1.5)+
	 theme_bw()

new.y <- function(a) {
   result <- (a *(dim(Zdf)[1]-1) + 0.5)/dim(Zdf)[1]
   print(result)
	}

Zsc<-new.y(ZGI3$rZyg)
zyg_logit <- betareg(Zsc~ ZGI3$GI)
summary(zyg_logit) #higher greenness, less zygo, but driven by low zygo in really wet 2019 march
plot(zyg_logit )#surprisingly good
#################################
#sanity check, NDVI values consistent with seasonal expectations
ggplot(data=TdfWW, aes(x=Date, y=GI, color= Season))+ 
		geom_jitter()+
		 scale_x_date(date_labels="%y-%m",date_breaks  ="1 month")+
		theme(axis.text.x = element_text(angle = 60, hjust = 1))
####################
#compare TNW and TR by NDVI using surveys as measure instead of season
psSurveys<-merge_samples(psWW, "Date")
psSurveys
Surveys_df<-data.frame(Date=row.names(sample_data(psSurveys)), 
			estimate_richness(psSurveys, measures = c("Observed", "Shannon")),
	           	Season=sample_data(psSurveys)$Season,
		    	GI=sample_data(psSurveys)$GI
			)	
#add in count of number of samples per survey, 
#these were calculated in dfpoints, in same order, just add col
Surveys_df$n<-dfpoints$n

#patterns due to high NDVI in Mar19? partially, TR is sig if outlier excluded
Surveys_df2<- Surveys_df[which(Surveys_df$GI<2.5),]

pdf("Obs_byGI_poplevel_legend.pdf", width=4, height=4)
ggplot(data=Surveys_df, aes(x=GI, y=Observed)) + 
	scale_color_manual(values=c("#076612", "#D1B305", "#914223"))+
	geom_point(data=Surveys_df, aes(x=GI, y=Observed, size=n,
				color=as.character(Season)), alpha=0.6)+
	geom_smooth(method = "glm", , se = T, 
        	method.args = list(family = "poisson"), 
		color="grey60", alpha=0.05,lty=2)+
	geom_smooth(data=Surveys_df2, aes(x=GI, y=Observed), method = "glm", , se = T, 
        	method.args = list(family = "poisson"), 
		color="grey40", alpha=0.05,lty=1)+
	scale_size(range = c(2,5)) +
	theme_classic()
	#theme(legend.position = "none")
dev.off()


obs_glm<-glm(Observed~GI+n,family=poisson, data=Surveys_df)#n does not improve model
obs_glm<-glm(Observed~GI+n,family=poisson, data=Surveys_df)
summary(obs_glm)

obs_glm<-glm(Observed~GI+n,family=poisson, data=Surveys_df2 )#n does not improve model
summary(obs_glm)
par(mfrow=c(2,2))
plot(obs_glm)

ggplot(data=Surveys_df, aes(x=GI, y=Shannon, color=Season)) + geom_point()
shan_lm<-lm(Shannon~GI,data=Surveys_df2)
summary(shan_lm)
plot(shan_lm)

###############################
#####calculate overlaps within trapping trip
survs<-data.frame(TdfWW  )%>% count(Date)#number of samps from each date (n=13), all >9 except 3/18/18

trips<-as.character(survs$Date)#[-3] getting rid of one with only 5 samples has no effect
#subset can't be used in loops, but prune can be

dfloop<-data.frame(Names=character(), PSI= as.numeric())
for (i in 1: 13) {
	date <- trips[i]
	remove_pops <- as.character(get_variable(psWW, "Date")) == trips[i]
	psSub<-prune_samples(remove_pops , psWW ) #
	Subotu<-data.frame(otu_table(psSub))
	Sub_in<-import.RInSp(filename=Subotu, col.header=TRUE, row.names = 0, data.type= "integer",
             print.messages=TRUE)
	Sub_PSI<-PSicalc(Sub_in, pop.diet = "sum", exclude = FALSE, replicates=99, precision = 1e-3)
	dfSub<-data.frame( Names=row.names(sample_data(psSub)), PSI=Sub_PSI$PSi)
	dfloop<-rbind(dfloop, dfSub)
	}

#makes dfloop, with rownames and within date sample PSI values

################
#get diet diversity metrics for individuals
WWdf<-data.frame(Names=row.names(sample_data(psWW)), 
			estimate_richness(psWW, measures = c("Observed", "Shannon")),
			Animal_ID=sample_data(psWW)$Animal_ID,
	           	Season=sample_data(psWW)$Season,
		    	GI=sample_data(psWW)$GI
			)	
WWdf2<-merge(WWdf, dfloop, by="Names")

#Pielou's evenness (aka Shannon or Shannon-Weaver/Wiener/Weiner evenness; H/ln(S). 
#H=shannon's diversity index, S is species richness
#species evenness ranges from zero (no evenness) to one (complete evenness)
WWdf2$PielE<-WWdf2$Shannon/log(WWdf2$Observed)

#####################################
#does mean diet richness increase with TNW? for surveys
Surveys_df #Shannon is the TNW of each survey 
df3<-data.frame(Names=row.names(sample_data(psWW)), 
			estimate_richness(psWW, measures = c("Observed")),
			Date=sample_data(psWW)$Date)

df4<-data.frame(df3 %>% group_by(Date) %>% 
				dplyr::summarize(Mean = mean(Observed), 
							sd = sd(Observed),
							n = n()))
str(df4)
Surveys_df$Date<-as.Date(Surveys_df$Date)

df5<-merge(Surveys_df, df4, by= "Date")
lm5<-lm(Shannon~Mean,data=df5)
par(mfrow=c(2,2))
plot(lm5) #check model diagnostics, looks ok

df5$Season<-as.factor(df5$Season)
df5$up<-df5$Mean+df5$sd
df5$down<-df5$Mean-df5$sd
pdf("MeanObs_byTNW.pdf", width=3, height=3)
ggplot() +
	geom_smooth(data=df5, aes(x=Shannon, y= Mean), method='lm', se=T, formula=y~x, color= "black")+
	scale_fill_manual(values=c("#076612", "#D1B305", "#914223"))+
	geom_errorbar(data=df5, aes(x=Shannon, ymin = down, ymax = up), color="dark grey")+
	geom_point(data=df5, aes(x=Shannon, y= Mean, fill=Season, size=n), shape=21)+
	scale_size(range = c(2,5)) +
	theme_classic()+
	xlab("TNW")+
	ylab("Mean diet richness")+
	theme(legend.position = "none")
dev.off()



#look at correlations between variables
WWdf3<-WWdf2[,-c(1,4,5, 6)] #keep obs, shan, Psi

library("Hmisc") #correlation matrix with p values
rcorr(as.matrix(WWdf3))

install.packages("corrplot")
library(corrplot)
mydata.cor = cor(WWdf3, method = c("spearman"))
corrplot(mydata.cor)  #some sig correlations (obs~shan, shan~ E, shan~PSI, PSI~E)
##################################################
#test individual diet patterns with resource abundance (GI)
#Richness
ggplot(data=WWdf2, aes(x=GI, y=Observed)) + geom_point()
obs_glm<-glm(Observed~GI,family=poisson, data=WWdf2)
summary(obs_glm)
par(mfrow=c(2,2))
plot(obs_glm) #looks pretty good even with Mar19

###make plot with regression line
newdat <- data.frame(GI = seq(-.75, 3.5, .1))
#pred <- predict(obs_glm, newdata = newdat, type = "response")
#plot(newdat$AGE, pred, type = "l")

ginv <- obs_glm$family$linkinv  ## inverse link function
prs <- predict(obs_glm, newdata = newdat, type = "link", se.fit=TRUE)
newdat$pred <- ginv(prs[[1]])
newdat$lo <- ginv(prs[[1]] - 1.96 * prs[[2]])
newdat$up <- ginv(prs[[1]] + 1.96 * prs[[2]])

pdf("Obs_byGI.pdf", width=3, height=3)
ggplot() +
	scale_color_manual(values=c("#076612", "#D1B305", "#914223"))+
	geom_ribbon(data=newdat, aes(x=GI, ymin=lo, ymax=up), alpha=.2)+
	geom_jitter(data=WWdf2, aes(x=GI, y=Observed, color=Season), width=.02)+
	geom_line(data=newdat, aes(x=GI, y=pred))+
	theme_classic()+
	theme(legend.position = "none")
dev.off()

#Shannon
ggplot(data=WWdf2, aes(x=GI, y=Shannon)) + geom_point()
shan_lm<-lm(Shannon~GI,data=WWdf2)
summary(shan_lm)
plot(shan_lm)#looks pretty good

##evenness
ggplot(data=WWdf2, aes(x=GI, y=PielE, color=Season)) + geom_point()
#PielE_lm<-lm(PielE~GI,data=WWdf2)
pi_logit <- betareg(PielE~ GI, data = WWdf2)
summary(pi_logit) #no sig effect of season on ind evenness
summary(PielE_lm)#no effect
#plot(PielE_lm)

#proportional overlap
ggplot(data=WWdf2, aes(x=GI, y=PSI, color = Season)) + geom_point()
#PSI_lm<-lm(PSI~GI,data=WWdf2)
psi_logit <- betareg(PSI~ GI, data = WWdf2)
summary(psi_logit)
plot(PSI_lm)

#Due to very high GI in one sampling trip? only marginally sig when this trip excluded
WWdf4<- WWdf2[which(WWdf2$GI< 2),]
ggplot(data=WWdf4, aes(x=GI, y=Observed)) + geom_point()
obs_glm<-glm(Observed~GI,family=poisson, data=WWdf4)
summary(obs_glm) #p = 0.0787  

PSI_lm<-lm(PSI~GI,data=WWdf4)
summary(PSI_lm) #p= 0.067 . 
###########################################
#go to diet_niche_density script to make geom_density niche plots
############################################
#######################################################
########################
##########

mean(R_PSI$PSi) #this is IS of the population, 
#proportional similarity always much lower than expected from null model (ie variance much higher than null)
#but is this a good null model? I don't think so.
R_PSI$PSi.montecarlo #gives matrix dim 48 (row=sample size) by 1001 (col=reps), 1st one is actual pop
R_sim_means<-apply(R_PSI$PSi.montecarlo,2,mean) #get column means (IS)
########################################################

#Season level measures
#merge each pop to get dominant plant, reads already equal between samples
merg = merge_samples(psWW, "Season")#turns seasons into rownames
sample_data(merg)$Season<-  factor(rownames(sample_data(merg)), 
					levels = c("R", "PR", "D"))

#total richness and niche width (shannon) in each season?
#total diversity of plants R:22, PR: 16, D:15
#TNW (shannon) by season R: 1.85, PR: 1.46, D 1.46
estimate_richness(merg, measures = c("Observed", "Shannon")) 

#what is the most commonly consumed plant in each season
#convert counts to relative abundance
psFr = transform_sample_counts(merg, function(x) x/sum(x))

plot_bar(psFr , x= "Season", fill="Family")+
	theme_classic() 

#get the most abundant plant OTU from each season
OT<-data.frame(otu_table(psFr))
plants<-colnames(OT)[max.col(OT,ties.method="first")]
#get max plant column index value (pick first if ties)
maxP<-max.col(OT, "first")
#get value from each of those columns
value <- OT[cbind(1:nrow(OT), maxP)]
res <- data.frame(Season=row.names(OT), plants, value)

#add plant family column to plant dataframe
tt<-tax_table(psFr)[,5]
res2<-merge(res, tt, by.x = "plants", by.y = 0)
#Most abundant plants; R: Zygo 32%, PR Kram 45%; D Zygo 50%
###############################################################

#get average, SD for individuals eating those plants each season
##add column with prop of each plant family to WWdf2
psWWprop = transform_sample_counts(psWW, function(x) x/sum(x))
ot2<-data.frame(otu_table(psWWprop))
tt2<-data.frame(tax_table(psWWprop))
colnames(ot2)<-tt2$Family

WWdf3<-merge(WWdf2, ot2, by.x = "Names", by.y=0)

#across year
AnMer = merge_samples(psWW, "Type")#turns into one sample
AnMerprop = transform_sample_counts(AnMer , function(x) x/sum(x))
ot2<-data.frame(otu_table(AnMerprop))
tt2<-data.frame(tax_table(AnMerprop))
colnames(ot2)<-tt2$Family

mean(WWdf3$Zygophyllaceae) 
sd(WWdf3$Zygophyllaceae)

mean(WWdf3$Krameriaceae) 
sd(WWdf3$Krameriaceae) 

#summary of individual diversity metrics by season

sumTab<-WWdf3%>%
  group_by(Season) %>%
  summarise(mean_zyg = mean(Zygophyllaceae), sd_zyg= sd(Zygophyllaceae),
		mean_kra = mean(Krameriaceae), sd_kra= sd(Krameriaceae),
		mean_obs = mean(Observed), sd_obs= sd(Observed),
		mean_shan = mean(Shannon), sd_shan= sd(Shannon),
		mean_PielE= mean(PielE), sd_PielE= sd(PielE),
		mean_PSI= mean(PSI), sd_PSI= sd(PSI)	)

#####################################################