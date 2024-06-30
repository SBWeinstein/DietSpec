#Mark recapture data from individuals resampled at White Water

library("phyloseq")
library("ggplot2")
library(reshape2) #melt function
library(glmmTMB)#for count data regression
library(dplyr) #summary info

## see leav_out script for additional code for diet component IS impacts

setwd("C:/Location/here") 

psRF<-readRDS(file = "ps_rarF_Bigtrnl_29Nov23.rds") 

#subset to only White Water
psWW<-subset_samples(psRF,Site == "White Water")

reps<-data.frame(sample_data(psWW))%>% count(Animal_ID) 

#list of animal IDs with n>2
reps3<-reps[ which(reps$n >2 ),]
dim(reps3)#27 rats
r3rats<-reps3$Animal_ID

#subset to animals with at least 3 samples
psWW3<-subset_samples(psWW,Animal_ID %in% r3rats)
data.frame(sample_data(psWW3))%>% count(Animal_ID)#check
psWW3b<-prune_taxa(taxa_sums(psWW3) > 0, psWW3)#get rid of taxa not in this dataset

######################
##basic stats about diet diversity per sample
dfWW3b<-data.frame(estimate_richness(psWW3b, measures = c("Observed")),
	            Animal_ID=sample_data(psWW3b)$Animal_ID,
	            Date=sample_data(psWW3b)$Date,
			Names=row.names(sample_data(psWW3b)))

mean(dfWW3b$Observed) 
sd(dfWW3b$Observed) 

#average diet diversity for each individual and prop of total diet diversity
df3<-dfWW3b%>%
  	group_by(Animal_ID) %>%
  	summarize(mean_Obs = mean(Observed), sd_Obs= sd(Observed))

#total diet diversity for each individual
#merge by animal ID, reads already equal between samples
mergID = merge_samples(psWW3b, "Animal_ID")#turns ID into rownames
head(sample_data(mergID))
df4<-data.frame(estimate_richness(mergID, measures = c("Observed")),
	            Names=row.names(sample_data(mergID))	)
df3df<-data.frame(df3)
df34<-merge(df3df, df4, by.x = "Animal_ID", by.y = 0)
df34$propObs<-df34$mean_Obs/df34$Observed

mean(df34$propObs) #0.57 proportion of total diversity
sd(df34$propObs) #0.12

mean(df34$Observed) #Total diet diversity 8.29
sd(df34$Observed) #2.32

#confirm total diet stabilizes w/ three samples
df34r<-merge(df34, reps3, by.y = "Animal_ID", by.x = "Names")

plot(df34r$Observed~df34r$n)

#count data
obs_glm<-glm(Observed~n,family=poisson, data=df34r)
summary(obs_glm) #not significant increase with reps 

#add plot of this to supplementary material
ggplot(df34r, aes(x=n, y=Observed ))+
	geom_smooth(method = glm, formula = y ~ x, color = "grey25",
              method.args = list(family = poisson)) +
	geom_point(size = 2)+	
	xlab("Number of samples from different surveys ")+
	ylab("Total individual diet richness")+
	#stat_regline_equation(label.x=30, label.y=310)+
 	theme_bw()

#############################################
#for all pairs, calculate bray-curtis (phyloseq, 0 is total overlap,
#note that 1-BC= Quantitative sorensen index and now 1 is complete overlap)
Bray<-distance(psWW3b, "bray") # Bray-Curtis, output is a dist object
BM<-as.matrix(Bray)  #converts to a matrix
BM[lower.tri(BM, diag=TRUE)] <- NA  #keeps only upper triangle, no diagonal
df <- melt(BM, varnames = c("samp1", "samp2"),value.name = "Bray", na.rm = T)

#Jaccard distance (phyloseq, jaccard distance is complementary to jaccard coefficient JD=1-JC. 
#for JC no overlap is 0. for JD total overlap is 0.
Jac<-distance(psWW3b, "jaccard", binary = TRUE) # vegdist jaccard
JM<-as.matrix(Jac)  #converts to a matrix
JM[lower.tri(JM, diag=TRUE)] <- NA  #keeps only upper triangle, no diagonal
df2 <- melt(JM, varnames = c("samp1", "samp2"), value.name = "Jac", na.rm = T)

######################

identical(df2$samp2,df$samp2)#identical lists, allows lazy data merge
identical(df2$samp1,df$samp1)

df4<-cbind(df, JacD=df2$Jac)

#sanity check-expected relationships
plot(df4$Bray~df4$JacD)
lm(df4$Bray~df4$Jac)

#Jac distance and Bray-curtis are dissimiarities, could convert to similarity indeces
#df4$QS<-(1-df4$Bray) #Quantitative Sorensen Index
#plot(df4$PS~df4$QS) #strongly correlated, some are identical, but not all
#df4$JacI<-(1-df4$JacD) #Jaccard Index, aka Jaccard similarity coefficient

#add columns for categorical self/non-self and time between samples
#test<-tail(df4)
md<-data.frame(sample_data(psWW3b))
md1<-md[,c("Animal_ID","Date")] #get animal Id,date for each sample (and rowname)

#add columns for animals IDs and sampling date for samp1 and samp2
t1<-merge(df4, md1, by.x = "samp1", by.y = 0)
t2<-merge(t1, md1, by.x = "samp2", by.y = 0)
names(t2)[names(t2) == "Animal_ID.x"] <- "Animal_ID.S1"
names(t2)[names(t2) == "Animal_ID.y"] <- "Animal_ID.S2"
names(t2)[names(t2) == "Date.x"] <- "Date.S1"
names(t2)[names(t2) == "Date.y"] <- "Date.S2"

#add column, if s1_animal=s2_animal, self; else not self (ns)
t2$self<-ifelse(t2$Animal_ID.S1 == t2$Animal_ID.S2, "self","ns")

# get days between sampling days, not including years.
#function for getting days between samples, ignoring year difference

d_diff <- function(x, y) {
    	x1 <- as.integer(format(x, "%j")) #get day of year
	y1 <- as.integer(format(y, "%j"))
   	v1 <-abs(x1-y1)  #absolute diff between dates
	v2 <-(365-max(x1,y1))+ min(x1,y1) #but end of year,eg jan is close to december
	min(v1,v2) #pick smaller of those two
	}


#make sure columns are recognized as dates
t2$Date.S1<-as.Date(t2$Date.S1, format = "%m/%d/%y")
t2$Date.S2<-as.Date(t2$Date.S2, format = "%m/%d/%y")

#I give up, the loop works, don't judge me. get days between each day of the year
for(i in 1:dim(t2)[1]) {
       t2$times[i]<- d_diff(t2$Date.S1[i], t2$Date.S2[i])
	}

t2[100:130, ]
t2$self <- as.factor(t2$self)

#Does self and/or time predict similarity?
#https://cran.r-project.org/web/packages/betareg/vignettes/betareg.pdf
#https://cran.r-project.org/web/packages/betareg/betareg.pdf
#use beta regression as similarity is a proportion 
# but there are a few 0 and 1 values
# if y has 0 and 1, transformation from Smithson and Verkuilen 2006 
#(y * (n − 1) + 0.5)/n    where n is the sample size: (dim(t2)[1]

new.y <- function(a) {
   result <- (a *(dim(t2)[1]-1) + 0.5)/dim(t2)[1]
   print(result)
}
#new.y(1)

library("betareg")

Bn<-new.y(t2$Bray)
Jacn<-new.y(t2$JacD)

t2$Jacn<-Jacn

WW_logit <- betareg(Bn ~ self + times, data = t2)
WW_logit2 <- betareg(Jacn~ self + times, data = t2)

summary(WW_logit)#ns
summary(WW_logit2)#sig  

#is this influenced by there being a lot more ns than ss comparisons?
t2 %>% count(self)  #7029 ns, 231 self

tself<-t2[which(t2$self=="self"),]
tns<- t2[which(t2$self=="ns"),]

#for tns, pick 231 rownames to keep
#get p-value for 1000 replicates

vecp<-numeric()
for (x in 1:1000){
	keeps<-sample(1:7029, 231, replace=F)
	tns2<-tns %>% slice(keeps)
	ts3<-rbind(tself, tns2)
	WW_logit3 <- betareg(Jacn~ self + times, data = ts3)
	pval<-summary(WW_logit3)$coefficients$mean[2,4] #for this model
	vecp<-c(vecp, pval)
	}

vecp #vector of all of the p values (1000)
length(vecp[vecp< 0.05])  #how many significant? 549/1000, so still sig in majority of subsamples


par(mfrow = c(2, 2))
plot(WW_logit2, which = 1:4, type = "pearson")

p<-ggplot(t2, aes(x=times, y=JacD, color=self, ))+
	#geom_point(aes(alpha=self), show.legend = F)+
	geom_line(aes(y = predict(WW_logit2, t2)), linewidth=1.5)+  #default link is logit
	xlab("Days between samples")+
	ylab("Jaccard Dissimilarity")+
	scale_y_continuous(position = "right")+
	labs(color = "")+
 	theme_bw() 
pdf("mark_recap_self.pdf", width=3.5, height=4)
p
dev.off()

# leave one out analysis:
#seperate script, uses outputs and packages loaded here
#using jaccard diss- which plants are driving significant self-self pattern

###################################


