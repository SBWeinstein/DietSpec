#test for individual spec in WW repeated surveys
#test if mean inter-individual variation (Psi) is sig less than expected
#from a Monte Carlo sim based null model
#build null model for proportional diet data

library("ggplot2")
library("dplyr")
library(phyloseq)
library("RInSp") #proportional similiarity, get average for each pop 

setwd("C:/Location/here") 

#file of just white water animals produced by White_water_multiyear.R script
psWW<-readRDS(file = "psWW_27Dec23.rds") 

#count of animals per sampling date
df<-data.frame(sample_data(psWW))
surv<-df %>% count(Date) #13 "pops"

###########################
######################################

#set up empty data frame for summary info for all pops
#set up empty list for list of mc vecofVs
list_mcVs <- list()
popVs<-data.frame(matrix(ncol=9,nrow=0, 
		dimnames=list(NULL, c("Date","n", "meanR", "totR", "TNW", "popV", "mean_mcV", "p", "sig"))))


for (pop in 1:nrow(surv)) {
	popVs[pop,1]<-surv$Date[pop]
	popVs[pop,2]<-surv$n[pop]

	#subset to one pop
	ps1pop<-subset_samples(psWW,Date == surv$Date[pop])
	ps1pop<-prune_taxa(taxa_sums(ps1pop) > 0, ps1pop)#reduce to plants there

	#for each animal in pop, get richness
	rich1<-estimate_richness(ps1pop, measures = c("Observed"))
	meanR<-mean(rich1$Observed) #mean richness
	popVs[pop,3]<-meanR
	
	#for the plants eaten, get relative abundance in pop diet
	psT<-merge_samples(ps1pop,"Date")
	popdietRS<-estimate_richness(psT, measures = c("Observed", "Shannon"))# pop level observed richness and shannon 
	popVs[pop,4]<-popdietRS$Observed
	popVs[pop,5]<-popdietRS$Shannon

	psTr  = transform_sample_counts(psT, function(x) x / sum(x) ) 
	OTU = data.frame(as(otu_table(psTr), "matrix"))
	propT<-c(OTU)

	Totu<-t(OTU)
	colnames(Totu)<-c("amount")
	plants<-c(colnames(OTU))  # list of plant taxa in pop

	## loop for the Null model:do it 1000 times for each pop
	#make an empty list for pop diet matrices and empty vector for pop var from RInSp
	listOfMCdiets <- list()
	vecofVs<- character()

		for (j in 1:1000){

		#Make empty dataframe for 1 MC pop diet, row names are plants eaten by pop 
		MCdiet<-data.frame(matrix(nrow = length(plants), ncol = 0))
		row.names(MCdiet) <- plants

		#for a population, 
		#loop through each animal and assign MC diet based on real diet richness
			for (i in 1:length(rich1$Observed)){
			pickplants<-sample(x=plants, size=rich1$Observed[i], prob=propT, replace=FALSE)
			picked<-data.frame(subset(Totu, (row.names(Totu) %in% pickplants)))
			##possible to never get one item, consider forcing all items, then rando the rest
			pickprops<-sample( x= row.names(picked), size=100, prob= picked$amount, replace=TRUE)
			diet1<-table(pickprops) #counts number of each item
			diet2<- data.frame(diet1)
			rownames(diet2) <- diet2[,1]  #
			diet2[,1] <- NULL
			colnames(diet2)[1]<-paste("Ind",dim(MCdiet)[2] +1 , sep = "")  #name individuals

			MCdiet<- merge(MCdiet, diet2 , by = 0, all=TRUE)
			MCdiet<- data.frame(MCdiet, row.names = 1)
			MCdiet[is.na(MCdiet)] <- 0 #replace NAs with 0
			#print(rich1$Observed[i])#testing loop
		}

		listOfMCdiets[[j]] <- MCdiet
		#get PSI and IS
		MCdietT<-t(MCdiet) #transpose; expects animals as rows, diet items as cols
		PCin<-import.RInSp(filename=MCdietT, col.header=TRUE, row.names = 0, data.type= "integer",
             print.messages=TRUE)
		PC_PSI<-PSicalc(PCin, pop.diet = "average", exclude = FALSE, replicates=3, precision = 1e-9)
		PC_PSI$IS  #get IS value 
		Var<-1-PC_PSI$IS #get V
		vecofVs[j]<-Var

	}
	list_mcVs[[pop]]<-vecofVs

	#test:is the Variance calculated for the real pop greater than seen in 1000 MC simulations?
	#P-value is the relative ranking of the test statistic among the sample values from the MC simulation

	#get the IS and Var from the real pop diets
	#ps1pop  the one population phyloseq object
	PCotu<-data.frame(otu_table(ps1pop))
	PCin<-import.RInSp(filename=PCotu, col.header=TRUE, row.names = 0, data.type= "integer",
             print.messages=TRUE)
	PC_PSI<-PSicalc(PCin, pop.diet = "average", exclude = FALSE, replicates=3, precision = 1e-9)
	PC_PSI$IS  #get IS value
	Var<-1-PC_PSI$IS #get real pop V
	popVs[pop,6]<-Var
	popVs[pop,7]<- mean(as.numeric(vecofVs))

	#P-value:the probability of obtaining the value of test statistic equal or more extreme than the observed one,
	# even if the null hypothesis is true. 
	# calculated as (the number of permutations where the value of test statistic is = observed test statistic + 1)
	# divided by (total number of permutations + 1). 
	#The +1 in the formula represents the observed value of the test statistic (on unpermuted data), 
	#which was added into the null distribution. 
	#Calculated P-value is then compared with a set of arbitrary thresholds (eg 0.05), 
	#to decide if significant at given level (e.g. P < 0.05), 

	out = ifelse(vecofVs< Var, "lessVar", "greaterVar")  #lessVar when MC V is smaller than real V

	#get number of reps when real pop exhibits more specialization than MC pop
	GV<-sum(out == "greaterVar")

	p<- (GV + 1)/(1000+1)  #
	popVs[pop,8]<-p
	sig <- ifelse(p < 0.05, "sig", "ns")
	popVs[pop,9]<-sig
}

popVs
saveRDS(popVs, "popVs_WW_11Jan24.rds") #
popVs<-readRDS( "popVs_WW_11Jan24.rds")

saveRDS(list_mcVs, "list_mcVs_11Jan24.rds")

plot(popVs$popV~popVs$TNW)
plot(popVs$mean_mcV~popVs$TNW)
####
#Figures
#adjust p-values for multiple comparisons
popVs$p.adj<-p.adjust(popVs$p, method = "BH")
popVs$adj.sig<- ifelse(popVs$p.adj>0.05, "ns", "sig")
popVs$Date<-surv$Date
surv<-df %>% count(Date)
dfGI<-Surveys_df[,c(1,4,5)] #this data frame is in wwseason code
dfGI$Date<-as.Date(dfGI$Date)#lost its date-y ness
df3<-merge(popVs, dfGI, by= "Date")
df3$Season<-as.character(df3$Season)

identical(dfGI$Date, popVs$Date)

pdf("WWnull_14feb24.pdf", width=3, height=3)
ggplot(df3, aes(x=TNW, y=popV)) +
	geom_point(data=df3, aes(x=TNW, y=mean_mcV, size=n, color=Season), shape=21)+
	geom_point(data=df3, aes(x=TNW, y=popV, size=n, fill=Season), shape =21)+
  	stat_smooth(method = "lm", 
              formula = y ~ x , se = FALSE, color= "black")+ 
	stat_smooth(data=popVs, aes(x=TNW, mean_mcV), 
		method = "lm", formula = y ~ x , se = FALSE, 
		color= "dark grey", linetype = "dashed")+
	scale_color_manual(values=c("#076612", "#D1B305", "#914223"))+
	scale_fill_manual(values=c("#076612", "#D1B305", "#914223"))+
	scale_size(range = c(2,5)) +
	theme_classic() +
	theme(legend.position = "none")
dev.off()

  	
# do the MC and real V values have diff slopes?
#gather dataframe; one column with popV and  mean_mcV valuep adjust , one with type of V
library(tidyr)
pop_reshape <- gather(popVs,variable, value, popV , mean_mcV)
mod1 <- aov(value~TNW*variable, data=pop_reshape) #see a sig interaction effect TNW:variable p<0.05
mod2 <- aov(value~TNW+variable, data=pop_reshape) #refit without interaction
anova(mod2,mod1)# removing interaction effect significantly reduces fit of model

#get equations for null and emperical models
emp<-lm(popVs$popV ~ popVs$TNW)
nul<-lm(popVs$mean_mcV ~ popVs$TNW)
summary(nul)
summary(emp)
