#diet specialization across all 57 pops

library("ggplot2")
library("dplyr")
library(phyloseq)
library("vegan")#species accumulation curves
library(ggpubr) #multipanel figure
library(viridis)
library("RInSp") #for Ind spec, V calculations
library(tidyr)
library("betareg")

setwd("PATH HERE") 

psRF<-readRDS(file = "ps_rarF_Bigtrnl_29Nov23.rds") 

#count of animals per pop
df<-data.frame(sample_data(psRF))
pops<-df %>% count(Phylo_Code)

#subset to pops with at least 3 animals, and have a Phylo_code
pop3<-pops[which(pops$n>2 & pops$Phylo_Code != ""),]
p3<-pop3$Phylo_Code
ps3p<-subset_samples(psRF, (Phylo_Code %in% p3))

#remove samples from multiple sampling rounds at same site
#WW, remove all except one bryanti batch
dfb<-df[which(df$Month == "November" & df$Year== 2019 & df$Phylo_Code == "B17"),]
b17<-row.names(dfb)
dfWW<-df[which(df$Site == "White Water"),] #all the white water animals
WWnames<-row.names(dfWW)
wwremove<-setdiff(WWnames, b17)  #other white water L, LB, B to remove
ps4<-subset_samples(ps3p,!(sample_names(ps3p) %in% wwremove))#remove extra wwbryanti

####for floridana, use 2022 Key largo batch
dfb1<-df[which(df$County == "Monroe" & df$Year== 2017),]
F17<-row.names(dfb1) #samples to remove
ps5<-subset_samples(ps4,!(sample_names(ps4) %in% F17))

#check the list
df5<-data.frame(sample_data(ps5)) 
df5%>% count(Spec) #13 rat species
cts5<-df5%>% count(Spec, Phylo_Code) #57 pops

cts6<-cts5 %>% count(Spec) #number of pops of each species
arrange(cts6, -n) #show in desending n order

#average individuals per pop
mean(cts5$n)
sd(cts5$n)
##############
#remove 2 highly uncertain samples from S40 population
#likley did not receive expected samples in loan
S40<-subset_samples(ps5, Phylo_Code == "S40")
S40 <-prune_taxa(taxa_sums(S40 ) > 0, S40 )#
plot_bar(S40, fill="Family", x="Animal_ID")
S40x<-c("MSB:Mamm:322261", "MSB:Mamm:322268") #samples to remove
ps5b<-subset_samples(ps5,!(Animal_ID %in% S40x))

#### #merge at pop level
ps6<-prune_taxa(taxa_sums(ps5b) > 0, ps5b)
merg = merge_samples(ps6, "Phylo_Code") #reads already equal between samples, phylo codes are row names

########################
#version with just lepida (for supplemental figure)
psLep<-subset_samples(ps5,Spec == "Neotoma lepida")
#merge at pop level
psLep<-prune_taxa(taxa_sums(psLep) > 0, psLep)#45 plant families
merg = merge_samples(psLep, "Phylo_Code")
####################

#Make figure 2
#identify most abundant diet component per pop
psFr = transform_sample_counts(merg, function(x) x/sum(x)) #convert counts to relative abundance

#data frame w/ most abundant OTU per pop
OT<-data.frame(otu_table(psFr))
plant<-colnames(OT)[max.col(OT,ties.method="first")] #get the most abundant plant from each sample
maxP<-max.col(OT, "first") #get max plant column index value (pick first if ties)
value <- OT[cbind(1:nrow(OT), maxP)]  #get value from each of those columns
res <- data.frame(row.names(OT), plant, value)   #dataframe: pop with dominant plant and its proportion
names(res)[1]<-"Phylo_Code" #fix name
TT1<-data.frame(tax_table(psFr)) #

resTT<-merge(res, TT1, by.x="plant", by.y = 0)
resTT2<-resTT[, c(1,2,3,8)]  #dominant plants by pop, one metric of spec

resTT2[order(resTT2$value),]

spec_amts<- resTT2[order(resTT2$value),]
spec_amts<-spec_amts[,-1]
spec_amts$value<-round(spec_amts$value, digits = 2)
spec_amts$SG<-ifelse(spec_amts$value < 0.6, "G", "S")
spec_amts%>% count(SG)  #how many are spec v generalists

# for each rat, get the amount of its pop's most abundant plant
psF_r  = transform_sample_counts(ps6, function(x) x / sum(x) )  
# Extract OTU relative abundance matrix from the phyloseq object
OTUf = data.frame(as(otu_table(psF_r), "matrix"))
P_code<-data.frame(sample_data(psF_r)[,10])
PC_OTUS<-cbind(OTUf, P_code)

#for each animal, put in the dominant diet item for the pop
OTUf3<-merge(x=PC_OTUS, y=resTT2, by = "Phylo_Code")  #drop animal IDs at this point

#make new column that has rat's value for that dominant plant
dom.plant.val <- vector("numeric", nrow(OTUf3))
for (row in 1:nrow(OTUf3)) {
	y<-as.vector(OTUf3$plant[row])
	dom.plant.val[row]<- OTUf3[row,y]
		}

OTUf4<-data.frame(OTUf3$Phylo_Code, OTUf3$plant, OTUf3$Family, dom.plant.val, OTUf3$value)
names(OTUf4)<-c("Phylo_Code","seq","plant.family", "dom.plant.value", "pop.avg.val")

#order phylo codes based on dominant plant abundance (specialization)
domO<- resTT2[order(resTT2$value),]
level_order_dom <- as.vector(domO$Phylo_Code)
OTUf4$Phylo_Code<-factor(OTUf4$Phylo_Code, levels =level_order_dom)

p<-ggplot(OTUf4, aes(x=Phylo_Code, y=dom.plant.value, fill=plant.family)) + #rev() to reverse order
  geom_boxplot(outlier.shape = NA)+
  scale_fill_viridis(discrete = TRUE)+
  geom_jitter(shape=16, position=position_jitter(width = 0.1), alpha=0.6)+
	geom_hline(yintercept=.6, lty=6, alpha=.8)+
	theme_bw() +
  	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0), 
		panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
	labs(x ="Population", 
		y = "Proportion of diet from dominant plant")
pdf( "diet_spec_26Apr24.pdf", width=10, height=6,useDingbats=FALSE)
#pdf( "diet_spec_Lepida_27May24.pdf", width=10, height=6,useDingbats=FALSE)
p
dev.off()

#calculate PNW and TR for each pop
popsdf<-data.frame(estimate_richness(merg , measures = c("Observed", "Shannon")))

df10<-merge(popsdf, spec_amts, by.x=0, by.y="Phylo_Code")

##add in mean richness for each pop
indR<-data.frame(estimate_richness(ps6 , measures = c("Observed")))
indR$Phylo_Code<-sample_data(ps6)$Phylo_Code
MR<-data.frame(indR%>% group_by(Phylo_Code) %>% summarize(MeanR= mean(Observed)))

df11<-merge(df10, MR, by.x=1, by.y="Phylo_Code")

#maximum TR for generalist pops
gens<-df11[which(df11$SG== "G"),]
max(gens$Observed)
mean(gens$Observed)
sd(gens$Observed)

#average individual richness in generalists
max(gens$MeanR)
mean(gens$MeanR)
sd(gens$MeanR)

#do specialists have lower TR, TNW and MeanRich?
###total richness
ggplot(data=df11, aes(x=SG, y=Observed))+ 
		geom_boxplot(outlier.shape=NA)+
		geom_jitter( aes( color = value), width = .1)
m1<-glm(Observed~SG, family="poisson",data=df11)#
summary(m1) #sig
df11%>% group_by(SG) %>% summarize(Mean= mean(Observed),SD =sd(MeanR) )

#Mean Richness
ggplot(data=df11, aes(x=SG, y=MeanR))+ 
		geom_boxplot(outlier.shape=NA)+
		geom_jitter( aes( color = value), width = .1)
m1<-lm(df11$MeanR~df11$SG)  # sig p =   0.000187
aov1 <- aov((df11$MeanR~df11$SG))  #  0.000187 **same
par(mfrow = c(2, 2))
plot(aov1)
summary(m1)

#Total niche width
ggplot(data=df11, aes(x=SG, y=Shannon))+ 
		geom_boxplot(outlier.shape=NA)+
		geom_jitter( aes( color = value), width = .1)
m1<-lm(df11$Shannon~df11$SG)  # sig p = 4.5e-10  ***
summary(m1)

#how much does TNW vary across pops?
max(df11$Shannon)/min(df11$Shannon) #28  #nearly 30 fold change 

#### niche expansion 
#does mean richness increase with PNW?
m1<-lm(df11$MeanR~df11$Shannon)  # sig p = 1.95e-08 ***
summary(m1)

##for figure 3
df11$sp<-as.character(substr(df11$Row.names, start = 1, stop = 1))
#add sample size of each pop using cts5
df12<-merge(df11, cts5, by.x=1, by.y="Phylo_Code")

spec.col<-c("#F7C908", #albigula yellow  A
		"#1B0121", #bryanti near black B
		"#B72467", #cinerea maroon C
		"#BFA8E2", #devia lilac D
            	"#A05B00", #fuscipes brown orange E
		"#841618", #floridana dark red F
		"#BFA211",  #magister yellow brown G
		"#4D04BF", #lepida purple L
		"#FF0000", #macrotis red M
		"#F998E4", #micropus light pink P
		"#EF8630", #stephensi orange S
		"#6B4349", #leucodon brown U
		"#E205E2") #mexicana pinkish X

P1<-ggplot(data=df12, aes(y=MeanR, x=Shannon))+ 
	geom_smooth(method='lm', se=T, formula=y~x, color= "black")+
	geom_point( aes(fill = sp, size=n), shape =21, alpha=0.8)+
	scale_size(range = c(2,5)) +
	scale_fill_manual(values=spec.col, guide="none")+
	theme_classic() +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
	ylab("Average individual diet richness")
	
######################################
###testing for specialization using monte carlo resampling null model 
#set up empty data frame for summary info for all pops
#set up empty list for list of mc vecofVs

list_mcVs <- list()
popVs<-data.frame(matrix(ncol=9,nrow=0, 
		dimnames=list(NULL, c("Phylo_Code","n", "meanR", "totR", "TNW", "popV", "mean_mcV", "p", "sig"))))

for (pop in 1:nrow(cts5)) {
	popVs[pop,1]<-cts5$Phylo_Code[pop]
	popVs[pop,2]<-cts5$n[pop]

	#subset to one pop
	ps1pop<-subset_samples(ps6,Phylo_Code== cts5$Phylo_Code[pop])
	ps1pop<-prune_taxa(taxa_sums(ps1pop) > 0, ps1pop)#reduce to plants there

	#for each animal in pop, get richness
	rich1<-estimate_richness(ps1pop, measures = c("Observed"))
	meanR<-mean(rich1$Observed) #mean richness
	popVs[pop,3]<-meanR
	
	#for the plants eaten, get relative abundance in pop diet
	psT<-merge_samples(ps1pop,"Phylo_Code")
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
	#The + 1 in the formula represents the observed value of the test statistic (on unpermuted data), 
	#which was added into the null distribution. 
	#Calculated P-value is then compared with a set of arbitrary thresholds (eg 0.05), 
	#to decide if significant at given level (e.g. P < 0.05), 

	out = ifelse(vecofVs< Var, "lessVar", "greaterVar")  #lessVar when MC V is smaller than real V

	#get number of reps when real pop exhibits more specialization than MC pop
	GV<-sum(out == "greaterVar")

	p<- (GV + 1)/(1000+1)  
	popVs[pop,8]<-p
	sig <- ifelse(p < 0.05, "sig", "ns")
	popVs[pop,9]<-sig
}


saveRDS(popVs, "popVs_26Apr24.rds") #
saveRDS(list_mcVs, "list_mcVs_26Apr24.rds")

plot(popVs$popV~popVs$TNW)
plot(popVs$mean_mcV~popVs$TNW)
####

popVs<-readRDS( "popVs_26Apr24.rds")
#Figures
#adjust p-values for multiple comparisons
popVs$p.adj<-p.adjust(popVs$p, method = "BH")
popVs$adj.sig<- ifelse(popVs$p.adj>0.05, "ns", "sig")

popVs$sp<-as.character(substr(popVs$Phylo_Code, start = 1, stop = 1))

#how many are sig v nonsig
popVs%>% count(adj.sig)

#who are the pops w/o IS?
ISpops<-data.frame(popVs[,c(1,10,11)])
df12<-merge(df11, ISpops, by.x=1, by.y=1)
arrange(df12, adj.sig) #they are the most specialized pops

#average amount of diet dominant plant in these pops
mean(0.76,0.99,0.94,0.97,0.96, 0.76)

#########Plot for figure three
P2<-ggplot(popVs, aes(x=TNW, y=popV)) +
	geom_point(data=popVs, aes(x=TNW, y=mean_mcV, size=n, color= sp), shape=21)+  #simulated
	geom_point(data=popVs, aes(x=TNW, y=popV, size=n,  fill=sp, shape =adj.sig), alpha =.8)+  #data
	scale_size(range = c(2,5)) +
	scale_color_manual(values=spec.col, guide="none")+
	scale_fill_manual(values=spec.col, guide="none")+
	 scale_shape_manual(values=22:21)+
  	stat_smooth(method = "lm", 
              formula = y ~ x , se = FALSE, color= "black")+ 
	stat_smooth(data=popVs, aes(x=TNW, mean_mcV), 
		method = "lm", formula = y ~ x , se = FALSE, 
		color= "dark grey", linetype = "dashed")+
	theme_classic() +
	ylab("Inter-individual diet variation (V")

pdf("allpopsV_TNW_26Apr24.pdf", width=4, height=7)
ggarrange(P1, P2, nrow = 2, align = c("v"), common.legend = TRUE)
dev.off()

####### ######
#########################################
# do the MC and real V values have diff slopes?
#gather dataframe; one column with popV and  mean_mcV value, one with type of V
pop_reshape <- gather(popVs,variable, value, popV , mean_mcV)
mod1 <- aov(value~TNW*variable, data=pop_reshape) #see a sig interaction effect TNW:variable p<0.05
mod2 <- aov(value~TNW+variable, data=pop_reshape) #refit without interaction
anova(mod2,mod1)# removing interaction effect significantly reduces fit of model
#slopes are significantly different

#get equations for null and emperical models
emp<-lm(popVs$popV ~ popVs$TNW) #y=0.272x+0.016
nul<-lm(popVs$mean_mcV ~ popVs$TNW)#y=0.135x-0.018
summary(emp)

#####################
##test whether mean richness correlates with V
#include in supplemental material

PR_logit <- betareg(popV~ meanR, data = popVs)
PR_lm<-lm(meanR~popV, data=popVs)#same results, significant positive correlation

# predict with new data
dfPR <- data.frame(popV = predict(PR_logit, data.frame(meanR = seq(1, 8, 0.1))),
                 meanR = seq(1, 8, 0.1))

PR<-ggplot(popVs, aes(x=meanR, y=popV, color= sp)) +
	geom_point()+ 
	scale_color_manual(values=spec.col, guide="none")+
		geom_line(data = dfPR, aes(y = popV, x = meanR), col="red") +	
	theme_classic() +
	ylab("Inter-individual diet variation (V)")

######
#Does annual precip predict PNW?
precip<-read.csv("diet_pops_30yrnorm_ppt_1980to2010_prism.csv")
precip2<-data.frame(Phylo_Code=precip$Phylo_Code, precip=precip$ppt_mm_1980to2010mean_prism)

popVs2<-merge(popVs, precip2, by="Phylo_Code")

#add column for obligate spec populations (stephensi and macrotis S, M)based on literature

sp<-substr(popVs2$Phylo_Code,1,1)
ifelse(sp[1]=="S"| sp[1]=="M", "Sp", "Ge")
popVs2$SG<-ifelse(sp=="S"| sp=="M", "Sp", "Ge")

PNW<-ggplot(popVs2, aes(x=precip, y=TNW, color= sp)) +
	geom_point()+ 
	scale_color_manual(values=spec.col, guide="none")+	
	theme_classic() +
	ylab("PNW")+ xlab("Precipitation (mm)")

lm_pcp<-lm(TNW~precip+SG, data=popVs2)#
summary(lm_pcp) #not sig

par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm_pcp, las = 1)  

MR<-ggplot(popVs2, aes(x=precip, y=meanR)) +
	#geom_smooth( method = lm)+	
	geom_point(data=popVs2, aes(x=precip, y=meanR, color=sp))+ 
	scale_color_manual(values=spec.col, guide="none")+
	theme_classic() +
	ylab("Mean Richness")+ xlab("Precipitation (mm)")

lm_mr<-lm(meanR~precip+SG, data=popVs2)#
summary(lm_mr) #not sig
plot(lm_mr) 

TR<-ggplot(popVs2, aes(x=precip, y=totR)) +
	#geom_smooth( method = lm)+	
	geom_smooth(method = "glm", , se = T, 
		method.args = list(family = "poisson"))+
	geom_point(data=popVs2, aes(x=precip, y=totR, color=sp))+ 
	scale_color_manual(values=spec.col, guide="none")+
	theme_classic() +
	ylab("Total Richness")+xlab("Precipitation (mm)")

#count data, poisson model
lm_tr<-glm(totR~precip+SG, family="poisson",data=popVs2)#
summary(lm_tr) #this one is signficant
plot(lm_tr) 


ggarrange(PNW + rremove("xlab")+rremove( "x.text"),
		 MR + rremove("xlab")+rremove( "x.text"), 
		TR, nrow = 3, ncol=1, 
		heights = c(2, 2, 2.3), align = c("v"), common.legend = TRUE)

#### Does county plant family count predict stuff?
famct<-read.csv("Diet_pops_factors_23mar24.csv")
famct2<-data.frame(Phylo_Code=famct$Phylo_Code, fams=famct$Plant_fams)

popVs3<-merge(popVs2, famct2, by="Phylo_Code")

#does precip and plant ct correlate
ggplot(popVs3, aes(x=precip, y=fams)) +
	geom_point()+ 
	theme_classic()

lm_rf<-glm(fams~precip, family ="poisson", data=popVs3)#
summary(lm_rf) #Not correlated
plot(lm_rf) 

PNWF<-ggplot(popVs3, aes(x=fams, y=TNW, color= sp)) +
	geom_point()+ 
	scale_color_manual(values=spec.col, guide="none")+	
	theme_classic() +
	ylab("PNW")+xlab("Plant Families")

lm_tr<-lm(TNW~fams+SG, data=popVs3)#
summary(lm_tr) #not sig
plot(lm_tr) 

MRF<-ggplot(popVs3, aes(x=fams, y=meanR, color= sp)) +
	geom_point()+ 
	scale_color_manual(values=spec.col, guide="none")+	
	theme_classic() +
	ylab("Mean Richness")+xlab("Plant Families")

lm_mr<-lm(meanR~fams+SG, data=popVs3)#
summary(lm_mr) #not sig
plot(lm_mr)

TRF<-ggplot(popVs3, aes(x=fams, y=totR, color= sp)) +
	geom_point()+ 
	scale_color_manual(values=spec.col, guide="none")+	
	theme_classic() +
	ylab("Total Richness")+xlab("Plant Families")

#count data, poisson model
lm_tr<-glm(totR~fams+SG, family="poisson",data=popVs3)#
summary(lm_tr) #not sig
plot(lm_mr)


pdf("dietfactors_22May24.pdf", width=4, height=7)
ggarrange(PNW + rremove("xlab")+rremove( "x.text"),
	PNWF+ rremove("xlab")+rremove( "x.text")+ rremove("ylab")+rremove( "y.text"),
	MR+ rremove("xlab")+rremove( "x.text") , 
	MRF + rremove("xlab")+rremove( "x.text")+ rremove("ylab")+rremove( "y.text"),
	TR, 
	TRF+ rremove("ylab")+rremove( "y.text"),
	nrow = 3, ncol=2, 
	heights = c(2, 2, 2.3),widths =c(2, 2.2),
	align=c("hv"), common.legend = TRUE, labels = "AUTO")
dev.off()