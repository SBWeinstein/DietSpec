#How does sampling a subset of pop impact estimates of specialization?

# use large diet data sets from different populations at different 
#seasons to examine how smaller sample sizes influence estimated diet parameters,
#over a range of levels of dietary specialization 

#data sets:
#whitewater: 3 Seasons : Wet(49), EarlyDry (66), LateDry (64)
#Whitney wells Brynti from hills: 2 Seasons (no winter): Spring (43), Summer/Fall (59)
#Whitney wells Lepida from flats: 2 Seasons (no winter): Spring (46), Summer/Fall (67)
#key largo (16 animals, all "winter" December/Feb)
#####################################

library("ggplot2")
library("dplyr")
library(phyloseq)
library("vegan")#species accumulation curves
library(ggpubr)#arrange plots
library(RColorBrewer)

setwd("PATH HERE") 

#for Whitewater and Key largo data
psRF<-readRDS(file = "ps_rarF_Bigtrnl_29Nov23.rds") 
#for Whitney wells data
psRFMatocq<-readRDS(file = "ps_rarF_MatocqBL_26Oct24.rds") 

##################################
# Function to produce summary statistics (mean and bootstrap 95CI)
#on subsampling violin plots
data_summary <- function(x) {
   	m <- mean(x)
	#ymin <- m-sd(x)
   	#ymax <- m+sd(x)
   	CIL<- sort(x)[0.025*length(x)]
	CIH<-sort(x)[0.975*length(x)]
   return(c(y=m,ymin=CIL,ymax=CIH))
}

######################
#########################
###White water animals####
##########################
psWW<-subset_samples(psRF,Site == "White Water")#179 samples
psWW<-prune_taxa(taxa_sums(psWW) > 0, psWW)#32 Taxa
df<-sample_data(psWW)
df$Season<- ifelse(df$Month == "July" | df$Month == "June" | df$Month == "May" , "ED","")
df$Season<- ifelse(df$Month == "January" | df$Month == "February" | df$Month == "March" , "W",df$Season)
df$Season<- ifelse(df$Month == "October" | df$Month == "November" | df$Month == "September"| df$Month == "August", "LD",df$Season)
sample_data(psWW)<-df  #replace sample data
data.frame(data.frame(df) %>% count(Season)) #no april or dec

#Create seperate objects for each season
psWWw<-subset_samples(psWW,Season == "W")
psWWw<-prune_taxa(taxa_sums(psWWw) > 0, psWWw)#23 taxa and 49 samples

psWWed<-subset_samples(psWW,Season == "ED")
psWWed<-prune_taxa(taxa_sums(psWWed) > 0, psWWed)#22 taxa and 66 samples

psWWld<-subset_samples(psWW,Season == "LD")
psWWld<-prune_taxa(taxa_sums(psWWld) > 0, psWWld)#17 taxa and 64 samples

###########################
#Key Largo animals (F23)
###########################
psKL<-subset_samples(psRF,Site == "Key Largo")#16 samples
psKL<-prune_taxa(taxa_sums(psKL) > 0, psKL)#21 taxa
#add a Season column to psKL so format matches others
df<-sample_data(psKL)
df$Season<-rep("win", dim(df)[1])
sample_data(psKL)<-df  #replace sample data

#######################
#WhitneyWells animals
#######################
#get Bryanti from flats
psB48<-subset_samples(psRFMatocq,Genotype == "N. bryanti")#177 samples
psB48<-subset_samples(psB48,Habitat == "hill")#102 samples
psB48<-prune_taxa(taxa_sums(psB48) > 0, psB48)#33 Taxa

Zdf<-data.frame(sample_data(psB48))
Zdf %>% count(Zdf$Season) #spring 43, S/F= 49 

psL48<-subset_samples(psRFMatocq,Genotype == "N. lepida")#125 samples
psL48<-subset_samples(psL48,Habitat == "flats")#113 samples
psL48<-prune_taxa(taxa_sums(psL48) > 0, psL48)#29 Taxa

Zdf<-data.frame(sample_data(psL48))
Zdf %>% count(Zdf$Season) #spring 46, S/F= 67 

#Create seperate objects for each season, each pop
psL48sp<-subset_samples(psL48,Season == "Spring")
psL48sp<-prune_taxa(taxa_sums(psL48sp) > 0, psL48sp)#27 taxa and 46 samples

psL48sf<-subset_samples(psL48,Season == "Summer/Fall")
psL48sf<-prune_taxa(taxa_sums(psL48sf) > 0, psL48sf)# 21 taxa and 67 samples

psB48sp<-subset_samples(psB48,Season == "Spring")
psB48sp<-prune_taxa(taxa_sums(psB48sp) > 0, psB48sp)#28 taxa and 43 samples

psB48sf<-subset_samples(psB48,Season == "Summer/Fall")
psB48sf<-prune_taxa(taxa_sums(psB48sf) > 0, psB48sf)#25 taxa and 59 samples 

#########################
### for a population in a season, how many samples does it take to characterize diet?
#make species accumulation curves for each population in the relevant seasons
#Species accumulation curves:
#classic method is "random" which finds the mean SAC
# and its standard deviation from random permutations of the data, 
# or subsampling without replacement 
#The "exact" method finds the expected SAC using sample-based rarefaction method

#list of season pops:
pops<-c(psL48sf,psL48sp, psB48sf, psWWed,psWWld,psB48sp, psKL, psWWw)
names<-c("L48sf", "L48sp", "B48sf", "WWed", "WWld", "B48sp", "KL", "WWw")

n=length(pops) #8 pops
spec_plot_list <- vector("list", length = n)

for (i in 1:n){
	Popotu<-data.frame(otu_table(pops[[i]]))
	pop_acc<-specaccum(Popotu, method="exact")#
	#creating a dataframe for ggplot2
	data <- data.frame(Individuals=pop_acc$sites, Richness=pop_acc$richness, SD=pop_acc$sd)
	
	#add second curve for only "common" diet items, atleast 20% in one diet
	pops2<- transform_sample_counts(pops[[i]], function(x) replace(x, x<(.2*sum(x)),0) )  
	pops2<-prune_taxa(taxa_sums(pops2) > 0, pops2) 
	Popotu2<-data.frame(otu_table(pops2))
	pop_acc2<-specaccum(Popotu2, method="exact")#
	data2 <- data.frame(Individuals=pop_acc2$sites, Richness=pop_acc2$richness, SD=pop_acc2$sd)

	plot<-ggplot() +
  		geom_point(data=data, aes(x=Individuals, y=Richness)) +
  		geom_line(data=data, aes(x=Individuals, y=Richness)) +
  		geom_ribbon(data=data ,aes(x=Individuals, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2)+

		geom_point(data=data2, aes(x=Individuals, y=Richness), color="red") +
  		geom_line(data=data2, aes(x=Individuals, y=Richness), color="red") +
  		geom_ribbon(data=data2 ,aes(x=Individuals, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill="red")+

		labs(x="Individuals")+
		annotate("label", x = 5, y = max(data$Richness)-1, label = names[i], size=3) +
		theme_classic()
	spec_plot_list[[i]]<-plot
	}

##see all the spec accum plots:
ggarrange(spec_plot_list[[1]],
	spec_plot_list[[2]],
	spec_plot_list[[3]],
	spec_plot_list[[4]],
	spec_plot_list[[5]],
	spec_plot_list[[6]],
	spec_plot_list[[7]],
	spec_plot_list[[8]],
		 nrow = 4, ncol=2, 
		widths = c(2, 2), align = c("hv"))

###
#how well do we capture the dominant diet component with lower sample numbers?
#expect this to depend on degree of specialization
#for a pop/season full dataset: what is the dominant diet item? average pop diet?
#first make a barplot for each pop/season using merged average diet of all individuals

n=length(pops) #8 pops
barplot_list <- vector("list", length = n)

for (i in 1:n){
	psPop<-pops[[i]]
	merg = merge_samples(psPop, "Season")#merge into 1 sample
	#identify most abundant diet component per pop
	psFr = transform_sample_counts(merg, function(x) x/sum(x)) #convert counts to relative abundance
	#plot_bar(psFr , fill="Family")

	##plot with rare taxa as "other"
	y4 <- psmelt(psFr) # create dataframe from phyloseq object
	y4$Family <- as.character(y4$Family) #convert to character
	y4$Family[y4$Abundance < 0.01] <- "< 1% of total" #rename Families with < 1% abundance

	#set order so that most abundant is always at the bottom
	y4 <- y4 %>% arrange(Abundance)
	y4$Family <- factor(y4$Family, levels = unique(y4$Family))

	#set color palette to accommodate the number of genera
	colourCount = length(unique(y4$Family))
	getPalette = colorRampPalette(brewer.pal(9, "Set1"))

	p <- ggplot(data=y4, aes(x=Sample, y=Abundance, fill=Family))
	plot<-p + geom_bar(aes(), stat="identity", position="stack") + 
		scale_fill_manual(values=rev(getPalette(colourCount))) +
		theme_classic()+
		theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
				legend.text = element_text(size=12))
	#plot
	barplot_list[[i]]<-plot
	}

##see all the bar plots:
ggarrange(barplot_list[[1]],
	barplot_list[[2]],
	barplot_list[[3]],
	barplot_list[[4]],
	barplot_list[[5]],
	barplot_list[[6]],
	barplot_list[[7]],
	barplot_list[[8]],
		 nrow = 4, ncol=2, 
		widths = c(2, 2), align = c("hv"))

###############################################
###How do subsampling values compare to the true dominant plant in the full diet?
#How well does a random subset of 20/10/5/3 animals from a pop and season capture this? 
#can we put 95% confidence intervals on estimate of the amount of the true dom plant in the pop diet?

#function for CI in table
data_summary2 <- function(x) {
   	m <- mean(x)
	sd <- sd(x)
   	CIL<- sort(x)[0.025*length(x)]
	CIH<-sort(x)[0.975*length(x)]
   return(c(prop=m,sd=sd, CIL=CIL,CIH=CIH))
	}

pops<-c(psL48sf,psL48sp, psB48sf, psWWed,psWWld,psB48sp, psKL, psWWw)
names<-c("L48sf", "L48sp", "B48sf", "WWed", "WWld", "B48sp", "KL", "WWw")

#loop thru the 8 pops
MRTR_list<-	list()    #MR and TR values from sim
domdf_list<-list()  #counts of dom plant from each iteration
domplantCI_list<-list()
plotforpop_list<-list()

for (j in 1:8){

testpop<-pops[[j]]
testpopname<-names[j]

	#For each pop, get "true" TR and MR
	TR<-dim(otu_table(testpop))[2]
	obs<-data.frame(suppressWarnings(estimate_richness(testpop, measures="Observed")))
	MR<-mean(obs$Observed)
	MRTR_true<-c(dim(otu_table(testpop))[1], MR, NA, TR, NA) #values match up with subsampling values

	#for each pop: get the "true" dominant plant from max amount of samples 
	mergS = merge_samples(testpop, "Season")#all into one
	psFr = transform_sample_counts(mergS, function(x) x/sum(x)) #convert counts to relative abundance
	OT<-data.frame(otu_table(psFr))
	plant<-colnames(OT)[max.col(OT,ties.method="first")] #get the most abundant plant for the pop
	maxP<-max.col(OT, "first") #get max plant column index value (pick first if ties)
	value <- OT[cbind(1:nrow(OT), maxP)]  #get value from each of those columns
	res <- c( plant, value, "NA", dim(otu_table(testpop))[1], 1)   #vector: dominant plant seq, proportion, other stuff that lines up subsample data #

	TT1<-data.frame(tax_table(psFr))[,4:5] #for identifying plant Families later

	#Set up for subsampling at each level within each pop

	inds<-row.names(sample_data(testpop))
	list20 <- list()
	list10 <- list()
	list5 <- list()
	list3 <- list()

	MRTR20 <- list()
	MRTR10 <- list()
	MRTR5 <- list()
	MRTR3 <- list()

	for (i in 1:1000){
		sub20 <- ifelse(dim(otu_table(testpop))[1]<20, 15,20)
		indsSub<-sample(inds, sub20)#do this part 1000x

		psSub<-subset_samples(testpop, row.names(sample_data(testpop)) %in% indsSub)
		psSub<-prune_taxa(taxa_sums(psSub) > 0, psSub) 
		
		#get TR and MR for this subsample
		TR<-dim(otu_table(psSub))[2]
		obs<-data.frame(suppressWarnings(estimate_richness(psSub, measures="Observed")))
		MR<-mean(obs$Observed)
		MRTR<-c(MR, TR, 20)
		MRTR20 [[i]]  <- MRTR

		mergS = merge_samples(psSub, "Season")#all into one
		psFr = transform_sample_counts(mergS, function(x) x/sum(x)) #convert counts to relative abundance
		OT<-data.frame(t(otu_table(psFr))) #convert to dataframe and transform so that each plant is a row
		OT <- cbind(rownames(OT), data.frame(OT, row.names=NULL)) #make sequences a column
		colnames(OT) <- c("seq", "prop")
	
		OT$it<-rep(i, dim(OT)[1]) #add column to show iteration
		OT$n<-rep(20, dim(OT)[1]) #add column to show sample size
		OT <- OT[order(OT$prop, decreasing = TRUE),]
		OT$rank<-seq(from= 1, to =dim(OT)[1], by=1)
		
		list20[[i]] <- OT #put dataframe into list20
		}
	df20 <- do.call("rbind",list20) #combine all vectors into a dataframe
	MRTR20_df<- do.call("rbind",MRTR20)#combine all of the MR and TR values

	for (i in 1:1000){
		indsSub<-sample(inds, 10)#

		psSub<-subset_samples(testpop, row.names(sample_data(testpop)) %in% indsSub)
		psSub<-prune_taxa(taxa_sums(psSub) > 0, psSub) 
		
		TR<-dim(otu_table(psSub))[2]
		obs<-data.frame(suppressWarnings(estimate_richness(psSub, measures="Observed")))
		MR<-mean(obs$Observed)
		MRTR<-c(MR, TR, 10)
		MRTR10 [[i]]  <- MRTR

		mergS = merge_samples(psSub, "Season")#all the same genotype, 1 sample
		psFr = transform_sample_counts(mergS, function(x) x/sum(x)) #convert counts to relative abundance
		OT<-data.frame(t(otu_table(psFr))) #convert to dataframe and transform so that each plant is a row
		OT <- cbind(rownames(OT), data.frame(OT, row.names=NULL))
		colnames(OT) <- c("seq", "prop")
	
		OT$it<-rep(i, dim(OT)[1]) #add column to show iteration
		OT$n<-rep(10, dim(OT)[1]) #add column to show sample size
		OT <- OT[order(OT$prop, decreasing = TRUE),]
		OT$rank<-seq(from= 1, to =dim(OT)[1], by=1)
		
		list10[[i]] <- OT #put dataframe into list10
		}
	df10 <- do.call("rbind",list10) #combine all vectors into a matrix
	MRTR10_df<- do.call("rbind",MRTR10)#combine all of the MR and TR values

	for (i in 1:1000){
		indsSub<-sample(inds, 5)#

		psSub<-subset_samples(testpop, row.names(sample_data(testpop)) %in% indsSub)
		psSub<-prune_taxa(taxa_sums(psSub) > 0, psSub) 
		
		TR<-dim(otu_table(psSub))[2]
		obs<-data.frame(suppressWarnings(estimate_richness(psSub, measures="Observed")))
		MR<-mean(obs$Observed)
		MRTR<-c(MR, TR, 5)
		MRTR5 [[i]]  <- MRTR

		mergS = merge_samples(psSub, "Season")#all the same genotype, 1 sample
		psFr = transform_sample_counts(mergS, function(x) x/sum(x)) #convert counts to relative abundance
		OT<-data.frame(t(otu_table(psFr))) #convert to dataframe and transform so that each plant is a row
		OT <- cbind(rownames(OT), data.frame(OT, row.names=NULL)) #make sequences a column
		colnames(OT) <- c("seq", "prop")
	
		OT$it<-rep(i, dim(OT)[1]) #add column to show iteration
		OT$n<-rep(5, dim(OT)[1]) #add column to show sample size
		OT <- OT[order(OT$prop, decreasing = TRUE),]
		OT$rank<-seq(from= 1, to =dim(OT)[1], by=1)
		
		list5[[i]] <- OT #put dataframe into list5
		}
	df5 <- do.call("rbind",list5) #combine all vectors into a matrix
	MRTR5_df<- do.call("rbind",MRTR5)#combine all of the MR and TR values

	for (i in 1:1000){
		indsSub<-sample(inds, 3)#

		psSub<-subset_samples(testpop, row.names(sample_data(testpop)) %in% indsSub)
		psSub<-prune_taxa(taxa_sums(psSub) > 0, psSub) 
		
		TR<-dim(otu_table(psSub))[2]
		obs<-data.frame(suppressWarnings(estimate_richness(psSub, measures="Observed")))
		MR<-mean(obs$Observed)
		MRTR<-c(MR, TR, 3)
		MRTR3 [[i]]  <- MRTR

		mergS = merge_samples(psSub, "Season")#all the same genotype, 1 sample
		psFr = transform_sample_counts(mergS, function(x) x/sum(x)) #convert counts to relative abundance
		OT<-data.frame(t(otu_table(psFr))) #convert to dataframe and transform so that each plant is a row
		OT <- cbind(rownames(OT), data.frame(OT, row.names=NULL)) #make sequences a column
		colnames(OT) <- c("seq", "prop")
	
		OT$it<-rep(i, dim(OT)[1]) #add column to show iteration
		OT$n<-rep(3, dim(OT)[1]) #add column to show sample size
		OT <- OT[order(OT$prop, decreasing = TRUE),]
		OT$rank<-seq(from= 1, to =dim(OT)[1], by=1)
		
		list3[[i]] <- OT #put dataframe into list3
		}
	df3 <- do.call("rbind",list3) #combine all vectors into a matrix
	MRTR3_df<- do.call("rbind",MRTR3)#combine all of the MR and TR values
	
	MRTR_all<-data.frame(rbind( MRTR20_df, MRTR10_df,MRTR5_df, MRTR3_df))
	names(MRTR_all)<-c("MR", "TR", "n")
	##MRTR_all$n <- factor(MRTR_all$n, levels=c("3", "5", "10", "20"))
	MRTR_sum<-MRTR_all%>% group_by(n) %>% summarize( mean_MR= mean(MR), sd_MR=sd(MR), mean_TR=mean(TR), sd_TR=sd(TR))
	MRTR_sum<-data.frame(MRTR_sum)
	MRTR_sum<-rbind(MRTR_sum, MRTR_true)####broken
	MRTR_sum$name<-rep(names[j], length(MRTR_sum$n))##forgot to replace 1 with j!!in big run

###save the summary of MR and TR for this pop
MRTR_list[[j]]<-	MRTR_sum

	dfall<-data.frame(rbind(res, df20, df10,df5, df3))
	dfall$n <- factor(dfall$n, levels=c(res[4], "20", "10", "5", "3"))
	dfallT<-merge(dfall, TT1, by.x="seq", by.y = 0) 

	##get the dominant plant from each iteration
	df1<-dfallT[which(dfallT$rank == 1 & dfallT$it != "NA" ),]
	df1$prop<-as.numeric(df1$prop)#make sure the proportions are numeric
	domplants<-df1 %>% group_by(n,Family) %>% summarize( Count = n() )
	domdf<-data.frame(domplants)
	domdf<-domdf[order(domdf$n, -domdf$Count),]
	domdf$name<-rep(names[j], length(domdf$n))

###get the summary on dominant plant counts for this pop
domdf_list[[j]]<-domdf
	
	#make barplot of distribution of IDed dom plants
	bplot<-ggplot(domdf, aes(x=n, y=Count, fill=Family)) + 
  		geom_bar(stat="identity")+
		theme_classic()+
		theme(legend.position="top")

	###########################
	#Plot estimated abundance of true dom plant in diet
	#for plotting, subset to seqs that match true dominant plant
	df1ps<-dfallT[which(dfallT$seq == res[1]& dfallT$it != "NA"),]
	df1ps$prop<-as.numeric(df1ps$prop)#make sure the proportions are numeric

	#summary data on dominant plant proportion calculations
	ds3<-c(data_summary2(df1ps[which(df1ps$n == 3 ),]$prop), 3, names[j])  #needs minimum # iterations to generate CIL
	ds5<-c(data_summary2(df1ps[which(df1ps$n == 5 ),]$prop), 5, names[j])
	ds10<-c(data_summary2(df1ps[which(df1ps$n == 10 ),]$prop), 10, names[j])
	ds20<-c(data_summary2(df1ps[which(df1ps$n == 20 ),]$prop), 20, names[j])

	dstrue<-c(res[2], NA,NA,NA,dim(otu_table(testpop))[1], names[j] ) #add the "true value" of the dom plant
	domplantCI<-data.frame(rbind( ds3,ds5,ds10,ds20, dstrue))
	domplantCI$V5<-as.numeric(domplantCI$V5)
	domplantCI<-domplantCI[order( -domplantCI$V5),] #reorder to match others ds3 on the bottom

##get dom plant average prop, sd, upper lower CI
domplantCI_list[[j]]<-domplantCI

	vplot<-ggplot(df1ps, aes(x=n, y=prop)) + 
 		geom_hline(yintercept = as.numeric(res[2]), color="red", linetype="dashed") +
		geom_violin(trim=FALSE,color="transparent", fill="grey40", alpha=0.8)+
		stat_summary(fun.data=data_summary, color="red")+
		theme_classic()

	resam_plots<-ggarrange(bplot + rremove("xlab") + rremove("x.text"),
			vplot + rremove("xlab") ,
		 	nrow = 2, ncol=1, 
			heights = c(1.5, 2), align = c("v"))
	
plotforpop_list[[j]]<-resam_plots
}


#make figure that has all the combo plots for the 8 pops

comboplot_list<-list()

for (k in 1:8){
	testpop<-pops[[k]]
	testpopname<-names[k]
	
	#combo plot for 4 pop figures
	combopop<-ggarrange(barplot_list[[k]],
			spec_plot_list[[k]], 
			plotforpop_list[[k]],
		 	 ncol=3)
comboplot_list[[k]]<-combopop
}


all_plots1.4<-ggarrange(comboplot_list[[1]],
			comboplot_list[[2]] ,
			comboplot_list[[3]] ,
			comboplot_list[[4]] ,
			comboplot_list[[5]] ,
			comboplot_list[[6]] ,
			comboplot_list[[7]] ,
			comboplot_list[[8]] ,
		 	nrow = 4, ncol=1, align = c("hv"))

#######
#combine summary data from all 8 pops for supplementary table:
MRTR_all8<- do.call("rbind",MRTR_list) #combine all dataframes into a big one
write.csv(MRTR_all8, "MRTR_all8_7Nov24.csv")

domdf_all8<- do.call("rbind",domdf_list) #combine all dataframes into a big one
write.csv(domdf_all8, "domdf_all8_7Nov24.csv")

domplantCI_all8<- do.call("rbind",domplantCI_list) #combine all dataframes into a big one
write.csv(domplantCI_all8, "domplantCI_all8_7Nov24.csv")

