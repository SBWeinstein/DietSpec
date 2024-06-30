#test which plants contribute to within individual diet similarity in plant presence
#starting with psWW3b phyloseq object
#and using packages and outputs from WWresamp script

#make sure dates recognized as dates
sample_data(psWW3b)$Date<- as.Date(sample_data(psWW3b)$Date,
  format = "%m/%d/%y")

##used in loop
#function for getting days between samples, ignoring year difference
d_diff <- function(x, y) {
    	x1 <- as.integer(format(x, "%j")) #get day of year
	y1 <- as.integer(format(y, "%j"))
   	v1 <-abs(x1-y1)  #absolute diff between dates
	v2 <-(365-max(x1,y1))+ min(x1,y1) #but end of year,eg jan is close to december
	min(v1,v2) #pick smaller of those two
	}

#function to transform to remove zeros or ones
new.y <- function(a) {
   result <- (a *(dim(t2)[1]-1) + 0.5)/dim(t2)[1]
   print(result)
}


#columns for categorical self/non-self and time since last sample
md<-data.frame(sample_data(psWW3b))
md1<-md[,c("Animal_ID","Date")] #get animal Id and date for each sample (and rowname)

###################
#is it a single plant?
#get list of all plants in animal diets (28 fams)
wwplants<-rownames(tax_table(psWW3b))

#empty data frame
leav<-data.frame(matrix(ncol=3,nrow=0, 
		dimnames=list(NULL, c("plant","estimate", "p"))))

for (i in 1:length(wwplants)){
	ditch<-wwplants[i]
	print(ditch)
	leav[i,1]<-ditch
	sub3 <- subset_taxa(psWW3b, rownames(tax_table(psWW3b)) != ditch )

	Jac<-distance(sub3, "jaccard", binary = TRUE) #Jac distance matrix
	JM<-as.matrix(Jac)  #converts to a matrix
	JM[lower.tri(JM, diag=TRUE)] <- NA  #keeps only upper triangle, no diagonal
	df2 <- melt(JM, varnames = c("samp1", "samp2"), value.name = "Jac", na.rm = T)

	#df2$JacI<-(1-df2$Jac) #convert to Jaccard Index, aka Jaccard similarity coefficient

	#add columns for animals IDs and sampling date for samp1 and samp2
	t1<-merge(df2, md1, by.x = "samp1", by.y = 0)
	t2<-merge(t1, md1, by.x = "samp2", by.y = 0)
	names(t2)[names(t2) == "Animal_ID.x"] <- "Animal_ID.S1"
	names(t2)[names(t2) == "Animal_ID.y"] <- "Animal_ID.S2"
	names(t2)[names(t2) == "Date.x"] <- "Date.S1"
	names(t2)[names(t2) == "Date.y"] <- "Date.S2"

	#add column, if s1_animal=s2_animal, self; else not self (ns)
	t2$self<-ifelse(t2$Animal_ID.S1 == t2$Animal_ID.S2, "self","ns")

	#make sure columns are recognized as dates
	t2$Date.S1<-as.Date(t2$Date.S1)
	t2$Date.S2<-as.Date(t2$Date.S2)

	#calculate days apart
		for(j in 1:dim(t2)[1]) {
     		  t2$times[j]<- d_diff(t2$Date.S1[j], t2$Date.S2[j])
			}

	t2$self <- as.factor(t2$self)

	#JacIn<-new.y(t2$JacI)
	Jacn<-new.y(t2$Jac)
	WW_logit2 <- betareg(Jacn~ self + times, data = t2) #beta regression
	#WW_logit2 <- lm(JacI~ self + times, data = t2) #linear regression


	#get estimate for selfself
	leav[i,2]<-summary(WW_logit2)$coefficients$mean[2,1] #beta
	#leav[i,2]<-summary(WW_logit2)$coefficients[2,1] 	#lm

	#get the p-value for selfself
	leav[i,3]<-summary(WW_logit2)$coefficients$mean[2,4] #beta
	#leav[i,3]<-summary(WW_logit2)$coefficients[2,4]		#lm
	}

leav_logit<-leav  #using beta regression

############3
#is it core plants (always in diet in dry, rainy, post rain)
#or ephemeral plants that are only there in rainy or PR season

#add extra column to meta_data for season
df<-sample_data(psWW3b)
df$Season<- ifelse(df$Month == "July" | df$Month == "June" | df$Month == "May" , "PR","")
df$Season<- ifelse(df$Month == "January" | df$Month == "February" | df$Month == "March" , "R",df$Season)
df$Season<- ifelse(df$Month == "October" | df$Month == "November" | df$Month == "September"| df$Month == "August" , "D",df$Season)

sample_data(psWW3b)<-df  #replace sample data
psWWs<-subset_samples(psWW3b,Season != "")#remove samples w/o Season
psWWs<-prune_taxa(taxa_sums(psWWs) > 0, psWWs)#get rid of taxa not in this dataset

merg1 = merge_samples(psWWs, "Season") #Just want family lists for seasons

OT<-data.frame(otu_table(merg1 ))
TOT<-data.frame(t(OT)) #transpose for euler package
TOT2<-TOT %>% mutate_all(as.logical) #uses dplyr

library(eulerr)
fit2 <- euler(TOT2)
plot(fit2, quantities = TRUE, fill= c("cornflowerblue", "green4", "orangered3"))

#get plant list present in all seasons (core)
core<- TOT2[which(TOT2$D==TRUE & TOT2$PR==TRUE & TOT2$R==TRUE),]
coreP<-rownames(core)

#ephemerals, not present in dry, but present in rainy
eph<- TOT2[which(TOT2$D==FALSE  & TOT2$R==TRUE),]
ephP<-rownames(eph)

#ephrpr<- TOT2[which(TOT2$D==FALSE  & TOT2$R==TRUE | TOT2$D==FALSE  & TOT2tunr$R==TRUE ),] #the same as eph

##test these groups: ephP and coreP
sub3 <- subset_taxa(psWW3b, rownames(tax_table(psWW3b)) %in% coreP ) #this keeps core taxa
subE <- subset_taxa(psWW3b, rownames(tax_table(psWW3b)) %in% ephP) #this keeps eph taxa

prunecore <- prune_taxa(!rownames(tax_table(psWW3b)) %in% coreP , psWW3b)#remove core taxa, not meaningful, as most diets go to zero
pruneeph <- prune_taxa(!rownames(tax_table(psWW3b)) %in% ephP , psWW3b)#remove eph taxa
#sub<-prunecore
sub<-pruneeph

	Jac<-distance(sub, "jaccard", binary = TRUE) #Jac distance matrix
	JM<-as.matrix(Jac)  #converts to a matrix
	JM[lower.tri(JM, diag=TRUE)] <- NA  #keeps only upper triangle, no diagonal
	df2 <- melt(JM, varnames = c("samp1", "samp2"), value.name = "Jac", na.rm = T)

	#add columns for animals IDs and sampling date for samp1 and samp2
	t1<-merge(df2, md1, by.x = "samp1", by.y = 0)
	t2<-merge(t1, md1, by.x = "samp2", by.y = 0)
	names(t2)[names(t2) == "Animal_ID.x"] <- "Animal_ID.S1"
	names(t2)[names(t2) == "Animal_ID.y"] <- "Animal_ID.S2"
	names(t2)[names(t2) == "Date.x"] <- "Date.S1"
	names(t2)[names(t2) == "Date.y"] <- "Date.S2"

	#add column, if s1_animal=s2_animal, self; else not self (ns)
	t2$self<-ifelse(t2$Animal_ID.S1 == t2$Animal_ID.S2, "self","ns")

	#make sure columns are recognized as dates
	t2$Date.S1<-as.Date(t2$Date.S1)
	t2$Date.S2<-as.Date(t2$Date.S2)

	#calculate days apart
		for(j in 1:dim(t2)[1]) {
     		  t2$times[j]<- d_diff(t2$Date.S1[j], t2$Date.S2[j])
			}

	t2$self <- as.factor(t2$self)

	JacIn<-new.y(t2$Jac)
	WW_logit2 <- betareg(JacIn~ self + times, data = t2) #beta regression, same qualitative pattern as lm
	summary(WW_logit2)	


#self-self similarity driven by core plants that are present in year round diets, 
#can remove all of the ephemeral palnts and sig effect remains, remove core plants and it disappears
