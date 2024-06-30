#diet diversity across species range
#using 5 species with atleast 4 sampled pops
#hows does richness and tnw vary between pops in a species? 

library("ggplot2")
library("dplyr")
library(phyloseq)
library("vegan") #species accumulation curves
library(reshape2) #melt function
library(multcomp) #for posthoc comparison

setwd("C:/location/here") 
 
psRF<-readRDS(file = "ps_rarF_Bigtrnl_29Nov23.rds") 

#set colors for the 5 species
lepcol<- "#4D04BF" #lepida purple
maccol<- "#FF0000" #macrotis red
flocol<-"#841618" #floridana dark red
albcol<-"#F7C908" #albigula yellow
brycol<- "#1B0121" #bryanti near black

sp5col<-c(lepcol,albcol, brycol, flocol, maccol)

#count of animals per pop
df<-data.frame(sample_data(psRF))
pops<-df %>% count(Phylo_Code)

#clean up data for pop analyses
#subset to pops with at least 3 animals, and have a Phylo_code
pop4<-pops[which(pops$n>2 & pops$Phylo_Code != ""),]
p4<-pop4$Phylo_Code
ps4p<-subset_samples(psRF, (Phylo_Code %in% p4))

data.frame(sample_data(ps4p))%>% count(Phylo_Code)

#restrict to species with multiple pops sampled
spec<-c("Neotoma lepida", "Neotoma bryanti", 
		"Neotoma macrotis", "Neotoma albigula", "Neotoma floridana")
psSp<-subset_samples(ps4p, (Spec %in% spec))

#remove samples from multiple sampling rounds at same site
#WW, remove all except one bryanti batch
dfb<-df[which(df$Month == "November" & df$Year== 2019 & df$Phylo_Code == "B17"),]
b17<-row.names(dfb)
dfWW<-df[which(df$Site == "White Water"),] #all the white water animals
WWnames<-row.names(dfWW)
wwremove<-setdiff(WWnames, b17)  #other white water L, LB, B to remove
psSp2<-subset_samples(psSp,!(sample_names(psSp) %in% wwremove))#remove extra wwbryanti

####for floridana, use 2022 Key largo batch
dfb1<-df[which(df$County == "Monroe" & df$Year== 2017),]
F17<-row.names(dfb1) #samples to remove
psSp3<-subset_samples(psSp2,!(sample_names(psSp2) %in% F17))

#check the list
df5<-data.frame(sample_data(psSp3)) 
cts5<-df5 %>% count(Phylo_Code) #looks reasonable

unique(df5$Site)

#average individuals per pop
mean(cts5$n)
sd(cts5$n)

#metrics for each pop 
pspops = merge_samples(psSp3, "Phylo_Code")#turns phylo codes to rownames
pspops <-prune_taxa(taxa_sums(pspops ) > 0, pspops )#
popsdf<-data.frame(estimate_richness(pspops , measures = c("Observed", "Shannon")))	
popsdf$sp<-substr(row.names(popsdf), 1, 1) #add back in species to dataframe as letter
popsdf$sp<-factor(popsdf$sp, levels = c("L", "A", "B", "F", "M"))
popsdf$phylo_code<-row.names(popsdf)

#add number of individuals per pop, use to set point size
popsdf2<-merge(popsdf, cts5, by.x="phylo_code", by.y= "Phylo_Code")

#how many of each species?
popsdf2 %>% count(sp) 

pops_long<-melt(popsdf2, id.vars=c("phylo_code",  "sp",  "n"))

#make plots showing variation in richness and TNW across pops in same species
pdf("spec5obs_tnw_d2.pdf", height=3, width=3)
ggplot(pops_long, aes(x=sp, y=value)) +
  	geom_jitter(aes(col=sp, size = n), shape=20,
			position=position_jitter(width = 0.1), alpha=0.6) + 
	facet_wrap(vars(variable), nrow = 2, scales = "free_y")+ scale_size_area()+
	scale_color_manual(values = sp5col)+
	theme_bw() + theme(legend.position = "none", strip.text.x = element_blank())
dev.off()

#Summary info
sumpops<-popsdf%>% group_by(sp) %>% 
  		summarise(mean_Obs=mean(Observed),
            sd_Obs= sd(Observed),
		max_Obs=min(Observed),
            min_Obs= max(Observed),
		mean_shan=mean(Shannon),
            sd_shan= sd(Shannon),
		max_shan=min(Shannon),
            min_shan= max(Shannon)	)

#merge by species to get species level richnes
psSpec<-merge_samples(psSp3, "Spec")#turns species to rownames

Specdf<-data.frame(estimate_richness(psSpec , measures = c("Observed", "Shannon")))	

#get average proportion of total richness captured in pop diet by species
Specdf$sp<-c("A", "B", "F", "L", "M")

df4<-merge(popsdf, Specdf, by = "sp")
df4$propObs<-df4$Observed.x/df4$Observed.y

props<-df4%>% group_by(sp) %>% 
  		summarise(mean_props=mean(propObs),
		sd_props= sd(propObs)	)

#pie charts for figure inset
pdf("pies.pdf", width=6, height=3)
par(mfrow = c(1, 5))
pie(c(0.223, 1-0.223), labels =NA)
pie(c(0.348, 1-0.348), labels =NA)
pie(c(0.327, 1-0.327), labels =NA)
pie(c(0.364, 1-0.364), labels =NA)
pie(c(0.429, 1-0.429), labels =NA)
dev.off()
###################
#pop species accumulation curve for figure
#sp5col<-c(lepcol,albcol, brycol, flocol, maccol)
spec5<-c( "Neotoma lepida", "Neotoma albigula", "Neotoma bryanti",  
		"Neotoma floridana", "Neotoma macrotis")

##plot the 5 species' diet accumulation curves
#plist = list() #not working?
pdf("diet_accum.pdf", height=3, width=11)
par(mfrow = c(1, 5)) 
for (i in 1:5){
	#print(spec5[i])
	psLep<-subset_samples(psSp3,Spec == spec5[i])
	psLep<-prune_taxa(taxa_sums(psLep) > 0, psLep)
	psLpop = merge_samples(psLep, "Phylo_Code")
	LepPopotu<-data.frame(otu_table(psLpop))

	lep_acc<-specaccum(LepPopotu, method="exact")#smoother curve than random
	plist[[i]] <- plot(lep_acc, ci.type="poly", col="black", lwd=2, ci.lty=0, 
			ci.col=alpha(sp5col[i],0.5),xlab="Populations",ylab="Diet items",  main= spec5[i],
			ylim=c(0,48))
	}
dev.off()

#for one species:
psLep<-subset_samples(psSp3,Spec == Species)
psLep<-prune_taxa(taxa_sums(psLep) > 0, psLep)
psLpop = merge_samples(psLep, "Phylo_Code")
LepPopotu<-data.frame(otu_table(psLpop))

#Species accumulation curves:
#classic method is "random" which finds the mean SAC
# and its standard deviation from random permutations of the data, 
# or subsampling without replacement 
#The "exact" method finds the expected SAC using sample-based rarefaction method
lep_acc<-specaccum(LepPopotu, method="exact")#smoother curve than random
plot(lep_acc, ci.type="poly", col="black", lwd=2, ci.lty=0, 
	ci.col="grey50", xlab="Populations", ylab="Diet Items",  main= Species)

##could fit models to this and extract information...
mod1 <- fitspecaccum(lep_acc, "lomolino")
coef(mod1)
fitted(mod1)
plot(mod1,add=TRUE, col="red")#adds this to the active plot
specslope(lep_acc, at=3.5)#for slope at a number of "plots"

####################
######## Stats
#Test extent to which diet richness and TNW for pop a
# function of species identity, latitude, sample size

#use popsdf2 data frame- add lat 
lats<-sample_data(pspops)[,12]

popsdf3<-merge(popsdf2, lats, by.x="phylo_code", by.y=0)

#shannon
ggplot(data=popsdf3, aes(x=Latitude, y=Shannon, color= sp, size = n)) + geom_point()
Sh_lm<-lm(Shannon~Latitude + sp + n ,data=popsdf3)
summary(Sh_lm)
#sequentially remove non-sig factors (-lat)
Sh_lm<-lm(Shannon~ sp + n ,data=popsdf3)
#-sample size
Sh_lm<-lm(Shannon~ sp,data=popsdf3)
##only species sig- post hoc tukey test
ph1<-glht(Sh_lm,linfct=mcp(sp = "Tukey")) 
summary(ph1)  #only macrotis and floridana sig differ

par(mfrow=c(2,2))
plot(Sh_lm) #inspect model, ok-ish

#richness
obs_glm<-glm(Observed~Latitude + sp + n ,family=poisson, data=popsdf3)
summary(obs_glm)
#sequentially remove non-sig factors (-lat)
obs_glm<-glm(Observed~ sp + n ,family=poisson, data=popsdf3)
#-sample size
obs_glm<-glm(Observed~ sp ,family=poisson, data=popsdf3)
##only species sig- post hoc tukey test
ph2<-glht(obs_glm,linfct=mcp(sp = "Tukey")) 
summary(ph2)  #only macrotis-floridana  and floridana-lepida sig differ

par(mfrow=c(2,2))
plot(obs_glm)

##############################
#TEST: at pop level is diet similarity higher within species 
	#and with increasing geographic proximity
#Could do a PERMANOVA for species effect, then Mantel test for geographic effect, 
# but that tests a slightly different question
pspops #41 pops in these 5 species
#Jaccard distance complementary to jaccard coefficient JD=1-JC. 
Jac<-distance(pspops, "jaccard", binary = TRUE) # vegdist jaccard
JM<-as.matrix(Jac)  #converts to a matrix
JM[lower.tri(JM, diag=TRUE)] <- NA  #keeps only upper triangle, no diagonal
dfJ <- melt(JM, varnames = c("samp1", "samp2"), value.name = "Jac", na.rm = T)

#Bray curtis distance (dissimilarity index)
Bray<-distance(pspops, "bray", binary = TRUE) # vegdist bray curtis
BM<-as.matrix(Bray)  #converts to a matrix
BM[lower.tri(BM, diag=TRUE)] <- NA  #keeps only upper triangle, no diagonal
dfB <- melt(BM, varnames = c("samp1", "samp2"), value.name = "Bray", na.rm = T)

identical(dfJ$samp2,dfB$samp2)#
identical(dfJ$samp1,dfB$samp1)#identical lists, allows lazy data merge

dfPJ<-cbind(dfJ, Bray=dfB$Bray)

######################################

#add species letters and if samesp of diffsp
dfPJ$sp1<-substr(dfPJ$samp1, 1, 1) 
dfPJ$sp2<-substr(dfPJ$samp2, 1, 1) 
dfPJ$same<-ifelse(dfPJ$sp2 == dfPJ$sp1, "sameSp","difSp")

#add lat and long, calculate difference between points
coords<-sample_data(pspops)[,11:12]
dfPJ2<-merge(dfPJ, coords, by.x="samp1", by.y=0)
dfPJ3<-merge(dfPJ2, coords, by.x="samp2", by.y=0)

library(geosphere)
dist<-0 #empty vector
for (i in 1:nrow(dfPJ3)) {
		P<- distm (c(dfPJ3$Longitude.x[i], dfPJ3$Latitude.x[i]), 
			c(dfPJ3$Longitude.y[i], dfPJ3$Latitude.y[i]), fun = distHaversine)
		dist[i]<-P
		}
dfPJ3$dist<-dist/1000  #in km
#dfPJ3$JacI<- 1-dfPJ3$Jac
##################
#function to transform data to remove zeros, from  Smithson & Verkuilen (2006)DOI: 10.1037/1082-989X.11.1.54
new.y <- function(a) {
   result <- (a *(dim(dfPJ3)[1]-1) + 0.5)/dim(dfPJ3)[1]
   print(result)
}
new.y(0)

library("betareg")

max(dfPJ3$Bray) #1- problem for beta regression- transform
min(dfPJ3$Bray)#0.2

max(dfPJ3$JacI) # 0.6666667
min(dfPJ3$JacI)#zero

Jacn<-new.y(dfPJ3$Jac)
Brayn<-new.y(dfPJ3$Bray)

J_logit <- betareg(Jacn~ same + dist, data = dfPJ3)
summary(J_logit)
plot(J_logit)

B_logit <- betareg(Brayn~ same + dist, data = dfPJ3)
summary(B_logit)
plot(J_logit)


sample_data(pspops)$sp<-substr(row.names(sample_data(pspops)), 1, 1) 
sample_data(pspops)$sp<-factor(sample_data(pspops)$sp, levels = c("L", "A", "B", "F", "M"))

##For FIGURE 4
orduB2 = ordinate(pspops, "PCoA", "bray")
evals <- orduB2$values$Eigenvalues
p3<-plot_ordination(pspops, orduB2, color="sp")

pdf("PCOA_spec5.pdf", height=4, width=4)
p3+theme_classic()+
	ggtitle("Diet by pop (BC)")+
	geom_point(size = 4, alpha=0.6)+
	scale_color_manual(values = sp5col)+
  	coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

orduJ2 = ordinate(pspops, "PCoA", "jaccard")
evals <- orduJ2$values$Eigenvalues
p4<-plot_ordination(pspops, orduJ2, color="sp")
p4+theme_classic()+
	ggtitle("Diet by pop (J)")+
  	coord_fixed(sqrt(evals[2] / evals[1]))

# Difference in diet composition between species? 
#Calculate bray curtis distance matrix
sp_bray <- phyloseq::distance(pspops, method = "bray")
sp_bray <- phyloseq::distance(pspops, method = "jaccard")

# make a data frame from the sample_data
sampledf2 <- data.frame(sample_data(pspops))

# Adonis test: Permutational Multivariate Analysis of Variance Using Distance Matrices
Br_ad<-adonis2(sp_bray ~ sample_data(pspops)$sp, data = sampledf2) #sig diff between species

install.packages("pairwiseAdonis")#not available for this R version?
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#default is 999 permutations
res<-pairwiseAdonis::pairwise.adonis(sp_bray, sampledf2[,"sp"], p.adjust.m = "BH",perm =10000)