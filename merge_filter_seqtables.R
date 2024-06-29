#Merge Sequence Tables
#R version (4.2.2)

#code based on:
#Dada2: https://benjjneb.github.io/dada2/tutorial.html
#For multiple run modifications: https://benjjneb.github.io/dada2/bigdata.html

###Load packages
.cran_packages <- c("ggplot2", "gridExtra", "knitr")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn", "BiocStyle")
# Load packages into session
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

setwd("WORKING DIRECTORY HERE") 

#load in all sequence tables

core<- readRDS("seqtab_trnL_core_17July23.rds")
Dylan2<- readRDS("seqtab_trnL_Dylan2_18July23.rds")
Dylan3<- readRDS("seqtab_trnL_Dylan3_18July23.rds")
Dylan4<- readRDS("seqtab_trnL_Dylan4_17July23.rds")
PhyloRerun<- readRDS("seqtab_trnL_phyloRerun_18July23.rds")
Plate1<- readRDS("seqtab_trnL_plate1_18July23.rds")
Plate2<- readRDS("seqtab_trnL_plate2_17July23.rds")
Plate4<- readRDS("seqtab_trnL_plate4_7April23.rds")
Plate435c<- readRDS("seqtab_trnL_plate435c_17July23.rds")
Plate3<-readRDS("seqtab_trnL_plate3_1Sep23.rds")

#merge all sequence tables together 
tables1 <- list(core, Dylan2, Dylan3, Dylan4, PhyloRerun, Plate1,
		Plate2, Plate4, Plate435c, Plate3)
st.all <- mergeSequenceTables(tables=tables1)
table(nchar(getSequences(st.all)))
dim(st.all)

#####################################################
#################################################
#remove chimeras by comparing each inferred sequence to the others in the table,
# and removing those that can be reproduced by stitching together two more abundant 
#sequences. it is typical that chimeras comprise a substantial fraction of inferred 
#sequence variants, but only a small fraction of all reads. 

seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=FALSE)
dim(seqtab.nochim)

#save merged
#This is the OTU table: OTUs as x, sample names as y
saveRDS(seqtab.nochim, "seqtab_trnL_merged_1Sep23.rds") 
write.csv(seqtab.nochim, "seqtab_trnL_merged_1Sep23.csv")
head(seqtab.nochim)
####

###Taxonomy assignment using STAND pipeline
#https://github.com/robertgreenhalgh/stand
##these are command line commands
##conda activate stand
#./run the script stuff


###make phyloseq object
library("phyloseq")

seqtab.nochim<- readRDS("seqtab_trnL_merged_1Sep23.rds")

taxatab<-read.csv("trnL.Consensus_50.2022.09.01.csv")#Note that this date was actually 2023
taxa<-as.matrix(data.frame(taxatab,row.names=1))
head(taxatab)

#prep sample data. ***make sure row names match sequence table row names
#read in csv of sample data, already in same order as sequence table data
samdf <- read.csv("meta_data_2Sep23.csv", header=TRUE)
samdf2<- samdf[,-1] #remove first column that has sample order
samdf3 <- data.frame(samdf2, row.names = 1)

head(rownames(samdf3))
head(rownames(seqtab.nochim))

#make phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf3),tax_table(taxa))
saveRDS(ps, "ps_raw_Bigtrnl_2Sep23.rds") #includes trnl plate 3 from other project

sum(sample_sums(ps)) #404207194

#remove samples that belong to other projects
ps1<-subset_samples(ps, exclude != "otherproj")

#remove taxa not assigned to Kingdom Viridiplantae
psV<-subset_taxa(ps1, Kingdom== "Viridiplantae")

#remove obvious bait contaminants(apple, peanut, soy, oats), Keep NAs
psb <- subset_taxa(psV, 
	Genus!="Malus" & 
	Genus!="Avena" & 
	Genus!="Arachis" & 
	 Genus!="Glycine"| 
	is.na(Genus))

##################################
#look at blanks
ps_blanks<-subset_samples(psb,Type != "feces")
ps_blanks<-prune_taxa(taxa_sums(ps_blanks) > 0, ps_blanks)

p<-plot_bar(ps_blanks, x="Animal_ID", fill= "Family")
p+facet_grid( cols=vars(plate), scales="free"

#removing components less than 1% of diet will remove any contamination.
#####################################################

##for each sample: set taxa < 1% to zero
ps2<- transform_sample_counts(psb, function(x) replace(x, x<(.01*sum(x)),0) )  
ps3<-prune_taxa(taxa_sums(ps2) > 0, ps2) #remove anything no longer present in database, reduces to 181 taxa

sumdf<-data.frame(sample_sums(ps1),sample_sums(psV),sample_sums(psb), sample_sums(ps3))

psf<-subset_samples(ps3,Type == "feces")
mean(sample_sums(psf)) #668874
sd(sample_sums(psf)) #889913

#################### How to subsample?
#woodrat diets are not *that* diverse
plot_richness(psf, x= "Phylo_Code", measures= "Observed")

#sequencing depths vary substantially; control for different sampling effort
#examine rarefaction curves
library("vegan")
#OTUs should be columns for vegan
tab <- otu_table(psf)
class(tab) <- "matrix" # force it into a matrix
#tab <- t(tab) # transpose observations to rows
rare <- rarecurve(tab, step=200, lwd=1, ylab="OTU",  label=F, xlim=c(0,100000))
abline( v=5000, col="blue") #rarefy to 5,000; captures community even in most diverse diets
abline( v=10000, col="purple") 
abline( v=12000, col="red")
abline( v=15000, col="orange")
################################################

psR<-rarefy_even_depth(psf, 5000) 

sums<-data.frame(sample_sums(psf))
write.csv(sums, "sample_reads_2Sep23.csv")

saveRDS(psR, "ps_rar_Bigtrnl_2Sep23.rds") #376 taxa and 522 samples

#Calculated percent of mOTUs IDed to Class, Family, Genus
tt<-data.frame(tax_table(psR))
library("dplyr")
tt %>% count(Family)#all ided to fam level
tt %>% count(Genus) #41 NAs, 190 genera including NAs, 
#(376-47)/376 = 87.5% of OTUs ided to genus
tt %>% count(Species)#some species have multiple OTUs
sum(sample_sums(psR))#2610000 reads
gNA<-subset_taxa(psR, is.na(Genus))
sum(sample_sums(gNA))#365193reads, (2610000 -365193)/2610000 = 86% of reads ided to genus

#merge to family level for subsequent analyes
psRF<-tax_glom(psR, taxrank="Family", NArm=FALSE)

saveRDS(psRF, "ps_rarF_Bigtrnl_2Sep23.rds") #94 taxa (families) and 522 samples

#move to other scripts for subsequent analyses
##################
###################################
#sanity check: we are seeing cactus feeding where we expect to see it. YES!
dev.new()
As<-subset_samples(psR, Species=="Neotoma albigula")
AS2<-prune_taxa(taxa_sums(As) > 0, As)
plot_bar(AS2, x="Animal_ID", fill= "Family")

ASC<-subset_taxa(AS2, Family == "Cactaceae")
plot_bar(ASC, x="Animal_ID", fill= "Genus")
