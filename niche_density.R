#make density niche plots
#start with psWW and packages from WWseason script
#this is a clumsy way to do this, but it works
#For a survey in the wet, early dry, and late dry seasons

###########late dry plot example, 19 rats, 11 plants
ps1S<-subset_samples(psWW,Date == "2019-11-01")
ps1S<-prune_taxa(taxa_sums(ps1S) > 0, ps1S)#get rid of taxa not in this dataset
ps1S= transform_sample_counts(ps1S, function(x) x/sum(x))

otu<-data.frame(otu_table(ps1S))
tt<-data.frame(tax_table(ps1S))
colnames(otu)<-tt$Family

cm<-data.frame(family=tt$Family, val=colMeans(otu)) 
cm<-cm[order(-cm$val),]

#reorder otu columns by plant abundance in means
otu_idx <- match(row.names(cm), colnames(otu))
otur  <- otu[ , otu_idx]
oturt<-data.frame(t(otur))
oturt$pos<-c(0,-1,1,-2,2,-3,3,-4,4,-5,5)
colSums(oturt)#all sum to 5000 reads
xnames<-substr(row.names(oturt[order(oturt$pos),]), start = 1, stop = 2)

pLD<-ggplot()+
	geom_density(data=oturt, aes(x= pos, y = oturt[,1]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,2]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,3]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,4]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,5]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,6]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,7]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,8]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,9]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,10]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,11]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,12]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,13]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,14]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,15]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,16]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,17]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,18]), stat= "identity", fill="#914223", alpha=.2)+
	geom_density(data=oturt, aes(x= pos, y = oturt[,19]), stat= "identity", fill="#914223", alpha=.2)+
	scale_x_continuous( breaks = seq(-5,5,1), labels = xnames)+
	theme_classic()+
	theme(axis.title.x=element_blank(),	axis.title.y=element_blank())
	
#####################
#wet season plot #2020-03-09 14 rats, 14 taxa
ps1S<-subset_samples(psWW,Date == "2020-03-09")
ps1S<-prune_taxa(taxa_sums(ps1S) > 0, ps1S)#get rid of taxa not in this dataset
ps1S= transform_sample_counts(ps1S, function(x) x/sum(x))

otu<-data.frame(otu_table(ps1S))
tt<-data.frame(tax_table(ps1S))
colnames(otu)<-tt$Family

cm<-data.frame(family=tt$Family, val=colMeans(otu)) 
cm<-cm[order(-cm$val),]

#reorder otu columns by plant abundance in means
otu_idx <- match(row.names(cm), colnames(otu))
otur  <- otu[ , otu_idx]
oturt2<-data.frame(t(otur))
oturt2$pos<-c(0,-1,1,-2,2,-3,3,-4,4,-5,5,-6,6,-7)
colSums(oturt2)#all sum to 1
xnames<-substr(row.names(oturt2[order(oturt2$pos),]), start = 1, stop = 2)

pW<-ggplot()+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,1]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,2]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,3]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,4]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,5]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,6]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,7]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,8]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,9]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,10]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,11]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,12]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,13]), stat= "identity", fill="#076612", alpha=.2)+
	geom_density(data=oturt2, aes(x= pos, y = oturt2[,14]), stat= "identity", fill="#076612", alpha=.2)+
	scale_x_continuous( breaks = seq(-7,6,1), labels = xnames)+
	theme_classic()+
	theme(axis.title.x=element_blank(),	axis.title.y=element_blank())

#####################
#Early Dry season plot ##2020-07-14 16 rats, 10 tax
ps1S<-subset_samples(psWW,Date == "2020-07-14")
ps1S<-prune_taxa(taxa_sums(ps1S) > 0, ps1S)#get rid of taxa not in this dataset
ps1S= transform_sample_counts(ps1S, function(x) x/sum(x))

otu<-data.frame(otu_table(ps1S))
tt<-data.frame(tax_table(ps1S))
colnames(otu)<-tt$Family

cm<-data.frame(family=tt$Family, val=colMeans(otu)) 
cm<-cm[order(-cm$val),]

#reorder otu columns by plant abundance in means
otu_idx <- match(row.names(cm), colnames(otu))
otur  <- otu[ , otu_idx]
oturt3<-data.frame(t(otur))
oturt3$pos<-c(0,-1,1,-2,2,-3,3,-4,4,-5)
colSums(oturt3)#all sum to 1
xnames<-substr(row.names(oturt3[order(oturt3$pos),]), start = 1, stop = 2)

pED<-ggplot()+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,1]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,2]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,3]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,4]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,5]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,6]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,7]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,8]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,9]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,10]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,11]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,12]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,13]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,14]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,15]), stat= "identity", fill="#D1B305", alpha=.2)+
	geom_density(data=oturt3, aes(x= pos, y = oturt3[,16]), stat= "identity", fill="#D1B305", alpha=.2)+	
	scale_x_continuous( breaks = seq(-5,4,1), labels = xnames)+
	theme_classic()+
	theme(axis.title.x=element_blank(),	axis.title.y=element_blank())
#####################

library(cowplot)
pdf("trip_density_TNW.pdf", width=3, height=6)
plot_grid(pW, pED, pLD, ncol = 1, align = 'v')
dev.off()
