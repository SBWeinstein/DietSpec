#range maps for figure 4
#R version 4.2.2

setwd("C:/location/here")
samples<-read.csv("Diet_pops_20Mar24.csv")

#############
library("ggplot2")
library("rnaturalearth")
library(ggrepel)#label repel
library(geojsonio)
library("sf")



Nalb<-st_read("neotoma_5spp_individual_epsg4326/n_albigula_epsg4326.geojson")
Nbry<-st_read("neotoma_5spp_individual_epsg4326/n_bryanti_epsg4326.geojson")
Nflor<-st_read("neotoma_5spp_individual_epsg4326/n_floridiana_epsg4326.geojson")
Nlep<-st_read("neotoma_5spp_individual_epsg4326/n_lepida_epsg4326.geojson")
Nmac<-st_read("neotoma_5spp_individual_epsg4326/n_macrotis_epsg4326.geojson")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
MainStates <- map_data("state")
theme_set(theme_classic())

#make points bigger and use grepel to adjust spacing
options(ggrepel.max.overlaps = Inf)

##add species to dataframe, extract from Phylo_Code
samples$sp<-as.character(substr(samples$Phylo_Code, start = 1, stop = 1))

flor_samp<-samples[which(samples$sp=='F'),]

sub_samp<-samples[ samples$sp %in% c("A","B", "L", "F", "M"),]

spec.col<-c("#F7C908", #albigula yellow  A
		"#1B0121", #bryanti near black B
		"#841618", #floridana dark red F
		"#4D04BF", #lepida purple L
		"#FF0000") #macrotis red M

#big plot for floridana and ranges

pdf("flor_plot.pdf", width= 5, height=4, useDingbats=FALSE)
ggplot(data = world) +
 	geom_polygon( data=MainStates, aes(x=long, y=lat, group=group),
                color="grey", fill= "white")+
	geom_sf(fill= NA)+
	geom_sf(data = Nalb, colour = NA, fill = "#F7C908", alpha =.2)+
	geom_sf(data = Nbry, colour = NA, fill = "#1B0121", alpha =.2)+
	geom_sf(data = Nflor, colour = NA, fill = "#841618", alpha =.2)+
	geom_sf(data = Nlep, colour = NA, fill = "#4D04BF", alpha =.2)+
	geom_sf(data = Nmac, colour = NA, fill = "#FF0000", alpha =.2)+
	coord_sf(xlim = c(-128, -67), ylim = c(24, 50), expand = FALSE)+
	scale_color_manual(values=spec.col, guide="none")+
	#geom_point(data=sub_samp, 
	#		aes(x=Longitude, y=Latitude, color = sp ), size = 1.5, alpha = 0.8)+
	geom_point(data=flor_samp, 
			aes(x=Longitude, y=Latitude ),color= "#841618", size = 2, alpha = 0.8)+
	#geom_text_repel(data=sub_samp, 
	#		aes(x=Longitude, y=Latitude, label=Phylo_Code, color =sp)) +
	xlab("Longitude")+
	xlab("Longitude")+
	ylab("Latitude")+
	theme(legend.position="none")
dev.off()

#smaller region for A and L spp
sub_sampLA<-samples[ samples$sp %in% c("A","L"),]
spec.colLA<-c("#F7C908", #albigula yellow  A
		"#4D04BF") #lep

pdf("LA_plot.pdf", width= 4, height=4)
ggplot(data = world) +
 	geom_polygon( data=MainStates, aes(x=long, y=lat, group=group),
                color="grey", fill= "white")+
	geom_sf(fill= NA)+
	geom_sf(data = Nalb, colour = NA, fill = "#F7C908", alpha =.1)+
	geom_sf(data = Nbry, colour = NA, fill = "#1B0121", alpha =.1)+
	geom_sf(data = Nflor, colour = NA, fill = "#841618", alpha =.1)+
	geom_sf(data = Nlep, colour = NA, fill = "#4D04BF", alpha =.1)+
	geom_sf(data = Nmac, colour = NA, fill = "#FF0000", alpha =.1)+
	coord_sf(xlim = c(-121, -106), ylim = c(31, 41), expand = FALSE)+
	scale_color_manual(values=spec.colLA, guide="none")+
	geom_point(data=sub_sampLA, 
			aes(x=Longitude, y=Latitude, color = sp ), size = 4, 
			alpha = 0.6)+
	xlab("Longitude")+
	xlab("Longitude")+
	ylab("Latitude")+
	theme(legend.position="none")
dev.off()

#even smaller region for B and M spp
sub_sampBM<-samples[ samples$sp %in% c("B","M"),]
spec.colBM<-c("#1B0121", #B
		"#FF0000") #M
pdf("BM_plot.pdf", width= 4, height=4)
ggplot(data = world) +
 	geom_polygon( data=MainStates, aes(x=long, y=lat, group=group),
                color="grey", fill= "white")+
	geom_sf(fill= NA)+
	#geom_sf(data = Nalb, colour = NA, fill = "#F7C908", alpha =.1)+
	geom_sf(data = Nbry, colour = NA, fill = "#1B0121", alpha =.1)+
	#geom_sf(data = Nflor, colour = NA, fill = "#841618", alpha =.1)+
	#geom_sf(data = Nlep, colour = NA, fill = "#4D04BF", alpha =.1)+
	geom_sf(data = Nmac, colour = NA, fill = "#FF0000", alpha =.1)+
	coord_sf(xlim = c(-121, -115), ylim = c(32.5, 37), expand = FALSE)+
	scale_color_manual(values=spec.colBM, guide="none")+
	geom_point(data=sub_sampBM, 
			aes(x=Longitude, y=Latitude, color = sp ), size = 4, 
			alpha = 0.6)+
	#geom_text_repel(data=sub_sampBM, 
	#		aes(x=Longitude, y=Latitude, label=Phylo_Code, color =sp)) +
	xlab("Longitude")+
	xlab("Longitude")+
	ylab("Latitude")+
	theme(legend.position="none")
dev.off()
		
