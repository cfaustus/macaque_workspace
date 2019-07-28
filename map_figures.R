#SYMPATRIC MACAQUE SPECIES
#getting shapefiles
library(sp)
library(rgdal)
# data <- readOGR("/Users/cfaust/Dropbox/MAMMTERR", "MAMMTERR") 
# names(data)
# unique <- unique(data@data$BINOMIAL)
# sympmac<-c("Macaca fascicularis", "Macaca nemestrina", "Macaca arctoides")
# sympmacID<-match(sympmac, data@data$BINOMIAL)
# 
# macs<-match(sympmac, unique)
# macs<-sort(macs)
# 
# for (i in 3291:3319) {
#   tmp <- data[data$BINOMIAL == unique[i], ] 
#   writeOGR(tmp, dsn="/Users/cfaust/Dropbox/Princeton/Projects/SimianMalaria/AlphaGlobin/SHP", 
#            unique[i], driver="ESRI Shapefile",overwrite_layer=TRUE)
# }


## Libraries
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(rgdal)

setwd("~/Dropbox/Projects/SimianMalaria/AlphaGlobinOX/macaque_workspace")

#Import shapefiles
fascicularis_range=readOGR("data/Macaca fascicularis.shp")
nemestrina_range=readOGR("data/Macaca nemestrina.shp")
arctoides_range=readOGR("data/Macaca arctoides.shp")

#SEAsia map
map('world', xlim=c(90,150),ylim=c(-12,30))
#map.axes() #adds axis of lat/long to map

plot(fascicularis_range, add=TRUE, col=alpha("red",0.5), border=FALSE)
plot(nemestrina_range, add=TRUE, col=alpha("red4",0.5), border=FALSE) #plots peni range
plot(arctoides_range, add=TRUE, col=alpha("black",0.5), border=FALSE) #plots peni range

#title("Macaque Ranges with AQ Phenotype")
par(font=1) #sets font back to italics
fasci=expression(italic(M.)~italic(fascicularis)) #make only some text italicized
nemi=expression(italic(M.)~italic(nemistrina)) #make only some text italicized
arc=expression(italic(M.)~italic(arctoides))
legend(125, 20,bty='n', c(fasci,nemi, arc),pch=19, cex=0.7, col=c(alpha("red",0.5),alpha("red4",0.5),alpha("black",0.5)))
