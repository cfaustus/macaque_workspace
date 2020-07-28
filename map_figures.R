#MAPS

## Libraries
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(rgdal)
library(sp)
library(mapplots)
library(maps)
library(plyr)
library(rgeos)
library(raster)

#########
#getting shapefiles
# data = readOGR("/Users/cfaust/Dropbox/MAMMTERR", "MAMMTERR") 
# names(data)
# unique = unique(data@data$BINOMIAL)
# sympmac = c("Macaca fascicularis", "Macaca nemestrina", "Macaca arctoides")
# sympmacID = match(sympmac, data@data$BINOMIAL)
# 
# macs<-match(sympmac, unique)
# macs<-sort(macs)
# 
# for (i in 3291:3319) {
#   tmp <- data[data$BINOMIAL == unique[i], ] 
#   writeOGR(tmp, dsn="/Users/cfaust/Dropbox/Princeton/Projects/SimianMalaria/AlphaGlobin/SHP", 
#            unique[i], driver="ESRI Shapefile",overwrite_layer=TRUE)
# }

### Phenotype Map
dat = read.csv("data/phenotype_ltmacaque_data.csv", header = T)
head(dat)
fascicularis_range=readOGR("data/macaques/Macaca fascicularis.shp")

# c('A','Q','P','AQ','AP','AQP',
#   'QP','AX','AQX','APX','QPX'), 
col_pal = c('#2166ac','#fddbc7','#b2182b','#67a9cf','#ef8a62','#d1e5f0')
#col_pal=c('#8c510a','#c7eae5','#01665e','#d8b365','#5ab4ac','#f6e8c3')
#col_pal=c('#291620','#2A2F3F','#184A50','#28644A','#607738','#A78138','#E98465')
# col_pal = c('white','gold3','#335B8E','gold',
#             '#B7DBDB','steelblue3', 'steelblue4','khaki3',
#             '#B5B867','#6CA18F','olivedrab4')
# map ####
map("world", xlim=c(93,130), ylim=c(-11.5,23), 
    col="gray90", fill=TRUE)

plot(fascicularis_range, add=TRUE, col=alpha("grey10",0.3), border=FALSE)
#points(dat$long, dat$lat, pch=19, col="red", cex=0.5,)  #plot sites to check locations
# phenotypes A	AQ	Q	P	AP	QP	AQP	AX	AQX	APX	QPX
segments(103.8,1.36,105.5,0.97, lwd=1)
for(i in 1:nrow(dat)){
  add.pie(z=c(dat$A[i],dat$Q[i], dat$AQ[i], dat$AX[i], dat$AQX[i],
              (dat$P[i]+dat$AP[i]+dat$AQP[i]+dat$QP[i]+dat$APX[i]+dat$QPX[i])), 
          radius=sqrt(dat$total[i])/8.2, 
          x=dat$long[i], y=dat$lat[i], lty = 0,
          col= col_pal,
          #c(alpha(col_pal, 0.6))
          labels="")
}
# legend ####
legend(121, 4.3, 
       title = 'phenotypes ',
       c('A','Q','AQ','AX','AQX',
         'other'), 
       fill = col_pal,
       bty = 'o', box.col = 'white',cex = 0.8,
       bg = alpha('grey80',0.75))
segments(97.5,-1,97.5,-10, lwd=2) #vertical bar
segments(97.5,-1,97.8,-1, lwd=2)
segments(97.5,-10,97.8,-10, lwd=2)
segments(127,17,127,8, lwd=2)
segments(127,17,126.7,17, lwd=2)
segments(127,8,126.7,8, lwd=2)
#export as 5x5 PDF

###########
### Virulent Map
cambodia = readOGR("data/KHM_adm/KHM_adm0.shp") #cambodia
singapore = readOGR("data/SGP_adm/SGP_adm0.shp")
thailand = readOGR("data/THA_adm/THA_adm0.shp")
indonesia = readOGR("data/IDN_adm/IDN_adm1.shp")
bali = indonesia[indonesia$NAME_1 == 'Bali',]
java_names = c('Jakarta Raya','Jawa Tengah', 'Jawa Barat', 'Jawa Timur',
               'Banten', 'Yogyakarta')
java = indonesia[indonesia$NAME_1 %in% java_names,]
malaysia = readOGR('data/MYS_adm/MYS_adm1.shp')
borneo_names = c('Sabah', 'Sarawak')
pen_malay = malaysia[!malaysia$NAME_1 %in% borneo_names,]

#intersection
vir_countries = bind(cambodia, singapore, bali, thailand, java, pen_malay)
fascicularis_range=readOGR("data/macaques/Macaca fascicularis.shp")
monkey_vir = intersect(fascicularis_range,vir_countries)

map("world", xlim=c(93,130), ylim=c(-11.5,23), 
    col="gray90", fill=TRUE)
plot(fascicularis_range, add=TRUE, col=alpha("grey10",0.3), border=FALSE)
plot(monkey_vir, add=TRUE, col=alpha("darkred",0.5), border=FALSE)
legend(110, 16.3, 
       title = 'ranges ',
       c('M. fascicularis', 'malaria surveys'), 
       fill = c('grey70', 'darkred'),
       bty = 'o', box.col = 'white',cex = 0.7)

malaria = read.csv("data/malaria_clean_maintext.csv", header =T)
head(malaria)
malaria$region
mal_tab = malaria[,4:7]
dimnames(mal_tab)[[1]] <- c("Cambodia",'Thailand','PeninsularMalaysia', 
                            'Singapore', "Java", "Bali")
mosaicplot(mal_tab, main ='',
           col= c('white','gray80', 'darkred', 'gray50'), 
           las = 4)

### Range overlap

setwd("~/Dropbox/Projects/SimianMalaria/AlphaGlobinOX/macaque_workspace")

#Import shapefiles
fascicularis_range=readOGR("data/macaques/Macaca fascicularis.shp")
nemestrina_range=readOGR("data/macaques/Macaca nemestrina.shp")
arctoides_range=readOGR("data/macaques/Macaca arctoides.shp")

#SEAsia map
map("world", xlim=c(85,145), ylim=c(-11.5,30), 
    col="gray90", fill=TRUE)
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
