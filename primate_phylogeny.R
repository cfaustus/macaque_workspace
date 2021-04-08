#creating a phylogeny of primates with alpha globin stats
library(ape)
library(phytools)

# importing mammalian supertree from Bininda-Emonds et al. 2007
# https://doi.org/10.1038/nature05634
mammalsbest <- read.nexus("data/Bininda-emonds_2007_mammals_best.nex")

macaques <- c("Macaca_arctoides","Macaca_assamensis","Macaca_thibetana",
              "Macaca_radiata","Macaca_sinica","Macaca_cyclopis",
              "Macaca_mulatta","Macaca_fuscata","Macaca_fascicularis",
              "Macaca_maura","Macaca_nigra","Macaca_ochreata",
              "Macaca_tonkeana","Macaca_nemestrina","Macaca_silenus",
              "Macaca_sylvanus")
macaques.pruned <- drop.tip(mammalsbest,mammalsbest$tip.label[-match(macaques, mammalsbest$tip.label)])

plot(macaques.pruned)
nodelabels()
edgelabels()
tiplabels()
edgecol <- rep('black', length(macaques.pruned$edge.length))
edgecol[c(17,26,9,3, 10, 6)] <- 'red3'
# change lines for macaques that we do not have hemoglobin types
edgelty <- rep(1, length(macaques.pruned$edge.length) )
edgelty[c(28,27,7)] <- 3
# change line with to reflect the log (base10) of sample size
edgewidth<- rep(1.5,length(macaques.pruned$edge.length))
edgewidth[c(26,22,17,16,15,14,24,23,20,9,3)] <- 
  c(log10(255),#nemestrina
    log10(75),#nigra
    log10(2371),#fascicularis
    log10(1271), #fuscata
    log10(533), #mulatta
    log10(152), #cyclopis
    log10(49), #tonkeana
    log10(18), #ochreata
    log10(66), #maura
    log10(43), #radiata
    log10(247)) #arctoides
tipcol<- rep('black',length(macaques.pruned$tip.label))
#tipcol[c(2,10,15)]<-'red3'
tipcol[c(2,3,10,12,13,15,16)]<-'gray62'
plot(macaques.pruned,type="phylogram", 
     edge.color=edgecol, edge.width=edgewidth,edge.lty=edgelty,
     label.offset=0, font=3, cex=0.7,
     tip.color=tipcol, adj=1)

