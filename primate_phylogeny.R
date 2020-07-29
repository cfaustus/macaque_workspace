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
edgecol[c(17,26,9,3, 10,7)] <- 'red3'
# change lines for macaques that we do not have hemoglobin types
edgelty <- rep(1, length(macaques.pruned$edge.length) )
edgelty[c(28,27,10,7,6)] <- 3
# change line with to reflect the log (base10) of th
edgewidth<- rep(1.5,length(macaques.pruned$edge.length))
edgewidth[c(26,22,17,16,15,14,24,23,20,9,3)] <- 
  c(log10(510),log10(53),log10(2371),log10(2539),
    log10(633),log10(301),log10(49),log10(18),log10(66),log10(38),log10(256))
tipcol<- rep('black',length(macaques.pruned$tip.label))
#tipcol[c(2,10,15)]<-'red3'
tipcol[c(2,3,10,12,13,15,16)]<-'gray62'
plot(macaques.pruned,type="phylogram", 
     edge.color=edgecol, edge.width=edgewidth,edge.lty=edgelty,
     label.offset=0, font=3, cex=0.7,
     tip.color=tipcol, adj=1)

