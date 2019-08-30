#creating a phylogeny of primates with alpha globin stats
setwd("/Users/cfaust/Dropbox/Princeton/Projects/SimianMalaria/AlphaGlobinOX/macaquephylogeny")
library(ape)
library(phytools)

#import binida edmonds mammalian super tree
mammalsbest <- read.nexus("Bininda-emonds_2007_mammals_best.nex")
#plot(mammalsbest)
#mammalsbest$tip.label
macaques<- c("Macaca_arctoides","Macaca_assamensis","Macaca_thibetana","Macaca_radiata","Macaca_sinica","Macaca_cyclopis","Macaca_mulatta","Macaca_fuscata",
            "Macaca_fascicularis","Macaca_maura","Macaca_nigra","Macaca_ochreata","Macaca_tonkeana","Macaca_nemestrina","Macaca_silenus","Macaca_sylvanus")
macaques.pruned<-drop.tip(mammalsbest,mammalsbest$tip.label[-match(macaques, mammalsbest$tip.label)])
plot(macaques.pruned)
nodelabels()
edgelabels()
tiplabels()
edgecol <- rep('black', length(macaques.pruned$edge.length))
edgecol[c(17,26,9,3)] <- 'red3'
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



#10K has more species

macaques10k<- read.nexus("consensusTree_10kTrees_macaca_down_vs2.nex")
                         consensusTree_10kTrees_macaca_down_12Jan2014.nex")
plot(macaques10k)
#plot with different edge colors for macaques with AQ
