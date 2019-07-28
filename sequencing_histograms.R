###############
#### histograms for gene freq
library(ggplot2)
library(reshape2)
library(car)

##### seq richness #####
uni = read.csv("data/unique_seqs.csv", header = T)
hist(uni$unique)
h = ggplot(uni, aes(unique)) +
  geom_histogram(binwidth = 1, 
                 col="black", fill = 'lightcyan3',
                 size=.5) +
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  theme_classic(base_size = 13) +
  theme( panel.background = element_rect(fill = "transparent"), # bg of the panel, 
         panel.grid.major = element_blank() ,# get rid of major grid
         panel.grid.minor = element_blank(), # get rid of minor grid
         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
         legend.box.background = element_rect(fill = "transparent") )+# get rid of legend panel bg
  facet_wrap(~gene, ncol = 1)
ggsave(h, filename = "plots/tr_tst2.png",  bg = "transparent")

##### protein richness #####
prot = read.csv("data/protein_pivot2.csv",header=T)
head(prot)
prot_long <- melt(prot, id.vars = c("loci", "aa_combo", "gel"))
prot_long$number <- recode(prot_long$variable, "'X0' = 0; 'X1' = 1; 'X2' = 2")

hba1_aa = prot_long[prot_long$loci == "HBA1",]
hba2_aa = prot_long[prot_long$loci == "HBA2",]

ggplot(hba1_aa, aes(x = number, y = value, fill = gel))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  scale_fill_manual("gel phenotypes", values = c("A" = "orange2", 
                                         "Q" = "seagreen3", "Q?" = "seagreen3", 
                                         "X" = "sienna3", "X?" = "sienna1")) +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo)

ggplot(hba2_aa, aes(x = number, y = value, fill = gel))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  scale_fill_manual("gel phenotypes", values = c("A" = "orange2", 
                                                 "Q" = "seagreen3", "Q?" = "seagreen3", 
                                                 "X" = "sienna3", "X?" = "sienna1")) +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo)


ggplot(hba1_aa, aes(x = number, y = value, fill = gel))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  scale_fill_manual("gel phenotypes", values = c("A" = "grey50", 
                                                "Q" = "grey50", "Q?" = "grey50", 
                                                 "X" = "grey50", "X?" = "grey50")) +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo) +
  theme(legend.position="none")

ggplot(hba2_aa, aes(x = number, y = value, fill = gel))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  scale_fill_manual("gel phenotypes", values = c("A" = "grey50", 
                                                  "Q" = "grey50", "Q?" = "grey50", 
                                                  "X" = "grey50", "X?" = "grey50")) +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo)+
  theme(legend.position="none")

