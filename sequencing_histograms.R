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


ggplot(hba1_aa, aes(x = number, y = value, fill = 'grey50'))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo) +
  theme(legend.position="none")

ggplot(hba2_aa, aes(x = number, y = value))+
  geom_bar(stat = 'identity', color = 'black', fill = 'grey50')+  
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo)+
  theme(legend.position="none")



##replace

###############
#### histograms for gene freq
library(ggplot2)
library(reshape2)
library(car)

##### seq richness #####
uni = read.csv("data/unique_seqs.csv", header = T)
hist(uni$unique)
ggplot(uni, aes(unique)) +
  geom_histogram(binwidth = 1, 
                 col="black", fill = 'grey50',
                 size=.5) +
  labs(x = 'number of unique sequences', y = 'number of macaques') +
  theme_classic(base_size = 13) +
  facet_wrap(~gene, ncol = 1)

##### protein richness #####
prot = read.csv("data/protein_pivot3.csv",header=T)
head(prot)
prot_long <- melt(prot, id.vars = c("loci", "aa_combo", "gel"))
prot_long$number <- recode(prot_long$variable, "'X0' = 0; 'X1' = 1; 'X2' = 2")
#prot_long$gel[prot_long$variable == 'X0'] = 'zero' 
hba1_aa = prot_long[prot_long$loci == "HBA1",]
hba2_aa = prot_long[prot_long$loci == "HBA2",]

# c('A','Q','P','AQ','AP','AQP',
#   'QP','AX','AQX','APX','QPX'), 

col_pal = c('goldenrod4','gold3','#335B8E','gold',
            '#B7DBDB','steelblue3', 'steelblue4','khaki3',
            '#B5B867','#6CA18F','olivedrab4')
ggplot(hba1_aa, aes(x = number, y = value, fill = gel))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of HBA1 unique sequences', y = 'number of macaques') +
  scale_fill_manual("gel phenotypes", values = c("A" = "goldenrod4", 
                                                 "Q" = "gold", "Q?" = "gold2", 
                                                 "X" = "#B5B867", "X?" = "#6CA18F")) +
  theme_classic(base_size = 13) +
  facet_grid(~aa_combo)

ggplot(hba2_aa, aes(x = number, y = value, fill = gel))+
  geom_bar(stat = 'identity', color = 'black')+  
  labs(x = 'number of HBA2 unique sequences', y = 'number of macaques') +
  scale_fill_manual("gel phenotypes", values = c("A" = "goldenrod4", 
                                                 "Q" = "gold1", "Q?" = "gold2", 
                                                 "X" = "#B5B867", "X?" = "#6CA18F")) +
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

f# A; AP, AQ, AQP AQX
pie_col = rev(c('goldenrod4', '#B7DBDB','gold','steelblue3','#B5B867'))

# c('A','Q','P','AQ','AP','AQP',
#   'QP','AX','AQX','APX','QPX'), 

col_pal = c('goldenrod4','gold3','#335B8E','gold',
            '#B7DBDB','steelblue3', 'steelblue4','khaki3',
            '#B5B867','#6CA18F','olivedrab4')
exp_obs = read.csv("data/phenotype_indonesia_data.csv", header =T)
#exp_obs$phenotype <- factor(exp_obs$phenotype, levels = c("AQX", "AQP", "AQ", "AP", "A"))
head(exp_obs)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

ggplot(exp_obs, aes(x = "",y = proportion, fill = phenotype)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + 
  facet_wrap(~study, ncol =1) +
  scale_fill_manual(values=pie_col) +
  blank_theme +
  theme(axis.text.x=element_blank()) 

indo = exp_obs[exp_obs$study == 'ours',] 
indo$phenotype <- factor(indo$phenotype, levels = rev(c('A','AQ','AX','AQX')))
indo$phenotype <- factor(indo$phenotype, levels = c('A','AQ','AX','AQX'))

indo_pal = c('#2166ac','#b2182b','#67a9cf','#ef8a62')
ggplot(indo, aes(x = "",y = proportion, fill = phenotype)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y",direction = 180) + 
  #facet_wrap(~study, ncol =1) +
  scale_fill_manual(values = c("A" = "#2166ac", 
                               "AQ" = "#b2182b", 
                               "AX" = '#67a9cf',
                               "AQX" = "#ef8a62")) +
  blank_theme +
  theme(legend.text=element_text(size=20)) +
  theme(axis.text.x=element_blank()) +
  theme(legend.title=element_text(size=22))
