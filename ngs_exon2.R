############################################################
## alpha 1 and alpha 2 processing - exon 2
############################################################

# libraries
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyr)

# load barcode key
keyDT <- fread('data/sample_barcode_relative.csv',check.names=TRUE)
keyDT[,c('barcode_seq','relative'):=NULL]
setnames(keyDT,c('Sample.ID','barcode_ID'),c('subject_id','subject'))

#########
# HBA1 cleaning#
# load alpha-1 data
# DT1 <- fread(file.path('data/alpha-1 pooled.csv'))
# unique(DT1$PCR)

# write.csv(DT1,'data/mm2f_alpha1_pooled.csv', row.names = F)
DT1 = fread(file.path('data/alpha-1 pooled.csv'))
DT1 <- DT1[PCR=='MM2F_S3R']
head(DT1)
DT1 = DT1[keyDT,on='subject',nomatch=0]
p = ggplot(DT1[,.(value=sum(value)),by=subject],aes(value))
p
p+geom_histogram()
p+geom_histogram(bins=40) + xlim(0,2000)
#there are two low count samples, one with <500 reads, and one with probably not much over 100.
setkey(DT1,value)
DT1[,.(value=sum(value)),by=subject_id][value<1000]
head(DT1)
# subject_id value
# 1:     Mf1_42   426
# 2:     Mf1_51    29
DT1_counts <- DT1[freq>2]
DT1_counts = DT1_counts[,.(count=sum(value)),by=subject_id]
mean(DT1_counts$count)
#eliminate the those less than 100
DT1_trun = DT1[subject_id!='Mf1_51']
#DT1_trun = DT1_trun[subject_id!='Mf1_42']

#length(unique(DT1_trun$subject_id))
p = ggplot(DT1_trun,aes(freq))
p+geom_histogram()
p+geom_histogram() + xlim(0,2) +ylim(0,20) #there's a cluster under ~1.7%, with a fairly normal distribution

setkey(DT1_trun,freq)
DT1_trun[value>0,freq] #check by eyeballing - a 2% cut is fine (the gap is between 1.7% and 2.6%)
fDT1 <- DT1_trun[freq>2]
fDT1[,sum(freq)]/DT1_trun[,sum(freq)]#how many reads were removed?
# 0.9867369 #<2% removed
HBA1 = read.csv("data/seq_conversion_HBA1.csv", header = T)
fDT1 = fDT1[HBA1,on='uniqueseq',nomatch=0] # match barcode to macaque id 
fDT1$converted_name = ordered(fDT1$converted_name, 
                              levels = c("HBA1.1", "HBA1.2", "HBA1.3",
                                         "HBA1.4", "HBA1.5", "HBA1.6",
                                         "HBA1.7", "HBA1.8", "HBA1.9",
                                         "HBA1.10", "HBA1.11", "HBA1.12",
                                         "HBA1.13"))
head(fDT1)
fDT1 = fDT1[DT1_counts,on='subject_id',nomatch=0]
fDT1$freq_wo2 = fDT1$value/fDT1$count

# frequency table 
head(fDT1)
fDT1_wide = spread(fDT1[,c('subject_id','converted_name','freq_wo2')],converted_name, freq_wo2)
write.csv(fDT1_wide, 'output/alpha1_exon2_freqtable_less2andlowcountremoved.csv', row.names = F)
length(fDT1[,unique(uniqueseq)])
# 13 number of alleles left in cleaned

#########
#alpha-2#
#########
#load alpha-2 data
# DT2 <- fread(file.path('data/alpha-2 pooled.csv'))
# DT2 <- DT2[PCR=='MM2F_S3R']
# write.csv(DT2,'data/mm2f_alpha2_pooled.csv', row.names = F)

DT2 = fread(file.path('data/alpha-2 pooled.csv'))
DT2 <- DT2[PCR=='MM2F_S3R']
DT2 = DT2[keyDT,on='subject',nomatch=0] # match barcode to macaque id 
head(DT2)

#examine reads/sample
p <- ggplot(DT2[,.(count=sum(count)),by=subject],aes(count))
p+geom_histogram()
p+geom_histogram(bins=20) + xlim(0,2500)
p+geom_histogram(bins=20) + xlim(0,1000)
#there is only one very low sample :
DT2[,.(count=sum(count)),by=subject_id][count<500]
# subject_id count
# 1:     Mf1_42   320
# 2:     Mf1_51     9
# 3:     Mf1_59   349

DT2_counts <- DT2[freq>2]
DT2_counts = DT2_counts[,.(count=sum(count)),by=subject_id]
mean(DT2_counts$count)
#let's eliminate the very lowest
head(DT2)
DT2_trun = DT2[subject_id!='Mf1_51']
#DT2_trun = DT2_trun[subject_id!='Mf1_42']
#DT2_trun = DT2_trun[subject_id!='Mf1_59']
#DT2_trun = DT2_trun[subject_id!='Mf1_42']
#DT2_trun = DT2_trun[subject_id!='Mf1_59']
#DT2_trun = DT2_trun[subject_id!='Mf1_61']

#frequency of individual sequences
p <- ggplot(DT2_trun,aes(freq))
p+geom_histogram()
p+geom_histogram() + xlim(0,5) +ylim(0,200) 
#eliminate everything under 2%
setkey(DT2_trun,freq)
DT2_trun[count>0,freq] #check by eyeballing - a 2% cut is fine (the gap is between 1.7% and 2.6%)
fDT2 = DT2_trun[freq>2]

#how many reads were removed?
fDT2[,sum(freq)]/DT2_trun[,sum(freq)]# [1] 0.9787715#<3% removed
length(fDT2[,unique(uniqueseq)]) #12 unique sequences left
length(fDT2[,unique(subject)])# [1] 75 macaques left
fDT2[,.(freq=sum(freq)/length(fDT2[,unique(subject)])),by=uniqueseq] #mean freq
fDT2[,.(freq=max(freq)),by=uniqueseq] #max freq
#all look real (lowest max is 23%)
#by contrast:
DT2[,.(freq=max(freq)),by=uniqueseq]

HBA2 = read.csv("data/seq_conversion_HBA2.csv", header = T)
fDT2 = fDT2[HBA2,on='uniqueseq',nomatch=0] # match barcode to macaque id
unique(fDT2$converted_name)
fDT2$converted_name = ordered(fDT2$converted_name, 
                              levels = c("HBA2.1", "HBA2.2", "HBA2.3",
                                         "HBA2.4", "HBA2.5", "HBA2.6",
                                         "HBA2.7", "HBA2.8", "HBA2.9",
                                         "HBA2.10", "HBA2.11", "HBA2.12"))

head(fDT2)
fDT2 = fDT2[DT2_counts,on='subject_id',nomatch=0]
fDT2$freq_wo2 = fDT2$count/fDT2$i.count

fDT2_wide = spread(fDT2[,c('converted_name','freq_wo2','subject_id')],converted_name, freq_wo2)
head(fDT2_wide)
write.csv(fDT2_wide, 'output/alpha2_exon2_freqtable_less2andlowcountremoved.csv', row.names = F)

#########
# visualising trimmed data
#stacked barcharts, with ablines at 25% intervals (solid) and  16.7% intervals (dashed)
p = ggplot(fDT1,aes(x=as.factor(subject_id),y=converted_name,fill=freq))
p1 = p + geom_tile(colour = "black") + scale_fill_gradient(low = "white",high = "steelblue") #heatmap
#3; 5; 10
blue13 = c('grey95', 'grey90', '#0570b0',
           'grey80', '#74a9cf','grey75',
           'grey70', 'grey65', 'grey60', '#bdc9e1',
           'grey55', 'grey50', 'grey45')
fDT1$freq_wo2 = fDT1$freq_wo2*100
fDT1$converted_name
p = ggplot(fDT1,aes(x=factor(1),y=freq_wo2,fill=converted_name))
p1 = p + 
  geom_bar(stat = "identity",  position = "stack")   +
  scale_fill_manual("unique sequences",values=blue13) +
  geom_hline(yintercept=25,linetype='solid') +
  geom_hline(yintercept=50,linetype='solid') +
  geom_hline(yintercept=75,linetype='solid') +
  geom_hline(yintercept=16.7,linetype='dashed') +
  geom_hline(yintercept=33.3,linetype='dashed') +
  geom_hline(yintercept=66.7,linetype='dashed') +
  geom_hline(yintercept=83.3,linetype='dashed') 

p1 +#theme(axis.text.x=element_text(angle=90,vjust=0.45)) +
  facet_wrap(~subject_id,nrow=6) + 
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
## HBA2



head(fDT2)
fDT2$freq_wo2 = fDT2$freq_wo2*100
p = ggplot(fDT2,aes(x=factor(1),y=freq_wo2,fill=converted_name))
red13 = c('#67001f', 'grey95',  'grey90','#b2182b',  'grey80', 'grey75', 'grey70', 'grey65', 'grey60', 'grey55', 'grey50', '#d6604d')
p1 = p + 
	geom_bar(stat = "identity",  position = "stack")   +
	scale_fill_manual("unique sequences",values=red13) +
	geom_hline(yintercept=25,linetype='solid') +
	geom_hline(yintercept=50,linetype='solid') +
	geom_hline(yintercept=75,linetype='solid') +
	geom_hline(yintercept=16.7,linetype='dashed') +
	geom_hline(yintercept=33.3,linetype='dashed') +
	geom_hline(yintercept=66.7,linetype='dashed') +
	geom_hline(yintercept=83.3,linetype='dashed') 

p1 +#theme(axis.text.x=element_text(angle=90,vjust=0.45)) +
	facet_wrap(~subject_id,nrow=6) + 
	theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
