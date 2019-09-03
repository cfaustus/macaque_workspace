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
#alpha-1#
# load alpha-1 data
# DT1 <- fread(file.path('data/alpha-1 pooled.csv'))
# unique(DT1$PCR)
# DT1 <- DT1[PCR=='MM2F_S3R']
# write.csv(DT1,'data/mm2f_alpha1_pooled.csv', row.names = F)
DT1 = fread(file.path('data/mm2f_alpha1_pooled.csv'))
head(DT1)
DT1 = DT1[keyDT,on='subject',nomatch=0]
p = ggplot(DT1[,.(value=sum(value)),by=subject],aes(value))
p
p+geom_histogram()
p+geom_histogram(bins=40) + xlim(0,2000)
#there are two low count samples, one with <500 reads, and one with probably not much over 100.
setkey(DT1,value)
DT1[,.(value=sum(value)),by=subject_id][value<1000]
# subject_id value
# 1:     Mf1_42   426
# 2:     Mf1_51    29
#eliminate the those less than 1000
DT1_trun = DT1[subject!='Mf1_42']
DT1_trun = DT1_trun[subject!='Mf1_51']

p = ggplot(DT1_trun,aes(freq))
p+geom_histogram()
p+geom_histogram() + xlim(0,2) +ylim(0,20) #there's a cluster under ~1.7%, with a fairly normal distribution

setkey(DT1_trun,freq)
DT1_trun[value>0,freq] #check by eyeballing - a 2% cut is fine (the gap is between 1.7% and 2.6%)
fDT1 <- DT1_trun[freq>2]
fDT1[,sum(freq)]/DT1_trun[,sum(freq)]#how many reads were removed?
# 0.986907 #<2% removed

# frequency table 
head(fDT1)
fDT1_wide = spread(fDT1[,c(3,5,6)],uniqueseq, freq)
write.csv(fDT1_wide, 'output/alpha1_exon2_freqtable_less2andlowcountremoved.csv', row.names = F)
length(fDT1[,unique(uniqueseq)])
fDT1$uniqueseq
# 13 number of alleles left in cleaned

#########
#alpha-2#
#########
#load alpha-2 data
# DT2 <- fread(file.path('data/alpha-2 pooled.csv'))
# DT2 <- DT2[PCR=='MM2F_S3R']
# write.csv(DT2,'data/mm2f_alpha2_pooled.csv', row.names = F)

DT2 = fread(file.path('data/mm2f_alpha2_pooled.csv'))
DT2 = DT2[keyDT,on='subject',nomatch=0] # match barcode to macaque id 
head(DT2)

#examine reads/sample
p <- ggplot(DT2[,.(count=sum(count)),by=subject],aes(count))
p+geom_histogram()
p+geom_histogram(bins=20) + xlim(0,2500)
p+geom_histogram(bins=20) + xlim(0,1000)
#there is only one very low sample :
DT2[,.(count=sum(count)),by=subject_id][count<1000]
## subject_id count
# 1:     Mf1_24   953
# 2:     Mf1_31   659
# 3:     Mf1_42   320
# 4:     Mf1_51     9
# 5:     Mf1_59   349
# 6:     Mf1_61   922

#let's eliminate the very lowest
DT2_trun = DT2[subject!='Mf1_24']
DT2_trun = DT2_trun[subject!='Mf1_31']
DT2_trun = DT2_trun[subject!='Mf1_42']
DT2_trun = DT2_trun[subject!='Mf1_51']
DT2_trun = DT2_trun[subject!='Mf1_59']
DT2_trun = DT2_trun[subject!='Mf1_61']
#frequency of individual sequences
p <- ggplot(DT2_trun,aes(freq))
p+geom_histogram()
p+geom_histogram() + xlim(0,5) +ylim(0,200) 
#eliminate everything under 2%
setkey(DT2_trun,freq)
DT2_trun[count>0,freq] #check by eyeballing - a 2% cut is fine (the gap is between 1.7% and 2.6%)
fDT2 <- DT2_trun[freq>2]

head(fDT2)
fDT2_wide = spread(fDT2[,c(3,5,6)],uniqueseq, freq)
write.csv(fDT2_wide, 'data/alpha2_exon2_freqtable_less2andlowcountremoved.csv', row.names = F)

#how many reads were removed?
fDT2[,sum(freq)]/DT2_trun[,sum(freq)]# [1] 0.9786663#<3% removed
length(fDT2[,unique(uniqueseq)]) #12 unique sequences left
length(fDT2[,unique(subject)])# [1] 78 macaques left
fDT2[,.(freq=sum(freq)/78),by=uniqueseq] #mean freq
fDT2[,.(freq=max(freq)),by=uniqueseq] #max freq
#all look real (lowest max is 23%)
#by contrast:
DT2[,.(freq=max(freq)),by=uniqueseq]


#########
#heatmaps
#stacked barcharts, with ablines at 25% intervals (solid) and  16.7% intervals (dashed)
p <- ggplot(fDT2,aes(x=factor(1),y=freq,fill=uniqueseq))
red13 = c('#67001f', 'grey95', '#b2182b', 'grey90', '#d6604d', 'grey80', 'grey75', 'grey70', 'grey65', 'grey60', 'grey55', 'grey50')
p1 <- p + 
	geom_bar(stat = "identity",  position = "stack")   +
	scale_fill_manual(values=red13) +
	geom_hline(yintercept=25,linetype='solid') +
	geom_hline(yintercept=50,linetype='solid') +
	geom_hline(yintercept=75,linetype='solid') +
	geom_hline(yintercept=16.7,linetype='dashed') +
	geom_hline(yintercept=33.3,linetype='dashed') +
	geom_hline(yintercept=66.7,linetype='dashed') +
	geom_hline(yintercept=83.3,linetype='dashed') +
  theme(axis.title.x = '')

p1 +
	#theme(axis.text.x=element_text(angle=90,vjust=0.45)) +
	ggtitle('alpha-2 MM2F_S3R') + 
	facet_wrap(~subject,nrow=6) + 
	theme_minimal()
