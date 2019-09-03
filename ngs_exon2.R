############################################################
## alpha 1 and alpha 2 analysis - exon 2
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

# load alpha-1 data
# DT1 <- fread(file.path('data/alpha-1 pooled.csv'))
# unique(DT1$PCR)
# DT1 <- DT1[PCR=='MM2F_S3R']
# write.csv(DT1,'data/mm2f_alpha1_pooled.csv', row.names = F)
DT1 = fread(file.path('data/mm2f_alpha1_pooled.csv'))
head(DT1)
DT1 = DT1[keyDT,on='subject',nomatch=0]

#load alpha-2 data
# DT2 <- fread(file.path('data/alpha-2 pooled.csv'))
# DT2 <- DT2[PCR=='MM2F_S3R']
# write.csv(DT2,'data/mm2f_alpha2_pooled.csv', row.names = F)

DT2 = fread(file.path('data/mm2f_alpha2_pooled.csv'))
DT2 = DT2[keyDT,on='subject',nomatch=0] # match barcode to macaque id 
head(DT2)
#########
#alpha-1#
#########
head(DT1)
#examine reads/sample
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
# 13 number of alleles left in cleaned


#########
#alpha-2#
#########

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

#heatmaps

#add missing zeros
plotDT2 <- fDT2[CJ(subject= fDT2[,unique(subject)], uniqueseq=fDT2[,unique(uniqueseq)]),on=c('subject','uniqueseq')]
plotDT2[is.na(freq),c('freq','count'):=0]

p <- ggplot(plotDT2,aes(x=as.factor(subject),y=uniqueseq,fill=freq))
p1 <- p + geom_tile(colour = "black") + scale_fill_gradient(low = "white",high = "steelblue")
p2 <- p1 +
	heatmap_theme +
	ggtitle('alpha-2 MM2F_S3R')
ggsave(file.path('plots','heatmaps','alpha-2 MM2F_S3R heatmap.pdf'),p2,width = 18, height = 4.8)

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
	geom_hline(yintercept=83.3,linetype='dashed')

p1 +
	#theme(axis.text.x=element_text(angle=90,vjust=0.45)) +
	ggtitle('alpha-2 MM2F_S3R') + 
	facet_wrap(~subject,nrow=6) + 
	blank_theme

#saved as 'alpha-2 frequency stacked barcharts.svg'
#Generally looks great, but a few require closer attention: Mf1_18, Mf1_24, Mf1_61, Mf1_62, Mf1_74
fullDT2 <- fread(file.path('data','alpha-2 pooled.csv'))
setnames(fullDT2,c('count','subject'),c('count','barcode'))
fullDT2 <- fullDT2[keyDT,on='barcode',nomatch=0]
setkey(fullDT2,subject)

fullDT2['Mf1_18'][freq>3]
#     barcode      PCR  uniqueseq count      freq subject
#  1:       8  F3F_S1R  F3F_S1R.3  1571 25.762545  Mf1_18
#  2:       8  F3F_S1R  F3F_S1R.4  1466 24.040669  Mf1_18
#  3:       8  F3F_S1R F3F_S1R.13  1565 25.664152  Mf1_18
#  4:       8  F3F_S1R F3F_S1R.14  1479 24.253854  Mf1_18
#  5:       8  F3F_S2R  F3F_S2R.1    59 22.097378  Mf1_18
#  6:       8  F3F_S2R  F3F_S2R.7   130 48.689139  Mf1_18
#  7:       8  F3F_S2R F3F_S2R.14     7  2.621723  Mf1_18 ? very low read count
#  8:       8  F3F_S2R F3F_S2R.24    65 24.344569  Mf1_18
#  9:       8 MM2F_S3R MM2F_S3R.1  1771 46.617531  Mf1_18
# 10:       8 MM2F_S3R MM2F_S3R.2  1700 44.748618  Mf1_18
# 11:       8 MM2F_S3R MM2F_S3R.5   267  7.028165  Mf1_18 ***
# 12:       8  S4F_S5R  S4F_S5R.3  1106 21.716081  Mf1_18
# 13:       8  S4F_S5R  S4F_S5R.4  1113 21.853524  Mf1_18
# 14:       8  S4F_S5R S4F_S5R.11  2692 52.856862  Mf1_18

fullDT2['Mf1_24'][freq>2]
#     barcode      PCR  uniqueseq count     freq subject
#  1:      17  F3F_S1R  F3F_S1R.1   394 25.27261  Mf1_24
#  2:      17  F3F_S1R  F3F_S1R.2   369 23.66902  Mf1_24
#  3:      17  F3F_S1R  F3F_S1R.7   407 26.10648  Mf1_24
#  4:      17  F3F_S1R  F3F_S1R.8   378 24.24631  Mf1_24
#  5:      17  F3F_S2R  F3F_S2R.1    41  4.85782  Mf1_24 ***
#  6:      17  F3F_S2R  F3F_S2R.5   774 91.70616  Mf1_24
#  7:      17 MM2F_S3R MM2F_S3R.1   534 56.03358  Mf1_24 ***
#  8:      17 MM2F_S3R MM2F_S3R.2   387 40.60860  Mf1_24
#  9:      17  S4F_S5R  S4F_S5R.3   565 48.37329  Mf1_24
# 10:      17  S4F_S5R  S4F_S5R.4   561 48.03082  Mf1_24

fullDT2['Mf1_61'][freq>2]
#     barcode      PCR  uniqueseq count      freq subject
#  1:      59  F3F_S1R  F3F_S1R.3   455 26.149425  Mf1_61
#  2:      59  F3F_S1R  F3F_S1R.4   407 23.390805  Mf1_61
#  3:      59  F3F_S1R F3F_S1R.11   431 24.770115  Mf1_61
#  4:      59  F3F_S1R F3F_S1R.12   437 25.114943  Mf1_61
#  5:      59  F3F_S2R  F3F_S2R.7   255 53.459119  Mf1_61
#  6:      59  F3F_S2R F3F_S2R.17   114 23.899371  Mf1_61
#  7:      59  F3F_S2R F3F_S2R.18    95 19.916143  Mf1_61 
#  8:      59 MM2F_S3R MM2F_S3R.1   446 48.373102  Mf1_61
#  9:      59 MM2F_S3R MM2F_S3R.2   396 42.950108  Mf1_61
# 10:      59 MM2F_S3R MM2F_S3R.5    62  6.724512  Mf1_61 ***
# 11:      59  S4F_S5R  S4F_S5R.3   324 23.841060  Mf1_61
# 12:      59  S4F_S5R  S4F_S5R.4   296 21.780721  Mf1_61
# 13:      59  S4F_S5R S4F_S5R.11   700 51.508462  Mf1_61

fullDT2['Mf1_62'][freq>2]
#     barcode      PCR  uniqueseq count      freq subject
#  1:      60  F3F_S1R  F3F_S1R.1  1453 25.388782  Mf1_62
#  2:      60  F3F_S1R  F3F_S1R.2  1460 25.511096  Mf1_62
#  3:      60  F3F_S1R  F3F_S1R.3  1448 25.301415  Mf1_62
#  4:      60  F3F_S1R  F3F_S1R.4  1360 23.763760  Mf1_62
#  5:      60 MM2F_S3R MM2F_S3R.1  1578 44.463229  Mf1_62
#  6:      60 MM2F_S3R MM2F_S3R.2  1586 44.688645  Mf1_62
#  7:      60 MM2F_S3R MM2F_S3R.5   311  8.763032  Mf1_62 ***
#  8:      60  S4F_S5R  S4F_S5R.3  1217 23.345482  Mf1_62
#  9:      60  S4F_S5R  S4F_S5R.4  1124 21.561481  Mf1_62
# 10:      60  S4F_S5R S4F_S5R.11  2686 51.525034  Mf1_62

fullDT2['Mf1_74'][freq>2]
#     barcode      PCR  uniqueseq count     freq subject
#  1:      74  F3F_S1R  F3F_S1R.1  1630 25.61685  Mf1_74
#  2:      74  F3F_S1R  F3F_S1R.2  1459 22.92944  Mf1_74
#  3:      74  F3F_S1R  F3F_S1R.5   769 12.08549  Mf1_74
#  4:      74  F3F_S1R  F3F_S1R.6   742 11.66117  Mf1_74
#  5:      74  F3F_S1R  F3F_S1R.7   896 14.08141  Mf1_74
#  6:      74  F3F_S1R  F3F_S1R.8   841 13.21704  Mf1_74
#  7:      74 MM2F_S3R MM2F_S3R.1   729 11.54029  Mf1_74  *** At least 8 copies!
#  8:      74 MM2F_S3R MM2F_S3R.2   753 11.92022  Mf1_74
#  9:      74 MM2F_S3R MM2F_S3R.5  3177 50.29286  Mf1_74
# 10:      74 MM2F_S3R MM2F_S3R.7  1543 24.42615  Mf1_74
# 11:      74  S4F_S5R  S4F_S5R.1  1266 24.77010  Mf1_74
# 12:      74  S4F_S5R  S4F_S5R.2  1218 23.83095  Mf1_74
# 13:      74  S4F_S5R  S4F_S5R.3   646 12.63941  Mf1_74
# 14:      74  S4F_S5R  S4F_S5R.4   617 12.07200  Mf1_74
# 15:      74  S4F_S5R S4F_S5R.24   623 12.18940  Mf1_74
# 16:      74  S4F_S5R S4F_S5R.25   630 12.32635  Mf1_74


############################################################
##
library(scales)

snp_a1 = c("MM2F_S3R.10", "MM2F_S3R.19","MM2F_S3R.5")
plotDT1$snp = ifelse(plotDT1$uniqueseq %in% snp_a1, 1, 0)
plotDT1_snp = plotDT1[!freq == 0,]
sum_a1 = ddply(plotDT1_snp, c("subject"), summarise,
               min_freq = min(freq)/100
)
plotDT1_snpb = plotDT1_snp[!snp == 0,]
sum_a1b = ddply(plotDT1_snpb, c("subject"), summarise,
               q1_freq = sum(freq)/100
)
a1 <- merge(sum_a1, sum_a1b, all = TRUE)
a1[is.na(a1)] <- 0
plot(1/a1$min_freq,a1$q_freq, xlim = c(1,4),
     bty = 'n', xlab = 'minimum copy number',
     ylab = 'proportion of sequences with pos 71 SNP')

snp_a2 = c("MM2F_S3R.1", "MM2F_S3R.11","MM2F_S3R.15")
plotDT2$snp = ifelse(plotDT2$uniqueseq %in% snp_a2, 1, 0)
plotDT2_snp = plotDT2[!freq == 0,]
sum_a2 = ddply(plotDT2_snp, c("subject"), summarise,
               min_freq = min(freq)/100
)
plotDT2_snpb = plotDT2_snp[!snp == 0,]
sum_a2b = ddply(plotDT2_snpb, c("subject"), summarise,
                q2_freq = sum(freq)/100
)
a2 <- merge(sum_a2, sum_a2b, all = TRUE)
a2[is.na(a2)] <- 0
a2$copy_n = round(1/a2$min_freq, 0)
a2$copy_n2 = ifelse(a2$copy_n %% 2 == 0,a2$copy_n,a2$copy_n*2)

plot(jitter(a2$copy_n2),a2$q_freq, xlim = c(0,14),
     bty = 'n', xlab = 'minimum copy number',
     ylab = 'proportion of sequences with pos 71 SNP',
     pch=19, col = alpha ('black', 0.4), cex = 2)

round(1/a2$min_freq, 0)
both = merge(a2,a1,by.x ='subject',by.y='subject',all=T)
plot(both$q1_freq,both$q2_freq,
  bty = 'n', xlab = 'proportion of alpha1 with pos 71 SNP',
  ylab = 'proportion of alpha2 with pos 71 SNP',
  pch=19, col = alpha ('black', 0.1), cex = 2)

DT2_2 = fullDT2[PCR=='MM2F_S3R']
DT2_2['Mf1_43'][freq>2]
