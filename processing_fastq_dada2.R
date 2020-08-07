# Processing Illumina fastq files
# Sequencing data in SRA: PRJNA639946
# Authored by Stephen Preston
# dada2 v1.12 used

library(dada2)
library(Biostrings)
library(DECIPHER)
library(data.table)

#a function to print the quality profiles of supplied fastq read files (using dada2)
checkqual <- function(fns){
    for(fn in fns) print(dada2::plotQualityProfile(fn))
}


#a function to call dada2 to trim and apply EE filtering, dereplicate, build an error model, and print the error plots to screen
run_dada2_pt1 <- function(fnFs,fnRs,sample.names,filt_path,truncLen,maxEE=c(2,2),truncQ=2,rdsdir=FALSE) {
    
    #check that filt_path doesn't exist, and create it
    ensuredir(filt_path)

    #construct full paths for filtered files
    filtFs <- file.path(filt_path, paste0(sample.names, "_1_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, "_2_filt.fastq.gz"))
    
    #trim and apply EE filter
    for(i in seq_along(fnFs)) { #iterate through file pairs
      dada2::fastqPairedFilter(
                                c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                                truncLen=truncLen, 
                                maxN=0, maxEE=maxEE, truncQ=truncQ, 
                                compress=TRUE, verbose=TRUE
                              )
    }
    
    #Dereplication
    derepFs <- dada2::derepFastq(filtFs, verbose=TRUE)
    derepRs <- dada2::derepFastq(filtRs, verbose=TRUE)
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    #Learn the Error Rates
    errF <- dada2::dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
    errR <- dada2::dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
    
    #look at error plots
    print(dada2::plotErrors(errF[[1]], nominalQ=TRUE))
    print(dada2::plotErrors(errR[[1]], nominalQ=TRUE))
    
    return(list(errF=errF,errR=errR,derepFs=derepFs,derepRs=derepRs))
}


#a function to call dada2 to run inference, merge pairs, and construct sequence table.  Optionally, save intermediate stages as RDS
run_dada2_pt2 <- function(derepFs,derepRs,errF,errR,sample.names,pool,rdsdir=FALSE,multithread=TRUE) {
    
    #if an rdsdir has been specified (to save the intermediate stages), check that it doesn't exist
    if(is.character(rdsdir)){
        ensuredir(rdsdir)
    }

    #Sample Inference
    dadaFs <- dada2::dada(derepFs, err=errF[[1]]$err_out, multithread=multithread,pool=pool)
    dadaRs <- dada2::dada(derepRs, err=errR[[1]]$err_out, multithread=multithread,pool=pool)

    #Merge paired reads
    mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

    #Construct the sequence table
    seqtab <- dada2::makeSequenceTable(mergers[names(mergers) != "Mock"])
    
    #save rds files for intermediate stages
    if(is.character(rdsdir)){
        saveRDS(errF,file.path(rdsdir,'errF.rds'))
        saveRDS(errR,file.path(rdsdir,'errR.rds'))
        saveRDS(dadaFs,file.path(rdsdir,'dadaFs.rds'))
        saveRDS(dadaRs,file.path(rdsdir,'dadaRs.rds'))
        saveRDS(mergers,file.path(rdsdir,'mergers.rds'))
    }
    
    return(seqtab)
}


#a simple function to create a directory, or to return an error if it already exists and is not empty.
ensuredir <- function(pth,pathID='directory'){
    if(dir.exists(pth)){
        f <- list.files(pth, all.files=TRUE, recursive=TRUE, full.names=TRUE)
        if(sum(file.info(f)$size)) stop(paste0('Specified filt_dir already exists, and is not empty.'))
    }else{
        dir.create(pth,recursive=TRUE)
    }
}


#a function to return a named vector of fastq file sizes, using the ShortRead package (which is required by dada2, anyway)
getsizes <- function(fns){
    sizes = ShortRead::qa(fns, type="fastq")[["readCounts"]]
    sizes_vector = sizes$read
    names(sizes_vector) = rownames(sizes)
    return(sizes_vector)
}


#a function to check that input files have pairs - returns a vector of the name-stem of those which do;
#prints a warning if any do not
checkpairings <- function(fns,seperator='_',readIDs=c('1','2'),suffix='.fastq'){   
    basenames = basename(fns) #get the filenames (seperate from the full path to the files)
    fns_with_sep <- grep(seperator,basenames,value=TRUE,fixed=TRUE) #exclude filenames which don't include the seperator
    
    #get the name stems - the names before the seperator
    #in case of the seperator appearing multiple times, split only on the final instance
    rootnames <- sapply(sapply(sapply(strsplit(basenames, seperator), function(x) x[1:length(x)-1]),paste,collapse=seperator),unlist)
    rootnames <- unique(rootnames)
    
    #check that both pair files are present for each name stem
    good_roots <- vector(mode='character')
    for(root in rootnames){
        present = vector(mode='logical')
        for(r in seq_along(readIDs)){
            #look for the pairfile
            #note that this relies on consistent formatting, and the correct (case sensitive) suffix
            if(paste0(root,seperator,readIDs[r],suffix) %in% basenames){
                present[r]=TRUE
            }else{
                present[r]=FALSE
                warning(paste0('Missing pair file: ',root,seperator,readIDs[r],suffix))
            }
        }
        if(all(present)) good_roots=c(good_roots,root) #add the root to the good_roots vector if all pairfiles are present
    }
    return(good_roots)
}


#####################################################
#constants for the experiment
seperator = '_'
readIDs=c('1','2')
suffix='.fastq'
maxEE=2

#directories and filenames
fastqdir <- 'F:/TCRSeq data/MiSeq23/19feb20/paired_barcodes_fastq'
analysisdir <- 'F:/TCRSeq data/MiSeq23/19feb20/for_paper/dada2'

libpaths = c('alpha1'='MiSeq23.11/CF_MM2F_CF_S3R','alpha2'='MiSeq23.12/CF_MM2F_CF_S3R')
seqtab_fns <- c('alpha1'='seqtab_alpha1_MM2F_S3R.csv','alpha2'='seqtab_alpha2_MM2F_S3R.csv')

#empty lists, to add to later
fnFs = vector(mode='list')
fnRs = vector(mode='list')
trunclengths = vector(mode='list')
inroots_keep = vector(mode='list')
filt_paths = vector(mode='list')
dada2_output = vector(mode='list')
rdsdir = vector(mode='list')

###########
#run dada2#
###########

#build paths and get input fastq file lists for each gene, excluding very small input files
for(gene in c('alpha1','alpha2')){
    #construct subdirectory paths
    indir <- file.path(fastqdir,libpaths[gene])
    filt_paths[[gene]] <- file.path(analysisdir,libpaths[gene],'filtered')
    rdsdir[[gene]] <- file.path(analysisdir,libpaths[gene],'rds')
    
    #get a list of root names of paired input files.  Exclude unpaired files, with a warning.  There shouldn't be any missing.
    inroots <- checkpairings(list.files(indir,pattern=".fastq$",full.names=TRUE),seperator,readIDs,suffix)
    
    #filter to get rid of very small files
    #as very small files tend to be junk, resulting in no reads after filtering, and a dada2 crash.
    
    infiles_1 <- paste0(inroots,seperator,readIDs[1],suffix)
    infiles_2 <- paste0(inroots,seperator,readIDs[2],suffix)
    sizes_1 <- getsizes(file.path(indir,infiles_1))#get read counts from read_1 files (should match read 2 files)
    
    plot(hist(sizes_1)) #print read length histogram to screen as a check
    
    reads_cutoff = 10 #the filter is set smaller than usual, as small numbers of reads may be informative here.
    keep_1 <- sizes_1 > reads_cutoff
    keep <- keep_1
    
    inroots_keep[[gene]] <- inroots[keep]
    fnFs[[gene]] = file.path(indir,paste0(inroots_keep[[gene]],seperator,readIDs[1],suffix))
    fnRs[[gene]] = file.path(indir,paste0(inroots_keep[[gene]],seperator,readIDs[2],suffix))
}

#print quality plots to screen (for the first 3 read_1 files, then the first 3 read _2 files)
#from these, select trimming lengths

checkqual(fnFs[['alpha1']][1:3])
checkqual(fnRs[['alpha1']][1:3])

checkqual(fnFs[['alpha2']][1:3])
checkqual(fnRs[['alpha2']][1:3])

#truncation lengths were obtained by running the code above
trunclengths[['alpha1']] = c(260,200)
trunclengths[['alpha2']] = c(250,200)

#run dada2 length truncation, EE filter, and error models, checking the error model plots for each gene
dada2_output[['alpha1']] <- run_dada2_pt1(
                              fnFs = fnFs[['alpha1']],
                              fnRs = fnRs[['alpha1']],
                              sample.names = inroots_keep[['alpha1']],
                              filt_path = filt_paths[['alpha1']],
                              truncLen = trunclengths[['alpha1']],
                              maxEE = maxEE,
                              rdsdir = rdsdir[['alpha1']]
                              )
                              
dada2_output[['alpha2']] <- run_dada2_pt1(
                              fnFs = fnFs[['alpha2']],
                              fnRs = fnRs[['alpha2']],
                              sample.names = inroots_keep[['alpha2']],
                              filt_path = filt_paths[['alpha2']],
                              truncLen = trunclengths[['alpha2']],
                              maxEE = maxEE,
                              rdsdir = rdsdir[['alpha2']]
                              )

#The error models look OK

#Run dada2 inference, merge pairs, remove chimeras, and save
for(gene in c('alpha1','alpha2')){
    seqtab = run_dada2_pt2(
                            derepFs = dada2_output[[gene]][['derepFs']],
                            derepRs = dada2_output[[gene]][['derepRs']],
                            errF = dada2_output[[gene]][['errF']],
                            errR = dada2_output[[gene]][['errR']],
                            sample.names = inroots_keep[[gene]],
                            rdsdir = rdsdir[[gene]],
                            pool=TRUE
                           )                       
    saveRDS(seqtab,file.path(rdsdir[[gene]],'seqtab_before_chimera_removal.rds')) #save "raw" sequence table as rds file  
    nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #remove chimeras
    saveRDS(nochim,file.path(rdsdir[[gene]],'seqtab_nochim.rds')) #save sequence table as rds file
    write.csv(nochim,file.path(analysisdir,seqtab_fns[gene])) #write sequence table as a csv file
    
    #free up memory
    rm(seqtab)
}
#Identified 5 bimeras out of 57 input sequences.
#Identified 12 bimeras out of 67 input sequences.

################################
#dada2 v1.12 doesn't combine sequences with their reverse complements.  So, to do so:
#a function to combine reverse complement sequences 
orientOTUtab.mat <- function(seqtab,threshold=0.05,processors=NULL){
    require(Biostrings)
    seqs = Biostrings::DNAStringSet(colnames(seqtab)) #convert sequences to a DNAStringSet
    rc_seqs = Biostrings::reverseComplement(seqs) #get reverse complements
    #now identify and combine duplicates    
    dup_seqs = seqs[seqs %in% rc_seqs] #identify sequences which have their reverse complement present
    for(n in seq_along(dup_seqs)){ #for each of these, add the reverse complement read count to the "parent" read count
        u_seqs = as.character(c(dup_seqs[n],Biostrings::reverseComplement(dup_seqs[n])))
        seqtab[,u_seqs[1]] = rowSums(seqtab[,u_seqs])
        seqtab[,u_seqs[2:length(u_seqs)]] = 0 #set the reverse complement read count to zero
    }
    storage.mode(seqtab) <- "integer"
    seqtab = seqtab[,colSums(seqtab)>0] #remove empty columns
    colnames(seqtab) = as.character(Biostrings::DNAStringSet(colnames(seqtab))) #restore names to sequence table
    return(seqtab)
}

#analysisdir <- 'C:/Users/Stephen Preston/Sync/Macaque/analysis_redone_for_paper/MM2_S3F_(for_paper)'
analysisdir <- 'F:/TCRSeq data/MiSeq23/19feb20/for_paper/dada2'

DTlist = vector(mode='list')
alleles_list = vector(mode='list')

for(i in seq_along(seqtab_fns)){ #process each gene in turn
    DTlist[[i]] <- data.table()
    alleles_list[[i]] <- vector(mode='character')
    temp_seqtab <- as.matrix(fread(file.path(analysisdir,seqtab_fns[i])),rownames=1) #load as matrix
    temp_oriented <- orientOTUtab.mat(temp_seqtab) #collapse reverse complements
    temp_alleles <- colnames(temp_oriented) #store sequences
    names(temp_alleles) <- paste0('HBA2',as.character(i),as.character(1:length(temp_alleles))) #rename alleles
    colnames(temp_oriented) <- names(temp_alleles) #rename columns in matrix to new sequence names
    tempDT <- as.data.table(temp_oriented,keep.rownames='subject') #convert to data.table
    tempDT[,PCR:=names(seqtab_fns)[i]] #create a column specifying the amplicon
    tempDTm <- melt(tempDT,id.vars=c('subject','PCR'),variable.name='allele',value.name='count') #melt
    tempDTm[,freq:=(count/sum(count))*100,by=subject] #calculate frequency column
    temp_alleles = as.character(DECIPHER::OrientNucleotides(DNAStringSet(temp_alleles),orientation='both')) #orient sequences
    DTlist[[i]] <- rbind(DTlist[[i]],tempDTm,fill=TRUE) #merge into master data.table for this gene
    alleles_list[[i]] = c(alleles_list[[i]],temp_alleles)#merge sequence IDs into master sequence vector for this gene
    rm(temp_seqtab,temp_oriented,temp_alleles,tempDT,tempDTm) #free up memory
    
    #######
    #tests#
    #######
    #check that no duplicates are present
    if(length(alleles_list[[i]])!=length(unique(alleles_list[[i]]))){stop('Duplicate sequences present.')}
    if(
       any(alleles_list[[i]] %in% as.character(
                                                    Biostrings::reverseComplement(
                                                                          Biostrings::DNAStringSet(alleles_list[[i]])
                                                                                        )
                                                )
           )
       ){
         stop('Reverse complements not correctly merged.')
        }
}

######
#save#
######
for(i in seq_along(seqtab_fns)){
    fwrite(DTlist[[i]],file.path(analysisdir,paste0('alpha-',as.character(i),'_MM2F_S3R_pooled.csv'))) #save sequence data.table
    fwrite(
            data.table(ID=names(alleles_list[[i]]),seq=alleles_list[[i]]),
            file.path(analysisdir,paste0('alpha-',as.character(i),'_MM2F_S3R_unique_seqs.csv'))
           ) #save sequences
}

