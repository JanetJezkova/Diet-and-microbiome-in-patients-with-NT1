library(dada2)

```{r message=FALSE, warning=FALSE}
dupl.concensus<-function(PHYLOS,NAMES){
  
  # exclude nonduplicated samples
  IDS<-as.character(data.frame(sample_data(PHYLOS))[,NAMES])
  IDS.dupl<-IDS[duplicated(IDS)]
  
  PHYLOSEQ<-prune_samples(IDS%in%IDS.dupl, PHYLOS)
  if(length(IDS.dupl)*2<length(IDS)) {NONUPLICATED<-prune_samples(!IDS%in%IDS.dupl, PHYLOS)
  print(paste("Following names are nonduplicated",sample_names(NONUPLICATED)))}
  
  CATS<-as.character(data.frame(sample_data(PHYLOSEQ))[,NAMES])
  CATS2<-levels(factor(CATS))
  OTU_TAB<-otu_table(PHYLOSEQ)
  rownames(OTU_TAB)<-CATS
  
  # i<-5
  for (i in 1:length(CATS2))
  {
    # print(CATS2[i])
    FILTER.act<-colSums(OTU_TAB[rownames(OTU_TAB)==CATS2[i],]>0)>1
    OTU_TAB[rownames(OTU_TAB)==CATS2[i],]
    OTU_TAB[rownames(OTU_TAB)==CATS2[i],]<-t(apply(OTU_TAB[rownames(OTU_TAB)==CATS2[i],],1,function(x) x*FILTER.act))
  }
  
  rownames(OTU_TAB)<-sample_names(PHYLOSEQ)
  otu_table(PHYLOSEQ)<-OTU_TAB
  PHYLOSEQ.clean<-prune_taxa(taxa_sums(PHYLOSEQ)>0,PHYLOSEQ)
  
  PHYLOSEQ.clean
}

merge.duplicates<-function(PHYLOSEQ,NAMES){
  CATS<-as.character(data.frame(sample_data(PHYLOSEQ))[,NAMES])
  sample_data(PHYLOSEQ)$duplic.id<-CATS
  SAMDAT<-sample_data(PHYLOSEQ)
  SAMDAT.sub<-subset(SAMDAT,duplicated(CATS)==F)
  FASTA<-refseq(PHYLOSEQ)
  rownames(SAMDAT.sub)<-SAMDAT.sub$duplic.id
  PHYLOSEQ.merge<-merge_samples(PHYLOSEQ,"duplic.id")
  sample_data(PHYLOSEQ.merge)<-SAMDAT.sub
  PHYLOSEQ.merge<-merge_phyloseq(PHYLOSEQ.merge,FASTA)
  PHYLOSEQ.merge
}

```


setwd("W:/R návody/Narko R/poslední skripty/ořezané sekvence+skripty122022/OREZANE_NEFILT/OREZANE_NEFILT")
PATH<-"W:/R návody/Narko R/poslední skripty/ořezané sekvence+skripty122022/OREZANE_NEFILT/OREZANE_NEFILT"
#List of forward and reverse reads
LIST<-list.files()
F_reads<-LIST[grep("-pair1.fastq.gz",LIST)]
R_reads<-LIST[grep("-pair2.fastq.gz",LIST)]
F_reads_TO<-paste0("S",F_reads)
R_reads_TO<-paste0("S",R_reads)

file.rename(from = F_reads,to=F_reads_TO)
file.rename(from = R_reads,to=R_reads_TO)

LIST<-list.files()
F_reads<-LIST[grep("-pair1.fastq.gz",LIST)]
R_reads<-LIST[grep("-pair2.gz",LIST)]

# #graphical representation of quality profiles
# system("zcat *_trus[12]R_trus-trimmed-pair1.fastq.gz > Merged.fastq")
# QP.f<-plotQualityProfile("Merged.fastq",aggregate = TRUE)+ggtitle("Forward reads")
# system("zcat *_trus[12]R_trus-trimmed-pair2.fastq.gz > Merged.fastq")
# QP.2<-plotQualityProfile("Merged.fastq",aggregate = TRUE)+ggtitle("Rewerse reads")
# system("rm Merged.fastq")
# 
# QP.f
# QP.2

# ggsave(QP.f,filename = "/media/kreising/DATA/data/Radka_Janet/Forward_reads.pdf")
# ggsave(QP.2,filename = "/media/kreising/DATA/data/Radka_Janet/Reverse_reads.pdf")

sample.names<-gsub("_trus-trimmed-pair1.fastq.gz","",F_reads)
sample.names<-gsub("-assigned-","",sample.names)
filtFs <- paste0(sample.names, "_READ1_filt.fastq.gz")
filtRs <- paste0(sample.names, "_READ2_filt.fastq.gz")

#Quality filtering
for(x in 1:length(F_reads)) {
  print(sample.names[x])
  fastqPairedFilter(c(F_reads[x], R_reads[x]), c(filtFs[x], filtRs[x]),
                    maxN=0, maxEE=2, minQ=2,truncQ=2,
                    compress=TRUE, verbose=TRUE,
                    minLen = c(270,190),truncLen = c(270,190))
}


###############################################
#DADA DENOISING################################
###############################################

#These commands denoise quality-filtered fastq files and build abundance matrix,
#(samples in rows, ASVs in columns)


#List of quality filtered fastq files
fns <- list.files()
fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) 

fnFs <- fastqs[grepl("_trus[12]R_READ1_filt.fastq.gz", fastqs)] 
fnRs <- fastqs[grepl("_trus[12]R_READ2_filt.fastq.gz", fastqs)] 
sample.names <- gsub("_READ1_filt.fastq.gz","",fnFs)

#fastq dereplication
derepFs <- derepFastq(fnFs,n = 1e+05, verbose=T)
derepRs <- derepFastq(fnRs,n = 1e+05, verbose=T)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#deoising
dadaFs <- dada(derepFs, selfConsist = TRUE,MAX_CONSIST=20)
dadaRs <- dada(derepRs, selfConsist = TRUE,MAX_CONSIST=20)

#merge denoised forward and reverse ASVs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 10,maxMismatch=1,justConcatenate=F)

#abundance matrix
seqtab <- makeSequenceTable(mergers)

save(seqtab,file = "W:/R návody/Narko R/poslední skripty/ořezané sekvence+skripty122022/OREZANE_NEFILT/FILTROVANÉ")

