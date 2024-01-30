library(dada2)
library(phyloseq)
library(ShortRead)

load("/media/kreising/DATA/data/Radka_Janet/DADA_PHYLOSEQ/seqtab.R")
setwd("/media/kreising/DATA/data/Radka_Janet/DADA_PHYLOSEQ/")
#########################################
#Chimeric sequences######################
#########################################

#extraxt ASVs fasta from abundance matrix
FASTA<-DNAStringSet(colnames(seqtab))
names(FASTA)<-colnames(seqtab)
writeFasta(FASTA,"haplo.fasta")

# #elimination of chimeric sequences by uchime (Terminal command)
# system("usearch8.0.1517_i86linux32 -uchime_ref haplo.fasta -db ~/DB/gold.fasta -nonchimeras haplo.uchime.fasta -strand plus")

#uchime nebode fungovat. Ale muze se pouzit neco takoveho:
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

############################################
#TAXONOMY###################################
############################################
# databaze se da stahnout tady: https://benjjneb.github.io/dada2/training.html

FASTA<-readDNAStringSet("haplo.uchime.fasta")
taxa <- assignTaxonomy(as.character(FASTA), 
                       refFasta="~/DB/DADA2/silva_nr99_v138_train_set.fa.gz", 
                       multithread=8,minBoot = 80)

############################################
#Create phyloseq object#####################
############################################

#OTU TABLE
seqtab<-otu_table(seqtab,taxa_are_rows = F)

#HAPLO
HAPLO<-readDNAStringSet("haplo.fasta")

#TAXO
TAXO<-tax_table(taxa)

PHYLOSEQ<-merge_phyloseq(seqtab,TAXO,HAPLO)
PHYLOSEQ_dupl<-PHYLOSEQ
sample_names(PHYLOSEQ)


#Provizorni sample metadata
SN<-sample_names(PHYLOSEQ_dupl)
SN.mod<-gsub("_F_trus[12]R","",SN)

DF<-data.frame(ID_D=SN,ID=SN.mod)
DF<-sample_data(DF)
sample_names(DF)<-SN

sample_data(PHYLOSEQ_dupl)<-DF
save(PHYLOSEQ_dupl,file = "/media/kreising/DATA/data/Radka_Janet/DADA_PHYLOSEQ/PHYLOSEQ_dupl.R")

# load("/media/kreising/DATA/data/Radka_Janet/DADA_PHYLOSEQ/PHYLOSEQ_dupl.R")
#write.table
