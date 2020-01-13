library(ape)
library(Biostrings)
library(seqinr)
library(annotate)
library(XML)

#set working directory
setwd("~/2. UCL/2 semestre/Bioinformatics/Project 2")

#reading the sequence
AEBO <- read.fasta("AEB0.fasta", seqtype = "AA")
AEBO <- toupper(paste((getSequence(AEBO[[1]])), collapse = ""))
N<-length(AEBO)

#finding 1000 hits for sequence query
hits = blastSequences(x = paste(">ID-1\n",AEBO,sep = ""), database = "nr", program = "blastp", hitListSize = 1000, filter = "mL", timeout = 40, as="XML")

definitions = sapply(hits["//Hit_def"], xmlValue)

#matching hits to family message
contain <- c()
proportion <- 0
for (i in 1:1000){
  contain[i] <- grepl("G-protein coupled receptor", definitions[i], fixed=TRUE)
  if (contain[i] == TRUE){
    proportion = proportion + 1
  }
}

#reading the receptors
list <- read.fasta("receptors.fasta", seqtype = "AA")

#finding 1000 hits for each receptor query
contain_rec <- c()
proportion_rec <- rep(0,10)
for (i in 1:10){
  receptor <- toupper(paste((getSequence(list[[i]])), collapse = ""))
  N_r<-length(receptor)
  
  hits_rec = blastSequences(x = paste(">ID-1\n",receptor,sep = ""), database = "nr", program = "blastp", hitListSize = 1000, filter = "mL", timeout = 150, as="XML")
  
  definitions_rec = sapply(hits_rec["//Hit_def"], xmlValue)
  
  #matching hits to family message
  for (j in 1:1000){
    contain_rec[j] <- grepl("G-protein coupled receptor", definitions_rec[j], fixed=TRUE)
    if (contain_rec[j] == TRUE){
      proportion_rec[i] = proportion_rec[i] + 1
    }
  }
}
  
