# Franco Chiesa Docampo

library(ape)
library(Biostrings)
library(seqinr)

#set working directory
#loading functions from pHMM.R
source("/Users/marcelochiesa/Desktop/EXAMS\ LLN/LGBIO2010/Bioinformatics/Assignment2/pHMM.R")

#reading the sequence
list <- read.fasta("/Users/marcelochiesa/Desktop/EXAMS\ LLN/LGBIO2010/Bioinformatics/Assignment2/AEB0.fasta", seqtype = "AA")
AEBO <- toupper(getSequence(list[[1]]))
N<-length(AEBO)

#converting the file into a readable one
convertHMMER3("/Users/marcelochiesa/Desktop/EXAMS\ LLN/LGBIO2010/Bioinformatics/Assignment2/7tm_1.hmm","/Users/marcelochiesa/Desktop/EXAMS\ LLN/LGBIO2010/Bioinformatics/Assignment2/pHMM")
pHMM <- readHMMFile("/Users/marcelochiesa/Desktop/EXAMS\ LLN/LGBIO2010/Bioinformatics/Assignment2/pHMM")

#performing viterbi algorithm
viterbi_1 <-viterbi(pHMM,AEBO)
sc_ps <- rep(0,10)
counter <- 0;
for (i in 1:10){
  permu_AEBO <- sample(AEBO,length(AEBO)) 
  viterbi_2 <-viterbi(pHMM,permu_AEBO)
  sc_ps[[i]] <- viterbi_2[["score"]][["log.values"]]
  counter <- counter + 1;
  print(counter)
}

perc_scoreXY <- 0
p_value <- 0

sort_sc <- rep(0,10)
sort_sc <- sort(sc_ps)
percentile <- ecdf(sort_sc)
perc_scoreXY <- percentile(viterbi_1[["score"]][["log.values"]])
p_value <- 1-perc_scoreXY


count_simil <- 0;
if (p_value < 0.05){ # p-value threshold 
  count_simil <- count_simil + 1;
}
