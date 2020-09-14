# Franco Chiesa Docampo

library(ape)
library(Biostrings)

#reading the sequence
list <- read.GenBank("NC_011725", as.character=TRUE)
seq <- list[[1]]
N<-length(seq)

#------------------------------------------------------------------------------#
#FUNCTION TO FIND ORF BASED ON DNA STRING 

ORF_finder <- function(DNA_string){
  #match with each start codon
  match_S1 <- matchPattern("ttg",DNA_string)
  match_S2 <- matchPattern("ctg",DNA_string)
  match_S3 <- matchPattern("ata",DNA_string)
  match_S4 <- matchPattern("att",DNA_string)
  match_S5 <- matchPattern("atc",DNA_string)
  match_S6 <- matchPattern("atg",DNA_string)
  match_S7 <- matchPattern("gtg",DNA_string)
  
  #match with each stop codon
  match_E1 <- matchPattern("tga",DNA_string)
  match_E2 <- matchPattern("taa",DNA_string)
  match_E3 <- matchPattern("tag",DNA_string)
  
  #merging and sorting the start position of all start codons
  start_pos <- c()
  start_pos <- c(start_pos, start(match_S1), start(match_S2), start(match_S3), start(match_S4), start(match_S5), start(match_S6), start(match_S7))
  start_pos <- sort(start_pos)
  
  #merging and sorting the start position of all stop codons
  end_pos <- c()
  end_pos <- c(end_pos, start(match_E1), start(match_E2), start(match_E3))
  end_pos <- sort(end_pos)
  
  #classification of start codons in each frame (F1, F2, F3)
  start_F1 <- c()
  start_F2 <- c()
  start_F3 <- c()
  
  for (i in 1:length(start_pos)) {
    if (start_pos[i]%%3 == 1) {
      # frame <- 1
      start_F1 <- c(start_F1, start_pos[i])
    }
    else if (start_pos[i]%%3 == 0) {
      #frame <- 2
      start_F2 <- c(start_F2, start_pos[i])
    }
    else {
      #frame <- 3
      start_F3 <- c(start_F3, start_pos[i])
    }
  }
  
  #classification of stop codons in each frame (F1, F2, F3)
  end_F1 <- c()
  end_F2 <- c()
  end_F3 <- c()
  
  for (i in 1:length(end_pos)) {
    if (end_pos[i]%%3 == 1) {
      # frame <- 1
      end_F1 <- c(end_F1, end_pos[i])
    }
    else if (end_pos[i]%%3 == 0) {
      #frame <- 2
      end_F2 <- c(end_F2, end_pos[i])
    }
    else {
      #frame <- 3
      end_F3 <- c(end_F3, end_pos[i])
    }
  }
  
  #getting the start and length of ORF in frame 1
  ORF_1_length  <- c()
  ORF_1_start <- c()
  i <- 1
  j <- 1
  flag <- 0
  
  while (i <= length(start_F1) && j <= length(end_F1)){
    if (end_F1[j] > start_F1[i] && flag == 0){
      ORF_1_length <- c(ORF_1_length, (end_F1[j] - start_F1[i]))
      ORF_1_start <- c(ORF_1_start, start_F1[i])
      i <- i + 1
      flag <- 1
    }
    else if (end_F1[j] > start_F1[i] && flag == 1){
      i <- i + 1
    }
    else {
      j <- j + 1
      flag <- 0
    }
  }
  
  #getting the start and length of ORF in frame 2
  ORF_2_length  <- c()
  ORF_2_start <- c()
  i <- 1
  j <- 1
  flag <- 0
  
  while (i <= length(start_F2) && j <= length(end_F2)){
    if (end_F2[j] > start_F2[i] && flag == 0){
      ORF_2_length <- c(ORF_2_length, (end_F2[j] - start_F2[i]))
      ORF_2_start <- c(ORF_2_start, start_F2[i])
      i <- i + 1
      flag <- 1
    }
    else if (end_F2[j] > start_F2[i] && flag == 1){
      i <- i + 1
    }
    else {
      j <- j + 1
      flag <- 0
    }
  }
  
  #getting the start and length of ORF in frame 3
  ORF_3_length  <- c()
  ORF_3_start <- c()
  i <- 1
  j <- 1
  flag <- 0
  
  while (i <= length(start_F3) && j <= length(end_F3)){
    if (end_F3[j] > start_F3[i] && flag == 0){
      ORF_3_length <- c(ORF_3_length, (end_F3[j] - start_F3[i]))
      ORF_3_start <- c(ORF_3_start, start_F3[i])
      i <- i + 1
      flag <- 1
    }
    else if (end_F3[j] > start_F3[i] && flag == 1){
      i <- i + 1
    }
    else {
      j <- j + 1
      flag <- 0
    }
  }
  
  #merging the lengths of the ORF
  ORF_length <- c(ORF_1_length, ORF_2_length, ORF_3_length)
  ORF_start <- c(ORF_1_start, ORF_2_start, ORF_3_start)
  
  #sorting the length
  ORF_length_sort <- sort(ORF_length)
  
  #order the length
  ORF_length_order <- order(ORF_length)
  
  return(list(ORF_length_sort,ORF_length_order,ORF_start))
}
#------------------------------------------------------------------------------#
#READING DNA 

#reading as a string
DNA <- BString(paste(seq, collapse = ""))

#obtaining the reverse complement
DNA_rev_comp <- reverse(chartr("atcg","tagc",DNA))

#function with DNA
list_ORF_DNA <- ORF_finder(DNA)
ORF_DNA <- list_ORF_DNA[[1]] #lengths
ORF_DNA_order <- list_ORF_DNA[[2]] #orders
ORF_DNA_start <- list_ORF_DNA[[3]] #starts

#function with DNA reverse complement
list_ORF_DNA_rev_comp <- ORF_finder(DNA_rev_comp)
ORF_DNA_rev_comp <- list_ORF_DNA_rev_comp[[1]] #lengths
ORF_DNA_rev_comp_order <- list_ORF_DNA_rev_comp[[2]] #orders
ORF_DNA_rev_comp_start <- list_ORF_DNA_rev_comp[[3]] #Starts

#merging ORF starts from DNA and reverse complement
ORF_DNA_Total_start <- c(ORF_DNA_start,ORF_DNA_rev_comp_start)

#merging and sorting ORF lengths from DNA and reverse complement
ORF_DNA_Total <- c(ORF_DNA,ORF_DNA_rev_comp)
ORF_DNA_T <- sort(ORF_DNA_Total)

#order of the ORF 
ORF_DNA_T_order <- order(ORF_DNA_Total)

#------------------------------------------------------------------------------#
#READING PERMUTATION OF DNA

#random permutation of DNA
DNA_perm <- sample(DNA, length(DNA))

#obtaining the reverse complement of the permutation
DNA_rev_comp_perm <- reverse(chartr("atcg","tagc",DNA_perm))

#function with DNA permutation
list_ORF_DNA_perm <- ORF_finder(DNA_perm)
ORF_DNA_perm <- list_ORF_DNA_perm[[1]] #lengths

#function with DNA reverse complement permutation
list_ORF_DNA_rev_comp_perm <- ORF_finder(DNA_rev_comp_perm)
ORF_DNA_rev_comp_perm <- list_ORF_DNA_rev_comp_perm[[1]] #lengths

#merging and sorting ORF from DNA and reverse complement permutation
ORF_DNA_Total_perm <- c(ORF_DNA_perm,ORF_DNA_rev_comp_perm)
ORF_DNA_T_perm <- sort(ORF_DNA_Total_perm)

#------------------------------------------------------------------------------#
#COUNTING LENGTH OF ORF

#counting number of ORF with specific length for DNA
counter10 <- 0
counter50 <- 0
counter100 <- 0
counter300 <- 0
counter500 <- 0
counter_plus <- 0

for (i in 1:length(ORF_DNA_T)){
  if (ORF_DNA_T[i] >= 10){
    counter10 <- counter10 + 1
  }
  if (ORF_DNA_T[i] >= 50){
    counter50 <- counter50 + 1
  }
  if (ORF_DNA_T[i] >= 100){
    counter100 <- counter100 + 1
  }
  if (ORF_DNA_T[i] >= 300){
    counter300 <- counter300 + 1
  }
  if (ORF_DNA_T[i] >= 500){
    counter500 <- counter500 + 1
  }
  else {
    counter_plus <- counter_plus + 1
  }
}

#counting number of ORF with specific length for DNA permutation
counter10_p <- 0
counter50_p <- 0
counter100_p <- 0
counter300_p <- 0
counter500_p <- 0
counter_plus_p <- 0

for (i in 1:length(ORF_DNA_T_perm)){
  if (ORF_DNA_T_perm[i] >= 10){
    counter10_p <- counter10_p + 1
  }
  if (ORF_DNA_T_perm[i] >= 50){
    counter50_p <- counter50_p + 1
  }
  if (ORF_DNA_T_perm[i] >= 100){
    counter100_p <- counter100_p + 1
  }
  if (ORF_DNA_T_perm[i] >= 300){
    counter300_p <- counter300_p + 1
  }
  if (ORF_DNA_T_perm[i] >= 500){
    counter500_p <- counter500_p + 1
  }
  else {
    counter_plus_p <- counter_plus_p + 1
  }
}

#------------------------------------------------------------------------------#

#find max length in DNA permutation
max_length_perm <- ORF_DNA_T_perm[length(ORF_DNA_T_perm)]

#find how many ORF in DNA are longer than this
counter_ORF <- 0

for (i in 1:length(ORF_DNA_T)) {
  if (ORF_DNA_T[i] > max_length_perm){
    counter_ORF <- counter_ORF + 1
  }
}

#find percentil 99% in DNA permutation
percent_99 <- ORF_DNA_T_perm[length(ORF_DNA_T_perm)*0.99] 

#find how many ORF in DNA are longer than this
counter_ORF_99 <- 0

for (i in 1:length(ORF_DNA_T)) {
  if (ORF_DNA_T[i] > percent_99){
    counter_ORF_99 <- counter_ORF_99 + 1
  }
}

#------------------------------------------------------------------------------#

#longest gene candidate
pos <- ORF_DNA_T_order[length(ORF_DNA_T)]
general_pos <- ORF_DNA_Total_start[pos]

#find string: normal = 1 -> DNA, normal = 0 -> reverse DNA
normal <- 0

for (i in 1:length(ORF_DNA)) {
  if (ORF_DNA_start[i] == general_pos){
    normal <- i
  }
  else {
    normal <- 0
  }
}

#find length and position of gene 
gene_length <- ORF_DNA_T[length(ORF_DNA_T)]
if (normal == 1){
  position <- ORF_DNA_order[length(ORF_DNA)]
  gene_position <- ORF_DNA_start[position]
}
if (normal == 0){
  position <- ORF_DNA_rev_comp_order[length(ORF_DNA_rev_comp)]
  gene_position <- ORF_DNA_rev_comp_start[position]
}

#find frame of correspondance 
frame <- 0
  
if (gene_position%%3 == 1) {
  frame <- 1
  }
if (gene_position%%3 == 0) {
  frame <- 2
}
if (gene_position%%3 == 2) {
  frame <- 3
}





