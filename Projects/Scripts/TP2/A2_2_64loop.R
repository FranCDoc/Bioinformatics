# Franco Chiesa Docampo

library(ape);library(Biostrings);library(seqinr);library(ape)

mat1 <- read.fasta("/Users/Franco/Desktop/Bioinformatics/Assignment2/AEB0.fasta",seqtype = "AA")
s1 <- getSequence(mat1[[1]]);s1 <- toupper(paste(s1,collapse = ""))
n1 <- length(getSequence(mat1[[1]]))

mat2 <- read.fasta("/Users/Franco/Desktop/Bioinformatics/Assignment2/receptors.fasta",seqtype = "AA")

nchar_pwa <- rep(0,64)
score_pwa <- rep(0,64)
sc_ps1 <- matrix(0, nrow = 64, ncol = 500)
sc_ps2 <- matrix(0, nrow = 64, ncol = 500)
counter <- 0;

for (j in 1:64){
  counter <- counter +1;
  print(counter)
  
  s2 <- getSequence(mat2[[j]]);s2 <- toupper(paste(s2,collapse = ""))
  n2 <- length(getSequence(mat2[[j]]))
  
  PAM300 <- read.table("/Users/Franco/Desktop/Bioinformatics/Assignment2/PAM300");PAM300 <- as.matrix(PAM300)
  colnames(PAM300)[colnames(PAM300) == "X."] <- "*"
  # pattern (s2) => FAMILY sequences; subject (s1) => AEB0 sequence
  l <-pairwiseAlignment(s2,s1,substitutionMatrix = PAM300, gapOpening = 11, gapExtension = 1, type = "local-global", scoreOnly = FALSE) # the lower the gapOpening the more "permission" has the algorithm to move between columns and rows in order to search for the max score. Naturally there will be more matches.
  nchar_pwa[[j]] <- nchar(l);
  score_pwa[[j]] <- score(l);
  print(type(l));print(l) # nchar gives the alignment size
  print('SCORE:');print(score(l))
  print('ALIGNED:');print(aligned(l)) # te arma el string de length(subject) (formado con aminoacidos del pattern) y lo pone en paralelo con el subject.
  print('PERCENTAGE OF ID:');print(pid(l))
  print('NMATCHES:');print(nmatch(l))
  print('NMISMATCH:');print(nmismatch(l))
  print('NCHAR:');print(nchar(l));
  print('///////////////////////////////////////////////////////////// - ///////////////////////////////////////////////////////////////')
  # Matches between an arbitrary pair of sequences may be due to chance (for example mutations, deletions or insertions) rather than homology ⇒ statistical significance required
  
  for (i in 1:500){
    permu_s1 <- sample(getSequence(mat1[[1]]),length(getSequence(mat1[[1]]))) 
    permu_s1 <- toupper(paste(permu_s1,collapse = ""))
    l1 <-pairwiseAlignment(s2,permu_s1,substitutionMatrix = PAM300, gapOpening = 11, gapExtension = 1, type = "local-global") # the lower the gapOpening the more "permission" has the algorithm to move between columns and rows in order to search for the max score. Naturally there will be more matches.
    sc_ps1[[j,i]] <- score(l1)/nchar(l1)
  }
  
  for (i in 1:500){
    permu_s2 <- sample(getSequence(mat2[[j]]),length(getSequence(mat2[[j]]))) 
    permu_s2 <- toupper(paste(permu_s2,collapse = ""))
    l2 <-pairwiseAlignment(permu_s2,s1,substitutionMatrix = PAM300, gapOpening = 11, gapExtension = 1, type = "local-global") # the lower the gapOpening the more "permission" has the algorithm to move between columns and rows in order to search for the max score. Naturally there will be more matches.
    sc_ps2[[j,i]] <- score(l2)/nchar(l2)
  }
}
sc_pstot <- cbind(sc_ps1,sc_ps2) # 64x1000 matrix

perc_scoreXY <- rep(0,64)
p_value <- rep(0,64)

counter <- 0;
for (k in 1:64){
  counter <- counter +1;
  print(counter)
  sort_sc <- 0;
  #for (l in 1:1000){
    sort_sc <- sort(sc_pstot[k,])
    percentile <- ecdf(sort_sc)
    perc_scoreXY[[k]] <- percentile(score_pwa[[k]]/nchar_pwa[[k]])
    p_value[[k]] <- 1-perc_scoreXY[[k]]
  #}
} 
count_simil <- 0;
for (k in 1:64){
  if (p_value[[k]] < 0.05){ # p-value threshold 
    count_simil <- count_simil + 1;
  }
}
