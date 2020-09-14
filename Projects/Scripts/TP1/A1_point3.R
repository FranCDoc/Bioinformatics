# Franco Chiesa Docampo

library(ape)
library(seqinr)
#reading the sequence 
list <- read.GenBank("NC_002662", as.character = TRUE); 
seq <- list[[1]]; 
N <- length(seq);
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
dimer_obs_freq <- (count(seq, word = 2))/N;
p <- (count(seq, word = 1))/N;
odd_ratios <- c(dimer_obs_freq[c(1)]/(p[c(1)]*p[c(1)]),dimer_obs_freq[c(2)]/(p[c(1)]*p[c(2)]),dimer_obs_freq[c(3)]/(p[c(1)]*p[c(3)]),dimer_obs_freq[c(4)]/(p[c(1)]*p[c(4)]),
                dimer_obs_freq[c(5)]/(p[c(2)]*p[c(1)]),dimer_obs_freq[c(6)]/(p[c(2)]*p[c(2)]),dimer_obs_freq[c(7)]/(p[c(2)]*p[c(3)]),dimer_obs_freq[c(8)]/(p[c(2)]*p[c(4)]),
                dimer_obs_freq[c(9)]/(p[c(3)]*p[c(1)]),dimer_obs_freq[c(10)]/(p[c(3)]*p[c(2)]),dimer_obs_freq[c(11)]/(p[c(3)]*p[c(3)]),dimer_obs_freq[c(12)]/(p[c(3)]*p[c(4)]),
                dimer_obs_freq[c(13)]/(p[c(4)]*p[c(1)]),dimer_obs_freq[c(14)]/(p[c(4)]*p[c(2)]),dimer_obs_freq[c(15)]/(p[c(4)]*p[c(3)]),dimer_obs_freq[c(16)]/(p[c(4)]*p[c(4)]))
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
exp_freqs <- c((p[c(1)]*p[c(1)]),(p[c(1)]*p[c(2)]),(p[c(1)]*p[c(3)]),(p[c(1)]*p[c(4)]),
               (p[c(2)]*p[c(1)]),(p[c(2)]*p[c(2)]),(p[c(2)]*p[c(3)]),(p[c(2)]*p[c(4)]),
               (p[c(3)]*p[c(1)]),(p[c(3)]*p[c(2)]),(p[c(3)]*p[c(3)]),(p[c(3)]*p[c(4)]),
                (p[c(4)]*p[c(1)]),(p[c(4)]*p[c(2)]),(p[c(4)]*p[c(3)]),(p[c(4)]*p[c(4)]))
