# Franco Chiesa Docampo

#reading the sequence

library(ape)
list <- read.GenBank("NC_002662", as.character = TRUE)
seq <- list[[1]]
#------------------------------------------------------------------------------------------------------------------------#
# 4 arguments to pass to function: sequence to analyze, nucleotide of interest, size of the sliding window, step of the sliding window
# t = 6 [s]
N <- length(seq);
kmer_content <- function(seq,nucleotide_1,nucleotide_2,window_s,step_s,graph_color,graph_title){
  
  positions <- 1:length(seq);
  count <- function(s,c){
    length(which(s==c))
  }
  
  start = 0;
  frame = 0;
  count_nucleotide = 0;
  num_frames <- floor((N-window_s)/step_s)+1;
  seq_window = rep(0L, num_frames);
  
  for (i in 0:num_frames){
    start <- i*step_s+1; 
    frame <- seq[c(start:(start+window_s-1))]
    count_nucleotide <- count(frame,nucleotide_1) + count(frame,nucleotide_2)
    seq_window[c(i+1)] <- count_nucleotide
    frame <- 0;
  }
  ratio_window = seq_window/window_s;
  positionsSeqWindow <- 1:length(seq_window); 
  
  plot(ratio_window ~ positionsSeqWindow,type="l",col=graph_color,ylab=graph_title)
  plot(ratio_window[c(250:350)] ~ positionsSeqWindow[c(250:350)],type="l",col=graph_color,ylab=graph_title) 
  plot(ratio_window[c(980:1180)] ~ positionsSeqWindow[c(980:1180)],type="l",col=graph_color,ylab=graph_title)
  
  return(ratio_window) 
}
#------------------------------------------------------------------------------------------------------------------------#
GC_content <- kmer_content(seq,"c","g",6000,2000,"red","GC content") 
AT_content <- kmer_content(seq,"a","t",6000,2000,"green","AT content")
#------------------------------------------------------------------------------------------------------------------------#



