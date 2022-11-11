library(ape)
 
args <- commandArgs(TRUE)
in_file  <- args[1]
out_file <- args[2]

# threshold, could be modified, if needed; 0.98 means that columns with gaps at frequency 98% and higher will be removed
threshold <- 0.98

seq_aln <- read.dna(in_file, format = "fasta")
trimmed_fas <- del.colgapsonly(seq_aln, threshold = threshold)
write.dna(trimmed_fas, file = out_file, format = "fasta", nbcol = -1, colsep = "")
