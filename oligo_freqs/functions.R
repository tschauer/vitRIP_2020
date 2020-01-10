

##################################################

keepBSgenomeSequences <- function(genome, seqnames)
{
      stopifnot(all(seqnames %in% seqnames(genome)))
      genome@user_seqnames <- setNames(seqnames, seqnames)
      genome@seqinfo <- genome@seqinfo[seqnames]
      genome
}


##################################################

getFreq <- function(my_ranges, my_width = 4, my_step = 1, prob = TRUE){
      
      my_seq = BSgenome::getSeq(my_genome, my_ranges)
      
      mcols(my_ranges) <- cbind(mcols(my_ranges), oligonucleotideFrequency(my_seq, width = my_width, step = my_step, as.prob = prob))
      
      return(my_ranges)
}


##################################################


exportOligofromFreqRanges <- function(my_freq_ranges, oligo = "TTTT"){
      
      my_freq_ranges_oligo <- my_freq_ranges[,oligo]
      colnames(mcols(my_freq_ranges_oligo)) <- "score"
      
      export.bw(my_freq_ranges_oligo, con = paste("ranges_bins_freq",oligo,".bw", sep=""))
      
}

##################################################