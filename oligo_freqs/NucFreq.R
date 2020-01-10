


rm(list = ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



#################################################################################################################################################################### 
#################################################################################################################################################################### 

########################################################  libraries ######################################################## 



library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(Biostrings)
library(GenomicRanges)
library(org.Dm.eg.db)

source("functions.R")


#################################################################################################################################################################### 
#################################################################################################################################################################### 

######################################################## dm6  annotation ######################################################## 



my_genome <- BSgenome.Dmelanogaster.UCSC.dm6

my_chromosomes <- c("chrX","chr2L","chr2R","chr3L","chr3R","chr4","chrM","chrY")

my_genome <- keepBSgenomeSequences(my_genome, my_chromosomes)



###########

my_tile_width <- 100


my_bins <- tileGenome(seqlengths = seqinfo(my_genome), 
                      tilewidth = my_tile_width, 
                      cut.last.tile.in.chrom = TRUE)



##########


my_gtf <- import.gff("../genome/dmel-all-r6.17.gtf")

#txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
#seqlevelsStyle(txdb) <- "Ensembl"

seqlevelsStyle(my_gtf) <- "UCSC"
seqlevels(my_gtf)[24] <- "chrM"

my_gtf <- keepSeqlevels(my_gtf, value = my_chromosomes, pruning.mode = "coarse")



my_genes    <-     my_gtf[grep("gene", my_gtf$type)]
ranges_exon <-     my_gtf[grep("exon",my_gtf$type)]
ranges_threeUTR <- my_gtf[grep("3UTR",my_gtf$type)]
ranges_fiveUTR <-  my_gtf[grep("5UTR",my_gtf$type)]


reduceRanges <- function(gr){
      grl <- reduce(split(gr, elementMetadata(gr)$gene_id))
      reducedRanges <- unlist(grl, use.names=FALSE)
      elementMetadata(reducedRanges)$gene_id <- rep(names(grl),elementNROWS(grl))
      reducedRanges$gene_symbol <- mapIds(org.Dm.eg.db, keys = reducedRanges$gene_id, keytype = "FLYBASE",column = "SYMBOL",multiVals = "first")
      return(reducedRanges)
}


my_genes <-        reduceRanges(my_genes)
ranges_exon <-     reduceRanges(ranges_exon)
ranges_threeUTR <- reduceRanges(ranges_threeUTR)
ranges_fiveUTR <-  reduceRanges(ranges_fiveUTR)


export.bed(my_genes,        con = "my_genes.bed")
export.bed(ranges_exon,     con = "ranges_exon.bed")
export.bed(ranges_threeUTR, con = "ranges_threeUTR.bed")
export.bed(ranges_fiveUTR,  con = "ranges_fiveUTR.bed")



#################################################################################################################################################################### 
#################################################################################################################################################################### 

######################################################## frequency ######################################################## 




ranges_bins_freq <- getFreq(my_bins, prob = FALSE)

seqlevels(ranges_bins_freq)[7] <- "chrmitochondrion_genome"


##################################################


exportOligofromFreqRanges(ranges_bins_freq, oligo = "TTTT")
exportOligofromFreqRanges(ranges_bins_freq, oligo = "AAAA")
exportOligofromFreqRanges(ranges_bins_freq, oligo = "GGGG")
exportOligofromFreqRanges(ranges_bins_freq, oligo = "CCCC")



##################################################



ranges_genes_freq <- getFreq(my_genes)
#ranges_exon_freq <-  getFreq(ranges_exon)
#ranges_fiveUTR_freq <-  getFreq(ranges_fiveUTR)


ranges_threeUTR_freq <-  getFreq(ranges_threeUTR)
ranges_threeUTR_freq <- ranges_threeUTR_freq[order(ranges_threeUTR_freq$TTTT, decreasing = TRUE),]
ranges_threeUTR_freq <- ranges_threeUTR_freq[!(duplicated(ranges_threeUTR_freq$gene_id)),]
ranges_threeUTR_freq[grep("Sxl",ranges_threeUTR_freq$gene_symbol), "TTTT"]


ranges_exon_freq <-  getFreq(ranges_exon)
ranges_exon_freq <- ranges_exon_freq[order(ranges_exon_freq$TTTT, decreasing = TRUE),]
ranges_exon_freq <- ranges_exon_freq[!(duplicated(ranges_exon_freq$gene_id)),]
ranges_exon_freq[grep("Sxl",ranges_exon_freq$gene_symbol), "TTTT"]



##################################################

#ranges_genes_4merfreq <- ranges_genes_freq[,c("gene_id","gene_symbol","TTTT","AAAA", "GGGG","CCCC")]

saveRDS(ranges_genes_freq, file = "ranges_genes_freq.rds")

saveRDS(ranges_threeUTR_freq, file = "ranges_threeUTR_freq.rds")

saveRDS(ranges_exon_freq, file = "ranges_exon_freq.rds")


#################################################################################################################################################################### 
#################################################################################################################################################################### 

######################################################## end ######################################################## 


