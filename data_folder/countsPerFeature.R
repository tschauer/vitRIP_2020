


rm(list = ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



########################################################  libraries ######################################################## 


library(GenomicFeatures)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)
library(RColorBrewer)
library(DESeq2)
library(affy)
library(pheatmap)
library(scales)
library(genefilter)
library(ggplot2)
library(sva)
library(grid)
library(GenomicFeatures)
library(gdata)
library(gridExtra)
library(gclus)
library(GenomicAlignments)
library(rtracklayer)






############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################  annotation ######################################################## 



dir.gtf <- "../genome/"
gtffile <- file.path(dir.gtf,"dmel-all-r6.17.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")


my_genes <- genes(txdb)


my_genes$symbol <- mapIds(org.Dm.eg.db, my_genes$gene_id, "SYMBOL", keytype="FLYBASE", multiVals="first")
my_genes <- my_genes[!(grepl("^His[1-4]",my_genes$symbol))]

my_genes$symbol <- gsub("lncRNA:", "", my_genes$symbol)



my_gtf <- import.gff(gtffile)



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


ranges_tRNA <-   reduce(my_gtf[grep("tRNA",my_gtf$gene_symbol)])
ranges_snRNA <-  reduce(my_gtf[grep("snRNA",my_gtf$gene_symbol)])
ranges_snoRNA <- reduce(my_gtf[grep("snoRNA",my_gtf$gene_symbol)])
ranges_rRNA <-   reduce(my_gtf[grep("rRNA",my_gtf$gene_symbol)])

ranges_roX1 <- my_genes[grep("roX1", my_genes$symbol)]
ranges_roX2 <- my_genes[grep("roX2", my_genes$symbol)]

ranges_threeUTR <- reduce(my_gtf[grep("3UTR",my_gtf$type)])
ranges_fiveUTR <-  reduce(my_gtf[grep("5UTR",my_gtf$type)])

ranges_exon <- reduce(my_gtf[grep("exon",my_gtf$type)])

ranges_intron <- setdiff(my_genes, 
                         c(ranges_exon, ranges_rRNA, ranges_snoRNA,ranges_snRNA,ranges_tRNA,ranges_fiveUTR, ranges_threeUTR), ignore.strand=FALSE)

ranges_exon <- setdiff(ranges_exon, 
                       c(ranges_rRNA, ranges_snoRNA,ranges_snRNA,ranges_tRNA,ranges_fiveUTR, ranges_threeUTR), ignore.strand=FALSE)



export.bed(ranges_threeUTR, "../data_folder/ranges_threeUTR.bed")
export.bed(ranges_fiveUTR, "../data_folder/ranges_fiveUTR.bed")
export.bed(ranges_intron, "../data_folder/ranges_intron.bed")
export.bed(ranges_exon, "../data_folder/ranges_exon.bed")

export.bed(ranges_tRNA, "../data_folder/ranges_tRNA.bed")
export.bed(ranges_snoRNA, "../data_folder/ranges_snoRNA.bed")
export.bed(ranges_snRNA, "../data_folder/ranges_snRNA.bed")





############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 





bam_dir <- "../BAM/"
bam_files <- list.files(bam_dir, pattern = ".bam$")
bam_files_path <- file.path(bam_dir, bam_files)







############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 







my_ranges <-  c("ranges_snoRNA", "ranges_snRNA", "ranges_tRNA","ranges_exon", "ranges_intron", "ranges_fiveUTR", "ranges_threeUTR")

ranges_counts <- parallel::mclapply(seq_along(bam_files), mc.cores = 16, FUN = function(bidx){
      
      reads_bam <- readGAlignmentPairs(bam_files_path[bidx])
      
      ranges_sum_counts <- lapply(seq_along(my_ranges), FUN = function(ridx){
            
            my_range <- get(my_ranges[ridx])
            
            overlap_counts <- countOverlaps(query = my_range,
                                            subject = invertStrand(reads_bam),
                                            ignore.strand=FALSE)
            
            sum_counts <- sum(overlap_counts)
            
            return(sum_counts)
            
      })
      
      names(ranges_sum_counts) <- gsub("ranges_","", my_ranges)
      
      ranges_sum_counts <- unlist(ranges_sum_counts)
      
      return(ranges_sum_counts)
})



names(ranges_counts) <- gsub("_[G,A,T,C].*","", bam_files)


ranges_mat <- matrix(unlist(ranges_counts),ncol = length(ranges_counts))

colnames(ranges_mat) <- names(ranges_counts)
rownames(ranges_mat) <- names(ranges_counts[[1]])



save(list = "ranges_mat", file = "../data_folder/ranges_mat.rda")









############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 
















