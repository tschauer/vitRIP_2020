




rm(list = ls())


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
library(BSgenome.Dmelanogaster.UCSC.dm6)

source("functions/functions_DESeq.R")


cbPalette <- c("#999999", "#0072B2", "#CC79A7", "#009E73", "#E69F00", "#D55E00", "#56B4E9", "#F0E442")




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################  annotation ######################################################## 




my_genome <- BSgenome.Dmelanogaster.UCSC.dm6
my_chromosomes <- c("chrX","chr2L","chr2R","chr3L","chr3R","chr4","chrM","chrY")
my_genome <- keepBSgenomeSequences(my_genome, my_chromosomes)



dir.gtf <- "genome/"
gtffile <- file.path(dir.gtf,"dmel-all-r6.17.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")

seqlevelsStyle(txdb) <- "UCSC"
seqlevels(txdb) <- gsub("mitochondrion_genome","chrM", seqlevels(txdb))


########################################################  


my_genes <- genes(txdb)

my_genes <- keepSeqlevels(my_genes, my_chromosomes, pruning.mode = "coarse")

my_genes$symbol <- mapIds(org.Dm.eg.db, my_genes$gene_id, "SYMBOL", keytype="FLYBASE", multiVals="first")



########################################################  


my_seq_genes <-  BSgenome::getSeq(my_genome, my_genes)
identical(names(my_seq_genes), my_genes$gene_id)

my_seq_genes <- my_seq_genes[width(my_seq_genes) > 100]
my_seq_genes <- my_seq_genes[width(my_seq_genes) < 1e5]


########################################################  


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 







############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################    results    ######################################################


my_ress <- list.files(path = "analysis/",pattern = "^res")

i=1

for(i in seq_along(my_ress)){
        
        my_name <- gsub(".txt","",my_ress[i])
        
        my_res <- read.delim(file.path("analysis/", my_ress[i]), header = T, row.names = 1)
        
        assign(my_name, my_res)
        
}


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

############################################                functions               ##########################################



bplotFreqs <- function(my_seq_set, 
                       my_seq_notset,
                       my_title = "",
                       #ylim = 500,                       
                       tophits = 25,
                       isFreq = FALSE,
                       my_width = 4){
        
        my_freqs1 <- oligonucleotideFrequency(my_seq_set, width = my_width, step = 1, as.prob = isFreq)
        
        my_order <- order(apply(my_freqs1,2, median), decreasing = T)
        my_freqs1 <- my_freqs1[,my_order]
        
        
        my_freqs2 <- oligonucleotideFrequency(my_seq_notset, width = my_width, step = 1, as.prob = isFreq)
        my_freqs2 <- my_freqs2[,my_order]
        
        stopifnot(identical(colnames(my_freqs1),colnames(my_freqs2)))
        
        colnames(my_freqs1) <- gsub("T","U", colnames(my_freqs1))
        colnames(my_freqs2) <- gsub("T","U", colnames(my_freqs2))
        
        my_freqs <- c(as.list(data.frame(my_freqs1)),
                      as.list(data.frame(my_freqs2)))
        
        #my_freqs <- my_freqs[order(names(my_freqs), decreasing = F)]
        
        if(isFreq){
                my_ylab <- "Frequency"
        } else {
                my_ylab <- "Counts"
        }
        
        
        bp <-boxplot(my_freqs[rep(1:tophits, each=2)+c(0,ncol(my_freqs1))], 
                     col =c("darkred","darkgrey"),
                     main = my_title,
                     las=2, outline=F, #ylim = c(0,ylim), 
                     ylab = my_ylab)
        
        legend("topright", legend = c(paste0("Gene set: ", length(my_seq_set)),
                                      paste0("Not in set: ", length(my_seq_notset))),
               fill = c("darkred","darkgrey"))
}





############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

############################################       OligoFreqs in Enriched Set      ##########################################





my_seq_types <- ls(pattern = "my_seq_genes")
my_ress <- ls(pattern = "^res")
my_widths <- c(4,5,6)

r=15
ty=1
my_width = 4



for(isFreq in c(TRUE, FALSE)){
        
        for(ty in seq_along(my_seq_types)){
                
                my_type_names <- gsub("my_seq_","",my_seq_types[ty])
                
                my_seq_type <- get(my_seq_types[ty])
                
                
                if(isFreq){
                        file_name <- paste0("analysis_OligoFreqs/OligoFreqs_Enriched_",my_type_names,".pdf")
                } else {
                        file_name <- paste0("analysis_OligoFreqs/OligoCounts_Enriched_",my_type_names,".pdf")
                }
                
                pdf(file = file_name, width = 16, height = 8, useDingbats = FALSE)
                
                par(oma=c(2,2,2,2), mar=c(5,5,5,2), mgp=c(3,1,0), bg=NA,
                    cex.axis = 1, cex.main = 1.5, cex.lab=1.25)
                
                par(mfcol=c(3,2))
                
                for(r in seq_along(my_ress)){
                        
                        my_res <- get(my_ress[r])
                        
                        my_enriched_ids <- rownames(my_res[my_res$padj < 0.01 & my_res$log2FoldChange > 0,])
                        
                        my_seq_set <-    my_seq_type[  names(my_seq_type) %in% my_enriched_ids]
                        my_seq_notset <- my_seq_type[!(names(my_seq_type) %in% my_enriched_ids)]
                        
                        if(length(my_seq_set) < 10){next()}
                        
                        for(my_width in my_widths){
                                
                                bplotFreqs(my_seq_set = my_seq_set, 
                                           my_seq_notset = my_seq_notset, 
                                           my_title = paste0("Enriched in Comparison: ", gsub(".*\\.","",my_ress[r])),
                                           tophits = 25,
                                           isFreq = isFreq, 
                                           my_width = my_width)
                        }
                }
                
                dev.off()
                
        }  
        
}




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



