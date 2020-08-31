


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


source("../functions/functions_DESeq.R")


cbPalette <- c("#999999", "#0072B2", "#CC79A7", "#009E73", "#E69F00", "#D55E00", "#56B4E9", "#F0E442")




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

my_chrX_genes <- my_genes[  seqnames(my_genes) == "X"]
my_chr2L_genes <- my_genes[  seqnames(my_genes) == "2L"]
my_chrY_genes <- my_genes[  seqnames(my_genes) == "Y"]



my_gtf <- import.gff(gtffile)

my_3UTR <- my_gtf[my_gtf$type == "3UTR",]
#my_3UTR_length <- aggregate(width(my_3UTR), by = list(my_3UTR$gene_symbol), FUN = sum)
my_3UTR_length <- aggregate(width(my_3UTR), by = list(my_3UTR$gene_symbol), FUN = max)

my_long3UTR_genes <- my_3UTR_length[,1][  my_3UTR_length[,2] > quantile(my_3UTR_length[,2], 0.95)  ]
names(my_long3UTR_genes) <- mapIds(org.Dm.eg.db, my_long3UTR_genes, "FLYBASE", keytype="SYMBOL", multiVals="first")
my_long3UTR_genes <- my_long3UTR_genes[!(is.na(names(my_long3UTR_genes)))]

class.long3UTR <- names(my_long3UTR_genes)

my_short3UTR_genes <- my_3UTR_length[,1][  my_3UTR_length[,2] < quantile(my_3UTR_length[,2], 0.05)  ]
names(my_short3UTR_genes) <- mapIds(org.Dm.eg.db, my_short3UTR_genes, "FLYBASE", keytype="SYMBOL", multiVals="first")
my_short3UTR_genes <- my_short3UTR_genes[!(is.na(names(my_short3UTR_genes)))]


class.short3UTR <- names(my_short3UTR_genes)



########################################################  classes  ######################################################## 

editedRNAs_symbols <- read.delim("../data_folder/nsmb.2675-S2.txt", stringsAsFactor = FALSE)
editedRNAs_symbols <- editedRNAs_symbols[,1]
editedRNAs_symbols <- gsub("-R.*", "", editedRNAs_symbols[-1:-2])

class.editedRNAs <-  unique(my_gtf$gene_id[ my_gtf$gene_symbol %in% editedRNAs_symbols])

##################################

class.mito <- my_genes[seqnames(my_genes) == "mitochondrion_genome"]$gene_id

#############################

snoRNA_snopy <- read.table("../data_folder/snoRNA_snopy.txt", header = T, sep = "\t", stringsAsFactors = FALSE)

class.snoRNA_HACA <- mapIds(org.Dm.eg.db, keys = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "H/ACA"], keytype = "SYMBOL", column = "FLYBASE")

class.snoRNA_CD <-  mapIds(org.Dm.eg.db, keys = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "C/D"], keytype = "SYMBOL", column = "FLYBASE")


#############################

class.ribo <- my_genes[grep("^RpL|^RpS", my_genes$symbol)]$gene_id

############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################  SampleTable  ######################################################



SampleTable_Global <- read.xls("../SampleTable_Global.xlsx", sheet = 1, header = TRUE)



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################    results    ######################################################


my_ress <- list.files(path = "../analysis/",pattern = "^res")

i=1

for(i in seq_along(my_ress)){
      
      my_name <- gsub(".txt","",my_ress[i])
      
      my_res <- read.delim(file.path("../analysis/", my_ress[i]), header = T, row.names = 1)
      
      assign(my_name, my_res)
      
}



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

###################################################      Fisher tests     ###################################################


my_ress <- ls(pattern = "^res")

my_classs <- ls(pattern = "^class")


test.table <- data.frame(Experiment = NA,
                         Comparison = NA,
                         RNA.Class = NA,
                         Sign.inClass = NA,
                         NS.inClass = NA,
                         Sign.notClass = NA,
                         NS.notClass = NA,
                         P.Value = NA)
r=15
cl=2

for(r in seq_along(my_ress)){
        
        for(cl in seq_along(my_classs)){
                
                my_res <- get(my_ress[r])
                
                my_class <- get(my_classs[cl])
                
                Experiment <- gsub("\\..*","",gsub("res.","",my_ress[r]))
                Comparison <- gsub(".*\\.","",my_ress[r])
                RNA.Class <- gsub("class.","", my_classs[cl])
                
                Sign.inClass <-  nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj < 0.01  & my_res$log2FoldChange > -100,])
                NS.inClass <-    nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj >= 0.01 & my_res$log2FoldChange > -100,])
                Sign.notClass <- nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj < 0.01  & my_res$log2FoldChange > -100,])
                NS.notClass <-   nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj >= 0.01 & my_res$log2FoldChange > -100,])
                
                P.Value <- fisher.test(matrix(c(Sign.inClass,Sign.notClass,NS.inClass,NS.notClass),nrow = 2, byrow = T))$p.value
                P.Value <- format(P.Value, digits = 2, scientific = T)
                
                test.table <- rbind(test.table, 
                                    c(Experiment,Comparison,RNA.Class,Sign.inClass,NS.inClass,Sign.notClass,NS.notClass,P.Value))
        }
        
        
        
}


test.table.all <- test.table[complete.cases(test.table),]

write.table(test.table.all, file = "enrichment_table_all.txt", sep="\t", quote = F, row.names = F)




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

###############################################      Fisher tests version 2     #############################################


my_ress <- ls(pattern = "^res")

my_classs <- ls(pattern = "^class")


test.table <- data.frame(Experiment = NA,
                         Comparison = NA,
                         RNA.Class = NA,
                         Sign.inClass = NA,
                         NS.inClass = NA,
                         Sign.notClass = NA,
                         NS.notClass = NA,
                         P.Value = NA)
r=15
cl=2

for(r in seq_along(my_ress)){
        
        for(cl in seq_along(my_classs)){
                
                my_res <- get(my_ress[r])
                
                my_class <- get(my_classs[cl])
                
                Experiment <- gsub("\\..*","",gsub("res.","",my_ress[r]))
                Comparison <- gsub(".*\\.","",my_ress[r])
                RNA.Class <- gsub("class.","", my_classs[cl])
                
                Sign.inClass <-  nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj < 0.01  & my_res$log2FoldChange > 0,])
                NS.inClass <-    nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj >= 0.01 & my_res$log2FoldChange > 0,])
                Sign.notClass <- nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj < 0.01  & my_res$log2FoldChange > 0,])
                NS.notClass <-   nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj >= 0.01 & my_res$log2FoldChange > 0,])
                
                P.Value <- fisher.test(matrix(c(Sign.inClass,Sign.notClass,NS.inClass,NS.notClass),nrow = 2, byrow = T))$p.value
                P.Value <- format(P.Value, digits = 2, scientific = T)
                
                test.table <- rbind(test.table, 
                                    c(Experiment,Comparison,RNA.Class,Sign.inClass,NS.inClass,Sign.notClass,NS.notClass,P.Value))
        }
        
        
        
}


test.table.up <- test.table[complete.cases(test.table),]

write.table(test.table.up, file = "enrichment_table_up.txt", sep="\t", quote = F, row.names = F)





############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

###############################################      Fisher tests version 3     #############################################


my_ress <- ls(pattern = "^res")

my_classs <- ls(pattern = "^class")


test.table <- data.frame(Experiment = NA,
                         Comparison = NA,
                         RNA.Class = NA,
                         Sign.inClass = NA,
                         NS.inClass = NA,
                         Sign.notClass = NA,
                         NS.notClass = NA,
                         P.Value = NA)
r=15
cl=2

for(r in seq_along(my_ress)){
        
        for(cl in seq_along(my_classs)){
                
                my_res <- get(my_ress[r])
                
                my_class <- get(my_classs[cl])
                
                Experiment <- gsub("\\..*","",gsub("res.","",my_ress[r]))
                Comparison <- gsub(".*\\.","",my_ress[r])
                RNA.Class <- gsub("class.","", my_classs[cl])
                
                Sign.inClass <-  nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj < 0.01  & my_res$log2FoldChange < 0,])
                NS.inClass <-    nrow(my_res[(rownames(my_res) %in% my_class)  & my_res$padj >= 0.01 & my_res$log2FoldChange < 0,])
                Sign.notClass <- nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj < 0.01  & my_res$log2FoldChange < 0,])
                NS.notClass <-   nrow(my_res[!(rownames(my_res) %in% my_class) & my_res$padj >= 0.01 & my_res$log2FoldChange < 0,])
                
                P.Value <- fisher.test(matrix(c(Sign.inClass,Sign.notClass,NS.inClass,NS.notClass),nrow = 2, byrow = T))$p.value
                P.Value <- format(P.Value, digits = 2, scientific = T)
                
                test.table <- rbind(test.table, 
                                    c(Experiment,Comparison,RNA.Class,Sign.inClass,NS.inClass,Sign.notClass,NS.notClass,P.Value))
        }
        
        
        
}


test.table.down <- test.table[complete.cases(test.table),]

write.table(test.table.down, file = "enrichment_table_down.txt", sep="\t", quote = F, row.names = F)



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



