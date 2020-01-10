


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

my_short3UTR_genes <- my_3UTR_length[,1][  my_3UTR_length[,2] < quantile(my_3UTR_length[,2], 0.05)  ]
names(my_short3UTR_genes) <- mapIds(org.Dm.eg.db, my_short3UTR_genes, "FLYBASE", keytype="SYMBOL", multiVals="first")
my_short3UTR_genes <- my_short3UTR_genes[!(is.na(names(my_short3UTR_genes)))]





editedRNAs_symbols <- read.xls("../data_folder/nsmb.2675-S2.xlsx", sheet = 4, stringsAsFactor = FALSE)
editedRNAs_symbols <- unique(editedRNAs_symbols[,1])
editedRNAs_symbols <- unique(gsub("-R.*", "", editedRNAs_symbols[-1:-2]))


editedRNAs <-  unique(my_gtf$gene_id[ my_gtf$gene_symbol %in% editedRNAs_symbols])




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

########################################################  SampleTable  ######################################################



SampleTable_Global <- read.xls("../SampleTable_Global.xlsx", sheet = 1, header = TRUE)




########################################################  count data  ########################################################


myCount_Dir <- "../counts/"

myCount_Files <- list.files(path = "../counts/", pattern = ".out.tab")
myCount_File_Path <- file.path(myCount_Dir, myCount_Files)



######################################################## 

### create count matrix ###


my_counts <-  makeCountTable(count_files = myCount_Files, 
                             count_file_path = myCount_File_Path,
                             stranded = TRUE)


my_counts_Genes <- my_counts[grep("FBgn", rownames(my_counts)),]


# remove Y chromosome
my_counts_Genes <- my_counts_Genes[!(rownames(my_counts_Genes) %in% my_genes[seqnames(my_genes) == "Y"]$gene_id ),]


write.csv(my_counts_Genes, "../data_folder/my_counts_genes.csv")



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 






      
      
      
      
      
      
      
      
      
      
      
      
      
      

############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

      
      
############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################       Setup All DDS objects      ###############################################





my_comparisons <- unique(SampleTable_Global$comparison)



for(my_comparison in my_comparisons){
      
      dds <- setupDDS(my_comparison = my_comparison)
      
      
      assign(paste0("dds.", my_comparison) , dds)
      
      
}





############################################       Batch Corrected Counts      ###############################################


my_ddss <- ls(pattern = "^dds\\.")


for(my_dds in my_ddss){
      
      dds <- get(my_dds)
      
      log2_norm_counts <- correctCounts(dds)
      
      #log2_norm_counts <- log2(counts(dds, normalized = TRUE)+1)
      
      
      assign(gsub("dds.","log2_norm_counts.", my_dds)  , log2_norm_counts)
      
}




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



######################################################## 


load("../data_folder/ranges_mat.rda")


ranges1_mat_mean <- calc_mean_ranges_counts(ranges_mat = ranges_mat, 
                                            select_features = c("snoRNA", "snRNA", "tRNA"),
                                            comparison = "MLE_titration_ATP_S2_RNA")



ranges2_mat_mean <- calc_mean_ranges_counts(ranges_mat = ranges_mat, 
                                            select_features = c("exon", "intron", "fiveUTR", "threeUTR"),
                                            comparison = "MLE_titration_ATP_S2_RNA")



ranges3_mat_mean <- calc_mean_ranges_counts(ranges_mat = ranges_mat, 
                                            select_features = c("snoRNA", "snRNA", "tRNA"),
                                            comparison = "eMLE_lowMLE_S2_cells|MLE_IP_4MSL_ATP_S2_RNA")




ranges4_mat_mean <- calc_mean_ranges_counts(ranges_mat = ranges_mat, 
                                            select_features = c("exon", "intron", "fiveUTR", "threeUTR"),
                                            comparison = "eMLE_lowMLE_S2_cells|MLE_IP_4MSL_ATP_S2_RNA")



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 














############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

#######################################################   Figure 1   ########################################################



######################################################## 



pdf("Figure1.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)

par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)





######################################################## 


my_conditions <- colData(dds.MLE_ATP_clone8_RNA)$Sample

my_colors <- c("#8D9093", "#C199CE", "#DA0EA1", "#5026D9", "#90E7F8", "#1971A9")


par(fig = c(0.30,0.47,0.80,1.00),mar=c(2,3,2,0), new = FALSE)

my_favorite_gene <- "roX1"

ylims <- c(-0.5,6)

plotDots(my_data = log2_norm_counts.MLE_ATP_clone8_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("cl.8 RNA")), bty = "n", x.intersp = 0, cex = 0.9)

axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)


par(fig = c(0.45,0.62,0.80,1.00),mar=c(2,3,2,0), new = TRUE)

my_favorite_gene <- "roX2"


plotDots(my_data = log2_norm_counts.MLE_ATP_clone8_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("cl.8 RNA")), bty = "n", x.intersp = 0, cex = 0.9)



######################################################## 



my_conditions <- colData(dds.MLE_titration_ATP_S2_RNA)$Sample

my_conditions <- factor(my_conditions, levels = unique(my_conditions)[c(1,3,2,5,6,4)])


par(fig = c(0.60,0.80,0.80,1.00),mar=c(2,3,2,0), new = TRUE)

my_favorite_gene <- "roX2"


plotDots(my_data = log2_norm_counts.MLE_titration_ATP_S2_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)





######################################################## 


par(fig = c(0.80,1.00,0.80,1.00), mar=c(1,0,1,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)


legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")






######################################################## 
######################################################## 



par(fig = c(0.00,0.25,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2")

my_contrast <- c("Sample", "cl8.MLE_plusATP", "cl8.Input")

summarizeResults(dds = "dds.MLE_ATP_clone8_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 ylab = "log2FC [IP-Input]")

legend("topleft", legend = expression(italic("cl.8 RNA")), bty = "n", x.intersp = 0, cex = 0.9)


######################################################## 






par(fig = c(0.30,0.55,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "High.MLE_plusATP", "Input")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 xlab = "", 
                 ylab = "log2FC [IP-Input]")


legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)


######################################################## 

par(fig = c(0.50,0.75,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "Med.MLE_plusATP", "Input")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 ylab = "")

legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)


######################################################## 


par(fig = c(0.70,0.95,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "Low.MLE_plusATP", "Input")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 xlab = "", ylab = "")


legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)




######################################################## 
######################################################## 





par(fig = c(0.00,0.275,0.40, 0.60),mar=c(3,3,2,2), mgp=c(1.5,0.5,0),new = TRUE)

par(cex.lab=0.6)
plotLog2FC(dds = dds.MLE_titration_ATP_S2_RNA, 
           contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
           contrast2 = c("Sample", "Low.MLE_plusATP", "Input"), 
           favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
           favorite_genes2 = "roX2", 
           selection_name1 = "mito",
           selection_name2 = "roX2", 
           selection_color2 = rgb(0.9,0.6,0,1), 
           selection2_point_size = 1,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-10,10)
)


######################################################## 

par(fig = c(0.200,0.475,0.40, 0.60),mar=c(3,3,2,2), new = TRUE)


par(cex.lab=0.6)
plotLog2FC(dds = dds.MLE_titration_ATP_S2_RNA, 
           contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
           contrast2 = c("Sample", "Low.MLE_plusATP", "Input"), 
           favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
           favorite_genes2 = my_genes[grep("^RpL|^RpS", my_genes$symbol)]$symbol, 
           selection_name1 = "snoRNA",
           selection_name2 = "ribo", 
           selection_color2 = rgb(0.1,0.8,0,1), 
           y_label = FALSE, 
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-10,10)
)



######################################################## 



ranges_colors <- c("#000000", "#F2F2F2", "#737373")


par(fig = c(0.425,0.725,0.40,0.60), mar=c(3,5,3,1), cex.lab=0.9,  mgp=c(1.5,0.5,0), new=TRUE)

barplot(height = t(as.matrix(ranges1_mat_mean[c(3,4,5,1,2),])), 
        horiz = TRUE,
        col = ranges_colors, las=1, xlab="Fraction")



######################################################## 


par(fig = c(0.450,0.750,0.40,0.60), mar=c(3,1,1,1), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")


legend("top", legend = colnames(ranges1_mat_mean), fill = ranges_colors, border = "black",
       horiz = TRUE, cex=0.8, pt.cex = 1.25, text.width = c(0.275, 0.250, 0.225), x.intersp = 0.15, bty = "n")



######################################################## 


par(fig = c(0.725,1.00,0.40, 0.60),mar=c(3,3,2,2),new = TRUE)


snoRNA_snopy <- read.table("../data_folder/snoRNA_snopy.txt", header = T, sep = "\t", stringsAsFactors = FALSE)


par(cex.lab=0.6)
plotLog2FC(dds = dds.MLE_titration_ATP_S2_RNA, 
           contrast1 = c("Sample", "High.MLE_plusATP", "High.MLE_minusATP"), 
           contrast2 = c("Sample", "Low.MLE_plusATP", "Low.MLE_minusATP"), 
           favorite_genes1 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "H/ACA"], 
           favorite_genes2 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "C/D"], 
           selection_name1 = "H/ACA",
           selection_name2 = "C/D", 
           selection_color1 = rgb(0.7,0,0.7,1), 
           selection_color2 = rgb(0,0.7,0.7,1), 
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-10,10)
)



######################################################## 
######################################################## 




par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



my_conditions <- colData(dds.MLE_titration_ATP_S2_RNA)$Sample

my_conditions <- factor(my_conditions, levels = unique(my_conditions)[c(1,3,2,5,6,4)])


par(fig = c(0.57,0.77,0.20,0.40),mar=c(2,3,2,0), new = TRUE)

my_favorite_gene <- "ND1"


plotDots(my_data = log2_norm_counts.MLE_titration_ATP_S2_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)





######################################################## 

par(fig = c(0.75,0.95,0.20,0.40),mar=c(2,3,2,0), new = TRUE)

my_favorite_gene <- "ND3"


plotDots(my_data = log2_norm_counts.MLE_titration_ATP_S2_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)





######################################################## 






######################################################## 


dev.off()




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


































############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

##############################################       Figure 1 related     ###################################################



######################################################## 



pdf("Figure1related.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)

par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)




######################################################## 



ranges_colors <- c("#000000", "#F2F2F2", "#737373")

ranges5_mat_mean <- calc_mean_ranges_counts(ranges_mat = ranges_mat, 
                                            select_features = c("snoRNA", "snRNA", "tRNA"),
                                            comparison = "^MLE_ATP_clone8_RNA$")


par(fig = c(0.400,0.700,0.40,0.60), mar=c(3,5,3,1), cex.lab=0.9,  mgp=c(1.5,0.5,0), new=TRUE)

barplot(height = t(as.matrix(ranges5_mat_mean)), 
        horiz = TRUE,
        col = ranges_colors, las=1, xlab="Fraction")


######################################################## 


# par(fig = c(0.425,0.725,0.40,0.60), mar=c(3,1,1,1), new = TRUE)
# 
# plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
# 
# 
# legend("top", legend = colnames(ranges5_mat_mean), fill = ranges_colors, border = "black",
#        horiz = TRUE, cex=0.8, pt.cex = 1.25, text.width = c(0.275, 0.250, 0.225), x.intersp = 0.15, bty = "n")
# 


######################################################## 



ranges_colors <- c("#000000", "#F2F2F2", "#737373")

ranges6_mat_mean <- calc_mean_ranges_counts(ranges_mat = ranges_mat, 
                                            select_features = c("snoRNA", "snRNA", "tRNA"),
                                            comparison = "MLE_ATP_OrR_Male_RNA|MLE_ATP_OrR_Female_RNA")


par(fig = c(0.675,0.975,0.40,0.60), mar=c(3,5,3,1), cex.lab=0.9,  mgp=c(1.5,0.5,0), new=TRUE)

barplot(height = t(as.matrix(ranges6_mat_mean)[c(4,2,3,1),]), 
        horiz = TRUE,
        col = ranges_colors, las=1, xlab="Fraction")


######################################################## 


par(fig = c(0.675,0.975,0.40,0.60), mar=c(3,1,1,1), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")


legend("top", legend = colnames(ranges6_mat_mean), fill = ranges_colors, border = "black",
       horiz = TRUE, cex=0.8, pt.cex = 1.25, text.width = c(0.275, 0.250, 0.225), x.intersp = 0.15, bty = "n")



######################################################## 

dev.off()


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 














############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

#######################################################   Figure S1   ########################################################


######################################################## 




pdf("FigureS1.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)

par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



######################################################## 

df <- getNumberOfDetectedReads(sample_type = "Input")


par(fig = c(0.325,0.575,0.80,1.00),mar=c(4,4,2,0))

df <- df[df$source != "Kc",]

set.seed(2)

stripchart(detected_gene_counts ~ source, data = df, 
           col = "black", cex =  0.7,
           ylim = c(0,6000), xlim = c(0.5,5.5),
           vertical = TRUE, jitter = 0.2, method = "jitter", pch=16, 
           ylab = "", xlab = "", las=2)

mtext(text = "detected genes \n per million reads", side = 2, line = 2.5, cex = 0.9)
mtext(text = "Inputs", side = 1, line = 2, cex = 0.9)


######################################################## 


par(fig = c(0.55,0.80,0.80,1.00),mar=c(4,4,2,0), new=TRUE)



plottingPCA(my_data = log2_norm_counts.MLE_titration_ATP_S2_RNA,
            color_palette = my_colors,
            conditions = my_conditions,
            quantiles = c(0,1),
            point_size = 0.7,
            show_labels = FALSE)


######################################################## 


par(fig = c(0.80,1.00,0.80,1.00), mar=c(1,0,1,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)

legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")



######################################################## 



par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)


######################################################## 

par(fig = c(0.50,0.75,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2", "roX1", "para")

my_contrast <- c("Sample", "High.MLE.OrR.M", "OrR.M")

summarizeResults(dds = "dds.MLE_ATP_OrR_Male_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 5*10^6), 
                 plotMA = TRUE, 
                 ylab = "log2FC [IP-Input]")

legend("topleft", legend = expression(italic("Male RNA")), bty = "n", x.intersp = 0, cex = 0.9)

######################################################## 

par(fig = c(0.70,0.95,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2", "roX1", "para")

my_contrast <- c("Sample", "High.MLE.OrR.F", "OrR.F")

summarizeResults(dds = "dds.MLE_ATP_OrR_Female_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 5*10^6), 
                 plotMA = TRUE, 
                 ylab = "", xlab="")

legend("topleft", legend = expression(italic("Female RNA")), bty = "n", x.intersp = 0, cex = 0.9)




######################################################## 
######################################################## 


par(fig = c(0.00,0.25,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "High.MLE_plusATP", "Low.MLE_plusATP")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 ylab = "log2FC",
                 my_main = "High - Low [+ATP]")




######################################################## 

par(fig = c(0.20,0.45,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "High.MLE_minusATP", "Low.MLE_minusATP")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 xlab = "", ylab = "", 
                 my_main = "High - Low [-ATP]")


######################################################## 


par(fig = c(0.50,0.75,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "High.MLE_plusATP", "High.MLE_minusATP")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 ylab = "log2FC",
                 my_main = "High MLE +/- ATP")




######################################################## 

par(fig = c(0.70,0.95,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX2")

my_contrast <- c("Sample", "Low.MLE_plusATP", "Low.MLE_minusATP")

summarizeResults(dds = "dds.MLE_titration_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 10^6), 
                 plotMA = TRUE, 
                 xlab = "", ylab = "", 
                 my_main = "Low MLE +/- ATP")





######################################################## 
######################################################## 



my_files <- list.files(pattern = "res.*.txt")

i=1

for(i in seq_along(my_files)){
      
      my_name <- gsub(".txt","", my_files[i])
      
      my_res <- read.table(my_files[i], sep="\t", row.names = 1, header = T)  
      
      assign(my_name, my_res)
}



######################################################## 


my_mat <- setupMatforHeatmap(res_name = "res.MLE_titration_ATP_S2_RNA.HighMLE+ATP-Input",
                             log2_norm_counts_name = "log2_norm_counts.MLE_titration_ATP_S2_RNA",
                             pval_cutoff = 0.01,
                             n_hits = 100,
                             scaling = TRUE)



par(oma=c(2,2,2,2), mar=c(5,4,3,4), mgp=c(2.5,1,0),
    cex.axis = 1.2, cex.main = 1.25, cex.lab=1.2, pch=19, cex=1)

par(fig = c(0.00,1.00,0.00,0.30),mar=c(4,4,6,4), new=TRUE)


min_value = -4
max_value = 4

col_dist <- dist(t(my_mat))
col_clust <- hclust(col_dist, method = "ward.D2")
col_clust <- reorder.hclust(col_clust, dis = col_dist)

row_clust <- hclust(dist((my_mat)))

plotHeatmap(my_mat = my_mat, 
            my_title = "",
            min_value = min_value, max_value = max_value, 
            my_row_order = row_clust$order, 
            my_col_order = col_clust$order,
            my_color_palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
            useRaster = FALSE)

axis(side = 1, at = seq(0, 1, length.out = ncol(my_mat)), 
     labels = colnames(my_mat)[col_clust$order],
     las=2, cex.axis=0.15, lwd=0, line = -1)


axis(side = 4, at = seq(0, 1, length.out = nrow(my_mat)), 
     labels = rownames(my_mat)[row_clust$order],
     las=2, cex.axis=0.15, lwd=0, line = -1)


par(fig = c(0.075,0.225,0.00,0.30),mar=c(3.5,0,5.5,4.5), new=TRUE, cex=0.01)

plot(as.dendrogram(row_clust),
     xlab="", ylab="", main="", yaxt="n", xaxt="n", horiz=TRUE)

par(cex=0.45) # no idea why to set to this size
par(fig = c(0,1.0,0.170,0.255),mar=c(4,8.75,4,8.75), new=TRUE, cex=0.01)

plot(as.dendrogram(col_clust),
     xlab="", ylab="", main="", yaxt="n", xaxt="n", horiz=FALSE)



par(oma=c(2,2,2,2), mar=c(5,4,3,4), mgp=c(2.5,1,0),
    cex.axis = 1.2, cex.main = 1.25, cex.lab=1.2, pch=19, cex=1)

par(fig = c(0.00,0.05,0.075,0.175), mar=c(0,0.5,0,1.00), mgp=c(2.5,0.5,0), new=TRUE)

plotHeatmapKey(my_mat = t(matrix(seq(min_value,max_value, length.out = 100))), my_title = "",
               min_value = min_value, max_value = max_value, 
               my_color_palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
               useRaster = TRUE, 
               show_yaxis = TRUE, cex.axis = 0.6)

axis(side = 4, at = 0.5, labels = "relative log2 counts", cex.axis=0.7, tck = -0.50)

######################################################## 






dev.off()




############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 











######################################################## 


my_files <- list.files(pattern = "res.*.txt")



i=1


for(i in seq_along(my_files)){
      
      my_name <- gsub(".txt","", my_files[i])
      
      my_res <- read.table(my_files[i], sep="\t", row.names = 1, header = T)  
      
      assign(my_name, my_res)
}


######################################################## 



ranges_genes_freq <- readRDS("../oligo_freqs/ranges_genes_freq.rds")


######################################################## 











############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

#######################################################   Figure 2   ########################################################



######################################################## 



pdf("Figure2.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



######################################################## 




par(fig = c(0.00,0.275,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  

######################################################## 




par(fig = c(0.50,0.775,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 




par(fig = c(0.700,0.975,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  




######################################################## 
######################################################## 



ranges_colors <- c("#D9D9D9", "#FFE78B", "#808080", "#B75858")


par(fig = c(0.05,0.40,0.60,0.80), mar=c(3,5,3,1), cex.lab=0.9,  mgp=c(1.5,0.5,0), new=TRUE)

barplot(height = t(as.matrix(ranges2_mat_mean[c(3,4,5,1,2),])), 
        horiz = TRUE,
        col = ranges_colors, las=1, xlab="Fraction")



######################################################## 


par(fig = c(0.25,0.55,0.60,0.80), mar=c(3,1,3,1), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")


legend("center", legend = colnames(ranges2_mat_mean), fill = ranges_colors, border = "black",
       horiz = FALSE, cex=0.9, pt.cex = 1.25, text.width = 0.05, bty = "n")





######################################################## 


par(fig = c(0.50,0.775,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_long3UTR_genes, 
              favorite_genes2 = "roX2", 
              selection_name1 = "long 3' UTR",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              selection2_point_size = 1,
              text_label2 = FALSE, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.700,0.975,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_long3UTR_genes, 
              favorite_genes2 = "roX2", 
              selection_name1 = "long 3' UTR",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              selection2_point_size = 1,
              text_label2 = FALSE, 
              y_label = FALSE,
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 
######################################################## 


par(fig = c(0.00,0.275,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "High.MLE_minusATP"),
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "H/ACA"], 
              favorite_genes2 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "C/D"], 
              favorite_genes3 = "roX2", 
              selection_name1 = "H/ACA",
              selection_name2 = "C/D", 
              selection_name3 = "roX2",
              selection_color1 = rgb(0.7,0,0.7,1), 
              selection_color2 = rgb(0,0.7,0.7,1), 
              selection_color3 = rgb(0.9,0.6,0,1), 
              y_label = FALSE,
              selection3_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "High.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "H/ACA"], 
              favorite_genes2 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "C/D"], 
              favorite_genes3 = "roX2", 
              selection_name1 = "H/ACA",
              selection_name2 = "C/D", 
              selection_name3 = "roX2",
              selection_color1 = rgb(0.7,0,0.7,1), 
              selection_color2 = rgb(0,0.7,0.7,1), 
              selection_color3 = rgb(0.9,0.6,0,1), 
              y_label = FALSE,
              selection3_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  






######################################################## 


par(fig = c(0.50,0.775,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Low.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "H/ACA"], 
              favorite_genes2 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "C/D"], 
              favorite_genes3 = "roX2", 
              selection_name1 = "H/ACA",
              selection_name2 = "C/D", 
              selection_name3 = "roX2",
              selection_color1 = rgb(0.7,0,0.7,1), 
              selection_color2 = rgb(0,0.7,0.7,1), 
              selection_color3 = rgb(0.9,0.6,0,1), 
              y_label = FALSE,
              selection3_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.700,0.975,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Low.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "H/ACA"], 
              favorite_genes2 = snoRNA_snopy$snoRNA.name[snoRNA_snopy$Box == "C/D"], 
              favorite_genes3 = "roX2", 
              selection_name1 = "H/ACA",
              selection_name2 = "C/D", 
              selection_name3 = "roX2",
              selection_color1 = rgb(0.7,0,0.7,1), 
              selection_color2 = rgb(0,0.7,0.7,1),
              selection_color3 = rgb(0.9,0.6,0,1), 
              y_label = FALSE,
              selection3_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  




######################################################## 

dev.off()

######################################################## 









############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 































############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

#################################################    Figure 2 related    ####################################################



######################################################## 



pdf("Figure2related.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



######################################################## 




par(fig = c(0.00,0.275,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_ATP_clone8_RNA, 
              contrast1 = c("Sample", "cl8.MLE_plusATP", "cl8.Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_ATP_clone8_RNA, 
              contrast1 = c("Sample", "cl8.MLE_plusATP", "cl8.Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 
######################################################## 




par(fig = c(0.00,0.275,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_ATP_OrR_Male_RNA, 
              contrast1 = c("Sample", "High.MLE.OrR.M", "OrR.M"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_ATP_OrR_Male_RNA, 
              contrast1 = c("Sample", "High.MLE.OrR.M", "OrR.M"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  




######################################################## 
######################################################## 




par(fig = c(0.00,0.275,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_ATP_OrR_Female_RNA, 
              contrast1 = c("Sample", "High.MLE.OrR.F", "OrR.F"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.MLE_ATP_OrR_Female_RNA, 
              contrast1 = c("Sample", "High.MLE.OrR.F", "OrR.F"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 
######################################################## 



par(fig = c(0.50,0.775,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)


my_res_name <- "res.MLE_titration_ATP_S2_RNA.HighMLE+ATP-Input"
my_res <- get(my_res_name)

my_subset <- my_res$symbol %in% editedRNAs_symbols

plot(density(my_res$log2FoldChange), 
     xlim = c(-4,4), ylim = c(0,0.8), lwd = 1.25,
     main = "", xlab="log2FC [IP-Input]")

lines(density(my_res$log2FoldChange[my_subset]), col="red", lwd=1.25)


my_res_name <- "res.MLE_ATP_clone8_RNA.cl8MLE+ATP-cl8Input"
my_res <- get(my_res_name)

my_subset <- my_res$symbol %in% editedRNAs_symbols

lines(density(my_res$log2FoldChange), col="black", lwd=1.25, lty=2)
lines(density(my_res$log2FoldChange[my_subset]), col="red", lwd=1.25, lty=2)

abline(v=0, col="grey32")

par(fig = c(0.700,0.975,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

plot.new()
legend("topleft", legend = c("all (S2)", "edited (S2)","all (cl.8)", "edited (cl.8)"), 
       col = c("black","red", "black", "red"), lwd=1.5, lty=c(1,1,2,2), cex=0.7)



######################################################## 


par(fig = c(0.50,0.775,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)


my_res_name <- "res.MLE_ATP_OrR_Male_RNA.HighMLEOrRM-OrRM"
my_res <- get(my_res_name)

my_subset <- my_res$symbol %in% editedRNAs_symbols

plot(density(my_res$log2FoldChange), 
     xlim = c(-4,4), ylim = c(0,0.8), lwd = 1.25,
     main = "", xlab="log2FC [IP-Input]")

lines(density(my_res$log2FoldChange[my_subset]), col="red", lwd=1.25)


my_res_name <- "res.MLE_ATP_OrR_Female_RNA.HighMLEOrRF-OrRF"
my_res <- get(my_res_name)

my_subset <- my_res$symbol %in% editedRNAs_symbols

lines(density(my_res$log2FoldChange), col="black", lwd=1.25, lty=2)
lines(density(my_res$log2FoldChange[my_subset]), col="red", lwd=1.25, lty=2)

abline(v=0, col="grey32")

par(fig = c(0.700,0.975,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

plot.new()
legend("topleft", legend = c("all (Male)", "edited (Male)","all (Female)", "edited (Female)"), 
       col = c("black","red", "black", "red"), lwd=1.5, lty=c(1,1,2,2), cex=0.7)








######################################################## 

dev.off()

######################################################## 









############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



























############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

######################################################   Figure S2   ########################################################



######################################################## 



pdf("FigureS2.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



######################################################## 
######################################################## 



par(fig = c(0.00,0.275,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.7)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "AAAA",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "GGGG",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  

######################################################## 




par(fig = c(0.50,0.775,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "AAAA",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 




par(fig = c(0.700,0.975,0.80, 1.00),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "GGGG",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  






######################################################## 
######################################################## 









par(fig = c(0.00,0.275,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  

######################################################## 




par(fig = c(0.50,0.775,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 




par(fig = c(0.700,0.975,0.60, 0.80),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "mito",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  







######################################################## 
######################################################## 




par(fig = c(0.00,0.275,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "AAAA",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "GGGG",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  

######################################################## 




par(fig = c(0.50,0.775,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "AAAA",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 




par(fig = c(0.700,0.975,0.40, 0.60),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "GGGG",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  






######################################################## 
######################################################## 









par(fig = c(0.00,0.275,0.20, 0.40),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.20, 0.40),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  

######################################################## 




par(fig = c(0.50,0.775,0.20, 0.40),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 




par(fig = c(0.700,0.975,0.20, 0.40),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_minusATP", "Input"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  




######################################################## 
######################################################## 








par(fig = c(0.00,0.275,0.00, 0.20),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "High.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  



######################################################## 




par(fig = c(0.200,0.475,0.00, 0.20),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "High.MLE_plusATP", "High.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  

######################################################## 




par(fig = c(0.50,0.775,0.00, 0.20),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Low.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  





######################################################## 




par(fig = c(0.700,0.975,0.00, 0.20),mar=c(3,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.6)

plotOligoFreq(dds = dds.MLE_titration_ATP_S2_RNA, 
              contrast1 = c("Sample", "Low.MLE_plusATP", "Low.MLE_minusATP"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "CCCC",
              favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
              favorite_genes2 = "roX2", 
              selection_name1 = "snoRNA",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = FALSE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  




######################################################## 
######################################################## 







######################################################## 

dev.off()

######################################################## 









############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



























############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

######################################################   Figure 3   ########################################################



######################################################## 



pdf("Figure3.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



######################################################## 
######################################################## 


par(fig = c(0.00,0.25,0.80,1.00),mar=c(2,4,2,0), new=FALSE)


my_favorite_gene <- c("roX1","roX2", "RpS29")

my_contrast <- c("Sample", "eMLE", "eInput")

summarizeResults(dds = "dds.eMLE_lowMLE_S2_cells", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "log2FC [IP-Input]")


######################################################## 

par(fig = c(0.25,0.425,0.80,1.00),mar=c(2,3,2,0), new = TRUE)

my_conditions <- colData(dds.eMLE_lowMLE_S2_cells)$Sample

my_conditions <- factor(my_conditions, levels = unique(my_conditions)[c(1,3,2,4)])

my_colors <- c("#8D9093", "#262626", "#1971A9", "#59A662")


my_favorite_gene <- "roX2"

ylims <- c(-1,9)

plotDots(my_data = log2_norm_counts.eMLE_lowMLE_S2_cells, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

#legend("topleft", legend = expression(italic("eMLE")), bty = "n", x.intersp = 0, cex = 0.9)


axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)


######################################################## 

par(fig = c(0.400,0.575,0.80,1.00),mar=c(2,3,2,0), new = TRUE)

my_conditions <- colData(dds.eMLE_lowMLE_S2_cells)$Sample

my_conditions <- factor(my_conditions, levels = unique(my_conditions)[c(1,3,2,4)])

my_colors <- c("#8D9093", "#262626", "#1971A9", "#59A662")


my_favorite_gene <- "RpS29"

ylims <- c(-1,9)

plotDots(my_data = log2_norm_counts.eMLE_lowMLE_S2_cells, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

#legend("topleft", legend = expression(italic("eMLE")), bty = "n", x.intersp = 0, cex = 0.9)



######################################################## 


par(fig = c(0.575,0.750,0.80,1.00), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)


legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")




######################################################## 




par(fig = c(0.725,1.000,0.80, 1.00),mar=c(2,4,2,1), mgp=c(2,0.5,0),new = TRUE)

par(cex.lab=0.8)

plotOligoFreq(dds = dds.eMLE_lowMLE_S2_cells, 
              contrast1 = c("Sample", "eMLE", "eInput"), 
              freq_ranges = ranges_genes_freq, 
              oligo = "TTTT",
              favorite_genes1 = my_long3UTR_genes, 
              favorite_genes2 = "roX2", 
              selection_name1 = "long 3´UTR",
              selection_name2 = "roX2",
              selection_color2 = rgb(0.9,0.6,0,1), 
              text_label2 = FALSE, 
              y_label = TRUE,
              selection2_point_size = 1, 
              pval_cutoff = 0.01, 
              lfc_cutoff = 0, 
              xlims = c(-10,10))  




######################################################## 
######################################################## 

par(mgp=c(1.5,0.5,0))

par(fig = c(0.00,0.25,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2", "RpS29")

my_contrast <- c("Sample", "Low.MLE", "Input")

summarizeResults(dds = "dds.MLE_IP_4MSL_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, xlab = "",
                 ylab = "log2FC [IP-Input]")


legend("top", legend = expression(italic("S2 RNA, anti-MLE IP")), bty = "n", x.intersp = 0, cex = 0.75)


######################################################## 


par(fig = c(0.195,0.445,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2", "RpS29")

my_contrast <- c("Sample", "Low.MLE_plus4MSL", "Input")

summarizeResults(dds = "dds.MLE_IP_4MSL_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "")


legend("top", legend = expression(italic("S2 RNA, anti-MLE IP")), bty = "n", x.intersp = 0, cex = 0.75)


######################################################## 


par(fig = c(0.410,0.660,0.60,0.80),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2", "RpS29", "rudhira")

my_contrast <- c("Sample", "Low.MLE_plus4MSL", "Low.MLE")

summarizeResults(dds = "dds.MLE_IP_4MSL_ATP_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "log2FC", xlab = "",
                 my_main = "Low MLE +/-DCC")

legend("top", legend = expression(italic("S2 RNA, anti-MLE IP")), bty = "n", x.intersp = 0, cex = 0.75)


######################################################## 



my_conditions <- colData(dds.MLE_IP_4MSL_ATP_S2_RNA)$Sample

my_conditions <- relevel(my_conditions, ref = "Input")

my_colors <- c("#8D9093", "#C199CE", "#DA0EA1", "#5026D9", "#90E7F8", "#1971A9")


par(fig = c(0.650,0.825,0.60,0.80),mar=c(2,3,2,0), new = TRUE)

my_favorite_gene <- "roX2"

ylims <- c(-1,6.5)

plotDots(my_data = log2_norm_counts.MLE_IP_4MSL_ATP_S2_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

#legend("topleft", legend = expression(italic("eMLE")), bty = "n", x.intersp = 0, cex = 0.9)


axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)




par(fig = c(0.825,1.00,0.60,0.80), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)


legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")




######################################################## 
######################################################## 

par(fig = c(0.650,0.825,0.40,0.60),mar=c(2,3,2,0), new = TRUE)

my_favorite_gene <- "RpS29"

ylims <- c(-1,6.5)

plotDots(my_data = log2_norm_counts.MLE_IP_4MSL_ATP_S2_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

#legend("topleft", legend = expression(italic("eMLE")), bty = "n", x.intersp = 0, cex = 0.9)


axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)




par(fig = c(0.825,1.00,0.60,0.80), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)


legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")




######################################################## 








######################################################## 

dev.off()

######################################################## 









############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 
























############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

######################################################   Figure S3   ########################################################



######################################################## 



pdf("FigureS3.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)





######################################################## 


par(fig = c(0.30,0.55,0.80,1.00),mar=c(4,4,2,0), new=TRUE)

my_conditions <- colData(dds.eMLE_lowMLE_S2_cells)$Sample

my_conditions <- factor(my_conditions, levels = unique(my_conditions)[c(3,1,4,2)])

my_colors <- c("#8D9093", "#262626", "#59A662", "#1971A9")

plottingPCA(my_data = log2_norm_counts.eMLE_lowMLE_S2_cells,
            color_palette = my_colors,
            conditions = my_conditions,
            quantiles = c(0,1),
            point_size = 0.7,
            show_labels = FALSE, 
            my_limits = c(-100,150))





######################################################## 

par(fig = c(0.55,0.775,0.80,1.00), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)


legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")




######################################################## 
######################################################## 



ranges_colors <- c("#000000", "#F2F2F2", "#737373")


par(fig = c(0.05,0.40,0.60,0.80), mar=c(5,5,1,1), cex.lab=0.9,  mgp=c(1.5,0.5,0), new=TRUE)

barplot(height = t(as.matrix(ranges3_mat_mean)[c(7,6,1,2,8,5),]), 
        horiz = TRUE,
        col = ranges_colors, las=1, xlab="Fraction")




######################################################## 


par(fig = c(0.25,0.55,0.60,0.80), mar=c(5,1,1,1), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")


legend("center", legend = colnames(ranges3_mat_mean), fill = ranges_colors, border = "black",
       horiz = FALSE, cex=0.9, pt.cex = 1.25, text.width = 0.05, bty = "n")





######################################################## 




ranges_colors <- c("#D9D9D9", "#FFE78B", "#808080", "#B75858")


par(fig = c(0.50,0.85,0.60,0.80), mar=c(5,5,1,1), cex.lab=0.9,  mgp=c(1.5,0.5,0), new=TRUE)

barplot(height = t(as.matrix(ranges4_mat_mean)[c(7,6,1,2,8,5),]), 
        horiz = TRUE,
        col = ranges_colors, las=1, xlab="Fraction")




######################################################## 




par(fig = c(0.75,0.95,0.60,0.80), mar=c(5,1,1,1), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")


legend("center", legend = colnames(ranges4_mat_mean), fill = ranges_colors, border = "black",
       horiz = FALSE, cex=0.9, pt.cex = 1.25, text.width = 0.05, bty = "n")



######################################################## 
######################################################## 



par(fig = c(0.00,0.25,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2")

my_contrast <- c("Sample", "FLAG.IP", "Input")

summarizeResults(dds = "dds.MLE_IP_FLAG_IP_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "log2FC [IP-Input]")


######################################################## 



par(fig = c(0.20,0.45,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2")

my_contrast <- c("Sample", "MLE.IP", "Input")

summarizeResults(dds = "dds.MLE_IP_FLAG_IP_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "")



######################################################## 



par(fig = c(0.45,0.70,0.40,0.60),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2")

my_contrast <- c("Sample", "FLAG.IP", "MLE.IP")

summarizeResults(dds = "dds.MLE_IP_FLAG_IP_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "log2FC", 
                 my_main = "FLAG - MLE")


######################################################## 


par(fig = c(0.70,0.85,0.40,0.60),mar=c(2,3,2,0), new = TRUE)

my_conditions <- colData(dds.MLE_IP_FLAG_IP_RNA)$Sample

my_conditions <- factor(my_conditions, levels = unique(my_conditions))

my_colors <- c("#8D9093", "#C199CE", "#DA0EA1", "#5026D9", "#90E7F8", "#1971A9")


my_favorite_gene <- "roX2"

ylims <- c(-1,9)

plotDots(my_data = log2_norm_counts.MLE_IP_FLAG_IP_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

#legend("topleft", legend = expression(italic("eMLE")), bty = "n", x.intersp = 0, cex = 0.9)


axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)




par(fig = c(0.85,0.975,0.40,0.60), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)


legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")




######################################################## 




######################################################## 

dev.off()

######################################################## 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 














############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

######################################################   Figure 4   ########################################################



######################################################## 



pdf("Figure4.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)



######################################################## 
######################################################## 


par(fig = c(0.00,0.25,0.80,1.00),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2")

my_contrast <- c("Sample", "4MSLplusMLE_plusATP", "4MSL_plusATP")

summarizeResults(dds = "dds.4MSL_MLE_WTalone_S2_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "log2FC \n [DCC+MLE - DCC-MLE]",
                 my_main = "")


legend("topleft", legend = expression(italic("S2 RNA")), bty = "n", x.intersp = 0, cex = 0.9)


######################################################## 


par(fig = c(0.200,0.450,0.80,1.00),mar=c(2,4,2,0), new=TRUE)


my_favorite_gene <- c("roX1","roX2")

my_contrast <- c("Sample", "cl8.4MSLplusMLE_plusATP", "cl8.4MSL_plusATP")

summarizeResults(dds = "dds.4MSL_MLE_ATP_clone8_RNA", 
                 contrast = my_contrast, 
                 favorite_gene = my_favorite_gene,
                 pval_cutoff = 0.01,
                 lfc_cutoff = 0,
                 ylims = c(-10,10),
                 xlims = c(1, 1e7), 
                 plotMA = TRUE, 
                 ylab = "", my_main = "DCC vitRIP +/- MLE")

legend("topleft", legend = expression(italic("cl.8 RNA")), bty = "n", x.intersp = 0, cex = 0.9)





######################################################## 
######################################################## 





par(fig = c(0.02,0.25,0.60,0.80),mar=c(2,3,2,0), new = TRUE)

my_conditions <- colData(dds.4MSL_MLE_WTmutants_ATP_S2_RNA)$Sample

my_conditions <- relevel(my_conditions, ref ="Input")

my_colors <- c("#8D9093", "#A25840", "#8F2C6A", "#544288", "#4C8B72", "#AD8900")


my_favorite_gene <- "roX2"

ylims <- c(-3,3)

plotDots(my_data = log2_norm_counts.4MSL_MLE_WTmutants_ATP_S2_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

#legend("topleft", legend = expression(italic("eMLE")), bty = "n", x.intersp = 0, cex = 0.9)


axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)



######################################################## 


par(fig = c(0.25,0.50,0.60,0.80), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)

my_legend <- gsub("4MSL\\+", "DCC \\+", my_legend)
my_legend <- gsub("\\+ATP", "", my_legend)

legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")


######################################################## 










######################################################## 

dev.off()

######################################################## 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

























############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 

######################################################   Figure S4   ########################################################



######################################################## 



pdf("FigureS4.pdf", width = 8.27, height = 11.69, useDingbats = FALSE)


par(oma=c(0,0,0,0), mar=c(5,4,3,4), mgp=c(1.5,0.5,0), bg=NA,
    cex.axis = 0.7, cex.main = 0.9, cex.lab=0.9, pch=19, cex=1)




######################################################## 
######################################################## 



par(fig = c(0.02,0.20,0.60,0.80),mar=c(2,3,2,0), new = TRUE)

my_conditions <- colData(dds.4MSL_MLE_ATP_clone8_RNA)$Sample

my_conditions <- relevel(my_conditions, ref ="cl8.Input")

my_colors <- c("#8D9093", "#BF6D40", "#40BF4F")


my_favorite_gene <- "roX2"

ylims <- c(-3,3)

plotDots(my_data = log2_norm_counts.4MSL_MLE_ATP_clone8_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("cl.8 RNA")), bty = "n", x.intersp = 0, cex = 0.9)

axis(2, at = mean(ylims), labels = "relative log2 counts", line = 1, lwd = 0, cex.axis = 0.8)




######################################################## 


par(fig = c(0.18,0.36,0.60,0.80),mar=c(2,3,2,0), new = TRUE)

my_conditions <- colData(dds.4MSL_MLE_ATP_clone8_RNA)$Sample

my_conditions <- relevel(my_conditions, ref ="cl8.Input")

my_colors <- c("#8D9093", "#BF6D40", "#40BF4F")


my_favorite_gene <- "roX1"

ylims <- c(-3,3)

plotDots(my_data = log2_norm_counts.4MSL_MLE_ATP_clone8_RNA, 
         my_favorite_gene = my_favorite_gene, 
         color_palette = my_colors, 
         color_groups = my_conditions, 
         conditions = my_conditions, 
         point_size = 0.7,
         x_label = rep("", length(levels(my_conditions))),
         ylims = ylims,
         xlims = c(0.5,length(levels(my_conditions))+0.5))

legend("topleft", legend = expression(italic("cl.8 RNA")), bty = "n", x.intersp = 0, cex = 0.9)


######################################################## 


par(fig = c(0.36,0.55,0.60,0.80), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)

my_legend <- gsub("4MSL\\+", "DCC \\+", my_legend)
my_legend <- gsub("\\+ATP", "", my_legend)

legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")


######################################################## 
######################################################## 





my_conditions <- colData(dds.4MSL_MLE_WTmutants_ATP_S2_RNA)$Sample

my_conditions <- relevel(my_conditions, ref ="Input")

my_colors <- c("#8D9093", "#A25840", "#8F2C6A", "#544288", "#4C8B72", "#AD8900")



par(fig = c(0.00,0.25,0.40,0.60),mar=c(3,3,2,0), new=TRUE)



plottingPCA(my_data = log2_norm_counts.4MSL_MLE_WTmutants_ATP_S2_RNA,
            color_palette = my_colors,
            conditions = my_conditions,
            quantiles = c(0,1),
            point_size = 0.7,
            show_labels = FALSE, 
            my_limits = c(-150,100))



######################################################## 


par(fig = c(0.25,0.50,0.40,0.60), mar=c(2,0,2,0), new = TRUE)

plot(0:1,0:1, type="n", main="", xlab="", ylab="", bty="n", xaxt="n", yaxt="n")

my_legend <- levels(my_conditions)
my_legend <- cleanLegend(my_legend)

my_legend <- gsub("4MSL\\+", "DCC \\+", my_legend)
my_legend <- gsub("\\+ATP", "", my_legend)

legend("topleft", legend = my_legend, col = my_colors, 
       pch=19, horiz = FALSE, cex=0.9, pt.cex = 1, text.width = 0.05, bty = "n")


######################################################## 
######################################################## 




par(fig = c(0.00,0.275,0.20, 0.40),mar=c(3,3,2,2), mgp=c(1.5,0.5,0),new = TRUE)

par(cex.lab=0.6)
plotLog2FC(dds = dds.4MSL_MLE_titration_S2_RNA, 
           contrast1 = c("Sample", "4MSLplusMLE_plusATP", "4MSL_plusATP"), 
           contrast2 = c("Sample", "Low.MLE_plusATP", "Input"), 
           favorite_genes1 = my_genes[seqnames(my_genes) == "mitochondrion_genome"]$symbol, 
           favorite_genes2 = "roX2", 
           selection_name1 = "mito",
           selection_name2 = "roX2", 
           selection_color2 = rgb(0.9,0.6,0,1), 
           selection2_point_size = 1,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-10,10)
)



######################################################## 



par(fig = c(0.225,0.500,0.20, 0.40),mar=c(3,3,2,2), mgp=c(1.5,0.5,0),new = TRUE)

par(cex.lab=0.6)
plotLog2FC(dds = dds.4MSL_MLE_titration_S2_RNA, 
           contrast1 = c("Sample", "4MSLplusMLE_plusATP", "4MSL_plusATP"), 
           contrast2 = c("Sample", "Low.MLE_plusATP", "Input"), 
           favorite_genes1 = my_genes[grep("snoRNA", my_genes$symbol)]$symbol, 
           favorite_genes2 = "roX2", 
           selection_name1 = "snoRNA",
           selection_name2 = "roX2", 
           selection_color2 = rgb(0.9,0.6,0,1), 
           selection2_point_size = 1,
           text_label2 = FALSE,
           pval_cutoff = 0.01, 
           lfc_cutoff = 0, 
           lims = c(-10,10)
)





########################################################
######################################################## 




par(fig = c(0.525,1.00,0.10,0.49),mar=c(4,4,6,4), new=TRUE)

min_value = -4
max_value = 4



my_mat <- setupMatforHeatmap(res_name = "",
                             log2_norm_counts_name = "log2_norm_counts.4MSL_MLE_titration_S2_RNA", 
                             my_selection = c("FBgn0019660", my_genes[seqnames(my_genes) == "mitochondrion_genome"]$gene_id),
                             scaling = TRUE)


col_dist <- dist(t(my_mat))
col_clust <- hclust(col_dist, method = "ward.D2")
col_clust <- reorder.hclust(col_clust, dis = col_dist)

row_clust <- hclust(dist( (my_mat)))

plotHeatmap(my_mat = t(my_mat), 
            my_title = "",
            min_value = min_value, max_value = max_value, 
            my_row_order = col_clust$order, 
            my_col_order = row_clust$order,
            my_color_palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
            useRaster = FALSE)

# axis(side = 2, at = seq(0, 1, length.out = ncol(my_mat)), 
#      labels = colnames(my_mat)[col_clust$order],
#      las=2, cex.axis=0.3, lwd=0, line = -0.9)

par(fig = c(0.55,0.70,0.10,0.49),mar=c(3.5,0,5.5,4.5), new=TRUE, cex=0.45)

plot(as.dendrogram(col_clust),
     xlab="", ylab="", main="", yaxt="n", xaxt="n", horiz=TRUE)


par(fig = c(0.525,1.00,0.355,0.445),mar=c(4,8,4,8), new=TRUE, cex=0.01)

plot(as.dendrogram(row_clust),
     xlab="", ylab="", main="", yaxt="n", xaxt="n", horiz=FALSE)



######################################################## 



par(oma=c(2,2,2,2), mar=c(5,4,3,4), mgp=c(2.5,1,0),
    cex.axis = 1.2, cex.main = 1.25, cex.lab=1.2, pch=19, cex=1)



par(fig = c(0.94,0.99,0.20,0.35), mar=c(0,0.75,0,0.75), mgp=c(2.5,0.5,0), new=TRUE)

plotHeatmapKey(my_mat = t(matrix(seq(min_value,max_value, length.out = 100))), my_title = "",
               min_value = min_value, max_value = max_value, 
               my_color_palette = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
               useRaster = TRUE, 
               show_yaxis = TRUE, cex.axis = 0.6)

axis(side = 4, at = 0.5, labels = "relative log2 counts", cex.axis=0.7, tck = -0.5)

######################################################## 



######################################################## 












######################################################## 

dev.off()

######################################################## 



############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 


############################################################################################################################# 
############################################################################################################################# 
############################################################################################################################# 





