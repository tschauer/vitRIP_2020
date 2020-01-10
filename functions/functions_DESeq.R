

library(GenomicFeatures)
library(AnnotationDbi)
library(org.Dm.eg.db)
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
library(lattice)
library(magick)



######################################################## make Count Table ######################################################## 



makeCountTable <- function(count_files, count_file_path, stranded = FALSE){
  
      if(stranded){
        cidx <- 4
      } else {
        cidx <- 2
      }
  
      for(i in seq_along(count_files)){
        
        tmp <- read.table(count_file_path[i])
        
        if(i == 1){
          my_counts <- tmp[,cidx] 
        } else {
          my_counts <- cbind(my_counts, tmp[,cidx])
        }
      }
      
      rownames(my_counts) <- tmp[,1]
      colnames(my_counts) <- gsub("_[G,A,T,C].*","", count_files)
      
      return(my_counts)
}





######################################################## make Kmer Table ######################################################## 



makeKmerTable <- function(kmer_files, kmer_file_path, sample_names){
      
      for(i in seq_along(kmer_files)){
            
            if(i == 1){
                  
                  tmp1 <- read.table(kmer_file_path[i])      
                  colnames(tmp1) <- c("id", sample_names[i])
                  
            } else {
                  
                  tmp2 <- read.table(kmer_file_path[i])
                  colnames(tmp2) <- c("id", sample_names[i])
                  
                  
                  tmp1 <- merge(data.table(tmp1), 
                                data.table(tmp2), by = "id")
            }
      }
      
      tmp1 <- as.data.frame(tmp1)
      rownames(tmp1) <- tmp1$id
      tmp1 <- tmp1[,-1]
      
      return(tmp1)
      
}







######################################################## Summarize Results and MAplot ######################################################## 



# magickPoints <- function(x, y, ...) {
#       w <- convertWidth(unit(1, "npc"), "in", valueOnly=TRUE)
#       h <- convertHeight(unit(1, "npc"), "in", valueOnly=TRUE)
#       cvp <- current.viewport()
#       dev <- dev.cur()
#       raster <- image_graph(width=w*72, height=h*72)
#       pushViewport(viewport(xscale=cvp$xscale, yscale=cvp$yscale))
#       panelDefault(x, y, ...)
#       dev.off()
#       dev.set(dev)
#       grid.raster(raster)
# }


summarizeResults <- function(dds, 
                             contrast, 
                             favorite_gene, 
                             pval_cutoff, lfc_cutoff = 0,
                             ylims = c(-5,5), xlims = c(1, 10^6),
                             plotMA = TRUE,
                             my_main = NULL,
                             xlab = "log10 mean counts", 
                             ylab = "log2 fold change", xaxt="n"
                             ){
  
      comparison_name <- gsub("dds.","", dds)
      
      dds <- get(dds)
      
    res <- results(dds, 
                   contrast = contrast, 
                   lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                   #, alpha = pval_cutoff
                   )
    
    #res <- lfcShrink(dds, contrast = contrast, res = res)
    
    
    res$chr <- mapIds(org.Dm.eg.db, rownames(res), "CHR", keytype="FLYBASE", multiVals="first")
    res$symbol <- mapIds(org.Dm.eg.db, rownames(res), "SYMBOL", keytype="FLYBASE", multiVals="first")
    
    res$symbol <- gsub("lncRNA:","",res$symbol)
    
    res$padj[is.na(res$padj)] <- 1
    res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
    res.sign <- subset(res, padj < pval_cutoff)
    
    if(is.null(my_main)){
          my_main <- paste(contrast[2], collapse = " - ")
          my_main <- cleanLegend(my_main)
    }
    



    if(plotMA){
          plot(res$baseMean, res$log2FoldChange, log="x",
              xlab = "", ylab = ylab, xaxt="n",
              col=rgb(0,0,0,0.1),  ylim=ylims, xlim=xlims, pch=19, cex = 0.25)
          
          axis(side = 1, at = c(1,10^2,10^4,10^6), labels = log10(c(1,10^2,10^4,10^6)))
          axis(1, at = exp(mean(log(xlims))), labels = xlab, line = 1, lwd = 0, cex.axis = 0.9)
          #axis(side = 1, at = c(100), labels = c(2))
          
          
          points(res.sign$baseMean, res.sign$log2FoldChange, col=rgb(0.8,0,0,0.5), pch=19, cex = 0.25)
          
          mtext(text = my_main, side = 3, line = 0.5, adj = 1, font=2)

          for(i in seq_along(favorite_gene)){
                
                my_favorite <- grep(paste0("^", favorite_gene[i],"$"), res$symbol)
                
                if(res$padj[my_favorite] < pval_cutoff){
                      
                      points(res$baseMean[my_favorite],
                             res$log2FoldChange[my_favorite],
                             col=rgb(0.9,0.6,0,1), pch=19)
                      
                      text(res$baseMean[my_favorite],
                           res$log2FoldChange[my_favorite],
                           labels = res$symbol[my_favorite], adj = c(0,-0.5),
                           col=rgb(0.9,0.6,0,1), cex = 0.9)
                } else {
                      
                      points(res$baseMean[my_favorite],
                             res$log2FoldChange[my_favorite],
                             col="grey", pch=19)
                      
                      text(res$baseMean[my_favorite],
                           res$log2FoldChange[my_favorite],
                           labels = res$symbol[my_favorite], adj = c(0,-0.5),
                           col="grey", cex = 0.9)
                }
          }
          
          
          
          abline(h=0, col="grey32")
          
          # legend("bottomright", legend =  c(paste("up", sum(res.sign$log2FoldChange > 0), sep = " - "),
          #                                   paste("down", sum(res.sign$log2FoldChange < 0), sep = " - ")),
          #        bg = NA, cex = 0.6, border = NA, bty = "n")
          
    }
    

    
    my_file_name <- gsub(" ","", cleanLegend(paste(contrast[2:3], collapse = " - ")))
    
    write.table(res[order(res$padj),], file = paste("res", comparison_name, my_file_name, "txt", sep="."), sep = "\t", row.names = T, col.names = NA)
    
}





######################################################## Summarize Kmer Results and MAplot ######################################################## 



summarizeKmerResults <- function(dds, contrast, favorite_gene, pval_cutoff, lfc_cutoff = 0, ylims = c(-5,5)){
      
      res <- results(dds, 
                     contrast = contrast, 
                     lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                     #, alpha = pval_cutoff
      )
      
      #res <- lfcShrink(dds, contrast = contrast, res = res)
      
      
      res$symbol <- rownames(res)
      
      
      res$padj[is.na(res$padj)] <- 1
      res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
      res.sign <- subset(res, padj < pval_cutoff)
      
      
      my_main <- paste(contrast[2:3], collapse = " - ")
      my_main <- gsub("plus", "+", my_main)
      my_main <- gsub("minus", "-", my_main)
      
      
      plot(res$baseMean, res$log2FoldChange, log="x", 
           main = my_main, xlab = "mean expression", ylab = "log2 fold change",
           col=rgb(0,0,0,0.1),  ylim=ylims, pch=19)
      points(res.sign$baseMean, res.sign$log2FoldChange, col=rgb(0.8,0,0,0.5), pch=19, cex = 0.5)
      
      
      for(i in seq_along(favorite_gene)){
            
            my_favorite <- grep(favorite_gene[i], res$symbol)
            
            if(res$padj[my_favorite] < pval_cutoff){
                  text(res$baseMean[my_favorite], 
                       res$log2FoldChange[my_favorite],
                       labels = res$symbol[my_favorite], adj = c(0,-0.5), 
                       col=rgb(0.8,0,0,1))
            } else {
                  text(res$baseMean[my_favorite], 
                       res$log2FoldChange[my_favorite],
                       labels = res$symbol[my_favorite], adj = c(0,-0.5), 
                       col="grey")
            }
      }
      
      
      
      abline(h=0, col="grey32")
      
      legend("bottomright", legend =  c(paste("up", sum(res.sign$log2FoldChange > 0), sep = " - "),
                                        paste("down", sum(res.sign$log2FoldChange < 0), sep = " - ")) )
      
      
      my_main <- gsub(" ","", my_main)
      
      write.csv(res[order(res$padj),], file = paste("res", my_main, "csv", sep="."))
      #write.csv(res.sign[order(res.sign$padj),], file = paste("res.sign", my_main, "csv", sep="."))
      
      
}









######################################################## log2FC vs log2FC ######################################################## 





plotLog2FC <- function(dds, 
                       contrast1,
                       contrast2,
                       favorite_genes1,
                       favorite_genes2,
                       pval_cutoff,
                       lfc_cutoff = 0,
                       lims = c(-5,5),
                       x_label = TRUE,
                       y_label = TRUE,
                       selection_color1 = rgb(0,0,0.9,1),
                       selection_color2 = rgb(0.8,0,0,1),
                       selection_name1 = "",
                       selection_name2 = "",
                       selection1_point_size = 0.5,
                       selection2_point_size = 0.5,
                       text_label1 = FALSE,
                       text_label2 = FALSE){
      
      
      ######
      
      res1 <- results(dds, 
                     contrast = contrast1, 
                     lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                     #, alpha = pval_cutoff
      )
      
      #res1 <- lfcShrink(dds, contrast = contrast1, res1 = res1,)
      
      #res1$chr <- mapIds(org.Dm.eg.db, rownames(res1), "CHR", keytype="FLYBASE", multiVals="first")
      res1$symbol <- mapIds(org.Dm.eg.db, rownames(res1), "SYMBOL", keytype="FLYBASE", multiVals="first")
      
      res1$symbol <- gsub("lncRNA:","",res1$symbol)
      
      res1$padj[is.na(res1$padj)] <- 1
      res1$log2FoldChange[is.na(res1$log2FoldChange)] <- 0
      res1.sign <- subset(res1, padj < pval_cutoff)
      
      if(x_label){
            my_label1 <- paste("log2FC \n", paste(contrast1[2:3], collapse = " - "))
            my_label1 <- cleanLegend(my_label1)
      } else {
            my_label1 <- ""
      }
      
      
      #######
      
      res2 <- results(dds, 
                      contrast = contrast2, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res2 <- lfcShrink(dds, contrast = contrast2, res2 = res2,)
      
      #res2$chr <- mapIds(org.Dm.eg.db, rownames(res2), "CHR", keytype="FLYBASE", multiVals="first")
      res2$symbol <- mapIds(org.Dm.eg.db, rownames(res2), "SYMBOL", keytype="FLYBASE", multiVals="first")
      
      res2$symbol <- gsub("lncRNA:","",res2$symbol)
      
      res2$padj[is.na(res2$padj)] <- 1
      res2$log2FoldChange[is.na(res2$log2FoldChange)] <- 0
      res2.sign <- subset(res2, padj < pval_cutoff)
      
      if(y_label){
            my_label2 <- paste("log2FC \n", paste(contrast2[2:3], collapse = " - "))
            my_label2 <- cleanLegend(my_label2)
      } else {
            my_label2 <- ""
      }
     
      
      #####
      
      stopifnot(
            identical(rownames(res1), rownames(res2))
      )
      
      #####
      
      plot(res1$log2FoldChange,
           res2$log2FoldChange,
           main = "", 
           xlab = my_label1, 
           ylab = my_label2,
           col=rgb(0,0,0,0.1),  ylim=lims, xlim=lims, pch=19, cex = 0.25)
      
      
      
      for(i in seq_along(favorite_genes1)){
            
            my_favorite <- favorite_genes1[i] == res1$symbol
            
            if(length(my_favorite) == 0){
                  next()
            }
            
            stopifnot(
                  identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
            )
            
      
                  points(res1$log2FoldChange[my_favorite], 
                       res2$log2FoldChange[my_favorite],
                       col=selection_color1, pch=19, cex = selection1_point_size)

                  if(text_label1){
                        text(res1$log2FoldChange[my_favorite],
                             res2$log2FoldChange[my_favorite],
                             labels = res1$symbol[my_favorite], adj = c(0,-0.5),
                             col=selection_color1, cex = 0.9)
                  }                  
                  
      }
      
      
      
      
      for(i in seq_along(favorite_genes2)){
            
            my_favorite <- favorite_genes2[i] == res1$symbol
            
            if(length(my_favorite) == 0){
                  next()
            }
            
            stopifnot(
                  identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
            )
            
            
            points(res1$log2FoldChange[my_favorite], 
                   res2$log2FoldChange[my_favorite],
                   col=selection_color2, pch=19, cex = selection2_point_size)
            
            if(text_label2){
                  text(res1$log2FoldChange[my_favorite],
                       res2$log2FoldChange[my_favorite],
                       labels = res1$symbol[my_favorite], adj = c(0,-0.5),
                       col=selection_color2, cex = 0.9)
            }                  
            
      }
      
      abline(h=0, v=0, col="grey32")
      abline(coef = c(0,1), col="grey32", lty=2)
      
            
      legend("topleft", legend = c(selection_name1, selection_name2), bg = "white",
             col = c(selection_color1, selection_color2), pch = 19, cex=0.7)
      
      
      
      

}




######################################################## Kmer log2FC vs log2FC ######################################################## 





plotKmerLog2FC <- function(dds, 
                       contrast1,
                       contrast2,
                       favorite_gene,
                       pval_cutoff,
                       lfc_cutoff = 0,
                       lims = c(-5,5),
                       selection_name = "",
                       text_label = FALSE){
      
      
      ######
      
      res1 <- results(dds, 
                      contrast = contrast1, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res1 <- lfcShrink(dds, contrast = contrast1, res1 = res1,)
      
      res1$symbol <- rownames(res1)
      
      res1$padj[is.na(res1$padj)] <- 1
      res1$log2FoldChange[is.na(res1$log2FoldChange)] <- 0
      res1.sign <- subset(res1, padj < pval_cutoff)
      
      
      my_label1 <- paste("log2FC", paste(contrast1[2:3], collapse = " - "))
      my_label1 <- gsub("plus", "+", my_label1)
      my_label1 <- gsub("minus", "-", my_label1)
      my_label1 <- gsub("\\.", " ", my_label1)
      my_label1 <- gsub("\\_", "", my_label1)
      
      #######
      
      res2 <- results(dds, 
                      contrast = contrast2, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res2 <- lfcShrink(dds, contrast = contrast2, res2 = res2,)
      
      res2$symbol <-rownames(res2)
      
      res2$padj[is.na(res2$padj)] <- 1
      res2$log2FoldChange[is.na(res2$log2FoldChange)] <- 0
      res2.sign <- subset(res2, padj < pval_cutoff)
      
      
      my_label2 <- paste("log2FC", paste(contrast2[2:3], collapse = " - "))
      my_label2 <- gsub("plus", "+", my_label2)
      my_label2 <- gsub("minus", "-", my_label2)
      my_label2 <- gsub("\\.", " ", my_label2)
      my_label2 <- gsub("\\_", "", my_label2)
      
      
      #####
      
      stopifnot(
            identical(rownames(res1), rownames(res2))
      )
      
      #####
      
      plot(res1$log2FoldChange,
           res2$log2FoldChange,
           main = "Comparison", 
           xlab = my_label1, 
           ylab = my_label2,
           col=rgb(0,0,0,0.1),  ylim=lims, xlim=lims, pch=19, cex = 0.5)
      
      
      # for(i in seq_along(res1$symbol)){
      #       
      #       my_favorite <- i
      #       
      #       stopifnot(
      #             identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
      #       )
      #       
      #       if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0.8,0,0,0.5), pch=19)
      #       } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0.8,0.5,0,0.5), pch=19)
      #       } else {
      #             points(res1$log2FoldChange[my_favorite], 
      #                  res2$log2FoldChange[my_favorite],
      #                  col=rgb(0,0,0,0), pch=19)
      #       }
      # }
      # 
      
      
      for(i in seq_along(favorite_gene)){
            
            my_favorite <- favorite_gene[i] == res1$symbol
            
            if(sum(my_favorite) == 0){
                  next()
            }
            
            stopifnot(
                  identical(res1$symbol[my_favorite], res2$symbol[my_favorite])
            )
            
            
            #if(res1$padj[my_favorite] < pval_cutoff & res2$padj[my_favorite] < pval_cutoff){
            
            points(res1$log2FoldChange[my_favorite], 
                   res2$log2FoldChange[my_favorite],
                   col=rgb(0.8,0,0,1), pch=19, cex = 0.75)
            
            if(text_label){
                  text(res1$log2FoldChange[my_favorite],
                       res2$log2FoldChange[my_favorite],
                       labels = res1$symbol[my_favorite], adj = c(0,-0.5),
                       col=rgb(0.8,0,0,1))
            }                  
            
            # } else if(res1$padj[my_favorite] < pval_cutoff | res2$padj[my_favorite] < pval_cutoff){
            #       
            #       points(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            col=rgb(0.8,0.5,0,1), pch =19)
            #       
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col=rgb(0.8,0.5,0,1))
            # } else {
            #       text(res1$log2FoldChange[my_favorite], 
            #            res2$log2FoldChange[my_favorite],
            #            labels = res1$symbol[my_favorite], adj = c(0,1), 
            #            col="grey")
            # }
      }
      
      abline(h=0, v=0, col="grey32")
      
      legend("topleft", legend = selection_name, bg = "white",
             #legend =  c("sign. in both", "sign. in one"),
             col = c(rgb(0.8,0,0,1)), pch = 19)
      
      
      
      
      
}





######################################################## log2FC vs oligo frequency ######################################################## 






plotOligoFreq <- function(dds, 
                       contrast1,
                       freq_ranges,
                       oligo = "TTTT",
                       favorite_genes1 = "",
                       favorite_genes2 = "",
                       favorite_genes3 = "",
                       pval_cutoff,
                       lfc_cutoff = 0,
                       xlims = c(-5,5),
                       ylims = c(0, 0.08),
                       selection_color1 = rgb(0,0,0.9,1),
                       selection_color2 = rgb(0.8,0,0,1),
                       selection_color3 = rgb(0,0.8,0,1),
                       selection_name1 = "",
                       selection_name2 = "",
                       selection_name3 = "",
                       selection1_point_size = 0.5,
                       selection2_point_size = 0.5,
                       selection3_point_size = 0.5,
                       x_label = TRUE,
                       y_label = TRUE,
                       text_label1 = FALSE,
                       text_label2 = FALSE,
                       text_label3 = FALSE,
                       show_cor = FALSE){
      
      
      ######
      
      res1 <- results(dds, 
                      contrast = contrast1, 
                      lfcThreshold = lfc_cutoff, independentFiltering = FALSE
                      #, alpha = pval_cutoff
      )
      
      #res1 <- lfcShrink(dds, contrast = contrast1, res1 = res1,)
      
      res1$symbol <- mapIds(org.Dm.eg.db, rownames(res1), "SYMBOL", keytype="FLYBASE", multiVals="first")
      res1 <- res1[!(is.na(res1$symbol)),]
      
      res1$symbol <- gsub("lncRNA:","",res1$symbol)
      
      res1$padj[is.na(res1$padj)] <- 1
      res1$log2FoldChange[is.na(res1$log2FoldChange)] <- 0
      res1.sign <- subset(res1, padj < pval_cutoff)
      
      #######
      
      if(x_label){
            my_label1 <- paste("log2FC \n", paste(contrast1[2:3], collapse = " - "))
            my_label1 <- cleanLegend(my_label1)
      } else {
            my_label1 <- ""
      }
      
      #######
      
      if(y_label){
            my_label2 <- ifelse(show_cor, paste("Frequency", gsub("TTTT","UUUU",oligo)), "Frequency")
      } else {
            my_label2 <- ""
      }


      
      #######
      
      
      res_merged <- merge(as.data.frame(res1), 
                          as.data.frame(freq_ranges), by.x = "row.names", by.y = "gene_id")
      
      #######
      
      plot(res_merged$log2FoldChange,
           res_merged[oligo][,1],
           main = "",
           xlim = xlims,
           ylim = ylims,
           xlab = "", 
           ylab = "",
           col=rgb(0,0,0,0.1), pch=19, cex = 0.25)
      
      axis(1, at = mean(xlims), labels = my_label1, line = 1.5, lwd = 0)
      
      mtext(text = my_label2, side = 2, line = 2)
      
      if(show_cor){
            title(paste("cor =", round(cor( res_merged$log2FoldChange,res_merged[oligo][,1], method = "pearson"),2)), line = 0.5, cex.main=1 )
            
      } else {
            title(gsub("TTTT","UUUU",oligo), line = 0.5, cex.main=1 )
      }
      
      #####
      
      
      
      
      for(i in seq_along(favorite_genes1)){
            
            my_favorite <- favorite_genes1[i] == res_merged$symbol
            
            if(sum(my_favorite) == 0){
                  next()
            }

            points(res_merged$log2FoldChange[my_favorite],
                   res_merged[oligo][,1][my_favorite],
                   col=selection_color1, pch=19, cex = selection1_point_size)
            
            if(text_label1){
                  text(res_merged$log2FoldChange[my_favorite],
                       res_merged[oligo][,1][my_favorite],
                       labels = res_merged$symbol[my_favorite], adj = c(0,-0.5),
                       col=selection_color1, cex = 0.9)
            }                  
            

      }
      
      for(i in seq_along(favorite_genes2)){
            
            my_favorite <- favorite_genes2[i] == res_merged$symbol
            
            if(sum(my_favorite) == 0){
                  next()
            }
            
            points(res_merged$log2FoldChange[my_favorite],
                   res_merged[oligo][,1][my_favorite],
                   col=selection_color2, pch=19, cex = selection2_point_size)
            
            if(text_label2){
                  text(res_merged$log2FoldChange[my_favorite],
                       res_merged[oligo][,1][my_favorite],
                       labels = res_merged$symbol[my_favorite], adj = c(0,-0.5),
                       col=selection_color2, cex = 0.9)
            }                  
            
            
      }
      
      for(i in seq_along(favorite_genes3)){
          
          my_favorite <- favorite_genes3[i] == res_merged$symbol
          
          if(sum(my_favorite) == 0){
              next()
          }
          
          points(res_merged$log2FoldChange[my_favorite],
                 res_merged[oligo][,1][my_favorite],
                 col=selection_color3, pch=19, cex = selection3_point_size)
          
          if(text_label2){
              text(res_merged$log2FoldChange[my_favorite],
                   res_merged[oligo][,1][my_favorite],
                   labels = res_merged$symbol[my_favorite], adj = c(0,-0.5),
                   col=selection_color3, cex = 0.9)
          }                  
          
          
      }
      
      #####
      
      abline(v=0, col="grey32")
      
      
      my_legend <- c(selection_name1, selection_name2, selection_name3)
      my_legend_color <- c(selection_color1, selection_color2, selection_color3)
      
      legend("topleft", legend = my_legend[my_legend != ""], bg = "white",
             #legend =  c("sign. in both", "sign. in one"),
             col = my_legend_color[my_legend != ""], pch = 19, cex=0.7)
      
  
}







######################################################## Heatmap ######################################################## 






createHeatmap <- function(my_res_name, data_counts, enriched = TRUE, n_genes = 200, favorite_genes = NULL, horizontal = FALSE){

      
      my_res <- get(my_res_name)
      
      if(nrow(my_res) == 1){next()}
      
      my_res <- my_res[order(my_res$padj),]

      if(enriched){
            my_res <- my_res[my_res$log2FoldChange > 0,]
      } else {
            my_res <- my_res[my_res$log2FoldChange < 0,]
      }
            
      
      if(!(is.null(favorite_genes))){
            my_res <- my_res[rownames(my_res) %in% favorite_genes,]
      }
      
      if(nrow(my_res) > n_genes){
            my_res <- my_res[1:n_genes,]
      }
      
      ############################
      
      mat <- data_counts[rownames(data_counts) %in% rownames(my_res),]
      mat <- mat - rowMeans(mat)
      my_gene_id <- rownames(mat)
      rownames(mat) <- mapIds(org.Dm.eg.db, rownames(mat), "SYMBOL", keytype="FLYBASE", multiVals="first")
      rownames(mat)[is.na(rownames(mat))] <-  my_gene_id[is.na(rownames(mat))] 
      
      ############################      
      
      my_mat_size <- nrow(mat)
      my_ht_size <- (1/(my_mat_size/100))*4.5
      if(my_ht_size > 12){my_ht_size <- 12}
      
      ############################
      
      if(horizontal){
            mat <- t(mat)
      }
      
      hm = pheatmap(mat,
                     cellwidth = my_ht_size, cellheight = my_ht_size, fontsize = 8,
                     fontsize_row = my_ht_size-1,  fontsize_col = my_ht_size-1, 
                     cluster_cols = TRUE, cluster_rows = TRUE,
                     border_color = NA, silent = TRUE, 
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                     breaks = seq(-1*max(abs(mat)), max(abs(mat)), length.out = 101))
      
      return(hm)
}








############################################################################################################################
############################################################################################################################
############################################################################################################################


#####################################################      PCA      ######################################################## 






plottingPCA <- function(my_data, 
                        xcomp = 1,
                        ycomp = 2,
                        color_palette,
                        conditions,
                        quantiles = c(0,1),
                        show_labels = TRUE,
                        point_size = 1.5,
                        my_limits = c(-100,100)){
      
      rv <- rowVars(my_data)
      
      selection <- (rv >  quantile(rv, quantiles[1])  & rv < quantile(rv, quantiles[2]))
      
      
      pca <- prcomp(t(my_data[selection, ]), scale. = TRUE)
      
      percentVar <- round(pca$sdev^2/sum(pca$sdev^2)*100,1)[1:10]
      
      
      
      plot(pca$x[, xcomp], pca$x[, ycomp]*-1, 
           col = color_palette[conditions], 
           pch=16, cex = point_size,
           xlab = paste("PC",xcomp," (", percentVar[xcomp], "%)", sep=""),
           ylab = paste("PC",ycomp," (", percentVar[ycomp], "%)", sep=""),
           xlim= my_limits, ylim=my_limits)
      
      title(main = "PCA", line = 0.5)
      
      
      if(show_labels){
            text(pca$x[, xcomp], pca$x[, ycomp]*-1, labels = rownames(pca$x), 
                 adj = -0.5, col = "gray32", cex=0.5)      
      }
      
      
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################       dot plots      ################################################### 



plotDots <- function(my_data, 
                     reference_level = "Input",
                     my_favorite_gene,
                     color_palette,
                     color_groups,
                     conditions,
                     point_size = 1.5,
                     ylims = c(0, 15),
                     xlims = NULL,
                     x_label = NULL){
      
      
      my_gene_id <- my_genes$gene_id[grep(paste0("^",my_favorite_gene,"$"),my_genes$symbol)]
      
      my_data_favorite <- my_data[rownames(my_data)==my_gene_id,]
      
      if(length(my_data_favorite) > 0){
            
            my_data_favorite_input <- mean(my_data_favorite[names(my_data_favorite)   %in% SampleTable_Global$id[SampleTable_Global$type == reference_level]])
            my_data_favorite_norma <- my_data_favorite - my_data_favorite_input
            
            
            grouped_data = data.frame(exp = as.numeric(my_data_favorite_norma),
                                      group =  as.numeric(conditions))
            
            plot(grouped_data$group + runif(length(grouped_data$group), -0.05, 0.05),
                 grouped_data$exp,
                 pch=19, cex = point_size, 
                 ylim =  ylims, xlim = xlims,
                 xaxt = "n", xlab = "", ylab = "",
                 col =  color_palette[color_groups])
            
            title(main = my_favorite_gene, line = 0.5)
            
            if(!(is.null(x_label))){
                  axis(side = 1, at = seq_along(unique(grouped_data$group)), labels = x_label)
                  
            } else {
                  axis(side = 1, labels = FALSE)
                  
            }
            
            
      }
      

      
      
}








############################################################################################################################
############################################################################################################################
############################################################################################################################

###################################################       Heatmaps       ################################################### 


plotHeatmap <- function(my_mat,
                        my_row_order,
                        my_col_order,
                        min_value = 0,
                        max_value = 15,
                        my_title,
                        my_color_palette,
                        show_xaxis = FALSE,
                        useRaster = TRUE,
                        cex.main = 1.25){
      
      my_mat[my_mat < min_value] <- min_value
      my_mat[my_mat > max_value] <- max_value
      
      image(t(my_mat[my_row_order, my_col_order]), 
            #main = my_title, cex.main = cex.main,
            col = my_color_palette, 
            breaks =  seq(min_value, max_value, length.out = 101),
            axes=FALSE, useRaster = useRaster)
      
      title(main = my_title, line = 1.5, cex.main = cex.main)
      
      if(show_xaxis){
            axis(side = 1, at = c(0,0.5,1), lwd = 0, lwd.ticks = 1, las=1,tck = -0.15, cex.axis = 0.8,
                 labels = c( paste(round(as.integer(colnames(my_mat))[1] / 1000), "kb", sep=""), "",
                             paste("+",round(as.integer(colnames(my_mat))[ncol(my_mat)] / 1000), "kb", sep="")
                 ))
      }
      
      
}





plotHeatmapKey <- function(my_mat,
                           min_value = 0,
                           max_value = 15,
                           my_title,
                           my_color_palette,
                           show_yaxis = FALSE,
                           show_xaxis = FALSE,
                           useRaster = TRUE,
                           cex.main = 1.25,
                           cex.axis = 0.6){
      
      my_mat[my_mat < min_value] <- min_value
      my_mat[my_mat > max_value] <- max_value
      
      image((my_mat), 
            main = my_title, cex.main = cex.main,
            col = my_color_palette, 
            breaks =  seq(min_value, max_value, length.out = 101),
            axes=FALSE, useRaster = useRaster)
      
      if(show_yaxis){
            axis(side = 4, at = c(0,1), lwd = 0, lwd.ticks = 1, las=1, cex.axis = cex.axis, tck = -0.25,
                 labels = round(c(min_value,max_value)))
      }
      
      if(show_xaxis){
            axis(side = 1, at = c(0,1), lwd = 0, lwd.ticks = 1, las=1, cex.axis = cex.axis, tck = -0.25,
                 labels = round(c(min_value,max_value)))
      }
      
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################






############################################################################################################################
############################################################################################################################
############################################################################################################################







setupDDS <- function(my_comparison = "eMLE_RIP_S2_cells" #,
                     #ref_level = "Input"
                     ){
      
      SampleTable <- SampleTable_Global[grep(paste0("^",my_comparison,"$"), SampleTable_Global$comparison),]
      SampleTable <- SampleTable[order(SampleTable$id),]
      
      
      
      my_counts_Genes_Comparison <- my_counts_Genes[, colnames(my_counts_Genes) %in%  SampleTable$id]
      my_counts_Genes_Comparison <- my_counts_Genes_Comparison[,colnames(my_counts_Genes_Comparison)]
      
      
      stopifnot(
            identical(colnames(my_counts_Genes_Comparison), as.character(SampleTable$id))
      )
      
      
      
      ######################################################## 
      
      
      #filter <- TRUE
      
      filter <- apply(my_counts_Genes_Comparison, 1, function(x) length(x[x>1]) >= ncol(my_counts_Genes_Comparison)/1.333 )
      
      my_counts_Filtered <- my_counts_Genes_Comparison[filter,]
      
      
      
      ######################################################## 
      
      
      
      
      my_conditions <- SampleTable$name
      my_conditions <- gsub("\\+","plus", my_conditions)
      my_conditions <- gsub("\\-","minus",my_conditions)
      my_conditions <- factor(my_conditions)
      
      #my_conditions <- relevel(my_conditions, ref = ref_level)
      
      my_colData <- DataFrame(Sample = my_conditions,
                              Batch = factor(SampleTable$replicate))
      
      rownames(my_colData) <- SampleTable$id
      
      
      ######################################################## 
      
      dds <- DESeqDataSetFromMatrix(countData = my_counts_Filtered, 
                                     colData = my_colData, 
                                     design = ~Batch+Sample)
      
      
      
      dds <- DESeq(dds)
      
      return(dds)
}












############################################################################################################################
############################################################################################################################
############################################################################################################################





correctCounts <- function(dds){
      
      #rld <- rlog(dds, blind=FALSE)
      
      batchVar <- colData(dds)$Batch
      
      modcombat <- model.matrix(~Sample, data = colData(dds))
      
      #BatchCorrect <- ComBat(dat = assay(rld), batch = batchVar, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
      BatchCorrect <- ComBat(dat = log2(counts(dds, normalized = TRUE)+1), batch = batchVar, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
      
      return(BatchCorrect)
    
      
}



############################################################################################################################
############################################################################################################################
############################################################################################################################


calc_mean_ranges_counts <- function(ranges_mat, 
                                    select_features = c("snoRNA", "snRNA", "tRNA"),
                                    comparison = "MLE_titration_ATP_S2_RNA"){
      

      ranges_mat <- ranges_mat[rownames(ranges_mat) %in% select_features,]
      
      ranges_mat <-  t( t(ranges_mat) / colSums(ranges_mat) )
      
      my_subset <- grepl(comparison, SampleTable_Global$comparison)
      
      

      ranges_matches <- data.frame(Sample = SampleTable_Global[my_subset,]$name,
                                   Id =     SampleTable_Global[my_subset,]$id)
      
      
      ranges_mat <- ranges_mat[,colnames(ranges_mat) %in% as.character(ranges_matches$Id)]
      
      
      ranges_mat_mean <- as.data.frame(apply(ranges_mat, 1, FUN = function(y){ 
            merged_data <- merge(y, as.data.frame(ranges_matches), by.x="row.names", by.y="Id")
            agg_data <- aggregate(x ~ Sample,  data=  merged_data, FUN=mean) 
      }))
      
      
      rownames(ranges_mat_mean) <- ranges_mat_mean[,1]
      ranges_mat_mean <- ranges_mat_mean[, seq(2, ncol(ranges_mat_mean),2)]
      
      colnames(ranges_mat_mean) <- gsub("\\.x","", colnames(ranges_mat_mean))
      colnames(ranges_mat_mean) <- gsub("five","5' ", colnames(ranges_mat_mean))
      colnames(ranges_mat_mean) <- gsub("three","3' ", colnames(ranges_mat_mean))
      
      rownames(ranges_mat_mean) <- gsub("_"," ", rownames(ranges_mat_mean))
      rownames(ranges_mat_mean) <- gsub("\\."," ", rownames(ranges_mat_mean))
      
      
      return(ranges_mat_mean)
      
}








############################################################################################################################
############################################################################################################################
############################################################################################################################





getNumberOfDetectedReads <- function(sample_type = "Input"){
      
      
      SampleTable <- SampleTable_Global[SampleTable_Global$type == sample_type,]
      SampleTable <- SampleTable[!duplicated(SampleTable$id),]
      
      
      my_counts_Genes_Inputs <- my_counts_Genes[,colnames(my_counts_Genes) %in% SampleTable$id]

      
      detected_gene_counts <- apply(my_counts_Genes_Inputs, 2, function(x){ sum(x > 0) / sum(x) * 10^6})
      
      
      df <- merge(SampleTable, data.frame(detected_gene_counts), by.x = "id", by.y = "row.names")
      df <- df[order(df$source, decreasing = TRUE),]
      
      df$source <- gsub("clone8", "cl.8", df$source)
      
      return(df)
}




############################################################################################################################
############################################################################################################################
############################################################################################################################







cleanLegend <- function(my_legend){
      
      my_legend <- gsub("\\_","", my_legend)
      my_legend <- gsub("\\."," ", my_legend)
      my_legend <- gsub("plus","+", my_legend)
      my_legend <- gsub("minus","-", my_legend)
      my_legend <- gsub("4MSL","DCC", my_legend)
      
      
      return(my_legend)
}





############################################################################################################################
############################################################################################################################
############################################################################################################################




setupMatforHeatmap <- function(res_name = "res.HighMLE+ATP-Input",
                               log2_norm_counts_name = "log2_norm_counts.MLE_titration_ATP_S2_RNA",
                               pval_cutoff = 0.01,
                               my_selection = NULL,
                               n_hits = 100,
                               scaling = TRUE){
      
      
      if(!(is.null(my_selection))){
            
            selected_gene_ids <-  my_selection
            
      } else {
            
            my_res <- get(res_name)
            my_res.sign <- my_res[my_res$padj < pval_cutoff,]
            
            selected_gene_ids <- rownames(head(my_res.sign[order(my_res.sign$log2FoldChange, decreasing = TRUE),], n_hits))
      }
      
      
      
      log2_norm_counts <- get(log2_norm_counts_name)
      
      my_mat <- log2_norm_counts[rownames(log2_norm_counts) %in% selected_gene_ids,]
      
      my_row_names <- rownames(my_mat)
      rownames(my_mat) <- mapIds(org.Dm.eg.db, my_row_names, "SYMBOL", keytype="FLYBASE", multiVals="first")
      rownames(my_mat) <- gsub("lncRNA:", "",rownames(my_mat))
      rownames(my_mat)[is.na( rownames(my_mat))] <- my_row_names[is.na( rownames(my_mat))]
      
      if(scaling){
            my_mat <- my_mat - rowMeans(my_mat)
      }
     
      
      return(t(my_mat))
}








############################################################################################################################
############################################################################################################################
############################################################################################################################









