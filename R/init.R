data("pbmc3k")

plotData <- as.data.frame(pbmc3k[["tsne"]]@cell.embeddings)
plotData$cluster <- pbmc3k$seurat_clusters


ggplot(plotData, aes(x = tSNE_1, y = tSNE_2, fill = cluster, color = cluster)) +
  stat_unchull(alpha = 0.3, size = 0.5) +
  geom_point(size = 0.5) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line()
  )



library(Seurat)
library(plotly)
library(shiny)

seurat_obj <- c()
set.seed(4422)

## functions here ##
# 

my_color = c("#0B90AA","#7dce94","#B1B336","#04384A","#66638B","#D74B4B","#FF652D","#F6AE2D","#AE8D65","#D8DEAE","#70AB8F")

getPalette = colorRampPalette(my_color)



#source("app.R")
xSelectCells <-function(input_seurat_obj,point_size = 2.5) {
  #seurat_obj <<- input_seurat_obj
  sc = sc_seurat_umap_inter_base(input_seurat_obj)
  dim(sc)
  ident = as.data.frame(sc@active.ident)
  colnames(ident) <- "ident"
  embeds = as.data.frame(Embeddings(sc[["umap"]]),col.names = T)
  embeds = cbind.data.frame(embeds, ident)
  embeds$nCount_RNA  = sc@meta.data[["nCount_RNA"]]
  embeds$nFeature_RNA = sc@meta.data[["nFeature_RNA"]]
  embeds$cluster <- sc@meta.data[["seurat_clusters"]]
  embeds$keys <- rownames(embeds)
  ax <- list(
    zeroline = FALSE
    # gridcolor = #bdc3c7,
    # gridwidth = 1
  )
  marker = list(size = point_size)
  plot_ly(sc,type="scatter", mode = "markers",
          x = ~UMAP_1, y = ~UMAP_2,
          key = ~keys,
          color = ~ident,
          colors = getPalette(length(levels(sc$ident))),
          marker = marker,
          hoverinfo = "all",
          hovertext = paste0(sc$keys,"\n","nCount:",sc$nCount_RNA,"\n","nFeature:",sc$nFeature_RNA)) %>% 
            layout(dragmode = "lasso", xaxis = ax, yaxis = ax)
  
}

sc_seurat_umap_inter_base <- function(sc){
  
  dim(sc)
  ident = as.data.frame(sc@active.ident)
  colnames(ident) <- "ident"
  embeds = as.data.frame(Embeddings(sc[["umap"]]),col.names = T)
  embeds = cbind.data.frame(embeds, ident)
  embeds$nCount_RNA  = sc@meta.data[["nCount_RNA"]]
  embeds$nFeature_RNA = sc@meta.data[["nFeature_RNA"]]
  embeds$cluster <- sc@meta.data[["seurat_clusters"]]
  embeds$keys <- rownames(embeds)
  
  return(embeds)
}

    d <- event_data("plotly_selected")
    if(is.null(d)) "Click and drag to select cells (double-click to clear)."
    else {
      cluster_brush_cells= d$key
      paste0(nrow(d)," barcodes have been selected.")
    }
  


