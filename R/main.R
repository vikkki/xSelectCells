#' xSelectCells
#' @name xSelectCells
#' @title Select cells from clusters
#' @description  A light weight shiny app tool to enable manually selection of cells from seurat object based on its umap reduction
#' @param input_seurat_obj a seurat obj with umap reduction
#' @export
#' @author Ashley(X) Zhao
#' @examples xSelectCells(seurat_obj)


library(Seurat)
library(plotly)
library(shiny)
# library(ggplot2)
# library(Cairo)   # For nicer ggplot2 output when deployed on Linux

# obj = readRDS("~/Documents/Lizard_tail/GEX/harmony/tl_harmony_fibsub.rds")
seurat_obj <- c()
set.seed(4422)

## functions here ##
# 


#source("app.R")
xSelectCells <-function(input_seurat_obj) {
  #seurat_obj <<- input_seurat_obj
  
  ui <- server <- NULL # avoid NOTE about undefined globals
  
  my_color = c("#0B90AA","#7dce94","#B1B336","#04384A","#66638B","#D74B4B","#FF652D","#F6AE2D","#AE8D65","#D8DEAE","#70AB8F")
  getPalette = colorRampPalette(my_color)
  
  ui <- fluidPage(
    plotlyOutput("umap_for_brush",
                 width = "auto",
                 height = "600px",),
    sliderInput("point_size",
                "Point size:",
                min = 0.1,  max = 20, value = 1),
    hr(),
    
    div(downloadButton("dl_select_cells", "Download selected cell barcodes"), align = "left"),
    hr(),
    h4("Brushed barcodes"),
    verbatimTextOutput("barcode_brush_info")
    
  )
  
  server <- function(input, output) {
    
    output$umap_for_brush <- renderPlotly(umap_for_brush())
    
    sc_seurat_umap_inter_base <- reactive({
      sc = input_seurat_obj
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
    })
    
    umap_for_brush <- reactive({
      
      sc = sc_seurat_umap_inter_base()
      # -- axits settings
      ax <- list(
        zeroline = FALSE
        # gridcolor = #bdc3c7,
        # gridwidth = 1
      )
      
      marker = list(size = input$point_size)
      
      plot_ly(sc,type="scatter", mode = "markers",
              x = ~UMAP_1, y = ~UMAP_2,
              key = ~keys,
              color = ~ident,
              colors = getPalette(length(levels(sc$ident))),
              marker = marker,
              hoverinfo = "all",
              hovertext = paste0(sc$keys,"\n","nCount:",sc$nCount_RNA,"\n","nFeature:",sc$nFeature_RNA)) %>% layout(dragmode = "lasso", xaxis = ax, yaxis = ax)
      
    })
    
    cluster_brush_cells <- reactiveVal(NULL)
    
    output$barcode_brush_info <- renderPrint({
      d <- event_data("plotly_selected")
      if(is.null(d)) "Click and drag to select cells (double-click to clear)."
      else {
        cluster_brush_cells(d$key)
        paste0(nrow(d)," barcodes have been selected.")
      }
    })
    
    output$dl_select_cells <- downloadHandler(
      filename = function() {
        paste0("Barcode_selected", ".csv")
      },
      content = function(file) {
        barcodes = cluster_brush_cells()
        write.table(barcodes,file, quote = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
      }
    )
  }
  #source(file_path, local = TRUE)
  server_env <- environment(server)
  
  # Here you add any variables that your server can find
  server_env$input_seurat_obj <- input_seurat_obj
  server_env$param <- param
  
  app <- shiny::shinyApp(ui, server)
  shiny::runApp(app, ...)

}

