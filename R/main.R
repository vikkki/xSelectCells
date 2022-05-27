#' xSelectCells
#' @name xSelectCells
#' @title Select cells from clusters
#' @description  A light weight shiny app tool to enable manually selection of cells from seurat object based on its umap reduction
#' @param input_seurat_obj a seurat obj with umap reduction
#' @export
#' @author Ashley(X) Zhao
#' @examples xSelectCells(seurat_obj)


library(Seurat)
# library(ggplot2)
# library(Cairo)   # For nicer ggplot2 output when deployed on Linux

# obj = readRDS("~/Documents/Lizard_tail/GEX/harmony/tl_harmony_fibsub.rds")
seurat_obj <- c()
set.seed(4422)

## functions here ##
my_color = c("#0B90AA","#7dce94","#B1B336","#04384A","#66638B","#D74B4B","#FF652D","#F6AE2D","#AE8D65","#D8DEAE","#70AB8F")
getPalette = colorRampPalette(my_color)


#source("app.R")
xSelectCells <-function(input_seurat_obj) {
  #seurat_obj <<- input_seurat_obj
  file_path <- system.file("app.R", package = "xSelectCells")
  if (!nzchar(file_path)) stop("Shiny app not found")
  ui <- server <- NULL # avoid NOTE about undefined globals
  source(file_path, local = TRUE)
  server_env <- environment(server)
  
  # Here you add any variables that your server can find
  server_env$input_seurat_obj <- input_seurat_obj
  server_env$param <- param
  
  app <- shiny::shinyApp(ui, server)
  shiny::runApp(app, ...)

}

