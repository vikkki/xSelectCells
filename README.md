# xSelectCells
A light weight shiny app tool to enable manually selection of cells from seurat object.

## Install
This mini app is based on Shiny, Seurat and plotly

```
# we'll need to install plotly, dplyr, seurat and shiny first if not already
devtools::install_github('vikkki/xSelectCells@main')
```
## Usage

To start, you can prepare your own data (normal single cell expression or spatial data) or use the example dataset included in the package by:

```
seurat_obj <- data(tiny_barin, package = "xSelectCells")
```
It's a sub-dataset of mouse anterior from 10X Genomics Visium Sagittal Mouse Brain [here](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Mouse_Brain_Sagittal_Anterior) and imported by ```stxBrain.SeuratData```.

### UMAP based selection

For the input seurat object, it should include umap reductions (sc[["umap"]]) and UMAP_ as key. Then to lauch shinyapp:

```
xSelectCells::xSelectCells(seurat_obj)
# or
barcodes <- xSelectCells::xSelectCells(seurat_obj)
```
![interface](https://raw.githubusercontent.com/vikkki/xSelectCells/main/img/xs1.png)
888
Select cells with lasso or rectangular tools, and download selected cells anytime by click "Download" button:

![select](https://raw.githubusercontent.com/vikkki/xSelectCells/main/img/xs2.png)
***

![download](https://raw.githubusercontent.com/vikkki/xSelectCells/main/img/xs3.png)
***

When the selection is done, click "Confirm" would end this shiny app session and return the final list of cells.

```R
Listening on http://127.0.0.1:3269
> barcodes
 [1] "GTAGCTAAGTATTGCC-7"  "AGGTGTTCATGGAATA-7"  "TTACCGCAGAGCAACC-10" "GCTTTCGAGGCAGGTT-1" 
 [5] "ATGCATGCACAAATCC-4"  "ACAACCAGTTACGGAG-5"  "CTTCTCTAGGCGATAC-5"  "GAGACTTTCACTGCTC-5" 
 [9] "AACAACCGTACCCACG-6"  "AGTGCCGAGGTCGCCT-6"  "ATCCTATCACGTTGGC-6"  "ATTCCCGCAAGAGTAT-6" 
...
```
### Spatial based selection
To use the information of coordinates:

```
xSelectCells::xSelectCells(xSelectCells(secB, type = "spatial", img = "secB")
```
Here we the image name under Seurat object's images slot as parameter ```img```, then shinyApp would plot spots associated to this image's coordinates information. (In the sample below it's "secB")

<img src="https://raw.githubusercontent.com/vikkki/xSelectCells/main/img/sx5.png" alt="drawing" width="50%"/>

Select spots:

![img](https://raw.githubusercontent.com/vikkki/xSelectCells/main/img/xs4.png)


### Then ...
Next, you can make a new seurat object by:
```
seurat_obj_sub <- subset(seurat_obj, cells = barcodes)
```

Highlight them in feature plot:
```
DimPlot(seurat_obj, cells.highlight = c(barcode))
```
And any further explorations.

## Future features
I'm still working on extracting coor info from spatial obj, expecially whenther is more than one slides. Really hope that this mini tool would save your a cup of coffe time during research! And any suggestion and idea would be more than welcome 🌶🌶🌶
