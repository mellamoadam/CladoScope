<details>
<summary>First we install and load necessary packages</summary>

```r
install.packages("stringr")
install.packages("smartsnp")
install.packages("adegenet")
install.packages("ggnewscale")
install.packages("zoo")
install.packages("dplyr")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("plotly")
install.packages("stats")
install.packages("pheatmap")
install.packages("ape")
install.packages("vegan")
install.packages("reshape2")
install.packages("cowplot")
install.packages("vcfR")
install.packages("ape")
install.packages("phangorn")
install.packages("ggtree")
install.packages("data.table")
install.packages("readxl")
install.packages("maps")
install.packages("sf")
install.packages("geojsonio")
install.packages("randomcoloR")
install.packages("gtools")
install.packages("rmarkdown")
install.packages("gtable")
install.packages("qpdf")
install.packages("scatterpie")
install.packages("remotes")
install.packages("scales")
install.packages('igraph', dependencies = TRUE)
install.packages('phytools')
install.packages('rnaturalearth')
remotes::install_github("liamrevell/phytools")
install.packages("devtools")
install_github("bbanbury/phrynomics") #For SNAPPER
devtools::install_github("TheWangLab/algatr")
install.packages("viridisLite")
devtools::install_github("TheWangLab/algatr")
install.packages("raster")
install.packages("geodata")
install.packages("rinat")
install.packages("Rcpp")
install.packages("terra")
install.packages("Biostrings")
install.packages("PipeMaster")
install.packages("phylotools")
install.packages("pdftools")
install.packages("magick")
system("brew install gdal")
install.packages("terra", type = "source", configure.args = "--with-proj-lib=$(brew --prefix)/lib/")
library(Biostrings)
library(PipeMaster)
devtools::install_github("DevonDeRaad/SNPfiltR")
install.packages("ggtree")
install.packages("rlang")
install.packages("patchwork")
if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("ggtree")

library(ggtree)
library(patchwork)  
library(scales)  
library(rinat)
library(SNPfiltR)
library(phylotools)
library(geodata)
library(Rcpp)
library(terra)
library(raster)
library(algatr)
library(viridisLite)
library(devtools)
library(rnaturalearth)
library(phrynomics)
library(phytools)
library(ggnewscale)
library(remotes)
library(scatterpie)
library(qpdf)
library(gtable)
library(tidyverse)
library(gtools)
library(randomcoloR)
library(geojsonio)
library(sf)
library(maps)
library(data.table)
library(phangorn)
library(vcfR)
library(cowplot)
library(vegan)
library(ape)
library(pheatmap)
library(stats)
library(tidyr)
library(plotly)
library(dplyr)
library(zoo)
library(adegenet) # If error then redo install.packages and say no instead of yes
library(rlang) # If error then redo install.packages and say no instead of yes
library(smartsnp)
library(stringr)
library(gridExtra)
library(grid)
library(reshape2)
library(readxl)
library(pdftools)
library(magick)
</details>




<details>
<summary>Now we can define and create paths to user folders and software execution files</summary>



# Main folders. Should end with a "/"
mainFolder="CladoScope/"
userPath = "/Users/adamaslam/"
folderPath=paste0(userPath, "Desktop/", mainFolder)
# folderPath="/Users/adamaslam/Desktop/Hypsiglena/Filtering/"
clusterFolderPath="/home/aaslam/" #If using cluster then this is relevant for some analyses that make files that reference other file locations


# Subfolders. Should end with a "/"
VCFPath="VCFs/"
TrawPath="Traw/"
PHYPath="PHYs/"
bedPath="Bed/"
admixturePath="Admixture/"
IQTreePath="Tree/"
ThreeDPath="3DPlots/"
pdfPath="OutputPDFs/"
nexusPath="Nexus/"
bppPath="BPP/"
testPath="Test/"
SVDQPath="SVDQ/"
binaryNexusPath="BinaryNexus/"
isolationPath="Isolation/"
gadmaPath="Gadma/"
DsuitePath="Dsuite/"
rawPath="Raw/"
tablePath="Tables/"
SoftwarePath="Software/"

# Software paths
phylipCoversionFolderPath=paste0(folderPath,SoftwarePath,"vcf2phylip-master")
pythonPath = "/opt/anaconda3/envs/tf/bin/python"
bppSWPathMac = "/Users/adamaslam/bpp-4.8.2-macos-aarch64/bin/bpp"
bppSWPathLinux = "bpp-4.8.2-linux-x86_64/bin/bpp"
gadmaSWPath = "/opt/anaconda3/bin/gadma"
gadmaClusterSWPath = "/home/aaslam/miniconda3/bin/gadma"
admixtureSWFolderPath=paste0(folderPath,SoftwarePath,"admixture_macosx-1.3.0/")
IQTreeSWPath="/opt/homebrew/bin/iqtree2"
easySFSPath="/Users/adamaslam/easySFS/easySFS.py"
DSuiteFBranchSWPath = "/Users/adamaslam/Dsuite/utils/dtools.py"
DtriosSWPath = "/Users/adamaslam/Dsuite/Build/Dsuite"

# Unfiltered raw VCF file name (e.g. ipyrad output)
rawFile="DeNovo.vcf"
rawString=str_sub(rawFile,end=-5)


knitr::opts_knit$set(root.dir = folderPath)


# Only run this at the very beginning of all analyses so results aren't deleted
# Commented out by default to avoid accidental deletion

# dir.create(file.path(folderPath), showWarnings = FALSE)
#
# dir.create(file.path(paste0(folderPath,VCFPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,TrawPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,PHYPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,bedPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,admixturePath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,IQTreePath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,ThreeDPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,nexusPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,bppPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,testPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,binaryNexusPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,isolationPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,gadmaPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,DsuitePath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,rawPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,tablePath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,SVDQPath)), showWarnings = FALSE)

# Create nested folders in pdfPath and admixturePath
# dir.create(file.path(paste0(folderPath,pdfPath,"FilterSelectionPlots")), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath,"FilterSelectionHeatmaps")), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath,"Trees")), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath,"AdmixtureBarPlots")), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath,"AdmixturePieMaps")), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath, isolationPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath, gadmaPath)), showWarnings = FALSE)
# dir.create(file.path(paste0(folderPath,pdfPath,"Tables/")), showWarnings = FALSE)

# dir.create(file.path(paste0(folderPath,admixturePath,"CVErrors")), showWarnings = FALSE)


</details>```
