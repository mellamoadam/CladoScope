## Overview
This tutorial provides a walkthrough of the CladoScope framework, developed to address the challenges of resolving population structure in taxa with deep divergence, high diversity, and complex biogeography such as the genus Hypsiglena. In systems like this, traditional population genomics approaches often struggle due to extensive and hierarchical structure, where multiple levels of divergence and cryptic lineages obscure clear genetic assignment. CladoScope implements an iterative, subset-based strategy that begins with a user-defined population map informed by prior knowledge. This map guides initial filtering and clustering, which are then refined through rounds of phylogenetic and clustering analyses. As population groupings are reassessed and updated, the workflow enables increasingly accurate and biologically meaningful subsetting. These refined subsets allow for higher-resolution genomic analyses, offering clearer insight into evolutionary relationships and population boundaries. The approach is especially valuable for researchers working with complex, deeply structured taxa and provides a replicable method for dissecting difficult population genomic patterns.
<br>
<details>
<summary>First we install and load necessary packages</summary>
<br>
          
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
```
</details>


<br>
Now we can define and create local paths to user folders and software execution files that keep the results and temporary files organized. Must be adjusted to your own system.
<details>
<summary>Define and create paths</summary>
<br>
          
```r

################################# USER INPUTS  #################################
# Main folders. Should end with a "/"
mainFolder="CladoScope/"
userPath = "/Users/adamaslam/"
folderPath=paste0(userPath, "Desktop/", mainFolder)
# If using cluster then this is relevant for some analyses that make files that reference other file locations
clusterFolderPath="/home/aaslam/" 
################################################################################

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
```
</details>


<br>
Pulling in raw data files including `coordinateFile` that contains sample coordinates, VCF file, map files (.shp and geojson for boundaries, .tiff for terrain). Includes a query for iNaturalist if specified that will later include citizen science observations as gray points in background of DAPC map.
<details>
<summary>Pull raw data</summary>
<br>
          
```r

################################# USER INPUTS  #################################
# Maximum allowable tolerance on GPS position for iNaturalist query
maxGPSAccuracyInat = 30000
iNatQuery = "Hypsiglena"
################################################################################

coordinateFile="Updated_Specimen_data_sheet_11272023.xlsx"
coordinatePath=paste0(folderPath,rawPath,coordinateFile)

vcfAll=read.vcfR(paste0(folderPath,rawPath,rawFile))
vcfSamples=colnames(vcfAll@gt)[-1]


# Maps. shp files corresponding shx file in same folder for sf_read
MexicoMap = paste0(folderPath,rawPath,"states.geojson")
terrainMapPath = paste0(folderPath, rawPath, "NE2_HR_LC_SR_W.tif")
internalBordersMapPath = paste0(folderPath, rawPath, "ne_10m_admin_1_states_provinces_lines.shp") # Need corresponding shx file in same folder for sf_read 
countryBordersMapPath = paste0(folderPath, rawPath, "ne_10m_admin_0_map_subunits.shp") # Need corresponding shx file in same folder for sf_read 


# Search iNaturalist for Hypsiglena observations (only if it hasn't been done yet. very time inefficient)
inatCoordsFile = paste0(folderPath, rawPath, "inatCoords.csv")

if (!file.exists(inatCoordsFile)) {
  inatData = get_inat_obs(query= iNatQuery, maxresults = 10000)
  inatCoords = inatData[inatData$public_positional_accuracy < maxGPSAccuracyInat, c("latitude", "longitude")]
  inatCoords = na.omit(inatCoords)
  write.csv(inatCoords, file=paste0(folderPath, rawPath, "inatCoords.csv"), row.names=FALSE)
}
inatCoords = read.csv(inatCoordsFile)


```
</details>



<br>
We must convert terrain raster into RGB colors for map representation along with formatting coordinates file and map bounding boxes for plotting compatibility.
<details>
<summary>Map related organization</summary>
<br>

```r

################################# USER INPUTS  #################################
# Change colnames and create df
coordinates=read_excel(coordinatePath)
names(coordinates)[names(coordinates) == 'Name_In_Seq_files']='ID' 
names(coordinates)[names(coordinates) == 'Lat']='Latitude' 
names(coordinates)[names(coordinates) == 'Long']='Longitude' 
coordinates=coordinates[,c('ID','Latitude','Longitude')]
coordinates=as.data.frame(coordinates)

# Estimating coordinates of outgroups from LSU museum specimen notes. Vague description on HWY 200 some # of km E of Guererro/Michoacan line. 
# Coordinates below likely within 5km radius
coordinates$Latitude[coordinates$ID == "LSUMZ_39534_Pseudoleptodeira_latifasciata"] = 18.02
coordinates$Longitude[coordinates$ID == "LSUMZ_39534_Pseudoleptodeira_latifasciata"] = -102.18
coordinates$Latitude[coordinates$ID=="LSUMZ_39571_latifasciata_latifasciata"]=17.99
coordinates$Longitude[coordinates$ID=="LSUMZ_39571_latifasciata_latifasciata"]=-102.15

# Change IDs to match VCF File IDs
coordinates$ID[coordinates$ID=="UTAR_52345_chlorophaea_chlorophaea"]="UTAR_52345F_chlorophaea_chlorophaea"
coordinates$ID[coordinates$ID=="UTAR_52350_jani_texana"]="UTAR_52350F_jani_texana"


# Define bounds for terrain map, with degrees of buffer in each direction
mapBoundsBuffer = 1.2

# If plotting is an issue due to map size, scale down (1 = no scaling)
downsampleTerrainMapFactor = 1

################################################################################

terrainMapBounds = c(min(coordinates$Longitude) - mapBoundsBuffer, max(coordinates$Longitude) + mapBoundsBuffer, min(coordinates$Latitude) - mapBoundsBuffer, max(coordinates$Latitude) + mapBoundsBuffer)  

terrainRaster = rast(terrainMapPath)
terrainRaster = crop(terrainRaster, terrainMapBounds)

if(downsampleTerrainMapFactor > 1){
  terrainRaster = aggregate(terrainRaster, fact = downsampleTerrainMapFactor, fun = mean) 
}

terrainDF = as.data.frame(terrainRaster, xy = TRUE)
colnames(terrainDF) = c("Longitude", "Latitude", "R", "G", "B")

# Convert RGB (0-255) to color. This map uses RGB columns that we convert to color here.
terrainDF$Color = rgb(terrainDF$R / 255, terrainDF$G / 255, terrainDF$B / 255)

# Now read in borders and crop. Internal are states/provinces
internalBordersMap = read_sf(internalBordersMapPath)
internalBordersBbox = st_bbox(c(xmin = terrainMapBounds[1], xmax = terrainMapBounds[2], 
                  ymin = terrainMapBounds[3], ymax = terrainMapBounds[4]), 
                crs = st_crs(internalBordersMap))  

internalBordersMap = st_crop(internalBordersMap, internalBordersBbox)  


countryBordersMap = read_sf(countryBordersMapPath)
countryBordersBbox = st_bbox(c(xmin = terrainMapBounds[1], xmax = terrainMapBounds[2], 
                  ymin = terrainMapBounds[3], ymax = terrainMapBounds[4]), 
                crs = st_crs(countryBordersMap))  

countryBordersMap = st_crop(countryBordersMap, countryBordersBbox)  

```
</details>



<br>
<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>





<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>





<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>





<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>





<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>





<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>





<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>













