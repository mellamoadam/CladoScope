# Overview
This tutorial provides a walkthrough of the CladoScope framework, developed to address the challenges of resolving population structure in taxa with deep divergence, high diversity, and complex biogeography such as the genus *Hypsiglena*. In systems like this, traditional population genomics approaches often struggle due to extensive and hierarchical structure, where multiple levels of divergence and cryptic lineages obscure clear genetic assignment. 
<br>
<br>
CladoScope implements an iterative, subset-based strategy that begins with a user-defined population map informed by prior knowledge. This map guides initial filtering and clustering, which are then refined through rounds of phylogenetic and clustering analyses. Performing the analysis in this iterative way is important because it provides a starting point to run the code, evaluate which filtering combinations yield the most optimal results, and then refine the population map accordingly before re-running the analysis for improved accuracy. 
<br>
<br>
As population groupings are reassessed and updated, the workflow enables increasingly accurate and biologically meaningful subsetting. These refined subsets allow for higher-resolution genomic analyses, offering clearer insight into evolutionary relationships and population boundaries. The approach is especially valuable for researchers working with complex, deeply structured taxa and provides a replicable method for dissecting difficult population genomic patterns.
<br>
<br>

# Data Preparation
## Load Packages
<details>
<summary>Package loading code</summary>
<br>
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

## Define and Create Paths
Now we can define and create local paths to user folders and software execution files that keep the results and temporary files organized. Must be adjusted to your own system.
<br>
<br>
<details>
<summary>Define and create paths code</summary>
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

## Pull Raw Data
Pulling in raw data files including `coordinateFile` that contains sample coordinates, VCF file, map files (.shp and geojson for boundaries, .tiff for terrain). Includes a query for iNaturalist if specified that will later include citizen science observations as gray points in background of DAPC map.
<br>
<br>
<details>
<summary>Pull raw data code</summary>
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

## Organize Map-Related Files
We must convert terrain raster into RGB colors for map representation along with formatting coordinates file and map bounding boxes for plotting compatibility. We also subset maps for the political boundaries necessary for plotting our samples.
<br>
<br>
<details>
<summary>Map-related organization code</summary>
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



# Now set the US states to include in the plot. No raw data file required. Only applicable for some species of course.
UsaStatesMap=c("california", "nevada", "arizona", "utah", "colorado", "new mexico", "texas", "oklahoma", "oregon", "washington", "idaho", "kansas")

OtherCountriesInMap=c("Mexico")
# Now set the Mx states to include in the plot. No raw data file required. Only applicable for some species of course.
MexicoStatesMap=c("Aguascalientes", "Baja California", "Baja California Sur", "Chihuahua", "Coahuila de Zaragoza", "Durango", "Guanajuato", "Hidalgo", "Jalisco", "Mexico", "Mexico City", "Michoacán", "Nayarit", "Nuevo León", "Querétaro", "San Luis Potosí", "Sinaloa", "Sonora", "Tamaulipas", "Zacatecas")

UsMap = map_data("state") %>% filter(region %in% UsaStatesMap)

mexicanStates = st_read(MexicoMap) %>% filter(state_name %in% MexicoStatesMap)
world = ne_countries(scale = "medium", returnclass = "sf")
mexico = world[world$name == OtherCountriesInMap, ]

```
</details>

## Define and Iterate Population Map
This is a section that will be iterated by the user after subsequent analyses. These variables can be empty during the initial run, or input based on some prior knowledge. Under the section USER INPUTS, the variables to change after subsequent analyses are the following:
          1. `samplesToRemove`: this variable can be manually changed if a fundamental error with a sample is observed.
          2. `hybridSamples`: this variable can be manually changed after running concatenated trees, PCA, DAPC, and ADMIXTURE to define a sample as being between observed populations.
          3. `popChanges`: this variable can be manually changed after running concatenated trees, PCA, and DAPC to reclassify a sample's population.
<br>
<br>
Note that the USER INPUTS section of the code below reflects population assignments after iterative assignment. Initially, `hybridSamples` and `popChanges` were input as empty vectors. After running concatenated trees, PCA, DAPC, and ADMIXTURE, samples that were better classified as being in other or between two populations became more clear. Rerun these steps until the population map has converged and remains unchanged. In the example below, this required three iterations, but may take a few more for more complex systems.
<br>
<br>
<details>
<summary>Population map code</summary>
<br>       
 
```r

################################# USER INPUTS  #################################
# Remove the following sample(s). Can be appended to in future analyses.
samplesToRemove = c("MVZ_230713_ochrorhyncha_nuchalata") # Museum specimen location error

outgroupPopulation = outgroupPopName = "Pseudoleptodeira_latifasciata"
outgroupSamples = c("LSUMZ_39571_latifasciata_latifasciata", "LSUMZ_39534_Pseudoleptodeira_latifasciata")


# Hybrids IDed from concatenated trees in future step
hybridSamples = c("KWS_245_spnov1", # Wilcox, Arizona. Concatenated tree places just outside H. jani.
                  "CAS_223504_chlorophaea_deserticola", # Anza Borrego area. In H. chlorophaea subsets, placed with E Sierran (deserticoloa West) samples.
                  "CAS_228973_ochrorhyncha_klauberi", # Anza Borrego area. In H. ochrorhyncha subset, placed with other klauberi samples.
                  "MVZ_236389_ochrorhyncha_klauberi", # Mid-Baja Norte. Admixture places about half ancestry from SoCal klauberi population, DAPC places there too.
                  "MVZ_180265_chlorophaea_loreala", # Four corners area. Sometimes placed outside H. jani, sometimes with H. jani (Western Clade), but definitely not H. chlorophaea loreala.
                  "AMNH_R_504774_jani_jani" # Central Mexico, Querétaro. Similar to KWS_245_spnov1, placed outside of H. jani
                  # "CAS_229952_chlorophaea_deserticola" # Might be a hybrid between deserticola (Eastern Clade) and deserticola (Western Clade)
                  )

# Some samples don't fit well with the ID in it's sample ID name based on analyses like concatenated trees or PCA. Rename them here using the format [sample ID]POP="[correct population]"
popChanges = c(
  "MVZ_226235_spnov1xtexana" = "H_jani_texana",
  "CAS_228935_chlorophaea_chlorophaea" = "sp_nov_1",
  "BYU_5630_chlorophaea_loreala" = "H_chlorophaea_deserticoloa",
  "JRO_694_chlorophaea_chlorophaea" = "sp_nov_2", # Fits better with Sonora samples in IQ-Tree but shares same admixture ancestry and DAPC group as sp_nov_2. Some IQ-Trees have it as a monophyletic group.
  "BYU_42832_chlorophaea_chlorophaea" = "sp_nov_2", # Fits better with Sonora samples in IQ-Tree but shares same admixture ancestry and DAPC group as sp_nov_2. Some IQ-Trees have it as a monophyletic group.
  "MF_21713_jani_jani" = "H_jani_texana", # Fits in well with H. jani (Eastern Clade)
  "MF_21732_jani_dunklei" = "H_jani_texana" # Fits in well with H. jani (Eastern Clade)
)


# If any sample ID contains a string that can be grabbed and used to define that group, but a different name is preffered, rename with popRename
popRename = c(
  "latifasciata" = "Pseudoleptodeira_latifasciata",
  "jani_texana" = "H_jani_texana",
  "jani_jani" = "H_jani_jani",
  "jani_dunklei" = "H_jani_dunklei",
  "chlorophaea_loreala" = "H_chlorophaea_loreala", 
  "chlorophaea_chlorophaea" = "H_chlorophaea_chlorophaea",  
  "chlorophaea_deserticola" = "H_chlorophaea_deserticola",  
  "ochrorhyncha_nuchalata" = "H_ochrorhyncha_nuchalata",  
  "ochrorhyncha_ochrorhyncha" = "H_ochrorhyncha_ochrorhyncha", 
  "ochrorhyncha_baueri" = "H_ochrorhyncha_baueri", 
  "ochrorhyncha_klauberi" = "H_ochrorhyncha_klauberi", 
  "unaocularus_unaocularus" = "H_unaocularus", 
  "torquata_torquata" = "H_torquata", 
  "affinis_affinis" = "H_affinis", 
  "catalinae_catalinae" = "H_catalinae", 
  "slevini" = "H_slevini", 
  "tanzeri_tanzeri" = "H_tanzeri", 
  "spnov1" = "sp_nov_1", 
  "spnov2" = "sp_nov_2"
)

# Create population map matchDF
# All populations (besides spnov1xtexana because it interferes with the for loops below since grepl finds spnov1 twice) for matching with raw data and future grouping in plots. popNames matches strings contained in the samples to the populations in the original names of popRename.
popNames = c("jani_texana",
           "jani_jani",
           "jani_dunklei",
           "tanzeri_tanzeri",
           "chlorophaea_chlorophaea",
           "chlorophaea_loreala",
           "chlorophaea_deserticola",
           "torquata_torquata",
           "ochrorhyncha_ochrorhyncha",
           "ochrorhyncha_baueri",
           "ochrorhyncha_klauberi",
           "ochrorhyncha_nuchalata",
           "unaocularus_unaocularus",
           "affinis_affinis",
           "catalinae_catalinae",
           "slevini",
           "latifasciata",
           "spnov1",
           "spnov2"
           ) 

################################################################################

# Initialize popmap with IDs
matchDF = data.frame(Sample = vcfSamples) 


# Match population names with the sample ID name
matchDF$Population = sapply(matchDF$Sample, function(m) {
  popNames[which(sapply(popNames, grepl, x=m, fixed=TRUE))[1]] 
})


# Match population names with the sample ID name for samples that don't belong in their assigned population using variable popChanges
matchDF$Population[matchDF$Sample %in% names(popChanges)] = popChanges[matchDF$Sample[matchDF$Sample %in% names(popChanges)]]


# Change population names using variable popRename
matchDF$Population[matchDF$Population %in% names(popRename)] = popRename[matchDF$Population[matchDF$Population %in% names(popRename)]]

# Change population names of any hybrids to "H_hybrid"
matchDF$Population = ifelse(matchDF$Sample %in% hybridSamples, "H_hybrid", matchDF$Population)


PopmapSummaryPDFPath = paste0(folderPath, pdfPath, tablePath, rawString, "PopmapSummary", ".pdf")
pdf(file = PopmapSummaryPDFPath)
grid.table(setNames(data.frame(table(matchDF$Population)), c("Population", "Count")) %>% arrange(desc(Count)), rows = NULL, theme = ttheme_minimal(base_size = 8, base_colour = "black"))
dev.off()

matchDF = matchDF[!(matchDF$Sample %in% samplesToRemove),]
rownames(matchDF)=NULL


matchDFTxtFilePath = file.path(folderPath,tablePath,paste0(rawString,"PopMap.txt"))
unlink(matchDFTxtFilePath,recursive=TRUE,force=TRUE)
write.table(
  matchDF, 
  file=matchDFTxtFilePath, sep="\t", quote=FALSE, 
  row.names=FALSE, col.names=FALSE
)

```
</details>


# Filter Triage
This section outlines how to perform VCF filtering using a grid-based approach, systematically varying thresholds across multiple parameters (i.e. MAC, MAF, individual missingness, site missingness, depth). By testing a wide range of combinations, we can identify filter sets that produce the most stable and informative outputs for downstream analysis. A critical component is defining `popsOfInterest`, which can start as a hypothesis based on monophyletic groups or closely related populations and then be refined iteratively as results from the concatenated tree, DAPC, and PCA are reviewed. 
<details>
<br>
<br>
<summary>Filtering triage code</summary>
<br>       
          
```r
################################# USER INPUTS  #################################
removeIndels = FALSE # Doesn't change VCF because all SNPs
replaceMultiallelic = TRUE # Used to change the genotype of each indiv at each locus if there is an allele that is not ref or alt present.
missingSitesThreshold = seq(0.6,0.9, by = 0.1) # Lower is more stringent. This is max missing sites threshold
missingIndivsThreshold = seq(0.8,0.95, by = 0.03) # This is the threshold used for all samples besides outgroupSamples. If any sample in outgroupSamples has more than missingIndivsThreshold, it is still included in subset. Only applicable when includeOutgroup = TRUE.
MACMinThreshold = 3
MAFMinThreshold = 0.05
MinMeanDPThreshold = 6
MaxMeanDPThreshold = 9
MinDPThreshold = 5
MaxDPThreshold = 12
thinningFilter = 200
includeOutgroup = c(TRUE,FALSE)


# List of populations of interest to examine
popsOfInterest = list(c("ochrorhyncha"),
  c("jani"),
  c("chlorophaea"),
  c("chlorophaea","sp_nov_1","sp_nov_2"),
  c("torquata","chlorophaea","sp_nov_1","sp_nov_2","unaocularus","catalinae"),
  c("torquata","ochrorhyncha","chlorophaea","jani","sp_nov","slevini","tanzeri","affinis","catalinae","unaocularus", "hybrid")
)

names(popsOfInterest) = c("Ochrorhyncha", 
                           "Jani", 
                           "Chlorophaea", 
                           "ChlorophaeaSpNov1SpNov2",
                           "ChlorophaeaSpNov1SpNov2TorquataUnaocularusCatalinae",
                           "All"
)

################################################################################         

# Define function that will append information about the VCF at each filter step to lists
updateCounts=function(filterStep) {
  samplesRemaining <<- c(samplesRemaining, ncol(vcfSubset@gt) - 1)
  SNPsRemaining <<- c(SNPsRemaining, nrow(vcfSubset@fix))
  filterValue <<- c(filterValue,filterStep)

    if (is.character(filterStep) && length(filterStep) == 1) {
    filterName <<- c(filterName,filterStep)  # Already a string, return as is
  } else {
    filterName <<- c(filterName,deparse(substitute(filterStep))) # Extract variable name as a string
  }
}


# Initialize lists for Gadma seqL
VCFFileNameList = seqLList = c()


pdf(paste0(folderPath,pdfPath,"FilterTriage.pdf"))

for (populations in 1:length(popsOfInterest)) {
  for (missingSite in missingSitesThreshold){
    for (OutgroupIncludedTF in includeOutgroup){

        pops = popsOfInterest[[populations]]
        currentPops = paste(pops, collapse = "")
        popName = names(popsOfInterest)[populations]


        popsDF = pops
        if (OutgroupIncludedTF) { # Include Pseudoleptodeira if specified to
          popsDF = c(popsDF,"latifasciata")
        }


        currentSubset = matchDF$Sample[grep(paste(popsDF,collapse="|"), matchDF$Population, ignore.case = TRUE)]

        # Subset the VCF to keep only the matching samples
        vcfSubset = vcfAll[, c(TRUE, colnames(vcfAll@gt)[-1] %in% currentSubset)] # Subset cols, appending TRUE to beginning to include 'FORMAT' header col

        samplesRemaining = c(ncol(vcfSubset@gt) - 1)
        SNPsRemaining = c(nrow(vcfSubset@fix))
        filterName = c("Pre-Filtering")
        filterValue = c("NA")

        #First remove indels. Only SNP data so no indels
        if (removeIndels) {
          vcfSubset = extract.indels(vcfSubset)
          updateCounts(removeIndels)
        }

          if (replaceMultiallelic) {
            # Find sites with multiple ALT alleles
            multiallelicSites = grepl(",", vcfSubset@fix[,"ALT"])

##################################### PLOT #####################################
              multiallelicSitesDF=data.frame(
              Category = c("Biallelic", "Multiallelic"),
              Count = as.numeric(table(multiallelicSites))
            )

            multiallelicSitesDF$Percent = round(multiallelicSitesDF$Count / sum(multiallelicSitesDF$Count) * 100, 1)
            numSNPsRemoved=multiallelicSitesDF$Count[multiallelicSitesDF=="Multiallelic"]
            percentSNPsRemoved=multiallelicSitesDF$Percent[multiallelicSitesDF=="Multiallelic"]

              p1 = ggplot(multiallelicSitesDF, aes(x = "", y = Count, fill = Category)) +
              geom_bar(stat = "identity", width = 1) +
                scale_fill_manual(values = c("Multiallelic" = "coral", "Biallelic" = "steelblue")) +
              coord_polar(theta = "y") +
              theme_void() +
              labs(title = "Proportion of Multiallelic vs Biallelic Sites",
                  subtitle = paste0("SNPs Before Filtering: ",last(SNPsRemaining),
                  "\nFilter Removes ",numSNPsRemoved," SNPs (",percentSNPsRemoved,"%)")) +
              geom_text(aes(label = paste0(Percent, "%")),
                        position = position_stack(vjust = 0.5), size = 4, fontface="bold") +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"))


################################################################################

          vcfSubset=vcfSubset[!multiallelicSites,] # Omit multiallelicSites sites

          updateCounts(replaceMultiallelic)
        }

        dp=extract.gt(vcfSubset, element="DP", as.numeric=TRUE)

        meanDepthCheck=apply(dp,1,function(x) mean(x[x>0],na.rm=TRUE)) #Mean depth for each site after dropping NAs and 0s

        # Identify sites where mean depth exceeds MaxMeanDPThreshold and replace with missing data. Don't do anything if all NAs in row. They will be deleted in a future step.
        highMeanDepthSites=!is.na(meanDepthCheck) & meanDepthCheck>=MaxMeanDPThreshold

        if (!is.na(MinMeanDPThreshold)) {
            lowMeanDepthSites=!is.na(meanDepthCheck) & meanDepthCheck<=MinMeanDPThreshold
        }else{
          lowMeanDepthSites=0
        }


##################################### PLOT #####################################

        lowerXLim=quantile(meanDepthCheck,.000001,na.rm=TRUE) # Cutoff for chart
        upperXLim=quantile(meanDepthCheck,.99,na.rm=TRUE) # Cutoff for chart

        numSNPsRemovedMaxMeanDp=sum(highMeanDepthSites)
        percentSNPsRemovedMaxMeanDp=round(sum(highMeanDepthSites)/last(SNPsRemaining)*100,1)

        numSNPsRemovedMinMeanDp=sum(lowMeanDepthSites)
        percentSNPsRemovedMinMeanDp=round(sum(lowMeanDepthSites)/last(SNPsRemaining)*100,1)

        p2 = ggplot(data.frame(meanDepthCheck), aes(x = meanDepthCheck)) +
        geom_histogram(bins = 50, fill = "steelblue", alpha = 0.5, color = "black") +
        geom_vline(xintercept = c(MinMeanDPThreshold,MaxMeanDPThreshold),
        col = "coral", linetype = "dashed", linewidth = 1) +
        labs(title = "Mean Depth Distribution",
        subtitle = paste0("SNPs Before Filtering: ",last(SNPsRemaining),
                          "\nMin Mean Depth Threshold: ", MinMeanDPThreshold,
                          "\nMin Mean Depth Removes ",numSNPsRemovedMinMeanDp," SNPs (",percentSNPsRemovedMinMeanDp,"%)",
                          "\nMax Mean Depth Threshold: ", MaxMeanDPThreshold,
                          "\nMax Mean Depth Removes ",numSNPsRemovedMaxMeanDp," SNPs (",percentSNPsRemovedMaxMeanDp,"%)"),
        x = "Mean Depth",y = "Count") +
        theme_minimal(base_size = 14) +
        xlim(lowerXLim, upperXLim) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        panel.grid.major = element_line(color = "grey90"))

################################################################################

        # Takes up a lot of space so remove
        rm(meanDepthCheck)

        # Update lists manually because applying both filters at the same time
        samplesRemaining=c(samplesRemaining,last(samplesRemaining))
        SNPsRemaining=c(SNPsRemaining,last(SNPsRemaining)-sum(highMeanDepthSites))
        filterName=c(filterName,"MaxMeanDPThreshold")
        filterValue=c(filterValue,MaxMeanDPThreshold)

        if (!is.na(MinMeanDPThreshold)) {
          # Update lists manually because applying both filters at the same time
          samplesRemaining=c(samplesRemaining,last(samplesRemaining))
          SNPsRemaining=c(SNPsRemaining,last(SNPsRemaining)-sum(lowMeanDepthSites))
          filterName=c(filterName,"MinMeanDPThreshold")
          filterValue=c(filterValue,MinMeanDPThreshold)
        }

        vcfSubset=vcfSubset[!highMeanDepthSites & !lowMeanDepthSites,] # Omit high and low Mean Depth Sites


        #Identify samples that fail depth thresholds. Ignore the 0s cause they are already "./.:0:0,0,0,0"
        depthViolations=(dp>0) & (dp<MinDPThreshold | dp>MaxDPThreshold)


##################################### PLOT #####################################

        depthDF=data.frame(Depth = as.vector(dp))
        downSampleVal=1e7
        depthLen=nrow(depthDF)
        # Randomly sample from depthDF because p3 can be too big to plot if not
        if (depthLen>downSampleVal){
          depthDFDownSampled = depthDF[sample(nrow(depthDF), size = 1e6), , drop = FALSE]
        }else{
          depthDFDownSampled=depthDF
        }
        lowerXLim=-1
        upperXLim=quantile(depthDFDownSampled,0.99,na.rm=TRUE)

        SNPxSamples=last(SNPsRemaining)*last(samplesRemaining)
        numSNPsRemovedMinDp=sum( (dp>0) & (dp<MinDPThreshold) )
        percentSNPsRemovedMinDp=round(numSNPsRemovedMinDp/SNPxSamples*100,2)

        numSNPsRemovedMaxDp=sum((dp>MaxDPThreshold) )
        percentSNPsRemovedMaxDp=round(numSNPsRemovedMaxDp/SNPxSamples*100,1)

        p3 = ggplot(depthDFDownSampled, aes(x = Depth)) +
        geom_histogram(bins = as.numeric(upperXLim-lowerXLim), fill="steelblue", alpha=0.5, color="black") +
        geom_vline(xintercept = c(MinDPThreshold,MaxDPThreshold),
        col = "coral", linetype = "dashed", linewidth = 1) +
        labs(title = paste0("Downsampled Depth Distribution \n(From ",formatC(depthLen, format = "e", digits = 2)," to ",downSampleVal,") Across All Sites & Samples"),
        subtitle = paste0("SNPs Across All Samples Before Filtering: ",SNPxSamples,
                          "\nMin Depth Threshold: ", MinDPThreshold,"; Min Depth Replaces ",numSNPsRemovedMinDp," SNPs (",percentSNPsRemovedMinDp,"%)\nWith Missing a Missing Genotype (Not Including Sites With dp=0)",
        "\nMax Depth Threshold: ", MaxDPThreshold,"; Max Depth Replaces ",numSNPsRemovedMaxDp," SNPs (",percentSNPsRemovedMaxDp,"%)\nWith Missing a Missing Genotype"),
        x = "Mean Depth",y = "Count") +
        theme_minimal(base_size = 14) +
        xlim(lowerXLim, upperXLim) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        panel.grid.major = element_line(color = "grey90"))

################################################################################

        # Takes up a lot of space so remove
        rm(depthDF)
        rm(dp)

        # Since we delete sites based on mean depth, depthViolations may have some sites that aren't present in vcfSubset anymore. To avoid the computationally expensive step of calculating dp, we subset depthViolations to only include the sites remaining after mean depth filter(s).
        validSites = vcfSubset@fix[, "ID"]  # Extract site IDs currently in vcfSubset
        depthViolations = depthViolations[validSites, , drop = FALSE]  # Subset depthViolations

        # Replace only failing samples with missing data
        vcfSubset@gt[,-1][depthViolations]="./.:0:0,0,0,0"

        # Takes up a lot of space so remove
        rm(depthViolations)

        # Nothing should change after these depth filters because it replaces by sample/site, not by site or individual
        updateCounts(MinDPThreshold)
        updateCounts(MaxDPThreshold)

        # Missing sites filter
        missingSNPs = apply(extract.gt(vcfSubset), 1, function(x) all(is.na(x))) # Identify SNPs with all missing data. Applies function over all rows to check if all cols are missing for each SNPs


##################################### PLOT #####################################
              AllMissingDF=data.frame(
                Category = c("Not All Missing","All Sites Missing"),
                Count = as.numeric(table(missingSNPs))
            )

            AllMissingDF$Percent = round(AllMissingDF$Count / sum(AllMissingDF$Count) * 100, 1)
            numSNPsRemoved=AllMissingDF$Count[AllMissingDF=="All Sites Missing"]
            percentSNPsRemoved=AllMissingDF$Percent[AllMissingDF=="All Sites Missing"]

              p4 = ggplot(AllMissingDF, aes(x = "", y = Count, fill = Category)) +
              geom_bar(stat = "identity", width = 1) +
                scale_fill_manual(values = c("All Sites Missing" = "coral", "Not All Missing" = "steelblue")) +
              coord_polar(theta = "y") +
              theme_void() +
              labs(title = "Proportion of Sites That Have All Missing Data",
                  subtitle = paste0("SNPs Before Filtering: ",last(SNPsRemaining),
                  "\nFilter Removes ",numSNPsRemoved," SNPs (",percentSNPsRemoved,"%)")) +
              geom_text(aes(label = paste0(Percent, "%")),
                        position = position_stack(vjust = 0.5), size = 4, fontface="bold") +
              theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"))


################################################################################

        vcfSubset = vcfSubset[!missingSNPs, ]


        gt = extract.gt(vcfSubset, element="GT", as.numeric=TRUE)

        #Calculate missingness per site
        missingRatePerSite = rowMeans(is.na(gt))



##################################### PLOT #####################################

        lowerXLim = min(missingRatePerSite)*.99
        upperXLim = max(missingRatePerSite)*1.01

        numSNPsRemovedSiteMissingness = sum(missingRatePerSite>missingSite)+sum(missingSNPs)
        percentSNPsRemovedMissingSite = round(numSNPsRemovedSiteMissingness/last(SNPsRemaining)*100,3)

        p5 = ggplot(data.frame(missingRatePerSite), aes(x = missingRatePerSite)) +
        geom_histogram(bins = 30, fill = "steelblue", alpha = 0.5, color = "black") +
        geom_vline(xintercept = missingSite,
        col = "coral", linetype = "dashed", linewidth = 1) +
        labs(title = "Site Missingness Distribution",
        subtitle = paste0("SNPs Before Filtering: ",last(SNPsRemaining),
                          "\nSite Missingness Threshold: ", missingSite,
                          "\nSite Missingness Threshold Removes ",numSNPsRemovedSiteMissingness," SNPs (",percentSNPsRemovedMissingSite,"%)"),
        x = "Site Missingness",y = "Count") +
        theme_minimal(base_size = 14) +
        xlim(lowerXLim, upperXLim) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        panel.grid.major = element_line(color = "grey90"))

################################################################################


        #Filter for site missingness (remove sites with missing site rate over missingSitesThreshold)
        vcfSubset=vcfSubset[missingRatePerSite<=missingSite,]


        updateCounts(missingSite)


        # Missing individuals filter
        missingRatePerInd=colMeans(is.na(gt))

        # Set missingness threshold (e.g., remove individuals with >91% missing genotypes)
        indivsToKeep=colnames(gt)[missingRatePerInd <= missingIndivsThreshold]

        # Add outgroups back in if this is an outgroup included file and they happened to be cut off by missingIndivsThreshold
         if (OutgroupIncludedTF){
           indivsToKeep = union(indivsToKeep, outgroupSamples)
         }


##################################### PLOT #####################################

        lowerXLim = min(missingRatePerInd)*.99
        upperXLim = max(missingRatePerInd)*1.01

        numSamplesRemovedIndMissingness = sum(missingRatePerInd>missingIndivsThreshold)
        percentSamplesRemovedMissingInd = round(numSamplesRemovedIndMissingness/last(samplesRemaining)*100,3)

        p6 = ggplot(data.frame(missingRatePerInd), aes(x = missingRatePerInd)) +
        geom_histogram(bins = 20, fill = "steelblue", alpha = 0.5, color = "black") +
        geom_vline(xintercept = missingIndivsThreshold,
        col = "coral", linetype = "dashed", linewidth = 1) +
        labs(title = "Individual Missingness Distribution",
        subtitle = paste0("Samples Before Filtering: ",last(samplesRemaining),
                          "\nIndividual Missingness Threshold: ", missingIndivsThreshold,
                          "\nIndividual Missingness Threshold Removes ",numSamplesRemovedIndMissingness," Samples (",percentSamplesRemovedMissingInd,"%)"),
        x = "Individual Missingness",y = "Count") +
        theme_minimal(base_size = 14) +
        xlim(lowerXLim, upperXLim) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
        panel.grid.major = element_line(color = "grey90"))

################################################################################


        vcfSubset = vcfSubset[, c(TRUE, colnames(vcfSubset@gt)[-1] %in% indivsToKeep)] #Subset cols, appending TRUE to beginning to include 'FORMAT' header col

        updateCounts(missingIndivsThreshold)

        # Now filter MAC/MAF
        gt = extract.gt(vcfSubset, element="GT", as.numeric=FALSE)  #Genotypes. Redo bc filtered out many in missing sites filter. May not be necessary though.

        # Copy over because to manipulate because it can't be converted to numeric
        gtNumeric = gt
        # Their use of as.numeric=TRUE calls 1/1 as a 1 but for mac/maf we want to count this as a 2
        gtNumeric[gtNumeric == "0/0"] = 0
        gtNumeric[gtNumeric == "0/1" | gtNumeric == "1/0"] = 1
        gtNumeric[gtNumeric == "1/1"] = 2

        gtNumeric = apply(gtNumeric, c(1, 2), as.numeric)

        # MAC/MAF. Remember, these are just sums across the loci. No hetero/homozygosity is taken into account.
        mac = rowSums(gtNumeric, na.rm = TRUE)
        maf = mac/(2*ncol(gt))  # x2 because diploid

        macMafIndicesToKeep = which(mac >= MACMinThreshold & maf >= MAFMinThreshold)

##################################### PLOT #####################################

        lowerXLimMAC = min(mac)
        upperXLimMAC = round(max(mac)*1.1,0)

        lowerXLimMAF = min(maf)
        upperXLimMAF = round(max(maf)*1.1,1)

        numSitesRemovedMACMAF=sum(mac<MACMinThreshold | maf<MAFMinThreshold)
        percentSitesRemovedMACMAF=round(numSitesRemovedMACMAF/last(SNPsRemaining)*100,1)

        p7 = ggplot(data.frame(mac), aes(x = mac)) +
          geom_histogram(bins = upperXLimMAC, fill = "steelblue", alpha = 0.5, color = "black") +
          geom_vline(xintercept = MACMinThreshold, col = "coral", linetype = "dashed", linewidth = 1) +
          labs(x = "MAC", y = "Count") +
          theme_minimal(base_size = 14) +
          xlim(lowerXLimMAC, upperXLimMAC) +
          theme(plot.title = element_blank(),
                panel.grid.major = element_line(color = "grey90"))

        p8 = ggplot(data.frame(maf), aes(x = maf)) +
          geom_histogram(bins = 30, fill = "steelblue", alpha = 0.5, color = "black") +
          geom_vline(xintercept = MAFMinThreshold, col = "coral", linetype = "dashed", linewidth = 1) +
          labs(x = "MAF", y = "Count") +
          theme_minimal(base_size = 14) +
          xlim(lowerXLimMAF, upperXLimMAF) +
          theme(plot.title = element_blank(),
                panel.grid.major = element_line(color = "grey90"))

        # Create shared title and subtitle
        titleTextMACMAF = "MAC and MAF Distributions"
        subtitleTextMACMAF = paste0("Sites Before Filtering: ", last(SNPsRemaining),
                               "\nMAC Threshold: ", MACMinThreshold, "; MAF Threshold: ", MAFMinThreshold,
                               "\nMAC and MAF Thresholds Remove ", numSitesRemovedMACMAF, " Sites (", percentSitesRemovedMACMAF, "%)")

################################################################################


        # Now branch vcfSubset to have one with MAC/MAF filters and one without
        vcfSubset = list(
          WithMACMAF = vcfSubset[macMafIndicesToKeep,],
          NoMACMAF = vcfSubset
        )

        # Now branch updateCounts variables to have one with MAC/MAF filters and one without
        samplesRemaining = list(
          WithMACMAF = c(samplesRemaining,ncol(vcfSubset$WithMACMAF@gt) - 1),
          NoMACMAF = samplesRemaining
        )
        SNPsRemaining = list(
          WithMACMAF = c(SNPsRemaining,nrow(vcfSubset$WithMACMAF@fix)),
          NoMACMAF = SNPsRemaining
        )
        filterValue = list(
          WithMACMAF = c(filterValue,paste0("MAC=",MACMinThreshold,"; MAF=",MAFMinThreshold)),
          NoMACMAF = filterValue
        )
        filterName = list(
          WithMACMAF = c(filterName,"MAC and MAF"),
          NoMACMAF = filterName
        )

        # Initialize variables that will be branched
        p9=p10=p11=c(ggplot())
        filterDFList=VCFFileName=list()
        SNPsRemainingBeforeThinning=SNPsRemainingAfterThinning=SNPsRemovedByThinning=PercentSNPsRemovedByThinning=c()

        # Branch and analyze each separately
        for (branch in names(vcfSubset)) {

          if(branch == "WithMACMAF"){
            VCFFileName[[branch]] = paste0(rawString, popName, "Outgroup", OutgroupIncludedTF, "RemoveIndels", removeIndels, "ReplaceMultiallelic", replaceMultiallelic, "MissingSites", missingSite, "MissingIndiv", missingIndivsThreshold, "MAC", MACMinThreshold, "MAF", MAFMinThreshold, "MinMeanDP", MinMeanDPThreshold, "MaxMeanDP", MaxMeanDPThreshold, "MinDP", MinDPThreshold, "MaxDP", MaxDPThreshold, "Thin", thinningFilter)
          }else{
            VCFFileName[[branch]] = paste0(rawString, popName, "Outgroup", OutgroupIncludedTF, "RemoveIndels", removeIndels, "ReplaceMultiallelic", replaceMultiallelic, "MissingSites", missingSite, "MissingIndiv", missingIndivsThreshold, "MACNoneMAFNone", "MinMeanDP", MinMeanDPThreshold, "MaxMeanDP", MaxMeanDPThreshold, "MinDP", MinDPThreshold, "MaxDP", MaxDPThreshold, "Thin", thinningFilter)
        }

          # Calculate seq length and save to txt file to be accessed in future analyses (e.g., Gadma)
          chromPos = data.frame(CHROM = vcfSubset[[branch]]@fix[, "CHROM"],  POS = as.numeric(vcfSubset[[branch]]@fix[, "POS"]))

          #Find the max position per CHROM (estimated CHROM length)
          maxPosPerChrom = aggregate(POS ~ CHROM, data = chromPos, max)
          colnames(maxPosPerChrom) = c("CHROM","MaxPos")
          #Count the number of variants per CHROM
          variantCounts = aggregate(POS ~ CHROM, data = chromPos, length)
          colnames(variantCounts) = c("CHROM","VariantCount")
          chromStats = merge(maxPosPerChrom, variantCounts, by = "CHROM")

          # Compute the variant density per CHROM
          chromStats$LengthEstimate = chromStats$MaxPos * (chromStats$VariantCount+1) / chromStats$VariantCount

          # Compute the total observed sequence length (sum of estimated chrom lengths)
          seqL = round(sum(chromStats$LengthEstimate))

          seqLList = c(seqLList, as.character(seqL))
          VCFFileNameList = c(VCFFileNameList, VCFFileName[[branch]])

          # Now thin after calculating seqL
          vcfSubset[[branch]] = distance_thin(vcfSubset[[branch]], min.distance = thinningFilter)

          samplesRemaining[[branch]] = c(samplesRemaining[[branch]],ncol(vcfSubset[[branch]]@gt) - 1)
          SNPsRemaining[[branch]] = c(SNPsRemaining[[branch]],nrow(vcfSubset[[branch]]@fix))
          filterValue[[branch]] = c(filterValue[[branch]],thinningFilter)
          filterName[[branch]] = c(filterName[[branch]],"Thinning")


##################################### PLOT #####################################
          SNPsRemainingBeforeThinning[[branch]] = last(SNPsRemaining[[branch]], 2)[1]
          SNPsRemainingAfterThinning[[branch]] = last(SNPsRemaining[[branch]])
          SNPsRemovedByThinning[[branch]] = SNPsRemainingBeforeThinning[[branch]] - SNPsRemainingAfterThinning[[branch]]

          PercentSNPsRemovedByThinning[[branch]] = round((SNPsRemovedByThinning[[branch]] / SNPsRemainingBeforeThinning[[branch]]) * 100, 1)

          thinDF=data.frame(
                    Thinning = c("SNPs Remaining", "SNPs Removed"),
                    Percent = c(100-PercentSNPsRemovedByThinning[[branch]],PercentSNPsRemovedByThinning[[branch]]) )

          if (branch=="NoMACMAF"){
            branchNote="Without MAC or MAF Filters"
          }else{
            branchNote="With MAC or MAF Filters"
          }

          p9[[branch]] = ggplot(thinDF, aes(x = "", y = Percent, fill = Thinning)) +
          geom_bar(stat = "identity", width = 1) +
          scale_fill_manual(values = setNames(c("steelblue", "coral"), levels(thinDF$Thinning))) +
          coord_polar(theta = "y") +
          theme_void() +
          labs(title = "Proportion of Sites Removed by Thinning",
              subtitle = paste0("SNPs Before Filtering: ",SNPsRemainingBeforeThinning[[branch]],
              "\nFilter Removes ",SNPsRemovedByThinning[[branch]]," SNPs (",PercentSNPsRemovedByThinning[[branch]],"%)","\n",branchNote)) +
          geom_text(aes(label = paste0(Percent, "%")),
                    position = position_stack(vjust = 0.5), size = 4, fontface="bold") +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"))


################################################################################


##################################### TABLE ####################################

      filterDFList[[branch]] = data.frame(
        Step = 1:length(filterName[[branch]]),
        Filter = filterName[[branch]],
        Value = filterValue[[branch]],
        Samples = samplesRemaining[[branch]],
        SNPs = SNPsRemaining[[branch]])

      rownames(filterDFList[[branch]]) = NULL

################################################################################


##################################### PLOT #####################################

        p10[[branch]] = ggplot(filterDFList[[branch]], aes(x = factor(Filter, levels = Filter), y = Samples,group=1)) +
          geom_line(color = "steelblue", linewidth = 1.2) +
          geom_point(color = "steelblue", size = 3) +
          labs(title = "Samples Remaining", y = "Count",x="") +
          theme_minimal(base_size = 14) +
          annotate("text", x = filterDFList[[branch]]$Filter[1], y = filterDFList[[branch]]$Samples[1],
          label = filterDFList[[branch]]$Samples[1], vjust = 2.5, fontface = "bold") +
          annotate("text", x = filterDFList[[branch]]$Filter[nrow(filterDFList[[branch]])], y = filterDFList[[branch]]$Samples[nrow(filterDFList[[branch]])],
          label = filterDFList[[branch]]$Samples[nrow(filterDFList[[branch]])], vjust = -2.5, fontface = "bold") +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

        # SNPs plot
        p11[[branch]] = ggplot(filterDFList[[branch]], aes(x = factor(Filter, levels = Filter), y = SNPs,group=2)) +
          geom_line(color = "coral", linewidth = 1.2) +
          geom_point(color = "coral", size = 3) +
          labs(title = "SNPs Remaining", y = "Count",x="") +
          theme_minimal(base_size = 14) +
          annotate("text", x = filterDFList[[branch]]$Filter[1], y = filterDFList[[branch]]$SNPs[1],
          label = scientific(filterDFList[[branch]]$SNPs[1]), vjust = 0.5,hjust = -.1,fontface = "bold") +
          annotate("text", x = filterDFList[[branch]]$Filter[nrow(filterDFList[[branch]])], y = filterDFList[[branch]]$SNPs[nrow(filterDFList[[branch]])],
          label = filterDFList[[branch]]$SNPs[nrow(filterDFList[[branch]])], vjust = -2.5, fontface = "bold") +
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


################################################################################


        VCFFilePath = paste0(folderPath, VCFPath, VCFFileName[[branch]], ".vcf")
        con = file(VCFFilePath, "w") #Open file
        writeLines(vcfSubset[[branch]]@meta,con) #Write header lines to the file
        bodyDF = as.data.frame(vcfSubset[[branch]]@fix)
        colnames(bodyDF)[colnames(bodyDF) == 'CHROM'] <- '#CHROM' # Needs to have # in front of this to be identified as the header column
        #Write the VCF body to the file with no header row and no row names
        write.table(cbind(bodyDF, vcfSubset[[branch]]@gt), file = con, sep= "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        close(con)


################################### ALL PLOTS ###################################
        # Create the table grob
        tableGrobObj = tableGrob(filterDFList[[branch]], theme = ttheme_minimal(base_size = 8, base_colour = "black"), rows = NULL)
        # Ensure popsDF is a character vector
        popsText = paste(popsDF, collapse = " ")

        # Dynamically adjust number of rows for pops text display
        numPops = length(popsDF)
        numCols = 4  # Max pops displayed per row before wrapping
        numRows = ceiling(numPops / numCols)

        # Split into rows to prevent overly long text lines
        popsRows = split(pops, ceiling(seq_along(pops) / numCols))
        popsTextList = lapply(popsRows, function(x) paste(x, collapse = " "))

        # Convert to grobs for multiple rows
        popsGrobs = lapply(popsTextList, function(text) textGrob(text, gp = gpar(fontsize = 10, col = "gray40")))

        if(OutgroupIncludedTF){
          outgroupText="Outgroup Samples Included"
        }else{
          outgroupText="Outgroup Samples Not Included"
        }

        grid.arrange(
          textGrob(popName, gp = gpar(fontsize = 18, fontface = "bold")),  # Main Title
          textGrob(outgroupText, gp = gpar(fontsize = 10, fontface = "bold")),
          do.call(arrangeGrob, c(popsGrobs, ncol = 1)),  # Dynamic pops text rows
          tableGrobObj,
          nrow = 4,
          heights = c(0.1,0.1, 0.1 * length(popsGrobs), 1.5)  # Adjust heights based on rows
        )

      }

      grid.arrange(
        arrangeGrob(p1, ncol = 1),
        arrangeGrob(p2, p3, ncol = 2),
        nrow = 2
      )

      grid.arrange(
        arrangeGrob(p4, ncol = 1),
        arrangeGrob(p5, p6, ncol = 2),
        nrow = 2
      )

      grid.arrange(
        textGrob(titleTextMACMAF, gp = gpar(fontsize = 16, fontface = "bold")),
        textGrob(subtitleTextMACMAF, gp = gpar(fontsize = 12, col = "gray40")),
        arrangeGrob(p7, p8, ncol = 2),
        do.call(arrangeGrob, c(p9[names(vcfSubset)], list(ncol = 2))),
        nrow = 4,
        heights = c(0.2, 0.2, 1, 1)
      )

      grid.arrange(do.call(arrangeGrob, c(p10[names(vcfSubset)], list(ncol = 2))),
                   do.call(arrangeGrob, c(p11[names(vcfSubset)], list(ncol = 2))))

################################################################################
      # Takes up a lot of space so remove
      rm(p1 ,p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)

    }
  }
}
dev.off()


# Write set name and  estimated sequence length to text file
seqLTxtFilePath = paste0(folderPath, gadmaPath, "seqL.txt")
writeLines(c("VCFFileName\tSeqL" , paste(VCFFileNameList, seqLList, sep = "\t")), seqLTxtFilePath)

```
</details>
For each combination of filter thresholds, a series of charts and analyses are generated to help determine the most appropriate filters for a given dataset. This process is difficult to fully automate because the optimal thresholds depend on the specific goals of downstream analyses, which can vary widely. Different biological systems may also prioritize different metrics. For instance, some taxa may require maximizing the number of SNPs, even if it means relaxing depth filters. The following charts are produced for each filter set to support this decision-making process.
<br>
<br>

![Filter Triage](https://github.com/mellamoadam/CladoScope/blob/main/Images/FilterTriage.png)

## Filter Set Selection
Using the charts above, we select the best filter sets for our dataset. 
<br>
<br>
<details>
<summary>Filter Selection</summary>
<br>      
          
```r

################################# USER INPUTS  #################################
# Set the values determined to be the "best" filter set based on filter comparison pdf. This block assumes that each subset (of samples) has the same "best" conditions with or without MAC/MAF filters and with or without outgroup samples.
removeIndelsBest = FALSE
replaceMultiallelicBest = TRUE
missingSitesThresholdBest = 0.6
missingSitesThresholdBestConcatenatedTree = 0.7
missingIndivsThresholdBest = 0.92
MACMinThresholdBest = 3
MAFMinThresholdBest = 0.05
MinMeanDPThresholdBest = 6
MaxMeanDPThresholdBest = 9
MinDPThresholdBest = 5
MaxDPThresholdBest = 12
thinningFilterBest = 200

################################################################################

bestFilterSetsOutgroupOmittedNoMACMAF = c()
bestFilterSetsOutgroupIncludedNoMACMAF = c()
bestFilterSetsOutgroupOmittedWithMACMAF = c()
bestFilterSetsOutgroupIncludedWithMACMAF = c()

for (populations in seq_along(popsOfInterest)) {
  for (OutgroupIncludedTF in includeOutgroup){
      popName = names(popsOfInterest)[populations]

      FilterStringWithMACMAF = paste0(rawString,popName,"Outgroup",OutgroupIncludedTF,"RemoveIndels",removeIndelsBest,"ReplaceMultiallelic",replaceMultiallelicBest,"MissingSites",missingSitesThresholdBest,"MissingIndiv",missingIndivsThresholdBest,"MAC",MACMinThreshold,"MAF",MAFMinThreshold,"MinMeanDP",MinMeanDPThreshold,"MaxMeanDP",MaxMeanDPThreshold,"MinDP",MinDPThreshold,"MaxDP",MaxDPThreshold,"Thin",thinningFilter)

      FilterStringNoMACMAF = paste0(rawString,popName,"Outgroup",OutgroupIncludedTF,"RemoveIndels",removeIndelsBest,"ReplaceMultiallelic",replaceMultiallelicBest,"MissingSites",missingSitesThresholdBest,"MissingIndiv",missingIndivsThresholdBest,"MACNoneMAFNone","MinMeanDP",MinMeanDPThreshold,"MaxMeanDP",MaxMeanDPThreshold,"MinDP",MinDPThreshold,"MaxDP",MaxDPThreshold,"Thin",thinningFilter)
      
      if (OutgroupIncludedTF) {
        bestFilterSetsOutgroupIncludedWithMACMAF = c(bestFilterSetsOutgroupIncludedWithMACMAF,FilterStringWithMACMAF)
        bestFilterSetsOutgroupIncludedNoMACMAF = c(bestFilterSetsOutgroupIncludedNoMACMAF,FilterStringNoMACMAF)
      }else{
        bestFilterSetsOutgroupOmittedWithMACMAF = c(bestFilterSetsOutgroupOmittedWithMACMAF,FilterStringWithMACMAF)
        bestFilterSetsOutgroupOmittedNoMACMAF = c(bestFilterSetsOutgroupOmittedNoMACMAF,FilterStringNoMACMAF)
      }
      
      
    }
}

bestFilters = c(bestFilterSetsOutgroupOmittedNoMACMAF,bestFilterSetsOutgroupIncludedNoMACMAF,bestFilterSetsOutgroupOmittedWithMACMAF,bestFilterSetsOutgroupIncludedWithMACMAF)



################################# USER INPUTS  #################################
# If any additional samples are desired to be removed from BED file, use the code below to save as outgroupSamplesPath and removed in the loop
RemoveSamplesBed=c()

# Files to convert to Bed/Fam/Bim
convertToBed=bestFilters

# Files to convert to PHY/Nexus/Bin.Nexus
convertToPhy_Nex_BinNex = list.files(path=paste0(folderPath, VCFPath),pattern=paste0(rawString,"AllOutgroupTRUE.+vcf"))


# Files to convert to Traw
convertToTraw = bestFilterSetsOutgroupOmittedWithMACMAF

################################################################################

# Convert VCF to Bed/Fam/Bim
if(length(RemoveSamplesBed)>0){
  RemoveSamplesFromBedTXTPath=paste0(folderPath,admixturePath, rawString,"RemoveSamplesBed.txt")
  header="#FID IID PID MID Sex Phenotype"
  fam_data=data.frame(FID = rep(0, length(RemoveSamplesBed)),IID = RemoveSamplesBed,PID = rep(0, length(RemoveSamplesBed)),MID = rep(0, length(RemoveSamplesBed)),Sex = rep(1, length(RemoveSamplesBed)),Phenotype = rep(-9, length(RemoveSamplesBed)))
  outputLines=c(header, apply(fam_data, 1, paste, collapse = " "))
  writeLines(outputLines, RemoveSamplesFromBedTXTPath)
}


for (filterSet in convertToBed){

  currentVCFPath=paste0(folderPath,VCFPath,filterSet,".vcf")

  bedSubsetPath=paste0(folderPath,bedPath,filterSet)

    #After reformatting the CHROM col, the following requires "--chr-set", which is specified as 95 and works with all files for some reason?
  system(paste0("plink --vcf ",currentVCFPath," --make-bed --out ",bedSubsetPath," --const-fid 0 --allow-extra-chr 0"))
  
  
  # Remove samples (if applicable) from bed file, save to temp bed/fam/bim files, then overwrite bedSubsetPath
  if(length(RemoveSamplesBed)>0){
    system(paste0("plink --bfile ", bedSubsetPath," --remove ", RemoveSamplesFromBedTXTPath," --make-bed ", "--out ", paste0(bedSubsetPath,"SamplesRemoved")))
    
    # Overwrite temp bed/fam/bim files to bedSubsetPath
    system(paste0("mv ",paste0(bedSubsetPath,"SamplesRemoved",".bed")," ",paste0(bedSubsetPath,".bed")))
    system(paste0("mv ",paste0(bedSubsetPath,"SamplesRemoved",".bim")," ",paste0(bedSubsetPath,".bim")))
    system(paste0("mv ",paste0(bedSubsetPath,"SamplesRemoved",".fam")," ",paste0(bedSubsetPath,".fam")))
    
    # Remove all temp files (.log/.nosex)
    file.remove(  list.files(path=paste0(folderPath,bedPath), pattern=paste0("\\.(nosex|log)$"), all.files=TRUE, full.names=TRUE)  )
  }

}


```
</details>

## Filter Set VCF Conversion
The code below will convert the selected filter sets to the necessary file types for future analyses (.nexus, .bin.nexus, .traw, .phy).  
<br>
<details>
<summary>VCF Conversion</summary>
<br>
          
```r
# Convert VCF to PHY, Nexus, and Bin.Nexus
for (filterSet in convertToPhy_Nex_BinNex){
    subsetAndFilterString=str_sub(filterSet,end=-5)
    system(paste0("cd ",phylipCoversionFolderPath, " && ", pythonPath," vcf2phylip.py --input ",folderPath, VCFPath, filterSet," -n -b"))

    #Automatic .phy naming system a bit strange so identify the file 
    currentPhyFile=list.files(path=phylipCoversionFolderPath,pattern=paste0(subsetAndFilterString,".+phy"))
    currentNexusFile=list.files(path=phylipCoversionFolderPath, pattern = paste0(subsetAndFilterString,".+nexus"))
    currentNexusFile=currentNexusFile[!grepl("\\.bin\\.nexus$", currentNexusFile)] # Drop .bin.nexus files
    currentBinaryNexusFile=list.files(path = phylipCoversionFolderPath, pattern=paste0(subsetAndFilterString, ".+.bin.nexus"))
    
    # Move to working folder
    system(paste0("mv ", phylipCoversionFolderPath,"/", currentPhyFile," ",folderPath, PHYPath, subsetAndFilterString, ".phy"))
    system(paste0("mv ", phylipCoversionFolderPath,"/",currentNexusFile," ",folderPath, nexusPath, subsetAndFilterString, ".nexus"))
    system(paste0("mv ", phylipCoversionFolderPath,"/", currentBinaryNexusFile," ",folderPath, binaryNexusPath, subsetAndFilterString, ".bin.nexus"))

}


# Convert VCF to Traw
for (filterSet in convertToTraw){
    system(paste0("plink --vcf ",folderPath,VCFPath,filterSet,".vcf --recode A-transpose --out ",folderPath,TrawPath,filterSet," --const-fid 0 --allow-extra-chr --chr-set 95"))
    file.remove(  list.files(path=paste0(folderPath,TrawPath), pattern=paste0("\\.(nosex|log)$"), all.files=TRUE, full.names=TRUE)  )
}


```
</details>

# Clustering
<details>
<summary>PCA, k-means, and DAPC</summary>
          
```r

################################# USER INPUTS  #################################
# Set list of VCF subset file paths to perform PCA and DAPC on (generally no outgroup and MAC/MAF filters on)
PCADAPCBestFilterList = bestFilterSetsOutgroupOmittedWithMACMAF

# If any additional samples are desired to be removed
RemoveSamplesPCADAPC=c()

# Take the number of PCs that explain this proportion of total variance for feeding into LDA
cumulativeSumPCACutoff = 0.8 
      
# Upper limit of number of populations to split subset into for K-means. For example, if we think the subset is 5 populations, we can test K=2:8 with expression below.
KLowerLimitKMeans=2
calculateKMeansSeq = function(populationLength){
  KLowerLimitKMeans:round(populationLength+sqrt(populationLength)*3)
}

# The code uses K-means to find the best number of clusters. To manually set for each subset, use the code below:
manualDAPCGroupNumSelection = FALSE
DAPCGroupNums=list()
DAPCGroupNums$Subset = str_extract(PCADAPCBestFilterList, paste0("(?<=",rawString,").*?(?=Outgroup)"))

DAPCGroupNums$GroupNum = c(3,3,3,4,4,5,13,17)

# Population map to use
PCADAPCPopMap = matchDF

KMeansSeed = 1

# Remove hybrids defined in hybridSamples. Already omitted from all subsets besides the subset "All"
removeHybrids = TRUE

################################################################################

# Initialize df that will contain all DAPC assignments for all subsets
DAPCAssignments = data.frame(Subset = character(),
                           Sample = character(),
                           DAPCGroup = character()
                           )

PCADAPCPDFPath = paste0(folderPath, pdfPath, rawString, "ManualDAPCGroupNumSelection", manualDAPCGroupNumSelection, "PCADAPC.pdf")
pdf(file = PCADAPCPDFPath)

for (currentBestFilterPCA in PCADAPCBestFilterList){
      currentSubset = str_extract(currentBestFilterPCA, paste0("(?<=",rawString,").*?(?=Outgroup)"))

      # Pull in Traw file that has 0's, 1's, 2's, and NA's representing distance from reference allele
      currentTrawDF = data.frame(read.table(paste0(folderPath, TrawPath, currentBestFilterPCA,".traw") , header=TRUE))[-c(1:6)]

      
      # Remove hybrid samples if specified. Hybrids would only be present in "All" subset since matchDF has a Population "H. hybrid" 
      if (removeHybrids) {
        currentTrawDF = currentTrawDF[, !(colnames(currentTrawDF) %in% paste("X0_", hybridSamples, sep = ""))]
      }
      

      # Edit currentTrawDF to take numeric  
      currentTrawDF %>% mutate_if(is.character,as.numeric)
      
      # Edit currentTrawDF to take transpose  
      currentTrawDF = t(currentTrawDF) 

      # Edit currentTrawDF to omit all NA samples
      currentTrawDF = currentTrawDF[,colSums(is.na(currentTrawDF))<nrow(currentTrawDF)] 
      
      # Edit currentTrawDF to mean impute
      currentTrawDF = na.aggregate(currentTrawDF) 
            
      if (length(RemoveSamplesPCADAPC)>0){
        currentTrawDF = currentTrawDF[!rownames(currentTrawDF) %in% RemoveSamplesPCADAPC, ]
      }      
      
      # Center data, perform PCA, and keep all PCs
      currentPCA = dudi.pca(currentTrawDF, nf = nrow(currentTrawDF), center=TRUE, scale=FALSE, scannf=FALSE) 

      # Get rid of the "X0_" that is in the rownames after PCA
      rownames(currentPCA$li) = sub("^X0_", "", rownames(currentPCA$li))
      
      # Match row names of currentPCA$li with PCADAPCPopMap$Sample and assign Population values
      currentPCA$li$pop = PCADAPCPopMap$Population[match(rownames(currentPCA$li), PCADAPCPopMap$Sample)]
      
      currentPCACumulativeSum = cumsum(currentPCA$eig)/sum(currentPCA$eig)

      currentScreeDF = data.frame(PC= 1:currentPCA$nf, Eigenvalues = currentPCA$eig, VariationExplained = currentPCA$eig/sum(currentPCA$eig), CumulativeSumVarExplained = cumsum(currentPCA$eig)/sum(currentPCA$eig))
      
      
##################################### PLOT #####################################  

      require(gridExtra)
      xBreaksSeq=seq(min(currentScreeDF$PC), max(currentScreeDF$PC)+2, round(max(currentScreeDF$PC)/10))
      
      p1 = ggplot(currentScreeDF, aes(PC, VariationExplained)) +
        geom_col(fill = "coral1") +
        scale_x_continuous(breaks = seq(1, max(currentScreeDF$PC), by = 1)) +
        theme(axis.text=element_text(size=12)) +
        scale_x_continuous(breaks=xBreaksSeq) +
        labs(x="Principal Component", y="% Variation Explained") + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      
      p2 = ggplot(currentScreeDF,aes(x = PC, y = CumulativeSumVarExplained)) +
        geom_line(linewidth=1.2,color = "brown") +
        geom_point(color="purple") +
        scale_x_continuous(breaks = seq(1, max(currentScreeDF$PC), by = 1)) +
        theme(axis.text=element_text(size=12)) +
        scale_x_continuous(breaks=xBreaksSeq) +
        scale_y_continuous(limits = c(0, 1.1)) +
        labs(x="Principal Component", y="Cumulative % Variation Explained") + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      
      grid.arrange(p1, p2, ncol=2,
                   bottom = textGrob(paste(substr(currentBestFilterPCA, 1, nchar(currentBestFilterPCA) %/% 2), substr(currentBestFilterPCA, nchar(currentBestFilterPCA) %/% 2 + 1, nchar(currentBestFilterPCA)), sep="\n"), x = 1, hjust = 1, gp = gpar(fontsize = 7.5)))
      
      
################################################################################ 
      

##################################### PLOT #####################################  

      p3 = s.class(currentPCA$li, fac = factor(currentPCA$li$pop),
                 col = transp(funky(length(unique(currentPCA$li$pop))),alpha=0.4),
                 axesel = FALSE, cstar = 0, cpoint = 3, pch = 18, csub = 0.55,
                 sub = paste(substr(currentBestFilterPCA, 1, nchar(currentBestFilterPCA) %/% 2), substr(currentBestFilterPCA, nchar(currentBestFilterPCA) %/% 2 + 1, nchar(currentBestFilterPCA)), sep="\n"))

################################################################################ 

##################################### PLOT #####################################  
      currentPCA$li$pop = as.factor(currentPCA$li$pop)
      
      p4 = ggplot(currentPCA$li, aes(x = Axis1, y = Axis2, color = pop, shape = pop)) +
        geom_point() +
        scale_shape_manual(values = 1:length(unique(currentPCA$li$pop))) +  
        theme(axis.text = element_text(size=12))+
        labs(x = "PC1", y = "PC2")+ 
        labs(caption = paste(substr(currentBestFilterPCA, 1, nchar(currentBestFilterPCA) %/% 2), substr(currentBestFilterPCA, nchar(currentBestFilterPCA) %/% 2 + 1, nchar(currentBestFilterPCA)), sep="\n"))+
        theme(plot.caption = element_text(size = 6,hjust = 0))
      grid.arrange(p4, ncol=1)
      
      
################################################################################ 

            
##################################### PLOT #####################################  
      p5 = plot_ly(x = currentPCA$li$Axis1, y = currentPCA$li$Axis2, z = currentPCA$li$Axis3,
                   type = "scatter3d",
                   mode = "markers",
                   color = factor(currentPCA$li$pop),
                   colors = RColorBrewer::brewer.pal(min(length(unique(currentPCA$li$pop)),9),"Set1"))
      
      htmlwidgets::saveWidget(as_widget(p5), 
                              paste0(folderPath, ThreeDPath, "PCA3D", currentBestFilterPCA, ".html"))
      
################################################################################ 
      
      # K-means on PCs, find best group # by AIC, DAPC (PCA-LDA) to visualize
      numAssumedPops = length(unique(currentPCA$li$pop))
      maxNClust = last(calculateKMeansSeq(numAssumedPops))  
      
      # Number of clusters must be less than the number of samples
      if ( maxNClust > nrow(currentTrawDF) ){
        maxNClust = nrow(currentTrawDF)-1
      } 
      
      set.seed(KMeansSeed)

      currentKMeans = find.clusters(currentTrawDF, n.clust=NULL, max.n.clust=maxNClust, n.pca = nrow(currentTrawDF), stat ="AIC", choose.n.clust=FALSE, criterion = "min") #goodfit #smoothNgoesup #goesup #min

##################################### PLOT ##################################### 
      
      plot(currentKMeans$Kstat, type = "o", xlab = "Number of Clusters (K)", ylab = "AIC",
           col = "blue", pch = 16, lwd = 2, cex = 1.2, xaxt = "n", bty = "l")
      axis(1, at = 1:length(currentKMeans$Kstat), labels = 1:length(currentKMeans$Kstat))
      grid(col="gray80")
      
################################################################################ 
 
      # Overwrite find.clusters if manualDAPCGroupNumSelection == TRUE
      if (manualDAPCGroupNumSelection){
        setNClust = DAPCGroupNums$GroupNum[DAPCGroupNums$Subset == currentSubset]
        currentKMeans = find.clusters(currentTrawDF, n.clust = setNClust, max.n.clust = maxNClust, n.pca = nrow(currentTrawDF), stat ="AIC", choose.n.clust = TRUE)
      }
      
      
      cumulativeSumPCACutoffPCs = min(which(currentPCACumulativeSum >= cumulativeSumPCACutoff))
      
      currentLDA = dapc(currentTrawDF, currentKMeans$grp, n.pca = cumulativeSumPCACutoffPCs, n.da = 100)

      
##################################### PLOT ##################################### 

      scatter(currentLDA, posi.da="topright", bg="white",cstar=0, scree.pca=TRUE, posi.pca = "bottomleft")

################################################################################ 
      
##################################### PLOT ##################################### 
      
      scatter(currentLDA, scree.da=FALSE, bg="white", pch=16, cell=0, cstar=0, solid=.4, cex=1, clab=0, leg = TRUE, txt.leg = paste("Cluster", 1:length(currentKMeans$size)))
      points(currentLDA$grp.coord[,1], currentLDA$grp.coord[,2], pch=4,cex=1, lwd=2, col="black")

      
################################################################################ 
      
      
      # Make DF that includes all IDs, popmap, PCA-LDA coordinates/group assignments
      grpDF = data.frame(currentKMeans$grp)
      grpDF$Sample = row.names(grpDF)
      grpDF$Sample = sub("^X0_", "", grpDF$Sample)
      grpDF = grpDF %>%
        rename(pop = currentKMeans.grp)

      row.names(grpDF) = NULL
      
      grpDFPopmapDF = inner_join(PCADAPCPopMap, grpDF, by = "Sample")
      
      # Append DAPC groupings to DAPCAssignments df
      grpDFPopmapDFAppend = grpDFPopmapDF
      grpDFPopmapDFAppend$Subset = currentSubset
      DAPCAssignments = rbind(DAPCAssignments, grpDFPopmapDFAppend)
      
      LDACoords = data.frame(currentLDA$ind.coord)
      LDACoords$Sample = row.names(LDACoords)
      LDACoords$Sample = sub("^X0_", "", LDACoords$Sample)
      row.names(LDACoords) = NULL
      
      
      grpDFPopmapLDACoordsDF = inner_join(grpDFPopmapDF, LDACoords)

      
##################################### PLOT ##################################### 
numPops = length(unique(grpDFPopmapLDACoordsDF$Population))

ggplot(grpDFPopmapLDACoordsDF, aes(LD1, LD2)) +
    geom_point(aes(color = factor(Population), shape = factor(Population)), size = 3, alpha = 0.8) +
    scale_color_manual(values = rep(RColorBrewer::brewer.pal(9, "Set1"), length.out = numPops )) + 
    scale_shape_manual(values = 1:length(unique(grpDFPopmapLDACoordsDF$Population))) +
    theme_classic(base_size = 14) +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")) +
    labs(x = "LD1", y = "LD2")

################################################################################ 
    
      grpDFPopmapCoordsDF = inner_join(grpDFPopmapDF, coordinates, by = c("Sample" = "ID"))
      grpDFPopmapCoordsDF$pop = as.character(grpDFPopmapCoordsDF$pop)

      # Save for making map later with color coordination between different analyses
      if (currentSubset == "All") {
        grpDFPopmapCoordsDFAll = grpDFPopmapCoordsDF
      }

      # For bounds on map plot
      latitudeLowerLimit = min(grpDFPopmapCoordsDF$Latitude,na.rm = TRUE) - mapBoundsBuffer
      latitudeUpperLimit = max(grpDFPopmapCoordsDF$Latitude,na.rm = TRUE) + mapBoundsBuffer
      longitudeLowerLimit = min(grpDFPopmapCoordsDF$Longitude,na.rm = TRUE) - mapBoundsBuffer
      longitudeUpperLimit = max(grpDFPopmapCoordsDF$Longitude,na.rm = TRUE) + mapBoundsBuffer
      
        
      p66 = ggplot() +
      # Add in terrain background
      geom_raster(data = terrainDF, aes(x = Longitude, y = Latitude, fill = Color)) +
      scale_fill_identity() + 
    
      # Add in borders
      geom_sf(data = countryBordersMap, fill = NA, color = "black", linewidth = 0.3) +
      geom_sf(data = internalBordersMap, fill = NA, color = "black", linewidth = 0.3, alpha = 0.9) +
    
      geom_point(data = grpDFPopmapCoordsDF, aes(x = Longitude, y = Latitude, fill = pop), size = 3, alpha= 0.7, shape = 21, color = "black", stroke = 0.5) +
    geom_text(data = grpDFPopmapCoordsDF, aes(x = Longitude, y = Latitude, label = pop), color = "white", size = 2, fontface = "bold", vjust = 0.5, hjust = 0.5) +
    labs(title = "PCA-LDA Group Assignments with Number of Clusters\n Chosen by K-means Minimum AIC", x = "", y = "") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank(), legend.position = "none") + 
    labs(caption=paste(substr(currentBestFilterPCA, 1, nchar(currentBestFilterPCA) %/% 2), substr(currentBestFilterPCA, nchar(currentBestFilterPCA) %/% 2 + 1, nchar(currentBestFilterPCA)), sep="\n"))+
    theme(plot.caption = element_text(size = 6,hjust = 0)) +
  
        coord_sf(xlim = c(longitudeLowerLimit, longitudeUpperLimit), 
             ylim = c(latitudeLowerLimit, latitudeUpperLimit), expand = FALSE)  
      
      grid.arrange(p66, ncol=1)
      
}

dev.off()


```
</details>

<br>
<br>

<br>
<br>

![PCA DAPC Kmeans](https://github.com/mellamoadam/CladoScope/blob/main/Images/Clustering.png)
<br>
<br>
The plots above illustrate results from several analyses, including PCA scree plots, AIC curves used to determine the optimal number of k-means clusters, 3D PCA visualizations, and Discriminant Analysis of Principal Components (DAPC). For the DAPC, we retain the principal components that collectively explain 80% of the total variance, and perform linear discriminant analysis using the number of clusters corresponding to the lowest AIC value from the k-means analysis.
<br>
<br>

# ADMIXTURE
<details>
<summary>Admixture</summary>
<br>
          
```r

################################# USER INPUTS  #################################
# Set list of VCF subset file paths to perform Admixture on (generally no outgroup and MAC/MAF filters on)
AdmixtureBestFilterList = bestFilterSetsOutgroupOmittedWithMACMAF

# Population map to use for calculating sequence of K-values to run
initializeAdmixturePopMap = PopMapDAPC

# Upper limit of number of populations to split subset into. For example, if we think the subset is 5 populations, we can test K = 2:8 with expression below.
KLowerLimit = 2
calculateKSeq = function(populationLength) {
  KLowerLimit:round(populationLength + sqrt(populationLength) * 1.5)
}

# Do not append hybridSamples to samplesToRemove since we want to see these samples in admixture 
samplesToRemoveAdmixture = c(samplesToRemove)

################################################################################


for (bestFilterBySubset in AdmixtureBestFilterList){

  currentVCFPath = paste0(folderPath, VCFPath, bestFilterBySubset, ".vcf")

  # Read in VCF to get 
  currentVCFAdmixture = read.vcfR(currentVCFPath)
      
  currentVCFSamples = colnames((currentVCFAdmixture)@gt)[-1]
  filteredVCFSamples = currentVCFSamples[!(currentVCFSamples %in% samplesToRemoveAdmixture)]
      
  admixturePopMap = initializeAdmixturePopMap
  admixturePopMap = admixturePopMap[admixturePopMap$Sample %in% filteredVCFSamples,]
  
  # Grab the name of the subset populations based on it's position in the file string (between rawString and "Outgroup")
  currentSubsetName = sub(paste0("^", rawString, "(.*?)(Outgroup).*"), "\\1", bestFilterBySubset)
  
  # Find the number of populations in the current subset. Complicated because in popsOfInterest, for example the popsOfInterest chlorophaea is actually 3 subspecies that are grabbed. Also works if popsOfInterest specified subspecies.
  currentSubsetPopulationLength = length(unique(admixturePopMap$Population))
  
  clusterSeq = calculateKSeq(currentSubsetPopulationLength)
  
  # Initialize plotlist and create plotlist iteration variable
  plotlist = list()
  iteration = 1
  CVErrorDF = data.frame(K=numeric(),Subset=character(),CVError=numeric())
  
  for (K in clusterSeq){
    # Run admixture
     system(paste0("cd ", admixtureSWFolderPath, " && ", admixtureSWFolderPath, "admixture ", folderPath, bedPath, bestFilterBySubset, ".bed ", K, " --cv  | tee ", folderPath, admixturePath, "CVErrors/CVLog",currentSubsetName,"_",K,".txt"), intern = TRUE)

    Sys.sleep(1) # Let files write before naming file in the next step

    CVError = system(paste0("grep -h CV ", folderPath, admixturePath, "CVErrors/CVLog", currentSubsetName, "_", K, ".txt | awk -F': ' '{print $2}'"),intern=TRUE)
    # Add a new row to the data frame
    CVErrorDF = rbind(CVErrorDF, data.frame(K=K,Subset=currentSubsetName,CVError=as.numeric(CVError)))
    
    # Move P and Q files to Admixture folder
    PPath = paste0(folderPath,admixturePath,bestFilterBySubset,".",K,".P") 
    QPath = paste0(folderPath,admixturePath,bestFilterBySubset,".",K,".Q")
    famPath = paste0(folderPath,bedPath,bestFilterBySubset,".fam")

    PFiles=list.files(path=paste0(admixtureSWFolderPath),pattern=paste0(bestFilterBySubset,".+P"),full.names=TRUE)
    QFiles=list.files(path=paste0(admixtureSWFolderPath),pattern=paste0(bestFilterBySubset,".+Q"),full.names=TRUE)
    system(paste0("mv ",PFiles," ",PPath))
    system(paste0("mv ",QFiles," ",QPath))

    
    P = read.table(PPath, header = FALSE)
    Q = read.table(QPath, header = FALSE)
    fam = read.table(famPath, header = FALSE)
    famDF = data.frame(fam)
    names(famDF)[names(famDF) == 'V2'] = 'Sample'
    
    # Match the samples to the population
    QJoined = inner_join(famDF, matchDF,by="Sample") 
    QJoined = cbind(QJoined[c("Sample","Population")],Q)
    QJoined = QJoined[order(QJoined$Population),]
    
    Q_long = melt(QJoined, id.vars = c("Sample", "Population"), variable.name = "AncestryComponent", value.name = "Proportion")
    
    Q_long$Sample = factor(Q_long$Sample, levels = unique(QJoined$Sample))
    
    p = ggplot(Q_long, aes(x = Sample, y = Proportion, fill = AncestryComponent)) + 
      geom_bar(stat = "identity") +
      labs(x = "Sample", y = "Ancestry", fill = "Ancestry Component") +
      scale_fill_manual(values = rainbow(ncol(QJoined) - 2)) +  
      theme_minimal() +
      theme(axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank())
    
    # Hide x-axis labels until the last plot (when K == length(clusterSeq))
    if (K < length(clusterSeq)) {
      p = p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    } else {
      p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
    }

    
    p = p + guides(fill = "none")

    
    plotlist[[iteration]] = p   

    
    
    # Join each sample with coordinates
    coords=coordinates
    colnames(coords)[colnames(coords) == "ID"]="Sample"
    coords$Sample = gsub("^X0_", "", coords$Sample)

    QMap=inner_join(QJoined, coords, by="Sample")
    
    # For bounds on map plot
    latitudeLowerLimit = min(QMap$Latitude,na.rm = TRUE) - mapBoundsBuffer
    latitudeUpperLimit = max(QMap$Latitude,na.rm = TRUE) + mapBoundsBuffer
    longitudeLowerLimit = min(QMap$Longitude,na.rm = TRUE) - mapBoundsBuffer
    longitudeUpperLimit = max(QMap$Longitude,na.rm = TRUE) + mapBoundsBuffer
    

    AdmixturePieMapsPDFPath=paste0(folderPath, pdfPath, 'AdmixturePieMaps/Pie', bestFilterBySubset, "K", ncol(Q), ".pdf")
    
    pdf(file=AdmixturePieMapsPDFPath)
  
    
    QMapAdmixtureColnames = colnames(QMap)[grepl("^V\\d+$", colnames(QMap))]
admixtureColors = setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(QMapAdmixtureColnames)), QMapAdmixtureColnames)

      
    p65 = ggplot() +
    # Add in terrain background
    geom_raster(data = terrainDF, aes(x = Longitude, y = Latitude, fill = Color)) +
    scale_fill_identity() + 
  
    # Add in borders
    geom_sf(data = countryBordersMap, fill = NA, color = "black", linewidth = 0.3) +
    geom_sf(data = internalBordersMap, fill = NA, color = "black", linewidth = 0.3, alpha = 0.9) +
    
    # Add pie charts on new scale so pies don't override raster
    new_scale_fill() +
    geom_scatterpie(aes(x = Longitude, y = Latitude, r = 0.24), 
                    data = QMap, 
                    cols = QMapAdmixtureColnames, 
                    color = "black", alpha = 0.8, linewidth = 0.1) +
    scale_fill_manual(values = admixtureColors) +  
    theme_void() +  
    theme(legend.position = "none") + 
      labs(caption = paste(substr(bestFilterBySubset, 1, nchar(bestFilterBySubset) %/% 2), substr(bestFilterBySubset, nchar(bestFilterBySubset) %/% 2 + 1, nchar(bestFilterBySubset)), sep="\n")) +
      theme(plot.caption = element_text(size = 6, hjust = 0)) +
      coord_sf(xlim = c(longitudeLowerLimit, longitudeUpperLimit), 
           ylim = c(latitudeLowerLimit, latitudeUpperLimit), 
           expand = FALSE)  
    
    grid.arrange(p65, ncol=1)

    dev.off()
    
    iteration=iteration+1
  }
  
  AdmixtureBarPlotsPDFPath=paste0(folderPath,pdfPath,'AdmixtureBarPlots/AdBar',bestFilterBySubset, ".pdf")
  pdf(file=AdmixtureBarPlotsPDFPath)
  footnote=textGrob(paste(substr(bestFilterBySubset, 1, nchar(bestFilterBySubset) %/% 2), substr(bestFilterBySubset, nchar(bestFilterBySubset) %/% 2 + 1, nchar(bestFilterBySubset)), sep="\n"), gp = gpar(fontsize = 6))

  grid.arrange(grobs=plotlist,bottom = footnote)
  dev.off()

  #Combine pdf for each K value into one and delete individuals
  AdmixturePiePDFList=mixedsort(list.files(path=paste0(folderPath,pdfPath,'AdmixturePieMaps'), pattern=paste0(bestFilterBySubset,".+pdf"), all.files=TRUE, full.names=TRUE))
  pdf_combine(input = paste0(AdmixturePiePDFList), output = paste0(folderPath,pdfPath,'AdmixturePieMaps/',bestFilterBySubset,"Subset.pdf"))
  file.remove(list.files(path=paste0(folderPath,pdfPath,"AdmixturePieMaps/"), pattern=paste0("Pie",".+K",".+pdf"), all.files=TRUE, full.names=TRUE)  )
 
  AdmixtureCVErrorPlotsPDFPath=paste0(folderPath,pdfPath,'CVErrorLinePlot', ".pdf")
  pdf(file=AdmixtureCVErrorPlotsPDFPath)
  ggplot(CVErrorDF,aes(x=K,y=CVError))+geom_line()+geom_point()+facet_wrap(~Subset,scales="free_y",ncol=1)+labs(x="K",y="CV Error")+theme_minimal()
  dev.off()
}

# Now combine all subset Admixture into one pdf and remove each individual subset
AdmixturePieSubsetsPDFList=list.files(path=paste0(folderPath,pdfPath,'AdmixturePieMaps/'), pattern="Subset.pdf", all.files=TRUE, full.names=TRUE)
pdf_combine(input = paste0(AdmixturePieSubsetsPDFList), output = paste0(folderPath,pdfPath,'AdmixturePieMaps/',"AdmixturePieMaps.pdf"))
file.remove(  AdmixturePieSubsetsPDFList )

AdmixtureBarSubsetsPDFList=list.files(path=paste0(folderPath,pdfPath,'AdmixtureBarPlots/'), pattern=paste0("AdBar",rawString,".+pdf"), all.files=TRUE, full.names=TRUE)
pdf_combine(input = paste0(AdmixtureBarSubsetsPDFList), output = paste0(folderPath,pdfPath,'AdmixtureBarPlots/',"AdmixtureBarPlots.pdf"))
file.remove(  AdmixtureBarSubsetsPDFList  )

```
</details>

<br>
<br>

![ADMIXTURE](https://github.com/mellamoadam/CladoScope/blob/main/Images/ADMIXTURE.png)
<br>
<br>
The plots above is an example of the results generated from ADMIXTURE analysis. Ancestral proportions are plotted as pie charts based on sampling locations. Various ancestral component values (K) are plotted for comparison purposes.  

<br>
<br>

# Concatenated Trees
We now run IQ-tree for all VCF files that do not subset out specific populations. 
<br>
<details>
<summary>IQ-Tree</summary>
<br>
 
```r

################################# USER INPUTS  #################################
# VCFTreePaths is the list of files to make IQ-Trees with. Conditions below pull VCFs with MAC/MAF filters, include the outgroup, and are subset: All
VCFTreePaths = list.files(path = paste0(folderPath, VCFPath), pattern = paste0(rawString, "AllOutgroupTRUE.*MAC[0-9]+MAF[0-9.]+.*\\.vcf$"))


# Don't plot if VCF has less than this many sites
minSitesForTree = 100

IQTreeModel="GTR+ASC"

bootstrapNum = 1000

# Popmap to use for tip colors
IQTreePopMap = matchDF

################################################################################

for (vcfTree in VCFTreePaths){
  subsetAndFilterString = str_sub(vcfTree,end=-5)
  # Only make tree if filtering protocol included 1 or more outgroup sample
  currentSampleList = system(paste0("awk 'BEGIN{FS=\"\t\"} /^#CHROM/ {for (i=10; i<=NF; i++) print $i}' ",folderPath,VCFPath,subsetAndFilterString,".vcf"),intern = TRUE)

  # Count lines that do not start with '#' (data lines)
  numSites = sum(!grepl("^#", readLines(paste0(folderPath,VCFPath,vcfTree))))

  # Only make tree if minumum site condition is met
  if (numSites >= minSitesForTree){
    # Attempt to run iqtree2. If invariant sites are present, iqtree will output a cleaned phy file with the same name, but with the suffix ".phy.varsites.phy" instead of the original ".phy" 
    system(paste0(IQTreeSWPath," -s ", folderPath, PHYPath, subsetAndFilterString, ".phy", " -m ", IQTreeModel, " -B ", bootstrapNum," -T AUTO -st DNA -redo"))
    
    Sys.sleep(1) # Let file write before naming file in the next step

    # Overwrite the original file and redo tree with the invariant sites removed if this step is required
    if(file.exists(paste0(folderPath, PHYPath, subsetAndFilterString, ".phy.varsites.phy"))) {
      system(paste0("mv ",folderPath,PHYPath,subsetAndFilterString,".phy.varsites.phy"," ",folderPath,PHYPath,subsetAndFilterString,".phy"))
      Sys.sleep(1) #Let file write before naming file in the next step
      #Rerun iqtree with invariant sites removed .phy file
      system(paste0("/opt/homebrew/bin/iqtree2 -s ",folderPath,PHYPath,subsetAndFilterString,".phy"," -m ", IQTreeModel," -B ",bootstrapNum," -T AUTO -st DNA"))
      Sys.sleep(1) #Let file write before naming file in the next step
    }
    
    currentTreeFile=list.files(path=paste0(folderPath,PHYPath),pattern=paste0(subsetAndFilterString,".+treefile"))
    
    # Create tree if the treefile was made successfully
    if (length(currentTreeFile)!=0){
      tree = read.tree(file=paste0(folderPath,PHYPath,currentTreeFile))
    }
    currentOutgroup = outgroupSamples[outgroupSamples %in% currentSampleList]
    #Root tree if any outgroup samples are present and if the treefile was made successfully
    if (length(currentOutgroup)!=0 & length(currentTreeFile)!=0){
      tree=root(tree, outgroup = currentOutgroup, resolve.root = TRUE)
    }
    
    # Create tree pdf if the treefile was made successfully
    if (length(currentTreeFile)!=0){
      #Create popmap for tree tip labels for color coordination to population names
      treeMatchDF=data.frame(Sample=tree$tip.label) # Initialize popmap with IDs
      treeMatchDF$Population = IQTreePopMap$Population[match(treeMatchDF$Sample, IQTreePopMap$Sample)]

      uniquePops=unique(treeMatchDF$Population)
      treeColors=distinctColorPalette(length(unique(treeMatchDF$Population)))
      named_colors=setNames(treeColors, uniquePops)
      treeMatchDF$Colors=named_colors[as.character(treeMatchDF$Population)]
      
      titleSplit = paste(substr(subsetAndFilterString, 1, nchar(subsetAndFilterString) %/% 2), substr(subsetAndFilterString, nchar(subsetAndFilterString) %/% 2 + 1, nchar(subsetAndFilterString)), sep="\n")
      pdf(file=paste0(folderPath,pdfPath,'Trees/',subsetAndFilterString, ".pdf"))
      plot.phylo(tree, main = titleSplit, cex=0.2, cex.main=0.6, show.node.label = TRUE, tip.color = treeMatchDF$Colors) 
      dev.off()
    }
    # Now move all tree-related files to separate folder
    allAllTreeFiles=list.files(path=paste0(folderPath,PHYPath),pattern="tree|mldist|split|bionj|log|ckp")
    
    file.rename(from = paste0(folderPath, PHYPath, allAllTreeFiles), to = paste0(folderPath, IQTreePath, allAllTreeFiles))
  }
}


treePDFList=mixedsort(list.files(path=paste0(folderPath,pdfPath,'Trees'), pattern=".pdf", all.files=TRUE, full.names=TRUE))
allTreesPath=paste0(folderPath,pdfPath,rawString,"Trees.pdf")
pdf_combine(input=c(treePDFList),output=allTreesPath)

```
</details>

# Color Coordinated Plots
After performing PCA, DAPC, and ADMIXTURE, we take our results and create plots that color coordinate the populations based on their respective DAPC populations. Here, the optimal ADMIXTURE ancestral proportion K-value is choosen as the minimum AIC value from k-means on each respective subset.
<br>
<br>
<details>
<summary>Color coordinated plots</summary>
<br>
          
```r
# In order to match the coloration of IQ Trees populations, DAPC Groups, and Admixture plots, we make TreeColPopMap that 

################################# USER INPUTS  #################################
# The following is for color coordinating DAPC, Admixture, and IQ-Trees
# We cannot just take the population assignments from the "All" subset in DAPC since because we want to grab the population assignments from each of the subsets, then manually name the assignments of the other populations after.
popAssignmentSubsets = c("Ochrorhyncha", "Jani", "ChlorophaeaSpNov1SpNov2TorquataUnaocularusCatalinae")

# Depending on how extensive the coverage of popAssignmentSubsets is, there may be some samples that don't have a named population. Define which population map to pull from for this assignment. Note that TreeColPopMap is created using DAPCAssignments, which doesn't contain any hybridSamples because it is based on the traw file that had those samples removed before analyses.
otherPopsPopMapMatchDF = matchDFJaniSplit

################################################################################


# Create the new dataframe for color coordination
TreeColPopMap = DAPCAssignments %>%
  filter(Subset %in% popAssignmentSubsets | Subset == "All") %>%
  distinct(Sample, .keep_all = TRUE) %>%
  mutate(Population = case_when(
    Subset %in% popAssignmentSubsets ~ paste(Subset, pop, sep = "_"),  # Concatenate Subset and pop for valid subsets
    TRUE ~ "Other"  # Assign "Other" for non-valid subsets
  )) 
TreeColPopMap = TreeColPopMap[, colnames(TreeColPopMap) %in% c("Sample", "Population", "Subset")]

# Match these "Other" samples to their otherPopsPopMapMatchDF Population
TreeColPopMap$Population[TreeColPopMap$Population == "Other"] = otherPopsPopMapMatchDF$Population[match(TreeColPopMap$Sample[TreeColPopMap$Population == "Other"], otherPopsPopMapMatchDF$Sample)]

numColors = length(unique(TreeColPopMap$Population))
colorPops = distinctColorPalette(length(unique(TreeColPopMap$Population)))

# If there are more than 12 populations, generate more distinct colors
if (numColors > 12) {
  colorPops = colorRampPalette(colorPops)(numColors)
}

# Create a vector with population names as keys and colors as values
popColorMap = setNames(colorPops, unique(TreeColPopMap$Population))

# Assign colors to the Color column in TreeColPopMap based on Population
TreeColPopMap$Color = popColorMap[TreeColPopMap$Population]


# Now TreeColPopMap contains the Sample, Population, Subset used for the population assignment, and Color. For some analyses, we might want to name these populations more clearly than by "subset_pop", so we create a new column PopulationDAPC and create grpDFPopmapCoordsDFAllColorCoordination to also include associated coordinates. This population map is called PopMapDAPC.

grpDFPopmapCoordsDFAllColorCoordination = merge(grpDFPopmapCoordsDFAll, TreeColPopMap[, c("Sample", "Color", "Subset")], by = "Sample", all.x = TRUE)

# Make a color map
colorMapping = data.frame(
  Color = unique(grpDFPopmapCoordsDFAllColorCoordination$Color),
  UniqueNumber = seq_along(unique(grpDFPopmapCoordsDFAllColorCoordination$Color))
)

# Merge the mapping back into the main data frame
grpDFPopmapCoordsDFAllColorCoordination = merge(grpDFPopmapCoordsDFAllColorCoordination, colorMapping, by = "Color", all.x = TRUE)

majorityMap = grpDFPopmapCoordsDFAllColorCoordination %>%
  group_by(UniqueNumber, Population) %>%
  tally() %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

majorityMap$Population = as.character(majorityMap$Population)
dupCounts = ave(majorityMap$Population, majorityMap$Population, FUN = seq_along)
duplicatedNames = duplicated(majorityMap$Population) | duplicated(majorityMap$Population, fromLast = TRUE)
majorityMap$Population[duplicatedNames] = paste0(majorityMap$Population[duplicatedNames], "_", dupCounts[duplicatedNames])

grpDFPopmapCoordsDFAllColorCoordination = grpDFPopmapCoordsDFAllColorCoordination %>%
  left_join(majorityMap[, c("UniqueNumber", "Population")], by = "UniqueNumber") %>%
  rename(PopulationDAPC = Population.y)

PopMapDAPC = grpDFPopmapCoordsDFAllColorCoordination[, c("Sample", "PopulationDAPC")]

```
</details>

<br>
<br>

![ColorCoordinated](https://github.com/mellamoadam/CladoScope/blob/main/Images/ColorCoordinated.png)
<br>
<br>
The plots above is an example of the results generated from ADMIXTURE analysis. Ancestral proportions are plotted as pie charts based on sampling locations. Various ancestral component values (K) are plotted for comparison purposes.  



<br>
<details>
<summary>First we install and load necessary packages</summary>
<br>
          
```r


```
</details>











