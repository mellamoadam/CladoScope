This repository provides CladoScope, a modular and iterative workflow for analyzing population structure using SNP data. It is designed for population genetic studies where group definitions are uncertain, and where initial assumptions about population membership may need refinement based on downstream analyses. This code will also run several analyses based on the final defined subsets, as well as prepare input text files necessary for several analyses. A full example with user defined parameters selected is shown in Hypsiglena.Rmd. The code is based in R but some bash commands are executed using system().

Inputs:
  1. VCF file
  2. loci file (ipyrad output; only necessary for BPP analysis)
  3. matchDF, a population map estimate for samples and their corresponding assumed population
  4. Map files for regions of interest (.shp file often used for borders; .tiff file often used for terrain)

Outputs:
  1. Converted files for subsequent analyses (.phy, .nexus, .bin.nexus, .traw, )
  2. Filter triage pdf (evaluation of various filtering schemes)
  3. Concatenated phylogenetic trees for all filter sets (for subset containing all samples)
  4. SVD Quartets phylogenetic species tree (for subset containing all samples)
  5. 2D and 3D PCA plots
  6. DAPC (PCA-LDA) plots with group ID number plotted on map
  7. Admixture pie chart maps pdf (for each K value, map with samples displayed as pie charts with ancestral proportions displayed)
  8. Color coordinated plots pdf that synchronizes color coordination between DAPC group IDs (for subset containing all samples), samples in concatenated tree, and ADMIXTURE pie chart maps and PCA plots with ADMIXTURE ancestral proportions plotted.
  9. Analaysis output files for isolation by distance and environment (collapses environmental variables into top PCA components)
  10. Analaysis output files for Gadma
  11. Analaysis output files for BPP
  12. Analaysis output files for DSuite (including heatmap grid of gene flow with Dtrios)

First, we use matchDF in the CREATE POPULATION MAP section based on our best guess for a population map. Then in the section VCF FILTERING, we filter each subset separately, but some of these population assignments might not be accurate. This allows us to define the best filter sets to use in DAPC, which gives us group assignments for each sample in their respective subset. Once that is done, we can rerun the CREATE POPULATION MAP section with a more accurate matchDF that we manually update, and then run VCF FILTERING, and DAPC. This is necessary because running DAPC on an entire dataset might not give meaningful results if there is a large amount of diversity. 

Simply assuming the subset a sample belongs to based on geography or morphological characteristics is antithetical to the whole reason for doing genetic analyses, however, it can give us a good starting point (just to determine decent filtering schemes and initial DAPC assignments). For example, we may have an idea of population assignments based on previous studies, geographical barriers, and morphology.

Workflow Summary: 
  Step 1. Make population assignment assumptions based on prior knowledge and save in section CREATE POPULATION MAP (data frame matchDF).
  Step 2. Run VCF FILTERING section and use output PDF to decide best filters, which are selected in section FILTER SELECTION.
  Step 3. Using "best" filter set, run IQ-TREE section. In this section, it is color coordinated with matchDF population assignments. 
  Step 4. This can be used to identify if any samples fall outside of their assumed population assignment. Update matchDF and hybridSamples accordingly.
  Step 5. Perform VCF FILTERING, FILTER SELECTION, DAPC etc. again with the updated subset definitions.
  Step 6. After running ADMIXTURE, update hybridSamples if neccesary.
