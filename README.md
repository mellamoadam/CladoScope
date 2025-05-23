This repository provides CladoScope, a modular and iterative workflow for analyzing population structure using SNP data. It is designed for population genetic studies where group definitions are uncertain, and where initial assumptions about population membership may need refinement based on downstream analyses. This code will also run several analyses based on the final defined subsets, as well as prepare input text files necessary for several analyses (Gadma, BPP, ADMIXTURE etc.).

Inputs:
  1. VCF file
  2. loci file (ipyrad output; only necessary for BPP analysis)
  3. matchDF, a population map estimate for samples and their corresponding assumed population

First, we use matchDF in the CREATE POPULATION MAP section based on our best guess for a population map. Then in the section VCF FILTERING, we filter each subset separately, but some of these population assignments might not be accurate. This allows us to define the best filter sets to use in DAPC, which gives us group assignments for each sample in their respective subset. Once that is done, we can rerun the CREATE POPULATION MAP section with a more accurate matchDF that we manually update, and then run VCF FILTERING, and DAPC. This is necessary because running DAPC on an entire dataset might not give meaningful results if there is a large amount of diversity. 

Simply assuming the subset a sample belongs to based on geography or morphological characteristics is antithetical to the whole reason for doing genetic analyses, however, it can give us a good starting point (just to determine decent filtering schemes and initial DAPC assignments). For example, we may have an idea of population assignments based on previous studies, geographical barriers, and morphology.

Workflow Summary: 
  Step 1. Make population assignment assumptions based on prior knowledge and save in section CREATE POPULATION MAP (data frame matchDF).
  Step 2. Run VCF FILTERING section and use output PDF to decide best filters, which are selected in section FILTER SELECTION.
  Step 3. Using "best" filter set, run IQ-TREE section. In this section, it is color coordinated with matchDF population assignments. 
  Step 4. This can be used to identify if any samples fall outside of their assumed population assignment. Update matchDF and hybridSamples accordingly.
  Step 5. Perform VCF FILTERING, FILTER SELECTION, DAPC etc. again with the updated subset definitions.
  Step 6. After running ADMIXTURE, update hybridSamples if neccesary.
