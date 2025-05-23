First, install necessary packages listed in requirements file along with raw data.

```{r CREATE POPULATION MAP, echo = FALSE}

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

janiGroups = list(
  "H_jani_texana_Pop1" = c("DBS_868_jani_texana", "DBS_900_jani_texana", "LJV_10635_jani_texana","MHP_10783_jani_texana", "MHP_12366_jani_texana", "MHP_13388_jani_texana","MHP_8260_jani_texana", "MF_21713_jani_jani", "MF_21732_jani_dunklei"),
  "H_jani_texana_Pop2" = c("AMNH_R_177719_jani_texana", "AMNH_R_504524_jani_texana", "AMNH_R_504527_jani_texana","CAS_228933_jani_texana", "CAS_228936_jani_texana", "CAS_228952_jani_texana","CAS_228965_jani_texana", "CAS_228965_jani_texana", "CAS_228966_jani_texana","CAS_228933_jani_texana", "CAS_228936_jani_texana", "CAS_228952_jani_texana","CAS_229281_jani_texana", "LSUMZ_84790_jani_texana","MSB_87826_jani_texana","UTAR_52351_jani_texana","MVZ_226235_spnov1xtexana"),
  "H_jani_texana_Pop3" = c("MVZ_236398_jani_texana", "ROM_5256_jani_texana", "ROM_15100_jani_texana","UTAR_52350F_jani_texana", "ROM_15094_jani_texana", "AMNH_R_504773_jani_texana","AMNH_R_177724_jani_texana", "AMNH_R_177715_jani_texana", "AMNH_R_177717_jani_texana","AMNH_R_177718_jani_texana", "AMNH_R_177720_jani_texana", "AMNH_R_504400_jani_texana","AMNH_R_504604_jani_texana", "AMNH_R_504615_jani_texana", "AMNH_R_504616_jani_texana","AMNH_R_504630_jani_texana", "AMNH_R_504631_jani_texana", "AMNH_R_504632_jani_texana","AMNH_R_504653_jani_texana","CAS_228960_jani_texana","CAS_228962_jani_texana","CAS_229229_jani_texana","CAS_229231_jani_texana","CAS_229281_jani_texana","CAS_229920_jani_texana","AMNH_R_504773_jani_texana","CAS_228960_jani_texana","CAS_228962_jani_texana","CAS_228965_jani_texana","CAS_228966_jani_texana","CAS_229229_jani_texana","CAS_229231_jani_texana","CAS_229920_jani_texana","CR_409137_jani_texana","GDC_13619_jani_texana","GDC_13715_jani_texana","MSB_87874_jani_texana","MSB_87876_jani_texana","MVZ_236398_jani_texana","ROM_15094_jani_texana","ROM_15100_jani_texana","ROM_5256_jani_texana","SM_CO2B_jani_texana","UTAR_34835_jani_texana","UTAR_51983_jani_texana","UTAR_52350F_jani_texana","UTEP_14082_jani_texana","UTEP_16307_jani_texana","UTEP_16309_jani_texana","UTEP_16317_jani_texana","UTEP_18484_jani_texana"),
   "H_jani_Other" = c("AMNH_R_504774_jani_jani","GDC_25928_jani_texana","MF_21713_jani_jani","MF_21732_jani_dunklei","UTEP_18438_jani_texana")
)
  

# Create population map matchDF
# All populations (besides spnov1xtexana because it interferes with the for loops below since grepl finds spnov1 twice) for matching with raw data and future grouping in plots
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


# Make matchDFUpdatedJaniSplit that splits jani texana into multiple populations 
matchDFJaniSplit = matchDF

for (pop in names(janiGroups)) {
  matchDFJaniSplit$Population[matchDFJaniSplit$Sample %in% janiGroups[[pop]]] = pop
}

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
