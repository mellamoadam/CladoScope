First, install necessary packages listed in requirements file along with raw data.

As a starting point for the population map, popNames will be used to pull the associated string from the sample ID name.
```{r CREATE POPULATION MAP, echo = FALSE}
# Define samples to exclude due to known issues like incorrect coordinates
samplesToRemove = c("sample1_Pop1")

# Define known hybrids for downstream exclusion
hybridSamples = c("sample2_Pop2", "sample3_Pop1")

# Rename populations based on analysis insights (e.g., PCA or phylogeny)
popChanges = c(
  "sample4_Pop1" = "sample4_Pop2",
  "sample5_Pop2" = "sample5_Pop3"
)

# Define population IDs
popNames = c("Pop1",
           "Pop2",
           "Pop3"
)

```







