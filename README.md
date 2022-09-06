# MergeION: Batch processing of LC-MS/MS data into a spectral library and several additional functions

Tandem mass spectrometry is a technique frequently used for small molecule identification. Automated structure elucidation is usually performed by spectral library search. Building a local high quality spectral library is an essentiel step thus often lacking in metabolomics and pharmaceutical laboratories. This is often due to the data confidentiality (e.g drug metadata) 

Our package fills these gaps and enables building local spectral libraries without sharing them in public domains. It works by extracting MS1 and MS2 scans from one or multiple raw chromatogram files according to m/z (and retention time) provided by users. They are then merged into a GNPS-style spectral library combining user-provided metadata. It is compatible with mzML/mzXML format acquired on Thermo, Water or Bruker instruments, in either DDA (Data-driven acquisition) or targeted MS/MS-mode.

In addition, several spectral search algorithms are available, allowing users to search an unknown spectrum in their local database or public databases (i.e. drug structures in GNPS, MASSBANK and DrugBANK).
  
# Preparation before library generation

## 1. Installation from Github in Rstudio

```R
# Install BiocManager if it has not been installed previously:
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# Install MergeION:
install.packages("remotes")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true") 
BiocManager::install("daniellyz/MergeION2")
```
# Tutorial 1: Forced degradation data
The tutorial is available: https://daniellyz.github.io/mergeion.github.io/forceddegradation.html
