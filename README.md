
# MergeION: Batch processing of LC-MS/MS data into a spectral library and several additional functions

Tandem mass spectrometry is a technique frequently used for small molecule identification. Automated structure elucidation is usually performed by spectral library search. Building a local high quality spectral library is an essentiel step thus often lacking in metabolomics and pharmaceutical laboratories. This is often due to the data confidentiality (e.g drug metadata) 

Our package fills these gaps and enables building local spectral libraries without sharing them in public domains. It works by extracting MS1 and MS2 scans from one or multiple raw chromatogram files according to m/z (and retention time) provided by users. They are then merged into a GNPS-style spectral library combining user-provided metadata. It is compatible with mzML/mzXML format converted from Thermo, Water or Bruker data files, in either DDA (Data-driven acquisition) or targeted MS/MS-mode.

Several library search algorithms are available, allowing users to search and annotate an unknown spectrum in their local database or public databases (i.e. drug structures in GNPS, MASSBANK and DrugBANK) - a quick start guide can be found below. Other functionalities such as local spectral library building, library lookup, analog search, and molecular networking can be found in the three walkthrough examples.
  
## Installation from Github in Rstudio

For MAC users, please first install XQuartz: https://www.xquartz.org/

For WINDOWS users, please use the newest 64 bit version of R, and first make sure JAVA is installed: https://www.java.com/download/ie_manual.jsp

```R
Sys.setenv(JAVA_HOME="C:/Program Files/Java/...)"
library(rJava)
```
Please check this forum post if you encounter errors: https://support.microsoft.com/en-us/topic/qa-when-i-try-to-load-the-rjava-package-using-the-library-command-i-get-an-error-531cb2e1-6ee1-5f5f-e4cf-40b819b5aaa3

```R
# Install BiocManager if it has not been installed previously:
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# Install MergeION:
options(timeout=9999999)
install.packages("remotes")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true") 
remotes::install_github("daqana/dqshiny", force = T)
remotes::install_github("schymane/RChemMass", force = T)
BiocManager::install("daniellyz/MergeION2")
```

## Quick start guide for compound annotation

We have pre-compiled a small molecule spectral database containing MS/MS spectra of 11,642 metabolites, natural products and drugs. The database combines public repositories such as GNPS, MassBank and reference spectra of in-house standards of approved drugs. Currently ESI-MS/MS spectra in our collection are all in positive ion mode. We are constantly populating this spectral database by combining various repositories after careful inspection of spectra and metadata quality. Users can annotate an unknown MS/MS spectrum with confidence by searching this spectral library.

This library should be first downloaded to your R working directory and loaded to your R environment:

```R
download.file("https://zenodo.org/record/7057435/files/GNPS_MASSBANK_PROCESSED_POS_CONSENSUS1.RData?download=1", "GNPS_MASSBANK_PROCESSED_POS_CONSENSUS1.RData")
load("GNPS_MASSBANK_PROCESSED_POS_CONSENSUS1.RData")
```

Now the library is loaded as _library1c_ in your R environment, now we can use this library to search and annotate an experimental spectrum. The MS/MS spectrum should be read as well into the R environment as a two-column matrix. In this example, we know the query spectrum corresponds to Cinnarizine. We use this example only for method validation:

```R
query.sp = read.csv("https://raw.githubusercontent.com/daniellyz/MESSAR/master/MESSAR_WEBSERVER_DEMO/example_cinnarizine.txt", header = F, sep = "\t")
head(query.sp)

```

Now it's time to collect query parameters into a R list. Don't be overwhelmed by the long list. Only important parameters to check are the _prec_mz_, which indicates the precursor mass, and _use_prec_, which forces precursor mass match in the search output by setting to TRUE (the classical approach for compound annotation):

```R
params.query.sp = list(prec_mz = 369.232, use_prec = T, polarity = "Positive", method = "Cosine", min_frag_match = 6, min_score = 0, reaction_type = "Metabolic")

```

One command-line to run spectral library search:

```R
search_result = library_query(input_library = library1c, query_spectrum = query.sp, params.query.sp = params.query.sp)

```

We can now print the candidate structure(s) found. In this example only one, and the Cosine spectral similarity is very high at 0.95: 

```R

print(search_result$consensus$metadata)

```

The INCHIKEY indicates we found the correct hit (Cinnarizine), we can inspect the spectra match by mirror plot:

```R
id_matched = search_result$consensus$metadata$ID
library_visualizer(library1c, id = id_matched, query_spectrum = query.sp)
```

![Alt text](https://github.com/daniellyz/MergeION2/raw/master/inst/mirror.png "Mirror Cinnarizine")

## Calling a GUI for compound annotation

Alternatively, you could perform spectral library by calling a webtool after successful installation of MergeION:

```R
library(MergeION) # First load MergeION
library(RChemMass) # Load the RChemMass package for structure visualization
runGUI()
```
A tutorial for the GUI is available at: https://github.com/daniellyz/mergeion.github.io/blob/gh-pages/Webtool_documentation.pdf

# Tutorial 1: Building a local database

The tutorial is available at: https://daniellyz.github.io/mergeion.github.io/Library_Generation.html

# Tutorial 2: Forced degradation data analysis

The tutorial is available at: https://daniellyz.github.io/mergeion.github.io/forceddegradation.html

# Tutorial 3: Urinary metabolomics data analysis

The tutorial is available at: https://daniellyz.github.io/mergeion.github.io/urinarymetabolomics.html

