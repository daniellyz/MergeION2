% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_lib2network.R
\name{process_lib2network}
\alias{process_lib2network}
\title{Create feature-based molecular networks from a conesnsus spectra library}
\usage{
process_lib2network(
  input_library,
  networking = T,
  polarity = c("Positive", "Negative"),
  params.screening = list(baseline = 1000, relative = 0.01, max_peaks = 200),
  params.search = list(mz_search = 0.005, ppm_search = 10),
  params.similarity = list(method = "Cosine", min.frag.match = 6, min.score = 0.6),
  params.network = list(topK = 10, reaction.type = "Metabolic", use.reaction = TRUE)
)
}
\arguments{
\item{input_library}{A list object. Must contain consensus library.}

\item{networking}{Boolean. TRUE to go through the entire molecular networking process. FALSE if only spectra alignment is performed, necessary for spectral library searching}

\item{polarity}{character. Either "Positive" or "Negative". Ion mode of the LC-MS/MS file.}

\item{params.screening}{Parameters for feature screening and MS2 spectra pre-processing from raw_data_file:
\itemize{
 \item{baseline:}{ Numeric. Absolute intensity threshold that is considered as a mass peak and written into the library.}
 \item{relative:}{ Numeric between 0 and 100. Relative intensity threshold of the highest peak in each spectrum, peaks above both absolute and relative thresholds are saved in the library.}
 \item{max_peaks:}{ Integer higher than 3. Maximum number of peaks kept per spectrum from the highest peak.}
}}

\item{params.search}{List of parameters for feature screening, combining precursor ions and fragments in the input file, as well as for searching in the spectral library. The list must contain following elements:
\itemize{
 \item{mz_search:}{ Numeric. Absolute mass tolerance in Da.}
 \item{ppm_search:}{ Numeric. Absolute mass tolerance in ppm.} 
}}

\item{params.network}{Parameters for post-filtering and annotation of network edges: based on feature correlation (if feature quantification table is provided) and mass difference
\itemize{
 \item{topK:}{ Integer higher than 0. For networking, the edge between two nodes are kept only if both nodes are within each other's TopK most similar nodes. For example, if this value is set at 20, then a single node may be connected to up to 20 other nodes. Keeping this value low makes very large networks (many nodes) much easier to visualize. We suggest keeping this value at 10.}
 \item{reaction.type:}{ Character. Either "Metabolic" and "Chemical". Type of transformation list used to annotate mass difference between connected features in molecular network.}
 \item{use.reaction:}{ Boolean. TRUE if keep only edges whose mass difference can be annotated to known metabolic or chemical reactions.}
}}
}
\description{
Function used by library_generator to create molecular networks
}