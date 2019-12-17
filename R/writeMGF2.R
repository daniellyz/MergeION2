#' Writing metadata and a list of spectra to mgf files
#'
#' The function writes a dataframe of metadata and a list of spectra to a mgf file
#'
#' @param library List that contains 2 elements.
#' \itemize{
#'   \item{sp:}{ List of all extracted spectra. Each spectrum is a data matrix with two columns: m/z and intensity}
#'   \item{metadata:}{ Data frame containing metadata of extracted scans. PEPMASS and RT are updated based on actually-detected scans. Following five columns are added: FILENAME, MSLEVEL, TIC, MASS_DEV, SCANNUMBER and SCANS}
#' }
#' @param con Name of the output library, the file extension must be mgf
#'
#' @examples
#'
#' data(DRUG_THERMO_LIBRARY)
#'
#' # Add new metadata "RESOLUTION = HIGH" to all scans:
#' library2$metadata$RESOLUTION = "HIGH"
#' # Write into a new mgf file:
#' writeMGF2(library2,"library_V2_bis.mgf")
#'
#' @importFrom tools file_ext
#' @export
#'
writeMGF2 <- function(library, con) {

  options(stringsAsFactors = FALSE)
  options(warn=-1)

  if (is.list(library)){
    if (length(library)==2 & "complete" %in% names(library)){
      library = library$complete
    }
    if (length(library)!=2 || (!is.list(library$sp)) || !is.data.frame(library$metadata)){
      stop("Please make sure your input library is a valid output of library_generator()!")
    }}

 if (is.character(con)){
    if (file_ext(con)!="mgf"){
      stop("The file extension of your input library must be mgf!")
    }} else {
      stop("The input must be the name of the mgf file!")}

  .cat <- function(..., file = con, sep = "", append = TRUE) {
    cat(..., file = file, sep = sep, append = append)
  }

  con <- file(description = con, open = "wt")
  on.exit(close(con))

  metadata = library$metadata
  splist = library$sp
  N=nrow(metadata)
  C=ncol(metadata)
  labels=colnames(metadata)
  for (i in 1:N) {
    .cat("\nBEGIN IONS\n")
    for (j in 1:C){
      .cat(labels[j],"=",as.character(metadata[i,j]),"\n")}
    sp=splist[[i]]
    .cat(paste(sp[,1],"\t",sp[,2], collapse = "\n"))
    .cat("\nEND IONS\n")
  }
}
