#' Loads protein-level meta data from a Perseus annotated phosphosite MaxQuant file.
#'
#' @description
#' This assumes that the file \code{pmetafile} has been previously generated using
#' Perseus by annotating a maxQuant PhosphoSite file (Phospho(STY).txt usually), with protein-level
#' and phosphosite-level annotations. Remove the data columns to avoid duplicate columns in downstream
#' merges with data column containing data frames or tibbles. The file must have the following columns
#'  'Unique identifier': key to merge with phosphosite data table.
#'  'Interpro.name': This and next two columns are used to annotate whether the protein is a kinase
#'  'Prosite.name'
#'  'Pfam'
#'
#' @param pmetafile character, name of protein meta data file
#' @param absfilepath character, file path
#' @param is_gzip_file logical, is the file gzipped?, default TRUE
#'
#' @return data frame (\code{tibble}) with the protein meta data and an additional column \code{proteinclass} annotating whether it is a kinase
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr as_tibble mutate if_else
#' @export
#'
#' @examples
#'
get_protein_meta_data <- function(pmetafile,
                                  absfilepath,
                                  is_gzip_file = T) {


  require(magrittr)
  require(dplyr)

  pmfile <- ""
  if (is.null(absfilepath)) {
    pmfile <- pmetafile # assume it is a full path if absfilepath null
  } else {
    pmfile <- file.path(absfilepath,
                        pmetafile)
  }


  if (is_gzip_file) {
    pmfile <- gzfile(pmfile)
  }

  pmetadf <- read.table(pmfile,
                      header = T,
                      fill = T,
                      quote = "\"",
                      stringsAsFactors = F,
                      na.strings = "NaN",
                      sep = "\t",
                      comment.char = "#")


  pmetadf <- dplyr::as_tibble(pmetadf)


  pmetadf %>% dplyr::mutate(proteinclass = dplyr::if_else( (grepl("kinase",Interpro.name, ignore.case = T) |
                                                grepl("kinase",Prosite.name, ignore.case = T)  |
                                                grepl("kinase",Pfam, ignore.case = T)), "kinase", "" ))
}
