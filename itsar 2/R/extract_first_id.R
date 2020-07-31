#' Extract the first string of each element containing delimiter-separated concatenated strings (in a character vector)
#'
#' @param concatstr character vector, each element of which is to be split and first element extracted from
#' @param splitchar character object to use for splitting by \code{strsplit}, default ";"
#'
#' @return character vector containing first string of each element of \code{concatstr}
#' @export
#'
#' @examples
#' library(SY351SILAC)
#' library(tidyverse)
#' psitestats$numgenes <- unlist(lapply(strsplit(psitestats$Gene.names, ";"), length))
#' table(psitestats$numgenes)
#' concatids <- psitestats %>% dplyr::filter(numgenes > 1) %>% dplyr::slice(1:20) %>% dplyr::select(Gene.names) %>% dplyr::pull()
#' cbind(extract_first_id(concatids, ";"),concatids)
#' concatids[6] <- "MYC"
#' concatids[2] <- ""
#' cbind(extract_first_id(concatids, ";"),concatids)
#'
extract_first_id <- function(concatstr,
                             splitchar = ";") {

  concat.split <- strsplit(concatstr, splitchar)

  gfirst <- unlist(lapply(concat.split, function(x) {
    if(length(x) == 0) {
      return("")
    } else {
      return(x[[1]])
    }
  }))

  gfirst[is.na(gfirst)] <- ""

  return(gfirst)

}
