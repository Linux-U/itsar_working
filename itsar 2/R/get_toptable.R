#' Table of statistics from an eBayes differential analysis object
#' @description
#' Extract a table containing the test statistics, p-values, adjusted p-values, etc. from
#' a linear model fit produced by lmFit and eBayes from the limma package. A suffix can be
#' provided to append to the column names for the statistics, to allow for downstream workflows.
#' @param fitobj list containing a linear model fit produced by lmFit and eBayes.
#'  fit should be an object of class MArrayLM as produced by lmFit and eBayes.
#' @param coef name of coefficient in \code{fitobj}
#' @param glist dataframe of annotation columns from the original data used to create \code{fitobj}
#' @param suffix suffix to append to the statistic columns
#' @return A data frame derived from a call to limma's topTable() function, with s2.prior and sigma from the lmFit and eBayes calls.
#' @importFrom magrittr %>%
#' @importFrom limma topTable
#' @importFrom tibble tibble
#' @importFrom dplyr as_tibble rename
#' @export
#' @examples
#'
get_toptable <- function(fitobj, coef, glist, suffix) {

  ttemp <- limma::topTable(fitobj, coef = coef,
                           genelist = glist,
                           number = Inf, p.value = 1,
                           confint = T,sort.by =  "none")

  ttemp$s2.post <-  fitobj$s2.post
  ttemp$sigma <- fitobj$sigma
  ttemp <- dplyr::as_tibble(ttemp)

  if( any(names(ttemp) == "Unique.identifier")) {
    ttemp <- ttemp %>% dplyr::rename( ID = Unique.identifier)
  }


  datinds <- which(!names(ttemp) %in% c("id", "ID"))

  names(ttemp)[datinds] <- paste(names(ttemp)[datinds], suffix, sep = "")
  ttemp
}

