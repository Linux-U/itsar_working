#' Create model matrix for a two label isotope-labeled SILAC design, given a targets data frame
#'
#'  Modeled after limma's \code{\link[limma]{modelMatrix}} function
#'  Make sure the the target data frame has two columns: Light and Heavy
#'  which indicates the sample's SILAC label.
#'  targets. e.g. for a label-swap design:
#'  samples   Light   Heavy
#'  1    null  DMSO  DMSO
#'  2    rep1  DMSO SY351
#'  3    rep2  DMSO SY351
#'  4  rep3LF SY351  DMSO
#' @param targets matrix or data frame with columns \code{Light} and \code{Heavy} describing which samples are labeled with which isotope
#' @param parameters Not used. matrix specifying contrasts between RNA samples which should correspond to regression coefficients. Row names should correspond to unique RNA sample names found in targets.
#' @param ref character string giving name of one of the RNA sources to be treated as reference. Exactly one argument of parameters or ref should be specified.
#' @param verbose logical, if TRUE then unique names found in targets will be printed to standard output
#'
#' @return produces a numeric design matrix with row names as in targets and column names as in parameters.
#' @export
#'
#' @examples
#' targets <- data.frame(samples = c('null','rep1','rep2','rep3LF'),
#'                       Light = c("DMSO","DMSO","DMSO","SY351"),
#'                       Heavy  = c("DMSO","SY351","SY351","DMSO"),
#'                       samplenames = c("Ratio.H.L.normalized.null",
#'                                      "Ratio.H.L.normalized.rep1",
#'                                      "Ratio.H.L.normalized.rep2",
#'                                      "Ratio.H.L.normalized.rep3LF"))
#' design <- silac_model_matrix(targets, ref = "DMSO")
#' design
#' design.isotopeffect <- cbind(rep(1,length(targets$samples)), design)
#' colnames(design.isotopeffect) <- c("isotope","SY351")
#' design.isotopeffect
#'
silac_model_matrix <- function (targets, parameters = NULL, ref = NULL, verbose = TRUE) {

  if (missing(targets))
    stop("targets is required argument")
  targets <- as.matrix(targets)
  if (!all(c("Light", "Heavy") %in% colnames(targets)))
    stop("targets should contain columns: Light and Heavy")
  if (is.null(parameters) == is.null(ref))
    stop("exactly one of the arguments parameters and ref should be specified")
  target.names <- sort(unique(as.vector(t(as.matrix(targets[,
                                                            c("Light", "Heavy")])))))
  if (verbose)
    cat("Found unique target names:\n", target.names,
        "\n")
  if (is.null(parameters)) {
    if (!(ref %in% target.names))
      stop(paste("\"", ref, "\" not among the target names found",
                 sep = ""))
    other.names <- setdiff(target.names, ref)
    target.names <- c(ref, other.names)
    ntargets <- length(target.names)
    parameters <- rbind(-1, diag(ntargets - 1))
    rownames(parameters) <- target.names
    colnames(parameters) <- other.names
  }
  else {
    parameters <- as.matrix(parameters)
    if (length(target.names) != nrow(parameters))
      stop("rows of parameters don't match unique target names")
    if (any(sort(target.names) != sort(rownames(parameters))))
      stop("rownames of parameters don't match unique target names")
    target.names <- rownames(parameters)
    ntargets <- nrow(parameters)
    if (ncol(parameters) != ntargets - 1)
      warning("number of parameters should be one less than number of targets")
  }
  narrays <- nrow(targets)
  J <- matrix(rep(target.names, narrays), ntargets, narrays)
  J <- t((t(J) == targets[, "Heavy"]) - (t(J) == targets[,
                                                       "Light"]))
  rownames(J) <- target.names
  colnames(J) <- rownames(targets)
  zapsmall(t(solve(crossprod(parameters), crossprod(parameters,
                                                    J))), 14)
}
