#' Read protein quantification table from MaxQuant output file proteinGroups.txt
#' @description
#' Read protein quantification table from a file derived from the MaxQuant output file proteinGroups.txt,
#' which may have been annotated using the Perseus program.
#' one row per protein, with samples represented in columns as specified in the
#' MaxQuant experimental design file summary.table. Be careful with SILAC  - FBS
#' contamination will result in unlabled peptides identical to human orthologs and
#' will inflate the light channel,and this function filters contaminants out. Note
#' that this function (read.table actually) will use make.names during the import of
#' the proteinGroups.txt file and will convert any column names in that file by using
#' make.names, and will for example convert spaces to "." in the resulting prot.table
#' @param proteinGroups.filename file name string, assuming it is located in current directory. Default "proteinGroups.txt"
#' @param datacolnames character array of the data column names to extract, which will ignore the datacolumnprefix parameter, optional
#' @param datacolumnprefix string used to identify the column name containing the quantification. Default "LFQ intensity "
#' @param uniprot.acc.column column with the uniprot accessions (protein identifiers), default "Majority protein IDs"
#' @param uniprot.id.column column with the uniprot entry name, (previously called uniprot IDs, default "UniProt names"
#' @param meta.colnms character array of additional columns to include
#' @param null.zeros logical, If TRUE, convert zero intensities to NA. Default TRUE
#' @param log2.transform logical, log2 transform the intensities? Default TRUE
#' @param is.gzip.file boolean indicating whether the file is gzipped (*.gz)
#' @return a list with the following elements: prot.table (the data table) and
#'     prot.colnms (the data column names indexed with datacolumnprefix parameter)
#' @examples
#' read_maxquant_proteingroupfile()
#' @export
read_maxquant_proteingroupfile <- function(proteinGroups.filename = "proteinGroups.txt",
                               datacolnames = NULL,
                               datacolumnprefix = "LFQ intensity ",
                               uniprot.acc.column = "Majority protein IDs",
                               uniprot.id.column = "UniProt names",
                               meta.colnms = c(
                                 "Protein.names", "Gene.names", "uniprot.ids",
                                 "uniprot.accs", "Isoform.ID", "Sequence.lengths"
                               ),
                               null.zeros = TRUE,
                               #  sumtable,
                               log2.transform = T,
                               is.gzip.file = F) {

  extract_first_id <- function(concatstr, splitchar = ";") {
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

  #  browser()
  #  groups <- sumtable$Experiment
  #  groups <- groups[groups != ""]

  protcols <- unlist(strsplit(readLines(proteinGroups.filename, n = 1), "\t"))

  if(is.gzip.file) {
    prot <- read.table(
      file = gzfile(proteinGroups.filename),
      header = T,
      sep = "\t",
      na.strings = c("NA","NaN","Infinite"),
      fill = T,
      stringsAsFactors = F,
      quote = "\""
    )

  } else {
    prot <- read.table(
      file = proteinGroups.filename,
      header = T,
      sep = "\t",
      na.strings = c("NA","NaN","Infinite"),
      fill = T,
      stringsAsFactors = F,
      quote = "\""
    )
  }

  if (!uniprot.acc.column %in% protcols) {
    stop(paste("\nError: column name ", uniprot.acc.column, " not in names of prot data frame:\n",
      paste(protcols, collapse = ";"), "\n",
      sep = ""
    ))
  }

  if (!uniprot.id.column %in% protcols) {
    stop(paste("Error: column name ", uniprot.id.column, " not in names of prot data frame:\n",
      paste(protcols, collapse = ";"), "\n",
      sep = ""
    ))
  }

  if (is.null(datacolnames) & datacolumnprefix != "") {
    intcolindices <- grep(datacolumnprefix, protcols, value = F)

    if (length(intcolindices) < 1) {
      stop(paste("Error: specified data column prefix string: ",
        datacolumnprefix,
        "\n\tnot found in proteinGroups.txt file columns:\n",
        protcols,
        collapse = ";"
      ),
      sep = ""
      )
    }
  } else if (length(datacolnames) > 1) {
    intcolindices <- which(protcols %in% datacolnames)

    if (length(intcolindices) != length(datacolnames)) {
      stop(paste("Error: specified column names in parameter datacolnames: ",
        paste(datacolnames, collapse = ";"),
        "\n\tnot found in proteinGroups.txt file:\n",
        paste(protcols, collapse = ";"),
        sep = ""
      ))
    }
  } else {
    stop("Error: Either datacolumnprefix or datacolnames must be specified correctly.")
  }



  #  contam <- prot[grep("CON__",prot[,uniprot.acc.column],invert = F),]

  names(prot)[intcolindices] <- gsub(datacolumnprefix, "", protcols[intcolindices])

  names(prot)[intcolindices] <- make.names(names(prot)[intcolindices])

  newdatacolnms <- names(prot[, intcolindices])

  if (null.zeros) {
    prot[, newdatacolnms][prot[, newdatacolnms] == 0] <- NA
    prot$numvalid <- apply(prot[, newdatacolnms], 1, function(x) sum(!is.na(x)))
  }



  uniprot.id.column <- make.names(uniprot.id.column)
  uniprot.acc.column <- make.names(uniprot.acc.column)


  prot <- prot[grep("REV__", prot[, uniprot.acc.column], invert = T), ]
  prot <- prot[grep("CON__", prot[, uniprot.acc.column], invert = T), ]

  #prot$Only.identified.by.site <- NULL
  #prot$Potential.contaminant <- NULL
  #prot$Reverse <- NULL


  if (log2.transform) {
    prot[, newdatacolnms] <- log2(prot[, newdatacolnms])
    prot[, newdatacolnms][is.na(prot[, newdatacolnms])] <- NA
  }

  prot$uniprot.id <- extract_first_id(prot[, uniprot.id.column])
  prot$uniprot.ids <- prot[, uniprot.id.column]
  prot[, uniprot.id.column] <- NULL

  prot$uniprot.acc <- extract_first_id(prot[, uniprot.acc.column])
  prot$uniprot.accs <- prot[, uniprot.acc.column]
  prot[, uniprot.acc.column] <- NULL

  prot$Gene.name <- extract_first_id(prot$Gene.name)

  nmsnew <- c("Gene.name", "uniprot.id", "uniprot.acc", newdatacolnms, meta.colnms)
  nmsnew <- nmsnew[nmsnew %in% names(prot)]
  nmsnew <- c(nmsnew, names(prot)[!names(prot) %in% nmsnew])
  prot <- prot[, nmsnew]

  list(
    prot.table = prot,
    prot.colnms = newdatacolnms
  )
}
