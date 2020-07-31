## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  fig.width=7, fig.height=7,
  fig.align = "center"
  )
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----loading-libs, message=FALSE----------------------------------------------
library(itsar)
library(tidyverse)
require(magrittr)
require(dplyr)
require(limma)
require(sva)

options(stringsAsFactors = FALSE)


## -----------------------------------------------------------------------------
#Read proteinGroups file
proteinGroupFile <- system.file("extdata",
                                       "proteinGroups_itsa_staurosporine_3temp_20190402KB.txt.gz",
                                       package = "itsar")
# read the first line of file to retrieve the column names, for properly setting up the sample names
# gzfile() is used because we are reading a gzipped tab delimited file
# 
pgfilecolnames <- as.character(
   read.table(
      file = gzfile(proteinGroupFile),
      header = F,
      nrows = 1,
      sep = "\t",
      na.strings = c("NA","NaN","Infinite"),
      fill = T,
      stringsAsFactors = F,
      quote = "\"")[1,])
 
### Filter data columns & define TMT data columns
# identifies column numbers that contain Label corrected intensity data (TMT lot-specific correction)
tmtdata_colnames <- grep("Reporter intensity corrected \\d+ 52C", pgfilecolnames, value = T)


## -----------------------------------------------------------------------------

plist <- itsar::read_maxquant_proteingroupfile(
  proteinGroups.filename = proteinGroupFile,
  datacolnames = tmtdata_colnames,
  uniprot.id.column = "Majority protein IDs",
  uniprot.acc.column = "Fasta headers",
  null.zeros = T,
  log2.transform = F,
  is.gzip.file = T)

prot.table <- plist$prot.table
tmtdata_colnames <- plist$prot.colnms

tmtdata_colnames
prot.table <- as_tibble(prot.table)


## -----------------------------------------------------------------------------
# Remove poorly quantified proteins with less than 2 unique peptides
prot.table <- dplyr::filter(prot.table, Reverse != "+", 
                            Potential.contaminant != "+", 
                            Only.identified.by.site != "+",
                             Reporter.intensity.count.3.52C > 2,
                            Unique.peptides > 1)




## -----------------------------------------------------------------------------

# Are there any duplicated gene names?

prot.table$Gene.names %>% duplicated() %>% any()

prot.table %>% dplyr::group_by(Gene.names) %>% 
  dplyr::summarize(frequency = n()) %>% 
  dplyr::arrange(desc(frequency)) %>% 
  dplyr::filter(frequency > 1)

#prot.table <- DEP::make_unique(prot.table, "Gene.names", "Protein.IDs", delim = ";")

# prot.table %>% dplyr::group_by(name) %>% 
#   dplyr::summarize(frequency = n()) %>% 
#   dplyr::arrange(desc(frequency)) %>% 
#   dplyr::filter(frequency > 1)


## -----------------------------------------------------------------------------
limma::plotDensities(prot.table[,tmtdata_colnames], legend = "topright")


## -----------------------------------------------------------------------------
prot.table[,tmtdata_colnames] <- log2(prot.table[,tmtdata_colnames] )

limma::plotDensities(prot.table[,tmtdata_colnames], legend = "topright")

prot.table[,tmtdata_colnames] <- limma::normalizeBetweenArrays(as.matrix(prot.table[,tmtdata_colnames]),
                                                        method="quantile" )



## -----------------------------------------------------------------------------


targets <- data.frame(samples = tmtdata_colnames,
                      treatment = factor(rep(c("Vehicle","Drug"), 5), levels = c("Vehicle","Drug")))
targets

origdesign <- model.matrix(~0+treatment, data = targets)
colnames(origdesign) <- c("Vehicle","Drug")  #V=Vehicle; D=Drug

origdesign

design <- model.matrix(~treatment, data = targets)
colnames(design) <- c("intercept","Drug") 



## -----------------------------------------------------------------------------
### Use limma for empirical bayes


fit <- limma::lmFit(prot.table[,tmtdata_colnames], origdesign)

fit <- limma::eBayes(fit) ##Apply empirical Bayes smoothing to the standard errors.

contrast <- limma::makeContrasts("Drug-Vehicle", levels=origdesign)

fit2 <- limma::contrasts.fit(fit, contrast)
fit2$Amean <- log2(prot.table$Intensity.52C)

fit2 <- limma::eBayes(fit2, trend = TRUE)

## -----------------------------------------------------------------------------
pstats <- itsar::get_toptable(fit2,
                              coef = 1, 
                              glist = prot.table[,c("id", "Gene.names", "uniprot.id")], 
                              suffix = ".uncor" )

pstats <- pstats %>% arrange(P.Value.uncor)

sum(pstats$P.Value.uncor <= 0.05)
sum(pstats$P.Value.uncor <= 0.01)
sum(pstats$P.Value.uncor<= 0.001)

pstats


## -----------------------------------------------------------------------------

mod0 <- model.matrix(~1, data = targets)
colnames(mod0) <- "intercept"
dmat <- as.matrix(prot.table[,tmtdata_colnames])

numsv <- sva::num.sv(dmat, design, method = "be", B = 100)


svobj <- sva::sva(dmat,design,mod0,n.sv = 3) 

modSv = cbind(design,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = sva::f.pvalue(dmat,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
sum(qValuesSv < 0.1, na.rm = T)


svafit <- lmFit(dmat , modSv)
contrast.matrix <- cbind("intercept" = c(1,0, rep(0,svobj$n.sv) ), "Drug"= c(0,1,rep(0,svobj$n.sv)))


fitContrasts = limma::contrasts.fit(svafit,contrast.matrix)

svafit$Amean <- log2(prot.table$Intensity.52C)

fit2_sva <- eBayes(fitContrasts, trend = TRUE)

## -----------------------------------------------------------------------------
pstats_sva <- itsar::get_toptable(fit2_sva, coef = "Drug", 
                                  glist = prot.table[,c("id", "Gene.names", "uniprot.id")], 
                                  suffix = ".sva" )

pstats_sva <- pstats_sva %>% arrange(P.Value.sva)

sum(pstats_sva$adj.P.Val.sva <= 0.05)
sum(pstats_sva$adj.P.Val.sva <= 0.01)
sum(pstats_sva$adj.P.Val.sva<= 0.001)

pstats_sva


## -----------------------------------------------------------------------------

pstats$logp <- -log10(pstats$P.Value.uncor)
pstats_sva$logp <- -log10(pstats_sva$P.Value.sva)
range(pstats$logp)
range(pstats$logFC.uncor)

range(pstats_sva$logp)
range(pstats_sva$logFC.sva)

volcanoplots_uncor <- itsar::ggvolcano(pstats, xvar = logFC.uncor,
                   yvar = logp,
                   labelvar = Gene.names.uncor,
                   adjpvalvar = adj.P.Val.uncor,
                   fdr.range.cut = c(0,0.01,0.1,1),
                   high.sig.threshold = 0.001,
                   annotgenes = T, xlimits = c(-2, 2), ylimits = c(0, 16))


volcanoplots_sva <- itsar::ggvolcano(pstats_sva, xvar = logFC.sva,
                   yvar = logp,
                   labelvar = Gene.names.sva,
                   adjpvalvar = adj.P.Val.sva,
                   fdr.range.cut = c(0,0.01,0.1,1),
                   high.sig.threshold = 0.001,
                   annotgenes = T, xlimits = c(-2, 2), ylimits = c(0, 16))




## ----fig.width = 5, fig.asp = 1.1---------------------------------------------

volcanoplots_uncor$ggplot_volcano


## ----fig.width = 5, fig.asp = 1.1---------------------------------------------

volcanoplots_sva$ggplot_volcano


## ----fig.width = 5, fig.asp = 1.2---------------------------------------------

volcanoplots_uncor$plotly_volcano


## ----fig.width = 5, fig.asp = 1.2---------------------------------------------

volcanoplots_sva$plotly_volcano


