#' Volcano plot of differential proteomic data with log2 ratios, p-values and adjusted p-values.
#' @description
#' \code{ggvolcano} returns a list with two volcano plot objects (ggplot2 and plotly).
#' The data frame (parameter tdf) contains the log2 ratios (xvar), p-values (yvar),
#' adjusted p-values (adjpvalvar), labeling variable (labelvar). The columns of tdf that
#' correspond to each of these variables are provided to the ggvolcano function as unquoted
#' column names (in the tidyverse fasion) for the corresponding variables in the tdf data frame.
#' @param tdf A data frame containing the log2 ratios (xvar), p-values (yvar), labeling variable (labelvar)
#' @param xvar x-axis bare variable name of numeric column in tdf data frame, usually log2 ratio
#' @param yvar y-axis bare variable name of numeric column in tdf data frame, usually log10(unadjusted p-value)
#' @param labelvar column in tdf (bare variable name) to use for labeling individual points
#' @param adjpvalvar column in tdf (bare variable name) to use for adjusted p-value
#' @param fdr.range.cut numeric vector, FDR intervals to indicate significance and for mapping color, sorted from smallest (zero) to largest (1).
#' @param high.sig.threshold numeric, all points with adj.p.value <= to this are annotated with text specified by labelvar
#' @param xlimits length-2 numeric, x-axis limits
#' @param ylimits length-2 numeric, y-axis limits
#' @param xtitle character, x-axis title
#' @param ytitle character, y-axis title
#' @param plottitle character, plot title
#' @param annotgenes logical, should genes with adjusted p-value < high.sig.threshold be explicitly labeled
#' @param revcols logical, reverse the color vector
#' @return a named list with two plot objects, both of which can be individually plotted by assigning
#'     to a new variable and typing variable name at command line:
#'       ggplot_volcano: a ggplot version of volcano plot
#'       plotly_volcano: an interactive plotly volcano plot.
#' @importFrom ggplot2 ggplot geom_point geom_vline geom_text scale_colour_manual theme_classic theme xlim ylim xlab ylab ggtitle aes_string element_blank
#' @importFrom rlang enquo quo_name
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr pull select filter pull
#' @importFrom plotly ggplotly add_annotations
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @export
#' @examples
#' library(SY351SILAC) # loads the psitestats data frame
#' library(tidyverse)
#' psitestats <- psitestats %>% dplyr::mutate(logp = -log10(P.Value.nullfiltisotope))
#'
#' mygglist <- ggvolcano(psitestats, xvar =  logFC.nullfiltisotope,
#'                       yvar = logp,
#'                       fdr.range.cut = c(0,0.01,1),
#'                       adjpvalvar = adj.P.Val.nullfiltisotope,
#'                       labelvar = gene.psite,
#'                       high.sig.threshold = 0.001,
#'                       xlimits = c(-3.5,3.5),
#'                       ylimits = c(0,8),
#'                       annotgenes = T,
#'                       xtitle = 'log2(SY351/DMSO)',
#'                       ytitle = 'log10(p-value)')
#'
#'
#' mygglist$ggplot_volcano
#'
#' mygglist$plotly_volcano
#'
ggvolcano <- function(tdf, # data frame containing the data to plot
                 xvar, # string name of x-variable
                 yvar,# string name of y-variable
                 labelvar, # string name of label variable, usually Gene.name
                 adjpvalvar,# string name of adjusted p-value to color points by significance
                 fdr.range.cut = c(0,0.01,0.05,0.1,1), # cuts used to categorize adj.p.vals
                 high.sig.threshold = 0.1, # all points with adj.p.value <= to this are annotated with gene symbols that must be named Gene.name,
                 xlimits = Inf,
                 ylimits = Inf,
                 xtitle = NULL,
                 ytitle = NULL,
                 plottitle = NULL,
                 annotgenes = T,
                 revcols = T) {

  if(length(fdr.range.cut) < 3) {
    stop("fdr.range.cut should be a numeric vector of monotonically increasing fdr cut values of length between 3 and 5.")
  }

  fdrcut <- function(adjpval, fdr.range.cut =  c(0,0.02,0.05,1)) {

    adjpval.ranges <- cut(adjpval, fdr.range.cut, include.lowest = T)
    fdr.range.cut <- sprintf("%.2f", fdr.range.cut)


    if(length(fdr.range.cut) == 3) {
      adjpval.ranges <- factor(adjpval.ranges,levels = rev(levels(adjpval.ranges)),ordered = T,
                               labels = rev( c( paste(fdr.range.cut[1],'< q <=',fdr.range.cut[2],sep = ''),
                                                paste(fdr.range.cut[2],'< q <=',fdr.range.cut[3],sep = ''))))
    } else if(length(fdr.range.cut) == 4) {
      adjpval.ranges <- factor(adjpval.ranges,levels = rev(levels(adjpval.ranges)),ordered = T,
                               labels = rev( c( paste(fdr.range.cut[1],'< q <=',fdr.range.cut[2],sep = ''),
                                                paste(fdr.range.cut[2],'< q <=',fdr.range.cut[3],sep = ''),
                                                paste(fdr.range.cut[3],'< q <=',fdr.range.cut[4],sep = ''))))
    } else if(length(fdr.range.cut) == 5) {
      adjpval.ranges <- factor(adjpval.ranges,levels = rev(levels(adjpval.ranges)),ordered = T,
                               labels = rev( c( paste(fdr.range.cut[1],'< q <=',fdr.range.cut[2],sep = ''),
                                                paste(fdr.range.cut[2],'< q <=',fdr.range.cut[3],sep = ''),
                                                paste(fdr.range.cut[3],'< q <=',fdr.range.cut[4],sep = ''),
                                                paste(fdr.range.cut[4],'< q <=',fdr.range.cut[5],sep = ''))))

    } else {
      stop("fdr.range.cut vector not of length between 3 and 5.")
    }
    return(adjpval.ranges)
  }
#  browser()
  xvar <- rlang::enquo(xvar)
  yvar <- rlang::enquo(yvar)
  adjpvalvar <- rlang::enquo(adjpvalvar)
  labelvar <- rlang::enquo(labelvar)

  adjpvalvarcut <- fdrcut(tdf %>% dplyr::select(!! adjpvalvar) %>% dplyr::pull(),
                          fdr.range.cut = fdr.range.cut)
  txtlabels <- tdf %>% dplyr::select(!! labelvar) %>% dplyr::pull()
  adjpvals <- tdf %>% dplyr::select(!! adjpvalvar) %>% dplyr::pull()
  txtlabels[adjpvals > high.sig.threshold] <- ""

  if(revcols) {
    rcols <- rev(RColorBrewer::brewer.pal(n = 7, name = "Reds"))[1:(length(fdr.range.cut) - 2)]
  } else {
    rcols <- RColorBrewer::brewer.pal(n = 7, name = "Reds")[1:(length(fdr.range.cut) - 2)]
  }

  gcol <- RColorBrewer::brewer.pal(n = 5, name = "Greys")[3]

  ggy <- tdf %>% ggplot2::ggplot(mapping = ggplot2::aes_string( rlang::quo_name(xvar),
                                                       rlang::quo_name(yvar))) +
    ggplot2::geom_point(ggplot2::aes(color = adjpvalvarcut), alpha = 0.5) +
                #   alpha = 0.2)) +  #text = quo_name(labelvar), #Gene.name,
    ggplot2::scale_colour_manual(values = c(gcol, rev(rcols))) +
    ggplot2::geom_vline(xintercept = 0,color = "grey") +
    ggplot2::geom_text(label = txtlabels, hjust = "outward", size = 3,
              check_overlap = T) +

    ggplot2::theme_classic(base_size = 15,
                           base_line_size = 2)  +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position= c(0.85,0.2))

  # xlab(xtitle) + ylab(ytitle)

  if ( !is.infinite(xlimits[1])) {
    ggy <- ggy + ggplot2::xlim(xlimits)
  }
  if ( !is.infinite(ylimits[1])) {
    ggy <- ggy + ggplot2::ylim(ylimits)
  }
  if(!is.null(xtitle) ) {
    ggy <- ggy + ggplot2::xlab(xtitle)
  }

  if(!is.null(ytitle) ) {
    ggy <- ggy + ggplot2::ylab(ytitle)
  }

  if(!is.null(plottitle) ) {
    ggy <- ggy +  ggplot2::ggtitle(plottitle)
  }

#  tdf.thresh10pct <- tdf %>% dplyr::filter(!! adjpvalvar <= 0.01 & !! adjpvalvar > high.sig.threshold)
  ggplotly <- tdf %>% ggplot2::ggplot(mapping = ggplot2::aes_string( rlang::quo_name(xvar),
                                                                rlang::quo_name(yvar))) +
    ggplot2::geom_point(ggplot2::aes(color = adjpvalvarcut), alpha = 0.5) +
    #   alpha = 0.2)) +  #text = quo_name(labelvar), #Gene.name,
    ggplot2::scale_colour_manual(values = c(gcol, rev(rcols))) +
    ggplot2::geom_vline(xintercept = 0,color = "grey") +

    ggplot2::theme_classic(base_size = 15,
                           base_line_size = 2)  +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position= c(0.85,0.2))

  # xlab(xtitle) + ylab(ytitle)

  if ( !is.infinite(xlimits[1])) {
    ggplotly <- ggplotly + ggplot2::xlim(xlimits)
  }
  if ( !is.infinite(ylimits[1])) {
    ggplotly <- ggplotly + ggplot2::ylim(ylimits)
  }
  if(!is.null(xtitle) ) {
    ggplotly <- ggplotly + ggplot2::xlab(xtitle)
  }

  if(!is.null(ytitle) ) {
    ggplotly <- ggplotly + ggplot2::ylab(ytitle)
  }

  if(!is.null(plottitle) ) {
    ggplotly <- ggplotly +  ggplot2::ggtitle(plottitle)
  }

  tdf_thresh_highsig <- tdf %>% dplyr::filter(!! adjpvalvar <= high.sig.threshold)

#  tdf.threshy <- tdf %>% dplyr::filter(!! yvar > 25 & !! xvar < 2 )

  ggplotly <- plotly::ggplotly(ggplotly, tooltip = "text")

  # ggply <- ggply %>% plotly::add_annotations(x = tdf.thresh10pct %>% dplyr::select(!! xvar) %>% dplyr::pull(),
  #                                            y  = tdf.thresh10pct %>% dplyr::select(!! yvar) %>% dplyr::pull(),
  #                                            text = tdf.thresh10pct %>% dplyr::select(!! labelvar) %>% dplyr::pull(),
  #                                            font = list(family = 'Arial',
  #                                                        size = 12,
  #                                                        color = 'black'),
  #                                            showarrow = T,
  #                                            visible = F,
  #                                            arrowhead = 4,
  #                                            arrowwidth = 0.5,
  #                                            arrowcolor = 'grey',
  #                                            arrowsize = .5,
  #                                            ax = 20,
  #                                            ay = -20,
  #                                            clicktoshow = "onoff")

  if (annotgenes) {
    ggplotly <- ggplotly %>% plotly::add_annotations(x = tdf_thresh_highsig %>% dplyr::select(!! xvar) %>% dplyr::pull(),
                                              y  = tdf_thresh_highsig %>% dplyr::select(!! yvar) %>% dplyr::pull(),
                                              text = tdf_thresh_highsig %>% dplyr::select(!! labelvar) %>% dplyr::pull(),
                                              font = list(family = 'Arial',
                                                          size = 13,
                                                          color = 'black'),
                                              showarrow = T,
                                              visible = T,
                                              arrowhead = 4,
                                              arrowwidth = 0.5,
                                              arrowcolor = 'grey',
                                              arrowsize = .5,
                                              ax = 20,
                                              ay = -20,
                                              clicktoshow = "onoff")
  }


  return(list(ggplot_volcano = ggy,
              plotly_volcano = ggplotly))


}

# Try:
# myggma <- tophits %>% ggmaplot(xvar = abcam, yvar = avgabcam, adjpvalvar = abcam.adj.P.Val)
# ggplotly(myggma, tooltip = "text")
