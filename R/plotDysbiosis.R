#' Plot Dysbiosis
#'
#' @name plotDysbiosis
#'
#' @details A wrapper for quick plot.
#'
#' @param df The data frame output from different dysbiosis calculators provided
#'           by dysbiosisR package.
#'
#' @param xvar Variable to plot on x-axis.
#'
#' @param yvar Variable to plot on y-axis.
#'
#' @param colors colors to use for plotting. Default is NULL.
#'
#' @param show_points Logical TRUE or FALSE.
#'
#' @return A ggplot2 object
#'
#' @examples
#' data("WirbelJ_2018")
#' # data are relative abundances summed to 100 or
#' # near 100 val
#' dysbiosis.oob <- dysbiosisOBB(WirbelJ_2018,
#'                               group_col = "disease",
#'                               control_label = "healthy",
#'                               case_label = "CRC",
#'                               seed_value = 1235,
#'                               add_tuneRF_params = list(ntreeTry=100,
#'                                                        stepFactor=1.5,
#'                                                        improve=0.01,
#'                                                        trace=TRUE,
#'                                                        dobest=FALSE),
#'                               ntree = 100, # increase for real data
#'                               plot_roc = TRUE)
#' plotDysbiosis(df=dysbiosis.oob,
#'               xvar="disease",
#'               yvar="oob.score",
#'               colors=c(CRC="brown3", healthy="steelblue"))
#' @author Sudarshan A. Shetty
#'
#' @export
NULL
plotDysbiosis <- function(df=NULL,
                          xvar=NULL,
                          yvar=NULL,
                          colors=NULL,
                          show_points= TRUE){


  if(!xvar %in% colnames(df)){
    stop(paste0("The input data does not have column '", xvar, "'"))
  }

  if(!yvar %in% colnames(df)){
    stop(paste0("The input data does not have column '", yvar, "'"))
  }

  if(is.null(colors)){
    pgrp <- df |>
      ggplot2::ggplot(ggplot2::aes_string(x=xvar, y=yvar, group=xvar)) +
      # ggplot2::geom_boxplot(width = 0.1,
      #                       alpha = 0.1,
      #                       outlier.shape = NA) +
      ggdist::stat_halfeye(color = "black",
                           size=1,
                           adjust = 0.5,
                           #justification = -0.1
                           )
    if(show_points){
      pgrp <- pgrp + ggplot2::geom_jitter(shape = 21,
                                          size = 2,
                                          alpha = 0.6,
                                          width = 0.05)
    }

    pgrp <- pgrp + ggplot2::theme(legend.position = "none") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())


  } else {
    pgrp <- df |>
      ggplot2::ggplot(ggplot2::aes_string(x=xvar,
                                          y=yvar,
                                          fill = xvar,
                                          group=xvar)) +
      # ggplot2::geom_boxplot(width = 0.1,
      #                       alpha = 0.1,
      #                       outlier.shape = NA) +
      ggdist::stat_halfeye(color = "black",
                           size=1,
                           adjust = 0.5,
                           #justification = -0.1
                           )

    if(show_points){
      pgrp <- pgrp + ggplot2::geom_jitter(shape = 21,
                                          size = 2,
                                          alpha = 0.6,
                                          width = 0.05)
    }
    pgrp <- pgrp +
      ggplot2::scale_color_manual(xvar,values = colors) +
      ggplot2::scale_fill_manual(xvar, values = colors) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  return(pgrp)
}
