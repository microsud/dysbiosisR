#' Dysbiosis Score Gradient Visualization
#'
#' @name plotDysbiosisGradient
#'
#' @details A wrapper for dysbiosis gradient visualization.
#'
#' @param df The data frame output from different dysbiosis calculators provided
#'           by dysbiosisR package.
#'
#' @param score Variable to plot on x-axis.
#'
#' @param colors colors to use for plotting. Default is NULL.
#'
#' @param show_points Logical TRUE or FALSE.
#'
#' @return A ggplot2 object
#'
#' @examples
#' data("WirbelJ_2018")
#' library(RColorBrewer)
#' dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
#' # get reference samples
#' ref.samples <- sample_names(subset_samples(WirbelJ_2018,
#'                                            disease == "healthy"))
#'
#' dysbiosis_1 <- dysbiosisMedianCLV(WirbelJ_2018,
#'                                   dist_mat = dist.mat,
#'                                   reference_samples = ref.samples)
#'
#' # get dysbiosis and normobiosis thresholds
#' dysbiosis_thres <- quantile(subset(dysbiosis_1, disease == "CRC")$score, 0.9)
#' normobiosis_thres <- quantile(subset(dysbiosis_1, disease == "CRC")$score, 0.1)
#'
#' plotDysbiosisGradient(df=dysbiosis_1,
#'                       score="score",
#'                       high_line = dysbiosis_thres,
#'                       low_line = normobiosis_thres,
#'                       group_var = "disease",
#'                       group_colors=c("healthy" = "steelblue", "CRC"= "brown3"),
#'                       point_size = 2,
#'                       bg_colors = rev(brewer.pal(9, "YlOrBr")))
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid rasterGrob unit
#' @importFrom grDevices colorRampPalette
#'
#' @author Sudarshan A. Shetty
#'
#' @export
NULL
plotDysbiosisGradient <- function(df=NULL,
                                  score=NULL,
                                  high_line = NULL,
                                  low_line = NULL,
                                  group_var = NULL,
                                  group_colors=NULL,
                                  point_size = 2,
                                  bg_colors = NULL){

  if(!score %in% colnames(df)){
    stop(paste0("The input data does not have column '", score, "'"))
  }

  if(!group_var %in% colnames(df)){
    stop(paste0("The input data does not have column '", group_var, "'"))
  }
  if(is.null(df) || is.null(score) || is.null(group_colors)){
    stop("All arguments must be specified")
  }

  df$sam.ind.x <- "samples"

  if(!is.null(bg_colors)){
    g <-  .make_gradient(deg = 90, n = 500, cols = bg_colors)
  } else{
    g <-  .make_gradient(deg = 90, n = 500, cols = rev(brewer.pal(9, "YlOrBr")))
  }

  p.res <- ggplot2::ggplot(df,
                           ggplot2::aes_string("sam.ind.x", "score")) +
    ggplot2::annotation_custom(
      grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    ggplot2::geom_hline(yintercept = high_line, lty ="dashed", color= "white")+
    ggplot2::geom_hline(yintercept = low_line, lty ="dashed") +
    ggplot2::geom_jitter(color="white",
                         ggplot2::aes_string(fill=group_var), shape=21, size=3,
                         position = ggplot2::position_jitter(seed = 42)) +
    ggplot2::scale_fill_manual(values = group_colors)+
    ggplot2::theme_minimal()+
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::labs(x="")
  return(p.res)

}

# https://stackoverflow.com/questions/30136725/plot-background-colour-in-gradient
.make_gradient <- function(deg = 90, n = 100, cols = NULL) {
  cols <- grDevices::colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(data = rep(seq(0, 1, length.out = n) * cos(rad), n),
                byrow = TRUE,
                ncol = n) +
    matrix(data = rep(seq(0, 1, length.out = n) * sin(rad), n),
           byrow = FALSE,
           ncol = n)
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(image = mat,
                   width = grid::unit(1, "npc"),
                   height = grid::unit(1, "npc"),
                   interpolate = TRUE)
}
