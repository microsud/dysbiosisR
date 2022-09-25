#' Dysbiosis Score Based on Median Community Level Variation
#'
#' @name dysbiosisMedianCLV
#'
#' @details Calculates median variation in a given sample compared to a reference
#'          sample group. Here \code{dysbiosisMedianCLV} will calculate median
#'          variation for a sample compared to a reference sample set. The user
#'          can provide a custom distance matrix. If the user provides a Bray-Curtis
#'          dissimilarity matrix, then the resulting score is comparable to the
#'          \code{dysbiosis score} reported in
#'          \code{Lloyd-Price J, Arze C, Ananthakrishnan AN et al. (2019)}.
#'
#' @param x A phyloseq object
#'
#' @param dist_mat A distance matrix. Can be output of \code{phyloseq::distance}
#'                 or \code{vegan::vegdist}.
#'
#' @param reference_samples Vector of samples to use as reference.
#'
#' @return A data frame with Median CLV values and sample information
#'
#' @examples
#' library(dysbiosisR)
#' # We use WirbelJ_2018 as test data
#' dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
#' ref.samples <- rownames(meta(subset_samples(WirbelJ_2018, disease == "healthy")))
#' db.1 <- dysbiosisMedianCLV(WirbelJ_2018,
#'                            dist_mat = dist.mat,
#'                            reference_samples = ref.samples)
#'
#' head(db.1)
#'
#' @references
#' \itemize{
#' \item{}{Lloyd-Price J, Arze C, Ananthakrishnan AN et al. (2019).
#' Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.
#' \emph{Nature}, 569(7758), pp.655-662.}
#'
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom stats median
#' @importFrom microbiome meta
#' @export
NULL
dysbiosisMedianCLV <- function(x = NULL,
                               dist_mat = NULL,
                               reference_samples = NULL){

  dist.mat <- Var1 <- Var2 <- score <- ids <- value <- NULL
  if(is.null(x) || is.null(dist_mat) || is.null(reference_samples)){
    stop("All arguments must be specified")
  }

  dist.mat <- dist_mat  |>
    as.matrix() |>
    as.data.frame.table(responseName = "value")

  dist.mat <- dist.mat |>
    dplyr::filter(Var1 != Var2) |>
    dplyr::filter(Var1 %in% reference_samples) |>
    dplyr::group_by(Var2) |>
    dplyr::summarise(score = median(value, na.rm = TRUE)) |>
    as.data.frame()
  colnames(dist.mat)[1] <- "ids"
  main.sam.df <- microbiome::meta(x)
  rownames(dist.mat) <- dist.mat$ids
  return(cbind(dplyr::select(dist.mat, -ids), main.sam.df))

}
