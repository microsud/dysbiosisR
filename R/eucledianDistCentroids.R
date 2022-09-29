#' Dysbiosis Score Based on Euclidean Distance to Group Centroids
#'
#' @name euclideanDistCentroids
#'
#' @details Calculates difference in \code{euclidean distance (ED)} for a sample
#'          to group centroids. For example, sample_1 to control centroid minus
#'          sample_1 to case centroid. The user can provide a custom distance
#'          matrix. This approach was used in \code{AlShawaqfeh MK et al. (2017)}.
#'
#' @param x A phyloseq object
#'
#' @param dist_mat A distance matrix. Can be output of \code{phyloseq::distance}
#'                 or \code{vegan::vegdist}.
#'
#' @param use_squared Logical. Default is FALSE. If TURE to the score is
#'                    calculated using the squared distance to group centroids.
#'                    see \code{usedist::dist_to_centroids}.
#'
#' @param group_col A column in \code{phyloseq::sample_data} with all control and
#'                  case labels.
#'
#' @param control_label A character string specifying control/healthy labels in
#'                      group_col.
#'
#' @param case_label A character string specifying case/disease labels in
#'                   group_col.
#'
#' @return A data frame with Centroid distance to each group and a score with
#'         sample information
#'
#' @examples
#' library(dysbiosisR)
#' # We use WirbelJ_2018 as test data
#' dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
#' db.1 <- euclideanDistCentroids(WirbelJ_2018,
#'                                dist_mat = dist.mat,
#'                                use_squared = TRUE,
#'                                group_col = "disease",
#'                                control_label = "healthy",
#'                                case_label = "CRC")
#'
#' head(db.1)
#'
#' @references
#' \itemize{
#' \item{}{AlShawaqfeh MK et al. (2017).
#' A dysbiosis index to assess microbial changes in fecal samples of dogs with
#' chronic inflammatory enteropathy.
#' \emph{FEMS microbiology ecology}, 93(11), p.fix136.}
#'
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom usedist dist_to_centroids
#' @importFrom microbiome meta
#' @export
NULL
euclideanDistCentroids <- function(x = NULL,
                                   dist_mat = NULL,
                                   use_squared = FALSE,
                                   group_col = NULL,
                                   control_label = NULL,
                                   case_label = NULL){

  main.sam.df <- cent.df <- cent.df.wide <- score <- CentroidDist <- Item <- NULL
  CentroidDist_Item <- CentroidDistance <- CentroidGroup <- NULL
  # main.sam.df <- cent.df <- cent.df.wide <- score <- CentroidDist <- Item <- NULL
  if(is.null(x) ||
     is.null(dist_mat) ||
     is.null(group_col) ||
     is.null(control_label) ||
     is.null(case_label)){
    stop("All arguments must be specified")
  }

  # check
  if(!group_col %in% colnames(meta(x))){
    stop(paste0("The specified Group column '", group_col, "' not found"))
  }

  # check
  if(!control_label %in% meta(x)[,group_col]){
    stop(paste0("The specified Control label '", control_label, "' not found"))
  }

  # check
  if(!case_label %in% meta(x)[,group_col]){
    stop(paste0("The specified Case label '", case_label, "' not found"))
  }

  main.sam.df <- meta(x)
  cent.df <- usedist::dist_to_centroids(dist_mat,
                                        main.sam.df[,group_col],
                                        squared = use_squared)
  cent.df.wide <- cent.df |>
    tidyr::pivot_wider(id_cols = Item,
                       names_from = CentroidGroup,
                       values_from = CentroidDistance) |>
    as.data.frame()

  cent.df.wide$score <- (cent.df.wide[,control_label]) - (cent.df.wide[,case_label])
  rownames(cent.df.wide) <- cent.df.wide$Item
  colnames(cent.df.wide) <- paste("CentroidDist", colnames(cent.df.wide), sep = "_")
  return(cbind(dplyr::select(cent.df.wide, -CentroidDist_Item), main.sam.df))
}
