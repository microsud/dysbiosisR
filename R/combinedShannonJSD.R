#' Combined Alpha Beta Diveristy Based Score
#'
#' @name combinedShannonJSD
#'
#' @details Calculates a dysbiosis score that combined Shannon
#'          diversity and Jensen–Shannon divergence.
#'
#' The combined alpha-beta diversity approach was used by
#' \code{Santiago M E et al. 2019}. This approach uses Shannon diversity as the
#' alpha diversity measure and Jensen–Shannon divergence as the beta diversity
#' measure. The score is mean difference of Shannon diversity between test sample
#' and all references samples multiplied by the mean JSD of the test sample to all
#' reference samples. When calculating this score for reference samples, the
#' sample being used is excluded from calculating means for alpha and beta diversity.
#'
#' @param x A phyloseq object
#'
#' @param reference_samples Vector of samples to use as reference.
#'
#' @return A data frame with combinedShannonJSD score and sample information
#'
#' @examples
#' data("WirbelJ_2018")
#' library(phyloseq)
#' ps <- WirbelJ_2018
#' # Define controls as reference samples
#' ref.samples <- sample_names(subset_samples(WirbelJ_2018,
#'                                            disease == "healthy"))
#' alpha.beta.dysbiosis <- combinedShannonJSD(ps,
#'                                            reference_samples = ref.samples)
#' head(alpha.beta.dysbiosis)
#'
#' @references
#' \itemize{
#' \item{}{Santiago M et al. (2019).
#' Microbiome predictors of dysbiosis and VRE decolonization in patients with
#' recurrent C. difficile infections in a multi-center retrospective study.
#' \emph{AIMS Microbiol}, 5:1–18.}
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom phyloseq distance sample_names
#' @importFrom microbiome meta diversity
#' @importFrom utils packageVersion stack
#' @export
NULL
combinedShannonJSD <- function(x = NULL,
                               reference_samples = NULL){
  jsd.mat <- shannon.df <- all.res.score <- ids <- NULL

  if(is.null(reference_samples)){
    stop("Please provide reference_samples")
  }

  jsd.mat <- phyloseq::distance(x, method = "jsd")
  shannon.df <- microbiome::diversity(x, index = "shannon")

  # First for reference samples.
  ref.res.df <- NULL
  for(r in reference_samples){
    # exclude r from ref samples
    #ref.samples <- ref.samples[!r %in% ref.samples]
    ref.samples.s <- reference_samples[reference_samples != r]
    ref.shannon <- shannon.df[ref.samples.s,]
    test.sam <- shannon.df[r,]
    dif.shannon.r <- lapply(ref.shannon,
                            function(x){test.sam - x})
    shan.res.r <- unlist(dif.shannon.r) |> mean(na.rm = TRUE)
    jsd.mat.r <- jsd.mat |>
      as.matrix()
    jsd.test.r <- jsd.mat.r[ref.samples.s,r] |> mean(na.rm = TRUE)
    score.res.r <- shan.res.r * jsd.test.r
    ref.res.df[[r]] <- score.res.r
  }
  test.samples <- sample_names(x)[!sample_names(x) %in% reference_samples]
  test.res.df <- NULL
  for(t in test.samples){
    ref.shannon <- shannon.df[reference_samples,]
    test.sam.t <- shannon.df[t,]
    dif.shannon <- lapply(ref.shannon,
                          function(x){test.sam.t - x})
    shan.res <- unlist(dif.shannon) |> mean(na.rm = TRUE)
    jsd.mat <- jsd.mat |>
      as.matrix()
    jsd.test <- jsd.mat[reference_samples,t] |> mean(na.rm = TRUE)
    score.res.t <- shan.res * jsd.test
    test.res.df[[t]] <- score.res.t
  }
  all.res.score <- rbind(stack(test.res.df), stack(ref.res.df))
  colnames(all.res.score) <- c("ShannonJSDScore", "ids")
  rownames(all.res.score) <- all.res.score$ids
  return(cbind(dplyr::select(all.res.score, -ids), meta(x)))
}
