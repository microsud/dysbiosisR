#' Cloud-based LOcally linear Unbiased Dysbiosis (CLOUD) test
#'
#' @name cloudStatistic
#'
#' @details Calculates CLOUD score.
#'
#' Cloud-based LOcally linear Unbiased Dysbiosis (CLOUD) test is a non-parametric
#' test and returns a measure of dysbiosis. The function was adapted from
#' the original article by \code{Montassier E et al. 2018}. Here, a user defines
#' a set of reference samples from which distance of every other sample is
#' calculated. When calculating the \code{CLOUD} stats the \emph{k} is an
#' important parameter specified by argument \code{k_num}. By default we use a
#' conservative 80 percent of the samples in each group.
#'
#' @param x A phyloseq object
#'
#' @param reference_samples Vector of samples to use as reference.
#'
#' @param dist_mat A distance matrix. Can be output of \code{phyloseq::distance}
#'                 or \code{vegan::vegdist}.
#'
#' @param k_num Neighbors to use. Default is 5 percent. User can
#'              define percent of samples. See \code{Montassier E et al. 2018}
#'              for more details.
#'
#' @param ndim Dimension of the space in which the data are to be represented.
#'             Default is -1. See \code{Montassier E et al. 2018}
#'
#' @return A data frame with CLOUD stats and sample information
#'
#' @examples
#' data("WirbelJ_2018")
#' library(phyloseq)
#' ps <- WirbelJ_2018
#' # Define controls as reference samples
#' ref.samples <- sample_names(subset_samples(WirbelJ_2018,
#'                                            disease == "healthy"))
#' dist.data <- phyloseq::distance(ps, "bray")
#' cloud.results <- cloudStatistic(x= ps,
#'                                 dist_mat = dist.data,
#'                                 reference_samples = ref.samples,
#'                                 ndim=-1,
#'                                 k_num=5)
#' head(cloud.results)
#'
#' @references
#' \itemize{
#'
#' \item{}{Montassier E et al. (2018). CLOUD: a non-parametric detection test
#' for microbiome outliers.
#' \emph{Microbiome}, 6(1), pp.1-10.}
#'
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom stats cmdscale dist
#' @importFrom microbiome meta
#' @export
NULL
cloudStatistic <- function(x =NULL,
                           dist_mat = NULL,
                           reference_samples = NULL,
                           ndim=-1,
                           k_num=5){

  if(is.null(x) || is.null(dist_mat) || is.null(reference_samples)){
    stop("All arguments must be specified")
  }

  # reference scores
  dist_mat.ref = as.matrix(dist_mat)[reference_samples,reference_samples]

  # get number of neighnour from percent defined
  k.num.ref <- length(reference_samples) * (k_num/100)

  ref.stat.df <- .cloud_function(d=dist_mat.ref,
                                 test.ix=reference_samples,
                                 k=k.num.ref,
                                 ndim=-1)
  ref.stat.df$log2Stats <- log2(ref.stat.df$stats)
  rownames(ref.stat.df) <- ref.stat.df$ids
  #ref.stat.df$group <- "reference"

  # test scores
  test.samples <- as.matrix(dist_mat)
  test.samples <- test.samples[!colnames(test.samples) %in% reference_samples,]
  test.samples <- rownames(test.samples)

  # if(length(test.samples)==length(test.samples)){
  #   k.num <- length(reference_samples)-1
  # } else{
  #   k.num <- length(reference_samples)
  # }

  # get number of neighbor from percent defined
  k.num.test <- length(test.samples) * (k_num/100)

  test.stat.df <- .cloud_function(d = as.matrix(dist_mat),
                                  test.ix = test.samples,
                                  k = k.num.test, #k.num,
                                  ndim = -1)
  test.stat.df$log2Stats <- log2(test.stat.df$stats)
  rownames(test.stat.df) <- test.stat.df$ids
  all.res <- rbind(ref.stat.df, test.stat.df)

  main.sam.df <- meta(x)

  main.sam.df$ids <- rownames(main.sam.df)
  all.res.df <- merge(all.res, main.sam.df, by="ids")
  rownames(all.res.df) <- all.res.df$ids
  return(all.res.df[,-1])
}

# code from  Montassier E et al.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6080375/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6080375/bin/40168_2018_514_MOESM1_ESM.docx
.cloud_function <- function(d, test.ix=NULL, k=NULL, ndim=-1){
  res.df <- test.dist <- NULL
  # if (!inherits(d, "matrix")) d <- as.matrix(d)
  #if(class(d) != 'matrix')
  stats <- numeric(length(test.ix))
  pvals <- numeric(length(test.ix))

  for(i in 1:length(test.ix)){
    ref.ix <- test.ix[-i]
    keep.ix <- c(test.ix[i], ref.ix)
    if(ndim > -1){
      pc <- cmdscale(d[keep.ix,keep.ix,drop=F],k=ndim)
      d.i <- as.matrix(dist(pc))
    } else {
      d.i <- d[keep.ix,keep.ix,drop=F]
    }
    test.dist <- mean(sort(d.i[1,-1])[1:k])
    ref.dists <- numeric(length(ref.ix))
    for(j in 1:length(ref.ix)){
      ref.dists[j] <- mean(sort(d.i[-1,-1][j,-j]))
    }

    stats[i] <- test.dist / mean(ref.dists)
    pvals[i] <- mean(test.dist < ref.dists)
  }

  res.df<- data.frame(stats=stats,
                      pvals=pvals,
                      ids = test.ix)

  return(res.df)

}
