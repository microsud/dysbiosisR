#' Score based on the distance to a healthy/reference plane
#'
#' @name distanceToReferencePlane
#'
#' @details Calculates a 'healthy' or 'reference plane'. The plane is calculated in a space
#'          derived from PCoA and based on user-defined distances between samples from healthy
#'          subjects. Briefly, a model is constructed using the samples from healthy subjects,
#'          which were fitted to a two-dimensional plane embedded in a three-dimensional space using
#'          the least-squares method. The plane is then restricted to only span the three-dimensional
#'          ranges of the healthy control samples. The plane is considered a proxy for the normal
#'          microbial variation in healthy subjects. \code{distanceToReferencePlane} calculates the
#'          Euclidean distance (within PCoA space) from each sample towards the plane, which can
#'          be considered a measure of (ab)normality. The user can provide a custom distance matrix.
#'          If the user provides a UniFrac distance matrix, then the resulting score is comparable to
#'          the \code{dysbiosis score} reported in \code{Halfvarson J, Brislawn CJ, Lamendella R et al.
#'          (2017)}.
#'
#' @param x A phyloseq object
#'
#' @param dist_mat A distance matrix. Can be output of \code{phyloseq::distance}
#'                 or \code{vegan::vegdist}.
#'
#' @param reference_samples Vector of samples to use as reference.
#'
#' @return A data frame with distance to reference plane ('dtrpScore') values and sample information.
#'
#' @examples
#' library(dysbiosisR)
#' # We use WirbelJ_2018 as test data
#' dist.mat <- phyloseq::distance(WirbelJ_2018, "bray")
#' ref.samples <- sample_names(subset_samples(WirbelJ_2018,
#'                                            disease == "healthy"))
#' dtrp.results <- distanceToReferencePlane(WirbelJ_2018,
#'                                          dist_mat = dist.mat,
#'                                          reference_samples = ref.samples)
#'
#' head(dtrp.results)
#'
#' @references
#' \itemize{
#' \item{}{Halfvarson J, Brislawn CJ, Lamendella R et al. (2017)
#' Dynamics of the human gut microbiome in inflammatory bowel disease.
#' \emph{Nature Microbiology}, 2, article number: 17004.}
#' \item{}{Vázquez-Baeza Y. (2017) Reference Plane, GitHub repository, https://github.com/ElDeveloper/reference-plane}
#' }
#'
#' @author Wouter A.A. de Steenhuijsen Piters
#'
#' @importFrom stats coef lm
#' @importFrom phyloseq sample_names ordinate
#' @importFrom microbiome meta diversity
#' @importFrom utils packageVersion stack
#' @export
NULL
distanceToReferencePlane <- function(x = NULL,
                                     dist_mat = NULL,
                                     reference_samples = NULL) {

  if(is.null(x) || is.null(dist_mat) || is.null(reference_samples)){
    stop("All arguments must be specified")
  }

  if(!all(sample_names(x) == rownames(as(dist_mat, "matrix"))) ||
     !all(sample_names(x) == colnames(as(dist_mat, "matrix")))) {
    stop("Phyloseq `sample_names()` and row- and/or colnames of the distance matrix are not in the same order")
  }

  ord <- ordinate(x, method = "PCoA", distance = dist_mat)
  reference <- ord$vectors[reference_samples, 1:3]
  abcd = .compute_coefficients(reference)

  dtrp.res <- apply(ord$vectors[, 1:3], 1, function(x) .point_to_segment_distance(abcd, x, reference)) |>
    stack()
  colnames(dtrp.res) <- c("dtrpScore", "ids")

  main.sam.df <- meta(x)
  main.sam.df$ids <- rownames(main.sam.df)

  all.res.df <- merge(dtrp.res, main.sam.df, by = "ids", sort = F)
  rownames(all.res.df) <- all.res.df$ids

  return(all.res.df[,-1])
}

.compute_coefficients <- function(xyz) {

  # Original Python-implementation: https://github.com/ElDeveloper/reference-plane
  # Fit a plane to the first three dimensions of a matrix
  #
  # xyz: The matrix of data to fit the plane to.
  x <- xyz[, 1]
  y <- xyz[, 2]
  z <- xyz[, 3]

  A <- cbind(x, y)

  mod <- lm(z ~. , data = data.frame(A))
  abcd <- c(coef(mod)[[2]], coef(mod)[[3]], -1, coef(mod)[[1]])
  # Returns a vector with coefficients `a`, `b`, `c` and `d` in the equation:
  #  a*x + b*y + c*z + d = 0.
  return(abcd)
}

.point_to_plane_distance <- function(abcd, point) {
  # Original Python-implementation: https://github.com/ElDeveloper/reference-plane
  # Calculates the euclidean distance from a point to a plane
  #
  # abcd: The four coefficients of an equation that defines a
  # plane of the form a*x + b*y + c*z + d = 0.
  #
  # point: The values for x, y and z for the point that you want
  # to calculate the distance to.

  abc = abcd[1:3]
  d = abcd[4]
  dist = abs((abc %*% point) + d)
  dist/norm(matrix(abc), "F")
  # Returs the distance from the point to the plane.
}

.point_to_segment_distance <- function(abcd, point, xyz) {
  # Original Python-implementation: https://github.com/ElDeveloper/reference-plane
  # Compute the distance from a point to a segment of a plane
  #
  # point: a point in PCoA-space.
  #
  # xyz: 3-dimensional matrix of sample coordinates.
  #
  # abcd: The four coefficients of an equation that defines a
  # plane of the form a*x + b*y + c*z + d = 0.
  plane <- function(abcd, xy) {
    a = abcd[1];b = abcd[2];c = abcd[3];d = abcd[4]
    x = xy[1];y = xy[2]
    return (a*x + b*y + d)/(-1*c)
  }

  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2)^2))

  a = abcd[1];b = abcd[2];c = abcd[3];d = abcd[4]
  p = point[1];q = point[2];r = point[3]
  l = ((d*-1) - p*a - b*q - c*r) / (a^2 + b^2 + c^2)
  extreme = c(p + l*a, q + l*b, r + l*c)

  for (i in 1:ncol(xyz)) {
    vector = xyz[,i]
    ranges = range(vector)
    if(extreme[i] < ranges[1]) { extreme[i] <- ranges[1]
    } else if(extreme[i] > ranges[2]) { extreme[i] <- ranges[2] } }
  extreme[length(extreme)] = plane(abcd, extreme[-length(extreme)])
  return(euc.dist(point, extreme))
  # returns the distance from the point to the segment of the plane that
  # spans the `xyz`-space.
}
