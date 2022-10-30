#'Overview of Dysbiosis Measures
#'
#' @name dysbiosisOverview
#'
#' @details Prints a list of dysbiosis score available in the
#'          \code{dysbiosisR} package.
#'
#'@author Sudarshan A. Shetty
#'
#' @export
NULL
dysbiosisOverview <- function(){

  cat("You are using dysbiosisR version: ")
  print(packageVersion('dysbiosisR'))

  cat("\nFollowing score (s) are currently supported\n")
  cat("1. dysbiosisMedianCLV: Lloyd-Price J, Arze C")
  cat("Ananthakrishnan AN et al. (2019)\n")
  cat("2. euclideanDistCentroids: AlShawaqfeh MK et al. (2017)\n")
  cat("3. cloudStatistic: Montassier E et al. (2018)\n")
  cat("4. combinedShannonJSD: Santiago M et al. (2019)\n")
  cat("5. dysbiosisOBB: Saffouri GB, Shields-Cutler R et al. (2019)\n")
  cat("6. distanceToReferencePlane: Halfvarson J, Brislawn CJ et al. (2017)\n")
}
