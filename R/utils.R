


#' @keywords Internal
.check_input_ref_samples <- function(dist.mat, ref.samples){

  in.sams <- intersect(colnames(as.matrix(dist.mat)), ref.samples)

  if(length(in.sams)!=length(ref.samples)){
    # length(ref.samples) - length(in.sams)
    stop("Distance matrix does not have same number of reference samples",
         "as input reference_samples",
         call. = FALSE)
  }

}
