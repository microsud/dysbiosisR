#' Random Forest Prediction Based Score
#'
#' @name dysbiosisOBB
#'
#' @details An out-of-bag probability of Random Forest of being classified as
#'          case or diseased group.
#'
#' The original article \code{Saffouri GB, Shields-Cutler R et al. 2019}
#' reported a Symptom Index abbreviated as SI. In this approach the feature
#' abundances are used for \code{Random Forest} classification and the resulting
#' out of bag (OOB) predicted probability of being classified in disease group
#' is considered as an SI or also dysbiosis index. The \code{dysbiosisOBB}
#' function in this package allows for calculating this measure with some level
#' of freedom of 'tuneRF' and 'randomForest' parameters via the \code{randomForest}
#' R package.
#'
#' @param x A phyloseq object
#'
#' @param group_col A column in \code{phyloseq::sample_data} with all control and
#'                  case labels.
#'
#' @param case_label A character string specifying case/disease labels in
#'                   group_col.
#'
#' @param seed_value The random seed of your session for reproducibility.
#'
#' @param add_tuneRF_params A list of arguments to pass for
#'                          \code{randomForest::tuneRF}
#'
#' @param ntree Number of trees. See \code{randomForest::randomForest}
#'
#' @param plot_roc Logical TRUE or FALSE to plot ROC curve.
#'
#' @param ... Additional arguments if necessary for
#'            \code{randomForest::randomForest}.
#'
#' @return A data frame with oob.score score and sample information
#'
#' @examples
#' data("WirbelJ_2018")
#' # data are relative abundances summed to 100 or
#' # near 100 val
#' dysbiosis.oob <- dysbiosisOBB(WirbelJ_2018,
#'                               group_col = "disease",
#'                               case_label = "CRC",
#'                               seed_value = 1235,
#'                               add_tuneRF_params = list(ntreeTry=100,
#'                                                        stepFactor=1.5,
#'                                                        improve=0.01,
#'                                                        trace=TRUE,
#'                                                        dobest=FALSE),
#'                               ntree = 100, # increase for real data
#'                               plot_roc = TRUE)
#' head(dysbiosis.oob)
#'
#' @references
#' \itemize{
#' \item{}{Saffouri GB, Shields-Cutler R et al. (2019).
#' Small intestinal microbial dysbiosis underlies symptoms associated with
#' functional gastrointestinal disorders.
#' \emph{Nature communications}, 10(1), pp.1-11.}
#' }
#'
#' @author Sudarshan A. Shetty
#'
#' @importFrom randomForest randomForest tuneRF
#' @importFrom microbiome meta abundances
#' @importFrom pROC roc ci
#' @export
NULL
dysbiosisOBB <- function(x= NULL,
                         group_col = NULL,
                         case_label = NULL,
                         seed_value = 1235,
                         add_tuneRF_params = list(ntreeTry=1000,
                                                  stepFactor=1.5,
                                                  improve=0.01,
                                                  trace=TRUE,
                                                  dobest=FALSE),
                         ntree = 1000,
                         plot_roc = TRUE,
                         ...){

  prop <- group.var <- arg.list.tuneRF.all <- bestmtry <- min.oob <- NULL
  OOBError <- forest.res <- case.prob <- roc2 <- mtry <- NULL
  if(is.null(x) ||
     is.null(group_col) ||
     is.null(case_label)){
    stop("All arguments must be specified")
  }
  message("The random seed of your session for reproducibility is: ", seed_value)

  # Get abundance
  prop <- abundances(x)

  if(is.factor(meta(x)[,group_col])){
    group.var <- meta(x)[,group_col]
  } else {
    group.var <- as.factor(meta(x)[,group_col])
  }

  arg.list.tuneRF.all <- c(list(x=t(prop),
         y=group.var,
         sampsize=table(group.var)), add_tuneRF_params)
  set.seed(seed_value)
  bestmtry <- do.call('tuneRF',arg.list.tuneRF.all)|>
    as.data.frame()

  min.oob <- bestmtry |>
    dplyr::filter(OOBError ==  min(bestmtry$OOBError)) |>
    dplyr::pull(mtry)

  message("The mtry used is : ", min.oob)

  forest.res <- randomForest(x=t(prop),
                             y=group.var,
                             sampsize=table(group.var),
                             importance = TRUE,
                             ntree = ntree,
                             mtry = min.oob,
                             ...)
  #
  case.prob <- stats::predict(forest.res,  type='prob')[, case_label]
  if(plot_roc){
    roc2 <- pROC::roc(group.var,
                      case.prob ,
                      plot=TRUE,
                      ci = TRUE,
                      auc.polygon=TRUE,
                      max.auc.polygon=TRUE,
                      grid=TRUE,
                      print.auc=TRUE)
    #message(ci(roc2))
  }

  res.df <- data.frame(oob.score = case.prob)
  return(cbind(res.df, meta(x)))

}








