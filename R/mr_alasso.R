#' Adaptive MR-Lasso
#'
#' Fits the adaptive MR-Lasso estimator for two-sample Mendelian randomization.
#' The bootstrap option implements a bagged version of the estimator.
#'
#' @param object An object of class \code{MRInput}.
#' @param distribution Distribution used for confidence intervals in the
#'   non-bootstrap estimator. Either \code{"normal"} or \code{"t-dist"}.
#' @param alpha Significance level for confidence intervals.
#' @param lambda Optional tuning parameter for the adaptive Lasso.
#' @param bootstrap Logical. If \code{FALSE}, the ordinary MR-ALasso estimator
#'   is used. If \code{TRUE}, the bootstrap/bagged estimator is used.
#' @param B Number of bootstrap replications when \code{bootstrap = TRUE}.
#' @param agg_threshold Selection-frequency threshold used to report the
#'   aggregated valid set when \code{bootstrap = TRUE}.
#'
#' @return An object of class \code{MRAlasso} or \code{MRAlassoB}.
#' @export
mr_alasso = function(object,
                     distribution = c("normal", "t-dist"),
                     alpha = 0.05,
                     lambda = numeric(0),
                     bootstrap = FALSE,
                     B = 200,
                     agg_threshold = 0.5) {

  distribution = match.arg(distribution)
  if (!is.logical(bootstrap) || length(bootstrap) != 1) {
    stop("'bootstrap' must be TRUE or FALSE.")
  }

  if (bootstrap) {
    return(
      .mr_alasso_bootstrap(inp = object,alpha_ci = alpha,B = B,
        agg_threshold = agg_threshold)
    )
  }

  .mr_alasso_core(object = object,distribution = distribution,alpha = alpha,
    lambda = lambda)
}
