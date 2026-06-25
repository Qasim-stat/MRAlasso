#' @export
print.MRAlasso = function(x, ...) {
  cat("\nMR-Adaptive-Lasso method\n\n")
  cat("Number of instruments :", x$SNPs, "\n")
  cat("Number of valid instruments :", x$Valid, "\n")
  cat("Median ratio estimate :", x$PilotTheta, "\n")
  cat("-----------------------------------------------------\n")
  cat(sprintf("Estimate  Std Error   95%% CI (LL,UL)    p-value\n"))

  if (is.finite(x$Estimate)) {
    cat(sprintf(" %6.3f %9.3f %10.3f, %0.3f   %12.3g\n",
                x$Estimate, x$StdError, x$CILower, x$CIUpper, x$Pvalue))
  } else {
    cat(" NA\n")
  }

  cat("-----------------------------------------------------\n")
  invisible(x)
}

#' @export
print.MRAlassoB = function(x, ...) {
  cat("\nMR-Bagged-Adaptive-Lasso method\n\n")
  cat("Number of instruments :", x$SNPs, "\n")
  cat("Number of valid instruments :", x$Valid, "\n")
  cat("Number of invalid instruments :", x$Invalid, "\n")
  cat("-----------------------------------------------------\n")
  cat(sprintf("Estimate  Std Error   95%% CI (LL,UL)    p-value\n"))

  if (is.finite(x$est)) {
    cat(sprintf(
      " %6.3f %9.3f %10.3f, %0.3f   %12.3g\n",
      x$est, x$se, x$lwr, x$upr, x$Pvalue
    ))
  } else {
    cat(" NA\n")
  }

  cat("-----------------------------------------------------\n")
  invisible(x)
}
