# Adaptive MR-Lasso (mr_aLasso)

.mr_alasso_core= function(object,
                      distribution = c("normal", "t-dist"),
                      alpha = 0.05, lambda = numeric(0)) {

  distribution = match.arg(distribution)

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required.")
  }

  #extract MRInput slots ----
  bx_raw = object@betaX
  by_raw = object@betaY
  bxse   = object@betaXse
  byse   = object@betaYse

  Bx = abs(bx_raw)
  By = by_raw * sign(bx_raw)

  nsnps = length(Bx)

  if (length(By) != nsnps || length(bxse) != nsnps || length(byse) != nsnps) {
    stop("Input lengths do not match.")
  }
  if (any(!is.finite(Bx)) || any(!is.finite(By)) ||
      any(!is.finite(bxse)) || any(!is.finite(byse))) {
    stop("Non-finite inputs detected.")
  }
  if (any(byse <= 0) || any(bxse <= 0)) {
    stop("Standard errors must be positive.")
  }

  if (!(distribution %in% c("normal", "t-dist"))) {
    stop("Distribution must be one of: 'normal', 't-dist'.")
  }

  #projection ----
  S = diag(byse^-2)
  S_half = diag(byse^-1)

  b = S_half %*% matrix(Bx, ncol = 1)
  Pb = b %*% solve(t(b) %*% b, t(b))

  xlas = (diag(nsnps) - Pb) %*% S_half
  ylas = (diag(nsnps) - Pb) %*% S_half %*% matrix(By, ncol = 1)


  #Adaptive weights ----
  ok = (Bx != 0) & is.finite(Bx) & is.finite(By)
  if (sum(ok) < 3) {
    stop("Too few SNPs with nonzero betaX for median-ratio pilot.")
  }

  pi_hat = By[ok] / Bx[ok]
  theta_pilot = stats::median(pi_hat, na.rm = TRUE)
  nu = 1; eps = 1e-6;
  alpha_tilde = By - theta_pilot * Bx
  omega = 1 / (abs(alpha_tilde) + eps)^nu

  finite_omega = omega[is.finite(omega)]
  if (length(finite_omega) == 0) {
    stop("All adaptive weights are non-finite.")
  }
  omega[!is.finite(omega)]  =max(finite_omega)
  # Fit adaptive lasso path ----
  if (length(lambda) != 0) {
    las_fit = glmnet::glmnet(x = xlas,y = ylas,intercept = FALSE,lambda = lambda,
                             penalty.factor = omega,standardize = FALSE)

    las_mod = list(
      fit = as.numeric(las_fit$beta[, 1]),lambda = lambda)
  } else {
    las_fit = glmnet::glmnet(x = xlas,y = ylas,intercept = FALSE,penalty.factor = omega,
                             standardize = FALSE)

    lamseq = sort(las_fit$lambda)
    lamlen = length(lamseq)
    rse = sapply(seq_len(lamlen), function(j) {
      col_idx = lamlen - j + 1
      av = which(las_fit$beta[, col_idx] == 0)
      #Need at least 2 "valid" SNPs to fit post model----
      if (length(av) <= 1) {
        return(c(NA_real_, length(av)))
      }

      mod = stats::lm(S_half[av, av] %*% By[av] ~ S_half[av, av] %*% Bx[av] - 1)
      c(sqrt(sum(mod$residuals^2) / mod$df.residual),length(av))
    })

    # rse[1, ] = residual SE, rse[2, ] = number selected valid
    rse_inc = rse[1, 2:lamlen] - rse[1, 1:(lamlen - 1)]

    het = which(is.finite(rse[1, 2:lamlen]) &
                  is.finite(rse_inc) &
                  is.finite(rse[2, 2:lamlen]) &
                  rse[1, 2:lamlen] > 1 &
                  rse_inc > (stats::qchisq(0.95, 1) / rse[2, 2:lamlen]))

    if (length(het) == 0) {
      lam_pos = lamlen
    } else {
      lam_pos =min(het)
    }

    num_valid = rse[2, ]
    feasible = which(num_valid > 1)

    if (length(feasible) == 0) {
      stop("No lambda in the path yields more than one selected valid instrument.")
    }

    min_lam_pos = min(feasible)
    if (lam_pos < min_lam_pos) {
      lam_pos = min_lam_pos}

    las_mod = list(fit = as.numeric(las_fit$beta[, (lamlen - lam_pos + 1)]),
                   lambda = lamseq[lam_pos])
  }

  a = las_mod$fit; e = By - a

  reg_est = as.numeric(solve(t(Bx) %*% S %*% Bx, t(Bx) %*% S %*% matrix(e, ncol = 1)))

  v = which(a == 0)

  if (length(v) > 1) {
    post_mod = summary(stats::lm(By[v] ~ Bx[v] - 1, weights = byse[v]^-2))
    post_est = as.numeric(post_mod$coef[, 1])
    post_se  = as.numeric(post_mod$coef[, 2]) / min(post_mod$sigma, 1)
  } else {
    post_est =NA_real_
    post_se  =NA_real_
    warning("Chosen lambda yields fewer than two valid instruments. Post-adaptive lasso method cannot be performed.")
  }

  #post-ALasso ----
  if (is.finite(post_est) && is.finite(post_se)) {
    if (distribution == "normal") {
      zcrit = stats::qnorm(1 - alpha / 2)
      ciLower = post_est - zcrit * post_se
      ciUpper = post_est + zcrit * post_se
      pvalue  = 2 * stats::pnorm(-abs(post_est / post_se))
    } else {
      df = length(v) - 1
      tcrit = stats::qt(1 - alpha / 2, df = df)
      ciLower = post_est - tcrit * post_se
      ciUpper = post_est + tcrit * post_se
      pvalue  = 2 * stats::pt(-abs(post_est / post_se), df = df)
    }
  } else {
    ciLower= NA_real_
    ciUpper= NA_real_
    pvalue= NA_real_
  }

  valid_snps = if (!is.null(object@snps)) as.character(object@snps[v]) else NULL

  out = list(
    Estimate = post_est,
    StdError = post_se,
    CILower = ciLower,
    CIUpper = ciUpper,
    Alpha = alpha,
    Pvalue = pvalue,
    SNPs = nsnps,
    RegEstimate = reg_est,
    RegIntercept = a,
    Valid = length(v),
    ValidIdx = v,
    ValidSNPs = valid_snps,
    Lambda = las_mod$lambda,
    PilotTheta = theta_pilot,
    AdaptiveWeights = omega
  )

  class(out) = "MRAlasso"
  out
}


