# mr_aLasso-b

bootstrap_snp_counts = function(J) {
  as.integer(stats::rmultinom(1, size = J, prob = rep(1 / J, J)))
}

make_bootstrap_input_from_counts = function(inp, counts) {
  idx_rep = rep(seq_along(counts), counts)
  if (length(idx_rep) < 2) return(NULL)

  inp_b = MendelianRandomization::mr_input(
    bx = inp@betaX[idx_rep],
    bxse = inp@betaXse[idx_rep],
    by = inp@betaY[idx_rep],
    byse = inp@betaYse[idx_rep],
    snps = inp@snps[idx_rep])

  list(input = inp_b,idx_rep = idx_rep)
}

extract_valid_alasso = function(fit, inp, idx_rep = NULL) {
  if (is.null(fit)) return(integer(0))

  #bootstrap samples----
  if (!is.null(idx_rep) && !is.null(fit$ValidIdx)) {
    idx_b = as.integer(fit$ValidIdx)
    idx_b = idx_b[is.finite(idx_b)]
    idx_b = idx_b[idx_b >= 1 & idx_b <= length(idx_rep)]
    return(sort(unique(idx_rep[idx_b])))
  }

  #Fallback: use SNP names if available----
  if (!is.null(fit$ValidSNPs) && !is.null(inp@snps)) {
    idx = match(unique(fit$ValidSNPs), inp@snps)
    idx = idx[is.finite(idx)]
    return(sort(unique(idx)))
  }

  #only for non-bootstrap fits----
  idx = as.integer(fit$ValidIdx)
  idx = idx[is.finite(idx)]
  sort(unique(idx))
}

post_ivw_from_valid_weighted = function(inp, valid_idx, weight_counts, alpha_ci = 0.05) {
  valid_idx = sort(unique(as.integer(valid_idx)))
  valid_idx = valid_idx[is.finite(valid_idx)]
  if (length(valid_idx) < 2) return(NULL)

  wc = weight_counts[valid_idx]
  keep = is.finite(wc) & wc > 0
  valid_idx = valid_idx[keep]
  wc = wc[keep]
  if (length(valid_idx) < 2) return(NULL)

  bx = inp@betaX[valid_idx]
  by = inp@betaY[valid_idx]
  byse = inp@betaYse[valid_idx]

  w = wc * (byse^-2)
  denom = sum(w * bx^2)
  if (!is.finite(denom) || denom <= 0) return(NULL)

  est = sum(w * bx * by) / denom
  resid = by - est * bx
  df = length(valid_idx) - 1
  if (df < 1) return(NULL)
  sigma2 = sum(w * resid^2) / df
  se_raw = sqrt(sigma2 / denom)
  se = se_raw / min(sqrt(sigma2), 1)
  z = stats::qnorm(1 - alpha_ci / 2)

  list(est = est,se = se,lwr = est - z * se,upr = est + z * se)
}

.mr_alasso_bootstrap = function(inp, alpha_ci = 0.05, B = 200, agg_threshold = 0.5) {
  J = length(inp@betaX)
  theta_b = rep(NA_real_, B)
  selected_list = vector("list", B)
  w_mat = matrix(0, nrow = J, ncol = B)

  for (b in seq_len(B)) {
    counts_b = bootstrap_snp_counts(J)
    w_mat[, b] = counts_b

    boot_obj = make_bootstrap_input_from_counts(inp, counts_b)
    if (is.null(boot_obj)) next
    inp_b = boot_obj$input
    idx_rep = boot_obj$idx_rep

    fit_b = tryCatch(
      .mr_alasso_core(inp_b, distribution = "normal", alpha = alpha_ci),
      error = function(e) NULL)
    if (is.null(fit_b)) next

    valid_idx_b = extract_valid_alasso(fit_b, inp, idx_rep = idx_rep)
    selected_list[[b]] = valid_idx_b

    if (length(valid_idx_b) < 2) next

    fit_theta_b = post_ivw_from_valid_weighted(inp,valid_idx_b,counts_b,
                                               alpha_ci = alpha_ci)
    if (is.null(fit_theta_b)) next

    theta_b[b] = fit_theta_b$est
  }

  ok = is.finite(theta_b)
  if (sum(ok) < 2) return(NULL)
  theta_tilde = mean(theta_b[ok])

  #variance estimate----
  w_centered = sweep(w_mat[, ok, drop = FALSE],1,
                     rowMeans(w_mat[, ok, drop = FALSE]),FUN = "-")
  theta_centered = theta_b[ok] - theta_tilde
  S_hat = as.vector(w_centered %*% theta_centered) / sum(ok)
  var_hat = sum(S_hat^2)
  se_hat = sqrt(max(var_hat, 0))
  z = stats::qnorm(1 - alpha_ci / 2)
  pvalue = if (is.finite(se_hat) && se_hat > 0) {
    2 * stats::pnorm(-abs(theta_tilde / se_hat))
  } else {
    NA_real_
  }

  #Aggregated selected set for reporting only----
  sel_freq = rep(0, J)
  for (b in which(ok)) {
    if (!is.null(selected_list[[b]]) && length(selected_list[[b]]) > 0) {
      sel_freq[selected_list[[b]]] = sel_freq[selected_list[[b]]] + 1
    }
  }

  sel_freq = sel_freq / sum(ok)
  agg_valid_idx = which(sel_freq >= agg_threshold)
  invalid_idx = setdiff(seq_len(J), agg_valid_idx)

  out = list(
    est = theta_tilde,
    se = se_hat,
    lwr = theta_tilde - z * se_hat,
    upr = theta_tilde + z * se_hat,
    Pvalue = pvalue,
    SNPs = J,
    Valid = length(agg_valid_idx),
    valid_idx = agg_valid_idx,
    Invalid = length(invalid_idx),
    ValidSNPs = inp@snps[agg_valid_idx],
    InvalidSNPs = inp@snps[invalid_idx],
    invalid_idx = invalid_idx,
    theta_boot = theta_b,
    sel_freq = sel_freq,
    n_boot_success = sum(ok)
  )

  class(out) = "MRAlassoB"
  out
}


