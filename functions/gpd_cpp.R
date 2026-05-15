#' Variance contribution from propensity score estimation (IPW version)
v_ps <- function(X_ps, V_ps, phi, D, ps, eff = eff, n = n) {
  w1 <- -(1 - ps) * (D / ps) * phi / sum(D / ps)
  A_gamma <- drop(crossprod(X_ps, w1) / n)
  M_ps <- crossprod(X_ps, X_ps * (ps * (1 - ps)))
  w2 <- (D - ps) * (D / ps / sum(D / ps) * phi - eff)
  Sigma <- drop(crossprod(X_ps, w2) / n)
  as.numeric(t(A_gamma) %*% V_ps %*% A_gamma - 2 * t(A_gamma) %*% solve(M_ps, Sigma))
}

#' Variance contribution from propensity score estimation (DR version)
v_ps_dr <- function(X_ps, V_ps, phi, D, ps, eff = eff, n = n) {
  w1 <- -(1 - ps) * (D / ps) * phi
  A_gamma <- drop(crossprod(X_ps, w1) / n)
  M_ps <- crossprod(X_ps, X_ps * (ps * (1 - ps)))
  w2 <- (D - ps) * (D / ps * phi - eff)
  Sigma <- drop(crossprod(X_ps, w2) / n)
  as.numeric(t(A_gamma) %*% V_ps %*% A_gamma - 2 * t(A_gamma) %*% solve(M_ps, Sigma))
}

#' Variance contribution from quantile regression threshold estimation
v_beta <- function(fit_rq, group = 1, Y, w, gpd_hat, thres, sigma, xi, phi, eff, u) {
  n <- length(Y)
  X_beta <- cbind(rep(1, n), group)
  sum_rq <- summary(fit_rq, covariance = TRUE, se = "nid")
  A_beta <- mean(-w * sapply(1:n, function(x)
    dgpd(gpd_hat, loc = 0, scale = sigma[x], shape = xi[x]))) * c(1, group)
  V_beta <- sum_rq$cov
  M_inv  <- sum_rq$Hinv
  phi_beta <- sweep(X_beta, 1, u - (Y <= X_beta %*% coef(fit_rq)), "*")
  center <- function(x) x - mean(x)
  Sigma_beta <- colMeans(center(phi_beta) * center(w * phi - eff))
  A_beta %*% V_beta %*% t(t(A_beta)) - 2 * A_beta %*% M_inv %*% t(t(Sigma_beta))
}

#' Variance contribution from GPD parameter estimation
v_gpd <- function(wf = w1 * f1,
                  X_sig = X1[, c(1, sig_idx + 1)],
                  X_xi  = X1[, c(1, shl_idx + 1)],
                  V_gpd, score, sigma, xi, tau, diffs, u,
                  g = w1 * phi1 - eff) {
  G <- G_gpd_cpp(tau, diffs, xi, sigma, u)
  dQ <- wf * cbind(G[, 2] * sigma * X_sig, G[, 1] * X_xi)
  A_gamma_gpd <- as.numeric(colMeans(dQ))
  M_inv <- V_gpd
  Sigma_gpd <- crossprod(score, g) / n
  as.numeric(A_gamma_gpd %*% V_gpd %*% t(t(A_gamma_gpd)) -
             2 * A_gamma_gpd %*% M_inv %*% Sigma_gpd)
}


#' TIEE: Tail-Informed Extreme Estimation for quantile treatment effects
#'
#' @param Y        Observed outcome vector
#' @param X        Covariate matrix
#' @param D        Binary treatment indicator
#' @param qn       Target quantile level (1 - p_n)
#' @param u        Threshold quantile for GPD tail model
#' @param R        Number of quasi-samples for KDE integration
#' @param tau      Grid for GPD quantile approximation
#' @param prop_scores  Optional pre-computed propensity scores
#' @param sig_idx  Column indices in X for GPD scale parameter
#' @param shl_idx  Column indices in X for GPD shape parameter
#' @param interval Search interval for theta optimization
#' @param heavy    Logical; heavy-tailed DGP
#' @param se       Logical; compute standard errors
tiee_cpp <- function(Y, X, D,
                     qn = 1 - 5 / n, u = 0.9, R = 1e6,
                     tau = tau, prop_scores = NULL,
                     sig_idx = 1:2, shl_idx = 1,
                     interval = c(0.5, 50), heavy = TRUE, se = TRUE) {

  n <- length(Y)

  # --- Quantile regression thresholds ---
  fit_rq <- quantreg::rq(Y ~ D, tau = u)
  thres  <- predict(fit_rq)
  thres1 <- predict(fit_rq, newdata = data.frame(X = X, D = 1))
  thres0 <- predict(fit_rq, newdata = data.frame(X = X, D = 0))

  # --- GPD fit on exceedances ---
  resid <- Y - thres
  idx   <- resid > 0
  fit_gpd <- ismev::gpd.fit(
    resid[idx], threshold = 0,
    ydat = cbind(D, X)[idx, ],
    sigl = sig_idx, shl = shl_idx,
    show = FALSE, siglink = exp  # log-sigma to ensure positivity
  )
  param <- fit_gpd$mle

  # --- Predicted scale and shape for D=1 and D=0 ---
  X1 <- as.matrix(cbind(1, D = 1, X))
  X0 <- as.matrix(cbind(1, D = 0, X))

  sigma1 <- exp(as.vector(X1[, c(1, sig_idx + 1)] %*% param[c(1, sig_idx + 1)]))
  sigma0 <- exp(as.vector(X0[, c(1, sig_idx + 1)] %*% param[c(1, sig_idx + 1)]))
  xi1    <- as.vector(X1[, c(1, shl_idx + 1)] %*% param[-c(1, sig_idx + 1)])
  xi0    <- as.vector(X0[, c(1, shl_idx + 1)] %*% param[-c(1, sig_idx + 1)])

  # --- GPD quantile estimates ---
  diffs <- diff(tau)
  tau   <- tau[-length(tau)]

  gpd_hat1 <- qgpd_matrix(tau, sigma1, xi1)
  gpd_hat0 <- qgpd_matrix(tau, sigma0, xi0)
  Yd1 <- gpd_hat1 + thres1
  Yd0 <- gpd_hat0 + thres0

  # --- Propensity scores ---
  if (is.null(prop_scores)) {
    hn <- floor(2 * n^(1 / 11))
    prop.fit <- glm(D ~ poly(X, hn), family = binomial)
    prop_scores <- fitted(prop.fit)
  }

  eff <- (qn - u) / (1 - u)
  w1  <- D / prop_scores
  w0  <- (1 - D) / (1 - prop_scores)

  # --- Point estimates ---
  result1 <- optimize(obj_cpp, interval, Yd = Yd1, w = w1,
                       diffs = diffs, eff = eff, R = R)
  result0 <- optimize(obj_cpp, interval, Yd = Yd0, w = w0,
                       diffs = diffs, eff = eff, R = R)
  theta1 <- result1$minimum
  theta0 <- result0$minimum
  eqte   <- theta1 - theta0

  if (!isTRUE(se)) {
    return(list(eqte = eqte, theta1 = theta1, theta0 = theta0))
  }

  # --- Standard errors via influence-function decomposition ---
  f1 <- dgpd_cpp(theta1 - thres1, scale = sigma1, shape = xi1) * (1 - u)
  f0 <- dgpd_cpp(theta0 - thres0, scale = sigma0, shape = xi0) * (1 - u)

  phi1 <- (Yd1 <= theta1) %*% diffs
  phi0 <- (Yd0 <= theta0) %*% diffs

  # Variance of weighted empirical process
  v_theta_1 <- var(w1 * phi1)
  v_theta_0 <- var(w0 * phi0)

  # PS variance contribution
  X_ps <- model.matrix(prop.fit)
  V_ps <- vcov(prop.fit)
  v_ps1 <- v_ps(X_ps, V_ps, phi1, D, prop_scores, eff, n)
  v_ps0 <- v_ps(X_ps, V_ps, phi0, 1 - D, 1 - prop_scores, eff, n)

  # GPD variance contribution
  score <- gpd_cpp(param, resid, cbind(1, D, X),
                   scale_cols = c(1, sig_idx + 1) - 1,  # 0-based indexing
                   shape_cols = c(1, shl_idx + 1) - 1)
  V_gpd <- fit_gpd$cov

  v_gpd1 <- v_gpd(wf = w1 * f1, X_sig = X1[, c(1, sig_idx + 1)],
                   X_xi = X1[, c(1, shl_idx + 1)],
                   V_gpd, score, sigma1, xi1, tau, diffs, u,
                   g = w1 * phi1 - eff)
  v_gpd0 <- v_gpd(wf = w0 * f0, X_sig = X0[, c(1, sig_idx + 1)],
                   X_xi = X0[, c(1, shl_idx + 1)],
                   V_gpd, score, sigma0, xi0, tau, diffs, u,
                   g = w0 * phi0 - eff)

  v1 <- v_theta_1 + v_ps1 + v_gpd1
  v0 <- v_theta_0 + v_ps0 + v_gpd0
  cov2 <- cov(w1 * phi1, w0 * phi0)

  # Density estimation via KDE on a grid of theta values
  tau2 <- seq(tau[1], tail(tau, 1), length.out = 100)
  sum_diffs <- sum(diffs)
  d1 <- diffs[1]
  delta <- 0.02 * diff(interval)
  tol <- 1e-4
  max_iter <- 200

  kde1 <- kde_path_cpp(tau2, interval, Yd1, w1, d1, sum(w1),
                        sum_diffs, R, delta, tol, max_iter)
  kde0 <- kde_path_cpp(tau2, interval, Yd0, w0, d1, sum(w0),
                        sum_diffs, R, delta, tol, max_iter)

  d_SJ1 <- density(kde1, bw = "SJ")
  d_SJ0 <- density(kde0, bw = "SJ")

  f1 <- approx(x = d_SJ1$x, y = d_SJ1$y, xout = theta1, rule = 2)$y
  f0 <- approx(x = d_SJ0$x, y = d_SJ0$y, xout = theta0, rule = 2)$y

  # Asymptotic standard error
  se <- sqrt((v1 / f1^2 + v0 / f0^2 - 2 * cov2 / f1 / f0) / n)
  lower <- eqte - se * qnorm(0.95)
  upper <- eqte + se * qnorm(0.95)

  list(
    eqte = eqte, theta1 = theta1, theta0 = theta0,
    f1 = f1, f0 = f0, se = se, lower = lower, upper = upper
  )
}
