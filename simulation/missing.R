# Misspecification study: PS model sensitivity (paper Table 3)
library(quantreg); library(ismev); library(actuar); library(evd); library(mev)
library(dplyr); library(tidyr); library(purrr); library(parallel); library(BMisc)
library(extRemes)

r_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
lapply(r_files, source)

# --- Configuration ---
h_fun_name <- "h3.fun"
h_fun      <- get(h_fun_name)
B          <- 1000
n_cores    <- 13
sig_idx    <- 1:2
shl_idx    <- 2
n          <- 1000
R          <- 1e6
interval   <- c(0.5, 100)
tau        <- seq(0.001, 0.9999, length.out = 1000)

pn_list  <- c(5 / n, 1 / n, 5 / (n * log(n)))
true_qte <- calculate_true_qtes(h_fun, pn_list)$qtes

# --- Main simulation loop ---
res_list <- parallel::mclapply(1:B, function(r) {

  tryCatch({
    simu <- h_fun(n)
    Y1 <- simu$Y1; Y0 <- simu$Y0
    D  <- simu$D;  X  <- simu$X
    Y  <- D * Y1 + (1 - D) * Y0

    mat <- vector("list", length(pn_list))
    for (j in seq_along(pn_list)) {
      pn    <- pn_list[j]
      k_seq <- n^(0.65)
      hn    <- floor(2 * n^(1 / 11))

      # Five misspecified PS models
      prop.fit <- glm(D ~ X, family = quasibinomial(link = "identity"))
      prop1 <- fitted(prop.fit)

      prop.fit <- glm(D ~ I(X^2), family = binomial(link = "logit"))
      prop2 <- fitted(prop.fit)

      Z <- rnorm(length(X))
      prop.fit <- glm(D ~ I(X^2) + Z + X:Z, family = binomial(link = "logit"))
      prop3 <- fitted(prop.fit)

      prop.fit <- glm(D ~ I(X^2), family = quasibinomial(link = "identity"))
      prop4 <- fitted(prop.fit)

      prop.fit <- glm(D ~ poly(X, hn), family = binomial)
      prop5 <- fitted(prop.fit)

      # TIEE with each misspecified PS
      tiee_prop1 <- tiee_cpp(Y, X, D, qn = 1 - pn, u = 1 - k_seq/n,
                              R = R, tau = tau, prop_scores = prop1,
                              sig_idx = sig_idx, shl_idx = shl_idx,
                              interval = interval, se = FALSE)
      tiee_prop2 <- tiee_cpp(Y, X, D, qn = 1 - pn, u = 1 - k_seq/n,
                              R = R, tau = tau, prop_scores = prop2,
                              sig_idx = sig_idx, shl_idx = shl_idx,
                              interval = interval, se = FALSE)
      tiee_prop3 <- tiee_cpp(Y, X, D, qn = 1 - pn, u = 1 - k_seq/n,
                              R = R, tau = tau, prop_scores = prop3,
                              sig_idx = sig_idx, shl_idx = shl_idx,
                              interval = interval, se = FALSE)
      tiee_prop4 <- tiee_cpp(Y, X, D, qn = 1 - pn, u = 1 - k_seq/n,
                              R = R, tau = tau, prop_scores = prop4,
                              sig_idx = sig_idx, shl_idx = shl_idx,
                              interval = interval, se = FALSE)
      tiee_prop5 <- tiee_cpp(Y, X, D, qn = 1 - pn, u = 1 - k_seq/n,
                              R = R, tau = tau, prop_scores = prop5,
                              sig_idx = sig_idx, shl_idx = shl_idx,
                              interval = interval, se = FALSE)

      mat[[j]] <- list(
        tiee_prop1 = tiee_prop1, tiee_prop2 = tiee_prop2,
        tiee_prop3 = tiee_prop3, tiee_prop4 = tiee_prop4,
        tiee_prop5 = tiee_prop5
      )
    }
    mat
  }, error = function(e) NA)
}, mc.cores = n_cores)

save(res_list, true_qte, pn_list, n, file = "scratch/missing.Rdata")

# --- Post-processing ---
extract_values <- function(model_obj, model_name) {
  if (all(is.na(model_obj))) return(NA)
  tryCatch(as.numeric(model_obj$eqte), error = function(e) NA)
}

cat("Processing results...\n")
df_raw <- map_dfr(seq_along(res_list), function(i) {
  iter_res <- res_list[[i]]
  if (is.null(iter_res) || all(is.na(iter_res))) return(NULL)
  map_dfr(seq_along(iter_res), function(j) {
    pn_res <- iter_res[[j]]
    map_dfr(names(pn_res), function(m_name) {
      data.frame(
        iter = i, pn_idx = j, pn_val = pn_list[j],
        model = m_name, true_qte = true_qte[j],
        est = extract_values(pn_res[[m_name]], m_name),
        row.names = NULL
      )
    })
  })
})

summary_stats <- df_raw %>%
  filter(!is.na(est)) %>%
  group_by(pn_idx, model) %>%
  summarise(
    pn_value = mean(pn_val),
    true_val = mean(true_qte),
    Bias    = mean(est) - mean(true_qte),
    AbsBias = abs(mean(est) - mean(true_qte)),
    SD      = sd(est),
    MSE     = mean((est - true_qte)^2),
    RMSE    = sqrt(mean((est - true_qte)^2)),
    n_valid = n(),
    .groups = "drop"
  ) %>%
  arrange(pn_idx, model)

print(as.data.frame(summary_stats))
