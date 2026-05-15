#!/usr/bin/env Rscript
library(quantreg); library(ismev); library(actuar); library(evd); library(mev)
library(dplyr); library(tidyr); library(parallel); library(BMisc)

# --- 1. Source helper functions ---
r_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_files, source))

# --- 2. Configuration ---
configs <- list(
  h_fun_name = "h1.fun",
  n          = 1000,
  R          = 1000,
  n_cores    = 14,
  tau_len    = 3000,
  sig_idx    = 1:2,
  shl_idx    = 1,
  interval   = c(0.5, 300),
  k_power    = 0.65
)

# Automatic path and file name generation
tag <- sprintf("%s_n%d_R%d", configs$h_fun_name, configs$n, configs$R)
output_dir <- file.path("scratch", configs$h_fun_name)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
log_file  <- file.path(output_dir, paste0(tag, ".txt"))
data_file <- file.path(output_dir, paste0(tag, ".Rdata"))

# --- 3. Initialise ---
h_fun     <- get(configs$h_fun_name)
pn_list   <- c(5 / configs$n, 1 / configs$n, 5 / (configs$n * log(configs$n)))
true_qte  <- calculate_true_qtes(h_fun, pn_list)$qtes

cat(sprintf("Simulation [%s] started at %s\n", tag, Sys.time()), file = log_file)

# --- 4. Parallel Monte Carlo loop ---
RNGkind("L'Ecuyer-CMRG")

res_list <- parallel::mclapply(1:configs$R, function(r) {
  if (r %% 50 == 0) {
    cat(sprintf("[Iter %d] PID %d at %s\n", r, Sys.getpid(), Sys.time()),
        file = log_file, append = TRUE)
  }

  tryCatch({
    simu <- h_fun(configs$n)
    Y1 <- simu$Y1; Y0 <- simu$Y0; D <- simu$D; X <- simu$X
    Y  <- D * Y1 + (1 - D) * Y0

    lapply(seq_along(pn_list), function(j) {
      pn    <- pn_list[j]
      k_seq <- configs$n^(configs$k_power)
      list(
        hill     = qte_extrapolation_hill(Y, X, D, pn, k_seq),
        pickland = qte_extrapolation_pickands(Y, X, D, pn, k_seq),
        zhang    = qte_firpo_zhang(Y, X, D, pn, N_bootstrap = 1000),
        tiee     = tiee_cpp(Y, matrix(X), D, qn = 1 - pn,
                            u = 1 - k_seq / configs$n,
                            R = 1e6, tau_length = configs$tau_len,
                            sig_idx = configs$sig_idx,
                            shl_idx = configs$shl_idx,
                            interval = configs$interval)
      )
    })
  }, error = function(e) {
    cat(sprintf("[ERROR Iter %d]: %s\n", r, e$message),
        file = log_file, append = TRUE)
    NA
  })
}, mc.cores = configs$n_cores, mc.set.seed = TRUE)

# --- 5. Save results ---
save(res_list, configs, pn_list, true_qte, file = data_file)
cat(sprintf("Simulation [%s] finished at %s\n", tag, Sys.time()),
    file = log_file, append = TRUE)

rm(res_list); gc()
