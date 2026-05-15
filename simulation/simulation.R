library(quantreg); library(ismev); library(actuar); library(evd); library(mev)
library(dplyr); library(tidyr); library(parallel); library(BMisc)

# Source all helper functions
r_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
lapply(r_files, source)

# ==============================================================================
# Configuration — edit these for each run
# ==============================================================================
h_fun_name <- "h3.fun"          # DGP function name
h_fun      <- get(h_fun_name)   # resolve function object

B        <- 1000                 # number of Monte Carlo iterations
n_cores  <- 62                   # parallel cores

sig_idx  <- 2                    # GPD scale covariate indices
shl_idx  <- 1:2                  # GPD shape covariate indices
n        <- 5000                 # sample size
R        <- 1e7                  # quasi-sample count for KDE
interval <- c(0.5, 300)          # optimisation search interval
tau      <- seq(0.0001, 0.9999, length.out = 2000)

# ==============================================================================
# Automatic path management
# ==============================================================================
file_tag   <- sprintf("%s_n%d_B%d", h_fun_name, n, B)
output_dir <- file.path("scratch", h_fun_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(sprintf("Created directory: %s", output_dir))
}
log_file  <- file.path(output_dir, paste0(file_tag, ".txt"))
save_file <- file.path(output_dir, paste0(file_tag, ".Rdata"))

# ==============================================================================
# Initialise log and compute true QTEs
# ==============================================================================
cat(paste0(">>> Simulation [", file_tag, "] started at ", Sys.time(), "\n"),
    file = log_file, append = FALSE)
cat(sprintf("Config: n=%d, B=%d, cores=%d\n", n, B, n_cores),
    file = log_file, append = TRUE)

pn_list  <- c(5 / n, 1 / n, 5 / (n * log(n)))
true_qte <- calculate_true_qtes(h_fun, pn_list)$qtes

# ==============================================================================
# Main parallel loop
# ==============================================================================
set.seed(12345)
RNGkind("L'Ecuyer-CMRG")

res_list <- parallel::mclapply(1:B, function(r) {

  if (r %% 50 == 0) {
    cat(sprintf("[Start] Iter %d on PID %d at %s\n", r, Sys.getpid(), Sys.time()),
        file = log_file, append = TRUE)
  }

  out <- tryCatch({
    simu <- h_fun(n)
    Y1 <- simu$Y1; Y0 <- simu$Y0
    D  <- simu$D;  X  <- simu$X
    Y  <- D * Y1 + (1 - D) * Y0

    mat <- vector("list", length(pn_list))
    for (j in seq_along(pn_list)) {
      pn    <- pn_list[j]
      k_seq <- n^(0.65)
      hill <- qte_extrapolation_hill(Y, X, D, pn, k_seq)
      pk   <- qte_extrapolation_pickands(Y, X, D, pn, k_seq)
      zg   <- qte_firpo_zhang(Y, X, D, pn, N_bootstrap = 1000)
      tiee <- tiee_cpp(Y, X, D, qn = 1 - pn, u = 1 - k_seq / n,
                       R = R, tau = tau,
                       sig_idx = sig_idx, shl_idx = shl_idx,
                       interval = interval)
      mat[[j]] <- list(hill = hill, pickland = pk, zhang = zg, tiee = tiee)
    }
    mat
  }, error = function(e) {
    cat(sprintf("[ERROR] Iter %d: %s\n", r, e$message),
        file = log_file, append = TRUE)
    return(NA)
  })

  if (r %% 50 == 0) {
    cat(sprintf("[Done ] Iter %d at %s\n", r, Sys.time()),
        file = log_file, append = TRUE)
  }

  out
}, mc.cores = n_cores, mc.set.seed = TRUE)

# ==============================================================================
# Save results
# ==============================================================================
save(res_list, true_qte, pn_list, n, B, file = save_file)
cat(sprintf("Results saved to: %s\n", save_file), file = log_file, append = TRUE)
