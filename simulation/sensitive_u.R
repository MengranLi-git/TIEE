# Sensitivity analysis: threshold u (paper Appendix B.2)
library(quantreg); library(ismev); library(actuar); library(evd)
library(dplyr); library(tidyr); library(parallel); library(BMisc)

# Source helper functions
r_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
lapply(r_files, source)

run_simulation <- function(h_fun, n = 1000, R = 1000, n_cores = 20,
                           log_file = "scratch/u_log.txt") {

  pn_list  <- c(5 / n, 1 / n, 5 / (n * log(n)))
  true_qte <- calculate_true_qtes(h_fun, pn_list)$qtes
  u_list   <- seq(0.85, 0.95, 0.01)

  res_list <- parallel::mclapply(1:R, function(r) {
    library(quantreg); library(ismev); library(actuar); library(evd)
    tryCatch({
      simu <- h_fun(n)
      Y1 <- simu$Y1; Y0 <- simu$Y0
      D  <- simu$D;  X  <- simu$X
      Y  <- D * Y1 + (1 - D) * Y0

      mat <- matrix(NA, nrow = length(pn_list), ncol = length(u_list),
                    dimnames = list(NULL, paste0("u_", u_list)))

      for (j in seq_along(pn_list)) {
        pn <- pn_list[j]
        for (i in seq_along(u_list)) {
          ipw <- tiee_cpp(Y, X, D, qn = 1 - pn, u = u_list[i], R = 1e6,
                          tau = seq(0.001, 0.999, length.out = 1000),
                          interval = c(0.5, 100), se = FALSE)
          mat[j, i] <- ipw$eqte
        }
      }

      if (r %% 100 == 0) {
        cat(sprintf("Progress: %d / %d at %s\n", r, R, Sys.time()),
            file = log_file, append = TRUE)
      }
      mat
    }, error = function(e) {
      cat(sprintf("ERROR at iteration %d: %s\n", r, e$message),
          file = log_file, append = TRUE)
      matrix(NA, nrow = length(pn_list), ncol = length(u_list))
    })
  }, mc.cores = n_cores)

  res_array <- array(
    unlist(res_list),
    dim = c(length(pn_list), length(u_list), R),
    dimnames = list(
      paste0("pn=", round(pn_list, length(u_list))),
      paste0("u_", u_list), NULL
    )
  )
  res_array <- aperm(res_array, c(3, 1, 2))

  cat(sprintf("Simulation finished at %s\n", Sys.time()),
      file = log_file, append = TRUE)

  list(res_array = res_array, pn_list = pn_list, true_qte = true_qte)
}

res_u <- run_simulation(h3.fun)
save.image("scratch/simu_u.Rdata")
