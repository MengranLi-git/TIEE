# Application: extreme quantile treatment effects on river stations
r_files <- list.files(path = "functions", pattern = "\\.R$", full.names = TRUE)
lapply(r_files, source)

load("data/eqte.Rdata")

library(parallel)
library(actuar)
library(quantreg); library(ismev); library(evd)
library(dplyr); library(tidyr); library(BMisc)
library(lubridate)

name <- unique(data$Name)
pn  <- 1 - 0.99
res <- list()

# --- Run TIEE, Hill, and Zhang for each station ---
for (i in (19:length(name))[-18]) {
  data_river <- data %>% filter(Name == name[i])
  Y    <- data_river$P
  Year <- year(data_river$Date)
  D    <- as.numeric(Year >= 1980)
  X    <- data_river %>% dplyr::select(z500, msl_hpa, wind, T, NAO, AO, AMO, ENSO, PDO)

  # Propensity score model
  hn <- 4
  prop.fit    <- glm(D ~ poly(X[[1]], 2) + poly(X[[2]], 2) +
                       X[[3]] + X[[5]] + poly(X[[6]], hn),
                     family = binomial)
  prop_scores <- pmax(pmin(fitted(prop.fit), 0.95), 0.05)

  u   <- 1 - n^0.65 / n
  tau <- seq(0.001, 0.999, length.out = 1000)
  X   <- as.matrix(X[, 1:6])

  tiee <- tiee_cpp(Y, X, D, qn = 1 - pn, u, R = 1e6, tau,
                   interval = c(0.5, 200), sig_idx = 1:2, shl_idx = 1,
                   se = TRUE, prop_scores = prop_scores, heavy = FALSE)
  hill <- qte_extrapolation_hill(Y, X, D, pn, n^(0.65),
                                  prop_scores = prop_scores)
  zg   <- qte_firpo_zhang(Y, X, D, pn, N_bootstrap = 100,
                           prop_scores = prop_scores)
  res[[i]] <- list(hill, zg, tiee)
  print(i)
}

save.image(file = "scratch/res.Rdata")

# --- Extract estimates ---
Hill  <- data.frame(
  eqte = sapply(res, function(x) x[[1]]$qte),
  lower = unlist(sapply(res, function(x) x[[1]]$CIlower)),
  upper = unlist(sapply(res, function(x) x[[1]]$CIupper))
)
Zhang <- data.frame(
  eqte = sapply(res, function(x) x[[2]]$qte),
  lower = unlist(sapply(res, function(x) x[[2]]$CIlower)),
  upper = unlist(sapply(res, function(x) x[[2]]$CIupper))
)
TIEE <- data.frame(
  eqte = sapply(res, function(x) x[[3]]$eqte),
  lower = unlist(sapply(res, function(x) x[[3]]$lower)),
  upper = unlist(sapply(res, function(x) x[[3]]$upper))
)

# Empirical QTE
emp <- numeric(23)
for (i in (1:length(name))[-18]) {
  data_river <- data %>% filter(Name == name[i])
  Y    <- data_river$P
  Year <- year(data_river$Date)
  D    <- as.numeric(Year >= 1980)
  emp[i] <- quantile(Y[D == 1], 0.99) - quantile(Y[D == 0], 0.99)
}

df_raw <- cbind(cbind(cbind(Hill, Zhang), TIEE), emp[-18])

# --- Format results table ---
library(kableExtra)

df <- as.data.frame(df_raw)
colnames(df) <- c(
  "Hill_est", "Hill_low", "Hill_upp",
  "Zhang_est", "Zhang_low", "Zhang_upp",
  "TIEE_est",  "TIEE_low",  "TIEE_upp",
  "Empirical"
)

# Number formatter: align positive and negative values via phantom minus
fnum <- function(x) {
  txt <- sprintf("%.2f", x)
  ifelse(x >= 0, paste0("\\phantom{-}", txt), txt)
}

# Combine estimate [lower, upper]
fmt_ci <- function(est, low, upp) {
  paste0(fnum(est), " [", fnum(low), ", ", fnum(upp), "]")
}

tab_disp <- data.frame(
  Station   = name[-18],
  Empirical = fnum(df$Empirical),
  Hill      = fmt_ci(df$Hill_est, df$Hill_low, df$Hill_upp),
  Zhang     = fmt_ci(df$Zhang_est, df$Zhang_low, df$Zhang_upp),
  TIEE      = fmt_ci(df$TIEE_est, df$TIEE_low, df$TIEE_upp)
)

# Bold TIEE estimates where CI excludes zero
is_sig <- (df$TIEE_low > 0) | (df$TIEE_upp < 0)
tab_disp$TIEE[is_sig] <- paste0("\\textbf{", tab_disp$TIEE[is_sig], "}")

kbl(tab_disp, booktabs = TRUE, format = "latex", escape = FALSE,
    align = "lllll",
    caption = paste0("Comparison of EQTE estimates. ",
                     "\\textbf{Bold} TIEE values indicate significance.")) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))
