# ==============================================================================
# Post-processing: compute metrics, generate tables and coverage plots
# ==============================================================================
library(dplyr); library(tidyr); library(purrr)
library(ggplot2); library(kableExtra)

# --- Shared helper functions ---------------------------------------------------

extract_values <- function(model_obj, model_name) {
  if (all(is.na(model_obj))) return(c(est = NA, lower = NA, upper = NA))
  tryCatch({
    if (model_name == "hill") {
      c(est = as.numeric(model_obj$qtes),
        lower = as.numeric(model_obj$CIlower),
        upper = as.numeric(model_obj$CIupper))
    } else if (model_name == "pickland") {
      c(est = as.numeric(model_obj$qte), lower = NA, upper = NA)
    } else if (model_name == "zhang") {
      c(est = as.numeric(model_obj$qte),
        lower = as.numeric(model_obj$CIlower),
        upper = as.numeric(model_obj$CIupper))
    } else if (model_name == "tiee") {
      c(est = as.numeric(model_obj$eqte),
        lower = as.numeric(model_obj$lower),
        upper = as.numeric(model_obj$upper))
    } else {
      c(est = NA, lower = NA, upper = NA)
    }
  }, error = function(e) c(est = NA, lower = NA, upper = NA))
}

# Convert raw res_list to a tidy data frame
flatten_results <- function(res_list, true_qte, pn_list) {
  map_dfr(seq_along(res_list), function(i) {
    iter_res <- res_list[[i]]
    if (is.null(iter_res) || all(is.na(iter_res))) return(NULL)
    map_dfr(seq_along(iter_res), function(j) {
      pn_res <- iter_res[[j]]
      map_dfr(names(pn_res), function(m_name) {
        vals <- extract_values(pn_res[[m_name]], m_name)
        data.frame(iter = i, pn_idx = j, pn_val = pn_list[j],
                   model = m_name, true_qte = true_qte[j],
                   est = vals["est"], lower = vals["lower"], upper = vals["upper"],
                   row.names = NULL)
      })
    })
  })
}

# Compute summary metrics from raw df
compute_metrics <- function(df) {
  df %>%
    filter(!is.na(est)) %>%
    group_by(pn_idx, model) %>%
    summarise(
      pn_value  = mean(pn_val),
      true_val  = mean(true_qte),
      Bias      = mean(est) - mean(true_qte),
      SD        = sd(est),
      MSE       = mean((est - true_qte)^2),
      Coverage  = mean(lower <= true_qte & upper >= true_qte, na.rm = TRUE),
      n_valid   = n(),
      .groups   = "drop"
    ) %>%
    arrange(pn_idx, model)
}

# --- Main pipeline: load 3 scenario files → metrics → table + plot -----------
#
# Usage example (run from the simulation/ directory):
#   source("analyze.R")
#   run_analysis(
#     files   = c("scratch/h1.fun_n1000_B1000.Rdata",
#                 "scratch/h2.fun_n1000_B1000.Rdata",
#                 "scratch/h3.fun_n1000_B1000.Rdata"),
#     tag     = "heavy1000",
#     tail    = "heavy",
#     scenario_labels = c('M[1]^"(H)"', 'M[2]^"(H)"', 'M[3]^"(H)"')
#   )

run_analysis <- function(files, tag, tail = c("heavy", "light"),
                         scenario_labels = NULL) {
  tail <- match.arg(tail)

  # 1. Load and flatten each scenario
  dfs <- lapply(files, function(f) {
    env <- new.env()
    load(f, envir = env)
    df <- flatten_results(env$res_list, env$true_qte, env$pn_list)
    compute_metrics(df)
  })

  # 2. Build table (Bias + MSE) — matches paper Table 1 / Appendix B.1
  process_scenario <- function(df, suffix) {
    df %>%
      ungroup() %>%
      select(pn_idx, model, Bias, MSE) %>%
      mutate(
        Bias_str = formatC(Bias, format = "f", digits = 3),
        MSE_str  = formatC(MSE,  format = "f", digits = 3),
        MSE_num  = MSE
      ) %>%
      select(pn_idx, model, Bias_str, MSE_str, MSE_num) %>%
      rename_with(~ paste0(., "_", suffix), -c(pn_idx, model))
  }

  d1 <- process_scenario(dfs[[1]], "S1")
  d2 <- process_scenario(dfs[[2]], "S2")
  d3 <- process_scenario(dfs[[3]], "S3")

  final_df <- d1 %>%
    left_join(d2, by = c("pn_idx", "model")) %>%
    left_join(d3, by = c("pn_idx", "model"))

  model_map <- c(hill = "Causal Hill", pickland = "Causal Pickands",
                 tiee = "TIEE-IPW", zhang = "Zhang-Firpo")
  pn_map <- c("3" = "$5/(n \\log n)$", "2" = "$1/n$", "1" = "$5/n$")

  df_table <- final_df %>%
    group_by(pn_idx) %>%
    mutate(
      MSE_str_S1 = ifelse(MSE_num_S1 == min(MSE_num_S1),
                          paste0("\\textbf{", MSE_str_S1, "}"), MSE_str_S1),
      MSE_str_S2 = ifelse(MSE_num_S2 == min(MSE_num_S2),
                          paste0("\\textbf{", MSE_str_S2, "}"), MSE_str_S2),
      MSE_str_S3 = ifelse(MSE_num_S3 == min(MSE_num_S3),
                          paste0("\\textbf{", MSE_str_S3, "}"), MSE_str_S3)
    ) %>%
    ungroup() %>%
    mutate(
      pn_label = pn_map[as.character(pn_idx)],
      model    = model_map[model]
    ) %>%
    mutate(pn_label = factor(pn_label,
                             levels = c("$5/(n \\log n)$", "$1/n$", "$5/n$"))) %>%
    arrange(pn_label, model) %>%
    select(pn_label, model,
           Bias_str_S1, MSE_str_S1,
           Bias_str_S2, MSE_str_S2,
           Bias_str_S3, MSE_str_S3)

  tail_label <- if (tail == "heavy") "Heavy" else "Light"
  caption <- sprintf("Estimation performance (Bias and MSE) under %s-tailed scenarios.", tolower(tail_label))

  kbl(df_table, format = "latex", booktabs = TRUE, escape = FALSE,
      col.names = c("$1-\\tau_n$", "Method",
                    "Bias", "MSE", "Bias", "MSE", "Bias", "MSE"),
      align = c("l", "l", rep("r", 6)),
      caption = caption) %>%
    add_header_above(c(" " = 2,
                       "Scenario $M_1$" = 2, "Scenario $M_2$" = 2, "Scenario $M_3$" = 2)) %>%
    collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    save_kable(file = sprintf("%s_table.tex", tag))

  cat(sprintf("Table saved: %s_table.tex\n", tag))

  # 3. Coverage plot
  if (is.null(scenario_labels)) {
    if (tail == "heavy") {
      scenario_labels <- c('M[1]^"(H)"', 'M[2]^"(H)"', 'M[3]^"(H)"')
    } else {
      scenario_labels <- c('M[1]^"(L)"', 'M[2]^"(L)"', 'M[3]^"(L)"')
    }
  }

  plot_data <- bind_rows(
    dfs[[1]] %>% mutate(scenario = scenario_labels[1]),
    dfs[[2]] %>% mutate(scenario = scenario_labels[2]),
    dfs[[3]] %>% mutate(scenario = scenario_labels[3])
  ) %>%
    filter(model != "pickland") %>%
    mutate(
      model_label = case_when(
        model == "zhang" ~ "Zhang",
        model == "hill"  ~ "Hill",
        model == "tiee"  ~ "TIEE"
      ),
      scenario = factor(scenario, levels = scenario_labels),
      pn_label = case_when(
        pn_idx == 1 ~ "5/n",
        pn_idx == 2 ~ "1/n",
        pn_idx == 3 ~ "5/(n~log~n)"
      ),
      pn_label = factor(pn_label, levels = c("5/n", "1/n", "5/(n~log~n)"))
    ) %>%
    mutate(
      se   = sqrt(Coverage * (1 - Coverage) / n_valid),
      ymin = pmax(0, Coverage - 1.96 * se),
      ymax = pmin(1, Coverage + 1.96 * se)
    )

  ggplot(plot_data, aes(x = pn_label, y = Coverage,
                        color = model_label, shape = model_label)) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.2, position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    facet_wrap(~ scenario, nrow = 1, labeller = label_parsed) +
    scale_color_manual(values = c(Zhang = "#D55E00", Hill = "#0072B2",
                                  TIEE = "#009E73")) +
    scale_x_discrete(labels = scales::parse_format()) +
    scale_y_continuous(limits = c(0.4, 1.05), breaks = seq(0.4, 1.0, 0.2)) +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_text(face = "bold"),
          axis.text.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold", size = 11)) +
    labs(x = TeX("Quantile level $1-\\tau_n$"),
         y = "Coverage Probability",
         color = "Method", shape = "Method")

  ggsave(sprintf("%s_coverage.pdf", tag), height = 4, width = 10)
  cat(sprintf("Plot saved: %s_coverage.pdf\n", tag))

  invisible(df_table)
}

# ==============================================================================
# Run all four analyses
# ==============================================================================

# Heavy-tailed, n = 1000
run_analysis(
  files = c("scratch/h1.fun_n1000_B1000.Rdata",
            "scratch/h2.fun_n1000_B1000.Rdata",
            "scratch/h3.fun_n1000_B1000.Rdata"),
  tag = "heavy1000", tail = "heavy"
)

# Heavy-tailed, n = 5000
run_analysis(
  files = c("scratch/h1.fun_n5000_B1000.Rdata",
            "scratch/h2.fun_n5000_B1000.Rdata",
            "scratch/h3.fun_n5000_B1000.Rdata"),
  tag = "heavy5000", tail = "heavy"
)

# Light-tailed, n = 1000
run_analysis(
  files = c("scratch/h1.light.fun_n1000_B1000.Rdata",
            "scratch/h2.light.fun_n1000_B1000.Rdata",
            "scratch/h3.light.fun_n1000_B1000.Rdata"),
  tag = "light1000", tail = "light"
)

# Light-tailed, n = 5000
run_analysis(
  files = c("scratch/h1.light.fun_n5000_B1000.Rdata",
            "scratch/h2.light.fun_n5000_B1000.Rdata",
            "scratch/h3.light.fun_n5000_B1000.Rdata"),
  tag = "light5000", tail = "light"
)

cat("All analyses complete.\n")
