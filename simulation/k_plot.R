# Plot: MSE vs grid size K (Appendix B.3)
load("scratch/simu_k.Rdata")
sim_result <- res_k
res_array  <- sim_result$res_array
pn_list    <- sim_result$pn_list
true_qte   <- sim_result$true_qte

library(ggplot2)
library(reshape2)

long_df <- reshape2::melt(res_array, varnames = c("rep", "pn_idx", "Method")) %>%
  rename(estimate = value) %>%
  mutate(
    pn   = pn_list[as.integer(pn_idx)],
    true = true_qte[as.integer(pn_idx)],
    err  = estimate - true
  )

long_df %>%
  filter(!is.na(estimate)) %>%
  group_by(pn, Method) %>%
  summarise(
    Bias = mean(abs(err), na.rm = TRUE),
    SD   = sd(estimate, na.rm = TRUE),
    MSE  = Bias^2 + SD^2,
    n    = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = Method, y = MSE)) +
  geom_line(aes(group = factor(pn), color = factor(pn),
                linetype = factor(pn)), linewidth = 1) +
  scale_x_discrete(labels = seq(100, 2000, 100)) +
  scale_linetype_manual(
    name   = TeX("$1-\\tau_n$"),
    values = c("solid", "dashed", "dotted"),
    labels = c(TeX("$5/\\log(n)$"), TeX("$1/n$"), TeX("$5/n$"))
  ) +
  scale_color_discrete(
    name   = TeX("$1-\\tau_n$"),
    labels = c(TeX("$5/\\log(n)$"), TeX("$1/n$"), TeX("$5/n$"))
  ) +
  labs(x = TeX("Grid size $K$")) +
  theme_bw()

ggsave("scratch/sensitive_k.pdf", width = 8, height = 3)
