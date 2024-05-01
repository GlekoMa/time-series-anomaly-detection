setwd_lw <- function(src_dir) {
  if (getwd() == src_dir) {
    print("Yes")
  } else {
    setwd(src_dir)
  }
}
setwd_lw(src_dir = "/home/gleko/Projects/lw")
source("outliers_functions.R")
asset_dir <- "assets/"


# Functions
# ---------
generate_outs_unif <- function(n) {
  # this func is only for ARMA(1, 1) ar=0.65, ma=0.35 whose avg abs is 4
  n1 <- sample(1:n, 1)
  n2 <- n - n1
  c(runif(n1, -15, -5), runif(n2, 5, 15))
}


test_data_synthetic_arma = function(pollution_rate, n = 1000, only_return_dfmm = TRUE) {
  # Simulate ARMA(1, 1) ar=0.65, ma=0.35
  set.seed(123)
  arma_sim = arima.sim(list(
    order = c(1, 0, 1),
    ar = 0.65,
    ma = 0.35
  ), n = n)
  # Random replace with uniform
  df_with_outs = generate_x_with_outs(x = arma_sim,
                            pollution_rate = pollution_rate,
                            outs_gen_func = generate_outs_unif)
  # Test & return
  df_metrics_methods = test_data(df_with_outs)
  df_metrics_methods$pollution_rate = pollution_rate
  if (only_return_dfmm) {
    return(df_metrics_methods)
  } else {
    l = list(
      x_ori = as.numeric(arma_sim),
      x_with_out = df_with_outs$value,
      df_with_outs = df_with_outs,
      df_metrics_methods = df_metrics_methods
    )
    return(l)
  }
}


# 1. ARMA
# --------
dfmm_multi_arma = do.call(rbind, lapply(seq(0.01, 0.2, 0.01), test_data_synthetic_arma))

df_plot_arma = dfmm_multi_arma
df_plot_arma$method = case_match(
  df_plot_arma$method,
  "Chang" ~ "Chang",
  "decomp+Chang" ~ "分解+Chang",
  "decomp+Chang+IQR" ~ "分解+Chang+IQR",
  "decomp+IQR" ~ "分解+IQR"
)
df_plot_arma = rename(
  df_plot_arma, 
  污染率 = pollution_rate,
  算法 = method
)
p_arma = ggplot(df_plot_arma, aes(x = 污染率, y = F1, color = 算法)) +
  geom_point() +
  geom_line()
ggsave(paste0(asset_dir, "compare-arma.png"),
  plot = p_arma, width = 8, height = 4
)

