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

p_arma = ggplot(dfmm_multi_arma, aes(x = pollution_rate, y = F1, color = method)) +
  geom_point() +
  geom_line()

ggsave(paste0(asset_dir, "compare-arma.png"),
  plot = p_arma, width = 8, height = 4
)
#recalls = l$recalls
#for (i in seq(0.02, 0.2, 0.01)) {
#  new_row = test_data_synthetic_arima(i)$recalls
#  recalls = add_row(recalls, new_row)
#}
#recalls = filter(recalls, !if_any(starts_with("recall"), is.na))
#
#recalls_long = pivot_longer(
#  recalls,
#  cols = starts_with("recall"),
#  names_to = "method",
#  values_to = "recall"
#)
#
#ggplot(recalls_long, aes(x = pollution_rate, y = recall, color = method)) +
#  geom_line() +
#  theme_light()
#
#
## Plot test
## ---------
#l = test_data_synthetic_arima(pollution_rate = 0.05, n = 50)
#x_ori = l$x_ori
#x_pollu = l$x
#
## Origin x plot
#ggplot(x_as_tibble(x_ori), aes(x = timestamp, y = value)) + 
#  geom_point(size = 0.4) + 
#  geom_line(linewidth = 0.3) + 
#  theme_light()
#
## Polluted x plot
#ggplot(x_as_tibble(x_pollu), aes(x = timestamp, y = value)) + 
#  geom_point(size = 0.01) + 
#  geom_line()
#
#
