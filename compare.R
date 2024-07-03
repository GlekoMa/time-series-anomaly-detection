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

generate_synthetic_arma <- function(pollution_rate, n = 1000) {
  # Simulate ARMA(1, 1) ar=0.65, ma=0.35
  set.seed(123)
  arma_sim <- arima.sim(list(
    order = c(1, 0, 1),
    ar = 0.65,
    ma = 0.35
  ), n = n)
  # Random replace with uniform
  df_with_outs <- generate_x_with_outs(x = arma_sim,
                                       pollution_rate = pollution_rate,
                                       outs_gen_func = generate_outs_unif)
  return(list(arma_sim = arma_sim, df_with_outs = df_with_outs))
}

test_data_synthetic_arma <- function(pollution_rate,
                                     n = 1000,
                                     only_return_dfmm = TRUE) {
  l <- generate_synthetic_arma(pollution_rate = pollution_rate, n = n)
  arma_sim <- l$arma_sim
  df_with_outs <- l$df_with_outs
  # Test & return
  df_metrics_methods <- test_data(df_with_outs) |>
    mutate(pollution_rate = pollution_rate)
  if (only_return_dfmm) {
    return(df_metrics_methods)
  } else {
    l <- list(
      x_ori = as.numeric(arma_sim),
      x_with_out = df_with_outs$value,
      df_with_outs = df_with_outs,
      df_metrics_methods = df_metrics_methods
    )
    return(l)
  }
}

str_split_select <- function(s, split = "-", ind = 1) {
  ssp <- strsplit(s, split = split)
  v <- purrr::map_chr(ssp, ~ .x[ind])
  return(v)
}


# 1. ARMA
# --------
df_metrics_methods_multi_arma <- do.call(rbind, lapply(seq(0.01, 0.2, 0.01), test_data_synthetic_arma))

df_plot_arma <- df_metrics_methods_multi_arma |>
  mutate(
    method = case_match(
      method,
      "Chang" ~ "Chang",
      "decomp+Chang" ~ "分解+Chang",
      "decomp+Chang+IQR" ~ "分解+Chang+IQR",
      "decomp+IQR" ~ "分解+IQR"
    )
  ) |>
  mutate(method = factor(method, levels = c(
    "分解+Chang+IQR", "分解+IQR", "分解+Chang", "Chang"
  ))) |>
  rename(污染率 = pollution_rate, 算法 = method)

p_arma <- ggplot(df_plot_arma, aes(x = 污染率, y = F1, color = 算法)) +
  geom_point() +
  geom_line()
ggsave(
  paste0(asset_dir, "compare-arma.png"),
  plot = p_arma,
  width = 8,
  height = 4
)


# 2. GutenTAG
# -----------
guten_csv <- "df_metrics_methods_multi_guten.csv"
if (file.exists(guten_csv)) {
  df_metrics_methods_multi_guten <- readr::read_csv(guten_csv)
} else {
  # this process need to run too long time
  df_metrics_methods_multi_guten <- test_data_based_csv_dir("../data/lw/synthetic")
  readr::write_csv(x = df_metrics_methods_multi_guten, file = guten_csv)
}


# 3. NAB
# ------
nab_csv <- "df_metrics_methods_multi_nab.csv"
if (file.exists(nab_csv)) {
  df_metrics_methods_multi_nab <- readr::read_csv(nab_csv)
} else {
  # this process need to run too long time (but shorter than guten)
  df_metrics_methods_multi_nab <- test_data_based_csv_dir("../data/lw/real")
  readr::write_csv(x = df_metrics_methods_multi_nab, file = nab_csv)
}


# 4. Box plot
# -----------
df_plot_box <- rbind(
  rename(df_metrics_methods_multi_arma, data = pollution_rate),
  df_metrics_methods_multi_guten,
  df_metrics_methods_multi_nab
) |>
  mutate(
    method = case_match(
      method,
      "Chang" ~ "Chang",
      "decomp+Chang" ~ "分解+Chang",
      "decomp+Chang+IQR" ~ "分解+Chang+IQR",
      "decomp+IQR" ~ "分解+IQR"
    )
  ) |>
  mutate(factor(method, levels = c(
    "分解+Chang+IQR", "分解+IQR", "分解+Chang", "Chang"
  ))) |>
  rename(精确率 = precision, 召回率 = recall, 算法 = method)

p_box_precision <- ggplot(df_plot_box, aes(x = 算法, y = 精确率)) +
  geom_boxplot(size = 0.3, outlier.shape = 1) +
  coord_flip()
p_box_recall <- ggplot(df_plot_box, aes(x = 算法, y = 召回率)) +
  geom_boxplot(size = 0.3, outlier.shape = 1) +
  coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())
p_box_f1 <- ggplot(df_plot_box, aes(x = 算法, y = F1)) +
  geom_boxplot(size = 0.3, outlier.shape = 1) +
  coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())
p_box_all <- p_box_precision + p_box_recall + p_box_f1
ggsave(
  paste0(asset_dir, "compare-box.png"),
  plot = p_box_all,
  width = 8,
  height = 2
)


# 5. Compare efficiency: comb_pro VS tsoutliers (test in ARMA)
# ------------------------------------------------------------
time_it <- function(func, n = 100) {
  arma_with_outs <- generate_synthetic_arma(0.01, n = n)$df_with_outs$value
  timings <- rep(0, 5)
  for (i in 1:5) {
    timing_temp <- system.time(func(arma_with_outs))
    timings[i] <- as.numeric(timing_temp[3])
  }
  timing <- mean(timings)
  return(timing)
}

time_it_multi_by_n <- function(func, ns) {
  timings <- rep(0, length(ns))
  for (i in 1:length(ns)) {
    timings[i] <- time_it(func = func, n = ns[i])
  }
  df_timings <- tibble(n = ns, timing = timings)
  return(df_timings)
}

df_timings_tso <- time_it_multi_by_n(tsoutliers::tso, ns = seq(500, 2000, 500))
df_timings_pro <- time_it_multi_by_n(get_outs_comb_pro, ns = seq(500, 2000, 500))

df_timings_tso$算法 <- "tsoutliers"
df_timings_pro$算法 <- "ours"
df_plot_timings <- rbind(df_timings_tso, df_timings_pro) |>
  rename(时间 = timing)
p_timings <- ggplot(df_plot_timings, aes(x = n, y = 时间, color = 算法)) +
  geom_point() +
  geom_line()
ggsave(
  paste0(asset_dir, "compare-timings.png"),
  plot = p_timings,
  width = 5,
  height = 2
)
