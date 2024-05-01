setwd_lw <- function(src_dir) {
  if (getwd() == "/home/gleko/lw/src") {
    print("Yes")
  } else {
    setwd("/home/gleko/lw/src")
  }
}
setwd_lw(src_dir = "/home/gleko/lw/src")
source("outliers_functions.R")
asset_dir <- "assets/"
set.seed(30)


# Functions
# ---------
generate_outs_unif <- function(n) {
  n1 <- sample(1:n, 1)
  n2 <- n - n1
  c(runif(n1, 15000, 17000), runif(n2, 38000, 40000))
}


# 1. Create formated data frame (original & polluted data)
# --------------------------------------------------------
start_date <- as.POSIXct("2000-6-4 00:30:00")
end_date <- as.POSIXct("2000-8-27 00:00:00")
date_sequence <- seq(start_date, end_date, by = "30 mins")
# note: taylor has 12 weeks, but we only use first 3 weeks
df_ori <- tibble(
  timestamp = date_sequence,
  value = as.numeric(taylor),
  is_anomaly = rep(0, length(taylor))
)[1:1008, ]

df_with_outs <- generate_x_with_outs(df_ori$value, 0.025, generate_outs_unif)

df_with_outs$timestamp <- df_ori$timestamp


# 2. Plot data: origin vs pollution
# ---------------------------------
df_plot_ori_pollu <- rbind(df_ori, df_with_outs)
df_plot_ori_pollu <- rename(df_plot_ori_pollu,
  时间 = timestamp,
  `电力需求(MW)` = value,
  值类型 = is_anomaly
)
df_plot_ori_pollu <- mutate(df_plot_ori_pollu,
  值类型 = ifelse(值类型 == 0, "正常", "异常"),
  序列类型 = c(rep("原始序列", 1008), rep("异常序列", 1008))
)
p_ori_pollu <- ggplot(df_plot_ori_pollu, aes(x = 时间, y = `电力需求(MW)`)) +
  geom_point(aes(color = 值类型), size = 0.6) +
  geom_line(linewidth = 0.4) + 
  scale_color_manual(values = c("正常" = "black", "异常" = "orange")) +
  facet_grid(rows = vars(序列类型))
ggsave(paste0(asset_dir, "example-taylor-ori-pollu.png"),
  plot = p_ori_pollu, width = 8, height = 4
)


# 3. Plot the decomposition of polluted data
# ------------------------------------------
x_msts <- msts(df_with_outs$value,
  seasonal.periods = c(2 * 24, 2 * 24 * 7),
  start = 2000 + 22 / 52
)
comps <- mstl(x_msts, robust = TRUE)

comps_plot <- comps
colnames(comps_plot) <- c("异常数据", "趋势项", "季节项(天)", "季节项(周)", "剩余项")

p_comps <- autoplot(comps_plot) +
  labs(x = "时间", y = "电力需求(MW)")
ggsave(paste0(asset_dir, "example-taylor-comps.png"),
  plot = p_comps, width = 8, height = 4
)


# 4. Get the decomposition process value
# --------------------------------------
# get remainder
rem <- remainder(comps)
# get x_adj (by seasadj: delete seasons of x)
detrend <- x_msts - trendcycle(comps)
strength <- 1 - var(rem) / var(detrend)
if (strength >= 0.6) {
  # the condition is TRUE in this example
  x_adj <- seasadj(comps)
} else {
  # no matter whether seasadj has been performed or not, we name x as x_adj.
  x_adj <- x_msts
}
# get smooth x_adj
mod <- supsmu(1:length(x_adj), x_adj)
# get residual
resid <- x_adj - mod$y


# 5. Plot decomposition process: value vs seasadj vs supsmu
# ---------------------------------------------------------
df_plot_val_seas_sups <- df_with_outs
df_plot_val_seas_sups$季节调整后 <- x_adj
df_plot_val_seas_sups$季节调整并平滑后 <- mod$y
df_plot_val_seas_sups <- rename(df_plot_val_seas_sups,
  时间 = timestamp,
  原始 = value
)
df_plot_val_seas_sups <- select(df_plot_val_seas_sups, -is_anomaly)

df_plot_val_seas_sups_long <- pivot_longer(
  df_plot_val_seas_sups,
  cols = c("原始", "季节调整后", "季节调整并平滑后"),
  names_to = "分解中间值",
  values_to = "电力需求(MW)"
)

p_val_seas_sups <- ggplot(
  df_plot_val_seas_sups_long,
  aes(x = 时间, y = `电力需求(MW)`, color = `分解中间值`)
) +
  geom_line() +
  scale_color_manual(values = c(
    "原始" = "#CCCCCC",
    "季节调整后" = "#888888",
    "季节调整并平滑后" = "#222222"
  ))
ggsave(paste0(asset_dir, "example-taylor-val-seas-sups.png"),
  plot = p_val_seas_sups, width = 8, height = 4
)


# 6. Plot decomposition process: remainder vs residual
# ----------------------------------------------------
df_plot_residual <- df_with_outs
df_plot_residual$剩余项 <- rem
df_plot_residual$重估后的剩余项 <- resid
df_plot_residual <- select(df_plot_residual, -value, -is_anomaly)
df_plot_residual <- rename(df_plot_residual, 时间 = timestamp)

df_plot_residual_long <- pivot_longer(
  df_plot_residual,
  cols = c("剩余项", "重估后的剩余项"),
  names_to = "剩余/重估后的剩余项",
  values_to = "电力需求(MW)"
)

p_rem_resid <- ggplot(
  df_plot_residual_long,
  aes(x = 时间, y = `电力需求(MW)`, color = `剩余/重估后的剩余项`)
) +
  geom_line() +
  labs(x = "时间", y = "电力需求(MW)") +
  scale_color_manual(values = c(
    "剩余项" = "#C0C0C0",
    "重估后的剩余项" = "#555555"
  ))
ggsave(paste0(asset_dir, "example-taylor-rem-resid.png"),
  plot = p_rem_resid, width = 8, height = 4
)


# 7. Detect outliers using TSA, IQR and union method
# --------------------------------------------------
# TSA method
m <- auto.arima(ts(resid))
IO <- detectIO(m)
AO <- detectAO(m)
outs_ind_TSA <- auto_recognize_AO_IO(AO, IO)$ind
# IQR method
resid.q <- quantile(resid, probs = c(0.25, 0.75))
iqr <- diff(resid.q)
limits <- resid.q + 3 * iqr * c(-1, 1)
outs_ind_IQR <- which((resid < limits[1]) | (resid > limits[2]))
# Union method
outs_ind_all <- union(outs_ind_TSA, outs_ind_IQR)


# 8. Create 'example_taylor_IO_AO_df.csv' to compare lambda
# ---------------------------------------------------------
df_IO <- tibble(lambda1 = IO$lambda1, ind = IO$ind)
df_AO <- tibble(lambda2 = AO$lambda2, ind = AO$ind)
df_IO_AO <- full_join(df_IO, df_AO, by = join_by(ind == ind))
df_IO_AO$labels <- auto_recognize_AO_IO(AO, IO)$labels
df_IO_AO <- select(arrange(df_IO_AO, ind), ind, lambda1, lambda2, labels)
readr::write_csv(df_IO_AO, file = paste0(asset_dir, "example_taylor_IO_AO_df.csv"))


# 9. Get metrics: Precision, Recall, F1
# -------------------------------------
is_anomaly_pred <- generate_is_anomaly(df_with_outs$value, outs_ind_all)
is_anomaly_true <- df_with_outs$is_anomaly
df_metrics <- get_all_metrics(is_anomaly_pred, is_anomaly_true)
print(df_metrics)


# 10. Plot data: pollu vs TSA vs IQR vs union
# ------------------------------------------
df_plot_outs_labeled <- df_with_outs
df_plot_outs_labeled$is_anomaly <- as.factor(df_plot_outs_labeled$is_anomaly)
df_plot_outs_labeled$异常检测来源 <- rep("-", nrow(df_plot_outs_labeled))
df_plot_outs_labeled$异常检测来源[outs_ind_TSA] <- "Chang"
df_plot_outs_labeled$异常检测来源[outs_ind_IQR] <- "IQR"
df_plot_outs_labeled$异常检测来源[intersect(outs_ind_TSA, outs_ind_IQR)] <- "Chang & IQR"
df_plot_outs_labeled <- rename(df_plot_outs_labeled,
  时间 = timestamp,
  `电力需求(MW)` = value,
  值类型 = is_anomaly
)
df_plot_outs_labeled <- mutate(df_plot_outs_labeled, 值类型 = ifelse(值类型 == 0, "正常", "异常"))

p_TSA_IQR <- ggplot(df_plot_outs_labeled, aes(x = 时间, y = `电力需求(MW)`)) +
  geom_point(mapping = aes(color = 值类型, shape = 异常检测来源), size = 2) +
  geom_line(linewidth = 0.4) +
  scale_color_manual(values = c("正常" = "black", "异常" = "orange")) +
  scale_shape_manual(values = c("-" = 20, "Chang" = 0, "IQR" = 2, "Chang & IQR" = 14)) +
  labs(y = "电力需求(MW)")
ggsave(paste0(asset_dir, "example-taylor-TSA-IQR.png"),
  plot = p_TSA_IQR, width = 8, height = 4
)
