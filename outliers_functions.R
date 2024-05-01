library(fpp3)
library(forecast)
library(TSA)
library(showtext)

font_add("SimSun", "/usr/share/fonts/windows-fonts/simsun.ttc")
font_add("Times New Roman", "/usr/share/fonts/windows-fonts/times.ttf")
showtext_auto()
theme_set(theme_light())
theme_update(text=element_text(size=26,  family="SimSun"))
theme_update(axis.text=element_text(size=24,  family="Times New Roman"))


# General
# -------
x_as_tibble = function(x) {
  df = tibble(timestamp = 1:length(x), value = x)
  return(df)
}

generate_is_anomaly = function(x, outliers_ind) {
  #' e.g.
  #' x = c(5, 4, 3, 2, 1); outliers_ind = c(2, 4)
  #' generate_is_anomaly(x, outliers_ind) = c(0, 1, 0, 1, 0)
  is_anomaly = rep(0, length(x))
  is_anomaly[outliers_ind] = 1
  return(is_anomaly)
}

generate_outliers_ind = function(is_anomaly) {
  return(which(as.logical(is_anomaly)))
}

get_precision = function(is_anomaly_pred, is_anomaly_true) {
  #' the format of is_anomaly_pred/true need be c(1, 0, 0, 1, 0, 0, 0) like
  # TP / (TP + FP)
  outs_pred = generate_outliers_ind(is_anomaly_pred)
  outs_true = generate_outliers_ind(is_anomaly_true)
  in_outs = intersect(outs_true, outs_pred)
  TP = length(in_outs)
  # Note: please check the examples of 'setdiff' to understand it's behavior
  FP = length(setdiff(outs_pred, outs_true))
  if (length(in_outs) == 0) {
    return(0)
  } else {
    return(TP / (TP + FP))
  }
}

get_recall = function(is_anomaly_pred, is_anomaly_true) {
  # TP / (TP + FN)
  outs_pred = generate_outliers_ind(is_anomaly_pred)
  outs_true = generate_outliers_ind(is_anomaly_true)
  in_outs = intersect(outs_true, outs_pred)
  if (length(in_outs) == 0) {
    return(0)
  } else {
    return(length(in_outs) / length(outs_true))
  }
}

get_F1 = function(is_anomaly_pred, is_anomaly_true) {
  precision = get_precision(is_anomaly_pred, is_anomaly_true)
  recall = get_recall(is_anomaly_pred, is_anomaly_true)
  if (precision == 0 & recall == 0) {
    return(0)
  } else {
    F1 = 2 * ((precision * recall) / (precision + recall))
    return(F1)
  }
}

get_all_metrics = function(is_anomaly_pred, is_anomaly_true) {
  df = tibble(
    precision = get_precision(is_anomaly_pred, is_anomaly_true),
    recall = get_recall(is_anomaly_pred, is_anomaly_true),
    F1 = get_F1(is_anomaly_pred, is_anomaly_true)
  )
  return(df)
}

generate_x_with_outs = function(x, pollution_rate, outs_gen_func, return_df = TRUE) {
  #' generate x with ouliers according to pollution_rate and outs_gen_func
  x_len = length(x)
  outs_len = round(x_len * pollution_rate)
  x_outs_index = sort(sample(1:x_len, outs_len))
  
  x_with_outs = x
  x_with_outs[x_outs_index] = outs_gen_func(outs_len)
  
  if (return_df == TRUE) {
    is_anomaly = rep(0, x_len)
    is_anomaly[x_outs_index] = 1
    df = tibble(timestamp = 1:x_len,
                value = x_with_outs,
                is_anomaly = is_anomaly)
    return(df)
  } else {
    return(x_with_outs)
  }
}

auto_recognize_AO_IO = function(AO, IO) {
  #' if lambda_IO > lambda_AO, recognize it as IO; otherwise, AO
  # get intersect, pure AO, pure IO index
  inter_ind = intersect(AO$ind, IO$ind)
  final_AO_ind = setdiff(AO$ind, inter_ind)
  final_IO_ind = setdiff(IO$ind, inter_ind)
  # get final AO IO labels
  AO_lam_inter = AO$lambda2[which(AO$ind %in% inter_ind)]
  IO_lam_inter = IO$lambda1[which(IO$ind %in% inter_ind)]
  inter_labels = ifelse(IO_lam_inter > AO_lam_inter, "IO", "AO")
  AO_labels = rep("AO", length(final_AO_ind))
  IO_labels = rep("IO", length(final_IO_ind))
  # return a data frame
  AO_IO_map = tibble(
    ind = c(inter_ind, final_AO_ind, final_IO_ind),
    labels = c(inter_labels, AO_labels, IO_labels)
  )
  return(arrange(AO_IO_map, ind))
}

get_outs_tsa = function(x) {
  x_ts = ts(x, frequency = findfrequency(x))
  # TSA's method have no extra effects in SARIMA,
  # and auto fit a SARIMA will spend tooooo much more time.
  # So we decide to ban it to fit seansonal.
  m = auto.arima(x_ts, seasonal = FALSE)
  AO = detectAO(m)
  IO = detectIO(m)
  outs_ind = auto_recognize_AO_IO(AO, IO)$ind
  is_anomaly = generate_is_anomaly(x, outs_ind)
  return(is_anomaly)
}

get_outs_comb = function(x, return_resid=FALSE) {
  if (is.constant(x)) {
    is_anomaly_comb = rep(0, length(x))
    return(is_anomaly_comb)
  }
  x_ts = ts(x, frequency = findfrequency(x))
  comps = mstl(x_ts, robust = TRUE)
  rem = remainder(comps)
  # if recognized season, a seasonally adjusted series is computed
  detrend = x_ts - trendcycle(comps)
  strength = 1 - var(rem) / var(detrend)
  if (strength >= 0.6) {
    x_ts = seasadj(comps)
  }
  # trend is estimated by applying Friedman’s super smoother
  mod = supsmu(1:length(x_ts), x_ts)
  # get resid
  resid = x_ts - mod$y
  if (return_resid) {
    return(resid)
  } else {
    is_anomaly_comb = get_outs_tsa(resid)
    return(is_anomaly_comb)
  }
}

get_outs_fore = function(x) {
  x_ts = ts(x, frequency = findfrequency(x))
  outs_ind = tsoutliers(x_ts, iterate = 1)$index
  is_anomaly = generate_is_anomaly(x, outs_ind)
  return(is_anomaly)
}

get_outs_comb_pro = function(x) {
  if (is.constant(x)) {
    is_anomaly_comb_pro = rep(0, length(x))
    return(is_anomaly_comb_pro)
  }
  resid = get_outs_comb(x, return_resid=TRUE)
  is_anomaly_comb = get_outs_tsa(resid)
  # union outliers detected by IQR
  resid.q = quantile(resid, probs = c(0.25, 0.75))
  iqr = diff(resid.q)
  limits = resid.q + 3 * iqr * c(-1, 1)
  outs_ind = which((resid < limits[1]) | (resid > limits[2]))
  is_anomaly = is_anomaly_comb
  is_anomaly[outs_ind] = 1
  return(is_anomaly)
}

get_outs_all_methods = function(x) {
  is_anomaly_fore = get_outs_fore(x)
  df = tibble(
    is_anomaly_tsa = get_outs_tsa(x),
    is_anomaly_comb = get_outs_comb(x),
    is_anomaly_fore = get_outs_fore(x),
    is_anomaly_comb_pro = get_outs_comb_pro(x)
  )
  return(df)
}

test_data = function(df) {
  #' df need to contains three col: 
  #' timestamp, value, is_anomaly
  outs_df = get_outs_all_methods(df$value)
  func_temp = function(s) { return(get_all_metrics(outs_df[[s]], df$is_anomaly)) }
  df_metrics_methods = do.call(rbind, lapply(names(outs_df), func_temp))
  df_metrics_methods$method = c("Chang", "decomp+Chang", "decomp+IQR", "decomp+Chang+IQR")
  return(df_metrics_methods)
}

test_data_based_csv_path = function(csv_path) {
  cat("[INFO] Start handle", csv_path, "\n")
  df = readr::read_csv(csv_path)
  df_metrics_methods = test_data(df)
  df_metrics_methods$data = strsplit(basename(csv_path), "\\.")[[1]][1]
  return(df_metrics_methods)
}

test_data_based_csv_dir = function(csv_dir) {
  csv_files = list.files(csv_dir,
                         pattern = "\\.csv$",
                         full.names = TRUE)
  df_metrics_methods_multi = do.call(rbind, lapply(csv_files, test_data_based_csv_path))
  return(df_metrics_methods_multi)
}
