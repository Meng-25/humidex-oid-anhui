# ==============================================================================
# Environmental Pollution - Reproducible Code Package (R)
# Project: Humidex-related risk of Other Infectious Diarrhea (OID) in Anhui, China
# Purpose: Data cleaning -> City-specific DLNM -> Meta-analysis/meta-regression -> Burden -> Figures/Tables
#
# How to run (from project root):
#   Rscript analysis_EP_checkable.R
#
# Folder structure expected:
#   ./all_cities.xlsx (placed in project root)
#   ./outputs/   (auto-created)
#
# Notes for journal checking:
# - No absolute paths (no setwd("C:/...")).
# - All outputs written to ./outputs with clear filenames.
# - If raw data cannot be shared, keep a de-identified/demo dataset (same structure) for reviewers,
#   or provide a Data Availability statement + code that runs once data are available.
# ==============================================================================

# -------------------------
# 0) Global options & paths
# -------------------------
options(stringsAsFactors = FALSE)
set.seed(123)

paths <- list(
  data_dir = ".",
  out_dir  = "outputs",
  fig_dir  = file.path("outputs", "figures"),
  tab_dir  = file.path("outputs", "tables"),
  log_dir  = file.path("outputs", "logs")
)

dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$log_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(paths$log_dir, "run_log.txt")
sink(log_file, split = TRUE)
cat("Run started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# -------------------------
# 1) Packages (with checks)
# -------------------------
pkgs <- c(
  "readxl", "dplyr", "tidyr", "lubridate",
  "splines", "dlnm", "mvmeta", "MASS",
  "ggplot2", "gridExtra"
)

missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nPlease install them, e.g. install.packages(c(...))")
}
invisible(lapply(pkgs, library, character.only = TRUE))

# (Optional) used in your original script
if (requireNamespace("car", quietly = TRUE)) library(car)
if (requireNamespace("corrplot", quietly = TRUE)) library(corrplot)

cat("R version:\n"); print(R.version.string); cat("\n")
cat("Session info:\n"); print(sessionInfo()); cat("\n")

# -------------------------
# 2) Helper functions
# -------------------------
calc_manual_i2_safe <- function(q_val, df) {
  if (length(q_val) > 1) q_val <- q_val[1]
  if (length(df) > 1) df <- df[1]
  if (is.na(q_val) || q_val == 0) return(NA_real_)
  max(0, (q_val - df) / q_val) * 100
}

get_q_stat_safe <- function(model) {
  q_obj <- summary(model)$qstat
  if (is.null(q_obj)) return(list(Q = NA, df = NA, p = NA))
  list(Q = q_obj$Q[1], df = q_obj$df[1], p = q_obj$pvalue[1])
}

calc_joint_wald <- function(model, var_name) {
  all_coef <- as.numeric(coef(model))
  all_vcov <- vcov(model)
  term_indices <- grep(var_name, rownames(all_vcov))
  if (length(term_indices) == 0) return(list(stat = NA, p = NA, df = NA))

  b_sub <- all_coef[term_indices]
  v_sub <- all_vcov[term_indices, term_indices]
  wald_stat <- tryCatch(as.numeric(t(b_sub) %*% solve(v_sub) %*% b_sub),
                        error = function(e) NA_real_)
  if (is.na(wald_stat)) return(list(stat = NA, p = NA, df = NA))
  p_val <- pchisq(wald_stat, df = length(b_sub), lower.tail = FALSE)
  list(stat = wald_stat, p = p_val, df = length(b_sub))
}

calc_stats_vec <- function(x, label) {
  x <- na.omit(as.numeric(x))
  data.frame(
    Variable = label,
    Mean = round(mean(x), 2),
    SD = round(sd(x), 2),
    Min = round(min(x), 2),
    P25 = round(quantile(x, 0.25), 2),
    Median = round(median(x), 2),
    P75 = round(quantile(x, 0.75), 2),
    Max = round(max(x), 2)
  )
}

qaic_quasipoisson <- function(model) {
  # QAIC = -2 * logLik / phi + 2k  (commonly used approximation)
  ll <- sum(dpois(model$y, model$fitted.values, log = TRUE))
  phi <- summary(model)$dispersion
  k <- length(coef(model))
  -2 * ll / phi + 2 * k
}

# -------------------------
# 3) Read and clean data
# -------------------------
data_file <- file.path(paths$data_dir, "all_cities.xlsx")
if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file,
       "\nExpected: ./all_cities.xlsx (project root)")
}
cat("Reading:", data_file, "\n")

raw_data <- readxl::read_excel(data_file) %>%
  mutate(
    date = as.Date(date),
    DOW  = as.factor(DOW),
    region_char = as.character(region_ah),
    # force numeric
    longitude = as.numeric(as.character(longitude)),
    latitude  = as.numeric(as.character(latitude)),
    green     = as.numeric(as.character(green))
  ) %>%
  mutate(
    Region_Label = dplyr::case_when(
      region_char == "1" ~ "Wanbei (North)",
      region_char == "2" ~ "Wanzhong (Central)",
      region_char == "3" ~ "Wannan (South)",
      TRUE ~ NA_character_
    )
  )

stopifnot(all(c("citys", "cases", "Humidex") %in% names(raw_data)))

# City-level metadata (averages)
city_meta <- raw_data %>%
  group_by(citys) %>%
  summarise(
    GDP = mean(GDP, na.rm = TRUE),
    hospital_beds_per_capita = mean(hospital_beds_per_capita, na.rm = TRUE),
    sewage_pipeline_coverage = mean(sewage_pipeline_coverage, na.rm = TRUE),
    urbanization_rate = mean(urbanization_rate, na.rm = TRUE),
    density = mean(density, na.rm = TRUE),
    longitude = mean(longitude, na.rm = TRUE),
    latitude  = mean(latitude, na.rm = TRUE),
    green     = mean(green, na.rm = TRUE),
    region_ah = first(region_ah),
    Region_Label = first(Region_Label),
    .groups = "drop"
  ) %>%
  arrange(citys)

# Standardize names used downstream
if ("density" %in% names(city_meta)) city_meta <- rename(city_meta, population_density = density)
if ("green" %in% names(city_meta))   city_meta <- rename(city_meta, green_coverage_rate = green)

# SES index (PCA on GDP + urbanization_rate)
vars_for_pca <- city_meta %>% dplyr::select(GDP, urbanization_rate)
pca <- prcomp(vars_for_pca, scale. = TRUE)
pc1 <- pca$x[, 1]
# align direction so higher = higher GDP
city_meta$SES_Index <- if (cor(city_meta$GDP, pc1, use = "complete.obs") < 0) -pc1 else pc1

cat("Data loaded. Cities:", nrow(city_meta), " | Records:", nrow(raw_data), "\n\n")

# -------------------------
# 4) Stage 1: City-specific DLNM
# -------------------------
cities <- sort(unique(raw_data$citys))

# model parameters (match manuscript baseline)
lag_max <- 14
df_var  <- 4
df_lag  <- 4

coef_mat <- matrix(NA_real_, nrow = length(cities), ncol = df_var,
                   dimnames = list(cities, paste0("b", 1:df_var)))
vcov_list  <- vector("list", length(cities)); names(vcov_list) <- cities
argvar_list <- vector("list", length(cities)); names(argvar_list) <- cities

for (i in seq_along(cities)) {
  city <- cities[i]
  sub <- raw_data %>%
    filter(citys == city) %>%
    arrange(date) %>%
    mutate(
      log_lag1 = log(dplyr::lag(cases, 1) + 1),
      log_lag2 = log(dplyr::lag(cases, 2) + 1)
    )

  cen_val <- median(sub$Humidex, na.rm = TRUE)

  cb <- crossbasis(
    sub$Humidex, lag = lag_max,
    argvar = list(fun = "ns", df = df_var, cen = cen_val),
    arglag = list(fun = "ns", df = df_lag)
  )

  n_years <- length(unique(lubridate::year(sub$date)))

  m <- glm(
    cases ~ cb + ns(date, df = 7 * n_years) + DOW + log_lag1 + log_lag2,
    family = quasipoisson,
    data = sub,
    na.action = na.exclude
  )

  red <- crossreduce(cb, m, type = "overall", cen = cen_val)
  coef_mat[i, ] <- coef(red)
  vcov_list[[i]] <- vcov(red)
  argvar_list[[city]] <- attributes(cb)$argvar
}

# Align metadata order to model order (important)
city_meta <- city_meta %>% arrange(factor(citys, levels = cities))

cat("Stage 1 DLNM finished.\n\n")

# -------------------------
# 5) Stage 2: Meta-analysis and meta-regression
# -------------------------
meta_base_reml <- mvmeta(coef_mat ~ 1, vcov_list, method = "reml")
q_base <- get_q_stat_safe(meta_base_reml)
i2_base <- calc_manual_i2_safe(q_base$Q, q_base$df)

cat("Baseline heterogeneity:\n")
cat("  Q =", round(q_base$Q, 2), " df =", q_base$df, " p =", q_base$p, "\n")
cat("  I2 =", round(i2_base, 1), "%\n\n")

# Univariate meta-regression screening (same set as manuscript)
covariates <- c(
  "GDP", "urbanization_rate", "SES_Index",
  "population_density", "hospital_beds_per_capita", "sewage_pipeline_coverage",
  "latitude", "longitude", "green_coverage_rate"
)

meta_base_ml <- mvmeta(coef_mat ~ 1, vcov_list, method = "ml")

stats_table <- data.frame(
  Variable = "Intercept (Baseline)",
  LRT_Chi2 = NA, LRT_P = NA,
  Wald_Chi2 = NA, Wald_P = NA,
  Residual_Q = round(q_base$Q, 1),
  Residual_Q_df = q_base$df,
  I2_Estimated = round(i2_base, 1),
  Heterogeneity_Explained = "Ref",
  stringsAsFactors = FALSE
)

for (v in covariates) {
  if (!v %in% names(city_meta)) next

  f_uni <- as.formula(paste("coef_mat ~", v))
  m_reml <- tryCatch(mvmeta(f_uni, data = city_meta, vcov_list, method = "reml"),
                     error = function(e) NULL)
  if (is.null(m_reml)) next

  wald <- calc_joint_wald(m_reml, v)
  q_v  <- get_q_stat_safe(m_reml)
  i2_v <- calc_manual_i2_safe(q_v$Q, q_v$df)
  expl <- max(0, (q_base$Q - q_v$Q) / q_base$Q * 100)

  # LRT based on ML fits
  m_ml <- tryCatch(mvmeta(f_uni, data = city_meta, vcov_list, method = "ml"),
                   error = function(e) NULL)
  if (!is.null(m_ml)) {
    lrt <- 2 * (as.numeric(logLik(m_ml)) - as.numeric(logLik(meta_base_ml)))
    if (lrt < 0) lrt <- 0
    df_diff <- attr(logLik(m_ml), "df") - attr(logLik(meta_base_ml), "df")
    lrt_p <- pchisq(lrt, df = df_diff, lower.tail = FALSE)
  } else {
    lrt <- NA; lrt_p <- NA
  }

  stats_table <- bind_rows(stats_table, data.frame(
    Variable = v,
    LRT_Chi2 = round(lrt, 2),
    LRT_P = ifelse(is.na(lrt_p), NA, sprintf("%.4f", lrt_p)),
    Wald_Chi2 = round(wald$stat, 2),
    Wald_P = ifelse(is.na(wald$p), NA, sprintf("%.4f", wald$p)),
    Residual_Q = round(q_v$Q, 1),
    Residual_Q_df = q_v$df,
    I2_Estimated = round(i2_v, 1),
    Heterogeneity_Explained = sprintf("%.1f%%", expl),
    stringsAsFactors = FALSE
  ))
}

write.csv(stats_table, file.path(paths$tab_dir, "Table_MetaRegression_Univariate.csv"), row.names = FALSE)
cat("Saved:", file.path(paths$tab_dir, "Table_MetaRegression_Univariate.csv"), "\n\n")

# Final prediction meta-regression model (as used in manuscript for BLUPs)
meta_pred <- mvmeta(coef_mat ~ SES_Index + latitude, data = city_meta, vcov_list, method = "reml")

# BLUPs
blup <- mvmeta::blup(meta_pred, vcov = TRUE)
blup_coef <- matrix(NA_real_, nrow = length(cities), ncol = ncol(coef_mat), dimnames = list(cities, colnames(coef_mat)))
blup_vcov <- vector("list", length(cities)); names(blup_vcov) <- cities

for (i in seq_along(cities)) {
  obj <- blup[[i]]
  blup_coef[i, ] <- if (!is.null(obj$blup)) obj$blup else obj[[1]]
  blup_vcov[[i]] <- if (!is.null(obj$vcov)) obj$vcov else obj[[2]]
}

cat("Meta-regression + BLUPs ready.\n\n")

# -------------------------
# 6) Burden calculation (Warm season, Humidex > MMT)
# -------------------------
warm_months <- 5:10
n_sim <- 1000
x_mmt <- seq(min(raw_data$Humidex, na.rm = TRUE),
             max(raw_data$Humidex, na.rm = TRUE),
             length.out = 200)

# MMT per city based on BLUP curve
mmt <- setNames(rep(NA_real_, length(cities)), cities)
for (i in seq_along(cities)) {
  city <- cities[i]
  b <- do.call(onebasis, c(list(x = x_mmt), argvar_list[[city]]))
  log_rr <- b %*% blup_coef[i, ]
  mmt[city] <- x_mmt[which.min(log_rr)]
}

burden <- data.frame()

for (i in seq_along(cities)) {
  city <- cities[i]
  sub <- raw_data %>% filter(citys == city) %>% arrange(date)

  mmt_city <- mmt[city]
  target <- (month(sub$date) %in% warm_months) & (sub$Humidex > mmt_city)

  if (sum(target, na.rm = TRUE) == 0) next

  argvar <- argvar_list[[city]]
  b_city <- do.call(onebasis, c(list(x = sub$Humidex), argvar))
  b_ref  <- do.call(onebasis, c(list(x = mmt_city), argvar))
  b_cen  <- sweep(b_city, 2, as.numeric(b_ref), "-")

  # simulate RR
  sim_b <- MASS::mvrnorm(n_sim, mu = blup_coef[i, ], Sigma = blup_vcov[[city]])
  log_rr_sim <- b_cen %*% t(sim_b)
  rr_sim <- exp(log_rr_sim)

  rr_target <- rr_sim[target, , drop = FALSE]
  cases_target <- sub$cases[target]

  calc_an <- function(rr_vec) {
    af <- (rr_vec - 1) / rr_vec
    af[af < 0] <- 0
    sum(af * cases_target, na.rm = TRUE)
  }

  an_dist <- apply(rr_target, 2, calc_an)

  total_cases <- sum(sub$cases, na.rm = TRUE)
  est <- mean(an_dist)
  ci  <- quantile(an_dist, c(0.025, 0.975))

  burden <- bind_rows(burden, data.frame(
    City = city,
    Total_Cases = total_cases,
    MMT = round(mmt_city, 1),
    AN = round(est, 1),
    AN_L = round(ci[1], 1),
    AN_H = round(ci[2], 1),
    AF = round(est / total_cases * 100, 3),
    AF_L = round(ci[1] / total_cases * 100, 3),
    AF_H = round(ci[2] / total_cases * 100, 3)
  ))
}

burden <- burden %>%
  left_join(city_meta, by = c("City" = "citys")) %>%
  arrange(desc(AF))

write.csv(burden, file.path(paths$tab_dir, "Table_Burden_WarmSeason.csv"), row.names = FALSE)
cat("Saved:", file.path(paths$tab_dir, "Table_Burden_WarmSeason.csv"), "\n\n")

# -------------------------
# 7) Descriptive table (Table 1)
# -------------------------
t1 <- bind_rows(
  data.frame(Variable = "--- Health Outcome (Daily) ---", Mean="", SD="", Min="", P25="", Median="", P75="", Max=""),
  calc_stats_vec(raw_data$cases, "Daily OID cases"),
  data.frame(Variable = "--- Exposure (Daily) ---", Mean="", SD="", Min="", P25="", Median="", P75="", Max=""),
  calc_stats_vec(raw_data$Humidex, "Humidex"),
  data.frame(Variable = "--- City-level (N=16) ---", Mean="", SD="", Min="", P25="", Median="", P75="", Max=""),
  bind_rows(
    calc_stats_vec(city_meta$GDP, "GDP"),
    calc_stats_vec(city_meta$urbanization_rate, "Urbanization rate"),
    calc_stats_vec(city_meta$SES_Index, "SES index (PCA)"),
    calc_stats_vec(city_meta$population_density, "Population density"),
    calc_stats_vec(city_meta$hospital_beds_per_capita, "Hospital beds per capita"),
    calc_stats_vec(city_meta$sewage_pipeline_coverage, "Sewage pipeline coverage"),
    calc_stats_vec(city_meta$green_coverage_rate, "Green coverage rate"),
    calc_stats_vec(city_meta$latitude, "Latitude"),
    calc_stats_vec(city_meta$longitude, "Longitude")
  )
)

write.csv(t1, file.path(paths$tab_dir, "Table1_Descriptive.csv"), row.names = FALSE)
cat("Saved:", file.path(paths$tab_dir, "Table1_Descriptive.csv"), "\n\n")

# -------------------------
# 8) Figures (core set)
# -------------------------

# Figure 1: Overall pooled exposure-response curve (baseline meta model)
cen_prov <- median(raw_data$Humidex, na.rm = TRUE)
x_all <- seq(min(raw_data$Humidex, na.rm = TRUE), max(raw_data$Humidex, na.rm = TRUE), length.out = 100)

b_all <- onebasis(x_all, fun = "ns", df = df_var)
log_rr <- b_all %*% coef(meta_base_reml)
rr_raw <- exp(log_rr)
# normalize at cen_prov
ref_i <- which.min(abs(x_all - cen_prov))
rr <- rr_raw / rr_raw[ref_i]

b_ref <- b_all[ref_i, , drop = FALSE]
b_cen <- sweep(b_all, 2, as.numeric(b_ref), "-")
se <- sqrt(diag(b_cen %*% vcov(meta_base_reml) %*% t(b_cen)))
df_fig1 <- data.frame(
  Humidex = x_all,
  RR = rr,
  Low = exp(log(rr) - 1.96 * se),
  High = exp(log(rr) + 1.96 * se)
)

p1 <- ggplot(df_fig1, aes(Humidex, RR)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cen_prov, linetype = "dotted") +
  geom_ribbon(aes(ymin = Low, ymax = High), alpha = 0.2) +
  geom_line(size = 1) +
  theme_bw() +
  labs(x = "Humidex", y = "Relative Risk", title = "Figure 1. Overall association (pooled)")

ggsave(file.path(paths$fig_dir, "Figure1_Overall.png"), p1, width = 7, height = 5, dpi = 300)

# Figure 5-like: Forest of AF
p_forest <- ggplot(burden, aes(y = reorder(City, AF), x = AF, color = Region_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = AF_L, xmax = AF_H), height = 0.3) +
  geom_point(aes(size = Total_Cases), alpha = 0.9) +
  theme_bw() +
  labs(x = "Attributable Fraction (%)", y = "", title = "Warm season burden by city") +
  theme(legend.position = "bottom")

ggsave(file.path(paths$fig_dir, "Figure_Forest_AF.png"), p_forest, width = 8, height = 6, dpi = 300)

cat("Figures saved in:", paths$fig_dir, "\n\n")

# -------------------------
# 9) (Optional) Sensitivity & subgroup analyses
# -------------------------
# For journal checking, you can keep extended analyses in separate scripts:
#   - analysis_sensitivity.R
#   - analysis_subgroup.R
# to keep the main pipeline short and easy to run.

cat("Run finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()
