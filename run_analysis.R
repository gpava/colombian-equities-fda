# Functional analysis of trading volume and returns in Colombian equities.
# Run from the repository root: Rscript --vanilla run_analysis.R
# Add --export-tables to write numerical summaries locally.
analysis_complete <- FALSE
company_files <- c(bancolombia = "Bancolombia.csv", bogota = "Bancodebogota.csv", bvc = "Bolsavaloresdecolombia.csv",
    ceargos = "CementosArgos.csv", celsia = "CelsiaSA.csv", ecopetrol = "Ecopetrol.csv", corpcolombia = "Corporacionfinancieradecolombia.csv",
    bolivar = "Grupobolivar.csv", mineros = "Mineros SA (MAS).csv", nutresa = "Nutresa (NCH).csv", terpel = "Organizacion Terpel SA (TPL).csv",
    promigas = "Promigas (PMG).csv", suramericana = "Suramericana (SIS).csv")
input_paths <- file.path("data", "raw", unname(company_files))
missing_files <- input_paths[!file.exists(input_paths)]
if (length(missing_files)) {
  stop(paste0("Some input files are missing. See DATA_ACCESS.md for access instructions and the CSV format.\nMissing files:\n",
              paste(paste0("  - ", missing_files), collapse = "\n")), call. = FALSE)
}
if (!all(file.exists(c("R/essential_functions.R", "R/reported_tables.R", "renv.lock")))) {
  stop("Run from the repository root containing run_analysis.R and DATA_ACCESS.md.", call. = FALSE)
}
if (as.character(getRversion()) != "4.4.3") {
  stop("This analysis requires R 4.4.3, as recorded in renv.lock.", call. = FALSE)
}
project_library <- file.path("renv", "library", paste0("R-", getRversion()), R.version$platform)
if (dir.exists(project_library)) .libPaths(c(normalizePath(project_library), .libPaths()))
required_packages <- c("dplyr", "tidyr", "purrr", "lubridate", "imputeTS", "fda", "zoo", "mgcv", "rugarch", "png", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), ". Restore renv.lock as described in README.md.", call. = FALSE)
}
lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
wrong_versions <- names(lock$Packages)[!vapply(lock$Packages, function(p) {
  requireNamespace(p$Package, quietly = TRUE) &&
    utils::packageVersion(p$Package) == package_version(p$Version)
}, logical(1))]
if (length(wrong_versions)) {
  stop("Missing or incompatible locked package versions: ", paste(wrong_versions, collapse = ", "),
       ". Restore renv.lock; see README.md.", call. = FALSE)
}
source("R/essential_functions.R", local = TRUE, encoding = "UTF-8")
source("R/reported_tables.R", local = TRUE, encoding = "UTF-8")
invisible(locate_ghostscript())
# Check input shape and identifiers without writing observations or echoing rows.
for (input_path in input_paths) {
  input_table <- utils::read.csv(input_path, sep = ";", stringsAsFactors = FALSE,
                                 check.names = FALSE, fileEncoding = "UTF-8-BOM")
  required_headers <- c("Fecha", "Último", "Apertura", "Máximo", "Mínimo", "Vol.")
  if (ncol(input_table) != 7L || !identical(names(input_table)[1:6], required_headers) ||
      !names(input_table)[7] %in% c("% var.", "%var.")) {
    stop("Unexpected CSV columns in ", input_path, ". See DATA_ACCESS.md.", call. = FALSE)
  }
  input_dates <- lubridate::dmy(input_table[[1]])
  if (anyNA(input_dates) || anyDuplicated(input_dates) ||
      any(input_dates < as.Date("2021-07-01") | input_dates > as.Date("2023-06-30"))) {
    stop("Invalid, duplicated or out-of-period dates in ", input_path, ". See DATA_ACCESS.md.", call. = FALSE)
  }
}
rm(input_table, input_dates)
suppressPackageStartupMessages(invisible(lapply(c("dplyr", "tidyr", "purrr", "lubridate", "fda", "zoo", "mgcv"), library, character.only = TRUE)))
analysis_started <- proc.time()[["elapsed"]]
message("Smoothing the curves, fitting the functional regression, and running diagnostics, leave-one-company-out refits and GARCH models.")

# Input parsing, imputation and univariate smoothing

convert_value <- function(s) {
    s <- as.character(s)
    s <- trimws(s)
    if (is.na(s) || s == "" || s == "-")
        return(NA_real_)
    s_up <- toupper(s)
    suffix <- substr(s_up, nchar(s_up), nchar(s_up))
    has_suffix <- suffix %in% c("K", "M", "B")
    num_str <- if (has_suffix) {
        substr(s_up, 1, nchar(s_up) - 1)
    }
    else {
        s_up
    }
    num_str <- gsub("\\.", "", num_str)
    num_str <- gsub(",", ".", num_str)
    value <- suppressWarnings(as.numeric(num_str))
    if (is.na(value))
        return(NA_real_)
    multiplier <- dplyr::case_when(!has_suffix ~ 1, suffix == "K" ~ 1000, suffix == "M" ~ 1e+06,
        suffix == "B" ~ 1e+09, TRUE ~ 1)
    value * multiplier
}

read_data <- function(path) {
    data <- read.csv(path, sep = ";", stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
    date <- gsub("/", "-", data[, 1])
    date <- lubridate::dmy(date)
    close <- as.numeric(gsub(",", ".", gsub("\\.", "", data[, 2])))
    volume <- lapply(data[, 6], convert_value)
    volume <- unlist(volume)
    data.frame(date = date, Close = close, Volume = volume) %>% dplyr::arrange(date)
}

create_base_dates <- function(company_series) {
    company_series %>% purrr::map(~dplyr::select(.x, date)) %>% dplyr::bind_rows() %>% dplyr::distinct(date) %>%
        dplyr::arrange(date)
}

assemble_wide_series <- function(company_series, variable, prefix) {
    df_wide <- create_base_dates(company_series)
    for (name in names(company_series)) {
        df_tmp <- company_series[[name]] %>% dplyr::select(date, dplyr::all_of(variable)) %>% dplyr::arrange(date) %>%
            dplyr::rename(`:=`(!!paste0(prefix, "_", name), dplyr::all_of(variable)))
        df_wide <- df_wide %>% dplyr::left_join(df_tmp, by = "date")
    }
    df_wide
}

compute_log_returns <- function(df_prices, impute_prices = TRUE) {
    df_base <- df_prices %>% dplyr::arrange(date)
    if (impute_prices) {
        df_base <- df_base %>% tidyr::fill(where(is.numeric), .direction = "down")
        df_base <- df_base %>% tidyr::fill(where(is.numeric), .direction = "up")
    }
    df_returns <- df_base %>% dplyr::mutate(dplyr::across(dplyr::starts_with("Close_"), ~c(NA_real_,
        diff(log(.x))), .names = "R_{sub('^Close_', '', .col)}")) %>% dplyr::select(date, dplyr::starts_with("R_")) %>%
        dplyr::slice(-1)
    rownames(df_returns) <- NULL
    df_returns
}

prepare_long_returns <- function(df_returns, company_labels) {
    df_returns %>% tidyr::pivot_longer(cols = -date, names_to = "Company", values_to = "Return") %>%
        dplyr::mutate(Company = gsub("^R_", "", Company), Company = dplyr::recode(Company, !!!company_labels),
            Company = factor(Company, levels = unname(company_labels)))
}

compute_metrics_fd <- function(Y_real, fd_smoothed, argvals) {
    Y_hat <- eval.fd(argvals, fd_smoothed)
    residuals <- Y_real - Y_hat
    RMS <- apply(residuals, 2, function(r) sqrt(mean(r^2, na.rm = TRUE)))
    sd_y <- apply(Y_real, 2, function(y) sd(y, na.rm = TRUE))
    range_y <- apply(Y_real, 2, function(y) diff(range(y, na.rm = TRUE)))
    RMS_rel <- ifelse(is.finite(sd_y) & sd_y > 0, RMS/sd_y, NA_real_)
    RMS_range <- ifelse(is.finite(range_y) & range_y > 0, RMS/range_y, NA_real_)
    data.frame(Series = colnames(Y_real), RMS = RMS, SD_Original = sd_y, RMS_Relative = RMS_rel,
        RMS_Range = RMS_range, row.names = NULL) %>% dplyr::arrange(dplyr::desc(RMS_Relative))
}

impute_volume_gam_log1p <- function(df, volume_col, k = NULL, method = "REML") {
    df <- df %>% dplyr::arrange(date)
    n <- nrow(df)
    if (is.null(k))
        k <- max(10, min(40, floor(n/12)))
    tnum <- as.numeric(df$date)
    dow <- lubridate::wday(df$date, label = TRUE, week_start = 1)
    y <- as.numeric(df[[volume_col]])
    z <- log1p(y)
    idx <- which(!is.na(z))
    if (length(idx) < 12) {
        z2 <- z
        z2 <- zoo::na.locf(z2, na.rm = FALSE)
        z2 <- zoo::na.locf(z2, fromLast = TRUE, na.rm = FALSE)
        z2[is.na(z2)] <- 0
        y_hat <- pmax(expm1(z2), 0)
        y_hat[!is.na(y)] <- y[!is.na(y)]
        return(y_hat)
    }
    fit <- mgcv::gam(z ~ s(tnum, k = k, bs = "cr") + dow, data = data.frame(z = z[idx], tnum = tnum[idx],
        dow = dow[idx]), method = method)
    z_pred <- predict(fit, newdata = data.frame(tnum = tnum, dow = dow))
    y_hat <- pmax(expm1(z_pred), 0)
    y_hat[!is.na(y)] <- y[!is.na(y)]
    y_hat
}

evaluate_gcv_fda <- function(Y, argvals, rangeval, nbasis, log10lambda, Lfdobj) {
    basis_obj <- fda::create.bspline.basis(rangeval = rangeval, nbasis = nbasis)
    lambda <- 10^log10lambda
    fdPar_obj <- fda::fdPar(basis_obj, Lfdobj, lambda)
    smooth_list <- fda::smooth.basis(argvals = argvals, y = Y, fdParobj = fdPar_obj)
    list(lambda = lambda, log10lambda = log10lambda, df_total = sum(smooth_list$df), gcv_total = sum(smooth_list$gcv),
        fd = smooth_list$fd, smooth = smooth_list)
}

compute_R2_t <- function(y_observed, y_estimated) {
    residuals <- y_observed - y_estimated
    SSE_t <- rowSums(residuals^2, na.rm = TRUE)
    SST_t <- rowSums(y_observed^2, na.rm = TRUE)
    R2_t <- 1 - SSE_t/SST_t
    R2_t[!is.finite(R2_t)] <- NA_real_
    R2_t
}

input_directory <- file.path("data", "raw")

company_labels <- c(bancolombia = "Bancolombia", bogota = "Banco de Bogotá", bvc = "Bolsa de Valores de Colombia",
    ceargos = "Cementos Argos", celsia = "Celsia S.A.", ecopetrol = "Ecopetrol", corpcolombia = "Corporación Financiera de Colombia",
    bolivar = "Grupo Bolívar", mineros = "Mineros", nutresa = "Nutresa", terpel = "Terpel", promigas = "Promigas",
    suramericana = "Suramericana")

company_series <- purrr::map(company_files, ~read_data(file.path(input_directory, .x)))

names(company_series) <- names(company_files)

dates_base <- create_base_dates(company_series)

df_prices <- assemble_wide_series(company_series = company_series, variable = "Close", prefix = "Close")

returns_before_price_imputation <- compute_log_returns(df_prices = df_prices, impute_prices = FALSE)

df_returns_interpolated <- returns_before_price_imputation

for (j in 2:ncol(df_returns_interpolated)) {
    df_returns_interpolated[[j]] <- imputeTS::na_interpolation(df_returns_interpolated[[j]],
        option = "spline")
}

rownames(df_returns_interpolated) <- NULL

df_return_interpolated_long <- prepare_long_returns(df_returns_interpolated, company_labels)

df_return_interpolated_long$date <- as.Date(df_return_interpolated_long$date)

dates_axis <- seq(from = as.Date("2021-09-01"), to = max(df_return_interpolated_long$date, na.rm = TRUE),
    by = "3 months")

axis_month_abbreviations <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov",
    "Dec")

labels_axis <- paste(axis_month_abbreviations[as.integer(format(dates_axis, "%m"))], format(dates_axis,
    "%Y"))

df_return_long <- prepare_long_returns(returns_before_price_imputation, company_labels)

df_return_long$date <- as.Date(df_return_long$date)

df_returns_final <- compute_log_returns(df_prices = df_prices, impute_prices = TRUE)

matrix_returns <- as.matrix(df_returns_final %>% dplyr::select(-date))

colnames(matrix_returns) <- gsub("^R_", "", colnames(matrix_returns))

n_days_return <- nrow(matrix_returns)

n_companies <- ncol(matrix_returns)

time_return <- seq_len(n_days_return)

range_return <- c(1, n_days_return)

nbasis_return_preliminary <- 140

base_return_preliminary <- fda::create.bspline.basis(rangeval = range_return, nbasis = nbasis_return_preliminary,
    norder = 4)

lambdas_log_return <- seq(-4, 6, by = 0.5)

gcv_return_preliminary <- rep(NA_real_, length(lambdas_log_return))

for (i in seq_along(lambdas_log_return)) {
    lambda_i <- 10^lambdas_log_return[i]
    fdPar_i <- fda::fdPar(base_return_preliminary, Lfdobj = 2, lambda = lambda_i)
    smooth_i <- fda::smooth.basis(time_return, matrix_returns, fdPar_i)
    gcv_return_preliminary[i] <- sum(smooth_i$gcv)
}

idx_lambda_return_preliminary <- which.min(gcv_return_preliminary)

lambda_return_preliminary <- 10^lambdas_log_return[idx_lambda_return_preliminary]

nbasis_grid_return <- c(80, 100, 125, 150)

lambdas_log_grid_return <- seq(-2, 6, by = 0.5)

results_return_opt <- data.frame()

for (nb in nbasis_grid_return) {
    base_temp <- fda::create.bspline.basis(range_return, nb, norder = 4)
    gcv_temp <- rep(NA_real_, length(lambdas_log_grid_return))
    for (i in seq_along(lambdas_log_grid_return)) {
        lambda_i <- 10^lambdas_log_grid_return[i]
        fdPar_i <- fda::fdPar(base_temp, Lfdobj = 2, lambda = lambda_i)
        smooth_i <- fda::smooth.basis(time_return, matrix_returns, fdPar_i)
        gcv_temp[i] <- sum(smooth_i$gcv)
    }
    idx_best <- which.min(gcv_temp)
    lambda_best_nb <- 10^lambdas_log_grid_return[idx_best]
    fdPar_best <- fda::fdPar(base_temp, Lfdobj = 2, lambda = lambda_best_nb)
    smooth_best <- fda::smooth.basis(time_return, matrix_returns, fdPar_best)
    y_hat_best <- fda::eval.fd(time_return, smooth_best$fd)
    residuals_best <- matrix_returns - y_hat_best
    rms_global <- sqrt(mean(residuals_best^2, na.rm = TRUE))
    sd_global <- sd(as.vector(matrix_returns), na.rm = TRUE)
    rms_rel_global <- rms_global/sd_global
    results_return_opt <- rbind(results_return_opt, data.frame(N_Bases = nb, Best_Lambda = lambda_best_nb,
        Log_Lambda = lambdas_log_grid_return[idx_best], Min_GCV = gcv_temp[idx_best], RMS = rms_global,
        RMS_Relative = rms_rel_global))
}

results_return_opt <- results_return_opt %>% dplyr::arrange(Min_GCV)

nbasis_sweep_return <- seq(40, 180, by = 10)

lambdas_log_sweep_return <- seq(-2, 6, by = 0.5)

results_K_return <- data.frame()

for (nb in nbasis_sweep_return) {
    base_temp <- fda::create.bspline.basis(range_return, nb, norder = 4)
    gcv_temp <- rep(NA_real_, length(lambdas_log_sweep_return))
    for (i in seq_along(lambdas_log_sweep_return)) {
        lambda_i <- 10^lambdas_log_sweep_return[i]
        fdPar_i <- fda::fdPar(base_temp, Lfdobj = 2, lambda = lambda_i)
        smooth_i <- fda::smooth.basis(time_return, matrix_returns, fdPar_i)
        gcv_temp[i] <- sum(smooth_i$gcv)
    }
    idx_best <- which.min(gcv_temp)
    results_K_return <- rbind(results_K_return, data.frame(K_Bases = nb, Min_GCV = gcv_temp[idx_best],
        Best_Lambda = 10^lambdas_log_sweep_return[idx_best]))
}

idx_min_K_return <- which.min(results_K_return$Min_GCV)

best_K_return_gcv <- results_K_return$K_Bases[idx_min_K_return]

min_gcv_K_return <- results_K_return$Min_GCV[idx_min_K_return]

row_K150 <- subset(results_K_return, K_Bases == 150)

gcv_K150 <- row_K150$Min_GCV


nbasis_return_final <- 150

base_return_final <- fda::create.bspline.basis(range_return, nbasis_return_final, norder = 4)

lambda_return_gcv <- 1e+05

lambda_return_mid <- 100

lambda_return_low <- 10

fdPar_return_gcv <- fda::fdPar(base_return_final, Lfdobj = 2, lambda = lambda_return_gcv)

fdPar_return_mid <- fda::fdPar(base_return_final, Lfdobj = 2, lambda = lambda_return_mid)

fdPar_return_low <- fda::fdPar(base_return_final, Lfdobj = 2, lambda = lambda_return_low)

fd_return_gcv <- fda::smooth.basis(time_return, matrix_returns, fdPar_return_gcv)$fd

fd_return_mid <- fda::smooth.basis(time_return, matrix_returns, fdPar_return_mid)$fd

fd_return_low <- fda::smooth.basis(time_return, matrix_returns, fdPar_return_low)$fd

company_idx <- 2

company_name <- colnames(matrix_returns)[company_idx]

lambda_return_final <- 100

fdPar_return_final <- fda::fdPar(base_return_final, Lfdobj = 2, lambda = lambda_return_final)

smooth_return_final <- fda::smooth.basis(time_return, matrix_returns, fdPar_return_final)

fd_returns <- smooth_return_final$fd

fd_returns$fdnames <- list(Time = "Days", Stock = colnames(matrix_returns), Return = "Log-return")

mean_returns <- fda::mean.fd(fd_returns)

table_summary_return <- data.frame(Parameter = c("Number of basis functions", "Lambda", "Total global GCV"), Value = c(as.character(nbasis_return_final),
    as.character(lambda_return_final), sprintf("%.8f", sum(smooth_return_final$gcv))))

table_performance_return <- compute_metrics_fd(Y_real = matrix_returns, fd_smoothed = fd_returns,
    argvals = time_return)

covariance_returns <- fda::var.fd(fd_returns)

range_covariance_return <- fd_returns$basis$rangeval

grid_covariance_return <- seq(range_covariance_return[1], range_covariance_return[2], length.out = 490)

matrix_covariance_return <- fda::eval.bifd(grid_covariance_return, grid_covariance_return, covariance_returns)

matrix_return_smoothed <- fda::eval.fd(time_return, fd_returns)

matrix_residuals_return <- matrix_returns - matrix_return_smoothed

indices_return_plot <- c(1, 2)

titles_return_left <- c("Bancolombia", "Banco de Bogotá")

titles_return_right <- c("Residuals | Bancolombia", "Residuals | Banco de Bogotá")

col_base <- "#CA717E"

col_points_row <- c("gray70", grDevices::adjustcolor(col_base, alpha.f = 0.35))

col_line_row <- c("black", col_base)

counter <- 0

for (i in indices_return_plot) {
    counter <- counter + 1
    y_obs <- matrix_returns[, i]
    y_fit <- matrix_return_smoothed[, i]
    y_res <- matrix_residuals_return[, i]
    rms_i <- sqrt(mean(y_res^2, na.rm = TRUE))
    m_obs <- max(abs(c(y_obs, y_fit)), na.rm = TRUE)
    if (!is.finite(m_obs) || m_obs == 0)
        m_obs <- 1
    ylim_obs <- c(-1, 1) * m_obs * 1.05
    m_res <- max(abs(y_res), na.rm = TRUE)
    if (!is.finite(m_res) || m_res == 0)
        m_res <- 1
    ylim_res <- c(-1, 1) * m_res * 1.05
}

lambda_return_fixed <- 100

return_basis_sizes <- c(150, 200, 250, 300, 365, 400)

results_bases_return <- data.frame()

for (nb in return_basis_sizes) {
    base_temp <- fda::create.bspline.basis(range_return, nb, norder = 4)
    fdPar_temp <- fda::fdPar(base_temp, Lfdobj = 2, lambda = lambda_return_fixed)
    smooth_temp <- fda::smooth.basis(time_return, matrix_returns, fdPar_temp)
    gcv_total <- sum(smooth_temp$gcv)
    y_hat_temp <- fda::eval.fd(time_return, smooth_temp$fd)
    residuals_temp <- matrix_returns - y_hat_temp
    rms_abs <- sqrt(colMeans(residuals_temp^2, na.rm = TRUE))
    sd_orig <- apply(matrix_returns, 2, sd, na.rm = TRUE)
    rms_rel_mean <- mean(rms_abs/sd_orig, na.rm = TRUE)
    results_bases_return <- rbind(results_bases_return, data.frame(N_Bases = nb, GCV_Total = gcv_total,
        RMS_Rel_Mean = rms_rel_mean))
}

df_volume <- assemble_wide_series(company_series = company_series, variable = "Volume", prefix = "Volume")

date_remove <- as.Date("2021-07-01")

df_volume <- df_volume %>% dplyr::filter(date != date_remove)

volume_cols <- grep("^Volume_", names(df_volume), value = TRUE)

df_volume_imputed <- df_volume

for (col in volume_cols) {
    df_volume_imputed[[col]] <- impute_volume_gam_log1p(df = df_volume_imputed, volume_col = col)
}

df_volume_long <- df_volume_imputed %>% dplyr::mutate(date = as.Date(date)) %>% dplyr::select(date, dplyr::all_of(volume_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(volume_cols), names_to = "Company", values_to = "Volume") %>%
    dplyr::mutate(Company = gsub("^Volume_", "", Company), Company = dplyr::recode(Company, !!!company_labels))

df_volume_log <- df_volume_imputed %>% dplyr::mutate(dplyr::across(dplyr::all_of(volume_cols), ~log1p(.x)))

volume_cols <- grep("^Volume_", names(df_volume_log), value = TRUE)

matrix_volume <- as.matrix(df_volume_log[, volume_cols])

colnames(matrix_volume) <- gsub("^Volume_", "", colnames(matrix_volume))

n_days_volume <- nrow(matrix_volume)

time_volume <- seq_len(n_days_volume)

range_volume <- c(1, n_days_volume)

Lfd_volume <- fda::int2Lfd(2)

nbasis_grid_volume <- seq(5, 25, by = 1)

log10lambda_grid_volume <- seq(1, 10, by = 1)

results_volume_list <- list()

for (nb in nbasis_grid_volume) {
    gcv_vec <- rep(NA_real_, length(log10lambda_grid_volume))
    df_vec <- rep(NA_real_, length(log10lambda_grid_volume))
    for (i in seq_along(log10lambda_grid_volume)) {
        lg <- log10lambda_grid_volume[i]
        tmp <- evaluate_gcv_fda(Y = matrix_volume, argvals = time_volume, rangeval = range_volume, nbasis = nb,
            log10lambda = lg, Lfdobj = Lfd_volume)
        gcv_vec[i] <- tmp$gcv_total
        df_vec[i] <- tmp$df_total
    }
    idx_best <- which.min(gcv_vec)
    best_log10lambda_nb <- log10lambda_grid_volume[idx_best]
    best_lambda_nb <- 10^best_log10lambda_nb
    results_volume_list[[as.character(nb)]] <- list(nbasis = nb, log10lambda_grid = log10lambda_grid_volume,
        gcv_vec = gcv_vec, df_vec = df_vec, best_log10lambda = best_log10lambda_nb, best_lambda = best_lambda_nb)
}

summary_volume <- do.call(rbind, lapply(results_volume_list, function(x) {
    idx_best <- which.min(x$gcv_vec)
    data.frame(nbasis = x$nbasis, best_log10lambda = x$log10lambda_grid[idx_best], best_lambda = 10^(x$log10lambda_grid[idx_best]),
        best_gcv = x$gcv_vec[idx_best], best_df_total = x$df_vec[idx_best])
}))

summary_volume <- summary_volume %>% dplyr::arrange(best_gcv)

best_nb_volume <- summary_volume$nbasis[1]

best_log10lambda_volume <- summary_volume$best_log10lambda[1]

best_lambda_volume <- summary_volume$best_lambda[1]

summary_volume_K <- summary_volume %>% dplyr::arrange(nbasis)

base_volume_final <- fda::create.bspline.basis(rangeval = range_volume, nbasis = best_nb_volume)

fdPar_volume_final <- fda::fdPar(base_volume_final, Lfd_volume, best_lambda_volume)

smooth_volume_final <- fda::smooth.basis(argvals = time_volume, y = matrix_volume, fdParobj = fdPar_volume_final)

fd_volume <- smooth_volume_final$fd

gcv_volume_final <- sum(smooth_volume_final$gcv)

fd_volume$fdnames <- list(Time = "Days", Stock = colnames(matrix_volume), Volume = "log(1 + trading volume)")

mean_volume <- fda::mean.fd(fd_volume)

covariance_volume <- fda::var.fd(fd_volume)

range_covariance_volume <- fd_volume$basis$rangeval

grid_covariance_volume <- seq(range_covariance_volume[1], range_covariance_volume[2], length.out = 490)

matrix_covariance_volume <- fda::eval.bifd(grid_covariance_volume, grid_covariance_volume, covariance_volume)

companies_volume_plot <- 1:min(2, ncol(matrix_volume))

titles_volume_left <- c("Bancolombia", "Banco de Bogotá")

titles_volume_right <- c("Residuals | Bancolombia", "Residuals | Banco de Bogotá")

matrix_volume_smoothed <- fda::eval.fd(time_volume, fd_volume)

cols <- seq_along(companies_volume_plot)

pt_col <- grDevices::adjustcolor(cols, alpha.f = 0.35)

ln_col <- grDevices::adjustcolor(cols, alpha.f = 0.95)

for (k in seq_along(companies_volume_plot)) {
    j <- companies_volume_plot[k]
    residual_volume_j <- matrix_volume[, j] - matrix_volume_smoothed[, j]
}

table_rms_volume <- compute_metrics_fd(Y_real = matrix_volume, fd_smoothed = fd_volume, argvals = time_volume)

compare_lambdas_volume <- function(lambdas) {
    out <- lapply(lambdas, function(lambda_i) {
        basis_i <- fda::create.bspline.basis(rangeval = range_volume, nbasis = best_nb_volume)
        fdPar_i <- fda::fdPar(basis_i, Lfd_volume, lambda_i)
        smooth_i <- fda::smooth.basis(time_volume, matrix_volume, fdPar_i)
        y_hat_i <- fda::eval.fd(time_volume, smooth_i$fd)
        resid_i <- matrix_volume - y_hat_i
        RMS_i <- apply(resid_i, 2, function(r) sqrt(mean(r^2, na.rm = TRUE)))
        sd_i <- apply(matrix_volume, 2, function(y) sd(y, na.rm = TRUE))
        RMS_rel_i <- ifelse(sd_i > 0, RMS_i/sd_i, NA_real_)
        data.frame(lambda = lambda_i, RMS_median = median(RMS_i, na.rm = TRUE), RMSrel_median = median(RMS_rel_i,
            na.rm = TRUE), RMSrel_max = max(RMS_rel_i, na.rm = TRUE), df = smooth_i$df)
    })
    do.call(rbind, out)
}

lambdas_volume_compare <- best_lambda_volume * c(0.1, 1, 10)


# Functional regression and penalty sensitivity
xfd_volume_centered <- fda::center.fd(fd_volume)

yfd_return_centered <- fda::center.fd(fd_returns)

x_basis <- xfd_volume_centered$basis

y_basis <- yfd_return_centered$basis

beta_coef_base <- fda::bifd(coef = matrix(0, x_basis$nbasis, y_basis$nbasis), sbasisobj = x_basis, tbasisobj = y_basis)

# We chose B after inspecting the three surfaces on the same 150-by-150 grid.
# The two penalties vary together across scenarios.
scenarios_sensitivity <- data.frame(
    scenario = c("A. Less smoothing", "B. Selected model", "C. More smoothing"),
    lambda_s = c(1e2, 1e3, 1e4),
    lambda_t = c(1e7, 1e8, 1e9)
)
selected_scenario <- match("B. Selected model", scenarios_sensitivity$scenario)

grid_s_sensitivity <- seq(from = min(x_basis$rangeval), to = max(x_basis$rangeval), length.out = 150)

grid_t_sensitivity <- seq(from = min(y_basis$rangeval), to = max(y_basis$rangeval), length.out = 150)

compute_scenario_sensitivity <- function(scenario_name, lambda_s, lambda_t) {
    beta_par_temp <- fda::bifdPar(beta_coef_base, fda::int2Lfd(2), fda::int2Lfd(2), lambda_s, lambda_t)
    beta_list_temp <- list(fda::fdPar(yfd_return_centered$basis, Lfdobj = 2, lambda = 1e+05), beta_par_temp)
    model_temp <- fda::linmod(xfd_volume_centered, yfd_return_centered, beta_list_temp)
    beta_mat_temp <- fda::eval.bifd(grid_s_sensitivity, grid_t_sensitivity, model_temp$beta1estbifd)
    y_estimated_temp <- fda::eval.fd(grid_t_sensitivity, model_temp$yhatfdobj)
    y_observed_temp <- fda::eval.fd(grid_t_sensitivity, yfd_return_centered)
    resid_temp <- y_observed_temp - y_estimated_temp
    var_y_temp <- apply(y_observed_temp, 1, var, na.rm = TRUE)
    var_res_temp <- apply(resid_temp, 1, var, na.rm = TRUE)
    r2_t_temp <- 1 - var_res_temp/var_y_temp
    r2_mean_temp <- mean(r2_t_temp[is.finite(r2_t_temp)], na.rm = TRUE)
    beta_mat_plot <- beta_mat_temp
    df_temp <- expand.grid(s = grid_s_sensitivity, t = grid_t_sensitivity) %>% mutate(beta = as.vector(beta_mat_plot),
        scenario = scenario_name, lambda_s = lambda_s, lambda_t = lambda_t, r2_mean = r2_mean_temp)
    return(df_temp)
}

df_sensitivity_plot <- do.call(rbind, lapply(seq_len(nrow(scenarios_sensitivity)), function(i) {
    compute_scenario_sensitivity(scenario_name = scenarios_sensitivity$scenario[i], lambda_s = scenarios_sensitivity$lambda_s[i],
        lambda_t = scenarios_sensitivity$lambda_t[i])
}))

format_power_10 <- function(x) {
    exponent <- round(log10(x))
    paste0("10^", exponent)
}

df_sensitivity_plot <- df_sensitivity_plot %>% mutate(lambda_s_lab = format_power_10(lambda_s),
    lambda_t_lab = format_power_10(lambda_t), scenario_label = paste0(scenario, "\n", "lambda_s = ",
        lambda_s_lab, "\n", "lambda_t = ", lambda_t_lab))

beta_lim <- range(df_sensitivity_plot$beta, na.rm = TRUE)

beta_breaks <- pretty(beta_lim, n = 7)

beta_labels <- function(x) {
    sprintf("%.6f", x)
}

jet_cols <- jet_publication

lambda_s_final <- scenarios_sensitivity$lambda_s[selected_scenario]
lambda_t_final <- scenarios_sensitivity$lambda_t[selected_scenario]

beta_par_final <- fda::bifdPar(bifdobj = beta_coef_base, Lfdobjs = fda::int2Lfd(2), Lfdobjt = fda::int2Lfd(2),
    lambdas = lambda_s_final, lambdat = lambda_t_final)

beta_list_final <- list(fda::fdPar(yfd_return_centered$basis, Lfdobj = 2, lambda = 1e+05), beta_par_final)

model_final <- fda::linmod(xfdobj = xfd_volume_centered, yfdobj = yfd_return_centered, beta_list_final)

grid_beta_s_final <- seq(x_basis$rangeval[1], x_basis$rangeval[2], length.out = 490)

grid_beta_t_final <- seq(y_basis$rangeval[1], y_basis$rangeval[2], length.out = 490)

beta_final_mat <- fda::eval.bifd(grid_beta_s_final, grid_beta_t_final, model_final$beta1estbifd)

limits_beta_final <- range(beta_final_mat, na.rm = TRUE)

dates_s <- as.Date(df_volume_log$date)

dates_t <- as.Date(df_returns_final$date)

stopifnot(length(dates_s) == nrow(beta_final_mat), length(dates_t) == ncol(beta_final_mat))

df_beta_final <- expand.grid(s = seq_along(dates_s), t = seq_along(dates_t))

df_beta_final$beta <- as.vector(beta_final_mat)

breaks_s <- unique(round(seq(1, length(dates_s), length.out = 9)))

breaks_t <- unique(round(seq(1, length(dates_t), length.out = 9)))

labels_s <- format(dates_s[breaks_s], "%b\n%Y")

labels_t <- format(dates_t[breaks_t], "%b\n%Y")

jet_cols <- jet_publication

beta_lim <- limits_beta_final

beta_breaks <- pretty(beta_lim, n = 7)

beta_breaks <- beta_breaks[beta_breaks >= beta_lim[1] & beta_breaks <= beta_lim[2]]

beta_labels <- function(x) {
    sprintf("%.6f", x)
}

summary_beta_s <- data.frame(s = seq_len(nrow(beta_final_mat)), date = dates_s, Beta_mean = rowMeans(beta_final_mat,
    na.rm = TRUE), Beta_RMS = sqrt(rowMeans(beta_final_mat^2, na.rm = TRUE)), Proportion_beta_positive = rowMeans(beta_final_mat >
    0, na.rm = TRUE))

if (nrow(beta_final_mat) == ncol(beta_final_mat)) {
    summary_beta_s$Beta_contemporaneous <- diag(beta_final_mat)
}

events_market <- data.frame(date = as.Date(c("2022-02-01", "2022-03-30", "2022-06-21", "2022-10-20",
    "2022-11-21", "2023-02-28", "2023-03-15", "2023-04-27", "2023-05-24")), Event = c("Resumption of Grupo SURA trading in the takeover-bid context",
    "Resumption of Banco de Bogotá trading after the BAC Holding spin-off", "First trading session after the presidential election",
    "Resumption of Nutresa trading after approval of the IHC takeover bid", "Nutresa decline after the unsuccessful IHC takeover bid",
    "Market reaction to cabinet changes and health reform", "Unusual BVC movement during regional exchange integration",
    "Change of finance minister and market reaction", "GEA-Gilinski memorandum concerning Nutresa and Grupo SURA"))

nearest_session_index <- function(date, session_dates) {
    which.min(abs(as.numeric(session_dates - date)))
}

events_market$Index_s <- sapply(events_market$date, nearest_session_index, session_dates = dates_s)

events_market$date_in_data <- dates_s[events_market$Index_s]

table_events_beta <- events_market %>% dplyr::left_join(summary_beta_s %>% dplyr::select(s, Beta_mean,
    Beta_RMS, Proportion_beta_positive, dplyr::any_of("Beta_contemporaneous")), by = c(Index_s = "s")) %>% dplyr::select(Index_s,
    date, date_in_data, Event, Beta_mean, Beta_RMS, Proportion_beta_positive, dplyr::any_of("Beta_contemporaneous"))

table_events_beta_formatted <- table_events_beta %>% dplyr::mutate(Beta_mean = round(Beta_mean,
    7), Beta_RMS = round(Beta_RMS, 7), Proportion_beta_positive = round(Proportion_beta_positive, 3), dplyr::across(dplyr::any_of("Beta_contemporaneous"),
    ~round(.x, 7)))

if (!identical(dates_s, dates_t)) {
    stop("Return and volume sessions do not align. Check the source files against DATA_ACCESS.md.", call. = FALSE)
}

df_beta_map <- expand.grid(s = seq_len(nrow(beta_final_mat)), t = seq_len(ncol(beta_final_mat)))

df_beta_map$beta <- as.vector(beta_final_mat)

events_plot <- table_events_beta %>% dplyr::mutate(Code = paste0("E", dplyr::row_number()))

breaks_dates <- unique(round(seq(1, length(dates_s), length.out = 9)))

month_abbreviations <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

labels_dates <- paste0(month_abbreviations[as.integer(format(dates_s[breaks_dates], "%m"))], "\n", format(dates_s[breaks_dates],
    "%Y"))

table_events_figure <- events_plot %>% dplyr::select(Code, date, Event, Beta_mean, Beta_RMS, Proportion_beta_positive,
    Beta_contemporaneous)

table_key_events <- events_plot %>% dplyr::select(Code, date, Event, Beta_mean, Beta_RMS, Proportion_beta_positive,
    Beta_contemporaneous) %>% dplyr::mutate(Beta_mean = sprintf("%.6f", Beta_mean), Beta_RMS = sprintf("%.6f",
    Beta_RMS), Proportion_beta_positive = sprintf("%.1f%%", 100 * Proportion_beta_positive), Beta_contemporaneous = sprintf("%.6f",
    Beta_contemporaneous))

breaks_map <- unique(round(seq(1, length(dates_s), length.out = 7)))

labels_map <- paste0(month_abbreviations[as.integer(format(dates_s[breaks_map], "%m"))], "\n", format(dates_s[breaks_map],
    "%Y"))

range_mean <- diff(range(summary_beta_s$Beta_mean, na.rm = TRUE))

events_plot$Label_mean <- events_plot$Beta_mean + rep(c(0.08, -0.11), length.out = nrow(events_plot)) *
    range_mean

range_rms <- diff(range(summary_beta_s$Beta_RMS, na.rm = TRUE))

events_plot$Label_RMS <- events_plot$Beta_RMS + rep(c(0.08, -0.1), length.out = nrow(events_plot)) *
    range_rms

grid_r2_final <- seq(y_basis$rangeval[1], y_basis$rangeval[2], length.out = 490)

y_estimated_final <- fda::eval.fd(grid_r2_final, model_final$yhatfdobj)

y_observed_final <- fda::eval.fd(grid_r2_final, yfd_return_centered)

residuals_final <- y_observed_final - y_estimated_final

R2_t_final <- compute_R2_t(y_observed = y_observed_final, y_estimated = y_estimated_final)

ylim_R2_final <- c(min(R2_t_final, 0, na.rm = TRUE), 1)

stopifnot(length(dates_t) == length(R2_t_final))

mean_temporal_R2_final <- mean(R2_t_final, na.rm = TRUE)

df_R2_final <- data.frame(date = dates_t, R2 = R2_t_final)

df_mean_R2_final <- data.frame(Mean_temporal = mean_temporal_R2_final)

breaks_R2 <- unique(round(seq(1, length(dates_t), length.out = 7)))

labels_R2 <- paste0(month_abbreviations[as.integer(format(dates_t[breaks_R2], "%m"))], "\n", format(dates_t[breaks_R2],
    "%Y"))

label_mean_R2 <- paste0("Temporal mean of R²(t) = ", sprintf("%.3f", mean_temporal_R2_final))


# Model diagnostics
fd_residual_centered <- yfd_return_centered - model_final$yhatfdobj

SSE_series_fd <- diag(fda::inprod(fd_residual_centered, fd_residual_centered))

MSV_series_fd <- diag(fda::inprod(yfd_return_centered, yfd_return_centered))

R2_series_fd <- 1 - SSE_series_fd/MSV_series_fd

names(R2_series_fd) <- colnames(y_observed_final)

table_R2_verification <- data.frame(Company = names(R2_series_fd), SSE_i = SSE_series_fd, MSV_i = MSV_series_fd,
    R2_i = R2_series_fd)

RMS_model <- colMeans((y_observed_final - y_estimated_final)^2, na.rm = TRUE)

MSV_model <- colMeans(y_observed_final^2, na.rm = TRUE)

fd_error_model <- yfd_return_centered - model_final$yhatfdobj

mean_error_model <- fda::mean.fd(fd_error_model)

sd_error_model <- fda::std.fd(fd_error_model)

grid_diag <- seq(yfd_return_centered$basis$rangeval[1], yfd_return_centered$basis$rangeval[2], length.out = 490)

matrix_observed_diag <- fda::eval.fd(grid_diag, yfd_return_centered)

matrix_predicted_diag <- fda::eval.fd(grid_diag, model_final$yhatfdobj)

matrix_residuals_diag <- matrix_observed_diag - matrix_predicted_diag

fd_error_model <- yfd_return_centered - model_final$yhatfdobj

stopifnot(all(dim(matrix_observed_diag) == c(490, 13)), all(dim(matrix_predicted_diag) == c(490, 13)),
    all(dim(matrix_residuals_diag) == c(490, 13)), length(grid_diag) == 490)

res_vec <- as.vector(matrix_residuals_diag)

fit_vec <- as.vector(matrix_predicted_diag)

df_res_fit <- data.frame(Fitted = fit_vec, Residual = res_vec) %>% dplyr::filter(is.finite(Fitted),
    is.finite(Residual))

mean_res_t <- rowMeans(matrix_residuals_diag, na.rm = TRUE)

sd_res_t <- apply(matrix_residuals_diag, 1, sd, na.rm = TRUE)

n_curves <- ncol(matrix_residuals_diag)

df_mean_res <- data.frame(Time = grid_diag, MeanResidual = mean_res_t, UpperBand = mean_res_t +
    1.96 * sd_res_t/sqrt(n_curves), LowerBand = mean_res_t - 1.96 * sd_res_t/sqrt(n_curves))

table_var_residual <- data.frame(Company = colnames(matrix_residuals_diag), RMS_residual = sqrt(colMeans(matrix_residuals_diag^2,
    na.rm = TRUE)), Variance_residual = apply(matrix_residuals_diag, 2, var, na.rm = TRUE)) %>% dplyr::arrange(dplyr::desc(RMS_residual))

covariance_residual <- fda::var.fd(fd_error_model)

grid_covariance_res <- seq(fd_error_model$basis$rangeval[1], fd_error_model$basis$rangeval[2], length.out = 150)

matrix_covariance_res <- fda::eval.bifd(grid_covariance_res, grid_covariance_res, covariance_residual)

cor_residual_companies <- cor(matrix_residuals_diag, use = "pairwise.complete.obs")

res_finite <- as.vector(matrix_residuals_diag)

res_finite <- res_finite[is.finite(res_finite)]

mean_res <- mean(res_finite)

sd_res <- sd(res_finite)

skewness_res <- mean(((res_finite - mean_res)/sd_res)^3)

kurtosis_res <- mean(((res_finite - mean_res)/sd_res)^4)

excess_kurtosis_res <- kurtosis_res - 3

table_normality_res <- data.frame(Mean = mean_res, SD = sd_res, Skewness = skewness_res,
    Kurtosis = kurtosis_res, Excess_Kurtosis = excess_kurtosis_res)

set.seed(123)

sample_shapiro <- sample(res_finite, size = min(5000, length(res_finite)))

shapiro_result <- shapiro.test(sample_shapiro)

alpha_hat <- as.vector(fda::eval.fd(time_return, model_final$beta0estfd))

rms_alpha <- sqrt(mean(alpha_hat^2, na.rm = TRUE))

Y_c_mat <- fda::eval.fd(time_return, yfd_return_centered)

rms_Y_c <- sqrt(mean(Y_c_mat^2, na.rm = TRUE))

Yhat_mat <- fda::eval.fd(time_return, model_final$yhatfdobj)

rms_Yhat <- sqrt(mean(Yhat_mat^2, na.rm = TRUE))


# Extreme-observation figure inputs
df_return_review <- df_returns_final

df_volume_review <- df_volume

df_return_review$date <- as.Date(df_return_review$date)

df_volume_review$date <- as.Date(df_volume_review$date)

return_cols_review <- grep("^R_", names(df_return_review), value = TRUE)

volume_cols_review <- grep("^Volume_", names(df_volume_review), value = TRUE)

if (length(return_cols_review) != 13) {
    stop("Expected 13 return columns; found ", length(return_cols_review), ".")
}

if (length(volume_cols_review) != 13) {
    stop("Expected 13 volume columns; found ", length(volume_cols_review), ".")
}

labels_review <- c(bancolombia = "Bancolombia", bogota = "Banco de Bogotá", bvc = "Bolsa de Valores de Colombia",
    ceargos = "Cementos Argos", celsia = "Celsia S.A.", ecopetrol = "Ecopetrol", corpcolombia = "CorpColombia",
    bolivar = "Grupo Bolívar", mineros = "Mineros", nutresa = "Nutresa", terpel = "Terpel", promigas = "Promigas",
    suramericana = "Suramericana")

keys_return <- sub("^R_", "", return_cols_review)

keys_volume <- sub("^Volume_", "", volume_cols_review)

unlabelled_keys <- setdiff(union(keys_return, keys_volume), names(labels_review))

if (length(unlabelled_keys) > 0) {
    stop("Missing labels for: ", paste(unlabelled_keys, collapse = ", "))
}

returns_long <- df_return_review %>% dplyr::select(date, dplyr::all_of(return_cols_review)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(return_cols_review), names_to = "Company_key", values_to = "Return") %>%
    dplyr::mutate(Company_key = sub("^R_", "", Company_key), Company = unname(labels_review[Company_key])) %>%
    dplyr::filter(!is.na(date), is.finite(Return), !is.na(Company)) %>% dplyr::select(date, Company,
    Return)

volume_observed_long <- df_volume_review %>% dplyr::select(date, dplyr::all_of(volume_cols_review)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(volume_cols_review), names_to = "Company_key", values_to = "Volume") %>%
    dplyr::mutate(Company_key = sub("^Volume_", "", Company_key), Company = unname(labels_review[Company_key])) %>%
    dplyr::filter(!is.na(date), !is.na(Company)) %>% dplyr::select(date, Company, Volume)

table_returns_extremes <- returns_long %>% dplyr::mutate(Magnitude_absolute = abs(Return)) %>%
    dplyr::group_by(Company) %>% dplyr::arrange(dplyr::desc(Magnitude_absolute), date, .by_group = TRUE) %>%
    dplyr::slice_head(n = 3) %>% dplyr::mutate(Rank = dplyr::row_number(), Direction = dplyr::if_else(Return >=
    0, "Increase", "Decrease")) %>% dplyr::ungroup() %>% dplyr::arrange(Company, Rank) %>% dplyr::transmute(date,
    Company, Rank, Return, Magnitude_absolute, Direction, Event_possible = NA_character_, Source_consulted = NA_character_,
    Link_source = NA_character_, Level_confidence = NA_character_, Decision = "Retain", Justification = NA_character_)

table_volume_high <- volume_observed_long %>% dplyr::filter(is.finite(Volume)) %>% dplyr::group_by(Company) %>%
    dplyr::arrange(dplyr::desc(Volume), date, .by_group = TRUE) %>% dplyr::slice_head(n = 3) %>% dplyr::mutate(Rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>% dplyr::arrange(Company, Rank) %>% dplyr::transmute(date, Company, Rank,
    Volume, Volume_log = log1p(Volume), Type_review = "High volume", Event_possible = NA_character_,
    Source_consulted = NA_character_, Link_source = NA_character_, Level_confidence = NA_character_,
    Decision = "Retain", Justification = NA_character_)

table_volume_low <- volume_observed_long %>% dplyr::filter(is.finite(Volume), Volume > 0) %>% dplyr::group_by(Company) %>%
    dplyr::arrange(Volume, date, .by_group = TRUE) %>% dplyr::slice_head(n = 3) %>% dplyr::mutate(Rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>% dplyr::arrange(Company, Rank) %>% dplyr::transmute(date, Company, Rank,
    Volume, Volume_log = log1p(Volume), Type_review = "Low positive volume", Event_possible = NA_character_,
    Source_consulted = NA_character_, Link_source = NA_character_, Level_confidence = NA_character_,
    Decision = "Retain", Justification = NA_character_)

table_summary_zeros <- volume_observed_long %>% dplyr::group_by(Company) %>% dplyr::summarise(Count_zeros = sum(Volume ==
    0, na.rm = TRUE), Count_observed = sum(!is.na(Volume)), .groups = "drop") %>% dplyr::mutate(Percentage_zeros = dplyr::if_else(Count_observed >
    0, 100 * Count_zeros/Count_observed, NA_real_)) %>% dplyr::arrange(Company) %>% dplyr::select(Company,
    Count_zeros, Percentage_zeros)

table_dates_zeros <- volume_observed_long %>% dplyr::filter(!is.na(Volume), Volume == 0) %>% dplyr::arrange(date,
    Company) %>% dplyr::select(date, Company, Volume)

summarise_dates_review <- function(table) {
    table %>% dplyr::group_by(date) %>% dplyr::summarise(Count_companies = dplyr::n_distinct(Company),
        Companies = paste(sort(unique(Company)), collapse = ", "), .groups = "drop") %>% dplyr::arrange(dplyr::desc(Count_companies),
        date)
}

table_dates_returns <- summarise_dates_review(table_returns_extremes)

table_dates_volume_high <- summarise_dates_review(table_volume_high)

table_dates_volume_low <- summarise_dates_review(table_volume_low)

count_per_company <- function(table) {
    table %>% dplyr::count(Company, name = "n") %>% dplyr::arrange(Company)
}

count_return <- count_per_company(table_returns_extremes)

count_v_high <- count_per_company(table_volume_high)

count_v_low <- count_per_company(table_volume_low)

stopifnot(nrow(table_returns_extremes) == 39, nrow(table_volume_high) == 39, nrow(table_volume_low) ==
    39, all(count_return$n == 3), all(count_v_high$n == 3), all(count_v_low$n == 3), all.equal(table_returns_extremes$Magnitude_absolute,
    abs(table_returns_extremes$Return), tolerance = 1e-12) == TRUE, all.equal(table_volume_high$Volume_log,
    log1p(table_volume_high$Volume), tolerance = 1e-12) == TRUE, all.equal(table_volume_low$Volume_log,
    log1p(table_volume_low$Volume), tolerance = 1e-12) == TRUE)


# Leave-one-company-out refits
refit_loco <- function(model_env) {
    e <- new.env(parent = model_env)
    evalq({
        lambda_s <- lambda_s_final
        lambda_t <- lambda_t_final
        lambda_alpha <- 1e+05
        if (!isTRUE(all.equal(lambda_s, 1000)) || !isTRUE(all.equal(lambda_t, 1e+08))) {
            stop("Final model penalties changed; review before running LOCO.")
        }
        companies <- as.character(fd_volume$fdnames[[2]])
        companies_y <- as.character(fd_returns$fdnames[[2]])
        if (!identical(companies, companies_y) || length(companies) != 13L || anyDuplicated(companies)) {
            stop("The 13 companies or their X/Y order do not match.")
        }
        nx <- ncol(as.matrix(xfd_volume_centered$coefs))
        ny <- ncol(as.matrix(yfd_return_centered$coefs))
        if (nx != 13L || ny != 13L)
            stop("Expected 13 curves in each fd object.")
        grid_s <- as.numeric(grid_beta_s_final)
        grid_t <- as.numeric(grid_beta_t_final)
        if (length(grid_s) != 490L || length(grid_t) != 490L || any(diff(grid_s) <= 0) || any(diff(grid_t) <=
            0)) {
            stop("Expected the 490 increasing grid positions of the final model.")
        }
        beta_full <- fda::eval.bifd(grid_s, grid_t, model_final$beta1estbifd)
        if (!identical(dim(beta_full), c(490L, 490L)) || any(!is.finite(beta_full)))
            stop("Invalid full coefficient surface.")
        if (exists("beta_final_mat", inherits = TRUE)) {
            if (!isTRUE(all.equal(beta_full, beta_final_mat, tolerance = 1e-10, check.attributes = FALSE))) {
                stop("The full coefficient surface does not match the final model.")
            }
        }
        weights_trapezoidal <- function(points) {
            steps <- diff(points)
            n <- length(points)
            if (n == 2L)
                return(rep(steps[1]/2, 2L))
            c(steps[1]/2, (steps[-length(steps)] + steps[-1L])/2, steps[length(steps)]/2)
        }
        weights <- outer(weights_trapezoidal(grid_s), weights_trapezoidal(grid_t))
        norm_l2 <- function(m) sqrt(sum(weights * m^2))
        norm_base <- norm_l2(beta_full)
        if (!is.finite(norm_base) || norm_base <= 0)
            stop("Invalid full-surface L2 norm.")
        max_beta <- max(abs(beta_full))
        thresholds <- c(0.01, 0.05, 0.1)
        surfaces_loco <- array(NA_real_, dim = c(490L, 490L, length(companies)), dimnames = list(NULL,
            NULL, companies))
        metrics <- vector("list", length(companies))
        sensitivity <- vector("list", length(companies))
        xcoef <- as.matrix(xfd_volume_centered$coefs)
        ycoef <- as.matrix(yfd_return_centered$coefs)
        for (i in seq_along(companies)) {
            retain <- setdiff(seq_along(companies), i)
            x_sub <- fda::fd(coef = xcoef[, retain, drop = FALSE], basisobj = xfd_volume_centered$basis)
            y_sub <- fda::fd(coef = ycoef[, retain, drop = FALSE], basisobj = yfd_return_centered$basis)
            x_sub <- fda::center.fd(x_sub)
            y_sub <- fda::center.fd(y_sub)
            beta_par <- fda::bifdPar(bifdobj = beta_coef_base, Lfdobjs = fda::int2Lfd(2), Lfdobjt = fda::int2Lfd(2),
                lambdas = lambda_s, lambdat = lambda_t)
            beta_list <- list(fda::fdPar(yfd_return_centered$basis, Lfdobj = 2, lambda = lambda_alpha),
                beta_par)
            fit <- fda::linmod(xfdobj = x_sub, yfdobj = y_sub, betaList = beta_list)
            beta_i <- fda::eval.bifd(grid_s, grid_t, fit$beta1estbifd)
            if (!identical(dim(beta_i), dim(beta_full)) || any(!is.finite(beta_i))) {
                stop("Invalid LOCO surface when omitting ", companies[i])
            }
            surfaces_loco[, , i] <- beta_i
            distance <- norm_l2(beta_i - beta_full)/norm_base
            agreement <- sign(beta_i) == sign(beta_full)
            agreement_global <- sum(weights * agreement)/sum(weights)
            idx_5 <- abs(beta_full) >= 0.05 * max_beta
            agreement_5 <- sum(weights[idx_5] * agreement[idx_5])/sum(weights[idx_5])
            metrics[[i]] <- data.frame(Company_omitted = companies[i], Distance_L2_relative = distance,
                Correlation_surfaces = stats::cor(as.vector(beta_full), as.vector(beta_i)), Agreement_sign_global = agreement_global,
                Agreement_sign_region_active_5pct = agreement_5, Max_difference_absolute = max(abs(beta_i -
                  beta_full)), stringsAsFactors = FALSE)
            sensitivity[[i]] <- do.call(rbind, lapply(thresholds, function(u) {
                active <- abs(beta_full) >= u * max_beta
                data.frame(Company_omitted = companies[i], Threshold_relative = u, Fraction_area_active = sum(weights[active])/sum(weights),
                  Agreement_sign_area_active = sum(weights[active] * agreement[active])/sum(weights[active]),
                  stringsAsFactors = FALSE)
            }))
        }
        metrics <- do.call(rbind, metrics)
        sensitivity_sign <- do.call(rbind, sensitivity)
        rownames(metrics) <- NULL
        rownames(sensitivity_sign) <- NULL
        results_loco <- list(companies = companies, grid_s = grid_s, grid_t = grid_t, beta_full = beta_full,
            surfaces_loco = surfaces_loco, metrics = metrics, sensitivity_sign = sensitivity_sign,
            lambda_s = lambda_s, lambda_t = lambda_t, lambda_alpha = lambda_alpha, reconstruction = "Omit one company, recenter the 12 remaining curve pairs, and refit with fixed bases and penalties")
        results_loco
    }, envir = e)
}


loco_results <- refit_loco(environment())

# GARCH comparison


set.seed(123)

df_volume_log_standardized <- df_volume_log %>% dplyr::mutate(dplyr::across(dplyr::all_of(volume_cols),
    ~as.numeric(scale(.x))))

names(df_returns_final) <- sub("^R_Close_", "R_", names(df_returns_final))

companies_garch <- sub("^R_", "", grep("^R_", names(df_returns_final), value = TRUE))

candidates_pq <- data.frame(p = c(1, 1, 2, 2), q = c(1, 2, 1, 2))

prepare_data_garch <- function(company, df_returns = df_returns_final, df_volume = df_volume_log_standardized) {
    column_return <- paste0("R_", company)
    column_volume <- paste0("Volume_", company)
    if (!column_return %in% names(df_returns)) {
        stop("Missing return column: ", column_return)
    }
    if (!column_volume %in% names(df_volume)) {
        stop("Missing volume column: ", column_volume)
    }
    data_r <- df_returns %>% dplyr::transmute(date = as.Date(date), Return = .data[[column_return]])
    data_v <- df_volume %>% dplyr::transmute(date = as.Date(date), Volume_log_z = .data[[column_volume]])
    data <- dplyr::inner_join(data_r, data_v, by = "date") %>% dplyr::arrange(date)
    if (anyNA(data)) {
        stop("The aligned data for ", company, " contain missing values.")
    }
    if (anyDuplicated(data$date)) {
        stop("Duplicate dates in the aligned data for ", company, ".")
    }
    data
}

table_grid_garch <- dplyr::bind_rows(lapply(companies_garch, function(company) {
    dat <- prepare_data_garch(company)
    data.frame(Company = company, N = nrow(dat), date_initial = min(dat$date), date_final = max(dat$date),
        row.names = NULL)
}))

fit_one_garch <- function(company, p, q, ar = 0, include_volume = TRUE, lag_lb = 10) {
    data <- prepare_data_garch(company)
    regressor <- NULL
    if (include_volume) {
        regressor <- matrix(data$Volume_log_z, ncol = 1)
        colnames(regressor) <- "Volume_log_z"
    }
    spec <- rugarch::ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(q, p)), mean.model = list(armaOrder = c(ar,
        0), include.mean = TRUE, external.regressors = regressor), distribution.model = "std")
    fit <- tryCatch(rugarch::ugarchfit(spec = spec, data = data$Return, solver = "hybrid", fit.control = list(stationarity = 1)),
        error = function(e) e)
    if (inherits(fit, "error")) {
        return(list(company = company, p = p, q = q, ar = ar, include_volume = include_volume, data = data,
            fit = NULL, error = conditionMessage(fit), lag_lb = lag_lb))
    }
    list(company = company, p = p, q = q, ar = ar, include_volume = include_volume, data = data,
        fit = fit, error = NA_character_, lag_lb = lag_lb)
}

extract_summary_garch <- function(obj) {
    if (is.null(obj$fit)) {
        return(data.frame(Company = obj$company, p = obj$p, q = obj$q, AR = obj$ar, N = nrow(obj$data),
            Mu = NA_real_, Delta_Volume = NA_real_, P_Volume_Robust = NA_real_, Omega = NA_real_,
            Persistence = NA_real_, Shape_t = NA_real_, LogLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
            P_LB_Residuals = NA_real_, P_LB_Squared = NA_real_, Convergence = NA_integer_, Error = obj$error,
            row.names = NULL))
    }
    fit <- obj$fit
    fit_slot <- methods::slot(fit, "fit")
    cf <- fit_slot[["coef"]]
    robust <- fit_slot[["robust.matcoef"]]
    residuals <- as.numeric(fit_slot[["residuals"]])
    sigma <- as.numeric(fit_slot[["sigma"]])
    llh <- as.numeric(fit_slot[["LLH"]])
    ipars <- fit_slot[["ipars"]]
    convergence <- as.integer(fit_slot[["convergence"]])
    if (is.null(cf) || length(cf) == 0) {
        stop("Could not extract coefficients from the fit slot for ", obj$company)
    }
    alpha_names <- grep("^alpha[0-9]+$", names(cf), value = TRUE)
    beta_names <- grep("^beta[0-9]+$", names(cf), value = TRUE)
    persistence <- sum(c(cf[alpha_names], cf[beta_names]), na.rm = TRUE)
    z <- residuals/sigma
    z <- z[is.finite(z)]
    lb_res <- stats::Box.test(z, lag = obj$lag_lb, type = "Ljung-Box", fitdf = obj$ar)
    lb_sq <- stats::Box.test(z^2, lag = obj$lag_lb, type = "Ljung-Box", fitdf = 0)
    get <- function(name) {
        if (name %in% names(cf))
            as.numeric(cf[name])
        else NA_real_
    }
    p_robust <- function(name) {
        if (!is.null(robust) && name %in% rownames(robust) && ncol(robust) >= 4) {
            as.numeric(robust[name, 4])
        }
        else {
            NA_real_
        }
    }
    n_obs <- length(residuals)
    n_par <- if (!is.null(ipars) && ncol(ipars) >= 4) {
        sum(ipars[, 4], na.rm = TRUE)
    }
    else {
        length(cf)
    }
    aic <- (-2 * llh)/n_obs + (2 * n_par)/n_obs
    bic <- (-2 * llh)/n_obs + (n_par * log(n_obs))/n_obs
    data.frame(Company = obj$company, p = obj$p, q = obj$q, AR = obj$ar, N = nrow(obj$data), Mu = get("mu"),
        Delta_Volume = get("mxreg1"), P_Volume_Robust = p_robust("mxreg1"), Omega = get("omega"),
        Persistence = persistence, Shape_t = get("shape"), LogLik = llh, AIC = aic, BIC = bic,
        P_LB_Residuals = as.numeric(lb_res$p.value), P_LB_Squared = as.numeric(lb_sq$p.value), Convergence = convergence,
        Error = NA_character_, row.names = NULL)
}

extract_parameters_variance <- function(obj) {
    if (is.null(obj$fit))
        return(NULL)
    fit_slot <- methods::slot(obj$fit, "fit")
    cf <- fit_slot[["coef"]]
    names <- c("omega", grep("^alpha[0-9]+$", names(cf), value = TRUE), grep("^beta[0-9]+$", names(cf),
        value = TRUE))
    names <- names[names %in% names(cf)]
    data.frame(Company = obj$company, p = obj$p, q = obj$q, AR = obj$ar, Parameter = names, Estimate = as.numeric(cf[names]),
        row.names = NULL)
}

select_garch_model <- function(company) {
    objects_ar0 <- lapply(seq_len(nrow(candidates_pq)), function(k) {
        fit_one_garch(company = company, p = candidates_pq$p[k], q = candidates_pq$q[k], ar = 0, include_volume = TRUE)
    })
    summary_ar0 <- dplyr::bind_rows(lapply(objects_ar0, extract_summary_garch)) %>% dplyr::mutate(Adequate = Convergence ==
        0 & is.finite(Persistence) & Persistence < 1 & P_LB_Residuals > 0.05 & P_LB_Squared > 0.05)
    if (any(summary_ar0$Adequate, na.rm = TRUE)) {
        idx <- which(summary_ar0$Adequate)[which.min(summary_ar0$BIC[summary_ar0$Adequate])]
        return(list(company = company, candidates = objects_ar0, table = summary_ar0, final = objects_ar0[[idx]],
            final_summary = summary_ar0[idx, , drop = FALSE] %>% dplyr::mutate(Status_final = "Adequate")))
    }
    objects_ar1 <- lapply(seq_len(nrow(candidates_pq)), function(k) {
        fit_one_garch(company = company, p = candidates_pq$p[k], q = candidates_pq$q[k], ar = 1, include_volume = TRUE)
    })
    summary_ar1 <- dplyr::bind_rows(lapply(objects_ar1, extract_summary_garch)) %>% dplyr::mutate(Adequate = Convergence ==
        0 & is.finite(Persistence) & Persistence < 1 & P_LB_Residuals > 0.05 & P_LB_Squared > 0.05)
    objects_all <- c(objects_ar0, objects_ar1)
    summary_all <- dplyr::bind_rows(summary_ar0, summary_ar1)
    if (any(summary_all$Adequate, na.rm = TRUE)) {
        adequate <- which(summary_all$Adequate)
        idx <- adequate[which.min(summary_all$BIC[adequate])]
        return(list(company = company, candidates = objects_all, table = summary_all, final = objects_all[[idx]],
            final_summary = summary_all[idx, , drop = FALSE] %>% dplyr::mutate(Status_final = "Adequate")))
    }
    eligible <- which(summary_all$Convergence == 0 & is.finite(summary_all$BIC) & is.finite(summary_all$Persistence) &
        summary_all$Persistence < 1)
    if (length(eligible) == 0) {
        stop("No admissible model converged for ", company)
    }
    idx <- eligible[which.min(summary_all$BIC[eligible])]
    list(company = company, candidates = objects_all, table = summary_all, final = objects_all[[idx]],
        final_summary = summary_all[idx, , drop = FALSE] %>% dplyr::mutate(Status_final = "Requires caution"))
}

company_test_garch <- companies_garch[1]

message("Checking the initial GARCH fit for ", company_test_garch)

obj_test_garch <- fit_one_garch(company = company_test_garch, p = 1, q = 1, ar = 0, include_volume = TRUE)

if (is.null(obj_test_garch$fit)) {
    stop("ugarchfit failed in the initial check: ", obj_test_garch$error)
}

message("Fit class: ", paste(class(obj_test_garch$fit), collapse = ", "))

summary_test_garch <- extract_summary_garch(obj_test_garch)

message("Initial GARCH check passed. Selecting models for all companies.")

results_selection <- setNames(lapply(companies_garch, select_garch_model), companies_garch)

table_candidates_garch <- dplyr::bind_rows(lapply(results_selection, function(x) x$table))

table_models_garch_final <- dplyr::bind_rows(lapply(results_selection, function(x) x$final_summary)) %>%
    dplyr::mutate(Company = dplyr::recode(Company, !!!company_labels), Model = paste0("AR(", AR,
        ")-GARCH(", p, ",", q, ")")) %>% dplyr::select(Company, Model, N, BIC, AIC, Persistence, Shape_t,
    P_LB_Residuals, P_LB_Squared, Status_final)

table_parameters_garch_final <- dplyr::bind_rows(lapply(results_selection, function(x) {
    extract_parameters_variance(x$final)
})) %>% dplyr::mutate(Company = dplyr::recode(Company, !!!company_labels), Model = paste0("AR(",
    AR, ")-GARCH(", p, ",", q, ")"))

compare_volume_model_final <- function(result_company) {
    final <- result_company$final
    with_volume <- final
    without_volume <- fit_one_garch(company = final$company, p = final$p, q = final$q, ar = final$ar, include_volume = FALSE)
    if (is.null(with_volume$fit) || is.null(without_volume$fit)) {
        return(data.frame(Company = final$company, p = final$p, q = final$q, AR = final$ar, Delta_Volume = NA_real_,
            P_Volume_Robust = NA_real_, Delta_AIC = NA_real_, Delta_BIC = NA_real_, LR = NA_real_,
            P_LR = NA_real_, row.names = NULL))
    }
    fit1 <- with_volume$fit
    fit0 <- without_volume$fit
    fit1_slot <- methods::slot(fit1, "fit")
    fit0_slot <- methods::slot(fit0, "fit")
    cf1 <- fit1_slot[["coef"]]
    robust <- fit1_slot[["robust.matcoef"]]
    delta <- if ("mxreg1" %in% names(cf1))
        as.numeric(cf1["mxreg1"])
    else NA_real_
    p_delta <- if (!is.null(robust) && "mxreg1" %in% rownames(robust) && ncol(robust) >= 4) {
        as.numeric(robust["mxreg1", 4])
    }
    else {
        NA_real_
    }
    ll0 <- as.numeric(fit0_slot[["LLH"]])
    ll1 <- as.numeric(fit1_slot[["LLH"]])
    res0 <- as.numeric(fit0_slot[["residuals"]])
    res1 <- as.numeric(fit1_slot[["residuals"]])
    ipars0 <- fit0_slot[["ipars"]]
    ipars1 <- fit1_slot[["ipars"]]
    n0 <- length(res0)
    n1 <- length(res1)
    k0 <- if (!is.null(ipars0) && ncol(ipars0) >= 4)
        sum(ipars0[, 4], na.rm = TRUE)
    else length(fit0_slot[["coef"]])
    k1 <- if (!is.null(ipars1) && ncol(ipars1) >= 4)
        sum(ipars1[, 4], na.rm = TRUE)
    else length(cf1)
    aic0 <- (-2 * ll0)/n0 + (2 * k0)/n0
    aic1 <- (-2 * ll1)/n1 + (2 * k1)/n1
    bic0 <- (-2 * ll0)/n0 + (k0 * log(n0))/n0
    bic1 <- (-2 * ll1)/n1 + (k1 * log(n1))/n1
    LR <- 2 * (ll1 - ll0)
    p_lr <- if (is.finite(LR) && LR >= 0) {
        stats::pchisq(LR, df = 1, lower.tail = FALSE)
    }
    else {
        NA_real_
    }
    data.frame(Company = final$company, p = final$p, q = final$q, AR = final$ar, Delta_Volume = delta,
        P_Volume_Robust = p_delta, Delta_AIC = aic1 - aic0, Delta_BIC = bic1 - bic0, LR = LR, P_LR = p_lr,
        row.names = NULL)
}

table_comparison_volume <- dplyr::bind_rows(lapply(results_selection, compare_volume_model_final)) %>%
    dplyr::mutate(Company = dplyr::recode(Company, !!!company_labels), Model = paste0("AR(", AR,
        ")-GARCH(", p, ",", q, ")"), Evidence_volume = dplyr::case_when(P_Volume_Robust < 0.05 &
        Delta_AIC < 0 & Delta_BIC < 0 ~ "Strong evidence", P_Volume_Robust < 0.05 & (Delta_AIC < 0 |
        Delta_BIC < 0) ~ "Partial evidence", P_Volume_Robust < 0.1 ~ "Weak evidence", TRUE ~ "No statistical evidence")) %>%
    dplyr::select(Company, Model, Delta_Volume, P_Volume_Robust, Delta_AIC, Delta_BIC, LR, P_LR,
        Evidence_volume)

extract_signbias <- function(result_company) {
    fit <- result_company$final$fit
    if (is.null(fit)) {
        return(NULL)
    }
    sb <- rugarch::signbias(fit)
    data.frame(Company = result_company$company, P_Sign_Bias = as.numeric(sb["Sign Bias", "prob"]),
        P_Negative_Sign_Bias = as.numeric(sb["Negative Sign Bias", "prob"]), P_Positive_Sign_Bias = as.numeric(sb["Positive Sign Bias",
            "prob"]), P_Joint_Effect = as.numeric(sb["Joint Effect", "prob"]), row.names = NULL)
}

table_symmetry_garch <- dplyr::bind_rows(lapply(results_selection, extract_signbias))

table_symmetry_garch <- table_symmetry_garch %>% dplyr::mutate(Company = dplyr::recode(Company, !!!company_labels),
    Conclusion_Symmetry = dplyr::case_when(P_Joint_Effect < 0.01 ~ "Strong evidence of asymmetry",
        P_Joint_Effect < 0.05 ~ "Evidence of asymmetry", P_Joint_Effect < 0.1 ~ "Weak evidence of asymmetry",
        TRUE ~ "No evidence of asymmetry"))


# Export EPS figures under descriptive names.
export_final_figures <- function(e, loco) {
  validate_numbers_final(list(e = e, loco = loco))
  validate_evaluation_fd(e)
  output_directory <- file.path("output", "figures")
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  draw_curve <- function(mat, meanfd, label) {
    base_style()
    matplot(seq_len(nrow(mat)), mat, type = "l", lty = 1, lwd = .65,
            col = col_companies(ncol(mat)), xlab = "Trading-session index", ylab = label)
    lines(seq_len(nrow(mat)), as.vector(evaluate_fd_stored(meanfd, seq_len(nrow(mat)))),
          lwd = 1.8, lty = 2, col = "black")
    legend("topright", "Cross-company mean", lty = 2, lwd = 1.8, bty = "n", cex = .8)
  }
  B <- e$beta_final_mat
  counts <- apply(sign(loco$surfaces_loco) == array(sign(B), dim = dim(loco$surfaces_loco)), c(1, 2), sum)
  draws <- list(
    smoothed_log_returns = function() draw_curve(e$matrix_return_smoothed, e$mean_returns, "Smoothed log-return"),
    smoothed_log_volume = function() draw_curve(e$matrix_volume_smoothed, e$mean_volume, "Smoothed log(1 + volume)"),
    coefficient_contour = function() draw_map(B, e$grid_beta_s_final, e$grid_beta_t_final, e$dates_s, e$events_plot, diagonal = TRUE),
    loco_influence_summary = function() draw_loco_compact(loco),
    extreme_log_returns = function() draw_extremes(e, "return"),
    highest_trading_volumes = function() draw_extremes(e, "high"),
    lowest_positive_trading_volumes = function() draw_extremes(e, "low"),
    volume_gcv_by_penalty = function() {
      d <- e$results_volume_list[["24"]]; base_style()
      plot(d$log10lambda_grid, d$gcv_vec, type = "b", pch = 16, cex = .65,
           xlab = expression(log[10](lambda)), ylab = "Total GCV")
      abline(v = d$best_log10lambda, lty = 2); mtext("K = 24", 3, line = .25, cex = .9)
    },
    volume_gcv_by_basis = function() {
      d <- e$summary_volume_K; base_style()
      plot(d$nbasis, d$best_gcv, type = "b", pch = 16, cex = .6,
           xlab = "Number of basis functions K", ylab = "Minimum total GCV over penalty grid")
      abline(v = e$best_nb_volume, lty = 2)
    },
    volume_smoothing_diagnostics = function() draw_smoothing_panels(e, TRUE),
    return_gcv_by_penalty = function() {
      base_style()
      plot(e$lambdas_log_return, e$gcv_return_preliminary, type = "b", pch = 16, cex = .6,
           xlab = expression(log[10](lambda)), ylab = "Total GCV")
      abline(v = e$lambdas_log_return[e$idx_lambda_return_preliminary], lty = 2)
    },
    return_gcv_by_basis = function() {
      d <- e$results_K_return; base_style()
      plot(d$K_Bases, (d$Min_GCV - .00834) * 1e6, type = "b", pch = 16, cex = .6,
           xlab = "Number of basis functions K", ylab = expression((GCV - .00834) %*% 10^6))
      abline(v = e$best_K_return_gcv, lty = 2); abline(v = 150, col = "#B2182B", lty = 3)
      legend("topright", c("Grid minimum", "Final K = 150"), lty = c(2, 3), col = c("black", "#B2182B"), bty = "n", cex = .8)
    },
    return_basis_comparison = function() {
      par(mfrow = c(1, 2))
      for (k in c("RMS_Rel_Mean", "GCV_Total")) {
        base_style(c(4.2, 5.8, 2, 1)); d <- e$results_bases_return
        plot(d$N_Bases, d[[k]], type = "b", pch = 16, cex = .6, xlab = "Number of basis functions", ylab = "")
        mtext(if (k == "GCV_Total") "Total GCV" else "Mean relative RMS", 2, line = 4.1, las = 0)
      }
    },
    return_smoothing_diagnostics = function() draw_smoothing_panels(e, FALSE),
    coefficient_penalty_sensitivity = function() draw_sensitivity_final(e),
    coefficient_surface = function() draw_surface(B, e$grid_beta_s_final, e$grid_beta_t_final),
    coefficient_profiles = function() {
      par(mfrow = c(2, 1)); draw_event_curve(e, "mean", "(a)"); draw_event_curve(e, "rms", "(b)")
    },
    pointwise_r_squared = function() draw_r2(as.numeric(e$dates_t), e$R2_t_final, e$dates_t),
    residuals_vs_fitted = function() draw_residual_fit_final(e),
    residual_correlations = function() draw_corr_final(e),
    residual_qq = function() {
      base_style()
      qqnorm(e$res_finite, main = "", xlab = "Theoretical normal quantile", ylab = "Sample residual quantile", pch = 16, cex = .4, col = "#526779")
      qqline(e$res_finite, col = "#B2182B", lwd = 1)
    },
    residual_histogram = function() {
      base_style()
      hist(e$res_finite, breaks = "FD", col = "#D8E1E8", border = "#444444", main = "", xlab = "Residual", ylab = "Frequency")
      abline(v = 0, lty = 2, lwd = .8)
    },
    loco_influence_metrics = function() draw_loco_full_final(loco),
    loco_sign_agreement_map = function() draw_map(counts, e$grid_beta_s_final, e$grid_beta_t_final,
       e$dates_s, count = TRUE, active = abs(B), label = "Refits\n(out of 13)")
  )
  dimensions <- rbind(
    c(143.4, 95),
    c(143.4, 95),
    c(143.4, 137),
    c(143.4, 103),
    c(143.4, 115),
    c(143.4, 115),
    c(143.4, 115),
    c(110, 85),
    c(120, 88),
    c(170, 150),
    c(110, 85),
    c(120, 88),
    c(170, 95),
    c(170, 150),
    c(143.4, 155),
    c(143.4, 112),
    c(143.4, 160),
    c(143.4, 90),
    c(143.4, 95),
    c(143.4, 137),
    c(110, 95),
    c(110, 85),
    c(143.4, 170),
    c(143.4, 137)
  )
  stopifnot(length(draws) == 24L)
  for (i in seq_along(draws)) {
    path <- file.path(output_directory, paste0(names(draws)[i], ".eps"))
    message("Exporting ", path)
    with_device(path, "eps", dimensions[i, 1], dimensions[i, 2], draws[[i]])
    # Remove identifying exporter metadata while preserving the drawing instructions.
    lines <- readLines(path, warn = FALSE)
    drop <- grepl("^%%(Creator|CreationDate|Title|For|Invocation):", lines)
    continuation <- FALSE
    for (j in seq_along(lines)) {
      if (grepl("^%%Invocation:", lines[j])) continuation <- TRUE
      else if (continuation && startsWith(lines[j], "%%+")) drop[j] <- TRUE
      else continuation <- FALSE
    }
    writeLines(lines[!drop], path, useBytes = TRUE)
    invisible(verify_eps(path, dimensions[i, 1], dimensions[i, 2]))
  }
  message("Verified 24 EPS figures in output/figures/. Keep generated outputs local.")
  invisible(NULL)
}
export_final_figures(environment(), loco_results)
analysis_complete <- TRUE
print_analysis_summary(environment())
if ("--export-tables" %in% commandArgs(trailingOnly = TRUE)) {
    export_reported_tables(environment())
} else {
    message("For local numerical summaries, call export_reported_tables() or rerun with --export-tables.")
}
message(sprintf("Finished in %.1f minutes. Computational workflow reproducible subject to third-party data access restrictions.",
                (proc.time()[["elapsed"]] - analysis_started) / 60))
