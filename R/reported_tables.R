# Export numerical summaries from the fitted analysis.

print_analysis_summary <- function(analysis) {
    with(analysis, {
        message(sprintf("Analysis summary: %d companies, %d aligned dates, %d return sessions.",
                        n_companies, nrow(dates_base), n_days_return))
        message(sprintf("Smoothing: volume K = %d, lambda = %g; returns K = %d, lambda = %g.",
                        best_nb_volume, best_lambda_volume, nbasis_return_final, lambda_return_final))
        message(sprintf("Functional fit: beta range [%+.6g, %+.6g]; mean pointwise R-squared = %.6f.",
                        limits_beta_final[1], limits_beta_final[2], mean_temporal_R2_final))
        message(sprintf("Matched event dates: %d/%d. Leave-one-company-out refits: %d. Adequate GARCH fits: %d/%d.",
                        sum(events_market$date == events_market$date_in_data), nrow(events_market),
                        nrow(loco_results$metrics),
                        sum(table_models_garch_final$Status_final == "Adequate"),
                        nrow(table_models_garch_final)))
    })
}

return_penalty_table <- function(analysis) {
    with(new.env(parent = analysis), {
        # Compare the low and selected penalties with the GCV minimum at K = 150.
        penalties <- c(lambda_return_low, lambda_return_final, row_K150$Best_Lambda)
        roughness_penalty <- fda::eval.penalty(base_return_final, fda::int2Lfd(2))
        do.call(rbind, lapply(penalties, function(lambda) {
            smooth <- if (lambda == lambda_return_final) {
                smooth_return_final
            } else {
                parameters <- fda::fdPar(base_return_final, Lfdobj = 2, lambda = lambda)
                fda::smooth.basis(time_return, matrix_returns, parameters)
            }
            metrics <- compute_metrics_fd(matrix_returns, smooth$fd, time_return)
            coefficients <- as.matrix(smooth$fd$coefs)
            roughness <- colSums(coefficients * (roughness_penalty %*% coefficients))
            data.frame(Lambda = lambda, Total_GCV = sum(smooth$gcv),
                       Effective_df = as.numeric(smooth$df)[1],
                       Mean_relative_RMS = mean(metrics$RMS_Relative, na.rm = TRUE),
                       Median_roughness = stats::median(roughness, na.rm = TRUE))
        }))
    })
}

collect_reported_tables <- function(analysis) {
    with(new.env(parent = analysis), {
        company_label <- function(keys) unname(company_labels[keys])
        company_ids <- names(company_labels)
        company_order <- unname(company_labels)

        original_dates <- vapply(company_series[company_ids], function(x) length(unique(x$date)), integer(1))
        coverage <- data.frame(
            Company = company_order,
            Original_dates = original_dates,
            Added_dates = nrow(dates_base) - original_dates,
            LOCF = colSums(is.na(df_prices[paste0("Close_", company_ids)])),
            GAM_percent = 100 * colSums(is.na(df_volume[paste0("Volume_", company_ids)]) &
                                           !is.na(df_volume_imputed[paste0("Volume_", company_ids)])) /
                nrow(df_volume)
        )

        event_columns <- events_plot[c("Code", "date", "date_in_data", "Index_s", "Event",
                                       "Beta_mean", "Beta_RMS", "Proportion_beta_positive",
                                       "Beta_contemporaneous")]
        names(event_columns)[1:4] <- c("Event_ID", "Event_date", "Matched_date", "Session_index")
        event_locations <- event_columns[c("Event_ID", "Event_date", "Matched_date", "Session_index", "Event")]
        event_locations$Exact_match <- event_locations$Event_date == event_locations$Matched_date
        # Match the selected events to the extreme-return and volume groups.
        context_ids <- c("suramericana", "bogota", "nutresa", "nutresa", "bvc")
        extreme_cases <- event_locations[c(1L, 2L, 4L, 5L, 7L), c("Event_ID", "Event_date", "Event")]
        extreme_cases$Company <- company_label(context_ids)
        is_selected <- function(group) {
            paste(extreme_cases$Event_date, unname(labels_review[context_ids])) %in%
                paste(group$date, group$Company)
        }
        extreme_cases$Top_absolute_return <- is_selected(table_returns_extremes)
        extreme_cases$Top_observed_volume <- is_selected(table_volume_high)
        extreme_cases$Bottom_positive_observed_volume <- is_selected(table_volume_low)

        volume_search <- utils::head(summary_volume[c("nbasis", "best_log10lambda", "best_lambda", "best_gcv")], 5)
        names(volume_search) <- c("Basis_size", "Log10_lambda", "Lambda", "Total_GCV")
        smoothing_metrics <- function(x, variable) {
            data.frame(Variable = variable, Company = company_label(x$Series),
                       RMS = x$RMS, SD = x$SD_Original, Relative_RMS = x$RMS_Relative,
                       Reconstruction_index = 1 - x$RMS_Relative^2)
        }
        smoothing <- rbind(smoothing_metrics(table_rms_volume, "Log-volume"),
                           smoothing_metrics(table_performance_return, "Log returns"))
        company_r2 <- table_R2_verification
        company_r2$Company <- company_label(company_r2$Company)
        residual_magnitude <- table_var_residual
        residual_magnitude$Company <- company_label(residual_magnitude$Company)
        residual_correlations <- cor_residual_companies[upper.tri(cor_residual_companies)]
        assumptions <- data.frame(
            Assumption = c("Linearity", "Zero mean", "Common covariance",
                           "Independence across curves", "Normality"),
            Diagnostic = c("Residuals versus fitted smoother", "Maximum absolute cross-company mean",
                           "Ratio of largest to smallest residual variance",
                           "Maximum absolute off-diagonal Pearson correlation",
                           "Pooled excess kurtosis"),
            Value = c(NA_real_, max(abs(mean_res_t)),
                      max(residual_magnitude$Variance_residual) / min(residual_magnitude$Variance_residual),
                      max(abs(residual_correlations)), table_normality_res$Excess_Kurtosis),
            Reference = c("output/figures/residuals_vs_fitted.eps",
                          "Cross-company means of evaluated residuals",
                          "output/tables/residual_magnitude.csv",
                          "output/figures/residual_correlations.eps",
                          paste("output/tables/residual_statistics.csv", "output/figures/residual_qq.eps",
                                "output/figures/residual_histogram.eps", sep = " | "))
        )

        candidates <- table_candidates_garch
        candidates$Company <- company_label(candidates$Company)
        orders <- candidates[candidates$AR == 0L, c("Company", "p", "q", "AIC", "BIC")]
        order_groups <- split(seq_len(nrow(orders)), orders$Company)
        orders$Minimum_BIC <- FALSE
        for (rows in order_groups) {
            finite <- rows[is.finite(orders$BIC[rows])]
            orders$Minimum_BIC[finite[which.min(orders$BIC[finite])]] <- TRUE
        }
        extensions <- candidates[candidates$AR == 1L,
                                 c("Company", "AR", "p", "q", "BIC", "Persistence",
                                   "P_LB_Residuals", "P_LB_Squared", "Adequate", "Convergence")]
        selected <- table_models_garch_final
        selected$Status_final[selected$Status_final == "Requires caution"] <- "Caution"
        parameters <- table_parameters_garch_final
        variance <- selected[c("Company", "Model", "Persistence")]
        for (parameter in unique(c("omega", "alpha1", "beta1", "beta2", parameters$Parameter))) {
            estimates <- parameters[parameters$Parameter == parameter, ]
            variance[[parameter]] <- estimates$Estimate[match(variance$Company, estimates$Company)]
        }
        variance <- variance[c("Company", "Model", setdiff(names(variance), c("Company", "Model", "Persistence")),
                               "Persistence")]
        volume <- table_comparison_volume
        volume$Evidence_volume <- sub(" evidence$", "", volume$Evidence_volume)
        volume$Evidence_volume[volume$Evidence_volume == "No statistical"] <- "No consistent evidence"
        evidence_order <- match(volume$Evidence_volume, c("Strong", "Partial", "Weak"))
        volume_evidence <- volume[which(!is.na(evidence_order)), c("Company", "Delta_Volume", "P_Volume_Robust",
                                                                   "Evidence_volume")]
        volume_evidence <- volume_evidence[order(evidence_order[!is.na(evidence_order)]), ]

        loco <- loco_results$metrics
        loco <- loco[order(loco$Distance_L2_relative, decreasing = TRUE), ]
        sensitivity <- loco_results$sensitivity_sign
        active_10 <- sensitivity[sensitivity$Threshold_relative == 0.10, ]
        agreement_10 <- active_10$Agreement_sign_area_active[match(loco$Company_omitted, active_10$Company_omitted)]
        loco_table <- data.frame(
            Omitted_company = company_label(loco$Company_omitted),
            Relative_L2 = loco$Distance_L2_relative,
            Correlation = loco$Correlation_surfaces,
            Global_sign_percent = 100 * loco$Agreement_sign_global,
            Active_5_percent = 100 * loco$Agreement_sign_region_active_5pct,
            Active_10_percent = 100 * agreement_10
        )
        ranges <- data.frame(
            Metric = c("Coefficient surface", "Pointwise R-squared", "Integrated company R-squared",
                       "LOCO relative L2", "LOCO surface correlation", "LOCO global sign (%)",
                       "LOCO active 5% sign (%)", "LOCO active 10% sign (%)", "GARCH persistence"),
            Minimum = c(limits_beta_final[1], min(R2_t_final, na.rm = TRUE), min(company_r2$R2_i),
                        min(loco_table$Relative_L2), min(loco_table$Correlation),
                        min(loco_table$Global_sign_percent), min(loco_table$Active_5_percent),
                        min(loco_table$Active_10_percent), min(selected$Persistence)),
            Maximum = c(limits_beta_final[2], max(R2_t_final, na.rm = TRUE), max(company_r2$R2_i),
                        max(loco_table$Relative_L2), max(loco_table$Correlation),
                        max(loco_table$Global_sign_percent), max(loco_table$Active_5_percent),
                        max(loco_table$Active_10_percent), max(selected$Persistence)),
            Mean = c(NA_real_, mean_temporal_R2_final, rep(NA_real_, 7)),
            Median = c(rep(NA_real_, 6), stats::median(loco_table$Active_5_percent),
                       stats::median(loco_table$Active_10_percent), NA_real_)
        )
        counts <- data.frame(
            Companies = n_companies, Common_dates = nrow(dates_base), Return_sessions = n_days_return,
            Exact_event_dates = sum(event_locations$Exact_match),
            Volume_K = best_nb_volume, Volume_lambda = best_lambda_volume,
            Return_K = nbasis_return_final, Return_lambda = lambda_return_final,
            Preliminary_return_GCV_lambda = lambda_return_preliminary,
            Beta_lambda_s = lambda_s_final, Beta_lambda_t = lambda_t_final,
            Company_R2_below_0_2 = sum(company_r2$R2_i < 0.2),
            Company_R2_above_0_5 = sum(company_r2$R2_i > 0.5),
            Extreme_returns = nrow(table_returns_extremes), High_volumes = nrow(table_volume_high),
            Low_positive_volumes = nrow(table_volume_low), Observed_zero_volumes = sum(table_summary_zeros$Count_zeros),
            Residual_evaluations = length(res_finite), LOCO_fits = nrow(loco_table),
            Adequate_GARCH_fits = sum(selected$Status_final == "Adequate"),
            Best_AR0_BIC_GARCH_1_1 = sum(orders$Minimum_BIC & orders$p == 1 & orders$q == 1)
        )
        list(
            event_coefficient_summary = event_columns[c("Event_ID", "Event_date", "Event", "Beta_mean", "Beta_RMS")],
            garch_volume_evidence = volume_evidence,
            data_coverage = coverage,
            extreme_event_matches = extreme_cases,
            top_volume_gcv_by_basis = volume_search,
            return_penalty_comparison = return_penalty_table(environment()),
            smoothing_metrics = smoothing,
            event_session_matches = event_locations,
            event_coefficient_metrics = event_columns[c("Event_ID", "Beta_mean", "Beta_RMS",
                                                         "Proportion_beta_positive", "Beta_contemporaneous")],
            company_r_squared = company_r2,
            residual_magnitude = residual_magnitude,
            residual_statistics = table_normality_res,
            functional_model_diagnostics = assumptions,
            garch_order_comparison = orders,
            garch_ar1_diagnostics = extensions,
            garch_selected_models = selected,
            garch_variance_parameters = variance,
            garch_sign_bias = table_symmetry_garch[c("Company", "P_Joint_Effect", "Conclusion_Symmetry")],
            garch_volume_comparison = volume,
            loco_influence_metrics = loco_table,
            analysis_counts_and_settings = counts,
            analysis_metric_ranges = ranges,
            return_basis_comparison = results_bases_return,
            regression_penalty_scenarios = unique(df_sensitivity_plot[c("scenario", "lambda_s", "lambda_t", "r2_mean")])
        )
    })
}

# Call after source("run_analysis.R") in the same session.
export_reported_tables <- function(analysis = parent.frame()) {
    if (!isTRUE(get0("analysis_complete", envir = analysis, inherits = FALSE))) {
        stop("Run run_analysis.R successfully before exporting tables.", call. = FALSE)
    }
    tables <- collect_reported_tables(analysis)
    directory <- file.path("output", "tables")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) {
        stop("Could not create output/tables/.", call. = FALSE)
    }
    paths <- file.path(directory, paste0(names(tables), ".csv"))
    for (i in seq_along(tables)) {
        utils::write.csv(tables[[i]], paths[i], row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
    }
    session_path <- file.path(directory, "sessionInfo.txt")
    writeLines(capture.output(utils::sessionInfo()), session_path)
    message(sprintf("Exported %d local CSV summaries and sessionInfo.txt to output/tables/.", length(paths)))
    invisible(c(paths, session_path))
}
