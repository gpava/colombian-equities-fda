# Temporal Heterogeneity in the Trading Volume–Return Relationship in Colombian Equities: A Function-on-Function Regression Approach

This is the R software for a study of 13 Colombian equities. It uses
B-splines to represent each company's log returns and log trading volume
as curves, then fits a function-on-function regression to examine how the
volume–return association varies across trading sessions. The workflow also
compares three penalty settings, refits the model after omitting each company
in turn, and fits company-level GARCH models.

**Computational workflow reproducible subject to third-party data access restrictions.**

## Repository files

| File | Purpose |
|---|---|
| `run_analysis.R` | Read the data, fit the models, check the results and export the figures |
| `R/essential_functions.R` | Draw the figures and check their numerical content, layout and EPS format |
| `R/reported_tables.R` | Print a short summary and export numerical results on request |
| `renv.lock` | Record the R version and 101 package versions |
| `DATA_ACCESS.md` | Describe the data source, access conditions, input format and preprocessing |
| `CITATION.cff`, `LICENSE`, `.gitignore` | Provide the software citation, license and rules for local files |

The release contains code and documentation. All figures and numerical
tables are generated on the reader's machine. Some graphs display observed
or transformed series. Keep the inputs and generated outputs in your local
working copy.

## Install the software

Use **R 4.4.3** and **Ghostscript**. The script looks for `gs` on macOS and
Linux and for `gswin64c` or `gswin32c` on Windows. You can also set
`R_GSCMD` to the Ghostscript executable. On macOS with Homebrew,
`brew install ghostscript` installs it. On Debian or Ubuntu, install the
`ghostscript` system package.

In **RStudio**, start a fresh R session and choose **Session → Set Working
Directory → Choose Directory…**. Select the folder that contains
`run_analysis.R`. Run the following commands once in the R console. The
package library stays inside your working copy.

```r
stopifnot(file.exists("run_analysis.R"), getRversion() == "4.4.3")
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
project_library <- file.path(
  getwd(), "renv", "library", paste0("R-", getRversion()), R.version$platform
)
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
renv::restore(project = getwd(), library = project_library,
              lockfile = "renv.lock", prompt = FALSE)
```

Restart R after restoring the packages and check that the working directory
still contains `run_analysis.R`. Restoration needs internet access. If your
platform installs packages from source, you may need C, C++ and Fortran
compilers and development libraries for individual packages. The script
checks the versions in the lockfile before fitting the models.

## Run in RStudio

Obtain the 13 CSVs under the provider's terms and place them in
`data/raw/`. The filenames, headers and locale conventions are listed in
[DATA_ACCESS.md](DATA_ACCESS.md).

```r
source("run_analysis.R", encoding = "UTF-8")
```

Before adding the inputs, you can run the same command as a quick check.
It lists the missing paths and directs you to `DATA_ACCESS.md` before
loading any packages. With the inputs present, it prints a summary of the
fit and reports `Verified 24 EPS figures` when all figures pass their checks.
The script creates `output/figures/` and replaces its own named files on
subsequent runs.

To export the numerical results, run this in the **same R session** after
the analysis finishes.

```r
export_reported_tables()
```

For a single, noninteractive run from the repository root, use the terminal.

```sh
Rscript --vanilla run_analysis.R --export-tables
```

Omit `--export-tables` to generate only the figures and console summary.
The optional export writes 24 CSVs and `sessionInfo.txt` to `output/tables/`.
The CSVs use UTF-8, comma-delimited fields, decimal points, ISO dates and
`NA` for unavailable values.

## Model settings and interpretation

The univariate curve-smoothing searches use generalized cross-validation
(GCV). For the functional regression, we compared coefficient surfaces for
three pairs of penalties, in the order (volume direction, return direction).
We evaluated A (10², 10⁷), B (10³, 10⁸) and C (10⁴, 10⁹) on a common
150 × 150 grid and selected B by inspection. The script reproduces this
choice for the final model and the 13 leave-one-company-out (LOCO) refits.
The final coefficient surface is evaluated on a 490 × 490 grid.

The mean pointwise R² for the final fit uses 490 return sessions. The
sensitivity table summarizes a variance ratio over 150 evaluation points.
The two summaries describe the fit on different grids and use different
calculations.

## Numerical outputs

Match rows by company, event, model order or smoothing parameter when
comparing runs. The session report records the platform, packages and
numerical libraries used.

| CSV in `output/tables/` | Contents |
|---|---|
| `analysis_counts_and_settings.csv` | Sample counts, selected penalties and basis sizes, and diagnostic counts |
| `analysis_metric_ranges.csv` | Ranges and central values for coefficients, fit, LOCO and GARCH persistence |
| `data_coverage.csv` | Data coverage, dates added and imputation by company |
| `extreme_event_matches.csv` | Selected event dates matched to extreme returns and volumes |
| `event_session_matches.csv` | Event dates, matched trading sessions and exact-match indicators |
| `event_coefficient_summary.csv` | Mean and RMS coefficients by event, with dates and descriptions |
| `event_coefficient_metrics.csv` | Mean, RMS, positive share and contemporaneous coefficient by event |
| `top_volume_gcv_by_basis.csv` | Five best volume basis sizes, each evaluated at its GCV-minimizing penalty |
| `return_penalty_comparison.csv` | Return smoothing diagnostics across penalties at K = 150 |
| `return_basis_comparison.csv` | Return smoothing diagnostics across basis sizes at a fixed penalty |
| `smoothing_metrics.csv` | Smoothing error and reconstruction index by variable and company |
| `regression_penalty_scenarios.csv` | Paired regression penalties and mean pointwise R² for each scenario |
| `company_r_squared.csv` | Company-level functional R² and sums of squares |
| `residual_magnitude.csv` | Residual magnitude and variance by company |
| `residual_statistics.csv` | Pooled residual statistics |
| `functional_model_diagnostics.csv` | Functional model checks and related outputs |
| `garch_order_comparison.csv` | AR(0) GARCH candidate AIC and BIC |
| `garch_ar1_diagnostics.csv` | Diagnostics for the AR(1) candidates evaluated |
| `garch_selected_models.csv` | Selected GARCH models and fit diagnostics |
| `garch_variance_parameters.csv` | Conditional variance parameters |
| `garch_sign_bias.csv` | Joint sign-bias tests |
| `garch_volume_comparison.csv` | GARCH volume contribution for all companies |
| `garch_volume_evidence.csv` | GARCH volume results classified as strong, partial or weak evidence |
| `loco_influence_metrics.csv` | Leave-one-company-out distances, correlations and sign agreement |

The reconstruction index is `1 - Relative_RMS^2`. GARCH(p,q) uses p GARCH
lags and q ARCH lags. AIC and BIC are per observation. A negative
likelihood-ratio statistic leaves its p-value unavailable. LOCO sign
agreement is reported in percentages.

## Figures generated locally

| Files | Content | Width × height (mm) |
|---|---|---|
| `smoothed_log_returns.eps` | Smoothed log returns and their cross-company mean | 143.4 × 95 |
| `smoothed_log_volume.eps` | Smoothed log volumes and their cross-company mean | 143.4 × 95 |
| `coefficient_contour.eps` | Coefficient contours and the nine event dates | 143.4 × 137 |
| `loco_influence_summary.eps` | Influence of each omitted company | 143.4 × 103 |
| `extreme_log_returns.eps` | Largest absolute log returns by company | 143.4 × 115 |
| `highest_trading_volumes.eps` | Highest recorded volumes by company | 143.4 × 115 |
| `lowest_positive_trading_volumes.eps` | Lowest positive recorded volumes by company | 143.4 × 115 |
| `volume_gcv_by_penalty.eps` | Volume GCV across penalties at 24 basis functions | 110 × 85 |
| `volume_gcv_by_basis.eps` | Best volume GCV at each basis size | 120 × 88 |
| `volume_smoothing_diagnostics.eps` | Log-volume curves and smoothing residuals for two companies | 170 × 150 |
| `return_gcv_by_penalty.eps` | Return GCV across penalties at 140 basis functions | 110 × 85 |
| `return_gcv_by_basis.eps` | Best return GCV at each basis size | 120 × 88 |
| `return_basis_comparison.eps` | Mean relative RMS and total return GCV across basis sizes | 170 × 95 |
| `return_smoothing_diagnostics.eps` | Log-return curves and smoothing residuals for two companies | 170 × 150 |
| `coefficient_penalty_sensitivity.eps` | Coefficient surfaces for the three penalty pairs | 143.4 × 155 |
| `coefficient_surface.eps` | Estimated volume–return coefficient surface | 143.4 × 112 |
| `coefficient_profiles.eps` | Mean and root-mean-square coefficients across return sessions | 143.4 × 160 |
| `pointwise_r_squared.eps` | Within-sample R² at each return session | 143.4 × 90 |
| `residuals_vs_fitted.eps` | Functional residuals versus fitted values | 143.4 × 95 |
| `residual_correlations.eps` | Correlations between company residual curves | 143.4 × 137 |
| `residual_qq.eps` | Normal Q–Q plot of pooled residuals | 110 × 95 |
| `residual_histogram.eps` | Histogram of pooled residuals | 110 × 85 |
| `loco_influence_metrics.eps` | Four measures of influence for each omitted company | 143.4 × 170 |
| `loco_sign_agreement_map.eps` | Number of refits that preserve the sign of the full coefficient surface | 143.4 × 137 |

Plot filenames describe the analysis and remain independent of manuscript
numbering. Ghostscript checks each EPS. The drawing code preserves the
original physical dimensions, colour scales and event annotations. The
compact LOCO plot also checks the clearance between the reference lines
and the value labels.

## Run time and numerical agreement

Allow a few minutes for the analysis once packages are installed. The
release took about two minutes with R 4.4.3 on Linux x86_64 and OpenBLAS
0.3.29. An earlier version took about 25 minutes on the author's Mac and
included an exploratory calculation that did not feed into the fitted
model. Hardware, numerical libraries and package installation affect
elapsed time.

On the same Linux system, the previous and revised code produced identical
files for all 24 EPS figures and all 24 numerical CSVs. The author’s macOS
run and the Linux run selected the same GARCH models and evidence categories.
Some GARCH estimates, diagnostic p-values and likelihood-ratio tests vary
across numerical environments. For Terpel, the macOS likelihood-ratio
statistic was negative and the associated p-value was unavailable. The
Linux fit returned a positive statistic. Compare the local
`sessionInfo.txt` files when investigating such differences.

## Data and code availability

The financial data were obtained from Investing.com. Redistribution of these
third-party inputs is restricted. Researchers can obtain the data directly
from the provider under the applicable terms. The corresponding author can
help with the expected format and reproduction workflow. Permission to
redistribute data rests with the relevant rights holder.

The analysis code used to produce the results and figures is available in
this repository. Reproduction of the empirical analysis requires the user
to supply the third-party input data described in [DATA_ACCESS.md](DATA_ACCESS.md).

## Authorship, citation and license

Gustavo A. Pava-Parra developed the software. David L. Ocampo R. and
Juan Carlos Riaño-Rojas reviewed and edited the accompanying manuscript.

Use `CITATION.cff` or GitHub's **Cite this repository** to cite the software.
When using an archived release, cite the DOI for that version in Zenodo.
The accompanying article has its own author list and citation.

The MIT license covers the R code and accompanying documentation. Financial
data from Investing.com are excluded. Dependencies retain their own licenses.
See [LICENSE](LICENSE).
