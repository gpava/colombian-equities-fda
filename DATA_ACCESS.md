# Accessing the data

## Source and permissions

The study and its corresponding author identify Investing.com as the source
of the historical prices and trading volumes. The original CSVs contain no
provider metadata. These are third-party data, and restrictions apply to their
redistribution. This repository grants no rights over them.

Obtain the data directly from the provider under the terms that apply to
your use, or use a copy you are authorized to access. Keep the files local.
Redistributing them through GitHub, Zenodo or another archive requires
permission from the relevant rights holder, including for a deposit with
restricted access. See the [provider's terms](https://www.investing.com/about-us/terms-and-conditions).

For help matching an export to the format below or reproducing the analysis,
contact the corresponding author, Gustavo Andres Pava Parra, using the
contact details associated with the study. He can advise on the workflow.
Contacting him does not guarantee that he can supply a copy. Access remains
subject to the provider's terms and permissions.

## Files and dates

The script expects the 13 filenames below. The abbreviations in parentheses
come from the original filenames. Check the provider or the author's
acquisition records for the exact ticker and share class before downloading
replacement data.

| Company | Required filename |
|---|---|
| Bancolombia | `Bancolombia.csv` |
| Banco de Bogotá | `Bancodebogota.csv` |
| Bolsa de Valores de Colombia (BVC) | `Bolsavaloresdecolombia.csv` |
| Cementos Argos | `CementosArgos.csv` |
| Celsia S.A. | `CelsiaSA.csv` |
| Ecopetrol | `Ecopetrol.csv` |
| Corporación Financiera de Colombia (Corficolombiana) | `Corporacionfinancieradecolombia.csv` |
| Grupo Bolívar | `Grupobolivar.csv` |
| Mineros (MAS in supplied filename) | `Mineros SA (MAS).csv` |
| Nutresa (NCH in supplied filename) | `Nutresa (NCH).csv` |
| Organización Terpel (TPL in supplied filename) | `Organizacion Terpel SA (TPL).csv` |
| Promigas (PMG in supplied filename) | `Promigas (PMG).csv` |
| Suramericana (SIS in supplied filename) | `Suramericana (SIS).csv` |

The study covers daily trading sessions from **1 July 2021 to 30 June 2023**.
Each company can have a different set of dates and its file may end before
30 June. The script combines the dates that actually appear across all 13
files. The original inputs give 491 distinct price dates. After dropping the
first date when taking log differences, the analysis uses 490 aligned trading
sessions for returns and volume. Keep the original calendar and any gaps.

Place your authorized copies in `data/raw/` under the repository root. Keep
the filenames exactly as shown, including spaces and parentheses. The script
reads these local files directly.

```text
data/
└── raw/
    ├── Bancolombia.csv
    ├── Bancodebogota.csv
    ├── Bolsavaloresdecolombia.csv
    ├── CementosArgos.csv
    ├── CelsiaSA.csv
    ├── Ecopetrol.csv
    ├── Corporacionfinancieradecolombia.csv
    ├── Grupobolivar.csv
    ├── Mineros SA (MAS).csv
    ├── Nutresa (NCH).csv
    ├── Organizacion Terpel SA (TPL).csv
    ├── Promigas (PMG).csv
    └── Suramericana (SIS).csv
```

## CSV format

Save the files as UTF-8, with or without a byte-order mark, and keep a single
header row. Use semicolons to separate columns, commas for decimal marks and
periods for thousands separators. The script reads the date, closing price
and trading volume from columns 1, 2 and 6. Keep all seven columns in their
original order, including the Spanish headers. The importer sorts the dates
and expects one row per company per date.

| Position | Original header | Input type | Interpretation |
|---|---|---|---|
| 1 | `Fecha` | Character date | Day, month and year. The reviewed files use `dd.mm.yyyy`. The parser also accepts `/` and `-` as date separators. |
| 2 | `Último` | Locale-formatted numeric text | Closing price in the provider's native quote units |
| 3 | `Apertura` | Locale-formatted numeric text | Opening price. Keep this column in the export. |
| 4 | `Máximo` | Locale-formatted numeric text | High price. Keep this column in the export. |
| 5 | `Mínimo` | Locale-formatted numeric text | Low price. Keep this column in the export. |
| 6 | `Vol.` | Numeric text, optionally suffixed | Observed trading volume in the provider's native units |
| 7 | `% var.` or `%var.` | Percentage text | Daily percentage change. Keep this column in the export. |

Leave missing cells blank or mark them `NA`. A single hyphen also counts as a
missing volume. An absent date becomes a missing cell when the script combines
the company series. A reported zero volume stays zero.

The volume parser trims spaces and expands `K`, `M` and `B` (in either case)
by 1,000, 1,000,000 and 1,000,000,000. It reads a period as a thousands
separator even when the value has a suffix. If you use a new export, check
its decimal format before converting it to the layout above.

Use prices and volumes in the same native units as the original exports.
The CSVs do not specify the currency, lot size, volume unit or treatment of corporate-action adjustments. The script
keeps the reported units and expands the volume suffixes. Check the source
records or provider for any instrument and unit details you need to match
when preparing replacement files.

## How the script uses the inputs

1. The script reads columns 1, 2 and 6 as `date`, `Close` and `Volume`. It sorts
   each company by date and combines the dates across companies.
2. To compute returns, it fills missing closing prices forward and fills any
   gaps at the beginning backward. It then takes differences of log closing
   prices and drops the first date. Returns are expressed as log differences
   without multiplying by 100.
3. It drops 1 July 2021 from the volume calendar before imputation. To fill
   missing volume values, it fits a GAM to `log1p(volume)` using a cubic
   regression spline of calendar date and a weekday factor, estimated by REML.
   It transforms the estimates back to the original scale, sets negative
   predictions to zero and keeps the observed values. With fewer than 12
   observed values, it fills forward and backward on the log scale and sets
   entries still missing to zero. It aligns the completed volumes with the
   return dates and applies `log1p` for the analysis.
4. The script also computes the exploratory spline interpolation of returns.
   The final analysis uses the returns calculated from the filled prices in
   step 2.
5. It smooths the curves with cubic B-splines and second-derivative penalties.
   The return curves use 150 basis functions and a penalty of 100. For
   volume, the GCV search selects 24 basis functions and a penalty of
   1,000 on the study inputs.
6. It centers the curves across companies and compares three functional
   regression penalty pairs on a common 150 × 150 grid. After comparing the
   surfaces, we selected scenario B, with penalties of 1,000 and 100,000,000
   in the volume and return directions. The intercept penalty is 100,000.
   The script evaluates the final coefficient surface on a 490 × 490 grid.
7. It runs 13 leave-one-company-out (LOCO) refits. Each refit omits one
   company’s pair of curves, recenters the remaining twelve pairs and
   uses the same preprocessing, univariate smooths,
   bases and penalties. The relative L2 distances and sign agreement use
   trapezoidal weights.
8. The GARCH comparison uses standardized log volume, the original candidate
   orders, the hybrid solver and seed 123.

The analysis keeps observations, fitted models and matrices in memory. It
writes the 24 EPS figures to `output/figures/`. The optional table
export described in [README.md](README.md) writes local summaries to
`output/tables/`. Keep those outputs with your local analysis files.
