# Real Data Assets

## UCI Communities and Crime

- Raw files:
  - `data/raw/communities/communities.data`
  - `data/raw/communities/communities.names`
- Processed benchmark-ready file:
  - `data/processed/communities_crime/grouped_regression.csv`
  - `data/processed/communities_crime/grouped_regression_meta.json`
  - `data/processed/communities_crime/grouped_regression_g_sparse_chain_edges.csv`
  - `data/processed/communities_crime/grouped_regression_g_dense_edges.csv` (when available)

### Processed CSV format

Primary CSV (`grouped_regression.csv`) columns:

- `y`: numeric response
- `group_id`: integer group label
- `group_label`: original group value
- feature columns (`X`)

Sidecar files:

- `grouped_regression_meta.json`: dataset metadata (feature names, group levels, source, etc.)
- `grouped_regression_g_sparse_chain_edges.csv`: chain graph edge list (`from,to,weight`)
- `grouped_regression_g_dense_edges.csv`: dense graph edge list when `k` is below dense limit

### Preprocessing applied

- Dropped identifier/target columns from predictors:
  - `state`, `county`, `community`, `communityname`, `fold`, `ViolentCrimesPerPop`
- Dropped groups with fewer than 2 rows (configurable in prep script)
- Kept numeric predictors only
- Dropped predictor columns with >20% missing values
- Median-imputed remaining missing values
- Dropped zero-variance predictors

### Rebuild

```bash
Rscript scripts/prepare_communities_crime.R
```

Optional 4th argument:
- `min_group_size` (default `2`)

## Multi-Dataset Grouped Regression Build

To download and process the broader grouped-regression benchmark suite:

```bash
Rscript scripts/prepare_grouped_regression_datasets.R all
```

Default suite currently covers:

- Communities and Crime
- UCI Wine Quality
- UCI Student Performance
- UCI Beijing Multi-Site Air Quality
- UCI ElectricityLoadDiagrams20112014
- OWID CO2
- World Bank WDI (API indicators)
- NYC TLC trip sample (Open Data API)
- NOAA GHCN daily sample (selected stations)
- CountyPlus (release asset)
- SARCOS (OpenML)
- School grouped regression fallback from `nlme` (`MathAchieve` + `MathAchSchool`)

This writes processed CSV-first outputs in per-dataset folders under `data/processed/`:

- `data/processed/<dataset>/grouped_regression.csv` (main tabular data)
- `data/processed/<dataset>/grouped_regression_meta.json` (metadata)
- `data/processed/<dataset>/grouped_regression_g_sparse_chain_edges.csv` (chain graph edges)
- `data/processed/<dataset>/grouped_regression_g_dense_edges.csv` (dense graph edges, when available)

A run summary is saved to:

- `data/processed/grouped_regression_datasets_summary.csv`
