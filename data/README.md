# Real Data Assets

## UCI Communities and Crime

- Raw files:
  - `data/raw/communities.data`
  - `data/raw/communities.names`
- Processed benchmark-ready file:
  - `data/processed/communities_crime_l2.rds`

### Processed object format (`communities_crime_l2.rds`)

A named list with:

- `X`: numeric matrix (`n x p`) predictors
- `y`: numeric response vector (`ViolentCrimesPerPop`)
- `groups`: integer group labels (derived from `state`)
- `group_state_id`: mapping from internal group index to original state id
- `feature_names`: predictor column names
- `community_name`: original community name field
- `G_dense`: dense all-pairs fusion matrix
- `G_sparse_chain`: chain fusion matrix
- `meta`: metadata

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
