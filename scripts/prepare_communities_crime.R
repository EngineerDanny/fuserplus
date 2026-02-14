#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
raw_data <- if (length(args) >= 1) args[[1]] else "data/raw/communities.data"
raw_names <- if (length(args) >= 2) args[[2]] else "data/raw/communities.names"
out_rds <- if (length(args) >= 3) args[[3]] else "data/processed/communities_crime_l2.rds"

if (!file.exists(raw_data)) stop("Missing file: ", raw_data)
if (!file.exists(raw_names)) stop("Missing file: ", raw_names)

name_lines <- readLines(raw_names, warn = FALSE)
attr_lines <- grep("^@attribute\\s+", name_lines, value = TRUE)
attr_names <- sub("^@attribute\\s+([^ ]+)\\s+.*$", "\\1", attr_lines)

d <- read.csv(raw_data, header = FALSE, na.strings = "?", stringsAsFactors = FALSE)
if (ncol(d) != length(attr_names)) {
  stop("Column count mismatch: data has ", ncol(d), " cols but names has ", length(attr_names))
}
colnames(d) <- attr_names

required <- c("state", "communityname", "ViolentCrimesPerPop")
missing_required <- setdiff(required, colnames(d))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

# Keep rows with valid response and grouping key.
d <- d[!is.na(d$state) & !is.na(d$ViolentCrimesPerPop), , drop = FALSE]

# Build groups from state ids (46 states in this dataset).
groups_raw <- as.integer(d$state)
groups <- as.integer(factor(groups_raw, levels = sort(unique(groups_raw))))

# Remove non-predictive identifiers and target.
drop_cols <- c("state", "county", "community", "communityname", "fold", "ViolentCrimesPerPop")
feature_cols <- setdiff(colnames(d), drop_cols)
Xdf <- d[, feature_cols, drop = FALSE]

# Retain numeric predictors only.
is_num <- vapply(Xdf, is.numeric, logical(1))
Xdf <- Xdf[, is_num, drop = FALSE]

# Drop very-missing columns; median-impute the rest.
na_frac <- colMeans(is.na(Xdf))
Xdf <- Xdf[, na_frac <= 0.20, drop = FALSE]

for (j in seq_len(ncol(Xdf))) {
  col <- Xdf[[j]]
  if (anyNA(col)) {
    med <- stats::median(col, na.rm = TRUE)
    if (is.na(med)) med <- 0
    col[is.na(col)] <- med
    Xdf[[j]] <- col
  }
}

# Remove zero-variance columns after imputation.
keep_var <- vapply(Xdf, function(z) stats::sd(z) > 0, logical(1))
Xdf <- Xdf[, keep_var, drop = FALSE]

X <- as.matrix(Xdf)
storage.mode(X) <- "double"
y <- as.numeric(d$ViolentCrimesPerPop)

k <- length(unique(groups))

G_dense <- matrix(1, k, k)
G_chain <- matrix(0, k, k)
if (k > 1) {
  for (i in seq_len(k - 1)) {
    G_chain[i, i + 1] <- 1
    G_chain[i + 1, i] <- 1
  }
}

out <- list(
  X = X,
  y = y,
  groups = groups,
  group_state_id = sort(unique(groups_raw)),
  feature_names = colnames(X),
  community_name = d$communityname,
  G_dense = G_dense,
  G_sparse_chain = G_chain,
  meta = list(
    source = "UCI Communities and Crime",
    n = nrow(X),
    p = ncol(X),
    k = k,
    dropped_feature_missing_threshold = 0.20,
    dropped_identifier_columns = drop_cols
  )
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, out_rds)

cat(sprintf("Saved %s\n", out_rds))
cat(sprintf("n=%d p=%d k=%d\n", nrow(X), ncol(X), k))
