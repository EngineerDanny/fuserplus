#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(foreign)
  library(Matrix)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

set.seed(1)

raw_root <- "data/raw"
processed_root <- "data/processed"
dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)

download_if_missing <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(dest)) {
    utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
  }
  dest
}

depmap_catalog <- function(cache_path) {
  download_if_missing("https://depmap.org/portal/data_page/api/data", cache_path)
  jsonlite::fromJSON(cache_path)
}

depmap_download_url <- function(catalog, file_name, release_name = NULL) {
  if (is.null(catalog$table) || !is.data.frame(catalog$table)) {
    stop("DepMap catalog does not contain a valid 'table' section.")
  }
  tab <- catalog$table
  idx <- which(tab$fileName == file_name & !is.na(tab$downloadUrl))
  if (!is.null(release_name)) {
    idx <- idx[tab$releaseName[idx] == release_name]
  }
  if (length(idx) == 0L) {
    msg <- if (is.null(release_name)) {
      sprintf("Could not find DepMap file '%s'.", file_name)
    } else {
      sprintf("Could not find DepMap file '%s' in release '%s'.", file_name, release_name)
    }
    stop(msg)
  }
  url <- tab$downloadUrl[idx[1L]]
  if (startsWith(url, "/")) url <- paste0("https://depmap.org", url)
  url
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

build_graphs <- function(k, dense_limit = 400L) {
  k <- as.integer(k)
  dense <- NULL
  if (k <= dense_limit) dense <- matrix(1, k, k)
  chain <- Matrix(0, k, k, sparse = TRUE)
  if (k > 1L) {
    idx <- seq_len(k - 1L)
    chain[cbind(idx, idx + 1L)] <- 1
    chain[cbind(idx + 1L, idx)] <- 1
  }
  list(G_dense = dense, G_sparse_chain = chain, dense_limit = dense_limit)
}

encode_mixed_df <- function(df) {
  out <- lapply(df, function(col) {
    if (is.numeric(col)) return(col)
    if (is.logical(col)) return(as.integer(col))
    as.integer(factor(col))
  })
  as.data.frame(out, stringsAsFactors = FALSE)
}

finalize_grouped <- function(
  df, y_col, group_col, feature_cols, dataset_name, source,
  min_group_size = 5L, max_rows = 200000L, dense_limit = 400L,
  extra_meta = list()
) {
  keep <- c(y_col, group_col, feature_cols)
  keep <- unique(keep[keep %in% names(df)])
  df <- df[, keep, drop = FALSE]
  names(df)[names(df) == y_col] <- ".y"
  names(df)[names(df) == group_col] <- ".group_raw"

  df <- df[!is.na(df$.y) & !is.na(df$.group_raw), , drop = FALSE]
  group_counts <- table(df$.group_raw)
  good_groups <- names(group_counts[group_counts >= min_group_size])
  df <- df[df$.group_raw %in% good_groups, , drop = FALSE]

  feat_df <- df[, setdiff(names(df), c(".y", ".group_raw")), drop = FALSE]
  feat_df <- encode_mixed_df(feat_df)

  # Keep numeric only and median-impute.
  is_num <- vapply(feat_df, is.numeric, logical(1))
  feat_df <- feat_df[, is_num, drop = FALSE]
  if (ncol(feat_df) == 0L) stop("No numeric features available after encoding.")

  for (j in seq_len(ncol(feat_df))) {
    v <- feat_df[[j]]
    if (anyNA(v)) {
      med <- stats::median(v, na.rm = TRUE)
      if (is.na(med)) med <- 0
      v[is.na(v)] <- med
      feat_df[[j]] <- v
    }
  }

  keep_var <- vapply(feat_df, function(z) stats::sd(z) > 0, logical(1))
  feat_df <- feat_df[, keep_var, drop = FALSE]
  if (ncol(feat_df) == 0L) stop("No non-constant features available.")

  n <- nrow(feat_df)
  if (!is.null(max_rows) && n > max_rows) {
    idx <- sort(sample.int(n, max_rows))
    feat_df <- feat_df[idx, , drop = FALSE]
    df <- df[idx, , drop = FALSE]
  }

  groups <- as.integer(factor(df$.group_raw))
  y <- as.numeric(df$.y)
  X <- as.matrix(feat_df)
  storage.mode(X) <- "double"

  k <- length(unique(groups))
  graphs <- build_graphs(k, dense_limit = dense_limit)

  list(
    X = X,
    y = y,
    groups = groups,
    group_levels = sort(unique(df$.group_raw)),
    feature_names = colnames(X),
    G_dense = graphs$G_dense,
    G_sparse_chain = graphs$G_sparse_chain,
    meta = c(list(
      dataset = dataset_name,
      source = source,
      n = nrow(X),
      p = ncol(X),
      k = k,
      y_col = y_col,
      group_col = group_col,
      dense_graph_limit = graphs$dense_limit
    ), extra_meta)
  )
}

save_grouped <- function(obj, slug) {
  dataset_dir <- file.path(processed_root, slug)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  data_csv <- file.path(dataset_dir, "grouped_regression.csv")
  data_rds <- file.path(dataset_dir, "grouped_regression.rds")
  meta_json <- file.path(dataset_dir, "grouped_regression_meta.json")
  chain_edges_csv <- file.path(dataset_dir, "grouped_regression_g_sparse_chain_edges.csv")
  dense_edges_csv <- file.path(dataset_dir, "grouped_regression_g_dense_edges.csv")

  # Primary tabular artifact for user control.
  group_levels <- obj$group_levels
  if (is.null(group_levels) && !is.null(obj$group_state_id)) {
    group_levels <- obj$group_state_id
  }
  if (is.null(group_levels)) {
    group_levels <- sort(unique(obj$groups))
  }
  if (length(group_levels) >= max(obj$groups)) {
    group_label <- as.character(group_levels[obj$groups])
  } else {
    group_label <- as.character(obj$groups)
  }
  dat <- data.frame(
    y = as.numeric(obj$y),
    group_id = as.integer(obj$groups),
    group_label = group_label,
    obj$X,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(obj$feature_names) && length(obj$feature_names) == ncol(obj$X)) {
    colnames(dat)[4:ncol(dat)] <- obj$feature_names
  }
  utils::write.csv(dat, data_csv, row.names = FALSE)

  matrix_to_edge_df <- function(G) {
    if (is.null(G)) return(data.frame(from = integer(0), to = integer(0), weight = numeric(0)))
    if (inherits(G, "sparseMatrix")) {
      trip <- Matrix::summary(G)
      keep <- trip$i < trip$j & trip$x != 0
      data.frame(
        from = trip$i[keep],
        to = trip$j[keep],
        weight = as.numeric(trip$x[keep]),
        row.names = NULL
      )
    } else {
      idx <- which(G != 0 & upper.tri(G), arr.ind = TRUE)
      data.frame(
        from = idx[, 1],
        to = idx[, 2],
        weight = as.numeric(G[idx]),
        row.names = NULL
      )
    }
  }

  chain_edges <- matrix_to_edge_df(obj$G_sparse_chain)
  utils::write.csv(chain_edges, chain_edges_csv, row.names = FALSE)

  if (!is.null(obj$G_dense)) {
    dense_edges <- matrix_to_edge_df(obj$G_dense)
    utils::write.csv(dense_edges, dense_edges_csv, row.names = FALSE)
  } else {
    dense_edges_csv <- NA_character_
  }

  # Benchmark-ready binary artifact.
  saveRDS(
    list(
      X = obj$X,
      y = obj$y,
      groups = obj$groups,
      group_levels = group_levels,
      feature_names = obj$feature_names,
      G_dense = obj$G_dense,
      G_sparse_chain = obj$G_sparse_chain,
      meta = obj$meta
    ),
    file = data_rds
  )

  meta <- obj$meta
  meta$feature_names <- as.character(obj$feature_names)
  meta$group_levels <- as.character(group_levels)
  meta$output_files <- list(
    data_csv = data_csv,
    data_rds = data_rds,
    sparse_chain_edges_csv = chain_edges_csv,
    dense_edges_csv = dense_edges_csv
  )
  writeLines(
    jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    con = meta_json
  )

  data_csv
}

download_github_raw_if_missing <- function(repo_path, dest) {
  # Use github.com/raw endpoint so Git LFS objects resolve to binary content.
  url <- sprintf("https://github.com/cBioPortal/datahub/raw/master/%s", repo_path)
  download_if_missing(url, dest)
}

read_cbio_clinical_patient <- function(path) {
  # cBioPortal clinical files use four metadata lines before the tabular header.
  utils::read.delim(
    path,
    sep = "\t",
    skip = 4L,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_cbio_expression <- function(path) {
  # Expression files are large and wide; fread is materially faster/more memory-friendly.
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(
      data.table::fread(path, sep = "\t", quote = "", data.table = FALSE, check.names = FALSE),
      stringsAsFactors = FALSE
    ))
  }
  utils::read.delim(
    path,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

prepare_communities <- function() {
  # Canonical location in a dataset-specific subfolder.
  dir.create(file.path(raw_root, "communities"), recursive = TRUE, showWarnings = FALSE)
  data_path <- file.path(raw_root, "communities", "communities.data")
  names_path <- file.path(raw_root, "communities", "communities.names")
  download_if_missing(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/communities/communities.data",
    data_path
  )
  download_if_missing(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/communities/communities.names",
    names_path
  )
  tmp_rds <- tempfile(pattern = "communities_tmp_", fileext = ".rds")
  on.exit(unlink(tmp_rds), add = TRUE)
  system2("Rscript", c("scripts/prepare_communities_crime.R", data_path, names_path, tmp_rds, "2"))
  obj <- readRDS(tmp_rds)
  out <- save_grouped(obj, slug = "communities_crime")
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_wine_quality <- function() {
  slug <- "wine_quality"
  zip_path <- file.path(raw_root, slug, "wine_quality.zip")
  unzip_dir <- file.path(raw_root, slug)
  download_if_missing("https://archive.ics.uci.edu/static/public/186/wine+quality.zip", zip_path)
  utils::unzip(zip_path, exdir = unzip_dir)
  red <- utils::read.csv(file.path(unzip_dir, "winequality-red.csv"), sep = ";", stringsAsFactors = FALSE)
  white <- utils::read.csv(file.path(unzip_dir, "winequality-white.csv"), sep = ";", stringsAsFactors = FALSE)
  red$wine_type <- "red"
  white$wine_type <- "white"
  d <- rbind(red, white)

  obj <- finalize_grouped(
    d, y_col = "quality", group_col = "wine_type",
    feature_cols = setdiff(names(d), c("quality", "wine_type")),
    dataset_name = "UCI Wine Quality",
    source = "https://archive.ics.uci.edu/dataset/186/wine+quality"
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_student_performance <- function() {
  slug <- "student_performance"
  zip_path <- file.path(raw_root, slug, "student_performance.zip")
  unzip_dir <- file.path(raw_root, slug)
  download_if_missing("https://archive.ics.uci.edu/static/public/320/student+performance.zip", zip_path)
  utils::unzip(zip_path, exdir = unzip_dir)
  nested_zip <- file.path(unzip_dir, "student.zip")
  if (file.exists(nested_zip)) utils::unzip(nested_zip, exdir = unzip_dir)
  mat <- utils::read.csv(file.path(unzip_dir, "student-mat.csv"), sep = ";", stringsAsFactors = FALSE)
  por <- utils::read.csv(file.path(unzip_dir, "student-por.csv"), sep = ";", stringsAsFactors = FALSE)
  mat$subject <- "math"
  por$subject <- "portuguese"
  d <- rbind(mat, por)

  obj <- finalize_grouped(
    d, y_col = "G3", group_col = "subject",
    feature_cols = setdiff(names(d), c("G3", "subject")),
    dataset_name = "UCI Student Performance",
    source = "https://archive.ics.uci.edu/dataset/320/student+performance",
    min_group_size = 2L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_beijing_air <- function() {
  slug <- "beijing_air_quality"
  zip_path <- file.path(raw_root, slug, "beijing_air_quality.zip")
  unzip_dir <- file.path(raw_root, slug)
  download_if_missing(
    "https://archive.ics.uci.edu/ml/machine-learning-databases/00501/PRSA2017_Data_20130301-20170228.zip",
    zip_path
  )
  utils::unzip(zip_path, exdir = unzip_dir)
  nested_zip <- file.path(unzip_dir, "PRSA2017_Data_20130301-20170228.zip")
  if (file.exists(nested_zip)) utils::unzip(nested_zip, exdir = unzip_dir)
  csv_files <- list.files(unzip_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  csv_files <- csv_files[grepl("^PRSA_Data_.*\\.csv$", basename(csv_files))]
  if (length(csv_files) == 0L) stop("No CSV files found in Beijing zip.")
  dlist <- lapply(csv_files, function(f) utils::read.csv(f, stringsAsFactors = FALSE))
  d <- do.call(rbind, dlist)

  # Robust naming across files.
  nms <- names(d)
  nms_low <- tolower(nms)
  y_col <- if ("PM2.5" %in% nms) "PM2.5" else if ("pm2.5" %in% nms_low) nms[which(nms_low == "pm2.5")[1]] else stop("No PM2.5 column.")
  group_col <- if ("station" %in% nms_low) nms[which(nms_low == "station")[1]] else stop("No station column.")

  feature_candidates <- c("PM10", "SO2", "NO2", "CO", "O3", "TEMP", "PRES", "DEWP", "RAIN", "WSPM", "year", "month", "day", "hour")
  feature_cols <- intersect(feature_candidates, nms)

  obj <- finalize_grouped(
    d, y_col = y_col, group_col = group_col, feature_cols = feature_cols,
    dataset_name = "UCI Beijing Multi-Site Air Quality",
    source = "https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data",
    max_rows = 180000L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_electricity <- function() {
  slug <- "electricity_load"
  zip_path <- file.path(raw_root, slug, "electricity_load.zip")
  unzip_dir <- file.path(raw_root, slug)
  download_if_missing(
    "https://archive.ics.uci.edu/static/public/321/electricityloaddiagrams20112014.zip",
    zip_path
  )
  utils::unzip(zip_path, exdir = unzip_dir)
  txt_path <- file.path(unzip_dir, "LD2011_2014.txt")
  d <- utils::read.csv(txt_path, sep = ";", dec = ",", check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) < 3L) stop("Unexpected ElectricityLoadDiagrams format.")

  n_timestamps_target <- 2500L
  n_clients_target <- 120L
  times_all <- as.POSIXct(d[[1]], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (all(is.na(times_all))) stop("Could not parse timestamp in ElectricityLoadDiagrams.")
  client_cols_all <- names(d)[-1L]

  # Pick active clients by nonzero usage and variability to avoid all-zero slices.
  client_nonzero <- vapply(client_cols_all, function(cl) {
    z <- suppressWarnings(as.numeric(d[[cl]]))
    sum(is.finite(z) & z != 0)
  }, numeric(1))
  client_sd <- vapply(client_cols_all, function(cl) {
    z <- suppressWarnings(as.numeric(d[[cl]]))
    stats::sd(z[is.finite(z)])
  }, numeric(1))
  active_clients <- client_cols_all[client_nonzero > 0 & is.finite(client_sd) & client_sd > 0]
  if (length(active_clients) == 0L) stop("No active electricity clients with nonzero and varying load found.")
  active_clients <- active_clients[order(client_nonzero[active_clients], decreasing = TRUE)]
  client_cols <- head(active_clients, min(n_clients_target, length(active_clients)))

  client_mat <- do.call(cbind, lapply(client_cols, function(cl) suppressWarnings(as.numeric(d[[cl]]))))
  if (!is.matrix(client_mat) || nrow(client_mat) == 0L) stop("Failed to build electricity client matrix.")

  # Pick a contiguous time window with the most active client measurements.
  row_activity <- rowSums(is.finite(client_mat) & client_mat != 0)
  n_timestamps <- min(n_timestamps_target, nrow(client_mat))
  if (n_timestamps < nrow(client_mat)) {
    csum <- c(0, cumsum(row_activity))
    window_sums <- csum[(n_timestamps + 1L):(nrow(client_mat) + 1L)] - csum[1L:(nrow(client_mat) - n_timestamps + 1L)]
    start_idx <- which.max(window_sums)
    row_idx <- seq.int(start_idx, start_idx + n_timestamps - 1L)
  } else {
    row_idx <- seq_len(n_timestamps)
  }
  d <- d[row_idx, c(1L, match(client_cols, names(d))), drop = FALSE]
  times <- times_all[row_idx]

  long <- do.call(rbind, lapply(client_cols, function(cl) {
    data.frame(group_raw = cl, time = times, load = as.numeric(d[[cl]]), stringsAsFactors = FALSE)
  }))
  long <- long[is.finite(long$load) & !is.na(long$time), , drop = FALSE]
  if (nrow(long) == 0L || !is.finite(stats::sd(long$load)) || stats::sd(long$load) <= 0) {
    stop("Electricity load is constant after active client/window selection; cannot build regression dataset.")
  }
  long <- long[order(long$group_raw, long$time), , drop = FALSE]
  long$lag1 <- ave(long$load, long$group_raw, FUN = function(z) c(NA_real_, head(z, -1L)))
  long$hour <- as.integer(format(long$time, "%H"))
  long$wday <- as.integer(format(long$time, "%u"))
  long$month <- as.integer(format(long$time, "%m"))
  long <- long[!is.na(long$lag1), , drop = FALSE]
  if (nrow(long) == 0L || !is.finite(stats::sd(long$load)) || stats::sd(long$load) <= 0) {
    stop("Electricity load became constant after lag feature construction; cannot build regression dataset.")
  }

  obj <- finalize_grouped(
    long, y_col = "load", group_col = "group_raw",
    feature_cols = c("lag1", "hour", "wday", "month"),
    dataset_name = "UCI ElectricityLoadDiagrams20112014",
    source = "https://archive.ics.uci.edu/dataset/321/electricityloaddiagrams20112014",
    max_rows = 200000L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_owid_co2 <- function() {
  slug <- "owid_co2"
  csv_path <- file.path(raw_root, slug, "owid-co2-data.csv")
  download_if_missing("https://raw.githubusercontent.com/owid/co2-data/master/owid-co2-data.csv", csv_path)
  d <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  keep <- c("country", "iso_code", "year", "co2_per_capita", "gdp", "population", "primary_energy_consumption")
  d <- d[, intersect(keep, names(d)), drop = FALSE]
  d <- d[nchar(d$iso_code) == 3, , drop = FALSE]
  d <- d[complete.cases(d[, c("co2_per_capita", "gdp", "population", "primary_energy_consumption")]), , drop = FALSE]

  obj <- finalize_grouped(
    d, y_col = "co2_per_capita", group_col = "country",
    feature_cols = c("year", "gdp", "population", "primary_energy_consumption"),
    dataset_name = "OWID CO2",
    source = "https://github.com/owid/co2-data",
    max_rows = 150000L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

fetch_wb_indicator <- function(indicator) {
  url <- sprintf(
    "https://api.worldbank.org/v2/country/all/indicator/%s?format=json&per_page=20000",
    utils::URLencode(indicator, reserved = TRUE)
  )
  js <- NULL
  last_err <- NULL
  for (attempt in 1:3) {
    out <- try(jsonlite::fromJSON(url), silent = TRUE)
    if (!inherits(out, "try-error")) {
      js <- out
      break
    }
    last_err <- out
    Sys.sleep(1.5 * attempt)
  }
  if (is.null(js)) {
    stop("World Bank API failed for ", indicator, ": ", as.character(last_err))
  }
  if (length(js) < 2L || is.null(js[[2]])) stop("World Bank API returned no data for ", indicator)
  d <- js[[2]]
  d <- d[, c("countryiso3code", "date", "value")]
  names(d) <- c("iso3c", "year", indicator)
  d$year <- as.integer(d$year)
  d
}

prepare_world_bank_wdi <- function() {
  slug <- "world_bank_wdi"
  indicators <- c("NY.GDP.PCAP.CD", "SP.POP.TOTL", "NE.EXP.GNFS.CD", "NE.IMP.GNFS.CD")
  dlist <- lapply(indicators, fetch_wb_indicator)
  d <- Reduce(function(a, b) merge(a, b, by = c("iso3c", "year"), all = FALSE), dlist)
  d <- d[nchar(d$iso3c) == 3, , drop = FALSE]
  d <- d[complete.cases(d[, indicators]), , drop = FALSE]

  obj <- finalize_grouped(
    d, y_col = "NY.GDP.PCAP.CD", group_col = "iso3c",
    feature_cols = c("year", "SP.POP.TOTL", "NE.EXP.GNFS.CD", "NE.IMP.GNFS.CD"),
    dataset_name = "World Bank WDI",
    source = "https://data.worldbank.org",
    max_rows = 200000L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_nyc_tlc <- function() {
  slug <- "nyc_tlc"
  csv_path <- file.path(raw_root, slug, "nyc_tlc_sample.csv")
  url <- "https://data.cityofnewyork.us/resource/pqfs-mqru.csv?$limit=120000"
  download_if_missing(url, csv_path)
  d <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  names(d) <- tolower(names(d))

  # Pick a usable group column.
  group_candidates <- c("pulocationid", "vendorid", "ratecodeid")
  group_col <- NULL
  for (gc in group_candidates) {
    if (gc %in% names(d)) {
      u <- unique(d[[gc]])
      u <- u[!is.na(u) & u != ""]
      if (length(u) >= 2L) {
        group_col <- gc
        break
      }
    }
  }
  if (is.null(group_col)) stop("No usable group column found in TLC sample.")

  y_col <- "total_amount"
  feature_candidates <- c("trip_distance", "fare_amount", "tip_amount", "tolls_amount", "passenger_count")
  feature_cols <- intersect(feature_candidates, names(d))
  for (nm in c(y_col, feature_cols)) d[[nm]] <- suppressWarnings(as.numeric(d[[nm]]))
  d <- d[complete.cases(d[, c(y_col, feature_cols, group_col)]), , drop = FALSE]

  obj <- finalize_grouped(
    d, y_col = y_col, group_col = group_col, feature_cols = feature_cols,
    dataset_name = "NYC TLC (sample)",
    source = "https://data.cityofnewyork.us",
    min_group_size = 20L,
    max_rows = 120000L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_noaa_ghcn <- function() {
  slug <- "noaa_ghcn"
  station_ids <- c("USW00094728", "USC00042319", "USW00023174")
  station_dir <- file.path(raw_root, slug)
  dir.create(station_dir, recursive = TRUE, showWarnings = FALSE)

  station_frames <- list()
  for (sid in station_ids) {
    gz_path <- file.path(station_dir, paste0(sid, ".csv.gz"))
    u <- sprintf("https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/%s.csv.gz", sid)
    ok <- TRUE
    tryCatch(download_if_missing(u, gz_path), error = function(e) ok <<- FALSE)
    if (!ok) next
    d <- utils::read.csv(gz_path, header = FALSE, stringsAsFactors = FALSE)
    names(d) <- c("station", "date", "element", "data_value", "m_flag", "q_flag", "s_flag", "obs_time")
    d <- d[d$element %in% c("TMAX", "TMIN", "PRCP", "SNOW"), c("date", "element", "data_value"), drop = FALSE]
    if (nrow(d) == 0L) next
    w <- reshape(d, idvar = "date", timevar = "element", direction = "wide")
    names(w) <- sub("^data_value\\.", "", names(w))
    w$station <- sid
    w$date <- as.Date(as.character(w$date), format = "%Y%m%d")
    w$year <- as.integer(format(w$date, "%Y"))
    w$doy <- as.integer(format(w$date, "%j"))
    # Unit normalization: TMAX/TMIN in tenths C, PRCP/SNOW in tenths mm.
    if ("TMAX" %in% names(w)) w$TMAX <- as.numeric(w$TMAX) / 10
    if ("TMIN" %in% names(w)) w$TMIN <- as.numeric(w$TMIN) / 10
    if ("PRCP" %in% names(w)) w$PRCP <- as.numeric(w$PRCP) / 10
    if ("SNOW" %in% names(w)) w$SNOW <- as.numeric(w$SNOW) / 10
    station_frames[[sid]] <- w
  }
  if (length(station_frames) == 0L) stop("No GHCN stations downloaded successfully.")
  d <- do.call(rbind, station_frames)
  d <- d[complete.cases(d[, c("TMAX", "TMIN")]), , drop = FALSE]
  d$PRCP[is.na(d$PRCP)] <- 0
  d$SNOW[is.na(d$SNOW)] <- 0

  obj <- finalize_grouped(
    d, y_col = "TMAX", group_col = "station",
    feature_cols = c("TMIN", "PRCP", "SNOW", "year", "doy"),
    dataset_name = "NOAA GHCN Daily (sampled stations)",
    source = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/",
    max_rows = 200000L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

choose_countyplus_target <- function(d, exclude_cols) {
  num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, exclude_cols)
  if (length(num_cols) == 0L) stop("No numeric columns available for CountyPlus target.")

  preferred <- c("hpi", "hpi_at_bdl", "unemployment_rate", "poverty_rate", "family_income", "median_income", "consumption")
  for (nm in preferred) {
    if (nm %in% num_cols) return(nm)
  }
  # Fallback: choose high-variance numeric column.
  sds <- vapply(d[num_cols], function(x) stats::sd(x, na.rm = TRUE), numeric(1))
  num_cols[which.max(sds)]
}

prepare_countyplus <- function() {
  slug <- "countyplus"
  dta_path <- file.path(raw_root, slug, "CountyPlus.dta")
  download_if_missing(
    "https://github.com/Clpr/CountyPlus/releases/download/v0.0.2/CountyPlus.dta",
    dta_path
  )
  d <- tryCatch(
    foreign::read.dta(dta_path, convert.factors = FALSE),
    error = function(e) {
      if (requireNamespace("readstata13", quietly = TRUE)) {
        as.data.frame(readstata13::read.dta13(dta_path))
      } else {
        stop("CountyPlus requires package 'readstata13' for modern .dta files. Install it and rerun.")
      }
    }
  )

  names_lower <- tolower(names(d))
  names(d) <- names_lower

  group_col <- intersect(c("fips", "fips_county", "county_fips", "countyfips"), names(d))
  if (length(group_col) == 0L) {
    id_like <- grep("fips", names(d), value = TRUE)
    if (length(id_like) > 0L) group_col <- id_like[1]
  } else {
    group_col <- group_col[1]
  }
  if (length(group_col) == 0L) stop("Could not identify CountyPlus group id (FIPS-like) column.")

  year_col <- if ("year" %in% names(d)) "year" else NULL
  target_col <- choose_countyplus_target(d, exclude_cols = c(group_col, year_col))

  # Numeric features with low missingness.
  num_cols <- names(d)[vapply(d, is.numeric, logical(1))]
  feat <- setdiff(num_cols, c(group_col, year_col, target_col))
  if (!is.null(year_col)) feat <- c(year_col, feat)
  if (length(feat) > 80L) {
    na_rate <- vapply(d[feat], function(x) mean(is.na(x)), numeric(1))
    feat <- names(sort(na_rate))[seq_len(80L)]
  }

  obj <- finalize_grouped(
    d, y_col = target_col, group_col = group_col, feature_cols = feat,
    dataset_name = "CountyPlus",
    source = "https://github.com/Clpr/CountyPlus/releases/tag/v0.0.2",
    min_group_size = 3L,
    max_rows = 220000L,
    dense_limit = 300L,
    extra_meta = list(target_selected = target_col)
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)), target = target_col)
}

prepare_sarcos_openml <- function() {
  slug <- "sarcos_openml"
  arff_path <- file.path(raw_root, slug, "sarcos.arff")
  download_if_missing("https://openml.org/data/v1/download/22102750/sarcos.arff", arff_path)
  d <- foreign::read.arff(arff_path)

  # OpenML default target is V22; V23-V28 are additional torques.
  target <- "V22"
  drop_extra_targets <- intersect(c("V23", "V24", "V25", "V26", "V27", "V28"), names(d))
  feature_cols <- setdiff(names(d), c(target, drop_extra_targets))

  # Synthetic grouping by decile of first joint position (V1) for grouped benchmarking.
  d$group_bin <- cut(d$V1, breaks = stats::quantile(d$V1, probs = seq(0, 1, by = 0.1), na.rm = TRUE), include.lowest = TRUE)

  obj <- finalize_grouped(
    d, y_col = target, group_col = "group_bin", feature_cols = feature_cols,
    dataset_name = "SARCOS (OpenML 43873)",
    source = "https://www.openml.org/d/43873",
    min_group_size = 50L
  )
  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_school_ilea <- function() {
  # Practical fallback: widely used school grouped-regression data shipped with nlme.
  data("MathAchieve", package = "nlme")
  data("MathAchSchool", package = "nlme")
  d <- merge(MathAchieve, MathAchSchool, by = "School", all.x = TRUE)

  obj <- finalize_grouped(
    d, y_col = "MathAch", group_col = "School",
    feature_cols = setdiff(names(d), c("MathAch", "School")),
    dataset_name = "School (nlme MathAchieve/MathAchSchool fallback)",
    source = "R package nlme datasets MathAchieve + MathAchSchool",
    min_group_size = 5L
  )
  out <- save_grouped(obj, "school_ilea")
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_nir_corn <- function() {
  slug <- "nir_corn_high_dim"
  mat_path <- file.path(raw_root, slug, "corn.mat")
  download_if_missing("https://www.eigenvector.com/data/Corn/corn.mat", mat_path)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("NIR Corn requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("m5spec", "mp5spec", "mp6spec", "propvals")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("NIR Corn .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_m5 <- as.matrix(m$m5spec$data)
  x_mp5 <- as.matrix(m$mp5spec$data)
  x_mp6 <- as.matrix(m$mp6spec$data)
  y_all <- as.matrix(m$propvals$data)
  if (ncol(y_all) < 1L) stop("NIR Corn propvals has no target columns.")

  # Default target is moisture (first property column in this data release).
  y_base <- as.numeric(y_all[, 1L])
  if (length(y_base) != nrow(x_m5) || nrow(x_m5) != nrow(x_mp5) || nrow(x_m5) != nrow(x_mp6)) {
    stop("NIR Corn dimensions are inconsistent across instruments/targets.")
  }

  wavelengths <- tryCatch({
    wl <- as.numeric(m$m5spec$axisscale[[2]][[1]][1, ])
    if (length(wl) != ncol(x_m5)) stop("bad_wl_len")
    wl
  }, error = function(e) {
    seq_len(ncol(x_m5))
  })

  X <- rbind(x_m5, x_mp5, x_mp6)
  y <- rep(y_base, times = 3L)
  grp <- c(rep("m5", nrow(x_m5)), rep("mp5", nrow(x_mp5)), rep("mp6", nrow(x_mp6)))

  wl_chr <- trimws(format(wavelengths, scientific = FALSE, trim = TRUE))
  feat_names <- paste0("wl_", gsub("[^0-9A-Za-z]+", "_", wl_chr))
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "NIR Corn (Eigenvector)",
    source = "https://www.eigenvector.com/data/Corn/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "moisture",
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_nir_tablets_high_dim <- function() {
  slug <- "nir_tablets_high_dim"
  mat_path <- file.path(raw_root, slug, "nir_shootout_2002.mat")
  download_if_missing("https://www.eigenvector.com/data/tablets/nir_shootout_2002.mat", mat_path)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("NIR Tablets requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c(
    "calibrate.1", "calibrate.2", "validate.1", "validate.2",
    "calibrate.Y", "validate.Y"
  )
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("NIR Tablets .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  get_data <- function(field_name) {
    x <- as.matrix(m[[field_name]][["data"]])
    if (is.null(dim(x)) || nrow(x) == 0L || ncol(x) == 0L) {
      stop("NIR Tablets field has empty data matrix: ", field_name)
    }
    x
  }

  x_cal_1 <- get_data("calibrate.1")
  x_cal_2 <- get_data("calibrate.2")
  x_val_1 <- get_data("validate.1")
  x_val_2 <- get_data("validate.2")
  y_cal <- get_data("calibrate.Y")
  y_val <- get_data("validate.Y")

  if (ncol(x_cal_1) != ncol(x_cal_2) || ncol(x_cal_1) != ncol(x_val_1) || ncol(x_cal_1) != ncol(x_val_2)) {
    stop("NIR Tablets spectra dimensions are inconsistent across splits/instruments.")
  }
  if (nrow(x_cal_1) != nrow(y_cal) || nrow(x_val_1) != nrow(y_val)) {
    stop("NIR Tablets spectra/response row counts are inconsistent.")
  }
  if (ncol(y_cal) < 3L || ncol(y_val) < 3L) {
    stop("NIR Tablets response matrix has fewer than 3 targets.")
  }

  # Use assay target (3rd response column), and calibration+validation only so p > n.
  target_idx <- 3L
  target_name <- "assay"
  try({
    y_labels <- trimws(as.character(m[["calibrate.Y"]][["label"]][[2]][[1]][, 1]))
    if (length(y_labels) >= target_idx && nzchar(y_labels[target_idx])) {
      target_name <- y_labels[target_idx]
    }
  }, silent = TRUE)

  x_1 <- rbind(x_cal_1, x_val_1)
  x_2 <- rbind(x_cal_2, x_val_2)
  y_1 <- as.numeric(c(y_cal[, target_idx], y_val[, target_idx]))
  y_2 <- as.numeric(c(y_cal[, target_idx], y_val[, target_idx]))

  wavelengths <- tryCatch({
    wl <- as.numeric(m[["calibrate.1"]][["axisscale"]][[2]][[1]][1, ])
    if (length(wl) != ncol(x_1)) stop("bad_wl_len")
    wl
  }, error = function(e) {
    seq_len(ncol(x_1))
  })

  wl_chr <- trimws(format(wavelengths, scientific = FALSE, trim = TRUE))
  feat_names <- paste0("wl_", gsub("[^0-9A-Za-z]+", "_", wl_chr))
  feat_names <- make.unique(feat_names, sep = "_dup")

  X <- rbind(x_1, x_2)
  y <- c(y_1, y_2)
  grp <- c(rep("spec_1", nrow(x_1)), rep("spec_2", nrow(x_2)))
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "NIR Tablets Shootout 2002 (Eigenvector)",
    source = "https://www.eigenvector.com/data/tablets/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = target_name,
      split_used = "calibrate+validate",
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_nir_tablets_hardness_high_dim <- function() {
  slug <- "nir_tablets_hardness_high_dim"
  mat_path <- file.path(raw_root, slug, "nir_shootout_2002.mat")
  download_if_missing("https://www.eigenvector.com/data/tablets/nir_shootout_2002.mat", mat_path)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("NIR Tablets requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c(
    "calibrate.1", "calibrate.2", "validate.1", "validate.2",
    "calibrate.Y", "validate.Y"
  )
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("NIR Tablets .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  get_data <- function(field_name) {
    x <- as.matrix(m[[field_name]][["data"]])
    if (is.null(dim(x)) || nrow(x) == 0L || ncol(x) == 0L) {
      stop("NIR Tablets field has empty data matrix: ", field_name)
    }
    x
  }

  x_cal_1 <- get_data("calibrate.1")
  x_cal_2 <- get_data("calibrate.2")
  x_val_1 <- get_data("validate.1")
  x_val_2 <- get_data("validate.2")
  y_cal <- get_data("calibrate.Y")
  y_val <- get_data("validate.Y")

  if (ncol(x_cal_1) != ncol(x_cal_2) || ncol(x_cal_1) != ncol(x_val_1) || ncol(x_cal_1) != ncol(x_val_2)) {
    stop("NIR Tablets spectra dimensions are inconsistent across splits/instruments.")
  }
  if (nrow(x_cal_1) != nrow(y_cal) || nrow(x_val_1) != nrow(y_val)) {
    stop("NIR Tablets spectra/response row counts are inconsistent.")
  }
  if (ncol(y_cal) < 2L || ncol(y_val) < 2L) {
    stop("NIR Tablets response matrix has fewer than 2 targets.")
  }

  # Use hardness target (2nd response column), and calibration+validation only so p > n.
  target_idx <- 2L
  target_name <- "hardness"
  try({
    y_labels <- trimws(as.character(m[["calibrate.Y"]][["label"]][[2]][[1]][, 1]))
    if (length(y_labels) >= target_idx && nzchar(y_labels[target_idx])) {
      target_name <- y_labels[target_idx]
    }
  }, silent = TRUE)

  x_1 <- rbind(x_cal_1, x_val_1)
  x_2 <- rbind(x_cal_2, x_val_2)
  y_1 <- as.numeric(c(y_cal[, target_idx], y_val[, target_idx]))
  y_2 <- as.numeric(c(y_cal[, target_idx], y_val[, target_idx]))

  wavelengths <- tryCatch({
    wl <- as.numeric(m[["calibrate.1"]][["axisscale"]][[2]][[1]][1, ])
    if (length(wl) != ncol(x_1)) stop("bad_wl_len")
    wl
  }, error = function(e) {
    seq_len(ncol(x_1))
  })

  wl_chr <- trimws(format(wavelengths, scientific = FALSE, trim = TRUE))
  feat_names <- paste0("wl_", gsub("[^0-9A-Za-z]+", "_", wl_chr))
  feat_names <- make.unique(feat_names, sep = "_dup")

  X <- rbind(x_1, x_2)
  y <- c(y_1, y_2)
  grp <- c(rep("spec_1", nrow(x_1)), rep("spec_2", nrow(x_2)))
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "NIR Tablets Shootout 2002 (Eigenvector, hardness target)",
    source = "https://www.eigenvector.com/data/tablets/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = target_name,
      split_used = "calibrate+validate",
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_swri_cn_high_dim <- function() {
  slug <- "swri_cn_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(raw_dir, "CNGATEST.ZIP")
  mat_path <- file.path(raw_dir, "CNGATEST.mat")
  download_if_missing("https://www.eigenvector.com/data/SWRI/CNGATEST.ZIP", zip_path)
  if (!file.exists(mat_path)) utils::unzip(zip_path, exdir = raw_dir)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("SWRI CN requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("cn.sd.hl", "cn.y.hl", "cn.sd.ll.a", "cn.sd.ll.b", "cn.y.ll.a", "cn.y.ll.b")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("SWRI CN .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_hl <- as.matrix(m$cn.sd.hl)
  y_hl <- as.numeric(m$cn.y.hl)
  x_ll_a <- as.matrix(m$cn.sd.ll.a)
  y_ll_a <- as.numeric(m$cn.y.ll.a)
  x_ll_b <- as.matrix(m$cn.sd.ll.b)
  y_ll_b <- as.numeric(m$cn.y.ll.b)

  if (ncol(x_hl) != ncol(x_ll_a) || ncol(x_hl) != ncol(x_ll_b)) {
    stop("SWRI CN spectra dimensions are inconsistent across subsets.")
  }
  if (nrow(x_hl) != length(y_hl) || nrow(x_ll_a) != length(y_ll_a) || nrow(x_ll_b) != length(y_ll_b)) {
    stop("SWRI CN spectra/response row counts are inconsistent.")
  }

  X <- rbind(x_ll_a, x_ll_b, x_hl)
  y <- c(y_ll_a, y_ll_b, y_hl)
  grp <- c(rep("ll_a", nrow(x_ll_a)), rep("ll_b", nrow(x_ll_b)), rep("hl", nrow(x_hl)))

  wl <- seq(750, by = 2, length.out = ncol(X))
  feat_names <- paste0("wl_", wl)
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "SWRI Diesel NIR (Cetane Number, CNGATEST)",
    source = "https://www.eigenvector.com/data/SWRI/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "CN",
      subset_labels = c("ll_a", "ll_b", "hl"),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_swri_bp50_high_dim <- function() {
  slug <- "swri_bp50_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(raw_dir, "bp50gatest.zip")
  mat_path <- file.path(raw_dir, "BP50GATEST.mat")
  download_if_missing("https://www.eigenvector.com/data/SWRI/bp50gatest.zip", zip_path)
  if (!file.exists(mat_path)) utils::unzip(zip_path, exdir = raw_dir)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("SWRI BP50 requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("bp50.s1d.hl", "bp50.y1.hl", "bp50.s1d.ll.a", "bp50.s1d.ll.b", "bp50.y1.ll.a", "bp50.y1.ll.b")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("SWRI BP50 .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_hl <- as.matrix(m$bp50.s1d.hl)
  y_hl <- as.numeric(m$bp50.y1.hl)
  x_ll_a <- as.matrix(m$bp50.s1d.ll.a)
  y_ll_a <- as.numeric(m$bp50.y1.ll.a)
  x_ll_b <- as.matrix(m$bp50.s1d.ll.b)
  y_ll_b <- as.numeric(m$bp50.y1.ll.b)

  if (ncol(x_hl) != ncol(x_ll_a) || ncol(x_hl) != ncol(x_ll_b)) {
    stop("SWRI BP50 spectra dimensions are inconsistent across subsets.")
  }
  if (nrow(x_hl) != length(y_hl) || nrow(x_ll_a) != length(y_ll_a) || nrow(x_ll_b) != length(y_ll_b)) {
    stop("SWRI BP50 spectra/response row counts are inconsistent.")
  }

  X <- rbind(x_ll_a, x_ll_b, x_hl)
  y <- c(y_ll_a, y_ll_b, y_hl)
  grp <- c(rep("ll_a", nrow(x_ll_a)), rep("ll_b", nrow(x_ll_b)), rep("hl", nrow(x_hl)))

  wl <- seq(750, by = 2, length.out = ncol(X))
  feat_names <- paste0("wl_", wl)
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "SWRI Diesel NIR (BP50, BP50GATEST)",
    source = "https://www.eigenvector.com/data/SWRI/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "BP50",
      subset_labels = c("ll_a", "ll_b", "hl"),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_swri_d4052_high_dim <- function() {
  slug <- "swri_d4052_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(raw_dir, "d4052gatest.zip")
  mat_path <- file.path(raw_dir, "D4052GATEST.mat")
  download_if_missing("https://www.eigenvector.com/data/SWRI/d4052gatest.zip", zip_path)
  if (!file.exists(mat_path)) utils::unzip(zip_path, exdir = raw_dir)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("SWRI D4052 requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("d.sd.hl", "d.y.hl", "d.sd.ll.a", "d.sd.ll.b", "d.y.ll.a", "d.y.ll.b")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("SWRI D4052 .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_hl <- as.matrix(m$d.sd.hl)
  y_hl <- as.numeric(m$d.y.hl)
  x_ll_a <- as.matrix(m$d.sd.ll.a)
  y_ll_a <- as.numeric(m$d.y.ll.a)
  x_ll_b <- as.matrix(m$d.sd.ll.b)
  y_ll_b <- as.numeric(m$d.y.ll.b)

  if (ncol(x_hl) != ncol(x_ll_a) || ncol(x_hl) != ncol(x_ll_b)) {
    stop("SWRI D4052 spectra dimensions are inconsistent across subsets.")
  }
  if (nrow(x_hl) != length(y_hl) || nrow(x_ll_a) != length(y_ll_a) || nrow(x_ll_b) != length(y_ll_b)) {
    stop("SWRI D4052 spectra/response row counts are inconsistent.")
  }

  X <- rbind(x_ll_a, x_ll_b, x_hl)
  y <- c(y_ll_a, y_ll_b, y_hl)
  grp <- c(rep("ll_a", nrow(x_ll_a)), rep("ll_b", nrow(x_ll_b)), rep("hl", nrow(x_hl)))

  wl <- seq(750, by = 2, length.out = ncol(X))
  feat_names <- paste0("wl_", wl)
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "SWRI Diesel NIR (D4052, D4052GATEST)",
    source = "https://www.eigenvector.com/data/SWRI/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "D4052",
      subset_labels = c("ll_a", "ll_b", "hl"),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_swri_freeze_high_dim <- function() {
  slug <- "swri_freeze_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(raw_dir, "freezegatest.zip")
  mat_path <- file.path(raw_dir, "FREEZEGATEST.mat")
  download_if_missing("https://www.eigenvector.com/data/SWRI/freezegatest.zip", zip_path)
  if (!file.exists(mat_path)) utils::unzip(zip_path, exdir = raw_dir)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("SWRI FREEZE requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("f.sd.hl", "f.y.hl", "f.sd.ll.a", "f.sd.ll.b", "f.y.ll.a", "f.y.ll.b")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("SWRI FREEZE .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_hl <- as.matrix(m$f.sd.hl)
  y_hl <- as.numeric(m$f.y.hl)
  x_ll_a <- as.matrix(m$f.sd.ll.a)
  y_ll_a <- as.numeric(m$f.y.ll.a)
  x_ll_b <- as.matrix(m$f.sd.ll.b)
  y_ll_b <- as.numeric(m$f.y.ll.b)

  if (ncol(x_hl) != ncol(x_ll_a) || ncol(x_hl) != ncol(x_ll_b)) {
    stop("SWRI FREEZE spectra dimensions are inconsistent across subsets.")
  }
  if (nrow(x_hl) != length(y_hl) || nrow(x_ll_a) != length(y_ll_a) || nrow(x_ll_b) != length(y_ll_b)) {
    stop("SWRI FREEZE spectra/response row counts are inconsistent.")
  }

  X <- rbind(x_ll_a, x_ll_b, x_hl)
  y <- c(y_ll_a, y_ll_b, y_hl)
  grp <- c(rep("ll_a", nrow(x_ll_a)), rep("ll_b", nrow(x_ll_b)), rep("hl", nrow(x_hl)))

  wl <- seq(750, by = 2, length.out = ncol(X))
  feat_names <- paste0("wl_", wl)
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "SWRI Diesel NIR (FREEZE, FREEZEGATEST)",
    source = "https://www.eigenvector.com/data/SWRI/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "FREEZE",
      subset_labels = c("ll_a", "ll_b", "hl"),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_swri_total_high_dim <- function() {
  slug <- "swri_total_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(raw_dir, "totalgatest.zip")
  mat_path <- file.path(raw_dir, "TOTALGATEST.mat")
  download_if_missing("https://www.eigenvector.com/data/SWRI/totalgatest.zip", zip_path)
  if (!file.exists(mat_path)) utils::unzip(zip_path, exdir = raw_dir)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("SWRI TOTAL requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("t.sd.hl", "t.y.hl", "t.sd.ll.a", "t.sd.ll.b", "t.y.ll.a", "t.y.ll.b")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("SWRI TOTAL .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_hl <- as.matrix(m$t.sd.hl)
  y_hl <- as.numeric(m$t.y.hl)
  x_ll_a <- as.matrix(m$t.sd.ll.a)
  y_ll_a <- as.numeric(m$t.y.ll.a)
  x_ll_b <- as.matrix(m$t.sd.ll.b)
  y_ll_b <- as.numeric(m$t.y.ll.b)

  if (ncol(x_hl) != ncol(x_ll_a) || ncol(x_hl) != ncol(x_ll_b)) {
    stop("SWRI TOTAL spectra dimensions are inconsistent across subsets.")
  }
  if (nrow(x_hl) != length(y_hl) || nrow(x_ll_a) != length(y_ll_a) || nrow(x_ll_b) != length(y_ll_b)) {
    stop("SWRI TOTAL spectra/response row counts are inconsistent.")
  }

  X <- rbind(x_ll_a, x_ll_b, x_hl)
  y <- c(y_ll_a, y_ll_b, y_hl)
  grp <- c(rep("ll_a", nrow(x_ll_a)), rep("ll_b", nrow(x_ll_b)), rep("hl", nrow(x_hl)))

  wl <- seq(750, by = 2, length.out = ncol(X))
  feat_names <- paste0("wl_", wl)
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "SWRI Diesel NIR (TOTAL, TOTALGATEST)",
    source = "https://www.eigenvector.com/data/SWRI/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "TOTAL",
      subset_labels = c("ll_a", "ll_b", "hl"),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_swri_visc_high_dim <- function() {
  slug <- "swri_visc_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  zip_path <- file.path(raw_dir, "viscgatest.zip")
  mat_path <- file.path(raw_dir, "VISCGATEST.mat")
  download_if_missing("https://www.eigenvector.com/data/SWRI/viscgatest.zip", zip_path)
  if (!file.exists(mat_path)) utils::unzip(zip_path, exdir = raw_dir)

  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("SWRI VISC requires package 'R.matlab'. Install it and rerun.")
  }

  m <- R.matlab::readMat(mat_path)
  required_fields <- c("v.sd.hl", "v.y.hl", "v.sd.ll.a", "v.sd.ll.b", "v.y.ll.a", "v.y.ll.b")
  missing_fields <- setdiff(required_fields, names(m))
  if (length(missing_fields) > 0L) {
    stop("SWRI VISC .mat missing expected fields: ", paste(missing_fields, collapse = ", "))
  }

  x_hl <- as.matrix(m$v.sd.hl)
  y_hl <- as.numeric(m$v.y.hl)
  x_ll_a <- as.matrix(m$v.sd.ll.a)
  y_ll_a <- as.numeric(m$v.y.ll.a)
  x_ll_b <- as.matrix(m$v.sd.ll.b)
  y_ll_b <- as.numeric(m$v.y.ll.b)

  if (ncol(x_hl) != ncol(x_ll_a) || ncol(x_hl) != ncol(x_ll_b)) {
    stop("SWRI VISC spectra dimensions are inconsistent across subsets.")
  }
  if (nrow(x_hl) != length(y_hl) || nrow(x_ll_a) != length(y_ll_a) || nrow(x_ll_b) != length(y_ll_b)) {
    stop("SWRI VISC spectra/response row counts are inconsistent.")
  }

  X <- rbind(x_ll_a, x_ll_b, x_hl)
  y <- c(y_ll_a, y_ll_b, y_hl)
  grp <- c(rep("ll_a", nrow(x_ll_a)), rep("ll_b", nrow(x_ll_b)), rep("hl", nrow(x_hl)))

  wl <- seq(750, by = 2, length.out = ncol(X))
  feat_names <- paste0("wl_", wl)
  d <- data.frame(group_raw = grp, y = y, X, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[3:ncol(d)] <- feat_names

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feat_names,
    dataset_name = "SWRI Diesel NIR (VISC, VISCGATEST)",
    source = "https://www.eigenvector.com/data/SWRI/",
    min_group_size = 2L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = "VISC",
      subset_labels = c("ll_a", "ll_b", "hl"),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_nci60_cellminer_high_dim <- function() {
  slug <- "nci60_cellminer_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  expr_zip <- file.path(raw_dir, "nci60_rna_seq_composite_expression.zip")
  drug_zip <- file.path(raw_dir, "dtp_nci60_zscore.zip")
  download_if_missing(
    "https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__RNA_seq_composite_expression.zip",
    expr_zip
  )
  download_if_missing(
    "https://discover.nci.nih.gov/cellminer/download/processeddataset/DTP_NCI60_ZSCORE.zip",
    drug_zip
  )

  expr_dir <- file.path(raw_dir, "nci60_rna_seq")
  drug_dir <- file.path(raw_dir, "dtp_zscore")
  if (!dir.exists(expr_dir)) utils::unzip(expr_zip, exdir = expr_dir)
  if (!dir.exists(drug_dir)) utils::unzip(drug_zip, exdir = drug_dir)

  expr_path <- file.path(expr_dir, "output", "RNA__RNA_seq_composite_expression.xls")
  drug_path <- file.path(drug_dir, "output", "DTP_NCI60_ZSCORE.xlsx")
  meta_path <- file.path(expr_dir, "output", "documentation", "NCI60_CELL_LINE_METADATA.xls")
  if (!all(file.exists(c(expr_path, drug_path, meta_path)))) {
    stop("Missing expected NCI-60 CellMiner files after unzip.")
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("NCI-60 CellMiner requires package 'readxl'. Install it and rerun.")
  }

  expr <- readxl::read_excel(expr_path, sheet = "Results", skip = 10)
  drug <- readxl::read_excel(drug_path, sheet = "all", skip = 8)
  meta <- readxl::read_excel(meta_path, sheet = "clc", skip = 7)

  normalize_cell <- function(x) {
    toupper(gsub("[^A-Za-z0-9]", "", x))
  }
  sanitize <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

  expr_cell_cols <- setdiff(names(expr), names(expr)[1:6])
  drug_cell_cols <- setdiff(names(drug), names(drug)[1:6])
  meta_cell_col <- "Cell Line Name"
  meta_group_col <- "tissue of origin a"
  if (!(meta_cell_col %in% names(meta)) || !(meta_group_col %in% names(meta))) {
    stop("NCI-60 metadata columns not found.")
  }

  expr_norm <- normalize_cell(expr_cell_cols)
  drug_norm <- normalize_cell(drug_cell_cols)
  meta_norm <- normalize_cell(meta[[meta_cell_col]])

  common_norm <- Reduce(intersect, list(expr_norm, drug_norm, meta_norm))
  if (length(common_norm) < 40L) {
    stop("Insufficient matched NCI-60 cell lines across expression/drug/metadata.")
  }

  idx_expr <- match(common_norm, expr_norm)
  idx_drug <- match(common_norm, drug_norm)
  idx_meta <- match(common_norm, meta_norm)
  common_cells <- expr_cell_cols[idx_expr]

  x_mat <- t(as.matrix(expr[, idx_expr + 6L, drop = FALSE]))
  storage.mode(x_mat) <- "double"

  gene_name <- as.character(expr[["Gene name d"]])
  gene_id <- suppressWarnings(as.integer(expr[["Entrez gene id e"]]))
  feat_names <- paste0(
    "g_",
    ifelse(!is.na(gene_id), gene_id, seq_len(nrow(expr))),
    "_",
    sanitize(ifelse(is.na(gene_name), "unknown", gene_name))
  )
  feat_names <- make.unique(feat_names, sep = "_dup")
  colnames(x_mat) <- feat_names
  rownames(x_mat) <- common_cells

  drug_df <- drug[, idx_drug + 6L, drop = FALSE]
  drug_num <- suppressWarnings(do.call(cbind, lapply(drug_df, as.numeric)))
  rownames(drug_num) <- seq_len(nrow(drug))
  n_non_na <- rowSums(is.finite(drug_num))
  sd_row <- apply(drug_num, 1L, function(z) stats::sd(z, na.rm = TRUE))
  drug_name <- as.character(drug[["Drug name"]])
  valid_name <- !is.na(drug_name) & nzchar(trimws(drug_name)) & trimws(drug_name) != "-"
  score <- ifelse(valid_name & n_non_na >= 45L, n_non_na + 1e-6 * sd_row, -Inf)
  best_idx <- which.max(score)
  if (!is.finite(score[best_idx])) {
    stop("Could not find a suitable drug response row with enough observations.")
  }
  y_all <- as.numeric(drug_num[best_idx, ])
  keep <- is.finite(y_all)
  if (sum(keep) < 40L) stop("Selected drug has too few finite response values.")

  group_raw <- as.character(meta[[meta_group_col]][idx_meta])[keep]
  y <- y_all[keep]
  X <- x_mat[keep, , drop = FALSE]

  group_counts <- table(group_raw)
  keep_groups <- names(group_counts[group_counts >= 3L])
  keep2 <- group_raw %in% keep_groups
  if (sum(keep2) < 30L) stop("Too few samples after group size filtering for NCI-60.")

  group_raw <- group_raw[keep2]
  y <- y[keep2]
  X <- X[keep2, , drop = FALSE]
  groups <- as.integer(factor(group_raw))

  graphs <- build_graphs(length(unique(groups)), dense_limit = 400L)
  obj <- list(
    X = X,
    y = y,
    groups = groups,
    group_levels = sort(unique(group_raw)),
    feature_names = colnames(X),
    G_dense = graphs$G_dense,
    G_sparse_chain = graphs$G_sparse_chain,
    meta = list(
      dataset = "NCI-60 CellMiner (RNA-seq + drug response)",
      source = "https://discover.nci.nih.gov/cellminer/loadDownload.do",
      n = nrow(X),
      p = ncol(X),
      k = length(unique(groups)),
      y_col = "drug_zscore",
      group_col = "tissue_of_origin",
      dense_graph_limit = graphs$dense_limit,
      target_selected = trimws(as.character(drug[["Drug name"]][best_idx])),
      target_nsc = as.character(drug[["NSC # b"]][best_idx]),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_ccle_depmap_gdsc_high_dim <- function() {
  slug <- "ccle_depmap_gdsc_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  catalog_path <- file.path(raw_dir, "depmap_data_catalog.json")
  catalog <- depmap_catalog(catalog_path)

  expr_release <- "DepMap Public 25Q3"
  gdsc_release <- "Sanger GDSC1 and GDSC2"
  expr_file <- "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
  model_file <- "Model.csv"
  gdsc_file <- "sanger-dose-response.csv"

  expr_url <- depmap_download_url(catalog, expr_file, expr_release)
  model_url <- depmap_download_url(catalog, model_file, expr_release)
  gdsc_url <- depmap_download_url(catalog, gdsc_file, gdsc_release)

  expr_path <- file.path(raw_dir, expr_file)
  model_path <- file.path(raw_dir, model_file)
  gdsc_path <- file.path(raw_dir, gdsc_file)
  download_if_missing(expr_url, expr_path)
  download_if_missing(model_url, model_path)
  download_if_missing(gdsc_url, gdsc_path)

  expr <- utils::read.csv(expr_path, check.names = FALSE, stringsAsFactors = FALSE)
  model <- utils::read.csv(model_path, check.names = FALSE, stringsAsFactors = FALSE)
  gdsc <- utils::read.csv(gdsc_path, check.names = FALSE, stringsAsFactors = FALSE)

  req_expr <- c("ModelID")
  req_model <- c("ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "TissueOrigin")
  req_gdsc <- c("ARXSPAN_ID", "DRUG_NAME", "auc")
  if (!all(req_expr %in% names(expr))) stop("DepMap expression file missing required columns.")
  if (!all(req_model %in% names(model))) stop("DepMap model file missing required columns.")
  if (!all(req_gdsc %in% names(gdsc))) stop("DepMap GDSC file missing required columns.")

  if ("IsDefaultEntryForModel" %in% names(expr)) {
    keep_default <- tolower(as.character(expr$IsDefaultEntryForModel)) %in% c("true", "t", "1")
    n_per_model <- ave(seq_len(nrow(expr)), expr$ModelID, FUN = seq_along)
    expr <- expr[keep_default | n_per_model == 1L, , drop = FALSE]
  }
  expr <- expr[!duplicated(expr$ModelID), , drop = FALSE]

  gdsc$ARXSPAN_ID <- as.character(gdsc$ARXSPAN_ID)
  gdsc$DRUG_NAME <- trimws(as.character(gdsc$DRUG_NAME))
  gdsc$auc <- suppressWarnings(as.numeric(gdsc$auc))
  gdsc <- gdsc[
    is.finite(gdsc$auc) & nzchar(gdsc$DRUG_NAME) & !is.na(gdsc$ARXSPAN_ID),
    c("ARXSPAN_ID", "DRUG_NAME", "auc"),
    drop = FALSE
  ]
  gdsc_agg <- stats::aggregate(auc ~ ARXSPAN_ID + DRUG_NAME, data = gdsc, FUN = mean)

  common_models <- intersect(as.character(expr$ModelID), as.character(gdsc_agg$ARXSPAN_ID))
  if (length(common_models) < 50L) {
    stop("Too few overlapping models between DepMap expression and GDSC.")
  }
  gdsc_agg <- gdsc_agg[gdsc_agg$ARXSPAN_ID %in% common_models, , drop = FALSE]

  drug_n_df <- stats::aggregate(auc ~ DRUG_NAME, data = gdsc_agg, FUN = length)
  names(drug_n_df)[2] <- "n"
  drug_sd_df <- stats::aggregate(auc ~ DRUG_NAME, data = gdsc_agg, FUN = stats::sd)
  names(drug_sd_df)[2] <- "sd"
  drug_stats <- merge(drug_n_df, drug_sd_df, by = "DRUG_NAME", all = TRUE)
  drug_n <- as.numeric(drug_stats$n)
  drug_sd <- as.numeric(drug_stats$sd)
  score <- ifelse(drug_n >= 200L & is.finite(drug_sd), drug_n + 1e-6 * drug_sd, -Inf)
  if (!any(is.finite(score))) {
    score <- ifelse(is.finite(drug_sd), drug_n + 1e-6 * drug_sd, -Inf)
  }
  best_idx <- which.max(score)
  if (!is.finite(score[best_idx])) {
    stop("Could not select a suitable GDSC drug target.")
  }
  target_selected <- as.character(drug_stats$DRUG_NAME[best_idx])

  y_df <- gdsc_agg[gdsc_agg$DRUG_NAME == target_selected, c("ARXSPAN_ID", "auc"), drop = FALSE]
  names(y_df) <- c("ModelID", "y")

  expr$ModelID <- as.character(expr$ModelID)
  d <- merge(y_df, expr, by = "ModelID", all = FALSE)
  if (nrow(d) < 100L) stop("Too few rows after joining selected GDSC target with expression.")

  model_sub <- unique(model[, req_model, drop = FALSE])
  model_sub$ModelID <- as.character(model_sub$ModelID)
  d <- merge(d, model_sub, by = "ModelID", all.x = TRUE)

  group_raw <- as.character(d$OncotreeLineage)
  miss <- is.na(group_raw) | !nzchar(group_raw)
  group_raw[miss] <- as.character(d$OncotreePrimaryDisease[miss])
  miss <- is.na(group_raw) | !nzchar(group_raw)
  group_raw[miss] <- as.character(d$TissueOrigin[miss])
  miss <- is.na(group_raw) | !nzchar(group_raw)
  group_raw[miss] <- "unknown"
  d$group_raw <- group_raw

  meta_cols <- c(
    "ModelID", "y", "group_raw",
    "SequencingID", "IsDefaultEntryForModel", "ModelConditionID", "IsDefaultEntryForMC",
    "OncotreeLineage", "OncotreePrimaryDisease", "TissueOrigin"
  )
  feature_cols <- setdiff(names(d), meta_cols)
  if (length(feature_cols) == 0L) stop("No expression feature columns found in CCLE/DepMap data.")

  obj <- finalize_grouped(
    d,
    y_col = "y",
    group_col = "group_raw",
    feature_cols = feature_cols,
    dataset_name = "CCLE + DepMap Expression + Sanger GDSC (single-drug)",
    source = "https://depmap.org/portal/download/all/",
    min_group_size = 5L,
    max_rows = NULL,
    extra_meta = list(
      target_selected = target_selected,
      expression_file = expr_file,
      expression_release = expr_release,
      response_file = gdsc_file,
      response_release = gdsc_release,
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_tcga_pan_cancer_high_dim <- function() {
  slug <- "tcga_pan_cancer_high_dim"
  raw_dir <- file.path(raw_root, slug)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  studies <- c(
    "brca_tcga_pan_can_atlas_2018",
    "coadread_tcga_pan_can_atlas_2018",
    "hnsc_tcga_pan_can_atlas_2018",
    "kirc_tcga_pan_can_atlas_2018",
    "luad_tcga_pan_can_atlas_2018",
    "lusc_tcga_pan_can_atlas_2018"
  )
  expr_file <- "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
  clin_file <- "data_clinical_patient.txt"
  max_per_group <- 200L

  X_list <- list()
  y_list <- list()
  g_list <- list()
  per_study_n <- integer(0)
  common_gene_names <- NULL

  for (study in studies) {
    study_dir <- file.path(raw_dir, study)
    dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

    expr_path <- file.path(study_dir, expr_file)
    clin_path <- file.path(study_dir, clin_file)
    download_github_raw_if_missing(file.path("public", study, expr_file), expr_path)
    download_github_raw_if_missing(file.path("public", study, clin_file), clin_path)

    clin <- read_cbio_clinical_patient(clin_path)
    if (!all(c("PATIENT_ID", "AGE") %in% names(clin))) {
      stop("TCGA clinical file missing PATIENT_ID/AGE for study: ", study)
    }
    clin$PATIENT_ID <- as.character(clin$PATIENT_ID)
    clin$AGE <- suppressWarnings(as.numeric(clin$AGE))
    clin <- clin[is.finite(clin$AGE) & !is.na(clin$PATIENT_ID) & nzchar(clin$PATIENT_ID), , drop = FALSE]
    if (nrow(clin) == 0L) next

    expr <- read_cbio_expression(expr_path)
    if (ncol(expr) < 3L) stop("Unexpected TCGA expression shape for study: ", study)

    symbol <- as.character(expr[[1L]])
    entrez <- suppressWarnings(as.integer(expr[[2L]]))
    gene_names <- paste0(
      "g_",
      ifelse(is.finite(entrez), entrez, seq_len(nrow(expr))),
      "_",
      sanitize_name(ifelse(is.na(symbol) | !nzchar(symbol), "unknown", symbol))
    )
    gene_names <- make.unique(gene_names, sep = "_dup")

    sample_ids <- names(expr)[-(1:2)]
    patient_ids <- substr(sample_ids, 1L, 12L)
    sample_type <- substr(sample_ids, 14L, 15L)
    idx_by_patient <- split(seq_along(patient_ids), patient_ids)
    chosen <- unlist(lapply(idx_by_patient, function(ix) {
      ix_tumor <- ix[sample_type[ix] == "01"]
      if (length(ix_tumor) > 0L) ix_tumor[1L] else ix[1L]
    }), use.names = FALSE)
    chosen <- sort(unique(chosen))

    chosen_patients <- patient_ids[chosen]
    age_map <- setNames(clin$AGE, clin$PATIENT_ID)
    valid <- !is.na(age_map[chosen_patients]) & is.finite(age_map[chosen_patients])
    chosen <- chosen[valid]
    chosen_patients <- chosen_patients[valid]
    if (length(chosen) < 20L) next

    x_block <- as.matrix(expr[, chosen + 2L, drop = FALSE])
    storage.mode(x_block) <- "double"
    X <- t(x_block)
    rm(x_block, expr)
    gc(verbose = FALSE)

    if (is.null(common_gene_names)) {
      common_gene_names <- gene_names
    } else if (!identical(common_gene_names, gene_names)) {
      stop("Gene feature ordering differs across TCGA studies; add alignment logic before continuing.")
    }
    colnames(X) <- common_gene_names

    y <- as.numeric(age_map[chosen_patients])
    group_raw <- if ("CANCER_TYPE_ACRONYM" %in% names(clin)) {
      as.character(setNames(clin$CANCER_TYPE_ACRONYM, clin$PATIENT_ID)[chosen_patients])
    } else {
      rep(toupper(sub("_tcga_pan_can_atlas_2018$", "", study)), length(chosen_patients))
    }
    miss_grp <- is.na(group_raw) | !nzchar(group_raw)
    group_raw[miss_grp] <- toupper(sub("_tcga_pan_can_atlas_2018$", "", study))

    if (nrow(X) > max_per_group) {
      keep_idx <- sort(sample.int(nrow(X), max_per_group))
      X <- X[keep_idx, , drop = FALSE]
      y <- y[keep_idx]
      group_raw <- group_raw[keep_idx]
    }

    X_list[[study]] <- X
    y_list[[study]] <- y
    g_list[[study]] <- group_raw
    per_study_n[study] <- nrow(X)
  }

  if (length(X_list) < 3L) {
    stop("Too few TCGA studies successfully processed.")
  }

  X_all <- do.call(rbind, X_list)
  y_all <- unlist(y_list, use.names = FALSE)
  g_all <- unlist(g_list, use.names = FALSE)

  if (anyNA(X_all)) {
    for (j in seq_len(ncol(X_all))) {
      v <- X_all[, j]
      if (anyNA(v)) {
        med <- stats::median(v, na.rm = TRUE)
        if (!is.finite(med)) med <- 0
        v[is.na(v)] <- med
        X_all[, j] <- v
      }
    }
  }

  keep_var <- apply(X_all, 2L, function(z) stats::sd(z) > 0)
  X_all <- X_all[, keep_var, drop = FALSE]
  if (ncol(X_all) == 0L) stop("No non-constant TCGA expression features after filtering.")

  groups <- as.integer(factor(g_all))
  group_levels <- levels(factor(g_all))
  graphs <- build_graphs(length(unique(groups)), dense_limit = 400L)
  obj <- list(
    X = X_all,
    y = y_all,
    groups = groups,
    group_levels = group_levels,
    feature_names = colnames(X_all),
    G_dense = graphs$G_dense,
    G_sparse_chain = graphs$G_sparse_chain,
    meta = list(
      dataset = "TCGA PanCancer Atlas (subset, expression z-scores)",
      source = "https://github.com/cBioPortal/datahub/tree/master/public",
      n = nrow(X_all),
      p = ncol(X_all),
      k = length(unique(groups)),
      y_col = "AGE",
      group_col = "CANCER_TYPE_ACRONYM",
      dense_graph_limit = graphs$dense_limit,
      target_selected = "AGE",
      tcga_studies = studies,
      max_per_group = max_per_group,
      per_study_n = as.list(per_study_n),
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

prepare_abide_connectome_high_dim <- function() {
  slug <- "abide_connectome_high_dim"
  raw_dir <- file.path(raw_root, slug)
  ts_dir <- file.path(raw_dir, "rois_cc200")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ts_dir, recursive = TRUE, showWarnings = FALSE)

  pheno_url <- "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Phenotypic_V1_0b_preprocessed1.csv"
  pheno_path <- file.path(raw_dir, "Phenotypic_V1_0b_preprocessed1.csv")
  download_if_missing(pheno_url, pheno_path)
  pheno <- utils::read.csv(pheno_path, stringsAsFactors = FALSE)

  required <- c("FILE_ID", "SITE_ID", "AGE_AT_SCAN")
  if (!all(required %in% names(pheno))) {
    stop("ABIDE phenotypic file missing required columns: ", paste(setdiff(required, names(pheno)), collapse = ", "))
  }

  pheno$FILE_ID <- as.character(pheno$FILE_ID)
  pheno$SITE_ID <- as.character(pheno$SITE_ID)
  pheno$AGE_AT_SCAN <- suppressWarnings(as.numeric(pheno$AGE_AT_SCAN))
  pheno <- pheno[
    !is.na(pheno$FILE_ID) &
      pheno$FILE_ID != "" &
      pheno$FILE_ID != "no_filename" &
      !is.na(pheno$SITE_ID) &
      pheno$SITE_ID != "" &
      is.finite(pheno$AGE_AT_SCAN),
    c("FILE_ID", "SITE_ID", "AGE_AT_SCAN"),
    drop = FALSE
  ]
  if (nrow(pheno) == 0L) stop("No usable ABIDE rows after phenotypic filtering.")

  min_group_size <- 20L
  max_sites <- 10L
  max_per_group <- 35L
  site_counts <- sort(table(pheno$SITE_ID), decreasing = TRUE)
  kept_sites <- names(site_counts[site_counts >= min_group_size])
  if (length(kept_sites) == 0L) stop("No ABIDE sites pass min_group_size filtering.")
  kept_sites <- head(kept_sites, max_sites)
  pheno <- pheno[pheno$SITE_ID %in% kept_sites, , drop = FALSE]

  split_by_site <- split(pheno, pheno$SITE_ID)
  for (site in names(split_by_site)) {
    d <- split_by_site[[site]]
    if (nrow(d) > max_per_group) {
      idx <- sort(sample.int(nrow(d), max_per_group))
      split_by_site[[site]] <- d[idx, , drop = FALSE]
    }
  }
  pheno_sub <- do.call(rbind, split_by_site)
  rownames(pheno_sub) <- NULL

  ts_url <- function(fid) {
    sprintf(
      "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/cpac/filt_global/rois_cc200/%s_rois_cc200.1D",
      fid
    )
  }

  get_connectome_vec <- function(fid) {
    local_path <- file.path(ts_dir, paste0(fid, "_rois_cc200.1D"))
    ok <- TRUE
    tryCatch(download_if_missing(ts_url(fid), local_path), error = function(e) ok <<- FALSE)
    if (!ok || !file.exists(local_path)) return(NULL)

    ts <- tryCatch(
      utils::read.table(local_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(ts) || ncol(ts) < 50L || nrow(ts) < 20L) return(NULL)
    M <- as.matrix(ts)
    storage.mode(M) <- "double"
    keep_var <- apply(M, 2L, function(z) stats::sd(z, na.rm = TRUE) > 0)
    M <- M[, keep_var, drop = FALSE]
    if (ncol(M) < 150L) return(NULL)

    C <- stats::cor(M, use = "pairwise.complete.obs")
    C[!is.finite(C)] <- 0
    diag(C) <- 0
    C[upper.tri(C, diag = FALSE)]
  }

  rows <- vector("list", nrow(pheno_sub))
  used <- 0L
  expected_p <- NULL
  feature_names <- NULL
  for (i in seq_len(nrow(pheno_sub))) {
    fid <- pheno_sub$FILE_ID[i]
    vec <- get_connectome_vec(fid)
    if (is.null(vec)) next
    if (is.null(expected_p)) {
      expected_p <- length(vec)
      feature_names <- paste0("conn_", seq_len(expected_p))
    }
    if (length(vec) != expected_p) next
    used <- used + 1L
    rows[[used]] <- list(
      y = as.numeric(pheno_sub$AGE_AT_SCAN[i]),
      group_raw = as.character(pheno_sub$SITE_ID[i]),
      x = as.numeric(vec)
    )
  }

  if (used < 200L) stop("Too few ABIDE subjects with usable connectome vectors.")
  rows <- rows[seq_len(used)]
  y <- vapply(rows, function(r) r$y, numeric(1))
  grp <- vapply(rows, function(r) r$group_raw, character(1))
  X <- do.call(rbind, lapply(rows, function(r) r$x))
  colnames(X) <- feature_names

  # Final group-size filter after availability checks.
  counts <- table(grp)
  keep_groups <- names(counts[counts >= min_group_size])
  keep <- grp %in% keep_groups
  X <- X[keep, , drop = FALSE]
  y <- y[keep]
  grp <- grp[keep]
  if (nrow(X) < 150L) stop("Too few ABIDE rows after final group-size filtering.")

  groups <- as.integer(factor(grp))
  group_levels <- levels(factor(grp))
  graphs <- build_graphs(length(unique(groups)), dense_limit = 400L)
  obj <- list(
    X = X,
    y = y,
    groups = groups,
    group_levels = group_levels,
    feature_names = colnames(X),
    G_dense = graphs$G_dense,
    G_sparse_chain = graphs$G_sparse_chain,
    meta = list(
      dataset = "ABIDE I Preprocessed Connectome (CC200, CPAC filt_global)",
      source = "https://preprocessed-connectomes-project.org/abide/",
      n = nrow(X),
      p = ncol(X),
      k = length(unique(groups)),
      y_col = "AGE_AT_SCAN",
      group_col = "SITE_ID",
      dense_graph_limit = graphs$dense_limit,
      target_selected = "AGE",
      pipeline = "cpac",
      strategy = "filt_global",
      derivative = "rois_cc200",
      min_group_size = min_group_size,
      max_sites = max_sites,
      max_per_group = max_per_group,
      p_gt_n = TRUE,
      regime = "real_high_dim"
    )
  )

  out <- save_grouped(obj, slug)
  list(path = out, n = nrow(obj$X), p = ncol(obj$X), k = length(unique(obj$groups)))
}

dataset_fns <- list(
  communities_crime = prepare_communities,
  wine_quality = prepare_wine_quality,
  student_performance = prepare_student_performance,
  beijing_air_quality = prepare_beijing_air,
  electricity_load = prepare_electricity,
  owid_co2 = prepare_owid_co2,
  world_bank_wdi = prepare_world_bank_wdi,
  nyc_tlc = prepare_nyc_tlc,
  noaa_ghcn = prepare_noaa_ghcn,
  countyplus = prepare_countyplus,
  sarcos_openml = prepare_sarcos_openml,
  school_ilea = prepare_school_ilea,
  nir_corn_high_dim = prepare_nir_corn,
  nir_tablets_high_dim = prepare_nir_tablets_high_dim,
  nir_tablets_hardness_high_dim = prepare_nir_tablets_hardness_high_dim,
  swri_cn_high_dim = prepare_swri_cn_high_dim,
  swri_bp50_high_dim = prepare_swri_bp50_high_dim,
  swri_d4052_high_dim = prepare_swri_d4052_high_dim,
  swri_freeze_high_dim = prepare_swri_freeze_high_dim,
  swri_total_high_dim = prepare_swri_total_high_dim,
  swri_visc_high_dim = prepare_swri_visc_high_dim,
  nci60_cellminer_high_dim = prepare_nci60_cellminer_high_dim,
  ccle_depmap_gdsc_high_dim = prepare_ccle_depmap_gdsc_high_dim,
  tcga_pan_cancer_high_dim = prepare_tcga_pan_cancer_high_dim,
  abide_connectome_high_dim = prepare_abide_connectome_high_dim
)

requested <- commandArgs(trailingOnly = TRUE)
requested <- if (length(requested) == 0L || identical(requested, "all")) names(dataset_fns) else requested
unknown <- setdiff(requested, names(dataset_fns))
if (length(unknown) > 0L) stop("Unknown dataset(s): ", paste(unknown, collapse = ", "))

results <- list()
for (nm in requested) {
  cat(sprintf("[%-22s] running...\n", nm))
  fn <- dataset_fns[[nm]]
  row <- list(
    dataset = nm,
    status = "ok",
    message = "",
    path = NA_character_,
    n = NA_integer_,
    p = NA_integer_,
    k = NA_integer_,
    p_gt_n = NA,
    regime = NA_character_
  )
  out <- tryCatch({
    out <- fn()
    list(ok = TRUE, out = out, msg = "")
  }, error = function(e) {
    list(ok = FALSE, out = NULL, msg = conditionMessage(e))
  })
  if (isTRUE(out$ok)) {
    row$path <- out$out$path %||% NA_character_
    row$n <- out$out$n %||% NA_integer_
    row$p <- out$out$p %||% NA_integer_
    row$k <- out$out$k %||% NA_integer_
    row$p_gt_n <- is.finite(row$p) && is.finite(row$n) && (row$p > row$n)
    row$regime <- if (isTRUE(row$p_gt_n)) "real_high_dim" else "real_standard"
    if (!is.null(out$out$target)) row$message <- paste("target:", out$out$target)
    cat(sprintf("[%-22s] ok  n=%s p=%s k=%s\n", nm, row$n, row$p, row$k))
  } else {
    row$status <- "failed"
    row$message <- out$msg
    cat(sprintf("[%-22s] failed: %s\n", nm, row$message))
  }
  results[[length(results) + 1L]] <- row
}

build_summary_from_processed <- function(processed_root) {
  license_record_for_dataset <- function(dataset_slug, source_url) {
    default <- list(
      license_status = "unknown",
      license_name = NA_character_,
      license_link = NA_character_,
      license_notes = NA_character_
    )

    uci_links <- c(
      communities_crime = "https://archive.ics.uci.edu/dataset/183/communities+and+crime",
      wine_quality = "https://archive.ics.uci.edu/dataset/186/wine+quality",
      student_performance = "https://archive.ics.uci.edu/dataset/320/student+performance",
      beijing_air_quality = "https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data",
      electricity_load = "https://archive.ics.uci.edu/dataset/321/electricityloaddiagrams20112014"
    )
    if (dataset_slug %in% names(uci_links)) {
      return(list(
        license_status = "open",
        license_name = "CC BY 4.0",
        license_link = uci_links[[dataset_slug]],
        license_notes = "UCI dataset page indicates Creative Commons Attribution 4.0."
      ))
    }

    if (dataset_slug == "owid_co2") {
      return(list(
        license_status = "mixed_open_with_exceptions",
        license_name = "CC BY 4.0 (OWID-produced components)",
        license_link = "https://github.com/owid/co2-data",
        license_notes = "OWID data are CC BY 4.0; some upstream source components may carry separate terms."
      ))
    }

    if (dataset_slug == "world_bank_wdi") {
      return(list(
        license_status = "mixed_open_with_exceptions",
        license_name = "CC BY 4.0 (default)",
        license_link = "https://www.worldbank.org/en/about/legal/terms-of-use-for-datasets",
        license_notes = "World Bank data are CC BY 4.0 by default, with possible third-party indicator exceptions."
      ))
    }

    if (dataset_slug == "countyplus") {
      return(list(
        license_status = "open",
        license_name = "MIT",
        license_link = "https://github.com/Clpr/CountyPlus/releases/tag/v0.0.2",
        license_notes = "Repository license file is MIT."
      ))
    }

    if (dataset_slug == "sarcos_openml") {
      return(list(
        license_status = "open",
        license_name = "Public (OpenML metadata)",
        license_link = "https://www.openml.org/api/v1/json/data/43873",
        license_notes = "OpenML metadata lists dataset licence as Public."
      ))
    }

    if (dataset_slug == "school_ilea") {
      return(list(
        license_status = "open",
        license_name = "GPL-2 | GPL-3",
        license_link = "https://cran.r-project.org/package=nlme",
        license_notes = "Derived from datasets shipped in CRAN nlme package."
      ))
    }

    if (dataset_slug == "nyc_tlc") {
      return(list(
        license_status = "open",
        license_name = "NYC Open Data terms",
        license_link = "https://opendata.cityofnewyork.us/overview/",
        license_notes = "Subject to NYC Open Data terms of use."
      ))
    }

    if (dataset_slug == "noaa_ghcn") {
      return(list(
        license_status = "open",
        license_name = "NOAA open/public data",
        license_link = "https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00861",
        license_notes = "Open government data with attribution/citation guidance."
      ))
    }

    if (dataset_slug == "tcga_pan_cancer_high_dim") {
      return(list(
        license_status = "open",
        license_name = "ODbL",
        license_link = "https://github.com/cBioPortal/datahub",
        license_notes = "cBioPortal datahub states ODbL."
      ))
    }

    if (dataset_slug == "abide_connectome_high_dim") {
      return(list(
        license_status = "restricted_noncommercial",
        license_name = "Attribution Non-Commercial",
        license_link = "https://www.nitrc.org/projects/fcon_1000",
        license_notes = "ABIDE/FCP terms indicate non-commercial restriction."
      ))
    }

    if (dataset_slug == "nci60_cellminer_high_dim") {
      return(list(
        license_status = "unclear_contact_required",
        license_name = "CellMiner terms/citation-required",
        license_link = "https://discover.nci.nih.gov/cellminer/loadDownload.do",
        license_notes = "NIH/NCI pages include citation and third-party rights caveats; confirm reuse terms before publication redistribution."
      ))
    }

    if (dataset_slug == "ccle_depmap_gdsc_high_dim") {
      return(list(
        license_status = "restricted_noncommercial",
        license_name = "Research/non-commercial terms",
        license_link = paste(
          "https://depmap.org/portal/terms/",
          "https://www.cancerrxgene.org/downloads",
          sep = "; "
        ),
        license_notes = "DepMap and GDSC data have research/non-commercial usage restrictions."
      ))
    }

    if (startsWith(dataset_slug, "nir_") || startsWith(dataset_slug, "swri_")) {
      eigenvector_license_pages <- c(
        nir_corn_high_dim = "https://eigenvector.com/resources/data-sets/nir-of-corn-samples-for-standardization-benchmarking/",
        nir_tablets_high_dim = "https://eigenvector.com/resources/data-sets/nir-spectra-of-pharmaceutical-tablets-from-shootout/",
        nir_tablets_hardness_high_dim = "https://eigenvector.com/resources/data-sets/nir-spectra-of-pharmaceutical-tablets-from-shootout/",
        swri_bp50_high_dim = "https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/",
        swri_cn_high_dim = "https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/",
        swri_d4052_high_dim = "https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/",
        swri_freeze_high_dim = "https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/",
        swri_total_high_dim = "https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/",
        swri_visc_high_dim = "https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/"
      )
      lic_link <- eigenvector_license_pages[[dataset_slug]] %||% as.character(source_url)
      return(list(
        license_status = "unclear_contact_required",
        license_name = "No explicit downstream reuse license",
        license_link = lic_link,
        license_notes = paste(
          "Pages provide provenance and indicate permission for Eigenvector to post/distribute,",
          "but do not publish a clear downstream open-data reuse license."
        )
      ))
    }

    default
  }

  meta_files <- list.files(
    processed_root,
    pattern = "grouped_regression_meta\\.json$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(meta_files) == 0L) {
    return(data.frame(
      dataset = character(0),
      status = character(0),
      message = character(0),
      path = character(0),
      n = integer(0),
      p = integer(0),
      k = integer(0),
      p_gt_n = logical(0),
      regime = character(0),
      source = character(0),
      license_status = character(0),
      license_name = character(0),
      license_link = character(0),
      license_notes = character(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(meta_files, function(meta_path) {
    meta <- jsonlite::fromJSON(meta_path)
    dataset_slug <- basename(dirname(meta_path))
    n <- as.integer(meta$n %||% NA_integer_)
    p <- as.integer(meta$p %||% NA_integer_)
    k <- as.integer(meta$k %||% NA_integer_)
    p_gt_n <- if (is.finite(n) && is.finite(p)) p > n else NA
    regime <- meta$regime %||% if (isTRUE(p_gt_n)) "real_high_dim" else "real_standard"
    msg <- if (!is.null(meta$target_selected)) paste("target:", meta$target_selected) else ""
    source <- meta$source %||% NA_character_
    license <- license_record_for_dataset(dataset_slug, source)
    data.frame(
      dataset = dataset_slug,
      status = "ok",
      message = msg,
      path = file.path(dirname(meta_path), "grouped_regression.csv"),
      n = n,
      p = p,
      k = k,
      p_gt_n = p_gt_n,
      regime = regime,
      source = as.character(source),
      license_status = as.character(license$license_status),
      license_name = as.character(license$license_name),
      license_link = as.character(license$license_link),
      license_notes = as.character(license$license_notes),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$dataset), , drop = FALSE]
}

summary_df <- build_summary_from_processed(processed_root)

run_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
failed_df <- run_df[run_df$status == "failed", c("dataset", "status", "message", "path", "n", "p", "k", "p_gt_n", "regime"), drop = FALSE]
if (nrow(failed_df) > 0L) {
  failed_df$source <- NA_character_
  failed_df$license_status <- "unknown"
  failed_df$license_name <- NA_character_
  failed_df$license_link <- NA_character_
  failed_df$license_notes <- NA_character_
  failed_df <- failed_df[, colnames(summary_df), drop = FALSE]
}
if (nrow(failed_df) > 0L) {
  missing_failed <- failed_df[!(failed_df$dataset %in% summary_df$dataset), , drop = FALSE]
  if (nrow(missing_failed) > 0L) {
    summary_df <- rbind(summary_df, missing_failed)
    summary_df <- summary_df[order(summary_df$dataset), , drop = FALSE]
  }
}

summary_path <- file.path(processed_root, "grouped_regression_datasets_summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)

cat("\n=== Summary ===\n")
print(summary_df, row.names = FALSE)
cat("\nSummary file:", summary_path, "\n")
