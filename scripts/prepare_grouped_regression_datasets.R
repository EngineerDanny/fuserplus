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

  meta <- obj$meta
  meta$feature_names <- as.character(obj$feature_names)
  meta$group_levels <- as.character(group_levels)
  meta$output_files <- list(
    data_csv = data_csv,
    sparse_chain_edges_csv = chain_edges_csv,
    dense_edges_csv = dense_edges_csv
  )
  writeLines(
    jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    con = meta_json
  )

  data_csv
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

  n_timestamps <- min(2500L, nrow(d))
  n_clients <- min(120L, ncol(d) - 1L)
  col_idx <- c(1L, seq.int(2L, n_clients + 1L))
  d <- d[seq_len(n_timestamps), col_idx, drop = FALSE]

  times <- as.POSIXct(d[[1]], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (all(is.na(times))) stop("Could not parse timestamp in ElectricityLoadDiagrams.")
  client_cols <- names(d)[-1L]

  long <- do.call(rbind, lapply(client_cols, function(cl) {
    data.frame(group_raw = cl, time = times, load = as.numeric(d[[cl]]), stringsAsFactors = FALSE)
  }))
  long <- long[is.finite(long$load) & !is.na(long$time), , drop = FALSE]
  long <- long[order(long$group_raw, long$time), , drop = FALSE]
  long$lag1 <- ave(long$load, long$group_raw, FUN = function(z) c(NA_real_, head(z, -1L)))
  long$hour <- as.integer(format(long$time, "%H"))
  long$wday <- as.integer(format(long$time, "%u"))
  long$month <- as.integer(format(long$time, "%m"))
  long <- long[!is.na(long$lag1), , drop = FALSE]

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
  school_ilea = prepare_school_ilea
)

requested <- commandArgs(trailingOnly = TRUE)
requested <- if (length(requested) == 0L || identical(requested, "all")) names(dataset_fns) else requested
unknown <- setdiff(requested, names(dataset_fns))
if (length(unknown) > 0L) stop("Unknown dataset(s): ", paste(unknown, collapse = ", "))

results <- list()
for (nm in requested) {
  cat(sprintf("[%-22s] running...\n", nm))
  fn <- dataset_fns[[nm]]
  row <- list(dataset = nm, status = "ok", message = "", path = NA_character_, n = NA_integer_, p = NA_integer_, k = NA_integer_)
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
    if (!is.null(out$out$target)) row$message <- paste("target:", out$out$target)
    cat(sprintf("[%-22s] ok  n=%s p=%s k=%s\n", nm, row$n, row$p, row$k))
  } else {
    row$status <- "failed"
    row$message <- out$msg
    cat(sprintf("[%-22s] failed: %s\n", nm, row$message))
  }
  results[[length(results) + 1L]] <- row
}

summary_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
summary_path <- file.path(processed_root, "grouped_regression_datasets_summary.csv")
utils::write.csv(summary_df, summary_path, row.names = FALSE)

cat("\n=== Summary ===\n")
print(summary_df, row.names = FALSE)
cat("\nSummary file:", summary_path, "\n")
