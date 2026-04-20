# =============================================================================
# I wrote this to drive the "My model" overlay on the SBE LiDAR composite carbon
# bar panel (05_sbe_lidar_maps_bars_composite.R): I read my modelled annual ACD
# change + uncertainty, map habitats onto the same display strata as the LiDAR
# bars, and I nudge/offset the points in the composite so they do not sit on top
# of the grey LiDAR error bars. I do not re-run the forest model from here; I
# just ingest exported CSVs from Inputs.
#
# I try `Inputs/carbon_recovery__quick_compare_model_vs_philipson.csv` first
# (I filter `source == "my_model"` and the annual ACD-increase metric when the
# table is in long form). If that is missing or empty, I fall back to
# `carbon_recovery__annual_acd_inc__summary.csv` and my older column-matching
# rules. I map a single "once-logged" model line onto three bar rows
# (enrichment, RIL, pooled) because I want the same my-model line under those
# LiDAR strata for an apples-to-apples read.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# I return the colour I use for my model (muted red) so I stay off the Okabe-Ito
# map / bar class colours in the same figure.
sbe_carbon_model_overlay_colour <- function() {
  "#924242"
}

# I nudge the model in Mg C ha-1 yr-1 to the *right* of the LiDAR CIs in x so the
# two intervals are visible side by side.
sbe_carbon_model_x_jitter <- function() {
  0.045
}

# I nudge the model slightly on the *discrete* y scale (position_nudge in the
# composite) so the red line sits a hair above the LiDAR error bars.
sbe_carbon_model_y_nudge <- function() {
  0.11
}

# I pick the first column in `candidates` that exists; I do this so I can
# re-use the same code when my CSVs use different naming.
sbe_pick_col <- function(df, candidates) {
  nn <- names(df)
  low <- tolower(nn)
  for (cand in candidates) {
    w <- which(low == tolower(cand))
    if (length(w) > 0L) {
      return(df[[nn[w[[1L]]]]])
    }
  }
  NULL
}

# I reduce messy habitat/labels in my input tables to a small set of keys I
# can expand to my composite bar class names below.
sbe_carbon_model_stratum_key <- function(x) {
  x0 <- tolower(stringr::str_squish(as.character(x)))
  dplyr::case_when(
    is.na(x0) | !nzchar(x0) ~ NA_character_,
    stringr::str_detect(x0, "twice") ~ "twice_logged",
    stringr::str_detect(x0, "restored|liana") ~ "restored",
    stringr::str_detect(x0, "primary|prim\\b|^prim") ~ "primary",
    stringr::str_detect(x0, "once") ~ "once_logged",
    TRUE ~ NA_character_
  )
}

# I turn one model row per coarse key (primary, once, twice, restored) into
# the rows ggplot needs—one y position per *bar* on the figure, with my
# once-logged value repeated for enrichment, RIL, and pooled.
sbe_expand_carbon_model_to_bar_classes <- function(d, bar_classes) {
  if (nrow(d) == 0L) {
    return(data.frame())
  }
  prim_disp <- "Primary (Danum)"
  twice_disp <- "Twice logged (SBE)"
  once_disps <- c(
    "Once-logged (enrichment)",
    "Once-logged (RIL)",
    "Once-logged (RIL and enrichment) pooled"
  )
  lev_cls <- rev(bar_classes)
  fac <- function(x) factor(x, levels = lev_cls)
  out <- list()
  if (any(d$key == "primary", na.rm = TRUE) && prim_disp %in% bar_classes) {
    row <- d[d$key == "primary", , drop = FALSE][1L, ]
    out[[length(out) + 1L]] <- data.frame(
      mean = row$mean,
      ymin = row$ymin,
      ymax = row$ymax,
      cls = fac(prim_disp),
      stringsAsFactors = FALSE
    )
  }
  if (any(d$key == "twice_logged", na.rm = TRUE) && twice_disp %in% bar_classes) {
    row <- d[d$key == "twice_logged", , drop = FALSE][1L, ]
    out[[length(out) + 1L]] <- data.frame(
      mean = row$mean,
      ymin = row$ymin,
      ymax = row$ymax,
      cls = fac(twice_disp),
      stringsAsFactors = FALSE
    )
  }
  if (any(d$key == "once_logged", na.rm = TRUE)) {
    row <- d[d$key == "once_logged", , drop = FALSE][1L, ]
    od <- intersect(once_disps, bar_classes)
    if (length(od) > 0L) {
      out[[length(out) + 1L]] <- data.frame(
        mean = rep(row$mean, length(od)),
        ymin = rep(row$ymin, length(od)),
        ymax = rep(row$ymax, length(od)),
        cls = fac(od),
        stringsAsFactors = FALSE
      )
    }
  }
  rest_disp <- "Restored (liana cut)"
  if (any(d$key == "restored", na.rm = TRUE) && rest_disp %in% bar_classes) {
    row <- d[d$key == "restored", , drop = FALSE][1L, ]
    out[[length(out) + 1L]] <- data.frame(
      mean = row$mean,
      ymin = row$ymin,
      ymax = row$ymax,
      cls = fac(rest_disp),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0L) {
    return(data.frame())
  }
  res <- dplyr::bind_rows(out)
  res$cls <- factor(as.character(res$cls), levels = lev_cls)
  res
}

# I add the same x-offset to point + interval so the whole "my model" mark moves
# together; the composite then applies my y nudge in separate layers.
sbe_carbon_model_apply_x_offset <- function(df, offset = NULL) {
  if (nrow(df) == 0L) {
    return(df)
  }
  if (is.null(offset)) {
    offset <- sbe_carbon_model_x_jitter()
  }
  df$mean_plot <- df$mean + offset
  df$ymin_plot <- ifelse(is.finite(df$ymin), df$ymin + offset, NA_real_)
  df$ymax_plot <- ifelse(is.finite(df$ymax), df$ymax + offset, NA_real_)
  df
}

# I accept either a .csv extension or no extension because I was inconsistent
# when I first dropped the file in Inputs.
sbe_carbon_quick_compare_path <- function() {
  b <- "carbon_recovery__quick_compare_model_vs_philipson"
  if (exists("nr2_input_path", inherits = TRUE)) {
    p_csv <- nr2_input_path(paste0(b, ".csv"))
    p_nc <- nr2_input_path(b)
  } else {
    p_csv <- file.path("Inputs", paste0(b, ".csv"))
    p_nc <- file.path("Inputs", b)
  }
  if (file.exists(p_csv)) {
    return(p_csv)
  }
  if (file.exists(p_nc)) {
    return(p_nc)
  }
  NA_character_
}

# I read my quick-compare table: the long form I use has `source`, `estimate`,
# `lwr95`, `upr95`, `habitat`—I filter to my annual ACD-increase rows. I keep a
# second branch for wide "my_model_*" column names in case I export that shape.
sbe_carbon_quick_compare_model_points_for_bar <- function(path, bar_classes) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(data.frame())
  }
  raw <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) {
      warning("Could not read quick-compare file: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0L) {
    return(data.frame())
  }

  nm <- names(raw)
  low <- tolower(nm)

  if ("source" %in% nm && "estimate" %in% nm) {
    src <- tolower(stringr::str_squish(as.character(raw$source)))
    keep_src <- src == "my_model"
    if ("metric" %in% nm) {
      met <- tolower(as.character(raw$metric))
      keep_src <- keep_src & stringr::str_detect(
        met,
        stringr::regex("annual.*acd.*increase|acd.*increase.*mg", ignore_case = TRUE)
      )
    }
    raw <- raw[keep_src, , drop = FALSE]
    if (nrow(raw) == 0L) {
      warning("Quick-compare: no `my_model` rows for annual ACD increase in: ", path)
      return(data.frame())
    }
    st <- sbe_pick_col(raw, c("habitat", "habitat_for_comparison", "stratum", "class"))
    if (is.null(st)) {
      warning("Quick-compare file: no habitat/stratum column found.")
      return(data.frame())
    }
    lo_nm <- if ("lwr95" %in% nm) "lwr95" else if ("lwr_95" %in% nm) "lwr_95" else NA_character_
    hi_nm <- if ("upr95" %in% nm) "upr95" else if ("upr_95" %in% nm) "upr_95" else NA_character_
    d <- data.frame(
      stratum = st,
      mean = suppressWarnings(as.numeric(raw$estimate)),
      ymin = if (!is.na(lo_nm) && nzchar(lo_nm)) suppressWarnings(as.numeric(raw[[lo_nm]])) else NA_real_,
      ymax = if (!is.na(hi_nm) && nzchar(hi_nm)) suppressWarnings(as.numeric(raw[[hi_nm]])) else NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    mm <- nm[grepl("my_model", low, fixed = TRUE)]
    if (length(mm) == 0L) {
      warning("Quick-compare file not in expected format (need source+estimate or my_model* columns): ", path)
      return(data.frame())
    }
    mlow <- tolower(mm)
    pick1 <- function(pat) {
      hit <- mm[grepl(pat, mlow, perl = TRUE)]
      if (length(hit) > 0L) hit[[1L]] else NA_character_
    }
    mean_col <- pick1("mean|estimate|median|pred|expect")
    if (is.na(mean_col) || !nzchar(mean_col)) {
      mean_col <- mm[[1L]]
    }
    lo_col <- pick1("lwr|lower|lo_|ci_lo|ymin|q2\\.5|q025")
    hi_col <- pick1("upr|upper|hi_|ci_hi|ymax|q97\\.5|q975")

    st <- sbe_pick_col(raw, c(
      "habitat_for_comparison", "habitat", "stratum", "class",
      "land_use", "intervention", "forest_type", "site_type"
    ))
    if (is.null(st)) {
      warning("Quick-compare file: no stratum/habitat column found.")
      return(data.frame())
    }

    # I use the same baseline / slope-1.0 filter I use for the annual ACD
    # summary when my table has trajectory and habitat columns.
    if ("trajectory_assumption" %in% names(raw) && "habitat_for_comparison" %in% names(raw)) {
      habl <- tolower(stringr::str_squish(as.character(raw$habitat_for_comparison)))
      tr <- as.character(raw$trajectory_assumption)
      keep <- (
        (habl == "twice-logged" & stringr::str_detect(tr, stringr::fixed("1.0"))) |
        (habl != "twice-logged" & stringr::str_detect(tr, stringr::regex("baseline", ignore_case = TRUE)))
      )
      if (any(keep, na.rm = TRUE)) {
        raw <- raw[keep, , drop = FALSE]
      }
    }

    d <- data.frame(
      stratum = st,
      mean = suppressWarnings(as.numeric(raw[[mean_col]])),
      ymin = if (!is.na(lo_col) && nzchar(lo_col)) suppressWarnings(as.numeric(raw[[lo_col]])) else NA_real_,
      ymax = if (!is.na(hi_col) && nzchar(hi_col)) suppressWarnings(as.numeric(raw[[hi_col]])) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  d <- d %>%
    dplyr::filter(is.finite(.data$mean)) %>%
    dplyr::mutate(key = sbe_carbon_model_stratum_key(.data$stratum)) %>%
    dplyr::filter(!is.na(.data$key)) %>%
    dplyr::distinct(.data$key, .keep_all = TRUE)

  if (nrow(d) == 0L) {
    warning("Quick-compare: no rows matched primary / once-logged / twice-logged / restored.")
    return(data.frame())
  }
  sbe_expand_carbon_model_to_bar_classes(d, bar_classes)
}

# I use this as a backup when the quick-compare file is not there: it matches
# my old annual summary layout (habitat_for_comparison, annual_acd_*, and I
# filter baselines the same way I do in the megatree pipeline for twice-logged).
sbe_carbon_model_points_for_bar <- function(path, bar_classes) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    message("Carbon model overlay skipped (file not found): ", path)
    return(data.frame())
  }
  raw <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) {
      warning("Could not read carbon model summary: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(raw) || nrow(raw) == 0L) {
    return(data.frame())
  }

  if ("habitat_for_comparison" %in% names(raw) && "trajectory_assumption" %in% names(raw)) {
    habl <- tolower(stringr::str_squish(as.character(raw$habitat_for_comparison)))
    tr <- as.character(raw$trajectory_assumption)
    keep <- (
      (habl == "twice-logged" & stringr::str_detect(tr, stringr::fixed("1.0"))) |
      (habl != "twice-logged" & stringr::str_detect(tr, stringr::regex("baseline", ignore_case = TRUE)))
    )
    if (any(keep, na.rm = TRUE)) {
      raw <- raw[keep, , drop = FALSE]
    }
  }

  st <- sbe_pick_col(raw, c(
    "habitat_for_comparison",
    "stratum", "habitat", "class", "land_use", "intervention",
    "model_stratum", "site_type", "forest_type"
  ))
  mn <- sbe_pick_col(raw, c(
    "annual_acd_increase_mean",
    "mean_mgc_ha_yr", "mean_acd_inc_mgc_ha_yr", "mean_annual_acd_inc",
    "mean", "estimate", "prediction",
    "annual_acd_increase_median", "median", "post_mean"
  ))
  lo <- sbe_pick_col(raw, c(
    "annual_acd_increase_lwr95",
    "lwr_95", "lwr", "q2_5", "q025", "conf_low", "lower", "ymin", "ci_low",
    "annual_acd_increase_lwr80"
  ))
  hi <- sbe_pick_col(raw, c(
    "annual_acd_increase_upr95",
    "upr_95", "upr", "q97_5", "q975", "conf_high", "upper", "ymax", "ci_high",
    "annual_acd_increase_upr80"
  ))
  se <- sbe_pick_col(raw, c("se", "std_err", "stderr", "sigma_mean"))

  if (is.null(st) || is.null(mn)) {
    warning(
      "carbon_recovery summary: need stratum-like + mean-like columns; found: ",
      paste(names(raw), collapse = ", ")
    )
    return(data.frame())
  }

  d <- data.frame(
    stratum = st,
    mean = as.numeric(mn),
    ymin = if (!is.null(lo)) as.numeric(lo) else NA_real_,
    ymax = if (!is.null(hi)) as.numeric(hi) else NA_real_,
    stringsAsFactors = FALSE
  )

  if (all(is.na(d$ymin) | is.na(d$ymax)) && !is.null(se)) {
    sev <- as.numeric(se)
    d$ymin <- d$mean - 1.96 * sev
    d$ymax <- d$mean + 1.96 * sev
  }

  d <- d %>%
    dplyr::filter(is.finite(.data$mean)) %>%
    dplyr::mutate(key = sbe_carbon_model_stratum_key(.data$stratum)) %>%
    dplyr::filter(!is.na(.data$key)) %>%
    dplyr::distinct(.data$key, .keep_all = TRUE)

  if (nrow(d) == 0L) {
    warning("carbon_recovery summary: no rows matched primary / once-logged / twice-logged.")
    return(data.frame())
  }
  sbe_expand_carbon_model_to_bar_classes(d, bar_classes)
}

# This is what 05_ sources: I try quick-compare first, and only then my annual
# file.
sbe_carbon_model_points_for_composite <- function(bar_classes) {
  pq <- sbe_carbon_quick_compare_path()
  if (!is.na(pq) && nzchar(pq)) {
    out <- sbe_carbon_quick_compare_model_points_for_bar(pq, bar_classes)
    if (nrow(out) > 0L) {
      message("Carbon model overlay from: ", pq)
      return(out)
    }
  }
  annual <- if (exists("nr2_input_path", inherits = TRUE)) {
    nr2_input_path("carbon_recovery__annual_acd_inc__summary.csv")
  } else {
    file.path("Inputs", "carbon_recovery__annual_acd_inc__summary.csv")
  }
  sbe_carbon_model_points_for_bar(annual, bar_classes)
}
