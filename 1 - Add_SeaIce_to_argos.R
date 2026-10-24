# ============================================================
# Extract Sea Ice concentration at animal tracking locations
# ============================================================
#
# Description:
#
# Extract sea-ice concentration (SIC) at ARGOS locations
# Data source: University of Bremen AMSR2 ASI daily sea-ice concentration grids
#
# Configure the paths and column names in the USER SETTINGS section.
#
# See readme file for more info on project pathin
#
# ---------------------------- USER SETTINGS ----------------------------

# Folder containing the input ARGOS CSV and where output files will be saved.
# Use "." to use the current working directory, or provide a path on your computer.
input_file <- file.path("data", "argos_tracks.csv")
output_file <- file.path("output", "argos_tracks_with_SIC.csv")
grid_file <- file.path("data", "LongitudeLatitudeGrid-s3125-Neumayer.hdf")

# ARGOS table column names (edit these if your CSV uses different names).
date_col <- "Date"
longitude_col <- "Lon"
latitude_col <- "Lat"
id_col <- "ID"  # used for reporting; set to NULL if there is no ID column

# Maximum distance allowed between a tracking location and the nearest grid cell.
max_distance_m <- 20000

# Search margin around the daily tracking locations, in degrees.
# Increase this if your locations cover a very large area or cross the dateline.
search_margin_deg <- 0.5

# Downloaded daily grids are stored here temporarily and reused if already present.
cache_dir <- file.path("data", "sea_ice_cache")

# Set TRUE to create monthly diagnostic plots in output_file's folder.
make_monthly_plots <- TRUE

# ----------------------------------------------------------------------

required_packages <- c("ncdf4", "geosphere")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the required packages first: ",
    paste(missing_packages, collapse = ", "),
    "\nRun: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

# Create output/cache folders when needed.
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input ARGOS CSV was not found: ", input_file)
}
if (!file.exists(grid_file)) {
  stop(
    "Static longitude/latitude grid was not found: ", grid_file,
    "\nDownload the matching University of Bremen grid file and set `grid_file`."
  )
}

# Read input data without changing the user's working directory.
argos <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c(date_col, longitude_col, latitude_col)
missing_columns <- setdiff(required_columns, names(argos))
if (length(missing_columns) > 0) {
  stop("Missing required column(s) in ARGOS CSV: ",
       paste(missing_columns, collapse = ", "))
}

# Parse dates robustly. Adjust this if your input date format is unusual.
argos[[date_col]] <- as.POSIXct(
  argos[[date_col]],
  tz = "UTC",
  tryFormats = c(
    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
    "%Y-%m-%d", "%d/%m/%Y %H:%M:%S",
    "%d/%m/%Y %H:%M", "%d/%m/%Y"
  )
)

if (all(is.na(argos[[date_col]]))) {
  stop(
    "No dates could be parsed in column `", date_col,
    "`. Check the date format in the input CSV."
  )
}

argos[[longitude_col]] <- suppressWarnings(as.numeric(argos[[longitude_col]]))
argos[[latitude_col]] <- suppressWarnings(as.numeric(argos[[latitude_col]]))

# Convert longitudes from 0–360 to -180–180 when necessary.
convert_lon360_to_180 <- function(lon) {
  ifelse(!is.na(lon) & lon > 180, lon - 360, lon)
}
argos[[longitude_col]] <- convert_lon360_to_180(argos[[longitude_col]])

# Initialize output column. Existing SIC values are intentionally overwritten.
argos$SI <- NA_real_

# Read the static grid once.
grid_nc <- ncdf4::nc_open(grid_file)

grid_vars <- names(grid_nc$var)
if (!("Longitudes" %in% grid_vars) || !("Latitudes" %in% grid_vars)) {
  stop(
    "The grid file must contain variables named `Longitudes` and `Latitudes`."
  )
}
grid_lon <- as.vector(ncdf4::ncvar_get(grid_nc, "Longitudes"))
grid_lat <- as.vector(ncdf4::ncvar_get(grid_nc, "Latitudes"))
ncdf4::nc_close(grid_nc)

if (length(grid_lon) != length(grid_lat)) {
  stop("Longitude and latitude grid variables have different lengths.")
}
grid_lon <- convert_lon360_to_180(grid_lon)

# Download a daily grid. Try the v5 filename first, then v5.4.
download_daily_grid <- function(date, cache_dir) {
  date_chr <- format(date, "%Y%m%d", tz = "UTC")
  year_chr <- format(date, "%Y", tz = "UTC")
  month_chr <- tolower(month.abb[as.integer(format(date, "%m", tz = "UTC"))])
  day_chr <- format(date, "%d", tz = "UTC")

  versions <- c("v5", "v5.4")
  last_error <- NULL

  for (version in versions) {
    filename <- sprintf("asi-AMSR2-s3125-%s-%s.hdf", date_chr, version)
    local_path <- file.path(cache_dir, filename)

    # Reuse a cached file if it exists and is non-empty.
    if (file.exists(local_path) && file.info(local_path)$size > 0) {
      return(local_path)
    }

    url <- sprintf(
      paste0(
        "https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/",
        "s3125/%s/%s/Neumayer/%s"
      ),
      year_chr, month_chr, filename
    )

    # Download to a temporary file so an incomplete download is not cached.
    temp_path <- tempfile(tmpdir = cache_dir, fileext = ".download")
    status <- tryCatch(
      utils::download.file(
        url, destfile = temp_path, mode = "wb", quiet = TRUE,
        method = "libcurl"
      ),
      error = function(e) {
        last_error <<- conditionMessage(e)
        return(1L)
      }
    )

    valid_download <- identical(as.integer(status), 0L) &&
      file.exists(temp_path) && file.info(temp_path)$size > 0

    if (valid_download) {
      if (file.rename(temp_path, local_path)) {
        return(local_path)
      }
      # Fallback if rename across filesystems fails.
      copied <- file.copy(temp_path, local_path, overwrite = TRUE)
      unlink(temp_path)
      if (copied && file.exists(local_path) && file.info(local_path)$size > 0) {
        return(local_path)
      }
    } else {
      unlink(temp_path)
    }
  }

  message(
    "Could not download SIC grid for ", date_chr,
    if (!is.null(last_error)) paste0(" (", last_error, ")") else ""
  )
  return(NULL)
}

# Find nearest grid cell for each point using a spherical nearest-neighbour
# search. This avoids creating a full point-by-grid distance matrix.
# The grid is first restricted to a geographic bounding box, then exact
# great-circle distances are calculated for the remaining candidates.
extract_sic_for_day <- function(points, ice_lon, ice_lat, ice_values,
                                margin_deg, max_distance_m) {
  result <- rep(NA_real_, nrow(points))

  good_points <- is.finite(points$lon) & is.finite(points$lat)
  good_grid <- is.finite(ice_lon) & is.finite(ice_lat) & is.finite(ice_values)
  if (!any(good_points) || !any(good_grid)) return(result)

  grid <- data.frame(
    lon = ice_lon[good_grid],
    lat = ice_lat[good_grid],
    ice = ice_values[good_grid]
  )

  # A common bounding box works well for regional tracks. For tracks close
  # to the dateline, fall back to the complete valid grid.
  lon_range <- range(points$lon[good_points], na.rm = TRUE)
  lat_range <- range(points$lat[good_points], na.rm = TRUE)
  crosses_dateline <- diff(lon_range) > 180

  if (!crosses_dateline) {
    candidates <- grid[
      grid$lon >= lon_range[1] - margin_deg &
      grid$lon <= lon_range[2] + margin_deg &
      grid$lat >= max(-90, lat_range[1] - margin_deg) &
      grid$lat <= min(90, lat_range[2] + margin_deg),
      , drop = FALSE
    ]
  } else {
    candidates <- grid[
      grid$lat >= max(-90, lat_range[1] - margin_deg) &
      grid$lat <= min(90, lat_range[2] + margin_deg),
      , drop = FALSE
    ]
  }

  if (nrow(candidates) == 0) {
    warning("No grid cells in the search window for this date.")
    return(result)
  }

  for (i in which(good_points)) {
    # Geodesic distance from one ARGOS point to all candidate grid cells.
    distances <- geosphere::distHaversine(
      matrix(c(points$lon[i], points$lat[i]), nrow = 1),
      as.matrix(candidates[, c("lon", "lat")])
    )
    nearest <- which.min(distances)

    if (length(nearest) == 1 && is.finite(distances[nearest]) &&
        distances[nearest] <= max_distance_m) {
      result[i] <- candidates$ice[nearest]
    }
  }
  result
}

# Process each unique calendar date.
valid_date <- !is.na(argos[[date_col]])
days_concerned <- sort(unique(as.Date(argos[[date_col]][valid_date], tz = "UTC")))

if (length(days_concerned) == 0) {
  stop("No valid dates are available to process.")
}

for (i in seq_along(days_concerned)) {
  day <- days_concerned[i]
  message(sprintf("Processing %s (%d/%d; %.1f%%)",
                  day, i, length(days_concerned),
                  100 * i / length(days_concerned)))

  row_idx <- which(valid_date &
                     as.Date(argos[[date_col]], tz = "UTC") == day)
  points <- data.frame(
    lon = argos[[longitude_col]][row_idx],
    lat = argos[[latitude_col]][row_idx]
  )
  usable_points <- is.finite(points$lon) & is.finite(points$lat)

  if (!any(usable_points)) {
    message("  No valid coordinates for this date; skipping.")
    next
  }

  daily_file <- download_daily_grid(day, cache_dir)
  if (is.null(daily_file)) {
    message("  SIC values remain NA for this date.")
    next
  }

  ice_nc <- tryCatch(
    ncdf4::nc_open(daily_file),
    error = function(e) {
      message("  Could not open daily grid: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ice_nc)) next

  ice_values <- tryCatch(
    ncdf4::ncvar_get(ice_nc, "ASI Ice Concentration"),
    error = function(e) {
      message("  SIC variable could not be read: ", conditionMessage(e))
      NULL
    },
    finally = ncdf4::nc_close(ice_nc)
  )
  if (is.null(ice_values)) next

  ice_values <- as.vector(ice_values)

  if (length(ice_values) != length(grid_lon)) {
    warning(
      "Daily SIC grid and static coordinate grid have different lengths on ",
      day, "; skipping this date."
    )
    next
  }

  # The source grid may encode missing data or special values. Convert
  # non-finite values to NA; inspect source documentation if special codes
  # appear as ordinary numeric values in your downloaded files.
  ice_values[!is.finite(ice_values)] <- NA_real_

  sic <- extract_sic_for_day(
    points = points,
    ice_lon = grid_lon,
    ice_lat = grid_lat,
    ice_values = ice_values,
    margin_deg = search_margin_deg,
    max_distance_m = max_distance_m
  )
  argos$SI[row_idx] <- sic
}

# Save the result.
utils::write.csv(argos, output_file, row.names = FALSE, na = "")
message("Saved output to: ", normalizePath(output_file, winslash = "/", mustWork = FALSE))
message("SIC extraction summary:")
print(summary(argos$SI))

# Optional monthly diagnostic plots.
if (isTRUE(make_monthly_plots)) {
  valid_plot <- !is.na(argos$SI) &
    is.finite(argos[[longitude_col]]) &
    is.finite(argos[[latitude_col]]) &
    !is.na(argos[[date_col]])

  if (any(valid_plot)) {
    plot_data <- argos[valid_plot, , drop = FALSE]
    plot_data$Month <- as.integer(format(plot_data[[date_col]], "%m", tz = "UTC"))

    plot_file <- file.path(dirname(output_file), "SIC_monthly_diagnostics.pdf")
    grDevices::pdf(plot_file, width = 8, height = 6)

    palette <- grDevices::colorRampPalette(c("blue", "cyan", "white"))(100)

    for (month_num in 1:12) {
      subset_month <- plot_data[plot_data$Month == month_num, , drop = FALSE]
      if (nrow(subset_month) == 0) {
        graphics::plot.new()
        graphics::title(main = month.name[month_num])
        graphics::text(0.5, 0.5, "No valid SIC observations")
        next
      }

      color_index <- cut(
        subset_month$SI,
        breaks = seq(0, 100, length.out = 101),
        labels = FALSE,
        include.lowest = TRUE
      )
      color_index <- pmax(1, pmin(100, color_index))

      graphics::plot(
        subset_month[[longitude_col]], subset_month[[latitude_col]],
        pch = 15, cex = 0.4, col = palette[color_index],
        xlab = "Longitude", ylab = "Latitude",
        main = month.name[month_num]
      )
    }

    grDevices::dev.off()
    message("Saved monthly diagnostic plots to: ", plot_file)
  } else {
    message("No valid SIC observations available for diagnostic plots.")
  }
}
