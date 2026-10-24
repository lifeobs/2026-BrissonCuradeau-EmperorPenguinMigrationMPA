
# ============================================================
# Extract mixed layer depth (MLD) at animal tracking locations
# ============================================================
#
# Description:
# Extracts daily mixed layer depth (mlotst) from annual NetCDF
# files at ARGOS tracking and dive locations.
#
# See readme file for more info on project pathing
#
#
# ============================================================


# 1. Install and load required packages -----------------------

packages <- c("lubridate", "ncdf4", "geosphere", "dplyr")

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))


# 2. User settings --------------------------------------------

# Main project folder (use "." for the current working directory)
project_dir <- "."

# Input CSV files
argos_file <- file.path(
  project_dir, "data", "input", "Argos Juv ALL + SI.csv"
)

dives_file <- file.path(
  project_dir, "data", "input", "Dives Juv ALL + SI.csv"
)

# Folder containing annual NetCDF files
environment_dir <- file.path(project_dir, "data", "environment")

# Folder for output CSV files
output_dir <- file.path(project_dir, "data", "output")

# Years to process
years_to_process <- c(2018, 2019, 2022, 2023, 2024)

# Maximum distance from the animal location to the candidate
# grid cells, in degrees of longitude/latitude
search_radius <- 0.5


# 3. Check files and prepare data -----------------------------

if (!file.exists(argos_file)) {
  stop("ARGOS CSV not found: ", argos_file)
}

if (!file.exists(dives_file)) {
  stop("Dives CSV not found: ", dives_file)
}

if (!dir.exists(environment_dir)) {
  stop("Environment folder not found: ", environment_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Argos <- read.csv(argos_file, stringsAsFactors = FALSE)
Dives <- read.csv(dives_file, stringsAsFactors = FALSE)

required_columns <- c("Date", "Lon", "Lat")

for (dataset_name in c("Argos", "Dives")) {
  dataset <- get(dataset_name)
  
  missing_columns <- setdiff(required_columns, names(dataset))
  
  if (length(missing_columns) > 0) {
    stop(
      dataset_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

# Convert dates and retain the calendar day for matching
Argos$Date <- as.Date(Argos$Date)
Dives$Date <- as.Date(Dives$Date)

if (anyNA(Argos$Date) || anyNA(Dives$Date)) {
  stop(
    "Some dates could not be parsed. Check the date format ",
    "in both input CSV files."
  )
}

Argos$MLD <- NA_real_
Dives$MLD <- NA_real_

Argos$Year <- lubridate::year(Argos$Date)
Dives$Year <- lubridate::year(Dives$Date)


# 4. Function to convert NetCDF time to dates -----------------

# Supports common CF time units such as:
#   seconds since 1970-01-01
#   hours since 1950-01-01 00:00:00
#   days since 1950-01-01

convert_netcdf_time <- function(nc, time_values) {
  
  time_units <- ncdf4::ncatt_get(nc, "time", "units")
  
  if (!isTRUE(time_units$hasatt)) {
    stop("The NetCDF time variable has no 'units' attribute.")
  }
  
  units_string <- time_units$value
  
  pattern <- paste0(
    "^\\s*(seconds|second|secs|sec|minutes|minute|mins|min|",
    "hours|hour|hrs|hr|days|day)\\s+since\\s+(.+?)\\s*$"
  )
  
  matches <- regexec(
    pattern, units_string,
    ignore.case = TRUE
  )
  
  parts <- regmatches(units_string, matches)[[1]]
  
  if (length(parts) != 3) {
    stop("Unsupported NetCDF time units: ", units_string)
  }
  
  unit <- tolower(parts[2])
  origin_string <- parts[3]
  
  origin <- as.POSIXct(
    origin_string,
    tz = "UTC",
    tryFormats = c(
      "%Y-%m-%d %H:%M:%OS",
      "%Y-%m-%dT%H:%M:%OS",
      "%Y-%m-%d",
      "%Y-%m-%d %H:%M"
    )
  )
  
  if (is.na(origin)) {
    stop("Could not parse NetCDF time origin: ", origin_string)
  }
  
  multiplier <- if (unit %in% c("seconds", "second", "secs", "sec")) {
    1
  } else if (unit %in% c("minutes", "minute", "mins", "min")) {
    60
  } else if (unit %in% c("hours", "hour", "hrs", "hr")) {
    3600
  } else if (unit %in% c("days", "day")) {
    86400
  } else {
    stop("Unsupported NetCDF time unit: ", unit)
  }
  
  as.Date(origin + time_values * multiplier, tz = "UTC")
}


# 5. Function to extract MLD for one location ------------------

extract_mld <- function(lon, lat, date, grid, max_radius) {
  
  if (is.na(lon) || is.na(lat) || is.na(date)) {
    return(NA_real_)
  }
  
  candidates <- grid[
    !is.na(grid$MLD) &
      abs(grid$lon - lon) < max_radius &
      abs(grid$lat - lat) < max_radius,
  ]
  
  if (nrow(candidates) == 0) {
    return(NA_real_)
  }
  
  distances <- geosphere::distGeo(
    p1 = cbind(lon, lat),
    p2 = as.matrix(candidates[, c("lon", "lat")])
  )
  
  candidates$distance <- distances
  
  candidates$MLD[which.min(candidates$distance)]
}


# 6. Process each year ----------------------------------------

for (yi in years_to_process) {
  
  message("Processing year: ", yi)
  
  nc_file <- file.path(environment_dir, paste0(yi, ".nc"))
  
  if (!file.exists(nc_file)) {
    warning("NetCDF file not found; skipping year: ", nc_file)
    next
  }
  
  Env <- ncdf4::nc_open(nc_file)
  
  tryCatch({
    
    MLD_array <- ncdf4::ncvar_get(Env, "mlotst")
    Lon <- ncdf4::ncvar_get(Env, "longitude")
    Lat <- ncdf4::ncvar_get(Env, "latitude")
    Time <- ncdf4::ncvar_get(Env, "time")
    
    Time_dates <- convert_netcdf_time(Env, Time)
    
    # Select animal locations from the year being processed
    argos_rows <- which(Argos$Year == yi)
    dives_rows <- which(Dives$Year == yi)
    
    target_dates <- sort(unique(c(
      Argos$Date[argos_rows],
      Dives$Date[dives_rows]
    )))
    
    target_dates <- target_dates[!is.na(target_dates)]
    
    for (target_date in target_dates) {
      
      time_index <- which(Time_dates == target_date)
      
      if (length(time_index) == 0) {
        next
      }
      
      # If multiple time steps occur on one calendar day,
      # use the first matching time step.
      time_index <- time_index[1]
      
      # Assumes the NetCDF array dimensions are longitude,
      # latitude, then time, as in the original script.
      MLD_matrix <- MLD_array[, , time_index]
      
      grid <- expand.grid(
        lon = Lon,
        lat = Lat
      )
      
      grid$MLD <- as.vector(MLD_matrix)
      
      grid <- grid[!is.na(grid$MLD), ]
      
      # Extract for ARGOS locations
      rows <- which(Argos$Date == target_date)
      
      if (length(rows) > 0) {
        for (j in rows) {
          Argos$MLD[j] <- extract_mld(
            lon = Argos$Lon[j],
            lat = Argos$Lat[j],
            date = target_date,
            grid = grid,
            max_radius = search_radius
          )
        }
      }
      
      # Extract for dive locations
      rows <- which(Dives$Date == target_date)
      
      if (length(rows) > 0) {
        for (j in rows) {
          Dives$MLD[j] <- extract_mld(
            lon = Dives$Lon[j],
            lat = Dives$Lat[j],
            date = target_date,
            grid = grid,
            max_radius = search_radius
          )
        }
      }
    }
    
  }, finally = {
    ncdf4::nc_close(Env)
  })
  
  message("Finished year: ", yi)
}


# 7. Save output files ----------------------------------------

write.csv(
  Argos,
  file.path(output_dir, "Argos Juv ALL + SI + MLD.csv"),
  row.names = FALSE
)

write.csv(
  Dives,
  file.path(output_dir, "Dives Juv ALL + SI + MLD.csv"),
  row.names = FALSE
)


# 8. Basic checks ---------------------------------------------

summary(Argos$MLD)
summary(Dives$MLD)

plot(
  Dives$MLD ~ Dives$Lat,
  xlab = "Latitude",
  ylab = "Mixed layer depth"
)

plot(
  Argos$MLD ~ Argos$Lat,
  xlab = "Latitude",
  ylab = "Mixed layer depth"
)

message("MLD extraction completed.")
