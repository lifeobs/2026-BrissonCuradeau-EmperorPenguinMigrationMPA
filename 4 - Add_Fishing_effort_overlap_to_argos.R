
# ============================================================
# Calculate fishing effort around animal tracking locations
# ============================================================
#
# Description:
# Assigns cumulative fishing effort to each animal location
# using fishing data from the current month and the previous
# 11 months.
#
# Fishing effort is summed within a specified radius around
# each location.
#
# ============================================================


# 1. Install and load required packages -----------------------

packages <- c("lubridate", "geosphere")

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))


# 2. User settings --------------------------------------------

# Main project folder; "." means the current working directory
project_dir <- "."

# Input files
argos_file <- file.path(
  project_dir,
  "data", "input",
  "Argos Juv ALL + SI + MLD + MPAs.csv"
)

fishing_file <- file.path(
  project_dir,
  "data", "environment",
  "Fishing Monthly.csv"
)

# Output folder
output_dir <- file.path(project_dir, "data", "output")

# Spatial search radius, in kilometres
search_radius_km <- 100

# Number of months to include, including the current month
months_to_include <- 12

# Years excluded from processing
excluded_years <- 2018


# 3. Check input files ----------------------------------------

if (!file.exists(argos_file)) {
  stop("Animal tracking file not found: ", argos_file)
}

if (!file.exists(fishing_file)) {
  stop("Fishing data file not found: ", fishing_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 4. Import data ----------------------------------------------

message("Reading input data...")

Argos <- read.csv(
  argos_file,
  stringsAsFactors = FALSE
)

Fish <- read.csv(
  fishing_file,
  stringsAsFactors = FALSE
)


# 5. Check required columns -----------------------------------

required_argos <- c("Date", "Year", "Lon", "Lat")

required_fish <- c("Time", "lon", "lat", "Effort")

missing_argos <- setdiff(
  required_argos,
  names(Argos)
)

missing_fish <- setdiff(
  required_fish,
  names(Fish)
)

if (length(missing_argos) > 0) {
  stop(
    "Animal tracking data is missing columns: ",
    paste(missing_argos, collapse = ", ")
  )
}

if (length(missing_fish) > 0) {
  stop(
    "Fishing data is missing columns: ",
    paste(missing_fish, collapse = ", ")
  )
}


# 6. Prepare dates and coordinates ----------------------------

# Convert animal tracking dates
Argos$Date <- as.Date(Argos$Date)

if (anyNA(Argos$Date)) {
  stop(
    "Some animal tracking dates could not be parsed. ",
    "Check the Date column."
  )
}

# Derive year from Date to avoid inconsistencies
Argos$Year <- lubridate::year(Argos$Date)

# Convert fishing month labels to YYYY-MM format
Fish$Time <- as.character(Fish$Time)

# Convert coordinates and effort to numeric
Argos$Lon <- as.numeric(Argos$Lon)
Argos$Lat <- as.numeric(Argos$Lat)

Fish$lon <- as.numeric(Fish$lon)
Fish$lat <- as.numeric(Fish$lat)
Fish$Effort <- as.numeric(Fish$Effort)

# Validate fishing month format
valid_month_format <- grepl(
  "^\\d{4}-\\d{2}$",
  Fish$Time
)

if (any(!valid_month_format & !is.na(Fish$Time))) {
  stop(
    "The Fishing data Time column must use YYYY-MM format, ",
    "for example 2024-03."
  )
}

# Check whether month labels represent valid calendar months
parsed_months <- as.Date(
  paste0(Fish$Time, "-01"),
  format = "%Y-%m-%d"
)

if (any(is.na(parsed_months) & !is.na(Fish$Time))) {
  stop("Some fishing month values are invalid.")
}

# Remove fishing records without usable coordinates or effort
valid_fish <- (
  !is.na(Fish$lon) &
    !is.na(Fish$lat) &
    !is.na(Fish$Effort) &
    is.finite(Fish$lon) &
    is.finite(Fish$lat) &
    is.finite(Fish$Effort) &
    Fish$lon >= -180 &
    Fish$lon <= 180 &
    Fish$lat >= -90 &
    Fish$lat <= 90
)

Fish <- Fish[valid_fish, ]


# 7. Initialize fishing effort column -------------------------

Argos$Fish <- NA_real_

# Identify locations to process
target_rows <- which(
  !Argos$Year %in% excluded_years &
    !is.na(Argos$Lon) &
    !is.na(Argos$Lat) &
    is.finite(Argos$Lon) &
    is.finite(Argos$Lat) &
    Argos$Lon >= -180 &
    Argos$Lon <= 180 &
    Argos$Lat >= -90 &
    Argos$Lat <= 90
)

# Set effort to zero for valid locations if no fishing records
# are found within the search window
Argos$Fish[target_rows] <- 0


# 8. Calculate fishing effort ---------------------------------

message(
  "Calculating fishing effort for ",
  length(target_rows),
  " locations..."
)

# Convert the search radius to metres
search_radius_m <- search_radius_km * 1000

for (counter in seq_along(target_rows)) {
  
  i <- target_rows[counter]
  
  # Current month of the animal location
  target_date <- Argos$Date[i]
  
  start_month <- as.Date(
    format(target_date, "%Y-%m-01")
  )
  
  # Current month and preceding months
  months_seq <- seq.Date(
    from = start_month,
    by = "-1 month",
    length.out = months_to_include
  )
  
  months_vector <- format(
    months_seq,
    "%Y-%m"
  )
  
  # Select fishing records from the relevant months
  subFish <- Fish[
    Fish$Time %in% months_vector,
  ]
  
  if (nrow(subFish) == 0) {
    next
  }
  
  # Calculate geographic distances in metres
  # Coordinate order: longitude, latitude
  target_coords <- c(
    Argos$Lon[i],
    Argos$Lat[i]
  )
  
  fishing_coords <- as.matrix(
    subFish[, c("lon", "lat")]
  )
  
  distances <- geosphere::distHaversine(
    p1 = fishing_coords,
    p2 = target_coords
  )
  
  # Keep fishing records within the search radius
  within_radius <- distances <= search_radius_m
  
  # Sum effort from qualifying fishing records
  Argos$Fish[i] <- sum(
    subFish$Effort[within_radius],
    na.rm = TRUE
  )
  
  # Progress message every 100 locations
  if (counter %% 100 == 0) {
    message(
      "Processed ",
      counter,
      " of ",
      length(target_rows),
      " locations."
    )
  }
}


# 9. Save results ---------------------------------------------

output_file <- file.path(
  output_dir,
  "Argos Juv ALL + SI + MLD + MPAs + Fish.csv"
)

write.csv(
  Argos,
  output_file,
  row.names = FALSE
)

message("Fishing effort calculation completed.")
message("Output saved to: ", output_file)


# 10. Summarize results ---------------------------------------

message("\nSummary of fishing effort:")

print(summary(Argos$Fish))

message(
  "Locations processed: ",
  length(target_rows)
)

message(
  "Locations with missing fishing effort: ",
  sum(is.na(Argos$Fish))
)
