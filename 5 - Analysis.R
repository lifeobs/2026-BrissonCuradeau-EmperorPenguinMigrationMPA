
# ============================================================
# ANIMAL TRACKING ANALYSES
# ============================================================
#
# Analyses:
#   1. Latitude by month and colony
#      - Compare mixed-effects models using AIC
#      - Estimate differences among colonies by month
#
#   2. Fishing effort
#      - Fit a negative-binomial generalized additive model
#      - Model seasonal patterns by colony
#
#   3. Fishing effort and MPA status
#      - Calculate the percentage of locations with positive
#        fishing effort per individual
#      - Fit a GAM to fishing effort by colony and MPA status
#
#   4. Distance from colony
#      - Calculate maximum distance travelled by each individual
#      - Summarize distances overall and by colony
#
#
# ============================================================


# 1. INSTALL AND LOAD PACKAGES -------------------------------

packages <- c(
  "lubridate",
  "lme4",
  "emmeans",
  "mgcv",
  "geosphere",
  "dplyr"
)

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(
  lapply(packages, library, character.only = TRUE)
)


# 2. USER SETTINGS -------------------------------------------

# Main project folder
# "." means the current working directory
project_dir <- "."

# Input files
tracking_file <- file.path(
  project_dir,
  "data", "input",
  "Argos Juv ALL + SI + MLD + MPAs.csv"
)

fishing_file <- file.path(
  project_dir,
  "data", "input",
  "Argos Juv ALL + SI + MLD + MPAs + Fish.csv"
)

# Output directory
output_dir <- file.path(
  project_dir,
  "data", "output"
)

# Years excluded from fishing and MPA analyses
excluded_fishing_years <- c(2018, 2020, 2025)

# Additional exclusion:
# exclude Atka colony locations from 2024 in the MPA analysis
exclude_atka_2024 <- TRUE

# Reference colony coordinates, in decimal degrees
# Longitude first, latitude second
colony_coordinates <- data.frame(
  Colony = c("Atka", "Jason", "Snow"),
  ref_lon = c(-8.133333, -56.587429, -52.417753),
  ref_lat = c(-69.800000, -66.025180, -64.840170)
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 3. HELPER FUNCTIONS ----------------------------------------

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
}

check_columns <- function(data, columns, dataset_name) {
  missing <- setdiff(columns, names(data))
  
  if (length(missing) > 0) {
    stop(
      dataset_name, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

save_model_summary <- function(model, filename) {
  capture.output(
    summary(model),
    file = file.path(output_dir, filename)
  )
}

# Parse dates and create month/year fields
prepare_dates <- function(data) {
  
  data$Date <- as.Date(data$Date)
  
  if (anyNA(data$Date)) {
    stop(
      "Some dates could not be parsed. ",
      "Check the Date column in the input data."
    )
  }
  
  data$Year <- lubridate::year(data$Date)
  data$Month_num <- lubridate::month(data$Date)
  
  data
}


# ============================================================
# ANALYSIS 1: LATITUDE BY MONTH AND COLONY
# ============================================================

message("\n--- Analysis 1: Latitude by month and colony ---")

check_file(tracking_file)

Argos_lat <- read.csv(
  tracking_file,
  stringsAsFactors = FALSE
)

check_columns(
  Argos_lat,
  c("Date", "Lat", "Colony", "ID", "Year"),
  "Tracking data"
)

Argos_lat <- prepare_dates(Argos_lat)

Argos_lat$Month <- factor(
  Argos_lat$Month_num,
  levels = 1:12
)

Argos_lat$Colony <- as.factor(Argos_lat$Colony)
Argos_lat$ID <- as.factor(Argos_lat$ID)
Argos_lat$Year <- as.factor(Argos_lat$Year)

# Remove incomplete records for this analysis
Argos_lat <- Argos_lat[
  complete.cases(
    Argos_lat[, c("Lat", "Month", "Colony", "ID", "Year")]
  ),
]

# Fit candidate models using maximum likelihood (REML = FALSE)
# This makes AIC comparisons between different fixed effects
# appropriate.

M1_lat <- lme4::lmer(
  Lat ~ Month + Colony + (1 | ID) + (1 | Year),
  data = Argos_lat,
  REML = FALSE
)

M2_lat <- lme4::lmer(
  Lat ~ Colony + (1 | ID) + (1 | Year),
  data = Argos_lat,
  REML = FALSE
)

M3_lat <- lme4::lmer(
  Lat ~ Month + (1 | ID) + (1 | Year),
  data = Argos_lat,
  REML = FALSE
)

M4_lat <- lme4::lmer(
  Lat ~ 1 + (1 | ID) + (1 | Year),
  data = Argos_lat,
  REML = FALSE
)

models_lat <- list(
  M1_Month_Colony = M1_lat,
  M2_Colony = M2_lat,
  M3_Month = M3_lat,
  M4_Intercept_only = M4_lat
)

# Compare models using AIC
aic_lat <- data.frame(
  Model = names(models_lat),
  AIC = sapply(models_lat, AIC)
)

aic_lat <- aic_lat[order(aic_lat$AIC), ]
aic_lat$Delta_AIC <- aic_lat$AIC - min(aic_lat$AIC)

write.csv(
  aic_lat,
  file.path(output_dir, "latitude_model_AIC.csv"),
  row.names = FALSE
)

print(aic_lat)

# Save the selected model summary
save_model_summary(
  M1_lat,
  "latitude_model_summary.txt"
)

# Estimate colony means by month and pairwise comparisons
emm_lat <- emmeans::emmeans(
  M1_lat,
  ~ Colony | Month
)

emm_lat_df <- as.data.frame(emm_lat)

write.csv(
  emm_lat_df,
  file.path(output_dir, "latitude_emmeans_by_month.csv"),
  row.names = FALSE
)

pairs_lat <- as.data.frame(
  summary(
    emmeans::pairs(emm_lat),
    infer = TRUE
  )
)

write.csv(
  pairs_lat,
  file.path(output_dir, "latitude_colony_pairwise_comparisons.csv"),
  row.names = FALSE
)

print(summary(M1_lat))
print(pairs_lat)


# ============================================================
# ANALYSIS 2: SEASONAL PATTERNS IN FISHING EFFORT
# ============================================================

message("\n--- Analysis 2: Seasonal fishing effort ---")

check_file(fishing_file)

Argos_fish <- read.csv(
  fishing_file,
  stringsAsFactors = FALSE
)

check_columns(
  Argos_fish,
  c("Date", "Year", "Month", "Colony", "ID", "Fish"),
  "Fishing data"
)

Argos_fish <- prepare_dates(Argos_fish)

# Exclude specified years
Argos_fish <- Argos_fish[
  !Argos_fish$Year %in% excluded_fishing_years,
]

Argos_fish$Colony <- as.factor(Argos_fish$Colony)
Argos_fish$ID <- as.factor(Argos_fish$ID)
Argos_fish$Year <- as.factor(Argos_fish$Year)

# Use numeric month for the seasonal smooth
Argos_fish$Month <- Argos_fish$Month_num

# Fishing effort is rounded to match the original analysis
Argos_fish$Fish <- round(as.numeric(Argos_fish$Fish))

Argos_fish <- Argos_fish[
  complete.cases(
    Argos_fish[, c("Fish", "Month", "Colony", "ID", "Year")]
  ),
]

if (any(Argos_fish$Fish < 0)) {
  stop("Fishing effort contains negative values.")
}

# A negative-binomial GAM with colony-specific seasonal patterns
M1_fish <- mgcv::gam(
  Fish ~ Colony +
    s(Month, by = Colony) +
    s(ID, bs = "re") +
    s(Year, bs = "re"),
  data = Argos_fish,
  family = mgcv::nb(),
  method = "REML"
)

print(summary(M1_fish))

save_model_summary(
  M1_fish,
  "fishing_effort_model_summary.txt"
)


# ============================================================
# ANALYSIS 3: FISHING EFFORT AND MPA STATUS
# ============================================================

message("\n--- Analysis 3: Fishing effort and MPA status ---")

Argos_mpa <- read.csv(
  fishing_file,
  stringsAsFactors = FALSE
)

check_columns(
  Argos_mpa,
  c(
    "Date", "Year", "Month", "Colony", "ID", "Fish",
    "MPA", "MPA_Est"
  ),
  "MPA analysis data"
)

Argos_mpa <- prepare_dates(Argos_mpa)

# Apply the original year exclusions
Argos_mpa <- Argos_mpa[
  !Argos_mpa$Year %in% excluded_fishing_years,
]

# Exclude Atka locations in 2024 if requested
if (exclude_atka_2024) {
  Argos_mpa <- Argos_mpa[
    !(Argos_mpa$Year == 2024 & Argos_mpa$Colony == "Atka"),
  ]
}

Argos_mpa$Colony <- as.factor(Argos_mpa$Colony)
Argos_mpa$ID <- as.factor(Argos_mpa$ID)
Argos_mpa$Year <- as.factor(Argos_mpa$Year)
Argos_mpa$Month <- Argos_mpa$Month_num

Argos_mpa$Fish <- round(as.numeric(Argos_mpa$Fish))

# Create mutually exclusive MPA categories:
# 0 = Outside any included MPA
# 1 = Proposed MPA only
# 2 = Established MPA
#
# Established status takes priority if a location intersects
# both established and proposed polygons.

Argos_mpa$MPA_status <- NA_character_

Argos_mpa$MPA_status[
  Argos_mpa$MPA == 0
] <- "Outside MPA"

Argos_mpa$MPA_status[
  Argos_mpa$MPA == 1
] <- "Proposed MPA"

Argos_mpa$MPA_status[
  Argos_mpa$MPA_Est == 1
] <- "Established MPA"

Argos_mpa$MPA_status <- factor(
  Argos_mpa$MPA_status,
  levels = c(
    "Outside MPA",
    "Proposed MPA",
    "Established MPA"
  )
)

# Calculate the percentage of locations with positive fishing
# effort for each individual within each colony and MPA status.
#
# This is a descriptive summary, not the response variable used
# in the GAM below.

occurrence_data <- Argos_mpa |>
  dplyr::filter(
    !is.na(Fish),
    !is.na(MPA_status)
  ) |>
  dplyr::group_by(Colony, MPA_status, ID) |>
  dplyr::summarise(
    n_locations = dplyr::n(),
    n_positive = sum(Fish > 0),
    percent_positive = 100 * mean(Fish > 0),
    .groups = "drop"
  )

write.csv(
  occurrence_data,
  file.path(output_dir, "fishing_occurrence_by_individual.csv"),
  row.names = FALSE
)

# Fit the original type of model to fishing effort at individual
# locations, rather than to the calculated percentages.
model_data <- Argos_mpa[
  complete.cases(
    Argos_mpa[, c(
      "Fish", "Colony", "MPA_status",
      "Month", "ID", "Year"
    )]
  ),
]

if (any(model_data$Fish < 0)) {
  stop("Fishing effort contains negative values.")
}

# Set the outside-MPA category as the reference where present
model_data$MPA_status <- stats::relevel(
  droplevels(model_data$MPA_status),
  ref = "Outside MPA"
)

M1_mpa <- mgcv::gam(
  Fish ~ Colony * MPA_status +
    s(Month, by = Colony) +
    s(ID, bs = "re") +
    s(Year, bs = "re"),
  data = model_data,
  family = mgcv::nb(),
  method = "REML"
)

print(summary(M1_mpa))

save_model_summary(
  M1_mpa,
  "fishing_effort_MPA_model_summary.txt"
)

# Save occurrence percentages summarized across individuals
# for each colony and MPA category.
occurrence_summary <- occurrence_data |>
  dplyr::group_by(Colony, MPA_status) |>
  dplyr::summarise(
    n_individuals = dplyr::n(),
    mean_percent_positive = mean(percent_positive, na.rm = TRUE),
    sd_percent_positive = sd(percent_positive, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  occurrence_summary,
  file.path(output_dir, "fishing_occurrence_summary.csv"),
  row.names = FALSE
)


# ============================================================
# ANALYSIS 4: MAXIMUM DISTANCE FROM COLONY
# ============================================================

message("\n--- Analysis 4: Maximum distance from colony ---")

# Use the full tracking dataset, without the exclusions used
# in the fishing and MPA analyses.

Argos_distance <- read.csv(
  tracking_file,
  stringsAsFactors = FALSE
)

check_columns(
  Argos_distance,
  c("Lon", "Lat", "Colony", "ID"),
  "Distance analysis data"
)

Argos_distance$Lon <- as.numeric(Argos_distance$Lon)
Argos_distance$Lat <- as.numeric(Argos_distance$Lat)

Argos_distance$Colony <- as.character(Argos_distance$Colony)
Argos_distance$ID <- as.factor(Argos_distance$ID)

# Add reference coordinates for each colony
Argos_distance <- dplyr::left_join(
  Argos_distance,
  colony_coordinates,
  by = "Colony"
)

unknown_colonies <- setdiff(
  unique(na.omit(Argos_distance$Colony)),
  colony_coordinates$Colony
)

if (length(unknown_colonies) > 0) {
  warning(
    "No reference coordinates provided for: ",
    paste(unknown_colonies, collapse = ", ")
  )
}

# Calculate distances only where animal and reference coordinates
# are available.
valid_distance <- (
  !is.na(Argos_distance$Lon) &
    !is.na(Argos_distance$Lat) &
    !is.na(Argos_distance$ref_lon) &
    !is.na(Argos_distance$ref_lat) &
    is.finite(Argos_distance$Lon) &
    is.finite(Argos_distance$Lat)
)

Argos_distance$distance_km <- NA_real_

if (any(valid_distance)) {
  
  Argos_distance$distance_km[valid_distance] <-
    geosphere::distHaversine(
      p1 = as.matrix(
        Argos_distance[valid_distance, c("Lon", "Lat")]
      ),
      p2 = as.matrix(
        Argos_distance[valid_distance, c("ref_lon", "ref_lat")]
      )
    ) / 1000
}

# Maximum distance for each individual within each colony
individual_max_distance <- Argos_distance |>
  dplyr::filter(!is.na(distance_km)) |>
  dplyr::group_by(ID, Colony) |>
  dplyr::summarise(
    max_distance_km = max(distance_km),
    .groups = "drop"
  )

write.csv(
  individual_max_distance,
  file.path(output_dir, "maximum_distance_by_individual.csv"),
  row.names = FALSE
)

# Overall summary
overall_distance_summary <- data.frame(
  n_individuals = nrow(individual_max_distance),
  mean_max_distance_km = mean(
    individual_max_distance$max_distance_km,
    na.rm = TRUE
  ),
  sd_max_distance_km = sd(
    individual_max_distance$max_distance_km,
    na.rm = TRUE
  )
)

write.csv(
  overall_distance_summary,
  file.path(output_dir, "maximum_distance_overall_summary.csv"),
  row.names = FALSE
)

print(overall_distance_summary)

# Summary by colony
distance_by_colony <- individual_max_distance |>
  dplyr::group_by(Colony) |>
  dplyr::summarise(
    n_individuals = dplyr::n(),
    mean_max_distance_km = mean(max_distance_km, na.rm = TRUE),
    sd_max_distance_km = sd(max_distance_km, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  distance_by_colony,
  file.path(output_dir, "maximum_distance_by_colony_summary.csv"),
  row.names = FALSE
)

print(distance_by_colony)


# ============================================================
# END OF ANALYSES
# ============================================================

message(
  "\nAll analyses completed. Results saved in: ",
  output_dir
)
