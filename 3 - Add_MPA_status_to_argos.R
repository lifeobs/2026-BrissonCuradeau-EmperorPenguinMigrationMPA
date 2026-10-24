# ============================================================
# Assign Marine Protected Area (MPA) status to animal locations
# ============================================================
#
# Description:
# Determines whether each animal tracking location falls
# within an established or proposed Marine Protected Area.
#
# ============================================================


# 1. Install and load required packages -----------------------

packages <- c("sf", "dplyr", "smoothr", "units")

missing_packages <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(packages, library, character.only = TRUE))


# 2. User settings --------------------------------------------

# Main project directory; "." means the current working directory
project_dir <- "."

# Animal location file
argos_file <- file.path(
  project_dir,
  "data", "input",
  "Argos Juv ALL + SI + MLD.csv"
)

# Output directory
output_dir <- file.path(project_dir, "data", "output")

# MPA shapefile paths.
# Each shapefile requires its associated files (.shp, .shx,
# .dbf, and .prj where available).

D1MPA_file <- file.path(
  project_dir, "data", "mpas",
  "D1MPA", "D1MPA v2 (2024).shp"
)

SGSI_file <- file.path(
  project_dir, "data", "mpas",
  "South Georgia and Sandwich Islands",
  "WDPA_WDOECM_Feb2026_Public_555547601_shp-polygons.shp"
)

WSMPA_phase1_file <- file.path(
  project_dir, "data", "mpas",
  "WSMPA_phase1", "all_polygons_densified_redrawn.shp"
)

WSMPA_phase2_file <- file.path(
  project_dir, "data", "mpas",
  "WSMPA_phase2", "MPA_tmp.shp"
)

established_MPA_file <- file.path(
  project_dir, "data", "mpas",
  "established", "CCAMLR_MPA_ESRI102020.shp"
)

# Threshold for filling small holes in the combined WSMPA polygon
hole_threshold_km2 <- 100

# Geographic bounding box used to crop established MPAs
# Longitude: 80 degrees W to 20 degrees W
# Latitude:   70 degrees S to 55 degrees S
crop_bbox <- c(
  xmin = -80,
  xmax = -20,
  ymin = -70,
  ymax = -55
)


# 3. Check input files ----------------------------------------

required_files <- c(
  argos_file,
  D1MPA_file,
  SGSI_file,
  WSMPA_phase1_file,
  WSMPA_phase2_file,
  established_MPA_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following files could not be found:\n",
    paste(missing_files, collapse = "\n"),
    "\n\nCheck the file paths in Section 2."
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 4. Read and prepare proposed MPAs ----------------------------

message("Reading proposed MPA shapefiles...")

# D1MPA
D1MPA <- sf::st_read(D1MPA_file, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

# South Georgia and South Sandwich Islands
SGSI <- sf::st_read(SGSI_file, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

# Weddell Sea MPA proposal: Phase 1
WSMPA1 <- sf::st_read(WSMPA_phase1_file, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

WSMPA1 <- sf::st_union(WSMPA1)

# Weddell Sea MPA proposal: Phase 2
WSMPA2 <- sf::st_read(WSMPA_phase2_file, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

WSMPA2 <- sf::st_union(WSMPA2)

# Combine both Weddell Sea proposals
WSMPA_combined <- sf::st_union(
  WSMPA1,
  WSMPA2
)

WSMPA_combined <- sf::st_make_valid(WSMPA_combined)

# Convert the geometry to an sf object
WSMPA_sf <- sf::st_sf(
  geometry = WSMPA_combined
)

# Fill small holes
WSMPA <- smoothr::fill_holes(
  WSMPA_sf,
  threshold = units::set_units(hole_threshold_km2, km^2)
)

WSMPA <- sf::st_make_valid(WSMPA)


# 5. Read and prepare established MPAs -------------------------

message("Reading established MPA shapefile...")

Established_raw <- sf::st_read(
  established_MPA_file,
  quiet = TRUE
)

# Convert to longitude/latitude
Established_raw <- sf::st_transform(
  Established_raw,
  4326
)

# Repair geometries
Established_raw <- sf::st_make_valid(Established_raw)

# Handle geometries crossing the antimeridian
Established_raw <- sf::st_wrap_dateline(
  Established_raw,
  options = c("WRAPDATELINE=YES"),
  quiet = TRUE
)

# Crop to the specified geographic area
bbox <- sf::st_bbox(
  crop_bbox,
  crs = sf::st_crs(Established_raw)
)

Established_crop <- sf::st_crop(
  Established_raw,
  bbox
)

Established_crop <- sf::st_make_valid(Established_crop)


# 6. Create combined MPA layers -------------------------------

message("Combining MPA layers...")

# Dissolve each category to simplify spatial intersection
WSMPA_union <- sf::st_union(WSMPA)
D1MPA_union <- sf::st_union(D1MPA)
SGSI_union <- sf::st_union(SGSI)
Established_union <- sf::st_union(Established_crop)

# Proposed MPAs
Proposed_MPA <- sf::st_union(
  WSMPA_union,
  D1MPA_union
)

# Established MPAs
Established_MPA <- sf::st_union(
  SGSI_union,
  Established_union
)

# All included MPAs
All_MPA <- sf::st_union(
  Proposed_MPA,
  Established_MPA
)

# Ensure geometries are valid
Proposed_MPA <- sf::st_make_valid(Proposed_MPA)
Established_MPA <- sf::st_make_valid(Established_MPA)
All_MPA <- sf::st_make_valid(All_MPA)


# 7. Read animal locations ------------------------------------

message("Reading animal tracking data...")

Argos <- read.csv(
  argos_file,
  stringsAsFactors = FALSE
)

required_columns <- c("Lon", "Lat")

missing_columns <- setdiff(
  required_columns,
  names(Argos)
)

if (length(missing_columns) > 0) {
  stop(
    "The animal location CSV is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Ensure coordinates are numeric
Argos$Lon <- as.numeric(Argos$Lon)
Argos$Lat <- as.numeric(Argos$Lat)

if (any(
  !is.na(Argos$Lon) &
  (Argos$Lon < -180 | Argos$Lon > 180)
)) {
  stop("Longitude values must be between -180 and 180.")
}

if (any(
  !is.na(Argos$Lat) &
  (Argos$Lat < -90 | Argos$Lat > 90)
)) {
  stop("Latitude values must be between -90 and 90.")
}

# Identify rows with usable coordinates
valid_coordinates <- (
  !is.na(Argos$Lon) &
    !is.na(Argos$Lat) &
    is.finite(Argos$Lon) &
    is.finite(Argos$Lat)
)


# 8. Convert locations to a spatial object ---------------------

message("Converting animal locations to spatial points...")

Argos_sf <- sf::st_as_sf(
  Argos,
  coords = c("Lon", "Lat"),
  crs = 4326,
  remove = FALSE,
  na.fail = FALSE
)

# Initialize indicators
Argos$MPA <- 0L
Argos$MPA_Est <- 0L
Argos$MPA_Prop <- 0L

# Rows without valid coordinates cannot be classified
Argos$MPA[!valid_coordinates] <- NA_integer_
Argos$MPA_Est[!valid_coordinates] <- NA_integer_
Argos$MPA_Prop[!valid_coordinates] <- NA_integer_


# 9. Assign MPA indicators ------------------------------------

message("Checking animal locations against MPAs...")

valid_rows <- which(valid_coordinates)

if (length(valid_rows) > 0) {
  
  valid_points <- Argos_sf[valid_rows, ]
  
  # Any included MPA
  inside_all <- lengths(
    sf::st_intersects(valid_points, All_MPA)
  ) > 0
  
  # Established MPAs
  inside_established <- lengths(
    sf::st_intersects(valid_points, Established_MPA)
  ) > 0
  
  # Proposed MPAs
  inside_proposed <- lengths(
    sf::st_intersects(valid_points, Proposed_MPA)
  ) > 0
  
  Argos$MPA[valid_rows] <- as.integer(inside_all)
  
  Argos$MPA_Est[valid_rows] <- as.integer(
    inside_established
  )
  
  Argos$MPA_Prop[valid_rows] <- as.integer(
    inside_proposed
  )
}


# 10. Save results --------------------------------------------

output_file <- file.path(
  output_dir,
  "Argos Juv ALL + SI + MLD + MPAs.csv"
)

write.csv(
  Argos,
  output_file,
  row.names = FALSE
)

message("MPA assignment completed.")
message("Output saved to: ", output_file)


# 11. Summarize results ---------------------------------------

message("\nMPA assignment summary:")

print(
  data.frame(
    Category = c(
      "Any MPA",
      "Established MPA",
      "Proposed MPA"
    ),
    Number_inside = c(
      sum(Argos$MPA == 1, na.rm = TRUE),
      sum(Argos$MPA_Est == 1, na.rm = TRUE),
      sum(Argos$MPA_Prop == 1, na.rm = TRUE)
    )
  )
)