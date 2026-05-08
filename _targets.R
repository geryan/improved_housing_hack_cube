
library(targets)
library(geotargets)
library(targets.utils)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "dplyr",
    "terra",
    "malariaAtlas",
    "geotargets"
  )
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_terra_rast(
    name = housing_2000,
    command = getRaster(
      dataset_id = "Explorer__2019_Nature_Africa_Housing_2000",
      file_path = "data/raster"
    )
  ),

  tar_terra_rast(
    name = housing_2015,
    command = getRaster(
      dataset_id = "Explorer__2019_Nature_Africa_Housing_2015",
      file_path = "data/raster"
    )
  ),


  tar_target(
    housing_model_fitted,
    fit_housing_model(
      housing_2000,
      housing_2015
    )
  ),

  tar_target(
    name = pointless_end_target,
    command = NULL
  )
)
