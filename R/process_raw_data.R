# Separate targets plan to process raw data that can't be uploaded to Dryad
#
# This should be run with:
# targets::tar_make(
#   script = "R/process_raw_data.R",
#   store = "_raw_data_cache"
# )

library(targets)
library(tarchetypes)

# Load packages, functions, tests
source("R/packages.R")
source("R/functions.R")
source("R/tests.R")

tar_plan(
  
  # Load "green list" from Ebihara and Nitta 2019
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  tar_file(ebihara_2019_zip_file, "data_raw/doi_10.5061_dryad.4362p32__v4.zip"),
  green_list = load_green_list(ebihara_2019_zip_file),
  
  # Load raw occurrence data of pteridophytes in Japan, excluding hybrids (717 taxa)
  tar_file(occ_point_data_raw_file, "data_raw/JP_pterid_excl_hyb200620.xlsx"),
  occ_point_data_raw = readxl::read_excel(
    occ_point_data_raw_file,
    col_types = c("text", "numeric", "numeric", "text", "text", "text", "text"),
    col_names = c("species", "longitude", "latitude", "date", "tns_barcode", "herbarium_code", "taxon_id"),
    skip = 1),
  
  # Standardize names to Green List
  occ_point_data = rename_taxa(occ_point_data_raw, green_list) %>%
    # check for missing data
    assert(not_na, longitude, latitude, taxon),
  
  # Load Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = load_ppgi(ebihara_2019_zip_file) %>% modify_ppgi,
  
  # Subset to just ferns (674 taxa)
  occ_point_data_ferns_unfiltered = subset_to_ferns(occ_point_data, ppgi),
  
  # Filter out duplicates, restrict to only points in second-degree mesh
  # Shape file downloaded from http://gis.biodic.go.jp/
  # http://gis.biodic.go.jp/BiodicWebGIS/Questionnaires?kind=mesh2&filename=mesh2.zip
  tar_file(mesh2_file, "data_raw/mesh2/mesh2.shp"),
  occ_point_data_ferns = filter_occ_points(
    occ_point_data = occ_point_data_ferns_unfiltered,
    shape_file = mesh2_file),
  
  # Calculate richness, abundance, and redundancy at four scales: 
  # 0.1, 0.2, 0.3, and 0.4 degree grid squares
  scales_to_test = c(0.1, 0.2, 0.3, 0.4),
  tar_target(
    comm_scaled_list,
    comm_from_points(
      species_coods = occ_point_data_ferns,
      resol = scales_to_test,
      lon = "longitude",
      lat = "latitude",
      species = "taxon"),
    pattern = map(scales_to_test)
  )
  
)
