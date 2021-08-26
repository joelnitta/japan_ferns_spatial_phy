# Separate targets plan to process raw data that can't be uploaded to Dryad
#
# This should be run with:
# targets::tar_make(script = "R/process_raw_data.R", store = "_raw_data_cache")

library(targets)
library(tarchetypes)

# Load packages, functions, tests
source("R/packages.R")
source("R/functions.R")
source("R/tests.R")

# Set parallel back-end
plan(callr)

# Summary: from raw occurrence data, test binning into grid cells at four scales,
# select the optimal scale, filter out poorly sampled grid cells, generate
# community data matrix (sites x species)
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
  
  # Calculate summary statistics for occurrence data, write out as CSV
  occ_point_data_summary = summarize_occ_data(
    occ_point_data_ferns_unfiltered = occ_point_data_ferns_unfiltered, 
    occ_point_data_ferns = occ_point_data_ferns),
  tar_file(
    occ_point_data_summary_out,
    write_csv_tar(occ_point_data_summary, "data/japan_ferns_occ_summary.csv")
  ),
  
  # Calculate latitudinal span by taxon, write out as CSV
  lat_span_summary = summarize_fern_lat_span(occ_point_data_ferns),
  tar_file(
    lat_span_summary_out,
    write_csv_tar(lat_span_summary, "data/japan_ferns_lat_span.csv")
  ),
  
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
  ),
  
  # After inspecting results (see SI), select 0.2 degrees as optimal scale.
  
  # Extract geographic shapes, richness, and number of specimens at 0.2 degree grid scale
  shape_ferns_full = shape_from_comm_scaled_list(comm_scaled_list, 0.2) %>%
    # calculate redundancy
    mutate(redundancy = 1 - (richness/abundance)),
  
  # Write out geographic shapes in GeoPackage format (https://www.geopackage.org/)
  tar_file(
    shape_ferns_full_out,
    st_write_tar(shape_ferns_full, "data/japan_ferns_shape_full.gpkg")
  ),
  
  # Extract community matrix from shapes
  comm_ferns_full = comm_from_comm_scaled_list(comm_scaled_list, 0.2),
  
  # Write out community matrix as CSV
  tar_file(
    comm_ferns_full_out,
    write_comm_to_csv(comm_ferns_full, "data/japan_ferns_comm_full.csv")
  ),
  
  # Subset geographic shapes to redundancy > 0.1
  shape_ferns = filter(shape_ferns_full, redundancy > 0.1),
  
  # Subset community matrix to communities with redundancy > 0.1
  comm_ferns = filter_comm_by_redun(
    comm = comm_ferns_full,
    shape = shape_ferns,
    cutoff = 0.1
  ),
  
  # Calculate redundancy across different grain sizes
  redundancy_by_res = calc_redundancy_by_res(comm_scaled_list),
  tar_file(
    redundancy_by_res_out,
    write_csv_tar(redundancy_by_res, "data/redundancy_by_res.csv")),
  
  # Assess sampling completeness with iNEXT
  inext_res = run_inext_on_ferns(occ_point_data_ferns),
  tar_file(
    inext_res_out,
    write_csv_tar(inext_res, "data/inext_results.csv")),
  
  # Clean raw lucid data (remove data in Japanese)
  tar_target(lucid_data_ja_raw, "data_raw/JpFernLucid_forJoel20200827.xlsx"),
  raw_trait_data = clean_lucid_traits(lucid_data_ja_raw),
  tar_file(
    raw_trait_data_out,
    write_csv_tar(raw_trait_data, "data/japan_ferns_traits_lucid.csv")
  ),
  
  # Load climate data downloaded from WorldClim database
  # The WorldClim data first needs to be downloaded with this:
  # raster::getData("worldclim", download = TRUE, var = "bio", res = 2.5, path = "data_raw/world_clim")
  tar_file(worldclim_dir, "data_raw/world_clim"),
  ja_climate_data = load_ja_worldclim_data(worldclim_dir),
  # Write out geographic shapes in GeoPackage format (https://www.geopackage.org/)
  tar_file(
    ja_climate_data_out,
    st_write_tar(ja_climate_data, "data/japan_climate.gpkg")
  )
  
)
