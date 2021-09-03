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
  
  # Load raw occurrence data of pteridophytes in Japan, excluding hybrids
  tar_file(occ_point_data_raw_file, "data_raw/JP_pterid_excl_hyb200620.xlsx"),
  occ_point_data_raw = readxl::read_excel(
    occ_point_data_raw_file,
    col_types = c("text", "numeric", "numeric", "text", "text", "text", "text"),
    col_names = c("species", "longitude", "latitude", "date", "tns_barcode", "herbarium_code", "taxon_id"),
    skip = 1),
  
  # Clean data and standardize names to Green List
  occ_point_data = clean_occ_point_data(occ_point_data_raw, green_list),
  
  # Load Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = load_ppgi(ebihara_2019_zip_file) %>% modify_ppgi,
  
  # Subset occurrence data to just ferns
  occ_point_data_ferns_unfiltered = subset_to_ferns(occ_point_data, ppgi),
  
  # Calculate summary statistics for occurrence data, write out as CSV
  occ_point_data_summary = summarize_occ_data(
    occ_point_data_ferns_unfiltered = occ_point_data_ferns_unfiltered, 
    occ_point_data_ferns = occ_point_data_ferns),
  tar_file(
    occ_point_data_summary_out,
    write_csv_tar(occ_point_data_summary, "data/japan_ferns_occ_summary.csv")
  ),
  
  # Filter out duplicates, restrict to only points in second-degree mesh
  # Shape file downloaded from http://gis.biodic.go.jp/
  # http://gis.biodic.go.jp/BiodicWebGIS/Questionnaires?kind=mesh2&filename=mesh2.zip
  tar_file(mesh2_zip_file, "data_raw/mesh2.zip"),
  
  japan_mesh2 = load_shape_from_zip(mesh2_zip_file, "mesh2.shp") %>% select(id = NAME),
  
  occ_point_data_ferns = filter_occ_points(
    occ_point_data = occ_point_data_ferns_unfiltered,
    mask = japan_mesh2),
  
  # get CRS (JGD2000) from japan_mesh2 for converting point data to community dataframe
  jgd2000_sp = get_sp_crs_from_sf(japan_mesh2),
  jgd2000_sf = st_crs(japan_mesh2),

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
      species = "taxon",
      crs = jgd2000_sp),
    pattern = map(scales_to_test)
  ),
  
  # After inspecting results (see SI), select 0.2 degrees as optimal scale.
  
  # Extract geographic shapes, richness, and number of specimens at 0.2 degree grid scale
  shape_ferns_full = shape_from_comm_scaled_list(comm_scaled_list, 0.2) %>%
    # calculate redundancy
    mutate(redundancy = 1 - (richness/abundance)) %>%
    # set CRS to be exactly same as JGD2000 in japan_mesh2
    st_transform(jgd2000_sf),
  
  # Write out geographic shapes in GeoPackage format (https://www.geopackage.org/)
  tar_file(
    shape_ferns_full_out,
    # Include manual time stamp so that sha is stable
    # https://github.com/r-spatial/rgeopackage
    st_write_tar(shape_ferns_full, "data/japan_ferns_shape_full.gpkg", time_stamp = as.Date("2021-09-03"))
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
    write_csv_tar(redundancy_by_res, "data/japan_ferns_redundancy_by_res.csv")),
  
  # Assess sampling completeness with iNEXT
  inext_res = run_inext_on_ferns(occ_point_data_ferns),
  tar_file(
    inext_res_out,
    write_csv_tar(inext_res, "data/japan_ferns_inext_results.csv")),
  
  # Clean raw lucid data (remove data in Japanese, correct names)
  tar_target(lucid_data_raw_file, "data_raw/JpFernLucid_forJoel20200827.xlsx"),
  raw_lucid_data = read_excel(lucid_data_raw_file, skip = 1), # skip top row, which is in Japanese
  
  # - raw data has some faulty names (synonyms, non-native taxa). Load data to clean these.
  tar_target(lucid_name_correction_file, "data_raw/lucid_taxa_name_correction.csv"),
  lucid_name_correction = read_csv(lucid_name_correction_file),
  
  lucid_data_ferns = clean_lucid_traits(
    raw_lucid_data = raw_lucid_data, lucid_name_correction = lucid_name_correction, 
    green_list = green_list, ppgi = ppgi),
  
  tar_file(
    lucid_data_ferns_out,
    write_csv_tar(lucid_data_ferns, "data/japan_ferns_traits_lucid.csv")
  ),
  
  # Load climate data downloaded from WorldClim database
  # The WorldClim data first needs to be downloaded with this:
  # raster::getData("worldclim", download = TRUE, var = "bio", res = 2.5, path = "data_raw/world_clim")
  tar_file(worldclim_dir, "data_raw/world_clim"),
  ja_climate_data = load_ja_worldclim_data(worldclim_dir, jgd2000_sf),
  # Write out geographic shapes in GeoPackage format (https://www.geopackage.org/)
  tar_file(
    ja_climate_data_out,
    st_write_tar(ja_climate_data, "data/japan_climate.gpkg", time_stamp = as.Date("2021-08-26"))
  ),

  ## Read in protected areas, assign protection levels following Kusumoto et al. 2017:
  # - high: no human activities allowed
  # - medium: permission required for economic activities
  # - low: protected area, but none of the above restrictions
  tar_file(protected_areas_zip_file, "data_raw/map17.zip"),
  protected_areas_other = load_protected_areas(protected_areas_zip_file),

  tar_file(protected_areas_forest_folder, "data_raw/forest_area_zip_files"),
  protected_areas_forest = load_protected_areas_forest(protected_areas_forest_folder) %>%
    sf::st_transform(jgd2000_sf),

  # Combine protected areas (high and medium only), write out
  protected_areas = combine_pa(protected_areas_other, protected_areas_forest),
  tar_file(
    protected_areas_out,
    st_write_tar(protected_areas, "data/japan_protected_areas.gpkg", time_stamp = as.Date("2021-08-30"))
  ),

  ## Read in deer distribution maps, write out
  tar_file(deer_range_zip_file, "data_raw/map14-1.zip"),
  japan_deer_range = load_deer_range(deer_range_zip_file),
  tar_file(
    deer_areas_out,
    st_write_tar(japan_deer_range, "data/japan_deer_range.gpkg", time_stamp = as.Date("2021-09-01"))
  ),

  ## Read in political map of Japan, set CRS, write out
  # downloaded from https://www.gsi.go.jp/kankyochiri/gm_japan_e.html
  tar_file(japan_pol_zip_file, "data_raw/gm-jpn-all_u_2_2.zip"),

  japan_shp = japan_pol_zip_file %>%
    load_shape_from_zip("polbnda_jpn.shp") %>%
    # Collapse all the political units down to just one shape for the country
    select(geometry) %>%
    summarize() %>%
    # set CRS to JGD2000
    sf::st_transform(jgd2000_sf),

  tar_file(
    japan_shp_out,
    st_write_tar(japan_shp, "data/japan_map.gpkg", time_stamp = as.Date("2021-09-03"))
  )
  
)
