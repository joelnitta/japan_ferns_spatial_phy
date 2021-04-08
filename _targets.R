library(targets)
library(tarchetypes)

# Load packages, functions
source("R/packages.R")
source("R/functions.R")

# Specify how to resolve futures for parallel tasks
plan(multicore)

tar_plan(
  
  # Load and process data ----
  
  # Unzip data files from Ebihara and Nitta 2019
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  tar_file(ebihara_2019_zip_file, "data_raw/doi_10.5061_dryad.4362p32__v4.zip"),
  
  ebihara_2019_data = unzip_ebihara_2019(
    dryad_zip_file = ebihara_2019_zip_file, 
    exdir = "data_raw/ebihara_2019"
  ),
  
  # Split out paths for each data file
  tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  tar_file(green_list_file, ebihara_2019_data[str_detect(ebihara_2019_data, "FernGreenListV1")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  
  # Load Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = read_csv(ppgi_file) %>% modify_ppgi,
  
  # Load Fern Green List (official taxonomy + conservation status for each species)
  green_list = read_excel(green_list_file) %>% tidy_japan_names(),
  
  # Load a map of Japan
  # downloaded from https://www.gsi.go.jp/kankyochiri/gm_japan_e.html
  # on 2020-08-26
  tar_file(japan_pol_file, "data_raw/gm-jpn-all_u_2_2/polbnda_jpn.shp"),
  jpn_pol = sf::st_read(japan_pol_file),
  
  # - collapse all the political units down to just one shape for the country
  japan_shp = jpn_pol %>%
    select(geometry) %>%
    summarize(),
  
  # Load points of interest for drawing a map of Japan
  tar_file(japan_points_raw_file, "data_raw/japan_points_raw.csv"),
  japan_points_raw = read_csv(japan_points_raw_file),
  
  # Load raw occurrence data of pteridophytes in Japan, excluding hybrids (717 taxa)
  tar_file(occ_point_data_raw_file, "data_raw/JP_pterid_excl_hyb200620.xlsx"),
  occ_point_data_raw = readxl::read_excel(
    occ_point_data_raw_file,
    col_types = c("text", "numeric", "numeric", "text", "text", "text", "text"),
    col_names = c("species", "longitude", "latitude", "date", "tns_barcode", "herbarium_code", "taxon_id"),
    skip = 1),
  
  # - standardize names to Green List
  occ_point_data = rename_taxa(occ_point_data_raw, green_list) %>%
    # check for missing data
    assert(not_na, longitude, latitude, taxon),
  
  # - subset to just ferns (674 taxa)
  occ_point_data_ferns_unfiltered = subset_to_ferns(occ_point_data, ppgi),
  
  # - filter out duplicates, restrict to only points in second-degree mesh
  # Shape file downloaded from http://gis.biodic.go.jp/
  # http://gis.biodic.go.jp/BiodicWebGIS/Questionnaires?kind=mesh2&filename=mesh2.zip
  tar_file(mesh2_file, "data_raw/mesh2/mesh2.shp"),
  occ_point_data_ferns = filter_occ_points(
    occ_point_data = occ_point_data_ferns_unfiltered,
    shape_file = mesh2_file),
  
  # Calculate richness, abundance, and redundancy at four scales: 
  # 0.1, 0.2, 0.3, and 0.4 degrees
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