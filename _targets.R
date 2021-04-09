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
  ),

  # - extract geographic shapes, richness, and number of specimens
  shape_ferns_full = shape_from_comm_scaled_list(comm_scaled_list, 0.2) %>%
    # calculate redundancy
    mutate(redundancy = 1 - (richness/abundance)),
  
  # - extract community matrix
  comm_ferns_full = comm_from_comm_scaled_list(comm_scaled_list, 0.2),
  
  # - subset geographic shapes to redundancy > 0.1
  shape_ferns = filter(shape_ferns_full, redundancy > 0.1),
  
  # - subset community matrix to communities with redundancy > 0.1
  comm_ferns = filter_comm_by_redun(
    comm = comm_ferns_full,
    shape = shape_ferns,
    cutoff = 0.1
  ),
  
  # - make community matrix subset to taxa endemic to Japan
  comm_ferns_endemic = subset_comm_to_endemic(
    comm = comm_ferns,
    green_list = green_list
  ),
  
  # Read in ultrametric phylogenetic tree of all pteridophytes,
  # not including hybrids (706 taxa total)
  tar_file(japan_pterido_tree_file, "data_raw/japan_pterido_tree_dated.tre"),
  japan_pterido_tree = ape::read.tree(japan_pterido_tree_file),
  
  # - subset to only ferns
  japan_fern_tree = subset_tree(
    phy = japan_pterido_tree, 
    ppgi = ppgi),
  
  # Format trait data, subset to ferns in tree
  tar_file(raw_trait_data_file, "data_raw/JpFernLUCID_forJoel.xlsx"),
  
  fern_traits = format_traits(
    path_to_lucid_traits = raw_trait_data_file,
    taxon_id_map = green_list,
    taxon_keep_list = japan_fern_tree$tip.label),
  
  # Transform continuous traits before making distance matrix
  traits_for_dist = transform_traits(
    fern_traits,
    trans_select = c("frond_width", "stipe_length", "number_pinna_pairs"),
    scale_select = c("frond_width", "stipe_length", "number_pinna_pairs")
  ) %>%
    # Make sure all traits are scaled within -1 to 1
    assert(within_bounds(-1,1), where(is.numeric)),
  
  # Make trait distance matrix using taxon IDs as labels
  trait_distance_matrix = make_trait_dist_matrix(traits_for_dist),
  
  # Conduct randomization tests of diversity metrics ----
  
  # - make list of communities for looping
  fern_comm_list = list(comm_ferns, comm_ferns_endemic),
  
  # - make list of biodiv metrics to calculate for each community
  fern_comm_metrics = list(
    c("pd", "rpd", "pe", "rpe"),
    c("pe", "rpe")
  ),
  
  # - specify data set names so we can filter results after looping
  fern_comm_names = c("ja_ferns", "ja_ferns_endemic"),
  
  # Conduct randomization tests of phylogeny-based metrics for all ferns and
  # ferns endemic to Japan only
  tar_target(rand_test_phy,
    run_rand_analysis(
      comm = fern_comm_list,
      phy = japan_fern_tree,
      null_model = "independentswap",
      n_reps = 999,
      n_iterations = 100000,
      metrics = fern_comm_metrics,
      dataset_name = fern_comm_names) %>%
      categorize_endemism,
    pattern = map(
      comm = fern_comm_list,
      metrics = fern_comm_metrics,
      dataset_name = fern_comm_names)
  )

)