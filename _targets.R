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
  tar_file(esm1_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ESM1")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  # tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  
  # Load Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = read_csv(ppgi_file) %>% modify_ppgi,
  
  # Load Fern Green List (official taxonomy + conservation status for each species)
  green_list = read_excel(green_list_file) %>% tidy_japan_names(),
  
  # Load reproductive mode data (Ebihara et al 2019 ESM1)
  repro_data_raw = read_csv(esm1_file),
  
  # Clean up reproductive mode data
  repro_data = process_repro_data(repro_data_raw, green_list),
  
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
      n_reps = 20,
      n_iterations = 100,
      metrics = fern_comm_metrics,
      dataset_name = fern_comm_names) %>%
      categorize_endemism,
    pattern = map(
      comm = fern_comm_list,
      metrics = fern_comm_metrics,
      dataset_name = fern_comm_names)
  ),
  
  # Separate out randomization results by dataset
  rand_test_phy_ferns = filter(rand_test_phy, dataset == "ja_ferns") %>% 
    select(-dataset),
  
  rand_test_phy_ferns_endemic = filter(rand_test_phy, dataset == "ja_ferns_endemic") %>% 
    select(-dataset),
  
  # Conduct randomization tests of trait-based metrics for all ferns
  rand_test_traits_ferns = run_rand_analysis(
    comm = comm_ferns,
    null_model = "independentswap",
    n_reps = 20,
    n_iterations = 100,
    metrics = c("fd", "rfd"),
    trait_distances = trait_distance_matrix,
    dataset_name = "ja_ferns") %>%
    select(-dataset), # don't need to include dataset name in the result
  
  # Analyze bioregions ----
  
  # - Assess optimal K-value for clustering by taxonomy
  k_taxonomy = find_k_taxonomy(comm_ferns),
  
  # - Assess optimal K-value for clustering by phylogeny
  k_phylogeny = find_k_phylogeny(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
  ),
  
  # - Cluster by taxonomy
  regions_taxonomy = cluster_taxonomic_regions(
    comm_df = comm_ferns,
    k = k_taxonomy[["optimal"]][["k"]]
  ),
  
  # - Cluster by phylogeny
  regions_phylogeny = cluster_phylo_regions(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
    k = k_phylogeny[["optimal"]][["k"]]
  ),
  
  # Traits ----
  
  # Analyze phylogenetic signal
  
  # - specify continuous traits
  cont_traits = c("frond_width", "stipe_length", "number_pinna_pairs"),
  
  # - analyze phy signal in continuous traits with K and lambda
  phy_sig_results = map_df(
    cont_traits,
    ~analyze_cont_phylosig(selected_trait = ., traits = fern_traits, phy = japan_fern_tree)
  ),
  
  # - analyze phy signal in binary traits with D
  fern_traits_binary = select(fern_traits, -any_of(cont_traits)),
  
  binary_sig_results = analyze_binary_phylosig(fern_traits_binary, japan_fern_tree),
  
  # Summarize traits
  traits_summary = make_trait_summary(fern_traits),
  
  # Reproductive mode ----
  # Calculate % apomictic species in each fern community
  percent_apo = calc_perc_apo(comm_ferns, repro_data),
  
  # Combine results ----
  
  # Combine spatial data, alpha diversity, and regions, add significance and endemism types

  # - All ferns
  biodiv_ferns_spatial =
    shape_ferns %>%
    left_join(rand_test_phy_ferns, by = c(grids = "site")) %>%
    left_join(rand_test_traits_ferns, by = c(grids = "site")) %>%
    left_join(percent_apo, by = c(grids = "site")) %>%
    left_join(regions_taxonomy %>% rename(taxonomic_cluster = cluster), by = "grids") %>%
    left_join(regions_phylogeny %>% rename(phylo_cluster = cluster), by = "grids") %>%
    categorize_endemism() %>%
    categorize_signif(),
  
  # - Japan endemics only
  biodiv_ferns_endemic_spatial =
    shape_ferns %>%
    left_join(rand_test_phy_ferns_endemic, by = c(grids = "site")) %>%
    categorize_endemism(),
  
  # Conservation analysis ----
  
  ## Read in protected areas (7 separate shape files corresponding to different kinds of areas)
  # Assign protection levels following Kusamoto et al. 2017
  # - high: no human activities allowed
  # - medium: permission required for economic activities
  # - low: protected area, but none of the above restrictions
  
  # 1: wilderness
  protected_1 = sf::st_read("data_raw/map17/原生自然環境保全地域_国指定自然環境保全地域.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝原生自然環境保全地域
        ZONE == 2 ~ "high", # 2＝特別地区
        ZONE == 3 ~ "high", # 3＝海中特別地区
        ZONE == 4 ~ "low" # 4＝普通地区
      )
    ),
  
  # 2: quasi-national parks
  protected_2 = sf::st_read("data_raw/map17/国定公園.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝特別保護地区
        ZONE == 20 ~ "high", # 20＝特別地域
        ZONE == 21 ~ "medium", # 21＝第1種特別地域
        ZONE == 22 ~ "medium", # 22＝第2種特別地域
        ZONE == 23 ~ "medium", # 23＝第3種特別地域
        ZONE == 3 ~ "low", # 3＝普通地区
        ZONE == 5 ~ "marine" #5＝海域公園地区
      )
    ),
  
  # 3: national wildlife protection areas
  protected_3 = sf::st_read("data_raw/map17/国指定鳥獣保護区.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "low", # 1＝鳥獣保護区（特別保護地区以外
        ZONE == 2 ~ "medium", # 2＝特別保護地区
      )
    ),
  
  # 4: national parks
  protected_4 = sf::st_read("data_raw/map17/国立公園.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝特別保護地区
        ZONE == 20 ~ "high", # 20＝特別地域
        ZONE == 21 ~ "medium", # 21＝第1種特別地域
        ZONE == 22 ~ "medium", # 22＝第2種特別地域
        ZONE == 23 ~ "medium", # 23＝第3種特別地域
        ZONE == 3 ~ "low", # 3＝普通地区
        ZONE == 5 ~ "marine" #5＝海域公園地区
      )
    ) %>%
    # Remove protected area in inland sea (marine)
    filter(NAME != "瀬戸内海"),
  
  # 5: prefectural wildlife protection areas
  protected_5 = sf::st_read("data_raw/map17/都道府県指定鳥獣保護区.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "low", # 1＝鳥獣保護区（特別保護地区以外
        ZONE == 2 ~ "medium", # 2＝特別保護地区
        ZONE == 3 ~ "low" # not specified, but assume no other special protection
      )
    ),
  
  # 6: prefectural natural parks
  protected_6 = sf::st_read("data_raw/map17/都道府県立自然公園.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝特別保護地区
        ZONE == 20 ~ "high", # 20＝特別地域
        ZONE == 3 ~ "low" # 3＝普通地区
      )
    ),
  
  # 7: prefectural protection areas
  protected_7 = sf::st_read("data_raw/map17/都道府県自然環境保全地域.shp") %>%
    mutate(
      status = case_when(
        ZONE == 0 ~ "high", # 0＝原生自然環境保全地域
        ZONE == 2 ~ "high", # 2＝特別地区
        ZONE == 4 ~ "low" # 4＝普通地区
      )
    ),
  
  # Combine protected areas into single dataframe
  protected_areas = combine_protected_areas(
    protected_1,
    protected_2,
    protected_3,
    protected_4,
    protected_5,
    protected_6,
    protected_7
  ),
  
  # Calculate percent protection for grid cells with significantly high biodiversity
  signif_cells_protected_area = calculate_protected_area(biodiv_ferns_spatial, protected_areas, japan_shp)
  
)