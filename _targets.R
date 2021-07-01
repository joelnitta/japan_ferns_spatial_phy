library(targets)
library(tarchetypes)

# Load packages, functions, tests
source("R/packages.R")
source("R/functions.R")
source("R/tests.R")

# Set parallel back-end
plan(callr)

# Specify how to resolve futures for parallel tasks
# e.g., cpr_rand_test()
plan(multicore, workers = 60) # IMPORTANT: Change the plan / num. workers as needed for your system!

tar_plan(
  
  # Run tests on custom functions ----
  # (result should be NULL, or will stop everything and issue an error)
  test_results = run_tests(),
  
  # Load and process various data ----
  
  # Unzip data files from Ebihara and Nitta 2019
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  tar_file(ebihara_2019_zip_file, "data_raw/doi_10.5061_dryad.4362p32__v4.zip"),
  
  ebihara_2019_data = unzip_ebihara_2019(
    dryad_zip_file = ebihara_2019_zip_file, 
    exdir = "data_raw/ebihara_2019"
  ),
  
  # Split out paths for each unzipped data file
  tar_file(ppgi_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ppgi")]),
  tar_file(green_list_file, ebihara_2019_data[str_detect(ebihara_2019_data, "FernGreenListV1")]),
  tar_file(esm1_file, ebihara_2019_data[str_detect(ebihara_2019_data, "ESM1")]),
  
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
  
  # Collapse all the political units down to just one shape for the country
  japan_shp = jpn_pol %>%
    select(geometry) %>%
    summarize(),
  
  # Load manually entered points of interest for drawing a map of Japan
  tar_file(japan_points_raw_file, "data_raw/japan_points_raw.csv"),
  japan_points_raw = read_csv(japan_points_raw_file),
  
  # Prepare occurrence data ----
  # Summary: from raw occurrence data, test binning into grid cells at four scales
  # select the optimal scale, filter out poorly sampled grid cells, generate
  # community data matrix (sites x species)
  
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
  ),
  
  # After inspecting results, 0.2 is optimal scale.
  
  # Extract geographic shapes, richness, and number of specimens at 0.2 degree grid scale
  shape_ferns_full = shape_from_comm_scaled_list(comm_scaled_list, 0.2) %>%
    # calculate redundancy
    mutate(redundancy = 1 - (richness/abundance)),
  
  # Extract community matrix
  comm_ferns_full = comm_from_comm_scaled_list(comm_scaled_list, 0.2),
  
  # Subset geographic shapes to redundancy > 0.1
  shape_ferns = filter(shape_ferns_full, redundancy > 0.1),
  
  # Subset community matrix to communities with redundancy > 0.1
  comm_ferns = filter_comm_by_redun(
    comm = comm_ferns_full,
    shape = shape_ferns,
    cutoff = 0.1
  ),
  
  # Make community matrix subset to taxa endemic to Japan
  comm_ferns_endemic = subset_comm_to_endemic(
    comm = comm_ferns,
    green_list = green_list
  ),
  
  # # Phylogenetic analysis ----
  # # Summary: combine Japan rbcL sequences with global sampling, infer global tree,
  # # estimate divergence times, trim tree to just Japan ferns
  # 
  # # Read in Japan rbcL alignment
  # japan_rbcL_raw = read_ja_rbcL_from_zip(ebihara_2019_zip_file),
  # 
  # # Rename Japan rbcL to use species names instead of taxon ID
  # japan_rbcL = rename_alignment(
  #   alignment = japan_rbcL_raw,
  #   taxon_id_map = green_list),
  # 
  # # Read in list of globally sampled genes including rbcL
  # tar_file(ftol_zip_file, "data_raw/ftol_data_release_v0.0.1.zip"),
  # 
  # ftol_data = unzip_ftol(
  #   zip_file = ftol_zip_file,
  #   exdir = "data_raw"
  # ),
  # 
  # # Split out paths for each data file
  # tar_file(ftol_plastid_concat, ftol_data[str_detect(ftol_data, "ftol_plastid_concat")]),
  # tar_file(ftol_plastid_parts, ftol_data[str_detect(ftol_data, "ftol_plastid_parts")]),
  # 
  # # Load list of aligned genes
  # broad_alignment_list = load_ftol_alignment(ftol_plastid_concat, ftol_plastid_parts),
  # 
  # # Read in calibration dates
  # tar_file(calibration_dates_file, "data_raw/testo_sundue_2016_calibrations.csv"),
  # plastome_calibration_dates = load_calibration_dates(calibration_dates_file),
  # 
  # # Combine the Japan rbcL data with global sampling
  # # - replaces any species with the same name with those from Japan
  # # - re-aligns rbcL
  # # - concatenates all genes into single alignment
  # plastome_alignment = combine_ja_rbcL_with_global(
  #   broad_alignment_list = broad_alignment_list,
  #   japan_rbcL = japan_rbcL),
  # 
  # # Infer tree
  # plastome_tree = jntools::iqtree(
  #   plastome_alignment,
  #   m = "GTR+I+G", bb = 1000, nt = "AUTO",
  #   redo = TRUE, echo = TRUE, wd = here::here("iqtree")),
  # 
  # # Root tree on bryophytes
  # plastome_tree_rooted = ape::root(
  #   plastome_tree,
  #   c("Anthoceros_angustus", "Marchantia_polymorpha", "Physcomitrium_patens")),
  # 
  # # Date tree
  # 
  # # Run initial treepl search to identify smoothing parameter
  # treepl_cv_results = run_treepl_cv(
  #   phy = plastome_tree_rooted,
  #   alignment = plastome_alignment,
  #   calibration_dates = plastome_calibration_dates,
  #   cvstart = "1000",
  #   cvstop = "0.000001",
  #   plsimaniter = "200000", # preliminary output suggested > 100000
  #   seed = 7167,
  #   thorough = TRUE,
  #   wd = here::here("treepl"),
  #   nthreads = 1,
  #   echo = TRUE
  # ),
  # 
  # # Run priming analysis to determine optimal states for other parameters
  # treepl_priming_results = run_treepl_prime(
  #   phy = plastome_tree_rooted,
  #   alignment = plastome_alignment,
  #   calibration_dates = plastome_calibration_dates,
  #   cv_results = treepl_cv_results,
  #   plsimaniter = "200000", # preliminary output suggested > 100000
  #   seed = 7167,
  #   thorough = TRUE,
  #   wd = here::here("treepl"),
  #   nthreads = 1,
  #   echo = TRUE
  # ),
  # 
  # # Run treePL dating analysis
  # treepl_dating_results = run_treepl(
  #   phy = plastome_tree_rooted,
  #   alignment = plastome_alignment,
  #   calibration_dates = plastome_calibration_dates,
  #   cv_results = treepl_cv_results,
  #   priming_results = treepl_priming_results,
  #   plsimaniter = "200000", # preliminary output suggested > 100000
  #   seed = 7167,
  #   thorough = TRUE,
  #   wd = here::here("treepl"),
  #   nthreads = 7,
  #   echo = TRUE
  # ),
  # 
  # # Subset to just pteridophytes in Japan
  # japan_pterido_tree = ape::keep.tip(treepl_dating_results, rownames(japan_rbcL)),
  
  japan_pterido_tree = ape::read.tree("data_raw/japan_pterido_tree_dated.tre"),
  
  # Subset to only ferns
  japan_fern_tree = subset_tree(
    phy = japan_pterido_tree, 
    ppgi = ppgi),
  
  # Format trait data ----
  
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
  
  # Make trait dendrogram
  japan_fern_trait_tree = make_trait_tree(trait_distance_matrix),
  
  # Climate data ----
  
  # Load climate data downloaded from WorldClim database
  tar_file(worldclim_dir, "data_raw/world_clim"),
  ja_climate_data = load_ja_worldclim_data(worldclim_dir),
  
  # Calculate mean climate values in each grid cell
  mean_climate = calc_mean_climate(shape_ferns, ja_climate_data),
  
  # Randomization tests of diversity metrics ----
  # Run these in the main workflow as they are parallelized with {future}
  
  # - all ferns, phylogenetic diversity and endemism
  tar_target(
    rand_test_phy_ferns,
    cpr_rand_test(
      comm = comm_ferns,
      phy = japan_fern_tree,
      null_model = "independentswap",
      n_reps = 999,
      n_iterations = 100000,
      metrics = c("pd", "rpd", "pe", "rpe")),
    deployment = "main"
  ),
  
  # - endemic ferns, phylogenetic endemism
  tar_target(
    rand_test_phy_ferns_endemic,
    cpr_rand_test(
      comm = comm_ferns_endemic,
      phy = japan_fern_tree,
      null_model = "independentswap",
      n_reps = 999,
      n_iterations = 100000,
      metrics = c("pe", "rpe")),
    deployment = "main"
  ),
  
  # - all ferns, trait-based diversity
  tar_target(
    rand_test_traits_ferns,
    cpr_rand_test(
      comm = comm_ferns,
      phy = japan_fern_trait_tree,
      null_model = "independentswap",
      n_reps = 999,
      n_iterations = 100000,
      metrics = c("pd", "rpd")) %>%
      # rename columns as "fd" instead of "pd"
      rename_with(~str_replace_all(., "pd_", "fd_")),
    deployment = "main"
  ),
  
  # Bioregions analysis ----
  
  # Assess optimal K-value for clustering by taxonomy
  k_taxonomy = find_k_taxonomy(comm_ferns),
  
  # Assess optimal K-value for clustering by phylogeny
  k_phylogeny = find_k_phylogeny(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
  ),
  
  # Cluster by taxonomy
  regions_taxonomy = cluster_taxonomic_regions(
    comm_df = comm_ferns,
    k = k_taxonomy[["optimal"]][["k"]]
  ),
  
  # Cluster by phylogeny
  regions_phylogeny = cluster_phylo_regions(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
    k = k_phylogeny[["optimal"]][["k"]]
  ),
  
  # Traits analysis ----
  
  # Analyze phylogenetic signal
  
  # Specify continuous traits
  cont_traits = c("frond_width", "stipe_length", "number_pinna_pairs"),
  
  # Analyze phy signal in continuous traits with K and lambda
  phy_sig_results = map_df(
    cont_traits,
    ~analyze_cont_phylosig(selected_trait = ., traits = fern_traits, phy = japan_fern_tree)
  ),
  
  # Analyze phy signal in binary traits with D
  fern_traits_binary = select(fern_traits, -any_of(cont_traits)),
  
  binary_sig_results = analyze_binary_phylosig(fern_traits_binary, japan_fern_tree),
  
  # Summarize traits
  traits_summary = make_trait_summary(fern_traits),
  
  # Reproductive mode analysis ----
  # Calculate % apomictic species in each fern community
  percent_apo = calc_perc_apo(comm_ferns, repro_data),
  
  # Combine results ----
  
  # Combine spatial data, alpha diversity, environmental data, and regions, add significance and endemism types
  
  # - All ferns
  biodiv_ferns_spatial =
    shape_ferns %>%
    left_join(format_cpr_res(rand_test_phy_ferns), by = "grids") %>%
    left_join(format_cpr_res(rand_test_traits_ferns), by = "grids") %>%
    left_join(regions_taxonomy %>% rename(taxonomic_cluster = cluster), by = "grids") %>%
    left_join(regions_phylogeny %>% rename(phylo_cluster = cluster), by = "grids") %>%
    # Add environmental data
    left_join(mean_climate, by = "grids") %>%
    # Add % apomictic taxa
    left_join(percent_apo) %>%
    # Classify endemism and significance of randomization tests
    cpr_classify_endem() %>%
    classify_signif("pd") %>%
    classify_signif("rpd") %>%
    classify_signif("fd") %>%
    classify_signif("rfd") %>%
    classify_signif("pe", one_sided = TRUE, upper = TRUE) %>%
    # Format factors
    mutate(taxonomic_cluster = as.factor(taxonomic_cluster) %>% fct_infreq %>% as.numeric %>% as.factor) %>%
    mutate(phylo_cluster = as.factor(phylo_cluster) %>% fct_infreq %>% as.numeric %>% as.factor) %>%
    # Add redundancy
    mutate(redundancy = 1 - (richness/abundance)),
  
  # - Japan endemics only
  biodiv_ferns_endemic_spatial =
    shape_ferns %>%
    left_join(format_cpr_res(rand_test_phy_ferns_endemic), by = "grids") %>%
    # Classify endemism and significance of randomization tests
    cpr_classify_endem() %>%
    classify_signif("pe", one_sided = TRUE, upper = TRUE),
  
  # Spatial modeling ----
  
  # Model the effects of environment and percent apomictic taxa on each biodiversity metric, 
  # while accounting for spatial autocorrelation
  
  ## Setup ----
  
  # Define variables for models
  # - response variables for environmental model
  resp_vars_env = c("richness", "pd_obs_z", "fd_obs_z", "rpd_obs_z", "rfd_obs_z"),
  # - response variables for reproductive model
  resp_vars_repro = c("pd_obs_z", "rpd_obs_z"),
  # - independent variables (both models)
  indep_vars = c("percent_apo", "temp", "precip", "precip_season", "area"),
  
  # Calculate area as rolling mean in 1 degree latitudinal windows
  lat_area_ja = calc_area_by_lat(japan_shp, lat_cut = 0.2, lat_window = 1),
  
  # Make biodiversity metrics dataframe with centroid of each site for models
  biodiv_ferns_cent = sf_to_centroids(biodiv_ferns_spatial) %>%
    # add area
    add_roll_area(lat_area_ja) %>%
    # drop single percent_apo outlier
    drop_apo_outlier %>%
    # keep only variables needed for model and only rows with zero missing data
    filter_data_for_model(c("grids", "lat", "long", resp_vars_env, indep_vars)),
  
  ## Correlation analysis ----
  
  # Check for correlation between independent variables in repro data
  t_test_results = run_mod_ttest_ja(
    sf_to_centroids(biodiv_ferns_spatial) %>%
      add_roll_area(lat_area_ja), 
    vars_select = c("temp", "temp_season", "precip", "precip_season", "percent_apo", "area")
  ),
  
  ## Analyze Moran's I ----
  
  # Make list of distances for run_moran_mc()
  dist_list = make_dist_list(biodiv_ferns_cent),
  
  # Prepare datasets for looping
  data_for_moran = tibble(
    var = c(resp_vars_env, indep_vars),
    data = list(biodiv_ferns_cent),
    dist_list = list(dist_list)
  ),
  
  # Loop over datasets and calculate Moran's I
  tar_target(
    morans_i,
    run_moran_mc(
      var_name = data_for_moran$var[[1]], 
      biodiv_data = data_for_moran$data[[1]], 
      listw = data_for_moran$dist_list[[1]], 
      nsim = 1000
    ),
    pattern = map(data_for_moran)
  ),
  
  ## Build spatial models ----
  
  # Prepare datasets for looping
  data_for_spamm = prepare_data_for_spamm(
    resp_var_env = resp_vars_env,
    resp_var_repro = resp_vars_repro,
    biodiv_ferns_cent = biodiv_ferns_cent
  ),
  
  # Loop across each formula and build a spatial model
  tar_target(
    spatial_models,
    run_spamm(
      formula = data_for_spamm$formula[[1]], 
      data = data_for_spamm$data[[1]], 
      resp_var = data_for_spamm$resp_var[[1]]
    ),
    pattern = map(data_for_spamm)
  ),
  
  ## Conduct LRTs ----
  
  # Prepare data for looping
  data_for_lrt = prepare_data_for_lrt(spatial_models, biodiv_ferns_cent),
  
  # Conduct LRTs in loop
  tar_target(
    lrt_results,
    run_spamm_lrt(
      null_formula = data_for_lrt$null_formula[[1]],
      full_formula = data_for_lrt$full_formula[[1]],
      data = data_for_lrt$data[[1]], 
      resp_var = data_for_lrt$resp_var[[1]],  
      comparison = data_for_lrt$comparison[[1]]
    ),
    pattern = map(data_for_lrt)
  ),
  
  ## Summarize spatial models ----
  
  # - model statistics
  model_stats = get_model_stats(spatial_models),
  # - model parameters (fixed effects)
  model_params = get_model_params(spatial_models),
  # - comparison between temp and reproductive models by cAIC
  aic_env_repro = compare_aic_env_repro(spatial_models),
  
  ## Predict model fits ---

  # Extract model fits in loop: use temperature for environmental models, and percent_apo for repro models
  tar_target(
    model_fits,
    predict_fit(spatial_models),
    pattern = map(spatial_models)
  ),
  
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
  signif_cells_protected_area = calculate_protected_area(biodiv_ferns_spatial, protected_areas, japan_shp),
  
  # Render manuscript
  # - doc
  tar_render(
    ms_doc,
    path = "ms/manuscript.Rmd",
    output_format = "bookdown::word_document2",
    output_file = here::here("results/manuscript.docx"),
    params = list(doc_type = "doc")
  ),
  # - pdf
  tar_render(
    ms_pdf,
    path = "ms/manuscript.Rmd",
    output_format = "bookdown::pdf_document2",
    output_file = here::here("results/manuscript.pdf"),
    params = list(doc_type = "pdf")
  ),
  # - SI pdf
  tar_render(
    si_pdf,
    path = "ms/SI.Rmd",
    output_format = "bookdown::pdf_document2",
    output_file = here::here("results/supp_info.pdf"),
    params = list(doc_type = "pdf")
  ),
  # - SI appendix on data exploration for models
  tar_render(
    si_data_exploration,
    path = "ms/data_exploration.Rmd",
    output_format = "rmarkdown::html_document",
    output_file = here::here("results/data_exploration.html"),
    params = list(knit_type = "targets")
  )
)
