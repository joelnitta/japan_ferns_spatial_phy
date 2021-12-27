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
# IMPORTANT: Change the plan / num. workers as needed for your system!
plan(multicore, workers = 30)

tar_plan(

  # Run tests on custom functions ----
  # (result should be NULL, or will stop everything and issue an error)
  test_results = run_tests(),

  # Load and process various data ----

  ## Data files from Ebihara and Nitta 2019 ----
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  tar_file(ebihara_2019_zip_file, "data/doi_10.5061_dryad.4362p32__v4.zip"),

  # - Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = load_ppgi(ebihara_2019_zip_file),

  # - Fern Green List (official taxonomy + conservation status for each species)
  green_list = load_green_list(ebihara_2019_zip_file),

  # - Raw reproductive mode data (Ebihara et al 2019 ESM1)
  repro_data_raw = load_repro_data(ebihara_2019_zip_file),

  # - Japan rbcL alignment
  japan_rbcL_raw = read_ja_rbcL_from_zip(ebihara_2019_zip_file),

  # - Summary of occurrence point data
  tar_file(occ_point_data_summary_file, "data/japan_ferns_occ_summary.csv"),
  occ_point_data_summary = read_csv(occ_point_data_summary_file),

  ## Occurrence data ----
  # Raw occurrence data processing was done in R/process_raw_data.R

  # Use Mollweide equal-area projection with longitude centered on Japan
  mollweide = "+proj=moll +lon_0=135",

  # Load grid cells, richness, abundance, redundancy at 20 km grid scale
  tar_file(japan_ferns_shape_full_file, "data/japan_ferns_shape_full.gpkg"),
  shape_ferns_full = load_jferns_shape(
    japan_ferns_shape_full_file, crs = mollweide),

  # Load community matrix at 20 km grid scale (not filtered by redundancy)
  tar_file(japan_ferns_comm_full_file, "data/japan_ferns_comm_full.csv"),
  comm_ferns_full = load_jferns_comm(japan_ferns_comm_full_file),

  # Load summaries of testing different grid cell sizes
  tar_file(redundancy_by_res_file, "data/japan_ferns_redundancy_by_res.csv"),
  redundancy_by_res = read_csv(
    redundancy_by_res_file,
    col_types = cols(res = col_character())
  ),

  # Load results of assessing sampling completeness with iNEXT
  tar_file(inext_res_file, "data/japan_ferns_inext_results.csv"),
  inext_res = read_csv(inext_res_file),

  # Subset geographic shapes to redundancy > 0.1
  shape_ferns = filter(shape_ferns_full, redundancy > 0.1),

  # Subset community matrix to communities with redundancy > 0.1,
  # drop any remaining empty species or sites
  comm_ferns = filter_comm_by_redun(
    comm = comm_ferns_full,
    shape = shape_ferns,
    cutoff = 0.1
  ) %>% drop_empty(),

  # Make community matrix subset to taxa endemic to Japan,
  # drop any remaining empty species or sites
  comm_ferns_endemic = subset_comm_to_endemic(
    comm = comm_ferns,
    green_list = green_list
  ) %>% drop_empty(),

  # Summarize latitudinal span by taxon
  lat_span_summary = summarize_fern_lat_span(comm_ferns, shape_ferns),

  ## Map data ----
  # Load (unlabled) map of Japan
  tar_file(japan_map_file, "data/japan_map.gpkg"),
  japan_shp = st_read(japan_map_file) %>% st_transform(mollweide),

  # Load manually entered points of interest for drawing a map of Japan
  tar_file(japan_map_points_file, "data/japan_map_points.csv"),
  japan_map_points = read_csv(japan_map_points_file),

  # Calculate latitudinal area (sq km) as rolling mean in 100 km latitudinal
  # windows across 20 km bands
  lat_area_ja = calc_area_by_lat(
    japan_shp, lat_cut = 20000, lat_window = 100000),

  # Phylogenetic analysis ----
  # Summary: combine Japan rbcL sequences with global sampling,
  # infer global tree, estimate divergence times, trim tree to just Japan ferns

  # Rename Japan rbcL to use species names instead of taxon ID
  japan_rbcL = rename_alignment(
    alignment = japan_rbcL_raw,
    taxon_id_map = green_list),

  # Read in list of globally sampled genes including rbcL
  tar_file(ftol_zip_file, "data/ftol_data_release_v0.0.1.zip"),
  broad_alignment_list = load_ftol_alignment(ftol_zip_file),

  # Read in calibration dates
  tar_file(calibration_dates_file, "data/testo_sundue_2016_calibrations.csv"),
  plastome_calibration_dates = load_calibration_dates(calibration_dates_file),

  # Combine the Japan rbcL data with global sampling
  # - replaces any species with the same name with those from Japan
  # - re-aligns rbcL
  # - concatenates all genes into single alignment
  plastome_alignment = combine_ja_rbcL_with_global(
    broad_alignment_list = broad_alignment_list,
    japan_rbcL = japan_rbcL),

  # Infer tree
  plastome_tree = jntools::iqtree(
    plastome_alignment,
    m = "GTR+I+G", bb = 1000, nt = "AUTO",
    redo = TRUE, echo = TRUE, wd = here::here("iqtree")),

  # Root tree on bryophytes
  plastome_tree_rooted = ape::root(
    plastome_tree,
    c("Anthoceros_angustus", "Marchantia_polymorpha", "Physcomitrium_patens")),

  # Date tree

  # Run priming analysis to determine optimal states for other parameters
  treepl_priming_results = run_treepl_prime(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    seed = 7167,
    thorough = TRUE,
    write_tree = TRUE,
    wd = here::here("treepl")
  ),

  # Run initial treepl search to identify smoothing parameter
  treepl_cv_results = run_treepl_cv(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    priming_results = treepl_priming_results,
    cvstart = "1000",
    cvstop =  "0.000001",
    seed = 7167,
    thorough = TRUE,
    write_tree = FALSE,
    wd = here::here("treepl"),
    nthreads = 20
  ),

  # Run treePL dating analysis
  treepl_dating_results = run_treepl(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cv_results = treepl_cv_results,
    priming_results = treepl_priming_results,
    plsimaniter = "1000000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 20
  ),

  # Subset to just pteridophytes in Japan
  japan_pterido_tree = ape::keep.tip(
    treepl_dating_results, rownames(japan_rbcL)),

  # Subset to only ferns, without collapsing any branches
  japan_fern_tree_uncollapsed = subset_tree(
    phy = japan_pterido_tree,
    ppgi = ppgi),

  # Collapse identical tips to zero-length polytomies
  japan_fern_tree = collapse_identical_tips(
    japan_fern_tree_uncollapsed, japan_rbcL),

  # Also make phylogram (not ultrametric tree) of ferns in Japan
  japan_fern_phylogram = ape::keep.tip(
    plastome_tree_rooted, japan_fern_tree$tip.label),

  # Format trait data ----

  # Format trait data, subset to ferns in tree
  tar_file(raw_trait_data_file, "data/japan_ferns_traits_lucid.csv"),

  fern_traits = format_traits(traits_lucid_path = raw_trait_data_file),

  # Transform continuous traits before making distance matrix
  traits_for_dist = transform_traits(
    fern_traits,
    trans_select = c("frond_width", "stipe_length", "number_pinna_pairs"),
    scale_select = c("frond_width", "stipe_length", "number_pinna_pairs")
  ) %>%
    # Make sure all traits are scaled within -1 to 1
    assert(within_bounds(-1, 1), where(is.numeric)),

  # Make trait distance matrix using taxon IDs as labels
  trait_distance_matrix = make_trait_dist_matrix(traits_for_dist),

  # Make trait dendrogram
  japan_fern_trait_tree = make_trait_tree(trait_distance_matrix),

  # Climate data ----

  # Load climate data from WorldClim database at 2.5 minute resolution
  tar_file(ja_climate_data_file, "data/japan_climate.gpkg"),
  ja_climate_data = sf::st_read(ja_climate_data_file) %>%
    st_transform(mollweide),

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
      null_model = "curveball",
      n_reps = 999,
      n_iterations = 100000,
      metrics = c("pd", "rpd", "pe", "rpe")),
    deployment = "main"
  ),

  # - all ferns, phylogenetic diversity, using phylogram (not ultrametric)
  tar_target(
    rand_test_phy_ferns_not_ult,
    cpr_rand_test(
      comm = comm_ferns,
      phy = japan_fern_phylogram,
      null_model = "curveball",
      n_reps = 999,
      n_iterations = 100000,
      metrics = c("pd", "rpd")),
    deployment = "main"
  ),

  # - endemic ferns, phylogenetic endemism
  tar_target(
    rand_test_phy_ferns_endemic,
    cpr_rand_test(
      comm = comm_ferns_endemic,
      phy = japan_fern_tree,
      null_model = "curveball",
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
      null_model = "curveball",
      n_reps = 999,
      n_iterations = 100000,
      metrics = c("pd", "rpd")) %>%
      # rename columns as "fd" instead of "pd"
      rename_with(~str_replace_all(., "pd_", "fd_")),
    deployment = "main"
  ),

  # Bioregions analysis ----

  # Assess optimal K-value for clustering by taxonomy
  k_taxonomy = find_k_taxonomy(comm_ferns, k_max = 10),

  # Assess optimal K-value for clustering by phylogeny
  k_phylogeny = find_k_phylogeny(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
    k_max = 10
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

  # Combine results
  # - original (default) IDs for bioregions
  bioregions_raw = combine_bioregions(regions_taxonomy, regions_phylogeny),
  # - bioregions relabeled in order of mean latitude
  bioregion_cutoff = 2, # later, lump bioregions into "other" if they have fewer than this number of grid-cells # nolint
  bioregions = relabel_bioregions_by_lat(
    bioregions_raw, shape_ferns, cutoff = bioregion_cutoff),

  # Traits analysis ----

  # Analyze phylogenetic signal

  # Specify continuous traits
  cont_traits = c("frond_width", "stipe_length", "number_pinna_pairs"),

  # Analyze phy signal in continuous traits with K and lambda
  phy_sig_results = map_df(
    cont_traits,
    ~analyze_cont_phylosig(
      selected_trait = .,
      traits = fern_traits,
      phy = japan_fern_tree_uncollapsed
    )
  ),

  # Analyze phy signal in binary traits with D
  fern_traits_binary = select(fern_traits, -any_of(cont_traits)),

  binary_sig_results = analyze_binary_phylosig(
    fern_traits_binary, japan_fern_tree_uncollapsed),

  # Summarize traits
  traits_summary = make_trait_summary(fern_traits),

  # Reproductive mode analysis ----

  # Clean up reproductive mode data
  repro_data = process_repro_data(repro_data_raw, green_list),

  # Calculate % apomictic species in each fern community
  percent_apo = calc_perc_apo(comm_ferns, repro_data),

  # Combine results ----

  # Combine spatial data, alpha diversity, environmental data, and regions,
  # add significance and endemism types

  # - All ferns
  biodiv_ferns_spatial =
    # Combine phylodiversity and bioregions
    shape_ferns %>%
    left_join(format_cpr_res(rand_test_phy_ferns), by = "grids") %>%
    left_join(format_cpr_res(rand_test_traits_ferns), by = "grids") %>%
    left_join(bioregions, by = "grids") %>%
    # Add environmental data
    add_roll_area(lat_area_ja) %>%
    left_join(mean_climate, by = "grids") %>%
    # Add % apomictic taxa
    left_join(percent_apo, by = "grids") %>%
    # Classify endemism and significance of randomization tests
    cpr_classify_endem() %>%
    classify_signif("pd") %>%
    classify_signif("rpd") %>%
    classify_signif("fd") %>%
    classify_signif("rfd") %>%
    classify_signif("pe", one_sided = TRUE, upper = TRUE) %>%
    # Add richness percent rank
    mutate(richness_obs_p_upper = dplyr::percent_rank(richness)),

  # - Japan endemics only
  biodiv_ferns_endemic_spatial =
    shape_ferns %>%
    left_join(format_cpr_res(rand_test_phy_ferns_endemic), by = "grids") %>%
    # Classify endemism and significance of randomization tests
    cpr_classify_endem() %>%
    classify_signif("pe", one_sided = TRUE, upper = TRUE),

  # - Phylogram (not ultrametric)
  biodiv_ferns_spatial_not_ult =
    shape_ferns %>%
    left_join(format_cpr_res(rand_test_phy_ferns_not_ult), by = "grids") %>%
    classify_signif("pd") %>%
    classify_signif("rpd") %>%
    # Add environmental data, % apomictic taxa
    add_roll_area(lat_area_ja) %>%
    left_join(mean_climate, by = "grids") %>%
    left_join(percent_apo),

  # Spatial modeling ----

  # Model the effects of environment and percent apomictic taxa on each
  # biodiversity metric, while accounting for spatial autocorrelation

  ## Setup ----

  # Define variables for models
  # - response variables for environmental model
  resp_vars_env = c(
    "richness", "pd_obs_z", "fd_obs_z", "rpd_obs_z", "rfd_obs_z"),
  # - response variables for reproductive model
  resp_vars_repro = c("pd_obs_z", "rpd_obs_z"),
  # - independent variables (both models)
  indep_vars = c("percent_apo", "temp", "precip", "precip_season", "lat_area"),

  # Make biodiversity metrics dataframe with centroid of each site for models
  # Drops single outlier site with extremely high % apomictic taxa and low
  # richness
  # and one site missing environmental data.
  # - ultrametric tree (full  analysis)
  biodiv_ferns_cent = spatial_to_cent_for_model(
    biodiv_ferns_spatial,
    c("grids", "lat", "long", resp_vars_env, indep_vars)),
  # - non ultrametric tree (% apo only)
  biodiv_ferns_not_ult = spatial_to_cent_for_model(
    biodiv_ferns_spatial_not_ult,
    c("grids", "lat", "long", resp_vars_repro, indep_vars)),

  ## Correlation analysis ----

  # Check for correlation between independent variables in repro data
  t_test_results = run_mod_ttest_ja(
    biodiv_ferns_spatial,
    vars_select = c(
      "temp", "temp_season", "precip",
      "precip_season", "percent_apo", "lat_area")
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
  # - ultrametric tree (full  analysis)
  data_for_spamm = prepare_data_for_spamm(
    resp_var_env = resp_vars_env,
    resp_var_repro = resp_vars_repro,
    biodiv_ferns_cent = biodiv_ferns_cent
  ),
  # - non ultrametric tree (% apo only)
  data_for_spamm_not_ult = prepare_data_for_spamm(
    resp_var_env = resp_vars_env,
    resp_var_repro = resp_vars_repro,
    biodiv_ferns_cent = biodiv_ferns_not_ult
  ) %>% filter(str_detect(formula, "percent_apo")),

  # Loop across each formula and build a spatial model
  # - ultrametric tree (full  analysis)
  tar_target(
    spatial_models,
    run_spamm(
      formula = data_for_spamm$formula[[1]],
      data = data_for_spamm$data[[1]],
      resp_var = data_for_spamm$resp_var[[1]]
    ),
    pattern = map(data_for_spamm)
  ),
  # - non ultrametric tree (% apo only)
  tar_target(
    spatial_models_not_ult,
    run_spamm(
      formula = data_for_spamm_not_ult$formula[[1]],
      data = data_for_spamm_not_ult$data[[1]],
      resp_var = data_for_spamm_not_ult$resp_var[[1]]
    ),
    pattern = map(data_for_spamm_not_ult)
  ),

  ## Conduct LRTs ----

  # Prepare data for looping (full analysis only)
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

  # Extract model fits in loop: use temperature for environmental models,
  # and percent_apo for repro models
  # - ultrametric tree (full  analysis)
  tar_target(
    model_fits,
    predict_fit(spatial_models),
    pattern = map(spatial_models)
  ),
  # - non ultrametric tree (% apo only)
  tar_target(
    model_fits_not_ult,
    predict_fit(spatial_models_not_ult),
    pattern = map(spatial_models_not_ult)
  ),

  # Conservation analysis ----

  ## Read in protected areas, assign protection levels following
  # Kusumoto et al. 2017:
  # - high: no human activities allowed
  # - medium: permission required for economic activities
  tar_file(protected_areas_file, "data/japan_protected_areas.gpkg"),
  protected_areas = st_read(protected_areas_file) %>% st_transform(mollweide),

  ## Read in deer distribution map
  tar_file(deer_range_file, "data/japan_deer_range.gpkg"),
  deer_range = st_read(deer_range_file) %>% st_transform(mollweide),

  # Crop siginficantly diverse grid cells to protected areas,
  # and to areas with deer
  protected_biodiv = crop_by_pa(
    protected_areas, biodiv_ferns_spatial, japan_shp),
  deer_danger_biodiv = crop_by_deer(
    deer_range, biodiv_ferns_spatial, japan_shp),

  # Calculate percentage of diverse cells in protected areas,
  # and in danger due to herbivory by deer
  signif_cells_protected_area = calculate_percent_protected(protected_biodiv),
  signif_cells_deer_danger = calculate_percent_deer_danger(deer_danger_biodiv),

  # Write out selected results files for dryad ----

  # Choose variables to include in biodiv results, convert to centroids
  biodiv_ferns_cent_dryad = biodiv_ferns_spatial_to_cent(
    biodiv_ferns_spatial,
    crs = 4612 # convert to JDG2000 CRS
    ),

  # japan_ferns_biodiv.csv
  tar_file(
    japan_ferns_biodiv_dryad_file,
    write_csv_tar(
      biodiv_ferns_cent_dryad, "results/dryad_files/japan_ferns_biodiv.csv")
  ),

  # japan_ferns_comm.csv
  tar_file(
    comm_ferns_dryad_file,
    comm_ferns %>% rownames_to_column("grids") %>%
      write_csv_tar("results/dryad_files/japan_ferns_comm.csv")
  ),

  # japan_ferns_shape.gpkg
  tar_file(
    shape_ferns_dryad_file,
    st_write_tar(
      shape_ferns,
      "results/dryad_files/japan_ferns_shape.gpkg",
      time_stamp = as.Date("2021-09-03")
    )
  ),

  # japan_ferns_traits.csv
  tar_file(
    fern_traits_dryad_file,
    write_csv_tar(fern_traits, "results/dryad_files/japan_ferns_traits.csv")
  ),

  # japan_ferns_tree.tre
  tar_file(
    japan_fern_phylogram_dryad_file,
    write_tree_tar(japan_fern_phylogram,
      "results/dryad_files/japan_ferns_tree.tre")
  ),

  # japan_ferns_tree_dated.tre
  tar_file(
    japan_fern_tree_dryad_file,
    write_tree_tar(japan_fern_tree,
      "results/dryad_files/japan_ferns_tree_dated.tre")
  ),

  # japan_ferns_tree_dated_uncollapsed.tre
  tar_file(
    japan_fern_tree_uncollapsed_dryad_file,
    write_tree_tar(japan_fern_tree_uncollapsed,
      "results/dryad_files/japan_ferns_tree_dated_uncollapsed.tre")
  ),

  # Etc ----
  # Lump regions with a small number of grid-cells
  # (fewer than `bioregion_cutoff`)
  biodiv_ferns_spatial_lumped = biodiv_ferns_spatial %>%
    mutate(taxonomic_cluster = fct_lump_min(
      taxonomic_cluster, bioregion_cutoff)
    ) %>%
    mutate(phylo_cluster = fct_lump_min(phylo_cluster, bioregion_cutoff)),

  # Render manuscript ----
  # Track ms files
  tar_file(refs_yaml, "ms/references.yaml"),
  tar_file(refs_other_yaml, "ms/references_other.yaml"),
  tar_file(template_file, "ms/template.docx"),
  tar_file(csl_file, "ms/apa-6th-edition.csl"),
  tar_file(ms_functions, "R/ms_functions.R"),

  # MS, docx format
  tar_render(
    ms_doc,
    knit_root_dir = here::here(),
    path = "ms/manuscript.Rmd",
    output_format = "bookdown::word_document2",
    output_file = here::here("results/manuscript.docx"),
    params = list(doc_type = "doc")
  ),
  # MS, pdf format
  tar_render(
    ms_pdf,
    knit_root_dir = here::here(),
    path = "ms/manuscript.Rmd",
    output_format = "bookdown::pdf_document2",
    output_file = here::here("results/manuscript.pdf"),
    params = list(doc_type = "pdf")
  ),
  # SI appendix 1
  tar_render(
    si_pdf,
    knit_root_dir = here::here(),
    path = "ms/SI.Rmd",
    output_format = "bookdown::pdf_document2",
    output_file = here::here("results/appendix_s1.pdf"),
    params = list(doc_type = "pdf")
  ),
  # SI appendix 2 on data exploration for models
  tar_render(
    si_data_exploration,
    knit_root_dir = here::here(),
    path = "ms/data_exploration.Rmd",
    output_format = "bookdown::pdf_document2",
    output_file = here::here("results/appendix_s2.pdf"),
    params = list(knit_type = "targets")
  ),
  tar_render(
    data_readme,
    knit_root_dir = here::here(),
    path = "ms/data_readme.Rmd",
    output_file = here::here("results/data_readme.rtf")
  )
)
