plan <- drake_plan(
  
  # Read in green list for standardizing names
  green_list = read_greenlist_from_zip("data_raw/doi_10.5061_dryad.4362p32__v4.zip") %>% 
    tidy_japan_names(),
  
  # Read in Japan rbcL alignment
  japan_rbcL_raw = read_ja_rbcL_from_zip("data_raw/doi_10.5061_dryad.4362p32__v4.zip") %>% as.matrix,
  
  # Rename Japan rbcL to use species names instead of taxon ID
  japan_rbcL = rename_alignment(
    alignment = japan_rbcL_raw, 
    taxon_id_map = green_list),
  
  # Read in list of 71 globally sampled genes including rbcL 
  broad_alignment_list = readRDS("data_raw/plastid_genes_aligned_trimmed_renamed.RDS"),
  
  # Read in calibration dates
  plastome_calibration_dates = load_calibration_dates(
    file_in("data_raw/testo_sundue_2016_calibrations.csv")),
  
  # Combine the Japan rbcL data with global sampling
  # - replaces any species with the same name with those from Japan
  # - re-aligns rbcL
  # - concatenates all genes into single alignment
  plastome_alignment = combine_ja_rbcL_with_global(
    broad_alignment_list = broad_alignment_list, 
    japan_rbcL = japan_rbcL),
  
  # Infer tree
  plastome_tree = jntools::fasttree(plastome_alignment),
  
  # Root tree on bryophytes
  plastome_tree_rooted = ape::root(
    plastome_tree,
    c("Anthoceros_angustus", "Marchantia_polymorpha", "Physcomitrella_patens")),
  
  # Date tree
  
  # Run initial treepl search to identify smoothing parameter
  treepl_cv_results = run_treepl_cv(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cvstart = "1000",
    cvstop = "0.000001",
    plsimaniter = "200000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 1,
    echo = TRUE
  ),
  
  # Run priming analysis to determine optimal states for other parameters
  treepl_priming_results = run_treepl_prime(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cv_results = treepl_cv_results,
    plsimaniter = "200000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 1,
    echo = TRUE
  ),
  
  # Run treePL dating analysis
  treepl_dating_results = run_treepl(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cv_results = treepl_cv_results,
    priming_results = treepl_priming_results,
    plsimaniter = "200000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 7,
    echo = TRUE
  )
  
)