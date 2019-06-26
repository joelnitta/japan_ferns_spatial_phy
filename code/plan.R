# Define analysis plan
plan <- drake_plan (
  
  # Load and process raw data ----
  
  # Pteridophyte Phylogeny Group I (PPGI) taxonomy
  # - original version
  ppgi_raw = read_csv(file_in("data/ppgi_taxonomy.csv")),
  # - modify slightly for Pteridophytes of Japan
  ppgi = modify_ppgi(ppgi_raw),
  
  # Reproductive mode data, with one row per species.
  repro_data_raw = read_excel(file_in("data/ESM1.xlsx")),
  repro_data = process_repro_data(repro_data_raw),
  
  # Occurrence data, with multiple rows per species.
  # Occurrences are presences in a set of 1km2 grid 
  # cells across Japan, not actual occurrence points of specimens.
  occ_data_raw = read_excel(
    file_in("data/ESM2.xlsx"),
    col_types = c("numeric", "text", "numeric", 
                  "numeric", "numeric", "text", 
                  "text", "text")
  ),
  
  occ_data = clean_names(occ_data_raw) %>%
    add_taxonomy(ppgi),
  
  # Phylogenetic tree of all non-hybrid pteridophyte
  # taxa built from rbcL gene.
  japan_pterido_tree_raw = read.nexus("data/PD170708Bayes2.nxs"),
  japan_pterido_tree = format_tip_labels(japan_pterido_tree_raw),
  
  # Basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # List of all 1km2 grid cells across Japan.
  all_cells = read_excel(file_in("data/2_grid_cells_all.xlsx")) %>%
    rename(longitude = x, latitude = y, secondary_grid_code = id),
  
  # Analyze basic statistics ----
  
  # Count species per grid cell excluding hybrids.
  species_per_cell = count_species_per_cell(occ_data, repro_data),
  
  # Count grid cells per species (CPS).
  cells_per_species = count_cells_per_species(occ_data),
  
  # Count CPS by reproductive mode.
  cps_by_repro = count_cells_per_species_by_repro(
    occ_data, repro_data
  ),
  
  # Calculate mean CPS by reproductive mode.
  cps_by_repro_means = avg_cells_per_species_by_repro(
    cps_by_repro
  ),
  
  # Run analysis of variance (AOV) on CPS by reproductive mode.
  cps_by_repro_model_summary = aov(
    n_grids ~ reproductive_mode, 
    data = cps_by_repro) %>% tidy,
  
  # Make tibble with latitudinal breadth and
  # reproductive mode for all non-hybrids.
  lat_by_repro = count_lat_by_repro(occ_data, repro_data),
  
  # Calculate mean latitudinal breadth per species by reproductive mode.
  lat_by_repro_means = avg_lat_by_repro(lat_by_repro),
  
  # Run analysis of variance (AOV) on 
  # latitudinal breadth by reproductive mode.
  lat_by_repro_model = aov(
    lat_breadth ~ reproductive_mode, 
    data = lat_by_repro),
  
  lat_by_repro_model_summary = tidy(lat_by_repro_model),
  
  # AOV of latidudinal breadth by reproductive mode 
  # showed a significant difference, so
  # run Tukey HSD test on results.
  lat_by_repro_tukey = TukeyHSD(lat_by_repro_model) %>% tidy,
  
  # Analyze community diversity ----
  
  # Make richness matrix (number of species per
  # 1km2 grid cell).
  richness = make_richness_matrix(occ_data),
  
  # Make community matrix (presence/absence of each species in
  # 1km2 grid cells), trim to only species in tree.
  comm_matrix = make_comm_matrix(occ_data) %>% 
    match_comm_and_tree(japan_pterido_tree, "comm"),
  
  ### Calculate phylogenetic diversity
  #
  # Calculating standard effect size (SES) involves generating
  # hundreds of null communities per each observed community, and
  # takes a long time (4-5 hours for 999 reps per community).
  #
  # So we will slice up the dataset into chunks and run each in parallel
  # to speed things up.
  #
  # First convert the community matrix to a dataframe
  # with species as rownames so it can be sliced up by columns and
  # keep the same rownames in each slice.
  comm_df = column_to_rownames(comm_matrix, "species"),
  
  # Split the dataset into 4 chunks (number of CPUs on my laptop)
  comm_split = target(
    drake_slice(comm_df, slices = 4, index = i, margin = 2),
    transform = map(i = !!seq_len(4))
  ),
  
  # Caclulate standard effect size of Faith's phylogenetic
  # diversity on each slice of the dataset.
  pd = target(
    ses_pd(
      # convert community matrix back to tibble
      comm = comm_split %>% rownames_to_column("species") %>% as_tibble, 
      phy = japan_pterido_tree, n_reps = 499),
    transform = map(i = !!seq_len(4), comm_split)
  ),
  
  # Combine sliced results into single dataframe.
  all_pd = target(
    bind_rows(pd),
    transform = combine(pd)
  ),
  
  # Combine PD and richness into single dataframe.
  alpha_div = merge_metrics(all_pd, richness),
  
  # Plots ----
  
  richness_map = make_diversity_map(
    div_data = all_alpha_div, 
    world_map = world_map, 
    occ_data = occ_data, 
    div_metric = "richness", 
    metric_title = "Richness"
  ),
  
  ses_pd_map = make_diversity_map(
    div_data = all_alpha_div, 
    world_map = world_map, 
    occ_data = occ_data, 
    div_metric = "ses_pd", 
    metric_title = "SES of PD"
  ),
  
  ses_pd_highlight_map = make_pd_highlight_map(
    div_data = alpha_div, 
    world_map = world_map, 
    occ_data = occ_data
  )
)
