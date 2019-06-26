# Define analysis plan
plan <- drake_plan (
  
  # Load and process raw data ----
  
  # Pteridophyte Phylogeny Group I (PPGI) taxonomy
  # - original version
  ppgi_raw = read_csv(
    file_in("data/ppgi_taxonomy.csv")),
  
  # - modify slightly for Pteridophytes of Japan
  ppgi = modify_ppgi(ppgi_raw),
  
  # Reproductive mode data, with one row per species.
  repro_data_raw = read_excel(
    file_in("data/ESM1.xlsx"),
    col_types = c("text", "text", "text", "text", "text", 
                  "numeric", "numeric", "numeric")),
  
  repro_data = process_repro_data(repro_data_raw),
  
  # Occurrence data, with multiple rows per species.
  # Occurrences are presences in a set of 1km2 grid 
  # cells across Japan, not actual occurrence points of specimens.
  occ_data_raw = read_excel(
    file_in("data/ESM2.xlsx"),
    col_types = c("text", "text", "text", 
                  "numeric", "numeric", "text", 
                  "text", "text")
  ),
  
  # - occurrence data including ferns and lycophytes
  occ_data_pteridos = clean_names(occ_data_raw) %>%
    add_taxonomy(ppgi),
  
  # - occurrence data including ferns only
  occ_data_ferns = 
    occ_data_pteridos %>%
    assert(not_na, class) %>%
    filter(class == "Polypodiopsida"),
  
  # Read in raw phylogenetic tree of all non-hybrid pteridophyte
  # taxa based on rbcL gene.
  japan_pterido_tree_raw = read.nexus("data/PD170708Bayes2.nxs"),
  
  # Process trees.
  # - tree including ferns and lycophtyes
  japan_pterido_tree = format_tip_labels(japan_pterido_tree_raw),
  
  # - tree including ferns only
  japan_fern_tree = drop.tip(
    japan_pterido_tree, 
    setdiff(japan_pterido_tree$tip.lab, occ_data_ferns$taxon_id)
  ),
  
  # Basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # List of all 1km2 grid cells across Japan.
  all_cells = read_excel(file_in("data/2_grid_cells_all.xlsx")) %>%
    rename(longitude = x, latitude = y, secondary_grid_code = id),
  
  # Analyze basic statistics ----
  
  # Count species per grid cell excluding hybrids.
  species_per_cell = count_species_per_cell(occ_data_pteridos, repro_data),
  
  # Count grid cells per species (CPS).
  cells_per_species = count_cells_per_species(occ_data_pteridos),
  
  # Count CPS by reproductive mode.
  cps_by_repro = count_cells_per_species_by_repro(
    occ_data_pteridos, repro_data
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
  lat_by_repro = count_lat_by_repro(occ_data_pteridos, repro_data),
  
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
  richness_pteridos = make_richness_matrix(occ_data_pteridos),
  
  richness_ferns = make_richness_matrix(occ_data_ferns),
  
  # Make community matrix (presence/absence of each species in
  # 1km2 grid cells), trim to only species in tree.
  comm_pteridos = make_comm_matrix(occ_data_pteridos) %>% 
    match_comm_and_tree(japan_pterido_tree, "comm"),
  
  comm_ferns = make_comm_matrix(occ_data_ferns) %>% 
    match_comm_and_tree(japan_fern_tree, "comm"),
  
  ### Calculate phylogenetic diversity
  #
  # Calculating standard effect size (SES) involves generating
  # hundreds of null communities per each observed community, and
  # takes a long time (4-5 hours for 999 reps per community).
  #
  # So we will slice up the dataset into chunks and run each in parallel
  # to speed things up.
  #
  # Do for both all pteridophytes together, and ferns only.
  
  # Split the dataset into 4 chunks (number of CPUs on my laptop).
  # (convert the community matrix to a dataframe
  # with species as rownames so it can be sliced up by columns and
  # keep the same rownames in each slice)
  comm_pteridos_split = target(
    drake_slice(column_to_rownames(comm_pteridos, "species"), 
                slices = 4, index = i, margin = 2),
    transform = map(i = !!seq_len(4))
  ),
  
  comm_ferns_split = target(
    drake_slice(column_to_rownames(comm_ferns, "species"), 
                slices = 4, index = i, margin = 2),
    transform = map(i = !!seq_len(4))
  ),
  
  # Caclulate standard effect size of Faith's phylogenetic
  # diversity on each slice of the dataset.
  pd_pteridos = target(
    ses_pd(
      # convert community matrix back to tibble
      comm = comm_pteridos_split %>% rownames_to_column("species") %>% as_tibble, 
      phy = japan_pterido_tree, n_reps = 999),
    transform = map(i = !!seq_len(4), comm_pteridos_split)
  ),
  
  pd_ferns = target(
    ses_pd(
      # convert community matrix back to tibble
      comm = comm_ferns_split %>% rownames_to_column("species") %>% as_tibble, 
      phy = japan_fern_tree, n_reps = 999),
    transform = map(i = !!seq_len(4), comm_ferns_split)
  ),
  
  # Combine sliced results into single dataframe.
  all_pd_pteridos = target(
    bind_rows(pd_pteridos),
    transform = combine(pd_pteridos)
  ),
  
  all_pd_ferns = target(
    bind_rows(pd_ferns),
    transform = combine(pd_ferns)
  ),
  
  # Combine PD and richness into single dataframe.
  alpha_div_pteridos = merge_metrics(all_pd_pteridos, richness_pteridos),
  alpha_div_ferns = merge_metrics(all_pd_ferns, richness_ferns),
  
  # Plots ----
  
  richness_pteridos_map = make_diversity_map(
    div_data = alpha_div_pteridos, 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "richness", 
    metric_title = "Richness"
  ),
  
  ses_pd_fern_map = make_diversity_map(
    div_data = alpha_div_ferns, 
    world_map = world_map, 
    occ_data = occ_data_pteridos, 
    div_metric = "ses_pd", 
    metric_title = "SES of PD"
  ),
  
  ses_pd_highlight_fern_map = make_pd_highlight_map(
    div_data = alpha_div_ferns, 
    world_map = world_map, 
    occ_data = occ_data_pteridos
  )
)
