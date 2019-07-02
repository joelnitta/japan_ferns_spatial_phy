# Define analysis plan
plan <- drake_plan (
  
  # Load and process raw data ----
  
  # Pteridophyte Phylogeny Group I (PPGI) taxonomy
  # - original version
  # IMPORTANT: check encoding on local machine
  ppgi_raw = data.table::fread(file_in("data/ppgi_taxonomy.csv"), encoding = "Latin-1"),
  
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
  
  # - occurrence data including ferns north of 30.1° lat
  # (Yakushima and northward) only. There are 16 genera that only
  # occur south of 30.1 lat, which I think may be skewing the
  # SES PD metric, so try it without them.
  occ_data_ferns_north = occ_data_ferns %>%
    filter(latitude > 30.1),
  
  occ_data_ferns_south = occ_data_ferns %>%
    filter(latitude < 30.1),
  
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
  
  # - tree including ferns north of 30.1 lat only.
  japan_fern_tree_north = drop.tip(
    japan_pterido_tree, 
    setdiff(japan_pterido_tree$tip.lab, occ_data_ferns_north$taxon_id)
  ),
  
  japan_fern_tree_south = drop.tip(
    japan_pterido_tree, 
    setdiff(japan_pterido_tree$tip.lab, occ_data_ferns_south$taxon_id)
  ),
  
  # Basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # List of all 1km2 grid cells across Japan with elevation.
  all_cells = read_csv(
    file_in("data/all_cells_el.csv"),
    col_types = "nncc?"),
  
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
  
  # Analyze percentage of sexual diploid ferns.
  percent_sex_dip_ferns =
    right_join(
      occ_data_ferns, 
      select(repro_data, -taxon_name),
      by = "taxon_id") %>%
    group_by(secondary_grid_code, latitude, longitude) %>%
    summarize(
      num_sex_dip = sum(sexual_diploid),
      num_total = n()
    ) %>%
    ungroup %>%
    mutate(percent_sex_dip = num_sex_dip / num_total),
  
  # Analyze percentage of sexual diploid pteridophytes.
  percent_sex_dip_pteridos =
    right_join(
      occ_data_pteridos, 
      select(repro_data, -taxon_name),
      by = "taxon_id") %>%
    group_by(secondary_grid_code, latitude, longitude) %>%
    summarize(
      num_sex_dip = sum(sexual_diploid),
      num_total = n()
    ) %>%
    ungroup %>%
    mutate(percent_sex_dip = num_sex_dip / num_total),
  
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
  
  comm_ferns_north = make_comm_matrix(occ_data_ferns_north) %>% 
    match_comm_and_tree(japan_fern_tree_north, "comm"),
  
  comm_ferns_south = make_comm_matrix(occ_data_ferns_south) %>% 
    match_comm_and_tree(japan_fern_tree_south, "comm"),
  
  ### Calculate phylogenetic diversity for ferns in N and S areas separately
  
  # Convert to dataframe with rows as communities, columns as species
  # for picante
  comm_pteridos_df = column_to_rownames(comm_pteridos, "species") %>% t(),
  comm_ferns_df = column_to_rownames(comm_ferns, "species") %>% t(),
  comm_ferns_north_df = column_to_rownames(comm_ferns_north, "species") %>% t(),
  comm_ferns_south_df = column_to_rownames(comm_ferns_south, "species") %>% t(),
  
  mpd_pteridos = picante::ses.mpd(
    samp = comm_pteridos_df, 
    dis = cophenetic(japan_pterido_tree),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  mpd_ferns = picante::ses.mpd(
    samp = comm_ferns_df, 
    dis = cophenetic(japan_fern_tree),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  mpd_ferns_north = picante::ses.mpd(
      samp = comm_ferns_north_df, 
      dis = cophenetic(japan_fern_tree_north),
      null.model = "independentswap",
      iterations = 10000,
      runs = 999),
  
  mpd_ferns_south = picante::ses.mpd(
    samp = comm_ferns_south_df, 
    dis = cophenetic(japan_fern_tree_south),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  mntd_pteridos = picante::ses.mntd(
    samp = comm_pteridos_df, 
    dis = cophenetic(japan_pterido_tree),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  mntd_ferns = picante::ses.mntd(
    samp = comm_ferns_df, 
    dis = cophenetic(japan_fern_tree),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  mntd_ferns_north = picante::ses.mntd(
    samp = comm_ferns_north_df, 
    dis = cophenetic(japan_fern_tree_north),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  mntd_ferns_south = picante::ses.mntd(
    samp = comm_ferns_south_df, 
    dis = cophenetic(japan_fern_tree_south),
    null.model = "independentswap",
    iterations = 10000,
    runs = 999),
  
  # Combine PD and richness into single dataframe.
  # Add elevation and lat/longs for all 1km2 grid cells, even for those
  # that didn't have any species.
  mpd_ferns_ns = rbind(mpd_ferns_north, mpd_ferns_south) %>% 
    rownames_to_column("secondary_grid_code") %>%
    as_tibble %>%
    clean_names %>%
    select(-ntaxa, -runs),
  
  mntd_ferns_ns = rbind(mntd_ferns_north, mntd_ferns_south) %>% 
    rownames_to_column("secondary_grid_code") %>%
    as_tibble %>%
    clean_names %>%
    select(-ntaxa, -runs),
  
  # Combine all diversity metrics into single df by site
  alpha_div_ferns_ns = select(all_cells, secondary_grid_code, latitude, longitude, elevation) %>%
    left_join(select(richness_ferns, secondary_grid_code, richness)) %>%
    left_join(mpd_ferns_ns) %>%
    left_join(mntd_ferns_ns) %>%
    left_join(select(percent_sex_dip_ferns, secondary_grid_code, percent_sex_dip)) %>%
    mutate(richness = replace_na(richness, 0)),
  
  alpha_div_ferns = select(all_cells, secondary_grid_code, latitude, longitude, elevation) %>%
    left_join(select(richness_ferns, secondary_grid_code, richness)) %>%
    left_join(mpd_ferns %>% rownames_to_column("secondary_grid_code") %>%
                as_tibble %>%
                clean_names %>%
                select(-ntaxa, -runs)) %>%
    left_join(mntd_ferns %>% rownames_to_column("secondary_grid_code") %>%
                as_tibble %>%
                clean_names %>%
                select(-ntaxa, -runs)) %>%
    left_join(select(percent_sex_dip_ferns, secondary_grid_code, percent_sex_dip)) %>%
    mutate(richness = replace_na(richness, 0)),
  
  alpha_div_pteridos = select(all_cells, secondary_grid_code, latitude, longitude, elevation) %>%
    left_join(select(richness_pteridos, secondary_grid_code, richness)) %>%
    left_join(mpd_pteridos %>% rownames_to_column("secondary_grid_code") %>%
                as_tibble %>%
                clean_names %>%
                select(-ntaxa, -runs)) %>%
    left_join(mntd_pteridos %>% rownames_to_column("secondary_grid_code") %>%
                as_tibble %>%
                clean_names %>%
                select(-ntaxa, -runs)) %>%
    left_join(select(percent_sex_dip_pteridos, secondary_grid_code, percent_sex_dip)) %>%
    mutate(richness = replace_na(richness, 0)),
  
  # Write out manuscript ----
  ms = rmarkdown::render(
    knitr_in(here::here("ms/japan_pteridos_biodiv.Rmd")),
    output_file = file_out(here::here("ms/japan_pteridos_biodiv.html")),
    quiet = TRUE)
)
