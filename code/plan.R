# Set vector of k-values to use for ecostructure
k_vals <- 2:10

# Define analysis plan
plan <- drake_plan (
  
  # Load and process raw data ----
  
  ### Taxonomy data ###
  
  # Pteridophyte Phylogeny Group I (PPGI) taxonomy
  # - original version
  ppgi_raw = read_csv(file_in("data_raw/ppgi_taxonomy.csv")),
  
  # - modify slightly for Pteridophytes of Japan
  ppgi = modify_ppgi(ppgi_raw),
  
  # Catalog of Life (COL) plants taxonomic data
  col_plants_raw = data.table::fread(here::here(
    "data_raw/archive-kingdom-plantae-phylum-tracheophyta-bl3/taxa.txt"
  ), encoding = "Latin-1") %>%
    as_tibble(),
  
  # Extract World Ferns database contained within COL, only use simple set of columns
  world_ferns = filter(
    col_plants_raw,
    str_detect(datasetName, "World Ferns")
  ) %>%
    select(c(
      "taxonID", "acceptedNameUsageID",
      "taxonomicStatus", "taxonRank",
      "scientificName", "genus", 
      "specificEpithet", "infraspecificEpithet"
    )),
  
  # Load Fern Green List, with conservation status for each species.
  green_list = read_excel(file_in("data_raw/FernGreenListV1.01.xls")) %>%
    select(taxon_id = ID20160331, scientific_name = `GreenList学名`,
           endemic = `固有`, conservation_status = `RL2012`) %>%
    mutate(taxon_id = as.character(taxon_id)),
  
  # Match fern and pteridophyte names to COL.
  gnr_results = match_with_gnr(green_list$scientific_name),
  
  # Resolve synonyms, drop any names that couldn't be resolved.
  resolved_names = resolve_synonyms(gnr_results, world_ferns) %>% 
    filter(!is.na(scientificName)),
  
  # Reproductive mode data, with one row per species.
  repro_data_raw = read_csv(
    file_in("data_raw/ESM1.csv"),
    col_types = "cccccnnnnn"),
  
  repro_data = process_repro_data(repro_data_raw),
  
  # Occurrence data, with multiple rows per species.
  # Occurrences are presences in a set of 1km2 grid 
  # cells across Japan, not actual occurrence points of specimens.
  occ_data_raw = read_csv(
    file_in("data_raw/ESM2.csv"),
    col_types = "cccnnccc"
  ),
  
  # - occurrence data including ferns and lycophytes
  occ_data_pteridos = clean_names(occ_data_raw) %>%
    add_taxonomy(ppgi) %>%
    # Verify that all secondary_grid_code values are in the all_cells data
    assert(in_set(all_cells$secondary_grid_code), secondary_grid_code),
  
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
  
  # GBIF data: cleaned points for all pteridophytes globally,
  # with names standardized to COL to species level (no infrasp. taxa)
  gbif_points_global = read_csv(file_in("data_raw/gbif_clean_no_obs.csv")),
  
  # Read in phylogenetic tree of all non-hybrid pteridophyte
  # taxa based on rbcL gene.
  japan_pterido_tree = read_nexus_in_zip(
    file_in("data_raw/japan_pterido_rbcl_cipres.zip"), 
    "infile.nex.con.tre")[[2]] %>%
    format_tip_labels,
  
  # Basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # List of all 10 km grid cells across Japan with elevation
  all_cells = read_csv(file_in("data/all_cells_el.csv")) %>%
    assert(is_uniq, secondary_grid_code),
  
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
    match_comm_and_tree(japan_pterido_tree, "comm"),
  
  comm_ferns_north = make_comm_matrix(occ_data_ferns_north) %>% 
    match_comm_and_tree(japan_pterido_tree, "comm"),
  
  comm_ferns_south = make_comm_matrix(occ_data_ferns_south) %>% 
    match_comm_and_tree(japan_pterido_tree, "comm"),
  
  ### Calculate phylogenetic diversity ###
  phy_mpd = target(
    ses_phy_mpd(
      comm,
      japan_pterido_tree,
      null.model = "independentswap",
      iterations = 10000,
      runs = 999
    ),
    transform = map(comm = c(comm_pteridos, comm_ferns, 
                             comm_ferns_north, comm_ferns_south))
  ),
  
  phy_mntd = target(
    ses_phy_mntd(
      comm,
      japan_pterido_tree,
      null.model = "independentswap",
      iterations = 10000,
      runs = 999
    ),
    transform = map(comm = c(comm_pteridos, comm_ferns, 
                             comm_ferns_north, comm_ferns_south))
  ),
  
  ### Analyze functional trait diversity ###
  
  # Make trait distance matrix
  # - first make tibble mapping taxon IDs to species names
  taxon_id_map = make_taxon_id_map(occ_data_pteridos),
  
  # - make trait distance matrix using taxon IDs as labels
  trait_distance_matrix = make_traits_dist_matrix(
    file_in("data_raw/JpFernLUCID_forJoel.xlsx"),
    taxon_id_map),
  
  func_mpd = target(
    ses_func_mpd(
      comm,
      trait_distance_matrix,
      null.model = "independentswap",
      iterations = 10000,
      runs = 999
    ),
    transform = map(comm = c(comm_pteridos, comm_ferns, 
                             comm_ferns_north, comm_ferns_south))
  ),
  
  func_mntd = target(
    ses_func_mntd(
      comm,
      trait_distance_matrix,
      null.model = "independentswap",
      iterations = 10000,
      runs = 999
    ),
    transform = map(comm = c(comm_pteridos, comm_ferns, 
                             comm_ferns_north, comm_ferns_south))
  ),
  
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
  
  # Ecostructure ----
  
  ### Global analysis (dispersion fields) ###
  
  # Rename pteridophyte species to match World Ferns (same as GBIF)
  occ_data_pteridos_renamed = 
    left_join(
      select(occ_data_pteridos, taxon_id, latitude, longitude), 
      select(green_list, taxon_id, scientific_name)
    ) %>%
    inner_join(
      resolved_names,
      by = c(scientific_name = "query")) %>%
    filter(is.na(genus)) %>%
    assert(not_na, genus) %>%
    assert(not_na, specificEpithet) %>%
    mutate(species = paste(genus, specificEpithet)) %>%
    mutate(site = paste(longitude, latitude, sep = "_")) %>%
    select(species, site),
  
  # Make community data matrix for renamed pteriphytes of Japan
  comm_pteridos_renamed = occ_data_pteridos_renamed %>%
    mutate(
      abundance = 1,
      site = as.character(site)) %>%
    spread(site, abundance) %>%
    mutate_at(vars(-species), ~replace_na(., 0)) %>%
    mutate(species = as.character(species)),
  
  # Crop global species records from GBIF to exclude Japan
  gbif_points_no_japan = exclude_japan_points(gbif_points_global, all_cells),
  
  # Make a global presence/absence matrix (excluding Japan)
  # based on GBFIF occurrence records
  comm_for_ecos_global_cropped = comm_from_points(
    species_coods = gbif_points_no_japan,
    resol = 1,
    rownames = TRUE
  ),
  
  # Make a combined presence/absence matrix
  # 1-degree grid cells outside of Japan,
  # 10-km grid cells inside Japan
  # NEED TO WRITE FUNCTION
  comm_for_ecos_global_combined = combine_presabs_mat(
    comm_for_ecos_global_cropped,
    comm_pteridos_renamed),
  
  # Make dispersion fields list
  dispersion_fields_list = pres_ab_to_disp(
    comm_for_ecos_global_combined
  ),
  
  # Make dispersion fields matrix
  dispersion_fields_matrix = dsp_to_matrix2(
    dispersion_fields_list,
    drop_zero = TRUE # drop all-zero columns, i.e., cells with no species
  ),
  
  # Keep only sites in Japan
  # NEED TO WRITE FUNCTION
  dispersion_fields_matrix_japan = select_disp_fields(
    dispersion_fields_matrix
  ),
  
  # Analyze geographical motifs using ecostructure
  geo_motifs_pteridos = target(
    ecostructure::ecos_fit(
      dispersion_fields_matrix_japan, 
      K = K, tol = 0.1, num_trials = 1),
    transform = map(K = !!k_vals)
  ),
  
  ### Regional analysis (species matrix) ###
  
  # Make input matrix. Species and geographic motifs don't
  # care about phylogeny, so use all pteridophytes together.
  comm_for_ecos_pteridos = make_ecos_matrix(comm_pteridos, all_cells),
  
  # Analyze motifs using ecostructure:
  # - species motifs
  species_motifs_pteridos = target(
    ecostructure::ecos_fit(
      comm_for_ecos_pteridos, 
      K = K, tol = 0.1, num_trials = 1),
    transform = map(K = !!k_vals)
  ),
  
  # - transposed species motifs
  species_motifs_trans_pteridos = target(
    ecostructure::ecos_fit(
      t(comm_for_ecos_pteridos), 
      K = K, tol = 0.1, num_trials = 1),
    transform = map(K = !!k_vals)
  ),
  
  # Write out manuscript ----
  ms = rmarkdown::render(
    knitr_in(here::here("ms/japan_pteridos_biodiv.Rmd")),
    output_file = file_out(here::here("ms/japan_pteridos_biodiv.html")),
    quiet = TRUE)
)
