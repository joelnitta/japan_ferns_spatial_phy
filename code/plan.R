# Set vector of k-values to use for ecostructure
k_vals <- 2:10

# Define analysis plan
plan <- drake_plan (

  # Load and process raw data ----

  ### Taxonomy data ###

  # Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = read_csv(file_in("data_raw/ppgi_taxonomy.csv")) %>%
    modify_ppgi,

  # Catalog of Life (COL) plants taxonomic data
  col_plants = read_col_plants(file_in("data_raw/archive-kingdom-plantae-phylum-tracheophyta-bl3/taxa.txt")),
  
  # Unzip data files from Ebihara and Nitta 2019
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  ebihara_2019_data = unzip_ebihara_2019(
    dryad_zip_file = file_in("data/doi_10.5061_dryad.4362p32__v4.zip"), 
    exdir = "data/ebihara_2019",
    produces_1 = file_out("data/ebihara_2019/FernGreenListV1.01E.xls"),
    produces_2 = file_out("data/ebihara_2019/ESM1.csv"),
    produces_3 = file_out("data/ebihara_2019/ESM2.csv"),
    produces_4 = file_out("data/ebihara_2019/japan_pterido_rbcl_cipres.zip"),
    produces_5 = file_out("data/ebihara_2019/2_grid_cells_all.csv")
  ),

  # Load Fern Green List, with conservation status for each species.
  green_list = read_excel(file_in("data/ebihara_2019/FernGreenListV1.01E.xls")) %>% 
    tidy_japan_names(),

  # Match fern and pteridophyte names to COL.
  resolved_names_auto = taxastand::resolve_fern_names(green_list$scientific_name, col_plants, resolve_to = "species"),

  # Read in list of manually resolved names.
  resolved_names_manual_fix = read_csv("data_raw/japan_pterido_names_manual_fix.csv") %>%
    filter(!is.na(scientificName)),

  # Update automatically resolved names with manually resolved names to get final list
  resolved_names = update_resolved_names(
    resolved_names_auto,
    resolved_names_manual_fix
  ),

  # Reproductive mode data, with one row per species.
  repro_data_raw = read_csv(
    file_in("data/ebihara_2019/ESM1.csv"),
    col_types = "cccccnnnnn"),

  repro_data = process_repro_data(repro_data_raw),
  
  # Raw environmental data
  ja_env_raw = read_csv(
    file_in("data/ja_env_data.csv"),
    col_types = "nnnnnnnnnnn"
  ),

  # List of all 10 km grid cells across Japan
  all_cells_raw = read_csv(
    file_in("data/ebihara_2019/2_grid_cells_all.csv"),
    col_types = "nnn") %>%
    rename(secondary_grid_code = id, latitude = y, longitude = x) %>%
    assert(is_uniq, secondary_grid_code),
  
  # Combine into list of 10km grid cells with environment
  all_cells = left_join(all_cells_raw, ja_env_raw, by = "secondary_grid_code"),

  # Occurrence data, with multiple rows per species.
  # Occurrences are presences in a set of 1km2 grid
  # cells across Japan, not actual occurrence points of specimens.
  occ_data_raw = read_csv(
    file_in("data/ebihara_2019/ESM2.csv"),
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
  # 10 km grid cell).
  richness_pteridos = make_richness_matrix(occ_data_pteridos),

  richness_ferns = make_richness_matrix(occ_data_ferns),

  # Make community matrix (presence/absence of each species in
  # 10 km grid cells), trim to only species in tree.
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

  # - format trait data
  traits_for_dist = format_traits(
    file_in("data_raw/JpFernLUCID_forJoel.xlsx"),
    taxon_id_map),

  # - make trait distance matrix using taxon IDs as labels
  trait_distance_matrix = make_trait_dist_matrix(
    traits_for_dist),

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

  ### Combine alpha diversity metrics ###

  # First combine metrics that were calculated separately for
  # ferns in north and south back together
  phy_mpd_comm_ferns_ns = rbind(phy_mpd_comm_ferns_north, phy_mpd_comm_ferns_south),

  phy_mntd_comm_ferns_ns = rbind(phy_mntd_comm_ferns_north, phy_mntd_comm_ferns_south),

  func_mpd_comm_ferns_ns = rbind(func_mpd_comm_ferns_north, func_mpd_comm_ferns_south),

  func_mntd_comm_ferns_ns = rbind(func_mntd_comm_ferns_north, func_mntd_comm_ferns_south),

  # Combine all measures of alpha diversity by data set.
  # Start with all_cells so that every grid cell is included,
  # even if there are 0 species there.
  alpha_div_pteridos =
    all_cells %>%
    left_join(select(richness_pteridos, secondary_grid_code, richness)) %>%
    left_join(clean_ses(phy_mpd_comm_pteridos)) %>%
    left_join(clean_ses(phy_mntd_comm_pteridos)) %>%
    left_join(clean_ses(func_mpd_comm_pteridos, "func_")) %>%
    left_join(clean_ses(func_mntd_comm_pteridos, "func_")) %>%
    left_join(select(percent_sex_dip_pteridos, secondary_grid_code, percent_sex_dip)) %>%
    mutate(richness = replace_na(richness, 0)),

  alpha_div_ferns =
    all_cells %>%
    left_join(select(richness_ferns, secondary_grid_code, richness)) %>%
    left_join(clean_ses(phy_mpd_comm_ferns)) %>%
    left_join(clean_ses(phy_mntd_comm_ferns)) %>%
    left_join(clean_ses(func_mpd_comm_ferns, "func_")) %>%
    left_join(clean_ses(func_mntd_comm_ferns, "func_")) %>%
    left_join(select(percent_sex_dip_ferns, secondary_grid_code, percent_sex_dip)) %>%
    mutate(richness = replace_na(richness, 0)),

  alpha_div_ferns_ns =
    all_cells %>%
    left_join(select(richness_ferns, secondary_grid_code, richness)) %>%
    left_join(clean_ses(phy_mpd_comm_ferns_ns)) %>%
    left_join(clean_ses(phy_mntd_comm_ferns_ns)) %>%
    left_join(clean_ses(func_mpd_comm_ferns_ns, "func_")) %>%
    left_join(clean_ses(func_mntd_comm_ferns_ns, "func_")) %>%
    left_join(select(percent_sex_dip_ferns, secondary_grid_code, percent_sex_dip)) %>%
    mutate(richness = replace_na(richness, 0))
  ,

  # Extra analysis ----

  # Haven't decided yet to include these in the MS

  ### Traits ###
  # Run NMDS on traits
  traits_nmds = vegan::metaMDS(trait_distance_matrix),

  # Plot species in NMDS trait space
  traits_nmds_plot = make_trait_nmds_plot(traits_nmds, ppgi, taxon_id_map),

  # Write out NMDS plot
  traits_nmds_plot_out = ggsave(
    plot = traits_nmds_plot,
    file = "results/traits_nmds_plot.pdf",
    width = 10,
    height = 6,
    units = "in"),

  # Plot and write out traits dendrogram
  traits_dendrogram = make_traits_dendrogram(trait_distance_matrix, taxon_id_map),

  ### Diversity by elevation + latitude scatter plots ###
  pterido_lat_el_plot = compose_lat_el_plots(alpha_div_pteridos, "Pteridophytes") %>%
    ggsave(
      plot = .,
      filename = "results/pterido_lat_el_plot.pdf",
      height = 9, width = 7, units = "in"),

  ferns_lat_el_plot = compose_lat_el_plots(alpha_div_ferns, "Ferns") %>%
    ggsave(
    plot = .,
    filename = "results/ferns_lat_el_plot.pdf",
    height = 9, width = 7, units = "in"),

  ferns_ns_lat_el_plot = compose_lat_el_plots(alpha_div_ferns_ns, "Ferns_NS") %>%
    ggsave(
    plot = .,
    filename = "results/ferns_ns_lat_el_plot.pdf",
    height = 9, width = 7, units = "in"),

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
    mutate(site = paste(longitude, latitude, sep = "_")) %>%
    select(species, site),

  # # Make community data matrix for renamed pteriphytes of Japan
  comm_pteridos_renamed = 
    occ_data_pteridos_renamed %>%
    # Have to use unique occurrences since some original names got 
    # collapsed to same resolved name (synonyms or varieties)
    unique() %>%
    mutate(species = str_replace_all(species, " ", "_")) %>%
    mutate(abundance = 1) %>%
    pivot_wider(names_from = species, values_from = abundance) %>%
    mutate_at(vars(-site), ~replace_na(., 0)),

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
  # 10-km grid cells inside Japan, 
  # limited to only species in common in both
  comm_for_ecos_global_combined = combine_presabs_mat(
    comm_for_ecos_global_cropped,
    comm_pteridos_renamed),
  
  # Split into datasets: 
  # combined at 1-degree resolution, 
  # global 1-degree plus JA 10km
  comm_ja = magrittr::extract(
    comm_for_ecos_global_combined, 
    rownames(comm_for_ecos_global_combined) %in% comm_pteridos_renamed$site, ),
  
  comm_global = magrittr::extract(
    comm_for_ecos_global_combined, 
    !rownames(comm_for_ecos_global_combined) %in% comm_pteridos_renamed$site, ),
  
  # Make dispersion fields list
  # - Combined global + Japan 1-degree
  dispersion_fields_list = pres_ab_to_disp(comm_for_ecos_global_combined),
  
  # - Global 1-degree
  dispersion_fields_list_global = pres_ab_to_disp(comm_global),
  
  # - Japan 10 km
  dispersion_fields_list_ja = pres_ab_to_disp(comm_ja, res = 1/111), # set resolution to ca. 10 km

  # Make dispersion fields matrix
  # - Combined global + Japan 1-degree
  dispersion_fields_matrix = dsp_to_matrix2(
    dispersion_fields_list,
    drop_zero = TRUE # drop all-zero columns, i.e., cells with no species
  ),
  
  # - Global (no Japan) Japan 1-degree
  dispersion_fields_matrix_global = dsp_to_matrix2(
    dispersion_fields_list_global,
    drop_zero = TRUE
  ),
  
  # - Japan only 10 km
  dispersion_fields_matrix_ja = dsp_to_matrix2(
    dispersion_fields_list_ja,
    drop_zero = TRUE
  ),

  # Keep only sites in Japan
  dispersion_fields_matrix_japan = dispersion_fields_matrix[rownames(dispersion_fields_matrix) %in% comm_pteridos_renamed$site,],
  
  # Write out dispersion fields matrix for running ecostructure analysis as separate plan
  dispersion_fields_matrix_japan_out = saveRDS(dispersion_fields_matrix_japan, "data/dispersion_fields_matrix_japan.RDS"),

  # Analyze geographical motifs using ecostructure
  geo_motifs_pteridos = target(
    ecostructure::ecos_fit(
      dispersion_fields_matrix_japan,
      K = K, tol = 0.1, num_trials = 1),
    transform = map(K = !!k_vals)
  )

  # ### Regional analysis (species matrix) ###
  #
  # # Make input matrix. Species and geographic motifs don't
  # # care about phylogeny, so use all pteridophytes together.
  # comm_for_ecos_pteridos = make_ecos_matrix(comm_pteridos, all_cells),
  #
  # # Analyze motifs using ecostructure:
  # # - species motifs
  # species_motifs_pteridos = target(
  #   ecostructure::ecos_fit(
  #     comm_for_ecos_pteridos,
  #     K = K, tol = 0.1, num_trials = 1),
  #   transform = map(K = !!k_vals)
  # ),
  #
  # # - transposed species motifs
  # species_motifs_trans_pteridos = target(
  #   ecostructure::ecos_fit(
  #     t(comm_for_ecos_pteridos),
  #     K = K, tol = 0.1, num_trials = 1),
  #   transform = map(K = !!k_vals)
  # ),

  # # Write out manuscript ----
  # ms = rmarkdown::render(
  #   knitr_in(here::here("ms/manuscript.Rmd")),
  #   output_file = file_out(here::here("ms/manuscript.pdf")),
  #   quiet = TRUE)
)
