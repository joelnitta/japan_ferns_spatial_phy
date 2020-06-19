# Set vector of k-values to use for ecostructure
k_vals <- 2:10

# Define analysis plan
plan <- drake_plan (

  # Load and process raw data ----
  
  # Load basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # Unzip data files from Ebihara and Nitta 2019
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  ebihara_2019_data = unzip_ebihara_2019(
    dryad_zip_file = file_in("data_raw/doi_10.5061_dryad.4362p32__v4.zip"), 
    exdir = "data_raw/ebihara_2019",
    produces_1 = file_out("data_raw/ebihara_2019/FernGreenListV1.01E.xls"),
    produces_2 = file_out("data_raw/ebihara_2019/ESM1.csv"),
    produces_3 = file_out("data_raw/ebihara_2019/ESM2.csv"),
    produces_4 = file_out("data_raw/ebihara_2019/japan_pterido_rbcl_cipres.zip"),
    produces_5 = file_out("data_raw/ebihara_2019/2_grid_cells_all.csv"),
    produces_6 = file_out("data_raw/ebihara_2019/ppgi_taxonomy.csv")
  ),
  
  # Lodd Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = read_csv(file_in("data_raw/ebihara_2019/ppgi_taxonomy.csv")) %>%
    modify_ppgi,
  
  # Load Fern Green List (conservation status for each species)
  green_list = read_excel(file_in("data_raw/ebihara_2019/FernGreenListV1.01E.xls")) %>% 
    tidy_japan_names(),

  # Load reproductive mode data (one row per species)
  repro_data_raw = read_csv(
    file_in("data_raw/ebihara_2019/ESM1.csv"),
    col_types = "cccccnnnnn"),

  repro_data = process_repro_data(repro_data_raw),
  
  # Load environmental data
  ja_env_raw_path = target("data/ja_env_data.csv", format = "file"),
  
  ja_env_raw = read_csv(
    ja_env_raw_path,
    col_types = "cnnnnnnnnnn"
  ) %>%
    rename(site = secondary_grid_code),

  # List of all 10 km grid cells across Japan
  all_cells_raw = read_csv(
    file_in("data_raw/ebihara_2019/2_grid_cells_all.csv"),
    col_types = "cnn") %>%
    rename(site = id, latitude = y, longitude = x) %>%
    assert(is_uniq, site),
  
  # Combine into list of 10 km grid cells with environment
  all_cells = left_join(all_cells_raw, ja_env_raw, by = "site"),
  
  # Load raw occurrence data of pteridophytes in Japan.
  # Occurrences are actual point data (one specimen per lat/long)
  # NOTE: this is a subset of all taxa (total 727 taxa, so probably excludes hybrids)
  occ_point_data_raw = readr::read_csv(
    file_in("data_raw/ja_fern_occs_raw.csv"),
    col_types = "lcnncccccc"),
  
  # Subset to just ferns (682 taxa)
  occ_point_data_ferns = subset_occ_point_data(
    occ_point_data_raw = occ_point_data_raw,
    ppgi = ppgi
  ),
  
  # Generate datasets at different spatial scales: 
  comm_04 = phyloregion::points2comm(
    dat = occ_point_data_ferns,
    res = 0.4,
    lon = "longitude",
    lat = "latitude",
    species = "taxon"
  ),
  
  # Load occurrence data of ferns and lycophytes, including hybrids (total 1089 taxa)
  # Occurrences are presences in a set of 10 x 10 km2 grid
  # cells across Japan, not actual occurrence points of specimens.
  occ_data_raw = read_csv(
    file_in("data_raw/ebihara_2019/ESM2.csv"),
    col_types = "cccnnccc"
  ),
  
  # - occurrence data of ferns and lycophytes, including hybrids (total 1089 taxa)
  pteridos = process_occ_data(
    occ_data_raw = occ_data_raw, 
    all_cells = all_cells, 
    ppgi = ppgi),

  # - occurrence data of ferns only, including hybrids (total 1043 taxa)
  ferns =
    pteridos %>%
    assert(not_na, class) %>%
    filter(class == "Polypodiopsida"),
  
  ferns_endemic = left_join(ferns, green_list, by = "taxon_id") %>%
    mutate(endemic = replace_na(endemic, "no")) %>%
    filter(endemic == "Endemic"),
  
  pteridos_endemic = left_join(pteridos, green_list, by = "taxon_id") %>%
    mutate(endemic = replace_na(endemic, "no")) %>%
    filter(endemic == "Endemic"),

  # Make tibble mapping taxon IDs to taxon names
  taxon_id_map = make_taxon_id_map(pteridos),

  # Read in phylogenetic tree of all pteridophyte
  # taxa based on rbcL gene, NOT including hybrids (706 taxa total)
  japan_pterido_tree = read_nexus_in_zip(
    file_in("data_raw/ebihara_2019/japan_pterido_rbcl_cipres.zip"),
    "japan_pterido_rbcl_cipres/infile.nex.con.tre")[[2]] %>%
    format_tip_labels %>%
    rename_tree(taxon_id_map),
  
  # Make community matrices, limited to taxa in tree (so, excluding hybrids).
  # Rows are sites, columns are species.
  comm = target(
    make_comm_matrix(occ_data, taxa_list = japan_pterido_tree$tip.label),
    transform = map(occ_data = c(pteridos, ferns, pteridos_endemic, ferns_endemic))
  ),

  # Format trait data
  raw_trait_data_path = target("data_raw/JpFernLUCID_forJoel.xlsx", format = "file"),
  
  traits_for_dist = format_traits(
    path_to_lucid_traits = raw_trait_data_path,
    taxon_id_map = taxon_id_map),

  # Make trait distance matrix using taxon IDs as labels
  trait_distance_matrix = make_trait_dist_matrix(traits_for_dist),
  
  # Analyze reproductive mode ----
  
  percent_sex_dip = target(
    calc_sex_dip(comm, repro_data = repro_data, taxon_id_map = taxon_id_map),
    transform = map(comm, .id = occ_data)
  ),

  # Analyze standard effect size (SES) of diversity metrics ----
  
  ses_div = target(
    run_ses_analysis(comm, japan_pterido_tree, n_reps = 999, metrics = c("mpd", "pd", "pe", "rpe")),
    transform = map(comm, .id = occ_data)
  )
  
  # 
  # # Make richness matrix (number of species per
  # # 10 km grid cell).
  # richness_pteridos = make_richness_matrix(occ_data_pteridos),
  # 
  # richness_ferns = make_richness_matrix(occ_data_ferns),
  # 
  # # Make community matrix (presence/absence of each species in
  # # 10 km grid cells), trim to only species in tree.
  # comm_pteridos = make_comm_matrix(occ_data_pteridos) %>%
  #   match_comm_and_tree(japan_pterido_tree, "comm"),
  # 
  # comm_ferns = make_comm_matrix(occ_data_ferns) %>%
  #   match_comm_and_tree(japan_pterido_tree, "comm"),
  # 
  # comm_ferns_north = make_comm_matrix(occ_data_ferns_north) %>%
  #   match_comm_and_tree(japan_pterido_tree, "comm"),
  # 
  # comm_ferns_south = make_comm_matrix(occ_data_ferns_south) %>%
  #   match_comm_and_tree(japan_pterido_tree, "comm"),
  # 
  # ### Calculate phylogenetic diversity ###
  # phy_mpd = target(
  #   ses_phy_mpd(
  #     comm,
  #     japan_pterido_tree,
  #     null.model = "independentswap",
  #     iterations = 10000,
  #     runs = 999
  #   ),
  #   transform = map(comm = c(comm_pteridos, comm_ferns,
  #                            comm_ferns_north, comm_ferns_south))
  # ),
  # 
  # phy_mntd = target(
  #   ses_phy_mntd(
  #     comm,
  #     japan_pterido_tree,
  #     null.model = "independentswap",
  #     iterations = 10000,
  #     runs = 999
  #   ),
  #   transform = map(comm = c(comm_pteridos, comm_ferns,
  #                            comm_ferns_north, comm_ferns_south))
  # ),
  # 
  # ### Analyze functional trait diversity ###
  # 
  # # - format trait data
  # traits_for_dist = format_traits(
  #   file_in("data_raw/JpFernLUCID_forJoel.xlsx"),
  #   taxon_id_map),
  # 
  # # - make trait distance matrix using taxon IDs as labels
  # trait_distance_matrix = make_trait_dist_matrix(
  #   traits_for_dist),
  # 
  # func_mpd = target(
  #   ses_func_mpd(
  #     comm,
  #     trait_distance_matrix,
  #     null.model = "independentswap",
  #     iterations = 10000,
  #     runs = 999
  #   ),
  #   transform = map(comm = c(comm_pteridos, comm_ferns,
  #                            comm_ferns_north, comm_ferns_south))
  # ),
  # 
  # func_mntd = target(
  #   ses_func_mntd(
  #     comm,
  #     trait_distance_matrix,
  #     null.model = "independentswap",
  #     iterations = 10000,
  #     runs = 999
  #   ),
  #   transform = map(comm = c(comm_pteridos, comm_ferns,
  #                            comm_ferns_north, comm_ferns_south))
  # ),
  # 
  # ### Combine alpha diversity metrics ###
  # 
  # # First combine metrics that were calculated separately for
  # # ferns in north and south back together
  # phy_mpd_comm_ferns_ns = rbind(phy_mpd_comm_ferns_north, phy_mpd_comm_ferns_south),
  # 
  # phy_mntd_comm_ferns_ns = rbind(phy_mntd_comm_ferns_north, phy_mntd_comm_ferns_south),
  # 
  # func_mpd_comm_ferns_ns = rbind(func_mpd_comm_ferns_north, func_mpd_comm_ferns_south),
  # 
  # func_mntd_comm_ferns_ns = rbind(func_mntd_comm_ferns_north, func_mntd_comm_ferns_south),
  # 
  # # Combine all measures of alpha diversity by data set.
  # # Start with all_cells so that every grid cell is included,
  # # even if there are 0 species there.
  # alpha_div_pteridos =
  #   all_cells %>%
  #   left_join(select(richness_pteridos, site, richness)) %>%
  #   left_join(clean_ses(phy_mpd_comm_pteridos)) %>%
  #   left_join(clean_ses(phy_mntd_comm_pteridos)) %>%
  #   left_join(clean_ses(func_mpd_comm_pteridos, "func_")) %>%
  #   left_join(clean_ses(func_mntd_comm_pteridos, "func_")) %>%
  #   left_join(select(percent_sex_dip_pteridos, site, percent_sex_dip)) %>%
  #   mutate(richness = replace_na(richness, 0)),
  # 
  # alpha_div_ferns =
  #   all_cells %>%
  #   left_join(select(richness_ferns, site, richness)) %>%
  #   left_join(clean_ses(phy_mpd_comm_ferns)) %>%
  #   left_join(clean_ses(phy_mntd_comm_ferns)) %>%
  #   left_join(clean_ses(func_mpd_comm_ferns, "func_")) %>%
  #   left_join(clean_ses(func_mntd_comm_ferns, "func_")) %>%
  #   left_join(select(percent_sex_dip_ferns, site, percent_sex_dip)) %>%
  #   mutate(richness = replace_na(richness, 0)),
  # 
  # alpha_div_ferns_ns =
  #   all_cells %>%
  #   left_join(select(richness_ferns, site, richness)) %>%
  #   left_join(clean_ses(phy_mpd_comm_ferns_ns)) %>%
  #   left_join(clean_ses(phy_mntd_comm_ferns_ns)) %>%
  #   left_join(clean_ses(func_mpd_comm_ferns_ns, "func_")) %>%
  #   left_join(clean_ses(func_mntd_comm_ferns_ns, "func_")) %>%
  #   left_join(select(percent_sex_dip_ferns, site, percent_sex_dip)) %>%
  #   mutate(richness = replace_na(richness, 0))
  # ,
  # 
  # # Extra analysis ----
  # 
  # # Haven't decided yet to include these in the MS
  # 
  # ### Traits ###
  # # Run NMDS on traits
  # traits_nmds = vegan::metaMDS(trait_distance_matrix),
  # 
  # # Plot species in NMDS trait space
  # # traits_nmds_plot = make_trait_nmds_plot(traits_nmds, ppgi, taxon_id_map),
  # 
  # # Write out NMDS plot
  # # traits_nmds_plot_out = ggsave(
  # #   plot = traits_nmds_plot,
  # #   file = "results/traits_nmds_plot.pdf",
  # #   width = 10,
  # #   height = 6,
  # #   units = "in"),
  # 
  # # Plot and write out traits dendrogram
  # # traits_dendrogram = make_traits_dendrogram(trait_distance_matrix, taxon_id_map),
  # 
  # ### Diversity by elevation + latitude scatter plots ###
  # # pterido_lat_el_plot = compose_lat_el_plots(alpha_div_pteridos, "Pteridophytes") %>%
  # #   ggsave(
  # #     plot = .,
  # #     filename = "results/pterido_lat_el_plot.pdf",
  # #     height = 9, width = 7, units = "in"),
  # # 
  # # ferns_lat_el_plot = compose_lat_el_plots(alpha_div_ferns, "Ferns") %>%
  # #   ggsave(
  # #   plot = .,
  # #   filename = "results/ferns_lat_el_plot.pdf",
  # #   height = 9, width = 7, units = "in"),
  # # 
  # # ferns_ns_lat_el_plot = compose_lat_el_plots(alpha_div_ferns_ns, "Ferns_NS") %>%
  # #   ggsave(
  # #   plot = .,
  # #   filename = "results/ferns_ns_lat_el_plot.pdf",
  # #   height = 9, width = 7, units = "in"),
  # 
  # # Ecostructure ----
  # 
  # ### Regional analysis (species matrix) ###
  # 
  # # Make input matrices for ecostructure
  # comm_for_ecos_pteridos = make_ecos_matrix(comm_pteridos, all_cells),
  # 
  # # remove communities with zero richness (after filtering to ferns)
  # comm_ferns_filtered =
  #   comm_ferns %>%
  #   pivot_longer(-species, values_to = "abun", names_to = "site") %>%
  #   group_by(site) %>%
  #   filter(sum(abun) > 0) %>%
  #   ungroup %>%
  #   pivot_wider(names_from = site, values_from = abun),
  # 
  # comm_for_ecos_ferns = make_ecos_matrix(comm_ferns_filtered, all_cells),
  # 
  # # Make dispersion files
  # # - dispersion fields list
  # dispersion_fields_ferns = pres_ab_to_disp(
  #   pres_ab_matrix = comm_for_ecos_ferns,
  #   xmin = pull(occ_data_pteridos, longitude) %>% min %>% floor, 
  #   xmax = pull(occ_data_pteridos, longitude) %>% max %>% ceiling,
  #   ymin = pull(occ_data_pteridos, latitude) %>% min %>% floor, 
  #   ymax = pull(occ_data_pteridos, latitude) %>% max %>% ceiling,
  #   res = 0.5
  # ),
  # 
  # # - dispersion fields matrix
  # dispersion_fields_matrix_ferns = dsp_to_matrix2(
  #   dispersion_fields_ferns,
  #   drop_zero = TRUE # drop all-zero columns, i.e., cells with no species
  # ),
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
  # species_motifs_ferns = target(
  #   ecostructure::ecos_fit(
  #     comm_for_ecos_ferns,
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
  # 
  # # dispersion fields
  # geo_motifs_ferns = target(
  #   ecostructure::ecos_fit(
  #     dispersion_fields_matrix_ferns,
  #     K = K, tol = 0.1, num_trials = 1),
  #   transform = map(K = !!k_vals)
  # ),
  # 
  # # Infomap Bioregions ----
  # 
  # # Prep data for Infomap Bioregions 
  # # Filter raw point occurrence data to ferns
  # occ_point_data_ferns =
  #   occ_point_data_raw %>%
  #   select(taxon_name = species, decimalLongitude, decimalLatitude) %>%
  #   add_taxonomy(ppgi) %>%
  #   assert(not_na, class) %>%
  #   filter(class == "Polypodiopsida"),
  # 
  # # Write data out for Infomap Bioregions
  # occ_point_data_ferns_out = readr::write_csv(
  #   occ_point_data_ferns, "data/ja_fern_occ.csv"),
  # 
  # # Run Infomap Bioregions via online server
  # #
  # # https://bioregions.mapequation.org/
  # # 
  # # Analysis settings:
  # # max cell size = 1
  # # min cell size = 0.5
  # # max cell capacity = 100
  # # min cell capacity = 10
  # # patch sparse grid cells = yes
  # # number of trials = 10
  # # number of cluster cost = 1
  # #
  # # Export resulting geojson file to "data" as "ja_fern_occ_bioregions.geojson"
  # 
  # # Read in results of running Infomap Bioregions externally
  # ja_bioregions_ferns = sf::st_read(file_in("data/ja_fern_occ_bioregions.geojson")) %>%
  #   # Join geometries within each bioregion
  #   mutate(geometry = sf::st_union(geometry, by_feature = TRUE)),
  # 
  # # Biodiverse ----
  # 
  # # Prep data for Biodiverse
  # comm_ferns_renamed = rename_comm(comm_ferns, taxon_id_map),
  # 
  # matrix_for_biodiverse_ferns = make_matrix_for_biodiverse(
  #   comm_ferns_renamed, all_cells_raw),
  # 
  # matrix_for_biodiverse_ferns_out = readr::write_csv(
  #   matrix_for_biodiverse_ferns, 
  #   file_out("data/matrix_for_biodiverse_ferns.csv")
  #   ),
  # 
  # tree_for_biodiverse_ferns = match_comm_and_tree(
  #   comm_ferns, japan_pterido_tree, "tree") %>%
  #   rename_tree(taxon_id_map),
  # 
  # tree_for_biodiverse_ferns_out = ape::write.tree(
  #   tree_for_biodiverse_ferns,
  #   file_out("data/tree_for_biodiverse_ferns.tre")
  # ),
  # 
  # # Also do for endemics only
  # comm_ferns_endemic = left_join(comm_ferns, green_list, by = c(species = "taxon_id")) %>%
  #   mutate(endemic = replace_na(endemic, "no")) %>%
  #   filter(endemic == "Endemic") %>%
  #   select(-endemic, -scientific_name, -conservation_status) %>%
  #   # Drop sites with 0 abundance
  #   pivot_longer(-species, names_to = "site", values_to = "abun") %>%
  #   filter(abun > 0) %>%
  #   pivot_wider(names_from = "site", values_from = "abun") %>%
  #   mutate_at(vars(-species), ~replace_na(., 0)),
  # 
  # matrix_for_biodiverse_ferns_endemic = make_matrix_for_biodiverse(
  #   comm_ferns_endemic, all_cells_raw),
  # 
  # matrix_for_biodiverse_ferns_endemic_out = readr::write_csv(
  #   matrix_for_biodiverse_ferns_endemic, 
  #   file_out("data/matrix_for_biodiverse_ferns_endemic.csv")
  # ),
  # 
  # tree_for_biodiverse_ferns_endemic = match_comm_and_tree(
  #   comm_ferns_endemic, japan_pterido_tree, "tree"),
  # 
  # tree_for_biodiverse_ferns_endemic_out = ape::write.tree(
  #   tree_for_biodiverse_ferns_endemic,
  #   file_out("data/tree_for_biodiverse_ferns_endemic.tre")
  # ),
  #
  # # # Write out manuscript ----
  # ms = rmarkdown::render(
  #   knitr_in(here::here("ms/manuscript.Rmd")),
  #   output_file = file_out(here::here("ms/manuscript.pdf")),
  #   quiet = TRUE)
)
