# Set vector of k-values to use for ecostructure
k_vals <- 2:20

# Define analysis plan
plan <- drake_plan (

  # Load and process data ----
  
  # Load basic world map.
  world_map = ggplot2::map_data("world") %>%
    rename(longitude = long, latitude = lat),
  
  # Unzip data files from Ebihara and Nitta 2019
  # This requires doi_10.5061_dryad.4362p32__v4.zip to be downloaded to data_raw/
  # from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32 first
  ebihara_2019_data = unzip_ebihara_2019(
    dryad_zip_file = file_in("data_raw/doi_10.5061_dryad.4362p32__v4.zip"), 
    exdir = "data_raw/ebihara_2019",
    produces_1 = file_out("data_raw/ebihara_2019/FernGreenListV1.01E.xls"), # Green List
    produces_2 = file_out("data_raw/ebihara_2019/ESM1.csv"), # Reproductive mode data
    produces_3 = file_out("data_raw/ebihara_2019/japan_pterido_rbcl_cipres.zip"), # rbcL tree
    produces_4 = file_out("data_raw/ebihara_2019/ppgi_taxonomy.csv") # PPGI taxonomy
  ),
  
  # Load Pteridophyte Phylogeny Group I (PPGI) taxonomy,
  # modified slightly for ferns of Japan
  ppgi = read_csv(file_in("data_raw/ebihara_2019/ppgi_taxonomy.csv")) %>%
    modify_ppgi,
  
  # Load Fern Green List (official taxonomy + conservation status for each species)
  green_list = read_excel(file_in("data_raw/ebihara_2019/FernGreenListV1.01E.xls")) %>% 
    tidy_japan_names(),
  
  # Load reproductive mode data (one row per species)
  repro_data_raw = read_csv(
    file_in("data_raw/ebihara_2019/ESM1.csv"),
    col_types = "cccccnnnnn"),

  repro_data = process_repro_data(repro_data_raw) %>%
    rename_taxa(green_list),
  
  # Load raw occurrence data of pteridophytes in Japan, excluding hybrids (717 taxa)
  occ_point_data_raw = readxl::read_excel(
    file_in("data_raw/JP_pterid_excl_hyb200620.xlsx"),
    col_types = c("text", "numeric", "numeric", "text", "text", "text", "text"),
    col_names = c("species", "longitude", "latitude", "date", "tns_barcode", "herbarium_code", "taxon_id"),
    skip = 1),

  # - standardize names to Green List
  occ_point_data = rename_taxa(occ_point_data_raw, green_list),
  
  # - subset to just ferns (674 taxa)
  occ_point_data_ferns = subset_to_ferns(occ_point_data, ppgi) %>%
    # check for missing data
    assert(not_na, longitude, latitude, taxon),
  
  # Calculate richness, abundance, and redundancy at three scales: 
  # 0.1, 0.2, 0.3, and 0.4 degrees
  richness_by_res = map_df(
    c("0.1" = 0.1, "0.2" = 0.2, "0.3" = 0.3, "0.4" = 0.4), 
    ~richness_from_points(occ_point_data_ferns, .), .id = "res"),
  
  # Decide that 0.2 scale is optimal, use this for downstream analyses
  comm_scaled_list_0.2 = phyloregion::points2comm(
    dat = occ_point_data_ferns,
    res = 0.2,
    lon = "longitude",
    lat = "latitude",
    species = "taxon"),
  
  # - extract community matrix
  comm_ferns = comm_from_points2comm(comm_scaled_list_0.2),
  
  # - extract geographic shapes, richness, and number of specimens
  shape_ferns = shape_from_points2comm(comm_scaled_list_0.2) %>%
    # calculate redundancy
    mutate(redundancy = 1 - (richness/abundance)),
    
  # - make community matrix subset to taxa endemic to Japan
  comm_ferns_endemic = subset_comm_to_endemic(
    comm = comm_ferns,
    green_list = green_list
  ),
  
  # Read in ultrametric phylogenetic tree of all pteridophytes,
  # not including hybrids (706 taxa total)
  japan_pterido_tree_path = target("data_raw/japan_pterido_tree_dated.tre", format = "file"),
  
  japan_pterido_tree = ape::read.tree(japan_pterido_tree_path),
  
  # - subset to only ferns
  japan_fern_tree = subset_tree(
    phy = japan_pterido_tree, 
    ppgi = ppgi),
  
  # Format trait data, subset to ferns
  raw_trait_data_path = target("data_raw/JpFernLUCID_forJoel.xlsx", format = "file"),
  
  traits_for_dist = format_traits(
    path_to_lucid_traits = raw_trait_data_path,
    taxon_id_map = green_list) %>%
    subset_to_ferns(ppgi),

  # Make trait distance matrix using taxon IDs as labels
  trait_distance_matrix = make_trait_dist_matrix(traits_for_dist),
  
  # Analyze reproductive mode ----
  
  percent_sex_dip = calc_sex_dip(
    comm = comm_ferns, 
    repro_data = repro_data),

  # Analyze standard effect size (SES) of diversity metrics ----
  
  ses_phy = target(
    run_ses_analysis(
      comm = comm, 
      phy = japan_fern_tree, 
      n_reps = 999, 
      metrics = metrics) %>% 
      categorize_endemism,
    transform = map(
      comm = c(comm_ferns, comm_ferns_endemic),
      metrics = c(
        c("pd", "rpd", "pe", "rpe"),
        c("pe", "rpe")
      ),
      .names = c("ses_phy_ferns", "ses_phy_ferns_endemic"))
  ),

  ses_traits_ferns = run_ses_analysis(
    comm = comm_ferns,
    n_reps = 999,
    metrics = c("fd", "rfd"),
    trait_distances = trait_distance_matrix),

  # Analyze phyloregions

  # - Assess optimal K-value for clustering by taxonomy
  # (plot this, then choose K manually)
  k_taxonomy = find_k_taxonomy(comm_ferns),

  # - Assess optimal K-value for clustering by phylogeny
  # (plot this, then choose K manually)
  k_phylogeny = find_k_phylogeny(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
  ),

  # - Cluster by taxonomy
  regions_taxonomy = cluster_taxonomic_regions(
    comm_df = comm_ferns,
    k = 8),

  # - Cluster by phylogeny
  regions_phylogeny = cluster_phylo_regions(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
    k = 6
  ),

  # Ecostructure ----

  species_motifs_ferns = target(
      ecostructure::ecos_fit(
        comm_ferns,
        K = K, tol = 0.1, num_trials = 1),
      transform = map(K = !!k_vals)
    ),
  
  species_motifs_ferns_combined = target(
    list(species_motifs_ferns) %>% set_names(k_vals),
    transform = combine(species_motifs_ferns)
  ),
  
  # Traits ----
  
  # Run NMDS on traits
  traits_nmds = vegan::metaMDS(trait_distance_matrix, k = 3),
  
  # Write out manuscript ----
  ms = rmarkdown::render(
    knitr_in("ms/manuscript.Rmd"),
    output_dir = here::here("results"),
    quiet = TRUE),
   
  si = rmarkdown::render(
    knitr_in("ms/SI.Rmd"),
    output_dir = here::here("results"),
    quiet = TRUE)

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
