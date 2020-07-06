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
  
  # Calculate richness, abundance, and redundancy at four scales: 
  # 0.1, 0.2, 0.3, and 0.4 degrees
  comm_scaled_list = target(
    comm_from_points(
      species_coods = occ_point_data_ferns,
      res = scale,
      lon = "longitude",
      lat = "latitude",
      species = "taxon"),
    transform = map(scale = c(0.1, 0.2, 0.3, 0.4))
  ),
  
  # Decide that 0.2 scale is optimal, use this for downstream analyses
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
  
)
