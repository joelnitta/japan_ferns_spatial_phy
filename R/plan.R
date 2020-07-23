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
  occ_point_data = rename_taxa(occ_point_data_raw, green_list) %>%
    # check for missing data
    assert(not_na, longitude, latitude, taxon),
  
  # - subset to just ferns (674 taxa)
  occ_point_data_ferns_unfiltered = subset_to_ferns(occ_point_data, ppgi),
  
  # - filter out duplicates, restrict to only points in second-degree mesh
  # Shape file downloaded from http://gis.biodic.go.jp/
  # http://gis.biodic.go.jp/BiodicWebGIS/Questionnaires?kind=mesh2&filename=mesh2.zip
  occ_point_data_ferns = filter_occ_points(
    occ_point_data = occ_point_data_ferns_unfiltered,
    shape_file = file_in("data_raw/mesh2/mesh2.shp")),
  
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
  
  # - extract geographic shapes, richness, and number of specimens
  shape_ferns_full = shape_from_points2comm(comm_scaled_list_0.2) %>%
    # calculate redundancy
    mutate(redundancy = 1 - (richness/abundance)),
  
  # - extract community matrix
  comm_ferns_full = comm_from_points2comm(comm_scaled_list_0.2),
  
  # - subset geographic shapes to redundancy > 0.1
  shape_ferns = filter(shape_ferns_full, redundancy > 0.1),
  
  # - subset community matrix to communities with redundancy > 0.1
  comm_ferns = filter_comm_by_redun(
    comm = comm_ferns_full,
    shape = shape_ferns,
    cutoff = 0.1
  ),
  
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
      null_model = "independentswap",
      n_reps = 999,
      n_iterations = 100000,
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
    null_model = "independentswap",
    n_reps = 999,
    n_iterations = 100000,
    metrics = c("fd", "rfd"),
    trait_distances = trait_distance_matrix),
  
  # Analyze bioregions ----
  
  # - Assess optimal K-value for clustering by taxonomy
  k_taxonomy = find_k_taxonomy(comm_ferns),
  
  # - Assess optimal K-value for clustering by phylogeny
  k_phylogeny = find_k_phylogeny(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
  ),
  
  # - Cluster by taxonomy
  regions_taxonomy = cluster_taxonomic_regions(
    comm_df = comm_ferns,
    k = k_taxonomy[["optimal"]][["k"]]
  ),
  
  # - Cluster by phylogeny
  regions_phylogeny = cluster_phylo_regions(
    comm_df = comm_ferns,
    phy = japan_fern_tree,
    k = k_phylogeny[["optimal"]][["k"]]
  ),
  
  # Traits ----
  
  # Run NMDS on traits
  traits_nmds = vegan::metaMDS(trait_distance_matrix, k = 3),
  
  # Combine results ----
  
  # Combine spatial data, alpha diversity, and regions
  
  # - All ferns
  ses_div_ferns_spatial =
    shape_ferns %>%
    left_join(ses_phy_ferns, by = c(grids = "site")) %>%
    left_join(ses_traits_ferns, by = c(grids = "site")) %>%
    left_join(regions_taxonomy %>% rename(taxonomic_cluster = cluster), by = "grids") %>%
    left_join(regions_phylogeny %>% rename(phylo_cluster = cluster), by = "grids"),
  
  # - Japan endemics only
  ses_div_ferns_endemic_spatial =
    shape_ferns %>%
    left_join(ses_phy_ferns_endemic, by = c(grids = "site")),
  
  # Conservation analysis ----
  
  ## Read in protected areas (7 separate shape files corresponding to different kinds of areas)
  # Assign protection levels following Kusamoto et al. 2017
  # - high: no human activities allowed
  # - medium: permission required for economic activities
  # - low: protected area, but none of the above restrictions
  
  # 1: wilderness
  protected_1 = sf::st_read("data_raw/map17/原生自然環境保全地域_国指定自然環境保全地域.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝原生自然環境保全地域
        ZONE == 2 ~ "high", # 2＝特別地区
        ZONE == 3 ~ "high", # 3＝海中特別地区
        ZONE == 4 ~ "low" # 4＝普通地区
      )
    ),
  
  # 2: quasi-national parks
  protected_2 = sf::st_read("data_raw/map17/国定公園.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝特別保護地区
        ZONE == 20 ~ "high", # 20＝特別地域
        ZONE == 21 ~ "medium", # 21＝第1種特別地域
        ZONE == 22 ~ "medium", # 22＝第2種特別地域
        ZONE == 23 ~ "medium", # 23＝第3種特別地域
        ZONE == 3 ~ "low", # 3＝普通地区
        ZONE == 5 ~ "marine" #5＝海域公園地区
      )
    ),
  
  # 3: national wildlife protection areas
  protected_3 = sf::st_read("data_raw/map17/国指定鳥獣保護区.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "low", # 1＝鳥獣保護区（特別保護地区以外
        ZONE == 2 ~ "medium", # 2＝特別保護地区
      )
    ),
  
  # 4: national parks
  protected_4 = sf::st_read("data_raw/map17/国立公園.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝特別保護地区
        ZONE == 20 ~ "high", # 20＝特別地域
        ZONE == 21 ~ "medium", # 21＝第1種特別地域
        ZONE == 22 ~ "medium", # 22＝第2種特別地域
        ZONE == 23 ~ "medium", # 23＝第3種特別地域
        ZONE == 3 ~ "low", # 3＝普通地区
        ZONE == 5 ~ "marine" #5＝海域公園地区
      )
    ) %>%
    # Remove protected area in inland sea (marine)
    filter(NAME != "瀬戸内海"),
  
  # 5: prefectural wildlife protection areas
  protected_5 = sf::st_read("data_raw/map17/都道府県指定鳥獣保護区.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "low", # 1＝鳥獣保護区（特別保護地区以外
        ZONE == 2 ~ "medium", # 2＝特別保護地区
        ZONE == 3 ~ "low" # not specified, but assume no other special protection
      )
    ),
  
  # 6: prefectural natural parks
  protected_6 = sf::st_read("data_raw/map17/都道府県立自然公園.shp") %>%
    mutate(
      status = case_when(
        ZONE == 1 ~ "high", # 1＝特別保護地区
        ZONE == 20 ~ "high", # 20＝特別地域
        ZONE == 3 ~ "low" # 3＝普通地区
      )
    ),
  
  # 7: prefectural protection areas
  protected_7 = sf::st_read("data_raw/map17/都道府県自然環境保全地域.shp") %>%
    mutate(
      status = case_when(
        ZONE == 0 ~ "high", # 0＝原生自然環境保全地域
        ZONE == 2 ~ "high", # 2＝特別地区
        ZONE == 4 ~ "low" # 4＝普通地区
      )
    ),
  
  # Combine protected areas into single dataframe
  protected_areas = combine_protected_areas(
    protected_1,
    protected_2,
    protected_3,
    protected_4,
    protected_5,
    protected_6,
    protected_7
  ),
  
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
