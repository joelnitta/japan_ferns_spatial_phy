# Define analysis plan
plan <- drake_plan (
  
  # Load data ----
  # - Electronic Supp. Mat. 1 contains data on reproductive
  # mode and GenBank accession number in one row per species.
  repro_data_raw = read_excel(file_in("data/ESM1.xlsx")),
  
  # - Electronic Supp. Mat. 2 contains occurrence data
  # with multiple rows per species. Occurrences are presences
  # in a set of 1km2 grid cells across Japan, not actual
  # occurrence points of specimens.
  occ_data_raw = read_excel(
    file_in("data/ESM2.xlsx"),
    col_types = c("numeric", "text", "numeric", 
                  "numeric", "numeric", "text", 
                  "text", "text")
  ),
  
  # Phylogenetic tree of all non-hybrid pteridophyte
  # taxa built from rbcL gene.
  japan_pterido_tree = read.nexus("data/PD170708Bayes2.nxs"),
  
  # Basic world map
  world_map = ggplot2::map_data("world"),
  
  # Process data ----
  repro_data = process_repro_data(repro_data_raw),

  occ_data = clean_names(occ_data_raw),
  
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
  
  # Run analysis of variance on CPS by reproductive mode.
  cps_by_repro_model_summary = aov(
    n_grids ~ reproductive_mode, 
    data = cps_by_repro) %>% tidy,
  
  # Make tibble with latitudinal breadth and
  # reproductive mode for all non-hybrids.
  lat_by_repro = count_lat_by_repro(occ_data, repro_data),
  
  # Calculate mean latitudinal breadth per species by reproductive mode.
  lat_by_repro_means = avg_lat_by_repro(lat_by_repro),
  
  # Run analysis of variance on 
  # latitudinal breadth by reproductive mode.
  lat_by_repro_model = aov(
    lat_breadth ~ reproductive_mode, 
    data = lat_by_repro),
  
  lat_by_repro_model_summary = tidy(lat_by_repro_model),
  
  # At least one mean is different, so run Tukey HSD test.
  lat_by_repro_tukey = TukeyHSD(lat_by_repro_model) %>% tidy,
  
  # Analyze community diversity ----
  
  # Make richness matrix (number of species per
  # 1km2 grid cell).
  richness = make_richness_matrix(occ_data),
  
  # Plots ----
  richness_map = make_richness_plot(
    richness, 
    occ_data, 
    world_map
  )
  
)