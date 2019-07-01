# Define analysis plan
test_swap_plan <- drake_plan (
  
  # Load and process raw data ----
  
  # Pteridophyte Phylogeny Group I (PPGI) taxonomy
  # - original version
  
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
  
  # Make community matrix (presence/absence of each species in
  # 1km2 grid cells), trim to only species in tree.
  comm_ferns = make_comm_matrix(occ_data_ferns) %>% 
    match_comm_and_tree(japan_fern_tree, "comm"),
  
  # Calculate phylogenetic diversity.
  # Species are rows.
  comm_ferns_df = column_to_rownames(comm_ferns, "species"),
  
  # MPD 
  # Transpose community matrix, keeping as dataframe.
  # Test on first community (now communities are rows) only
  mpd_ferns = 
    ses.mpd(
      samp = t(comm_ferns_df)[1:2,], 
      dis = cophenetic(japan_fern_tree),
      null.model = "independentswap",
      iterations = 100,
      runs = 100)
  
)
