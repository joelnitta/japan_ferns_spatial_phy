# Data processing ----

#' Clean up reproductive mode data for pteridophytes of Japan
#'
#' @param data Tibble. Raw data read-in from Electronic Supp. Mat. 1
#'
#' @return Tibble
process_repro_data <- function (data) {
  
  data %>%
  clean_names %>%
    mutate(
      reproductive_mode = case_when(
        reproductive_mode == 0 ~ "unknown",
        reproductive_mode == 1 ~ "sexual", 
        reproductive_mode == 2 ~ "apomictic",
        reproductive_mode == 3 ~ "sex_apo",
        TRUE ~ "unknown"
      ) %>% as.factor,
      sexual_diploid = case_when(
        sexual_diploid == 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      sexual_polyploid = case_when(
        sexual_polyploid == 1 ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    rename(rbcl_genbank_no = rbc_l_gen_bank_accession_no) %>%
    mutate(rbcl_genbank_no = str_remove_all(rbcl_genbank_no, "\\*"))
  
}

# Basic stats ----

#' Count species per grid cell excluding hybrids
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode mata, with
#' one row per taxon, excluding hybrids.
#' @return tibble
count_species_per_cell <- function (occ_data, repro_data) {
  occ_data %>%
  filter(taxon_id %in% repro_data$taxon_id) %>%
  group_by(secondary_grid_code) %>%
  count(sort = TRUE)
}

#' Count grid cells per species
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @return tibble
count_cells_per_species <- function (occ_data) {
  occ_data %>%
  group_by(taxon_name) %>%
  count(sort = TRUE)
}

#' Count number of grid cells per species by reproductive mode
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode mata, with
#' one row per taxon, excluding hybrids.
#' @return tibble
count_cells_per_species_by_repro <- function(occ_data, repro_data) {
  occ_data %>%
  group_by(taxon_id) %>%
  summarize(
    n_grids = n()
  ) %>%
  inner_join(select(repro_data, taxon_id, reproductive_mode)) %>%
  ungroup() %>%
  filter(reproductive_mode != "unknown")
}

#' Get mean number of grid cells per species by reproductive mode
#'
#' @param cells_per_species_by_repro 
#'
#' @return Tibble
avg_cells_per_species_by_repro <- function (cells_per_species_by_repro) {
  cells_per_species_by_repro %>%
  group_by(reproductive_mode) %>%
  summarize(
    mean_grids = mean(n_grids),
    n = n(),
    sd = sd(n_grids)
  )
}

#' Determine latitudinal breadth for each species,
#' including reproductive mode
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param repro_data Reproductive mode mata, with
#' one row per taxon, excluding hybrids.
#'
#' @return Tibble. One species per row with latitudinal
#' breadth (max - min) and reproductive mode, excluding
#' hybrids and repro mode unknown.
count_lat_by_repro <- function(occ_data, repro_data) {
  occ_data %>%
  group_by(taxon_id) %>%
  summarize(
    lat_breadth = max(latitude) - min(latitude)
  ) %>%
  inner_join(select(repro_data, taxon_id, reproductive_mode)) %>%
  ungroup() %>%
  filter(reproductive_mode != "unknown")
}

#' Get mean latitudinal breadth per species by reproductive mode
#'
#' @param cells_per_species_by_repro 
#'
#' @return Tibble
avg_lat_by_repro <- function(lat_by_repro) {
  lat_by_repro %>%
  group_by(reproductive_mode) %>%
  summarize(
    lat_breadth = mean(lat_breadth),
    n = n(),
    sd = sd(lat_breadth)
  )
}

# Taxonomy ----

#' Modify Pteridophyte Phylogeny Group I taxonomy
#' to match the version used for Pteridophytes of Japan
#'
#' @param ppgi Original PPGI taxonomy
#'
#' @return tibble
modify_ppgi <- function (ppgi) {
  
  # Use normal "e" for Isoetes
  ppgi_mod <- 
    ppgi %>%
    select_if(is.character) %>%
    mutate_all(
      ~snakecase::to_any_case(
        ., 
        case = "none", 
        transliterations = c("Latin-ASCII"), 
        parsing_option = 0)
    )
  
  # Add genera missing in PPGI that are included in Japan pteridophyte checklist
  # Use the sister (or encompassing) genus for each, so other higher-order
  # taxonomy will be correct
  anisocampium_dat <-
    ppgi_mod %>% filter(genus == "Athyrium") %>%
    mutate(genus = "Anisocampium")
  
  humata_dat <-
    ppgi_mod %>% filter(genus == "Davallia") %>%
    mutate(genus = "Humata")
  
  drynaria_dat <-
    ppgi_mod %>% filter(genus == "Aglaomorpha") %>%
    mutate(genus = "Drynaria")
  
  bind_rows(ppgi_mod, anisocampium_dat, humata_dat, drynaria_dat) %>%
    select(genus, family, order, class)
  
}

#' Add taxonomy data to occurrence data
#'
#' @param occ_data Occurence data, including
#' at least one column called "taxon_name" where
#' the first word separated by spaces is the genus name.
#' @param taxonomy_data Taxonomy data at the
#' genus level and higher
#'
#' @return tibble
add_taxonomy <- function(occ_data, taxonomy_data) {
  occ_data %>%
    mutate(genus = str_split(taxon_name, " ") %>% map_chr(1)) %>%
    filter(!is.na(genus)) %>%
    left_join(taxonomy_data)
}

#' Make a tibble for mapping taxon IDs to taxon names
#'
#' @param occ_data_pteridos All occurrence data of pteridophytes in Japan
#'
#' @return tibble
#' 
make_taxon_id_map <- function(occ_data_pteridos) {
  occ_data_pteridos %>% select(taxon_id, taxon = taxon_name) %>% unique %>%
  mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
  mutate(taxon = str_remove_all(taxon, "\\.")) %>%
  mutate(taxon = str_remove_all(taxon, "_var")) %>%
  mutate(taxon = str_remove_all(taxon, "_subsp")) %>%
  mutate(taxon = str_remove_all(taxon, "_x"))
}

#' Resolve synonyms in matched names
#' 
#' Synonyms in matched names will be matched to the taxonomic standard.
#' If any are synonyms these will be replaced with the accepted name.
#'
#' @param matched_names Dataframe of names matched to a taxonomic standard
#' (output of match_with_gnr).
#' @param taxonomic_standard Dataframe of standard names in Darwin Core format
#'
#' @return Tibble
#' 
resolve_synonyms <- function(matched_names, taxonomic_standard) {
  
  # Check format of names standard input
  checkr::check_data(taxonomic_standard, values = list(
    taxonID = 1L,
    acceptedNameUsageID = c(1L, NA),
    taxonomicStatus = "a",
    scientificName = "a",
    genus = "a",
    specificEpithet = c("a", NA),
    infraspecificEpithet = c("a", NA)),
    key = "taxonID"
  )
  
  checkr::check_data(matched_names, values = list(
    query = "a",
    matched_name = c("a", NA)
  )
  )
  
  merged <-
    left_join(
      matched_names,
      taxonomic_standard,
      by = c(matched_name = "scientificName"))
  
  not_resolved <- merged %>%
    mutate(taxonomicStatus = replace_na(taxonomicStatus, "no match")) %>%
    filter(taxonomicStatus != "accepted name") %>%
    filter(taxonomicStatus != "synonym") %>%
    select(
      query,
      scientificName = matched_name,
      genus,
      specificEpithet,
      infraspecificEpithet,
      fail_reason
    )
  
  accepted_names <- merged %>%
    filter(taxonomicStatus == "accepted name") %>%
    select(
      query,
      scientificName = matched_name,
      genus,
      specificEpithet,
      infraspecificEpithet
    )
  
  synonyms <- merged %>%
    filter(taxonomicStatus == "synonym") %>%
    select(query, matched_name, acceptedNameUsageID) %>%
    left_join(taxonomic_standard, by = c(acceptedNameUsageID = "taxonID")) %>%
    select(
      query,
      scientificName,
      genus,
      specificEpithet,
      infraspecificEpithet
    )
  
  bind_rows(accepted_names, synonyms, not_resolved)
}


# Traits ----

#' Get the mean value of a numeric trait from data formatted for lucid
#'
#' @param x Vector of values with numeric trait data formatted for lucid
#'
#' @return Means of each trait value
#' 
get_lucid_mean <- function (x) {
  if(isTRUE(is.na(x))) return (NA)
  assertthat::assert_that(is.character(x))
  num_colons <- stringr::str_count(x, ":")
  assertthat::assert_that(num_colons == 3)
  vals <- stringr::str_split(x, ":") %>%
    map(as.numeric) %>%
    unlist
  mean(c(vals[[2]], vals[[3]]))
}


#' Transform traits
#'
#' @param traits dataframe; trait data with one column per trait.
#' @param log_trans logical; should log-transform be applied?
#' @param scale_traits logical; should traits be scaled?
#' @param small_number Arbitrarily small number to use in place
#' of 0 before log-transform
#' @param trans_select Character vector of trait names to log
#' transform.
#' @param scale_select Character vector of trait names to rescale.
#'
#' @return dataframe
#' 
transform_traits <- function (traits, 
                              log_trans = TRUE, 
                              scale_traits = TRUE, 
                              small_number = 0.1, 
                              trans_select = c("dissection", "stipe", "length", "width", 
                                               "rhizome", "pinna"), 
                              scale_select = c("sla", "dissection", "stipe", "length", 
                                               "width", "rhizome", "pinna")
) {
  
  # Log-transform
  if (log_trans == TRUE) {
    traits <-
      traits %>%
      assertr::verify(trans_select %in% colnames(traits)) %>%
      assertr::assert(is.numeric, trans_select) %>%
      # Replace zeros with arbitrarily small number
      dplyr::mutate_at(trans_select, ~ifelse(. == 0, small_number, .)) %>%
      dplyr::mutate_at(trans_select, log)
  }
  
  # Rescale by dividing original value by range of 
  # that value (max - min) across the dataset
  if (scale_traits == TRUE) {
    traits <-
      traits %>% 
      assertr::verify(scale_select %in% colnames(traits)) %>%
      assertr::assert(is.numeric, scale_select) %>%
      dplyr::mutate_at(
        scale_select, ~ . / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
      )
  }
  
  traits
  
}

#' Make a traits distance matrix
#' 
#' Use raw fern and lycophyte trait data
#' formatted for lucid
#'
#' @param path_to_lucid_traits Path to raw trait data
#'
#' @return Distance matrix
#' 
make_traits_dist_matrix <- function(path_to_lucid_traits, taxon_id_map) {
  
  # Read in raw trait data for pteridophytes of Japan.
  # These were originally formatted for lucid dichotomous key software.
  # So they are mostly quantitative traits that have been converted to binary format,
  # or numeric traits. There are a lot of traits. One row per taxon.
  traits <- read_excel("data_raw/JpFernLUCID_forJoel.xlsx", skip = 1) %>%
    clean_names() %>%
    select(-x1, -x222) %>%
    rename(taxon = x2) %>%
    mutate(taxon = str_replace_all(taxon, ":", "_"))
  
  # Check for NA values: nope
  traits %>% map(is.na) %>% unlist %>% any
  
  # Separate out into numeric and binary traits
  # (numeric container "number" in name, assume binary otherwise)
  traits_numeric <- select(traits, taxon, contains("number"))
  traits_binary <-  select(traits, -contains("number"))
  
  ### Cleanup binary traits ###
  
  # Check what are the unique values in the "binary" traits:
  # Almost all 0 and 1s, a few other values
  traits_binary %>%
    select(-taxon) %>%
    map(~unique(.) %>% sort) %>%
    unlist %>%
    table
  
  # These are what the values mean according to the lucid manual
  # "common and misinterpreted" means that the user may incorrectly think the trait is
  # absent when it's actually frequently present.
  # 
  # 0=absent
  # 1=common
  # 2=rare
  # 3=uncertain
  # 4=common and misinterpreted
  # 5=rare and misinterpreted
  # 6=not scoped
  
  # Reformat so all traits are either present (1), absent (0), or NA
  traits_binary <-
    traits_binary %>%
    mutate_at(vars(-taxon), as.numeric) %>%
    mutate_if(is.numeric, ~case_when(
      . == 2 ~ 1,
      . == 3 ~ NaN,
      . == 4 ~ 1,
      . == 5 ~ 1,
      TRUE ~ .
    )) 
  
  traits_binary %>%
    select(-taxon) %>%
    map(~unique(.) %>% sort) %>%
    unlist %>%
    table
  
  traits_binary %>% map(is.na) %>% unlist %>% sum
  
  # Reformat presence/absence traits. These have one column each for "presence" (0 or 1),
  # "absence" (also 0 or 1), and sometimes another related state ("caducous" etc).
  # Combine these into a single "present" column.
  
  # Check names of presence/absence columns
  traits_binary %>% select(contains("abs"), contains("pres")) %>%
    colnames %>% sort
  
  traits_binary <- 
    traits_binary %>%
    # - pseudo_veinlet
    mutate(
      leaf_lamina_pseudo_veinlet_present = case_when(
        leaf_lamina_pseudo_veinlet_absent == 1 ~ 0,
        TRUE ~ leaf_lamina_pseudo_veinlet_present
      )) %>%
    select(-leaf_lamina_pseudo_veinlet_absent) %>%
    # - rhachis_adaxial_side_grooved
    mutate(
      leaf_lamina_rhachis_adaxial_side_grooved_present = case_when(
        leaf_lamina_rhachis_adaxial_side_grooved_absent == 1 ~ 0,
        leaf_lamina_rhachis_adaxial_side_grooved_present_continuous_to_costa_groove == 1 ~ 1,
        leaf_lamina_rhachis_adaxial_side_grooved_present_not_continuous_to_costa_groove == 1 ~ 1,
        TRUE ~ 0
      )) %>%
    select(-leaf_lamina_rhachis_adaxial_side_grooved_absent,
           -leaf_lamina_rhachis_adaxial_side_grooved_present_continuous_to_costa_groove,
           -leaf_lamina_rhachis_adaxial_side_grooved_present_not_continuous_to_costa_groove) %>%
    # - terminal_pinna
    mutate(
      leaf_lamina_terminal_pinna_present = case_when(
        leaf_lamina_terminal_pinna_absent == 1 ~ 0,
        leaf_lamina_terminal_pinna_absent == 0 ~ 1
      )
    ) %>%
    select(-leaf_lamina_terminal_pinna_absent) %>%
    # - indusium
    mutate(
      leaf_sorus_indusium_present = case_when(
        leaf_sorus_indusium_presence_absence_absent == 1 ~ 0,
        leaf_sorus_indusium_presence_absence_present_caducous == 1 ~ 1,
        TRUE ~ leaf_sorus_indusium_presence_absence_present
      )) %>%
    select(-contains("leaf_sorus_indusium_presence_absence")) %>%
    rename(leaf_sorus_false_indusium_present = leaf_sorus_false_indusium)
  
  # Check names of presence/absence columns
  traits_binary %>% select(contains("abs"), contains("pres")) %>%
    colnames %>% sort
  
  #### Clean up numeric traits ###
  # Replace missing (0 or 3) with NA,
  # take the mean of the range of normal values otherwise
  traits_numeric <-
    traits_numeric %>%
    # All numeric traits are preceded with '1:', but this doesn't mean anything
    mutate_at(vars(-taxon), ~str_remove(., "^1\\:")) %>%
    mutate_at(vars(-taxon), ~na_if(., "0")) %>%
    mutate_at(vars(-taxon), ~na_if(., "3")) %>%
    mutate_at(vars(-taxon), ~map_dbl(., get_lucid_mean))
  
  # Split numeric traits into those measured on sterile vs fertile plants
  # (i.e., with or without spores)
  
  traits_numeric_sterile <- select(traits_numeric, taxon, contains("sterile")) %>%
    rename_all(~str_remove(., "sterile_")) %>%
    gather(trait, value, -taxon) %>%
    mutate(type = "sterile")
  
  traits_numeric_fertile <- select(traits_numeric, taxon, contains("fertile")) %>%
    rename_all(~str_remove(., "fertile_")) %>%
    gather(trait, value, -taxon) %>%
    mutate(type = "fertile")
  
  # Combine these, taking the maximum value regardless of sterile or fertile
  traits_numeric_combined <-
    bind_rows(traits_numeric_sterile, traits_numeric_fertile) %>%
    group_by(taxon, trait) %>%
    summarize(
      value = max(value, na.rm = TRUE)
    ) %>%
    ungroup %>%
    mutate(value = na_if(value, -Inf)) %>%
    spread(trait, value) %>%
    rename_all(~str_remove(., "_number"))
  
  # Log-transform and scale numeric traits
  num_trait_names <- select(traits_numeric_combined, -taxon) %>% colnames()
  
  traits_numeric_combined_trans <- transform_traits(
    traits_numeric_combined,
    trans_select = num_trait_names,
    scale_select = num_trait_names
  )
  
  # Subset to only completely sampled species
  traits_numeric_combined_trans_complete <- 
    filter(traits_numeric_combined_trans, complete.cases(traits_numeric_combined_trans))
  
  ### Distance matrix
  
  # Combine all numeric and categorical traits
  traits_for_dist <- left_join(traits_numeric_combined_trans, traits_binary)
  
  # Convert species names to taxon id codes
  missing_taxon_id <-
    traits_for_dist %>% left_join(taxon_id_map) %>%
    select(taxon_id, everything()) %>%
    select(taxon_id, taxon) %>%
    filter(is.na(taxon_id))
  
  assertthat::validate_that(
    nrow(missing_taxon_id) == 0,
    msg = glue::glue("{nrow(missing_taxon_id)} taxa missing taxa IDs and dropped")
  )
  
  traits_for_dist <-
    traits_for_dist %>% inner_join(taxon_id_map) %>%
    select(-taxon) 
  
  # Set up weighting.
  trait_categories <-
    traits_for_dist %>%
    select(-taxon_id) %>%
    colnames %>%
    tibble(trait = .) %>%
    mutate(value_type = case_when(
      trait %in% colnames(traits_binary) ~ "binary",
      trait %in% colnames(traits_numeric_combined_trans) ~ "numeric"
    )) %>%
    mutate(comp_trait = case_when(
      str_detect(trait, "leaf_sorus_shape") ~ "sorus_shape",
      str_detect(trait, "leaf_sorus_indusium_shape") ~ "indusium_shape",
      str_detect(trait, "leaf_sorus_indusium_margin") ~ "indusium_margin",
      str_detect(trait, "leaf_sorus_false_indusium_present") ~ "false_indusium_present",
      str_detect(trait, "leaf_sorus_indusium_present") ~ "indusium_present",
      str_detect(trait, "leaf_lamina_shape_sterile_frond") ~ "shape_sterile_frond",
      str_detect(trait, "leaf_lamina_texture") ~ "leaf_lamina_texture",
      str_detect(trait, "leaf_lamina_color") ~ "leaf_lamina_color",
      str_detect(trait, "leaf_lamina_vennation") ~ "leaf_lamina_vennation",
      str_detect(trait, "leaf_lamina_pseudo_veinlet_present") ~ "pseudo_veinlet_present",
      str_detect(trait, "leaf_lamina_terminal_pinna") ~ "terminal_pinna",
      str_detect(trait, "leaf_lamina_lateral_pinna_shape") ~ "lateral_pinna_shape",
      str_detect(trait, "leaf_lamina_lateral_pinna_stalk") ~ "lateral_pinna_stalk",
      str_detect(trait, "leaf_lamina_margin") ~ "leaf_lamina_margin",
      str_detect(trait, "leaf_lamina_rhachis_adaxial_side_grooved") ~ "rhachis_adaxial_side_grooved",
      str_detect(trait, "leaf_lamina_terminal_pinna") ~ "terminal_pinna",
      TRUE ~ trait
    )) %>%
    add_count(comp_trait) %>%
    mutate(weight_by_comp_trait = 1 / n) %>%
    add_count(value_type) %>%
    mutate(weight_by_type = 1 / n) %>%
    mutate(final_weight = weight_by_comp_trait * weight_by_type)
  
  trait_categories$weight_by_comp_trait %>% sum
  trait_categories$weight_by_type %>% sum
  trait_categories$final_weight %>% sum
  
  # Make sure traits for calculating the distance matrix are in correct
  # order for weighting
  traits_for_dist <- select(traits_for_dist, taxon_id, trait_categories$trait)
  
  # Convert to dataframe for gowdis
  traits_df <- traits_for_dist %>%
    column_to_rownames("taxon_id")
  
  # Run gowdis with trait weightings
  dist_mat <- FD::gowdis(traits_df, w = trait_categories$final_weight) 
  
}

# Community diversity ----

#' Calculate species richness for all grid cells
#' 
#' The long / lats of the occurrence data
#' are already cell centroids, so it is trivial
#' to compute species richness.
#'
#' @param occ_data Occurrence data, where each row is the 
#' occurrence of a species in a grid cell
#' @return tibble
make_richness_matrix <- function (occ_data) {
  occ_data %>%
    group_by(secondary_grid_code) %>%
    summarize(
      richness = n()
    ) %>%
    left_join(
      select(occ_data, latitude, longitude, secondary_grid_code) %>% unique
    )
}


#' Format tip labels in Japanese pteridophytes rbcL tree
#'
#' The tips in original phylogeny file are coded with two numbers 
#' separated by an underscore,
#' e.g., 601_14.
#' The second part is the taxon_id in occurrence and reproductive 
#' mode data. Keep only this as the tip label.
#'
#' @param phy Phylogeny of Japanese pteridophytes with tip labels
#' formatted as two numbers separated by an underscore
#'
#' @return List of class "phylo"
format_tip_labels <- function (phy) {
  
  phy$tip.label <- str_split(phy$tip.label, "_") %>% map_chr(2)
  
  phy

}

#' Make a community matrix
#'
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' 
#' @return tibble. One column for species then the rest
#' of the columns are presence/absence of that species in 
#' each site, where a "site" is a 1km2 grid cell. Names
#' of sites are grid-cell codes.
#' Species are stored as taxon_id values.
make_comm_matrix <- function (occ_data) {
  occ_data %>%
    select(species = taxon_id, site = secondary_grid_code) %>%
    mutate(
      abundance = 1,
      site = as.character(site)) %>%
    spread(site, abundance) %>%
    mutate_at(vars(-species), ~replace_na(., 0)) %>%
    mutate(species = as.character(species))
}

#' Calculate Faith's PD for a single community
#' 
#' Does not include the root when summing branch lengths.
#' 
#' At least two taxa in the community must be present in the tree; otherwise,
#' communities with 0 species will return 0,
#' communities with 1 species will return NA.
#' 
#' @param single_comm Single community, formatted as tibble with
#' two columns: 'species' (character) and 'abundance' (numeric)
#' @param phy Phylogenetic tree
#' @param shuffle_tips Logical; should tips of the phylogeny be
#' randomized? Used for generating null distributions.
#' 
#' @return a number: phylogenetic diversity (sum of branch lengths)
#' for that community
calc_pd <- function(single_comm, phy, shuffle_tips = FALSE) {
  
  # Filter to only species present in the focal community
  single_comm <- filter(single_comm, abundance > 0)
  
  # Return 0 if there are zero species present
  if(nrow(single_comm) == 0) return (0)
  
  # The phylogenetic distance for a single species without using the root is
  # undefined.
  if(nrow(single_comm) == 1) return (NA)
  
  # Optionally shuffle tips when generating null distributions
  if(isTRUE(shuffle_tips)) phy <- picante::tipShuffle(phy)
  
  # Prune tree to only species present in the community
  phy <- ape::keep.tip(phy, single_comm$species)
  
  # Get sum of branches remaining
  sum(phy$edge.length)
}

#' Generate random values of phylogenetic diversity
#' for a single community.
#' 
#' Used for generating null distributions. Randomization
#' done by shuffling the tips of the tree.
#' 
#' @param n_reps Number of times to repeat the randomization.
#' @param single_comm Single community, formatted as tibble with
#' two columns: 'species' (character) and 'abundance' (numeric)
#' @param phy Phylogenetic tree
#' 
#' @return numeric vector: randomized phylogenetic diversity values 
#' for that community
run_pd_rnd <- function(n_reps, single_comm, phy) {
  map_dbl(1:n_reps, ~calc_pd(single_comm = single_comm, phy = phy, shuffle_tips = TRUE))
}

#' Get the rank of a value amongst other numbers
#'
#' @param value The value of interest
#' @param other_nums The other guys
#'
#' @return Number: the rank of our value of
#' interest compared to the other numbers
#'
#' @examples
#' # Should be near 20
#' get_rank(20, runif(100) * 100)
#' get_rank(NA, runif(100) * 100)
get_rank <- function(value, other_nums) {
  if(is.na(value) | is.null(value)) return(NA)
  if(any(is.na(other_nums)) | is.null(other_nums)) return(NA)
  assert_that(assertthat::is.number(value))
  assert_that(is.numeric(other_nums))
  combined <- c(value, other_nums) %>% set_names(c("value", other_nums))
  which(names(sort(combined)) == "value")
}

#' Calculate the standard effect size (SES) of Faith's PD across multiple communities.
#' 
#' Does not include the root when summing branch lengths. At least two taxa in each
#' community must be present in the tree. Null communities are simulated by randomly shuffling 
#' taxon names at the tips of the tree.
#' 
#' @param comm community matrix. One column must be named
#' 'species', and the rest should correspond to presence or absence of species
#' in communities (sites).
#' @param phy phylogenetic tree
#' @param n_reps number of null communities to simulate for each real community
#' 
#' @return a data frame -- observed Faith's PD (pd_obs), 
#' mean value of PD across all simulated communities (pd_rand_mean),
#' standard deviation of PD across all simulated communities (pd_rand_sd),
#' rank of the observed PD relative to simulated 
#' values (pd_obs_rank), standard effect size of PD (ses_pd), and 
#' probability of the observed value (pd_obs_p).
ses_pd <- function (comm, phy, n_reps) {
  
  assert_that(isTRUE(all.equal(
  sort(comm$species), sort(phy$tip.label) )),
  msg = "Species don't match exactly between 'comm' and 'phy'"
  )
  
  # Nest by community, then calculate PD for each
  comm %>%
    gather(site, abundance, -species) %>%
    nest(-site) %>%
    mutate(
      pd_obs = map_dbl(data, ~ calc_pd(., phy = phy)),
      pd_rnd = map(data, ~ run_pd_rnd(single_comm = ., phy = phy, n_reps = n_reps)),
      pd_rand_mean = map_dbl(pd_rnd, ~ mean(., na.rm = TRUE)),
      pd_rand_sd = map_dbl(pd_rnd, ~ sd(., na.rm = TRUE)),
      ses_pd = (pd_obs - pd_rand_mean) / pd_rand_sd,
      pd_obs_rank = map2_dbl(pd_obs, pd_rnd, ~ get_rank(.x, .y)),
      pd_obs_p = pd_obs_rank/(n_reps + 1)
      ) %>%
    select(-data, -pd_rnd)
  
}

#' Merge diversity metrics
#'
#' @param all_pd Phylogenetic diversity metrics including
#' observed PD, SES of PD, etc.
#' @param richness Species richness of communities (1km2 grid
#' cells)
#' @param all_cells Longitude, latitude, and elevation of all 
#' grid cells (including those not in `all_pd` or `richness`)
#'
#' @return Tibble including rows for all grid cells in 1km2 grid cells
#' of Japan dataset. NA (pd values) or 0 (richness) will be entered
#' for grid cells that weren't in `all_pd` or `richness`.
merge_metrics <- function (all_pd, richness, all_cells) {
  all_pd %>%
    rename(secondary_grid_code = site) %>%
    left_join(richness) %>%
    select(-latitude, -longitude) %>%
    right_join(all_cells) %>%
    mutate(richness = replace_na(richness, 0))
}

#' Match community data and tree
#' 
#' Order of species in comm will be rearranged to match the
#' phylogeny.
#'
#' @param comm Community data frame, with one column for sites and
#' the rest for species.
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_comm_and_tree <- function (comm, phy, return = c("comm", "tree")) {
  
  assert_that("species" %in% colnames(comm))
  
  # Keep only species in phylogeny
  comm <- comm %>%
    filter(species %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, comm$species))
  
  # Get comm in same order as tips
  comm <- left_join(
    tibble(species = phy$tip.label),
    comm
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))
  
  # Return comm or tree
  assert_that(return %in% c("tree", "comm"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (comm)
  }
  
}

# Ecostructure ----

#' Convert community tibble into matrix for ecostructure
#'
#' @param comm_pteridos Tibble of community data, with one
#' column for species names, the rest sites, and rows
#' as species.
#' @param all_cells Tibble with one column corresponding to
#' site names, one for longitude, and one for latitude
#'
#' @return Matrix with rows as sites and columns as species.
#' Rownames are longitude, latitude separated by underscore.
#' 
make_ecos_matrix <- function (comm_pteridos, all_cells) {
  
  # Get vector of site names
  pterido_sites <-
    comm_pteridos %>% 
    gather(secondary_grid_code, abundance, -species) %>%
    pull(secondary_grid_code) %>%
    unique
  
  # Make tibble of all longitudes and latitudes by
  # site name (here, secondary_grid_code) for renaming
  # rows in matrix
  long_lats <-
    all_cells %>%
    transmute(secondary_grid_code, long_lat = paste(longitude, latitude, sep = "_")) %>%
    unique %>%
    # Only use those actually in the community data so we
    # don't end up with a bunch of NA values after joining.
    filter(secondary_grid_code %in% pterido_sites)
  
  comm_pteridos %>% 
    gather(secondary_grid_code, abundance, -species) %>%
    spread(species, abundance) %>%
    # Remove grid cell flagged as missing until get fixed data
    filter(secondary_grid_code != "513613") %>%
    # Check for grid cells in comm data but missing from all_cells data
    verify(secondary_grid_code %in% long_lats$secondary_grid_code) %>%
    # Check that grid cells in comm data are unique
    assert(is_uniq, secondary_grid_code) %>%
    left_join(long_lats) %>%
    select(-secondary_grid_code) %>%
    column_to_rownames("long_lat") %>%
    as.matrix
  
}

# Slightly tweaked version of ecos_plot_pie that takes sf object as background map
# instead of reading in shape file.
ecos_plot_pie2 <- function (
  omega = NULL, coords = NULL, bgmap_path = NULL, adjust = FALSE, 
  thresh = 0.7, long_lim = c(-180, 180), lat_lim = c(-60, 90), 
  coastline_lwd = 10, intensity = 1, radius = 0.5, color = c("dodgerblue2", 
                                                             "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", 
                                                             "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", 
                                                             "gray70", "khaki2", "maroon", "orchid1", "deeppink1", 
                                                             "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", 
                                                             "yellow3", "darkorange4", "brown", "red", "cornflowerblue", 
                                                             "cyan", "brown4", "burlywood", "darkgoldenrod1", "azure4", 
                                                             "green", "deepskyblue", "yellow", "azure1"), pie_control = list(), 
  image_width = 1000, image_height = 800, path = "geostructure_plot.tiff") 
{
  if (is.null(coords)) {
    if (is.null(rownames(omega))) {
      stop("coords not provided, omega rownames do not have latitude longitude\n           information either")
    }
    latlong_chars <- rownames(omega)
    coords <- cbind.data.frame(as.numeric(sapply(latlong_chars, 
                                                 function(x) strsplit(x, "_")[[1]][1])), as.numeric(sapply(latlong_chars, 
                                                                                                           function(x) strsplit(x, "_")[[1]][2])))
    colnames(coords) <- c("lat", "long")
  }
  else {
    if (dim(coords)[1] != dim(omega)[1]) {
      stop("coords provided, but the number of rows in coords data does not\n           match the number of rows in omega matrix")
    }
  }
  pie_control_default <- list(edges = 200, clockwise = TRUE, 
                              init.angle = 90, density = NULL, angle = 45, border = NULL, 
                              lty = NULL, label.dist = 1.1)
  pie_control <- modifyList(pie_control_default, pie_control)
  if (is.null(bgmap_path)) {
    message("reading background map shapefile from inst/extdata/ne_110m_coastline \n            folder")
    GlobalCoast <- sf::st_read(system.file("extdata", "ne_110m_land", 
                                           "ne_110m_land.shp", package = "ecostructure"), quiet = T)
  }
  else {
    GlobalCoast <- bgmap_path
  }
  glob <- c(xmin = long_lim[1], xmax = long_lim[2], ymin = lat_lim[1], 
            ymax = lat_lim[2])
  glob <- sf::st_bbox(glob)
  glob <- structure(glob, crs = sf::st_crs(GlobalCoast))
  GlobalCoast <- suppressWarnings(suppressMessages(sf::st_intersection(GlobalCoast, 
                                                                       sf::st_as_sfc(glob))))
  if (adjust) {
    idx <- which(omega[, 1] > thresh)
    omega <- omega[-idx, ]
    coords <- coords[-idx, ]
    omega <- omega[, -1]
    omega <- t(apply(omega, 1, function(x) return(x/sum(x))))
  }
  output_type <- strsplit(path, "[.]")[[1]][2]
  if (output_type == "tiff") {
    tiff(path, width = image_width, height = image_height)
  }
  else if (output_type == "png") {
    png(path, width = image_width, height = image_height)
  }
  else if (output_type == "pdf") {
    pdf(path, width = image_width, height = image_height)
  }
  else {
    stop("the output image may either be of  tiff, png or pdf extension")
  }
  plot(sf::st_geometry(GlobalCoast), axes = T, main = "", reset = F, 
       xaxs = "i", yaxs = "i", lwd = coastline_lwd)
  par(lwd = 0.01)
  invisible(lapply(1:dim(omega)[1], function(r) do.call(mapplots::add.pie, 
                                                        append(list(z = as.integer(100 * omega[r, ]), x = coords[r, 
                                                                                                                 1], y = coords[r, 2], labels = c("", "", ""), radius = radius, 
                                                                    col = sapply(color, scales::alpha, intensity)), pie_control))))
  invisible(dev.off())
}

# Plotting ----

#' Get the lower, upper, or absolute maximum value
#' of a variable in a dataframe
#' 
#' For setting plotting limits manually
#'
#' @param data Dataframe
#' @param var Name of variable (column) in dataframe
#' @param digits Number of digits desired in output
#' @param type Type of limit to calculate: "min" is lower,
#' "max" is upper, and "abs" is the absolute greatest value.
#'
#' @return Number
#'
#' @examples
#' get_limit(mtcars, disp, "max")
get_limit <- function (data, var, type = c("min", "max", "abs"), digits = 2) {
  
  var <- enquo(var)
  
  switch(type,
         
         max = data %>% 
           pull(!!var) %>% 
           max(na.rm = TRUE) %>% 
           multiply_by(10^digits) %>% ceiling %>% divide_by(10^digits),
         
         min = data %>% 
           pull(!!var) %>% 
           min(na.rm = TRUE) %>% 
           multiply_by(10^digits) %>% floor %>% divide_by(10^digits),
         
         abs = c(
           data %>% pull(!!var) %>% max(na.rm = TRUE),
           data %>% pull(!!var) %>% min(na.rm = TRUE)) %>% 
           abs %>% max %>%
           multiply_by(10^digits) %>% ceiling %>% divide_by(10^digits)
  )
  
}

#' Make a plot showing species richness on a map of Japan
#'
#' @param richness Species richness in 1km grid cells, 
#' with longitude and latitude
#' @param occ_data Occurrence data used to generate species richness matrix
#' @param world_map Background world map
#'
#' @return ggplot object

#' Make a plot showing selected alpha diversity metric on a map of Japan
#'
#' @param div_data Alpha diversity matrix; rows are communities
#' (1km2 grid cells), and columns are various alpha diversity metrics.
#' @param world_map Background world mapp
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param div_metric Selected diversity metric to plot.
#' Must one of the column names of div_dat.
#' @param metric_title Character: title to use for legend for
#' diversity metric.
#'
#' @return ggplot object
make_diversity_map <- function (div_data, world_map, occ_data, div_metric, metric_title) {
  
  div_metric <- sym(div_metric)
  
  ggplot(world_map, aes(x = longitude, y = latitude)) +
    geom_polygon(aes(group = group), fill = "light grey") +
    geom_tile(data = div_data,
              aes(fill = !!div_metric),
              color = "black") + 
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor, 
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor, 
               pull(occ_data, latitude) %>% max %>% ceiling)
    ) +
    # scale_fill_viridis_c(na.value="transparent") +
    jntools::blank_x_theme() +
    jntools::blank_y_theme() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    labs(
      fill = metric_title
    )
  
}

#' Make a plot showing selected SES of PD on map of Japan with
#' only significant communities colored.
#'
#' @param div_data Alpha diversity matrix; rows are communities
#' (1km2 grid cells), and columns are various alpha diversity metrics.
#' @param world_map Background world mapp
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#'
#' @return ggplot object
make_highlight_map <- function (div_data, world_map, occ_data, div_metric, sig_metric, metric_title) {
  
  div_metric <- sym(div_metric)
  sig_metric <- sym(sig_metric)
  
  # gghighlight only "knows about" data in the most recent layer
  # So make plot with world map on top of highlighted SES of PD,
  # then rearrange layers.
  plot <-
    ggplot(div_data, aes(x = longitude, y = latitude)) +
    geom_tile(aes(fill = !!div_metric), color = "black") +
    gghighlight( (!!sig_metric > 0.975 | !!sig_metric < 0.025) & !is.na(!!sig_metric) ) +
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor, 
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor, 
               pull(occ_data, latitude) %>% max %>% ceiling)
    ) +
    geom_polygon(data = world_map, aes(group = group), fill = "light grey")
  
  # See
  # https://stackoverflow.com/questions/20249653/insert-layer-underneath-existing-layers-in-ggplot2-object
  plot$layers <- plot$layers[c(1,3,2)]
  
  plot +
    scale_fill_scico(palette = "vik", na.value="transparent") +
    jntools::blank_x_theme() +
    jntools::blank_y_theme() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    labs(
      fill = metric_title
    )
  
}


#' Make scatter plot with linear model
#'
#' @param data Input data
#' @param resp_var Name of response variable in data
#' @param indep_var Name of independent variable in data
#' @param resp_var_print Label to print on y axis
#' @param indep_var_print Label to print on x axis
#' @param digits_p Number of digits to output for p value
#' @param digits_r2 Number of digits to output for R-squared
#'
#' @return ggplot object
#' 
#' @examples
#' make_scatter_with_lm(mtcars, "mpg", "disp")
make_scatter_with_lm <- function(data, indep_var, resp_var,
                                 indep_var_print = indep_var,
                                 resp_var_print = resp_var,
                                 digits_p = 3, digits_r2 = 3) {
  
  resp_var_sym <- sym(resp_var)
  indep_var_sym <- sym(indep_var)
  
  model <- lm(formula(glue::glue("{resp_var} ~ {indep_var}")), data = data)
  
  p <- model %>% glance %>% pull(p.value) %>% round(digits_p)
  r2 <- model %>% glance %>% pull(r.squared) %>% round(digits_r2)
  
  plot <- ggplot(data, aes(!!indep_var_sym, !!resp_var_sym)) +
    geom_point(alpha = 0.2) +
    jntools::standard_theme() +
    labs(
      y = resp_var_print,
      x = indep_var_print
    )
  
  if(p < 0.05) {
    plot <- plot +
      geom_smooth(
        method = "lm",
        size = 0.5,
        fill = "pink"
      ) +
      annotate("text", 
               x = Inf, y = Inf, 
               label = glue::glue("italic(R) ^ 2 == {r2}"),
               parse = TRUE,
               hjust = 1.2,
               vjust = 1.2)
  }
  
  plot
  
}
