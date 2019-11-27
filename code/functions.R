# Data processing ----

#' Read in a nexus file contained in a zipped archive
#'
#' @param zip_folder Path to zip file
#' @param nexus_file Name of nexus file within zip file
#'
#' @return List
#'
read_nexus_in_zip <- function (zip_folder, nexus_file) {

  temp_dir <- tempdir()

  unzip(zip_folder, exdir = temp_dir)

  ape::read.nexus(fs::path(temp_dir, nexus_file))

}

#' Read in Catalog of Life plants data
#'
#' @param path_to_col Path to Catalog of Life plants data
#'
#' @return Tibble
#'
read_col_plants <- function (path_to_col) {

  read_tsv(path_to_col,
    col_types = rep("c", 31) %>% paste(collapse = "")
  ) %>%
    mutate(taxonID = as.numeric(taxonID), acceptedNameUsageID = as.numeric(acceptedNameUsageID))

}


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

#' Get a list of URLs to download data for a kokudosuuchi dataset
#'
#' @param dataset Name of the dataset
#'
#' @return Character vector; list of URLs for that dataset
#'
#' @examples
#' get_ksj_urls("標高・傾斜度3次メッシュ")
get_ksj_urls <- function (dataset) {
  # Get all KSJ metadata
  ksj_metadata <- kokudosuuchi::getKSJSummary()
  # Get identifier for the dataset we're interested in
  dataset_id <- dplyr::filter(ksj_metadata, stringr::str_detect(title, dataset)) %>%
    dplyr::pull(identifier)
  # Get a tibble of URLs to zip files with this data
  urls <- kokudosuuchi::getKSJURL(dataset_id)
  urls$zipFileUrl
}

#' Unzip Ebihara and Nitta 2017 Ecol Mono data zip file and 
#' extract needed data files.
#' 
#' The dryad data zip file should be downloaded from 
#' https://datadryad.org/stash/dataset/doi:10.5061/dryad.df59g
#' (click on "Download dataset")
#'
#' @param dryad_zip_file Path to the data zip file downloaded from Dryad.
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @param ... Extra arguments; not used by this function, but
#' meant for tracking with drake.
#' @return Unzipped data files:
#' - rbcL_clean_sporos.fasta: rbcL sequences of sporophytes from Moorea
#'
unzip_ebihara_2019 <- function (dryad_zip_file, exdir, ...) {
  
  # Unzip only the needed files
  unzip(dryad_zip_file, "FernGreenListV1.01E.xls", exdir = exdir)
  unzip(dryad_zip_file, "ESM1.csv", exdir = exdir)
  unzip(dryad_zip_file, "ESM2.csv", exdir = exdir)
  unzip(dryad_zip_file, "japan_pterido_rbcl_cipres.zip", exdir = exdir)
  unzip(dryad_zip_file, "2_grid_cells_all.csv", exdir = exdir)
  unzip(dryad_zip_file, "ppgi_taxonomy.csv", exdir = exdir)
  
}

#' Tidy taxonomic data of pteridophytes of Japan
#'
#' Data is from Japan Green list
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
tidy_japan_names <- function (data) {
  data %>%
    select(taxon_id = ID20160331, scientific_name = `GreenList Name`,
           endemic = Endemism, conservation_status = RL2012) %>%
    mutate(taxon_id = as.character(taxon_id)) %>%
    select(taxon_id, scientific_name)
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

#' Update automatically resolved names with manually fixed names
#'
#' @param resolved_names_auto Dataframe of automatically resolved fern names;
#' output of taxastand::resolved_fern_names()
#' @param resolved_names_manual_fix Dataframe of manually fixed names to add for
#' names that failed to be resolved automatically
#'
#' @return Tibble
#'
update_resolved_names <- function (resolved_names_auto, resolved_names_manual_fix) {

  resolved_names_auto %>%
    # Drop excluded names
    filter_at(vars(contains("exclude")), all_vars(. == FALSE)) %>%
    select(-contains("exclude")) %>%
    # Drop anything not resoved to species
    filter(!is.na(species)) %>%
    # Add manually fixed names
    select(query, species) %>%
    bind_rows(
      select(resolved_names_manual_fix, query, species)
    )

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

#' Format traits for further analysis
#'
#' Use raw fern and lycophyte trait data
#' formatted for lucid. Output taxon names will be IDs.
#'
#' @param path_to_lucid_traits Path to raw trait data
#' @param taxon_id_map Tibble mapping taxon IDs to taxon names
#'
#' @return Tibble
#'
format_traits <- function(path_to_lucid_traits, taxon_id_map) {

  # Read in raw trait data for pteridophytes of Japan.
  # These were originally formatted for lucid dichotomous key software.
  # So they are mostly quantitative traits that have been converted to binary format,
  # or numeric traits. There are a lot of traits. One row per taxon.
  traits <- read_excel("data_raw/JpFernLUCID_forJoel.xlsx", skip = 1) %>%
    clean_names() %>%
    select(-x1, -x222) %>%
    rename(taxon = x2) %>%
    mutate(taxon = str_replace_all(taxon, ":", "_")) %>%
    # Check for NA values
    assert(not_na, everything())

  # Separate out into numeric and binary traits
  # (numeric container "number" in name, assume binary otherwise)
  traits_numeric <- select(traits, taxon, contains("number"))
  traits_binary <-  select(traits, -contains("number"))

  ### Cleanup binary traits ###

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

  # Reformat presence/absence traits. These have one column each for "presence" (0 or 1),
  # "absence" (also 0 or 1), and sometimes another related state ("caducous" etc).
  # Combine these into a single "present" column.

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

  # Drop any traits that have only one state
  drop_single_state <-
    map_df(traits_for_dist, n_distinct) %>%
    gather(trait, n_states) %>%
    filter(n_states == 1) %>%
    pull(trait)

  traits_for_dist[,!colnames(traits_for_dist) %in% drop_single_state]
}

#' Categorize binary traits into their categorical versions
#'
#' Detects traits by name and groups them into categories
#' (combined traits)
#'
#' @param traits Tibble of traits in binary format
#'
#' @return Tibble
#'
categorize_traits <- function(traits) {
  traits %>%
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
    ))
}

#' Summarize traits for supplemental info
#'
#' @param traits_for_dist Formatted trait data
#'
#' @return Tibble
#'
make_trait_summary <- function (traits_for_dist) {

  # Count trait states to categorize traits.
  # Those with > 3 states should be numeric.
  # (binary can be 1, 0, or missing).
  states <- map_df(traits_for_dist, n_distinct) %>%
    gather(trait, n_states)

  # Summarize traits
  traits_for_dist %>%
    select(-taxon_id) %>%
    colnames %>%
    tibble(trait = .) %>%
    left_join(states) %>%
    mutate(trait_type = case_when(
      n_states > 3 ~ "continuous",
      TRUE ~ "binary"
    )) %>%
    # Collapse binarized traits into their
    # categorical version by name.
    categorize_traits %>%
    add_count(comp_trait) %>%
    select(comp_trait, trait_type, n) %>%
    unique() %>%
    mutate(trait_type = case_when(
      trait_type == "continuous" ~ trait_type,
      n > 1 ~ "qualitative",
      TRUE ~ trait_type
    )) %>%
    mutate(n_states = case_when(
      trait_type == "qualitative" ~ n
    )) %>%
    select(trait = comp_trait, trait_type, n_states) %>%
    arrange(trait_type, trait)
}

#' Make a traits distance matrix
#'
#' @param traits_for_dist Formatted trait data
#'
#' @return Distance matrix
#'
make_trait_dist_matrix <- function (traits_for_dist) {

  # Count trait states to categorize traits.
  # Those with > 3 states should be numeric.
  # (binary can be 1, 0, or missing).
  states <- map_df(traits_for_dist, n_distinct) %>%
    gather(trait, n_states)

  # Set up weighting.
  trait_categories <-
    traits_for_dist %>%
    select(-taxon_id) %>%
    colnames %>%
    tibble(trait = .) %>%
    left_join(states) %>%
    mutate(trait_type = case_when(
      n_states > 3 ~ "continuous",
      TRUE ~ "binary"
    )) %>%
    # Collapse binarized traits into their
    # categorical version by name.
    categorize_traits %>%
    add_count(comp_trait) %>%
    mutate(weight_by_comp_trait = 1 / n) %>%
    add_count(trait_type) %>%
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
  FD::gowdis(traits_df, w = trait_categories$final_weight)

}

#' Make NMDS plots for pteridophytes based on traits
#'
#' Makes two plots: one at order level, and another just for
#' Polypodiales.
#'
#' @param nmds Output of vegan::metaMDS() on traits matrix
#' @param ppgi PPGI taxonomy
#' @param taxon_id_map Dataframe mapping taxon IDs to PPGI species names
#'
#' @return GGPlot object
#'
make_trait_nmds_plot <- function (nmds, ppgi, taxon_id_map) {

  nmds_points = nmds[["points"]] %>%
    as.data.frame %>%
    rownames_to_column("taxon_id") %>%
    as_tibble %>%
    left_join(taxon_id_map) %>%
    mutate(genus = str_split(taxon, "_") %>% map_chr(1)) %>%
    left_join(ppgi) %>%
    mutate(
      order = factor(order),
      family = factor(family))

  a <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = order, shape = order)) +
    geom_point() +
    scale_shape_manual(values = 1:n_distinct(nmds_points$order)) +
    labs(title = "Pteridophytes")

  nmds_points_polypod <- filter(nmds_points, order == "Polypodiales")

  b <- ggplot(nmds_points_polypod, aes(x = MDS1, y = MDS2, color = family, shape = family)) +
    geom_point() +
    scale_shape_manual(values = 1:n_distinct(nmds_points_polypod$family)) +
    labs(title = "Polypodiales")

  a + b
}

#' Make trait dendrogram
#'
#' @param trait_distance_matrix Distance matrix based on traits
#' @param taxon_id_map  Dataframe mapping taxon IDs to PPGI species names
#'
#' @return Nothing; externally, the plot will be written to
#' results/traits_dendrogram.pdf
#'
make_traits_dendrogram <- function(trait_distance_matrix, taxon_id_map) {

  # Cluster species by trait distances
  h <- hclust(trait_distance_matrix)

  # Relabel from taxon id to species names
  h[["labels"]] <- taxon_id_map$taxon[match(h[["labels"]], taxon_id_map$taxon_id)]

  # Output plot
  pdf(height = 50, width = 8, file = "results/traits_dendrogram.pdf")
  plot(ape::as.phylo(h), cex = 0.4)
  dev.off()
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
#' phylogeny. Species missing from either will be dropped.
#'
#' @param comm Community data frame, with one column for species and
#' the rest for sites (i.e., rows as species and columns as sites).
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
#' library(picante)
#' data(phylocom)
#' phylo <- phylocom$phylo
#' comm_small <- phylocom$sample[,sample(1:ncol(phylocom$sample), 10, replace = FALSE)]
#' comm_small <- comm_small %>% t %>% as.data.frame %>% rownames_to_column("species")
#' phylo_small <- ape::drop.tip(phylo, phylo$tip.label[sample(1:length(phylo$tip.label), 10, replace = FALSE)])
#' match_comm_and_tree(comm_small, phylo_small, "comm")
#' match_comm_and_tree(comm_small, phylo_small, "tree")
#'
match_comm_and_tree <- function (comm, phy, return = c("comm", "tree")) {

  assert_that("species" %in% colnames(comm))

  comm_original <- comm
  phy_original <- phy

  # Keep only species in phylogeny
  comm <- comm %>%
    filter(species %in% phy$tip.label)

  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, comm$species))

  # Get comm in same order as tips
  comm <- left_join(
    tibble(species = phy$tip.label),
    comm,
    by = "species"
  )

  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, phy$tip.label)))

  if(nrow(comm_original) != nrow(comm)) {
    print(glue::glue("Dropped {nrow(comm_original) - nrow(comm)} species not in phy from comm"))
  }

  if(ape::Ntip(phy_original) != ape::Ntip(phy)) {
    print(glue::glue("Dropped {Ntip(phy_original) - Ntip(phy)} species not in comm from phy"))
  }

  # Return comm or tree
  assert_that(return %in% c("tree", "comm"))

  if(return == "tree") {
    return (phy)
  } else {
    return (comm)
  }

}

#' Match community data and traits
#'
#' Order of species in comm will be rearranged to match the
#' phylogeny. Species missing from either will be dropped.
#'
#' @param comm Community data frame, with one column for species and
#' the rest for sites (i.e., rows as species and columns as sites).
#' @param traits Traits distances matrix
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the community, pruned so that only species occurring in both datasets
#' are included.
#'
#' @examples
#' library(picante)
#' library(tidyverse)
#' data(phylocom)
#' comm_small <- phylocom$sample[,sample(1:ncol(phylocom$sample), 10, replace = FALSE)] %>%
#'   t %>% as.data.frame %>% rownames_to_column("species")
#' traits_small <- phylocom$traits[sample(1:nrow(phylocom$traits), 10, replace = FALSE) ,]
#' dist_mat <- FD::gowdis(traits_small)
#' match_comm_and_traits(comm_small, dist_mat, "comm")
#' match_comm_and_traits(comm_small, dist_mat, "traits")
match_comm_and_traits <- function (comm, traits, return = c("comm", "traits")) {

  assert_that("species" %in% colnames(comm))

  comm_original <- comm
  traits_original <- traits

  # Keep only species in traits
  comm <- comm %>%
    filter(species %in% attributes(traits)[["Labels"]])

  # Trim traits matrix to only species in comm
  in_both <- base::intersect(attributes(traits)[["Labels"]], comm$species)
  traits <- usedist::dist_subset(traits, in_both)

  # Get comm in same order as traits
  comm <- left_join(
    tibble(species = attributes(traits)[["Labels"]]),
    comm,
    by = "species"
  )

  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, attributes(traits)[["Labels"]])))

  if(nrow(comm_original) != nrow(comm)) {
    print(glue::glue("Dropped {nrow(comm_original) - nrow(comm)} species not in traits from comm"))
  }

  traits_labels_original <- attributes(traits_original)[["Labels"]]
  traits_labels <- attributes(traits)[["Labels"]]

  if(length(traits_labels_original) != length(traits_labels)) {
    print(glue::glue("Dropped {length(traits_labels_original) - length(traits_labels)} species not in comm from traits"))
  }

  # Return comm or tree
  assert_that(return %in% c("traits", "comm"))

  if(return == "traits") {
    return (traits)
  } else {
    return (comm)
  }

}

#' Analyze Standard Effect Size (SES) of mean phylogenetic distance (MPD)
#'
#' Although at least >3 species must match between comm and tree, any non-matching
#' species will be dropped before running the analysis.
#'
#' @param comm Community data frame, with one column for species and
#' the rest for sites (i.e., rows as species and columns as sites).
#' @param tree Phylogeny (list of class "phylo")
#' @param species_col Name of column with species names
#' @param null.model Type of null model to use for picante::ses.mpd
#' @param iterations Number of iterations to use for "independentswap"
#' null model
#' @param runs Number of times to perform the randomization
#'
#' @return Dataframe. Results of ape::ses.mpd
#'
ses_phy_mpd <- function(comm, tree, species_col = "species",
                    null.model = "independentswap",
                    iterations = 10000,
                    runs = 999) {

  assertthat::assert_that(
    species_col %in% colnames(comm),
    msg = "value for 'species_col' not one of the columns of comm")

  # Drop species not matching between community and tree
  comm <- match_comm_and_tree(comm, tree, "comm")
  tree <- match_comm_and_tree(comm, tree, "tree")

  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, tree$tip.label)))

  # Convert community to dataframe with rows as sites and columns as species
  comm_df <- tibble::column_to_rownames(comm, species_col) %>% t()

  # Run ses mpd
  picante::ses.mpd(
    samp = comm_df,
    dis = cophenetic(tree),
    null.model = null.model,
    iterations = iterations,
    runs = runs)

}

#' Analyze Standard Effect Size (SES) of mean nearest taxon distance (MNTD)
#'
#' Although at least >3 species must match between comm and tree, any non-matching
#' species will be dropped before running the analysis.
#'
#' @param comm Community data frame, with one column for species and
#' the rest for sites (i.e., rows as species and columns as sites).
#' @param tree Phylogeny (list of class "phylo")
#' @param species_col Name of column with species names
#' @param null.model Type of null model to use for picante::ses.mpd
#' @param iterations Number of iterations to use for "independentswap"
#' null model
#' @param runs Number of times to perform the randomization
#'
#' @return Dataframe. Results of ape::ses.mpd
#'
ses_phy_mntd <- function(comm, tree, species_col = "species",
                        null.model = "independentswap",
                        iterations = 10000,
                        runs = 999) {

  assertthat::assert_that(
    species_col %in% colnames(comm),
    msg = "value for 'species_col' not one of the columns of comm")

  # Drop species not matching between community and tree
  comm <- match_comm_and_tree(comm, tree, "comm")
  tree <- match_comm_and_tree(comm, tree, "tree")

  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, tree$tip.label)))

  # Convert community to dataframe with rows as sites and columns as species
  comm_df <- tibble::column_to_rownames(comm, species_col) %>% t()

  # Run ses mpd
  picante::ses.mntd(
    samp = comm_df,
    dis = cophenetic(tree),
    null.model = null.model,
    iterations = iterations,
    runs = runs)

}

#' Analyze Standard Effect Size (SES) of functional mean phylogenetic distance (MPD)
#'
#' Although at least >3 species must match between comm and tree, any non-matching
#' species will be dropped before running the analysis.
#'
#' @param comm Community data frame, with one column for species and
#' the rest for sites (i.e., rows as species and columns as sites).
#' @param traits Distance matrix of traits
#' @param species_col Name of column with species names
#' @param null.model Type of null model to use for picante::ses.mpd
#' @param iterations Number of iterations to use for "independentswap"
#' null model
#' @param runs Number of times to perform the randomization
#'
#' @return Dataframe. Results of ape::ses.mpd
#'
ses_func_mpd <- function(comm, traits, species_col = "species",
                         null.model = "independentswap",
                         iterations = 10000,
                         runs = 999) {

  assertthat::assert_that(
    species_col %in% colnames(comm),
    msg = "value for 'species_col' not one of the columns of comm")

  # Drop species not matching between community and traits
  comm <- match_comm_and_traits(comm, traits, "comm")
  traits <- match_comm_and_traits(comm, traits, "traits")

  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, attributes(traits)[["Labels"]])))

  # Convert community to dataframe with rows as sites and columns as species
  comm_df <- tibble::column_to_rownames(comm, species_col) %>% t()

  # Run ses mpd
  picante::ses.mpd(
    samp = comm_df,
    dis = traits,
    null.model = null.model,
    iterations = iterations,
    runs = runs)

}

#' Analyze Standard Effect Size (SES) of functional nearest mean taxonomic distance (MNTD)
#'
#' Although at least >3 species must match between comm and tree, any non-matching
#' species will be dropped before running the analysis.
#'
#' @param comm Community data frame, with one column for species and
#' the rest for sites (i.e., rows as species and columns as sites).
#' @param traits Distance matrix of traits
#' @param species_col Name of column with species names
#' @param null.model Type of null model to use for picante::ses.mpd
#' @param iterations Number of iterations to use for "independentswap"
#' null model
#' @param runs Number of times to perform the randomization
#'
#' @return Dataframe. Results of ape::ses.mpd
#'
ses_func_mntd <- function(comm, traits, species_col = "species",
                         null.model = "independentswap",
                         iterations = 10000,
                         runs = 999) {

  assertthat::assert_that(
    species_col %in% colnames(comm),
    msg = "value for 'species_col' not one of the columns of comm")

  # Drop species not matching between community and traits
  comm <- match_comm_and_traits(comm, traits, "comm")
  traits <- match_comm_and_traits(comm, traits, "traits")

  # Make sure that worked
  assert_that(isTRUE(all.equal(comm$species, attributes(traits)[["Labels"]])))

  # Convert community to dataframe with rows as sites and columns as species
  comm_df <- tibble::column_to_rownames(comm, species_col) %>% t()

  # Run ses mpd
  picante::ses.mntd(
    samp = comm_df,
    dis = traits,
    null.model = null.model,
    iterations = iterations,
    runs = runs)

}

#' Clean up output from ses.mpd and ses.mntd
#'
#' @param ses_mpd_results Dataframe; output of picante::ses.mpd or
#' picante::ses.mntd
#' @param id Name of column to assign for rownames
#' @param cols_keep Column names from original output to keep
#' @param prefix String to append to column names
#'
#' @return Tibble
#'
clean_ses <- function (
  ses_mpd_results,
  prefix = "phy_",
  id = "secondary_grid_code"
) {

  ses_mpd_results <-
    select(ses_mpd_results, contains("obs.z"))

  colnames(ses_mpd_results) <- paste0(prefix, colnames(ses_mpd_results))

  ses_mpd_results %>%
    rownames_to_column(id) %>%
    as_tibble %>%
    clean_names

}

# Geospatial ----

#' Exclude points in Japan from GBIF data
#'
#' @param gbif_points_global Dataframe; global occurrences of pteridophytes,
#' with species, latitude, and longitude
#' @param all_cells Dataframe; centroids of 10km grid cells in Japan
#'
#' @return Tibble; global occurrences of pteridophytes, excluding any that
#' fall in the 10km grid cells in Japan
#'
exclude_japan_points <- function (gbif_points_global, all_cells) {

  # Convert global occurrence points to sf object.
  gbif_points_global_sf <-
    purrr::map2(
      gbif_points_global$decimallongitude,
      gbif_points_global$decimallatitude,
      ~ sf::st_point(c(.x, .y))) %>%
    sf::st_sfc(crs = 4326) %>%
    sf::st_sf(gbif_points_global, .)

  # Convert Japan 10 km grid cell centroids to sf object
  all_cells_sf <- purrr::map2(
    all_cells$longitude,
    all_cells$latitude,
    ~ sf::st_point(c(.x, .y))) %>%
    sf::st_sfc(crs = 4326) %>%
    sf::st_sf(all_cells, .)

  # Find all points within 10 km of centers of grid cells in Japan
  points_within_japan <- sf::st_is_within_distance(
    gbif_points_global_sf, all_cells_sf, 10)

  # Make logical vector of all points outside of Japan
  is_outside_of_japan <- lengths(points_within_japan) == 0

  # Exclude points inside Japan
  dplyr::filter(gbif_points_global, is_outside_of_japan)

}

#' Combine the global occurrence matrix with local (Japan) matrix
#'
#' @param comm_for_ecos_global_cropped Global pteridophyte occurrence matrix
#' with Japan occurrences removed (cropped), of class matrix.
#' @param comm_pteridos_renamed Occurrences of pteridophytes in Japan as tibble.
#'
#' @return Matrix
#'
combine_presabs_mat <- function(comm_for_ecos_global_cropped, comm_pteridos_renamed) {
  
  # Convert the cropped (i.e., without Japan) global matrix to tibble
  comm_for_ecos_global_cropped_tibble <-
    comm_for_ecos_global_cropped %>% 
    as.data.frame %>%
    rownames_to_column("site") %>%
    as_tibble
  
  # Find common columns betwen the two ('site' and species)
  common_cols <- intersect(
    colnames(comm_for_ecos_global_cropped_tibble),
    colnames(comm_pteridos_renamed)
  )
  
  # Combine the cropped global dataset (1-degree grid cells)
  # with the Japan dataset (10km x 10 km grid cells)
  combined <-
  bind_rows(
    select(comm_for_ecos_global_cropped_tibble, common_cols),
    select(comm_pteridos_renamed, common_cols)
  ) %>%
    # Convert to matrix
    column_to_rownames("site") %>%
    as.matrix()
  
  # Only keep those cells with at least one species, and species with
  # at least on occurence
  combined[rowSums(combined) > 0, colSums(combined) > 0]
  
}

# Ecostructure ----

#' Make community matrix from species' occurrences
#'
#' A grid is defined according to xmn, xmx, ymn, ymx, and reso of square
#' cells where x and y values correspond to latitude and longitude. This
#' is then filled in with the species occurrences, and converted to a
#' community matrix where the rows are cells and columns are species.
#'
#' @param species_coods Dataframe of species occurrences. Must include
#' columns "species", "decimallongitude", and "decimallatitude".
#' No missing values allowed.
#' @param xmn Minimum longitude used to construct the presence-absence grid.
#' @param xmx Maximum longitude used to construct the presence-absence grid.
#' @param ymn Minimum latitude used to construct the presence-absence grid.
#' @param ymx Maximum latitude used to construct the presence-absence grid.
#' @param reso Degree of spatial resolution used to construct the presence-absence grid.
#' @param crs Character or object of class CRS. PROJ.4 type description of a
#' Coordinate Reference System (map projection) used to construct the
#' presence-absence grid.
#' @param rownames Logical; should the rownames of the matrix be the
#' coordinates of each grid cell? If FALSE, separate columns will be
#' created for grid cell coordinates. Warning: setting rownames on a large
#' matrix may require a large amount of memory and is not recommended.
#' @param abun Logical; should the cells of the matrix be species' abundances?
#' If false, values of the cell will indicate presence (1) or absence (0) of
#' each species.
#'
#' @return Matrix
#'
#' @examples
#' # Make test data set of 10 species with 100 occurrences total.
#' test_data <- data.frame(species = letters[1:10],
#'   decimallongitude = runif(100, -2, 2),
#'   decimallatitude = runif(100, -2, 2),
#'   stringsAsFactors = FALSE)
#' # Presence/absence
#' comm_from_points(test_data)
#' # Abundance
#' comm_from_points(test_data, abun = TRUE)

comm_from_points <- function(species_coods,
                             xmn = -180, xmx = 180,
                             ymn = -90, ymx = 90,
                             resol = 1,
                             crs = sp::CRS("+proj=longlat +datum=WGS84"),
                             rownames = FALSE,
                             abun = FALSE) {

  checkr::check_data(
    species_coods,
    values = list(
      species = "a",
      decimallongitude = 1,
      decimallatitude = 1
    )
  )

  assertr::verify(species_coods, decimallongitude <= 180, success_logical)

  assertr::verify(species_coods, decimallatitude <= 90, success_logical)

  assertthat::assert_that(!is.factor(species_coods$species))

  # Create empty raster
  r <- raster::raster(
    resolution = resol,
    xmn = xmn,
    xmx = xmx,
    ymn = ymn,
    ymx = ymx,
    crs = crs)

  ### Extract cell IDs from points by species

  # Coordinates must be long, lat for raster
  # (which assumes x, y position in that order)
  species_coods <- species_coods[,c("species", "decimallongitude", "decimallatitude")]

  # Extract cell IDs from points by species
  species_coods <- tidyr::nest(species_coods, data = c(decimallongitude, decimallatitude))

  # raster::cellFromXY needs data to be data.frame, not tibble
  species_coods <- dplyr::mutate(species_coods, data = purrr::map(data, as.data.frame))

  cells_occur <- purrr::map(
    species_coods$data,
    ~ raster::cellFromXY(xy = ., object = r)
  )

  names(cells_occur) <- species_coods$species

  rm(species_coods)

  # Get vector of non-empty cells
  non_empty_cells <- sort(unique(unlist(cells_occur)))

  # Convert list of cell occurrences from cell IDs to vector
  # the length of non_empty_cells with 1 (or number of times for abundance)
  # for each cell where that species occurs.

  if(isTRUE(abun)) {
    # abudance
    cells_occur <- purrr::map(cells_occur, count_abun, all_cells = non_empty_cells)
  } else {
    # presence/absence
    cells_occur <- purrr::map(cells_occur, ~as.numeric(non_empty_cells %in% .))
  }

  # Combine these into a matrix
  pres_abs_mat <- t(do.call(rbind, cells_occur))

  rm(cells_occur)

  # Add coordinates

  # Make vector of xy (lat/long) coordinates in raster,
  # named by cell ID.

  # When setting names, assumes cell ID naming order goes
  # from upper-left most cell
  # to bottom right-most cell (as in raster docs).
  cell_xy <- raster::xyFromCell(r, 1:raster::ncell(r))
  cell_xy <- cell_xy[non_empty_cells,]

  if(isTRUE(rownames)) {
    cell_xy <- paste(cell_xy[,"x"], cell_xy[,"y"], sep = "_")
    rownames(pres_abs_mat) <- cell_xy
  } else {
    colnames(cell_xy) <- c("Longitude(x)", "Latitude(y)")
    pres_abs_mat <- cbind(cell_xy, pres_abs_mat)
  }

  return(pres_abs_mat)

}


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

#' Convert a presence/absence matrix to a list of
#' dispersion fields.
#' 
#' A dispersion field is the global richness for the species from 
#' one particular site.
#'
#' @param pres_ab_matrix Presence/absence matrix. Must be formatted
#' with species as columns and sites as rows. Rownames must be formatted
#' as the longitude and latitude of the centroid of the site separated
#' by an underscore, e.g. '160_20'.
#' @param proj Coordinate reference system
#' @param xmin Minimum longitude to use when projecting the dispersion fields
#' @param xmax Maximum longitude to use when projecting the dispersion fields
#' @param ymin Minimum latitude to use when projecting the dispersion fields
#' @param ymax Maximum latitude to use when projecting the dispersion fields
#'
#' @return List, of length equal to the number of sites (rows) in the
#' input presence/absence matrix. 
#' Each item in the list is a raster -- the dispersion field for that site.
#' 
#' @examples
#' test_matrix <- matrix(data = rbinom(100, 1, 0.5), 10, 10)
#' longs <- sample(c(1:5), replace=TRUE, size=10)
#' lats <- sample(c(1:5), replace=TRUE, size=10)
#' rownames(test_matrix) <- paste(longs, lats, sep = "_")
#' colnames(test_matrix) <- letters[1:10]
#' disp_list <- pres_ab_to_disp(test_matrix, 1,5,1,5, res = 1)
#' test_matrix
#' plot(disp_list[[1]])
pres_ab_to_disp <- function (pres_ab_matrix,
                             xmin = -180, xmax = 180,
                             ymin = -90, ymax = 90,
                             proj = sp::CRS(' +proj=longlat +ellps=WGS84'),
                             res = 1) {
  
  assertthat::assert_that(is.matrix(pres_ab_matrix))
  
  # Make empty raster
  ras <- raster::raster(xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj, resolution = res)
  
  # Prepare loop to get cell IDs
  idx <- c()
  site_names <- rownames(pres_ab_matrix)
  
  # `idx` is a vector of cell IDs in the raster, only including cells in the
  # pres/abs matrix
  for (i in 1:length(site_names)){
    idx[i] <- raster::cellFromXY(ras, as.numeric(unlist(strsplit(site_names[i], "[_]"))))
  }
  
  # Prepare loop for species rasters
  species_rasters <- list()
  k <- length(colnames(pres_ab_matrix))
  
  # `species_rasters` is a list of rasters, one for each species in the
  # pres/abs matrix
  for (i in 1:k){
    ras <- raster::raster(xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj, resolution = res)
    ras[idx] <- pres_ab_matrix[,i]
    species_rasters[[i]] <- ras
  }
  
  names(species_rasters) <- colnames(pres_ab_matrix)
  
  # Prepare loop for dispersion fields
  disp_rasters <- list()
  k <- length(site_names)
  
  # `disp_rasters` is a list of rasters; each one is the dispersion field
  # of a site in the pres-abs matrix, i.e., the global richness of species
  # at that site.
  for (i in 1:k){
    spec_in_cell <- names(pres_ab_matrix[i,][pres_ab_matrix[i,] > 0])
    if (length(spec_in_cell) < 2){
      disp_rasters[[i]] <- species_rasters[spec_in_cell][[1]]
    } else {
      disp_rasters[[i]] <- sum(raster::stack(species_rasters[names(pres_ab_matrix[i,][pres_ab_matrix[i,] > 0])]))
    }
  }
  
  names(disp_rasters) <- rownames(pres_ab_matrix)
  
  return(disp_rasters)
  
}

#' Convert a list of dispersion fields to a matrix
#'
#' @param dispersion.field List of dispersion fields, 
#' where each dispersion field is a raster.
#' @param 
#'
#' @return Matrix of with number of rows equal to the length
#' of `dispersion.field` and number of columns equal to the 
#' number of rows x the number of columns of each dispersion
#' field in the list.
#' 
#' @examples
#' test_matrix <- matrix(data = rbinom(100, 1, 0.5), 10, 10)
#' longs <- sample(c(1:5), replace=TRUE, size=10)
#' lats <- sample(c(1:5), replace=TRUE, size=10)
#' rownames(test_matrix) <- paste(longs, lats, sep = "_")
#' colnames(test_matrix) <- letters[1:10]
#' disp_list <- pres_ab_to_disp(test_matrix, 1,5,1,5, res = 1)
#' dsp_to_matrix2(disp_list, drop_zero = TRUE)
dsp_to_matrix2 <- function (dispersion.field, drop_zero = FALSE) {
  
  # Set up results matrix: number of rows equal to length of dispersion
  # field list, with number of columns equal to the total number of cells
  # in each dispersion field.
  map_data <- matrix(
    ncol = dim(dispersion.field[[1]])[1] * dim(dispersion.field[[1]])[2], 
    nrow=length(dispersion.field))
  
  # Convert dispersion field list to matrix.
  for (l in 1:length(dispersion.field)) {
    temp_data <- dispersion.field[[l]]
    temp_data[is.na(temp_data)] <- 0
    map_data[l,] <- as.vector(temp_data)
  }
  
  # Set row and column names.
  rownames(map_data) <- names(dispersion.field)
  
  lat_long_names <- c()
  for (i in 1:dim(map_data)[2]){
    lat_long_names[i] <- paste0(raster::xyFromCell(temp_data,i), collapse = "_")
  }
  colnames(map_data) <- lat_long_names
  
  # Optionally drop columns that are all zeros
  if (isTRUE(drop_zero)) {
    map_data <- map_data[,colSums(map_data) != 0]
  }
  
  return(map_data)
}

#' Fit an ecostructure model and save the output as RDS
#' 
#' The output RDS file will be saved to `out_dir` and named
#' automatically according to motif, dataset, and other settings
#' used for ecostructure::ecos_fit().
#'
#' @param out_dir Directory to save the results
#' @param motif Type of motif used for ecostructure. Must be one of
#' "species", "species-trans", "geo", "phylo", or "traits"
#' @param dataset Name of dataset used
#' 
#' Other arguments as in ecostructure::ecos_fit()
#'
#' @return NULL
#' @examples
#' library(ecostructure)
#' data("himalayan_birds")
#' species_counts <- t(Biobase::exprs(himalayan_birds))
#' run_ecos(out_dir = "tmp", motif = "species", dataset = "birds",
#' dat = species_counts, K = 2, tol = 0.1, num_trials = 1)
#' 
run_ecos <- function (out_dir, motif, dataset, 
                      dat, max_dat = NULL, K, tol = 0.1, num_trials = 1, 
                      fit_control = list()) {
  
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(out_dir)))
  
  out_file = glue::glue("ecosfit_motif={motif}_dataset={dataset}_K={K}_tol={tol}_trials={num_trials}.RDS")
  
  out_path = fs::path(out_dir, out_file)
  
  results <- ecostructure::ecos_fit(dat = dat, max_dat = max_dat, 
                                    K = K, tol = tol, num_trials = num_trials, 
                                    fit_control = fit_control)
  
  saveRDS(object = results, file = out_path)
  
  results
  
}

# Plotting ----

#' Define ggplot theme
#'
#' BW theme with no gridlines, black axis text, main font size 11,
#' axis ticks size 9.
#'
standard_theme2 <- function () {
  ggplot2::theme_bw() +
    theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour="black"),
      axis.text.y = ggplot2::element_text(colour="black")
    )
}

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
make_diversity_map <- function (div_data, world_map, occ_data, div_metric, metric_title, label) {

  div_metric <- sym(div_metric)

  ggplot(world_map, aes(x = longitude, y = latitude)) +
    geom_polygon(aes(group = group), fill = "light grey") +
    geom_tile(data = div_data,
              aes(fill = !!div_metric)) +
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor,
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor,
               pull(occ_data, latitude) %>% max %>% ceiling)
    )  +
    labs(
      fill = metric_title
    ) +
    annotate(
      "text", label = label, size = 12/.pt, fontface = "bold",
      x = -Inf,
      y = Inf,
      vjust = 1.2, hjust = -0.5) +
    jntools::blank_x_theme() +
    jntools::blank_y_theme() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_text(size = 20/.pt),
      legend.text = element_text(size = 16/.pt),
      legend.justification=c(1,0),
      legend.position=c(1.1,0))
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
    standard_theme2() +
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


#' Plot an aspect of alpha diversity by latitude and elevation
#'
#' @param alpha_div Tibble of alpha diversity metrics including
#' elevatin and latitude
#' @param var Name of variable to plot
#' @param subtitle Subtitle to add to plot
#'
#' @return GGplot object
#'
make_lat_el_pd_plot <- function (alpha_div, var, subtitle) {

  var_sym <- rlang::sym(var)

  ggplot(alpha_div, aes(x = elevation, y = latitude, color = !!var_sym)) +
    geom_point(size = 0.7, alpha = 0.7, stroke=0) +
    scale_color_scico(
      palette = "vik",
      na.value="dark grey",
      limits = c(
        -get_limit(alpha_div, !!var_sym, "abs", 3),
        get_limit(alpha_div, !!var_sym, "abs", 3)
      )) +
    labs(
      color = "SES",
      subtitle = subtitle) +
    standard_theme2()

}

#' Plot an aspect of alpha diversity by latitude and elevation
#'
#' @param alpha_div Tibble of alpha diversity metrics including
#' elevation and latitude
#' @param var Name of variable to plot
#' @param legend_label Label to use for legend
#' @param subtitle Subtitle to add to plot
#'
#' @return GGplot object
#'
make_lat_el_rich_plot <- function (alpha_div, var, legend_label, subtitle) {

  var_sym <- rlang::sym(var)

  ggplot(alpha_div, aes(x = elevation, y = latitude, color = !!var_sym)) +
    geom_point(size = 0.7, alpha = 0.7, stroke=0) +
    scale_color_scico(palette = "bamako", na.value="grey") +
    labs(
      color = legend_label,
      subtitle = subtitle) +
    standard_theme2()

}

#' Creat a set of elevation by latitude plots for a given dataset
#'
#' @param alpha_div Tibble of alpha diversity metrics including
#' elevation and latitude
#' @param main_title Main title to use for plot
#'
#' @return GGplot object
#'
compose_lat_el_plots <- function (alpha_div, main_title) {

  # Make tibble of variables to feed into plotting funcs
  ses_plot_vars <- tibble(
    var = c("phy_mpd_obs_z", "phy_mntd_obs_z", "func_mpd_obs_z", "func_mntd_obs_z"),
    subtitle = c("MPDphy", "MNTDphy", "MPDfunc", "MNTDfunc")
  ) %>%
    mutate(alpha_div = list(alpha_div))

  richness_plot_vars <- tibble(
    var = c("richness", "percent_sex_dip"),
    subtitle = c("Richness", "Percent sex. dip."),
    legend_label = c("n_spp", "%")
  ) %>%
    mutate(alpha_div = list(alpha_div))

  ses_plots <- pmap(ses_plot_vars, make_lat_el_pd_plot)

  richness_plots <- pmap(richness_plot_vars, make_lat_el_rich_plot)

  ses_plots[[1]] <- ses_plots[[1]] +
    labs(title = main_title)

  wrap_plots(c(ses_plots,richness_plots), ncol = 2, nrow = 3)

}

# Drake ----

# Fill-in memory settings as named lists

#' Specify queue and memory requested for a drake plan
#' running on the hydra cluster
#'
#' @param plan Drake plan. Must include column for "mem_type", which
#' includes values "small", "medium", "large", or "very large".
#' The amount of memory and run queue will be set accordingly
#' (5gb, 50gb, 100gb, or 150gb per job).
#'
#' @return Tibble
#' 
specify_resources <- function (plan) {
  
  plan$resources <- ""
  
  for (i in 1:nrow(plan)) {
    # Default is 1gb, short-run time if not specified
    plan$resources[i] <- list(list(queue = "sThC.q", memory = "mres=1G"))
    # Change if mem_type is specified
    if (isTRUE(plan$mem_type[[i]] == "small")) {
      plan$resources[i] <- list(
        list(queue = "lThC.q", memory = "mres=5G,h_data=5G,h_vmem=5G") # 24 hours, 5Gb
      ) }
    if (isTRUE(plan$mem_type[[i]] == "medium")) {
      plan$resources[i] <- list(
        list(queue = "uThM.q -l lopri", memory = "mres=50G,h_data=50G,h_vmem=50G,himem") # unlimited, 50Gb
      ) }
    if (isTRUE(plan$mem_type[[i]] == "large")) {
      plan$resources[i] <- list(
        list(queue = "uThM.q -l lopri", memory = "mres=100G,h_data=100G,h_vmem=100G,himem") # unlimited, 100Gb
      ) }
    if (isTRUE(plan$mem_type[[i]] == "very large")) {
      plan$resources[i] <- list(
        list(queue = "uThM.q -l lopri", memory = "mres=150G,h_data=150G,h_vmem=150G,himem") # unlimited, 150Gb
      ) }
  }
  
  plan
  
}

