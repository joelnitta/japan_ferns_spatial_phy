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
    rename(
      rbcl_genbank_no = rbc_l_gen_bank_accession_no,
      species = taxon_name) %>%
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
    # Simplify taxon names, replace space with underscore
    taxastand::add_parsed_names(scientific_name, taxon) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
    mutate(taxon_id = as.character(taxon_id)) %>%
    assert(not_na, taxon, taxon_id, scientific_name) %>%
    assert(is_uniq, taxon, taxon_id, scientific_name)
}

#' Subset data to only ferns
#'
#' @param data Input dataframe with column "taxon"
#' @param ppgi Pteridophyte phylogeny group I taxonomy
#'
#' @return Dataframe subset to ferns only
#' 
subset_to_ferns <- function(data, ppgi) {
  
  cols_keep <- colnames(data)
  
  data %>%
    mutate(genus = str_split(taxon, " |_") %>% map_chr(1)) %>%
    assert(not_na, genus) %>%
    # Add higher-level taxonomy
    left_join(ppgi, by = "genus") %>%
    assert(not_na, class) %>%
    # Filter to only ferns
    filter(class == "Polypodiopsida") %>%
    select(all_of(cols_keep))
  
}

#' Subset occurrence point data to only ferns
#'
#' @param phy Phylogeny
#' @param ppgi Pteridophyte phylogeny group I taxonomy
#'
#' @return Phylogeny subset to ferns only
#' 
subset_tree <- function(phy, ppgi) {
  
  tips_keep <-
  tibble(tip = phy$tip.label) %>%
    mutate(genus = str_split(tip, "_") %>% map_chr(1)) %>%
    # Add higher-level taxonomy
    left_join(ppgi, by = "genus") %>%
    assert(not_na, class) %>%
    # Filter to only ferns
    filter(class == "Polypodiopsida") %>%
    # Check for missing data
    assert(not_na, tip) %>%
    assert(is_uniq, tip) %>%
    pull(tip)
  
  ape::keep.tip(phy, tips_keep)
  
}

#' Rename taxa in data by taxon ID code
#'
#' The `green_list` has the official taxon names. These
#' will be mapped by taxon ID code to the input data,
#' and offical names used in replacement of original species names.
#'
#' @param data Dataframe with columns `taxon_id` and `species`
#' @param green_list Dataframe of standard taxonomy with
#' columns `taxon_id` and `taxon`
#'
#' @return Dataframe
#' 
rename_taxa <- function (data, green_list) {
  
  data %>%
    left_join(select(green_list, taxon_id, taxon), by = "taxon_id") %>%
    assert(not_na, taxon) %>%
    select(-species, -taxon_id)
  
}

#' Subset a community to only endemic taxa
#'
#' @param comm Community matrix with colnames as species
#' @param green_list Dataframe with columns `taxon` and `endemic`
#'
#' @return Subsetted community matrix
#' 
subset_comm_to_endemic <- function (comm, green_list) {
  
  endemic_taxa <-
    green_list %>%
    filter(!is.na(endemic)) %>%
    verify(all(endemic == "Endemic")) %>%
    pull(taxon)
  
  comm[,colnames(comm) %in% endemic_taxa]
  
}

#' Extract community dataframe from points2comm()
#' 
#' Also converts the community dataframe to presence/absence
#'
#' @param data Output of phyloregion::points2comm()
#'
#' @return Dataframe
#' 
comm_from_points2comm <- function (data) {
  
  # Convert input to "dense" dataframe
  data[["comm_dat"]] %>%
    phyloregion::sparse2dense() %>% 
    as.data.frame() %>%
    # Temporarily store rownames in "site" column so tidyverse doesn't obliterate them
    rownames_to_column("site") %>%
    # Convert to presence-absence
    mutate_if(is.numeric, ~ifelse(. > 0, 1, 0)) %>%
    column_to_rownames("site") %>%
    # Make sure everything is 0-1
    assert(in_set(c(0,1)), everything()) %>%
    assert(not_na, everything())
  
}

#' Extract shape dataframe from points2comm()
#'
#' @param data Output of phyloregion::points2comm()
#'
#' @return Dataframe (sf object)
#' 
shape_from_points2comm <- function (data) {
  
  data[["poly_shp"]] %>% sf::st_as_sf()
  
}

#' Filter occurrence points using a polygon mask
#'
#' @param occ_point_data Occurrence points dataframe,
#' including columns for taxon, longitude, and latitude
#' @param shape_file Path to shape file to use as mask
#'
#' @return Filtered data points
filter_occ_points <- function(occ_point_data, shape_file) {
  
  # Remove duplicates (same species collected on same day at same site)
  occ_point_data <-
    occ_point_data %>%
    # only about 3,000 have date missing. conservatively consider these the same day.
    mutate(date = replace_na(date, "missing")) %>%
    assert(not_na, taxon, longitude, latitude, date) %>%
    group_by(taxon, longitude, latitude, date) %>%
    summarize(
      n = n(),
      .groups = "drop"
    )
  
  
  # Read in shape file to use for mask
  # (here, the second-degree mesh map of Japan)
  second_degree_mesh <- sf::st_read(shape_file) %>%
    select(id = NAME)
  
  # Convert point data to SF object,
  # with same projection as shape file
  occ_point_data <- sf::st_as_sf(occ_point_data, coords = c("longitude", "latitude"), crs = st_crs(second_degree_mesh))
  
  # Join point data to mesh map
  sf::st_join(occ_point_data, second_degree_mesh, join = sf::st_within) %>%
    # Filter to only those with a grid ID (so, those that are within the mesh map)
    filter(!is.na(id)) %>%
    # Convert to tibble with columns for taxon, longitude, and latitude
    as_tibble() %>%
    mutate(coords = sf::st_coordinates(geometry)) %>%
    select(taxon, coords) %>%
    mutate(
      coords = as.data.frame(coords),
      longitude = coords$X,
      latitude = coords$Y) %>%
    select(-coords)
  
}

#' Filter a community dataframe by sampling redundancy
#' 
#' Only those communities exceeding the minimum sampling redundancy
#' will be kept
#'
#' @param comm Community data matrix, with species as columns and
#' rows as sites (rownames are site names)
#' @param shape Spatial dataframe including columns for "grids" (community
#' site name) and redundancy
#' @param cutoff Redundancy cutoff value to use
#'
#' @return Dataframe
#' 
filter_comm_by_redun <- function (comm, shape, cutoff = 0.1) {
  
  grids_keep <-
    shape %>%
    filter(redundancy > cutoff) %>%
    pull(grids)
  
  comm[rownames(comm) %in% grids_keep,]
  
}

#' Combine shapes of protected areas in Japan
#'
#' @param ... Shape files read in with sf::st_read()
#'
#' @return Simple feature collection
#' 
combine_protected_areas <- function (...) {
  
  bind_rows(...) %>%
    assert(not_na, status) %>%
    # Remove all marine areas
    filter(status != "marine") %>%
    # Add area
    mutate(area = st_area(.) %>% units::set_units(km^2)) %>%
    # Convert status to ordered factor
    mutate(
      status = factor(status, ordered = TRUE, levels = c("low", "medium", "high"))
    )
  
}

# Reproductive mode ----

#' Calculate percent of sexual diploid taxa per grid cell (site)
#'
#' @param comm Community data matrix, with species as columns
#' and one column indicating "site"
#' @param repro_data Reproductive data, including columns 
#' for reproductive type (`sexual_diploid`, `sexual_polyploid`)
#'
#' @return Tibble with percentage of sexual diploids per cell
#' 
calc_sex_dip <- function (comm, repro_data) {
  
  # First subset comm and repro_data to the species in common
  taxa_keep <- intersect(colnames(comm), repro_data$taxon)
  
  n_taxa_to_drop_from_comm <- setdiff(colnames(comm), taxa_keep) %>% length()
  
  if (n_taxa_to_drop_from_comm > 0) message (glue::glue("Dropping {n_taxa_to_drop_from_comm} species in community missing from reproductive data"))
  
  comm <- comm[,taxa_keep]
  
  comm %>%
    rownames_to_column("site") %>%
    as_tibble() %>%
    # Convert comm to long
    pivot_longer(names_to = "taxon", values_to = "abundance", -site) %>% 
    filter(abundance > 0) %>%
    # Join reproductive data
    left_join(repro_data, by = "taxon") %>%
    # Calculate percent sexual diploid per site
    assert(not_na, sexual_diploid, taxon) %>%
    group_by(site) %>%
    summarize(
      num_sex_dip = sum(sexual_diploid), # `sexual_diploid` is logical; TRUE for sex diploids, FALSE otherwise
      num_total = n(),
      .groups = "drop"
    ) %>%
    ungroup %>%
    mutate(percent_sex_dip = num_sex_dip / num_total)
  
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
    mutate(genus = str_split(taxon_name, " |_") %>% map_chr(1)) %>%
    filter(!is.na(genus)) %>%
    left_join(taxonomy_data, by = "genus")
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
      assertr::assert(is.numeric, all_of(trans_select)) %>%
      # Replace zeros with arbitrarily small number
      dplyr::mutate_at(all_of(trans_select), ~ifelse(. == 0, small_number, .)) %>%
      dplyr::mutate_at(all_of(trans_select), log)
  }
  
  # Rescale by dividing original value by range of
  # that value (max - min) across the dataset
  if (scale_traits == TRUE) {
    traits <-
      traits %>%
      assertr::verify(scale_select %in% colnames(traits)) %>%
      assertr::assert(is.numeric, all_of(scale_select)) %>%
      dplyr::mutate_at(
        all_of(scale_select), ~ . / (max(., na.rm = TRUE) - min(., na.rm = TRUE))
      )
  }
  
  traits
  
}

#' Format traits for further analysis
#'
#' Use raw fern and lycophyte trait data
#' formatted for lucid.
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
  traits <- read_excel(path_to_lucid_traits, skip = 1) %>%
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
    )) %>%
    assert(in_set(c(0,1)), -taxon)
  
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
    group_by(taxon, trait) %T>% 
    # Normally this next step would issue a warning because a bunch of
    # values have NA for both sterile and fertile, resulting in -Inf.
    # But we fix that in the next step anyways.
    {options(warn=-1)} %>%
    summarize(
      value = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup %>%
    mutate(value = na_if(value, -Inf)) %>%
    spread(trait, value) %>%
    rename_all(~str_remove(., "_number"))
  
  # Log-transform and scale numeric traits
  num_trait_names <- select(traits_numeric_combined, -taxon) %>% colnames()
  
  traits_numeric_combined_trans <- transform_traits(
    traits_numeric_combined,
    trans_select = all_of(num_trait_names),
    scale_select = all_of(num_trait_names)
  )
  
  # Combine all numeric and categorical traits
  traits_for_dist <- left_join(traits_numeric_combined_trans, traits_binary, by = "taxon")
  
  # Fix some taxon names (synonyms)
  traits_for_dist <-
  traits_for_dist %>%
    mutate(taxon = case_when(
      # taxon == "Athyrium_nudum" ~ missing
      taxon == "Athyrium_opacum_opacum" ~ "Athyrium_opacum",
      # taxon == "Azolla_cristata" ~ , missing
      # taxon == "Dryopteris_intermedia" ~ , missing
      # taxon == "Pityrogramma_calomelanos" ~ , non-native?
      # taxon == "Psilotum_nudum" ~ , non-native?
      # taxon == "Selaginella_uncinata" ~ , missing
      taxon == "Stegnogramma_griffithii_wilfordii" ~ "Thelypteris_griffithii_wilfordii",
      taxon == "Stegnogramma_gymnocarpa_amabilis" ~ "Thelypteris_gymnocarpa_amabilis",
      taxon == "Stegnogramma_pozoi_mollissima" ~ "Thelypteris_pozoi_mollissima",
      taxon == "Thelypteris_aurita" ~ "Phegopteris_aurita",
      taxon == "Thelypteris_bukoensis" ~ "Phegopteris_bukoensis",
      taxon == "Thelypteris_decursivepinnata" ~ "Phegopteris_decursivepinnata",
      taxon == "Thelypteris_nipponica" ~ "Thelypteris_nipponica_nipponica",
      taxon == "Thelypteris_ogasawarensis" ~ "Macrothelypteris_ogasawarensis",
      taxon == "Thelypteris_phegopteris" ~ "Phegopteris_connectilis",
      taxon == "Thelypteris_subaurita" ~ "Phegopteris_subaurita",
      taxon == "Thelypteris_torresiana_calvata" ~ "Macrothelypteris_torresiana_calvata",
      taxon == "Thelypteris_torresiana_torresiana" ~ "Macrothelypteris_torresiana_torresiana",
      taxon == "Thelypteris_viridifrons" ~ "Macrothelypteris_viridifrons",
      TRUE ~ taxon
    ))
  
  # Convert species names to taxon id codes
  missing_taxon_id <-
    filter(traits_for_dist, !(taxon %in% taxon_id_map$taxon)) %>%
    pull(taxon)
  
  # Drop any missing names
  traits_for_dist <- filter(traits_for_dist, taxon %in% taxon_id_map$taxon)
    
  msg <- assertthat::validate_that(
    length(missing_taxon_id) == 0,
    msg = glue::glue("The following taxa in the trait data could not be verified in the taxa list and have been dropped: {paste(missing_taxon_id, collapse = ', ')}")
  )
  
  if(is.character(msg)) message(msg)
  
  ### Find correlated traits to drop ###
  
  # helper function to identify static traits
  get_static_traits <- function (data) {
    data %>%
      pivot_longer(names_to = "trait", values_to = "value", -taxon) %>%
      group_by(trait) %>%
      count(value) %>%
      ungroup %>%
      filter(!is.na(value)) %>%
      select(trait) %>%
      count(trait) %>%
      filter(n < 2) %>%
      pull(trait)
  }
  
  # helper function to identify low-variance traits
  # (less than two occurrences of a given trait state)
  get_low_var_traits <- function (data) {
    data %>% 
      pivot_longer(names_to = "trait", values_to = "value", -taxon) %>%
      group_by(trait) %>%
      count(value) %>%
      filter(!is.na(value)) %>%
      filter(n < 3) %>%
      pull(trait)
  }
  
  # To calculate Pearson's correlation co-efficient, 
  # can't allow any missing traits, so use a subset of the data
  traits_corr_test <- ggplot2::remove_missing(traits_for_dist)
  
  # First make vector of any static trait (only a single trait state)
  traits_corr_static <- 
    traits_corr_test %>%
    select(all_of(colnames(traits_binary))) %>% # consider binary traits only
    get_static_traits
  
  # Next make a vector of any trait with less than two occurrences of a given trait state
  traits_corr_low_var <-
    traits_corr_test %>%
    select(all_of(colnames(traits_binary))) %>% # consider binary traits only
    get_low_var_traits
  
  # Drop the low and non-varying traits, drop the "taxon" column
  traits_corr_test <- select(traits_corr_test, -any_of(c(traits_corr_static, traits_corr_low_var))) %>%
    select(-taxon)
  
  # Calculate correlations
  correlations <- cor(traits_corr_test)
  
  # Select columns to remove with absolute correlation > 0.6
  traits_correlated_to_drop <- caret::findCorrelation(correlations, cutoff = 0.6) %>%
    magrittr::extract(colnames(traits_corr_test), .)
  
  ### Drop zero and low-variance traits from full dataframe ###
  
  # Vector of any static trait (only a single trait state)
  traits_static <- 
    traits_for_dist %>%
    select(all_of(colnames(traits_binary))) %>% # consider binary traits only
    get_static_traits
  
  # Next make a vector of any trait with less than two occurrences of a given trait state
  traits_low_var <-
    traits_for_dist %>%
    select(all_of(colnames(traits_binary))) %>% # consider binary traits only
    get_low_var_traits
  
  # Drop correlated and non-varying traits from full dataframe
  traits_to_drop <- c(traits_static, traits_low_var, traits_correlated_to_drop) %>% unique()
  traits_for_dist <- select(traits_for_dist, -any_of(traits_to_drop))
  
  # There are still a few traits left once we use the full data that are correlated.
  # Manually remove these.
  still_correlated <- c(
    "leaf_sorus_shape_globular",
    "leaf_sorus_indusium_shape_globular",
    "leaf_lamina_vennation_unbrached_single",
    "leaf_sorus_shape_not_forming_distinc_sorus"
  )
  
  traits_for_dist %>%
    select(-taxon) %>%
    corrr::correlate() %>%
    pivot_longer(-rowname, names_to = "var2") %>%
    rename(var1 = rowname) %>%
    arrange(desc(value)) %>%
    filter(value > 0.6) %>%
    verify(all(var1 %in% still_correlated)) %>%
    verify(all(var2 %in% still_correlated), success_fun = success_logical)
  
  traits_for_dist <- select(traits_for_dist, -leaf_sorus_shape_globular, -leaf_lamina_vennation_unbrached_single)
  
  # Verify that observed correlations in final data are less than 0.65
  traits_for_dist %>%
    select(-taxon) %>%
    corrr::correlate() %>%
    pivot_longer(-rowname) %>%
    assert(within_bounds(-0.65, 0.65), value, success_fun = success_logical)
  
  traits_for_dist
  
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
    select(-taxon) %>%
    colnames %>%
    tibble(trait = .) %>%
    left_join(states, by = "trait") %>%
    mutate(trait_type = case_when(
      n_states > 3 ~ "continuous",
      TRUE ~ "binary"
    )) %>%
    # Collapse binarized traits into their
    # categorical version by name.
    categorize_traits %>%
    # First weigh by combined trait
    add_count(comp_trait) %>%
    mutate(weight_by_comp_trait = 1 / n) %>%
    select(-n) %>%
    # Then weigh by trait type
    add_count(trait_type) %>%
    mutate(weight_by_type = 1 / n) %>%
    select(-n) %>%
    # Then weigh by combination of combined trait and trait type
    mutate(final_weight = weight_by_comp_trait * weight_by_type)
  
  # Check trait weights
  trait_categories$weight_by_comp_trait %>% sum
  trait_categories$weight_by_type %>% sum
  trait_categories$final_weight %>% sum
  
  # Make sure traits for calculating the distance matrix are in correct
  # order for weighting
  traits_for_dist <- select(traits_for_dist, taxon, all_of(trait_categories$trait))
  
  # Convert to dataframe for gowdis
  traits_df <- traits_for_dist %>%
    column_to_rownames("taxon")
  
  # Run gowdis with trait weightings
  FD::gowdis(traits_df, w = trait_categories$final_weight)
  
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

#' Calculate sampling coverage at each site of a community matrix
#' 
#' Coverage (ratio of species present in sample to that in true assemblage) 
#' estimated following method of Chao et al 2014 (doi: 10.1111/2041-210X.12247)
#'
#' @param comm_list Output of phyloregion::points2comm(), a list of two items:
#' - `comm_dat`: community data frame
#' - `poly_shp`: shapefile of grid cells with the values per cell.
#' (this function extracts and uses `comm_dat`)
#'
#' @return Dataframe, with columns for "site" and "coverage"
#' 
calculate_coverage <- function (comm_list) {
  
  # Extract community matrix and convert to data frame
  comm <-
  magrittr::extract2(comm_list, "comm_dat") %>%
    phyloregion::sparse2dense() %>%
    as.data.frame()
  
  # Convert data frame to transposed version, with rows as species
  # and columns as sites (for iNEXT)
  comm_t <-
    comm %>%
    rownames_to_column("site") %>%
    pivot_longer(values_to = "abundance", names_to = "species", -site) %>%
    pivot_wider(values_from = "abundance", names_from = "site") %>%
    column_to_rownames("species")
  
  # Run iNEXT and extract sampling coverage results
  iNEXT::iNEXT(comm_t, datatype = "abundance") %>% 
    magrittr::extract2("DataInfo") %>%
    select(site, coverage = SC)
  
}

# Biodiverse ----

# Reformat community matrix (columns as sites and rows as species) to
# format for biodiverse (rows as sites, columns as species, with lat/long)
make_matrix_for_biodiverse <- function (comm, cell_xy) {
  
  comm %>%
    pivot_longer(-species, names_to = "site", values_to = "count") %>%
    pivot_wider(names_from = species, values_from = "count") %>%
    left_join(cell_xy, by = "site") %>%
    select(site, longitude, latitude, everything())
  
}

#' Rename species in community matrix from number codes to species names
#'
#' @param comm Community matrix
#' @param taxon_id_map Tibble mapping species names to number codes
#'
#' @return Dataframe
#' 
rename_comm <- function (comm, taxon_id_map) {
  comm %>%
    left_join(taxon_id_map, by = c(species = "taxon_id")) %>%
    assert(not_na, taxon) %>%
    select(-species) %>%
    select(species = taxon, everything())
}

#' Rename tips in tree from number codes to species names
#'
#' @param tree Phylogenetic tree with tips as number codes
#' @param taxon_id_map Tibble mapping species names to number codes
#'
#' @return Phylogenetic tree
rename_tree <- function (tree, taxon_id_map) {
  
  new_tip_labs <- tree$tip.label %>%
    tibble(taxon_id = .) %>%
    left_join(taxon_id_map, by = "taxon_id") %>%
    assert(not_na, taxon) %>%
    assert(is_uniq, taxon)
  
  tree$tip.label <- new_tip_labs$taxon
  
  tree
  
}

#' Calculate diversity metrics for a single random community
#' 
#' The independent swap method of Gotelli (2000) is used, which randomizes
#' the community matrix while maintaining species occurrence frequency and
#' sample species richness.
#' 
#' For description of metrics available, see run_ses_analysis()
#'
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names)
#' @param phy Input phylogeny with total branch length scaled to 1
#' @param phy_alt Alternative phylogeny where all branches are of equal length, scaled to 1
#' @param n_iterations Number of iterations to use when shuffling random community
#' @param metrics Names of metrics to calculate. Must one or more of
#' 'pd', 'rpd', 'fd', 'rfd', 'pe', or 'rpe'
#'
#' @return List of vectors. Each vector is a biodiversity metric measured on the
#' random community, in the same order as the rows in the input community.
#' 
calc_biodiv_random <- function (comm_df, phy = NULL, phy_alt = NULL, trait_tree = NULL, trait_tree_alt = NULL, n_iterations, metrics) {
  
  # Make sure selection of metrics is OK
  assert_that(is.character(metrics))
  assert_that(
    length(metrics) > 0,
    msg = "At least one biodiversity metric must be selected")
  assert_that(
    all(metrics %in% c("pd", "rpd", "fd", "rfd", "pe", "rpe")),
    msg = "Biodiversity metrics may only be selected from 'pd', 'rpd', 'fd', 'rfd', 'pe', or 'rpe'"
  )
  
  # Make sure names match between community and tree
  if (any(metrics %in% c("pd", "pe"))) assert_that(isTRUE(
    all.equal(sort(phy$tip.label), sort(colnames(comm_df)))
  ))
  
  if (any(metrics %in% c("rpd", "rpe"))) assert_that(isTRUE(
    all.equal(sort(phy_alt$tip.label), sort(colnames(comm_df)))
  ))
  
  if (any(metrics %in% c("fd"))) assert_that(isTRUE(
    all.equal(sort(trait_tree$tip.label), sort(colnames(comm_df)))
  ))
  
  if (any(metrics %in% c("rfd"))) assert_that(isTRUE(
    all.equal(sort(trait_tree_alt$tip.label), sort(colnames(comm_df)))
  ))
  
  # Make sure phylogeny has been rescaled to total branch length of 1 for RPE or RFD
  if (any(metrics %in% c("rpe", "rpd"))) assert_that(isTRUE(all.equal(sum(phy$edge.length), 1)))
  if (any(metrics %in% c("rpe", "rpd"))) assert_that(isTRUE(all.equal(sum(phy_alt$edge.length), 1)))
  if (any(metrics %in% c("rfd"))) assert_that(isTRUE(all.equal(sum(trait_tree$edge.length), 1)))
  if (any(metrics %in% c("rfd"))) assert_that(isTRUE(all.equal(sum(trait_tree_alt$edge.length), 1)))
  
  # Convert comm to sparse matrix format for phyloregions
  comm_sparse <- phyloregion::dense2sparse(comm_df)
  
  # Generate random community
  random_comm <- picante::randomizeMatrix(comm_df, null.model = "independentswap", iterations = n_iterations)
  random_comm_sparse <- phyloregion::dense2sparse(random_comm)
  
  # Calculate statistics for random community
  # - set up null vectors first
  pd <- NULL
  pd_alt <- NULL
  rpd <- NULL
  fd <- NULL
  fd_alt <- NULL
  rfd <- NULL
  pe <- NULL
  pe_alt <- NULL
  rpe <- NULL
  
  # - calculate selected metrics
  if ("pd" %in% metrics) pd <- phyloregion::PD(random_comm_sparse, phy)
  if ("rpd" %in% metrics) pd_alt <- phyloregion::PD(random_comm_sparse, phy_alt)
  if ("rpd" %in% metrics) rpd <- pd / pd_alt
  if ("fd" %in% metrics) fd <- phyloregion::PD(random_comm_sparse, trait_tree)
  if ("rfd" %in% metrics) fd_alt <- phyloregion::PD(random_comm_sparse, trait_tree_alt)
  if ("rfd" %in% metrics) rfd <- fd / fd_alt
  if ("pe" %in% metrics) pe <- phyloregion::phylo_endemism(random_comm_sparse, phy, weighted = TRUE)
  if ("rpe" %in% metrics) pe_alt <- phyloregion::phylo_endemism(random_comm_sparse, phy_alt, weighted = TRUE)
  if ("rpe" %in% metrics) rpe <- pe / pe_alt
  
  # Output results
  list(
    pd = pd,
    pd_alt = pd_alt,
    rpd = rpd,
    fd = fd,
    fd_alt = fd_alt,
    rfd = rfd,
    pe = pe,
    pe_alt = pe_alt,
    rpe = rpe
  ) %>%
    # Only keep non-NULL results
    purrr::compact()
  
}

#' Extract standard effect size (and other related statistics) for a single
#' diversity metric given random values and observed values of the metric
#'
#' @param random_vals List of list of vectors. Each list of vectors is a biodiversity metric measured on a
#' random community, in the same order as the rows in the input community.
#' @param obs_vals Observed values of the biodiversity metric
#' @param metric Name of the metric ("mpd", "mntd", "mpd_morph", "mntd_morph", "pd", "pe", or "rpe")
#'
#' @return Tibble
#' 
get_ses <- function (random_vals, obs_vals, metric) {
  
  assert_that(is.string(metric))
  
  assert_that(
    all(metric %in% c("pd", "pd_alt", "rpd", "fd", "fd_alt", "rfd", "pe", "pe_alt", "rpe")),
    msg = "Biodiversity metrics may only be selected from 'pd', 'rpd', 'fd', 'rfd', 'pe', or 'rpe'"
  )
  
  random_vals_trans <- transpose(random_vals)
  
  results <-
    tibble(
      random_values = transpose(random_vals_trans[[metric]]) %>% map(as_vector),
      obs_val = obs_vals
    ) %>%
    transmute(
      obs = obs_val,
      rand_mean = purrr::map_dbl(random_values, ~mean(., na.rm = TRUE)),
      rand_sd = purrr::map_dbl(random_values, ~sd(., na.rm = TRUE)),
      obs_rank = purrr::map2_dbl(.x = obs_val, .y = random_values, ~rank(c(.x, .y), ties.method = "average")[[1]]),
      obs_z = (obs_val - rand_mean) / rand_sd,
      obs_p = obs_rank/(length(random_vals) + 1)
    )
  
  colnames(results) <- paste(metric, colnames(results), sep = "_")
  
  results
}

#' Run standard effect size analysis for a set of biodiversity metrics
#' 
#' The biodiversity metrics analyzed include:
#'   - pd: Phylogenetic diversity (Faith 1992 https://doi.org/10.1016/0006-3207(92)91201-3)
#'   - rpd: Relative phylogenetic diversity
#'   - fd: Functional diversity
#'   - rfd: Relative functional diveristy
#'   - pe: Phylogenetic endemism (Rosauer 2009 https://doi.org/10.1111/j.1365-294x.2009.04311.x)\
#'   - rpe: Relative phylogenetic endemism (Mishler 2014 https://doi.org/10.1038/ncomms5473)
#'   
#' The independent swap method of Gotelli (2000) is used to generate random communities, which randomizes
#' the community matrix while maintaining species occurrence frequency and
#' sample species richness. Each random community is generated by shuffling observed
#' data 10,000 times.
#'   
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names)
#' @param phy Input phylogeny with total branch length scaled to 1
#' @param n_reps Number of random communities to replicate
#' @param metrics Names of metrics to calculate. Must one or more of
#' 'pd', 'rpd', 'fd', 'rfd', 'pe', or 'rpe'
#'
#' @return Tibble. For each of the biodiversity metrics, the observed value (_obs), 
#' mean of the random values (_rand_mean), SD of the random values (_rand_sd), 
#' rank of the observed value vs. the random values (_obs_rank), standard effect size
#' (_obs_z), and p-value (_obs_p) are given. Type of phylogenetic endemism (neo, paleo, or mixed)
#' is given as "endem_type".
#' 
run_ses_analysis <- function(comm_df, phy = NULL, trait_distances = NULL, n_reps, metrics) {
  
  # Make dummy phy_alt and trait_tree_alt in case one of these isn't being analyzed
  phy_alt <- NULL
  trait_tree <- NULL
  trait_tree_alt <- NULL
  
  # If running both phylogenetic and functional metrics,
  # subset to only taxa with represented in all datasets
  if (any(str_detect(metrics, "^pd$|^rpd$|^pe$|^rpe$")) & any(str_detect(metrics, "^fd$|^rfd$"))) {
    
    # Convert distances to matrix for subsetting
    trait_distances <- as.matrix(trait_distances)
    
    taxa_keep <- intersect(phy$tip.label, colnames(comm_df), colnames(trait_distances))
    phy <- ape::keep.tip(phy, taxa_keep)
    comm_df <- comm_df[,colnames(comm_df) %in% taxa_keep]
    trait_distances <- trait_distances[taxa_keep, taxa_keep] %>% as.dist()
    
  }
  
  # If analyzing phylogenetic diversity, 
  # first match tips of tree and column names of community data frame
  if (any(str_detect(metrics, "^pd$|^rpd$|^pe$|^rpe$"))) {
    
    # Use only taxa that are in common between phylogeny and community
    subsetted_data <- picante::match.phylo.comm(phy = phy, comm = comm_df)
    phy <- subsetted_data[["phy"]]
    comm_df <- subsetted_data[["comm"]]
    
    assert_that(
      isTRUE(
        all.equal(
          sort(colnames(comm_df)),
          sort(phy$tip.label)
        )
      ),
      msg = "Tip names don't match between community and phylogeny"
    )
    
    # Make alternative tree with equal branch lengths
    phy_alt <- phy
    phy_alt$edge.length <- rep(length(phy_alt$edge.length), 1)
    # rescale so total phy length is 1
    phy_alt$edge.length <- phy_alt$edge.length / sum(phy_alt$edge.length)
    # rescale original phy so total length is 1
    phy$edge.length <- phy$edge.length / sum(phy$edge.length)
    
  }
  
  # If analyzing function diversity, 
  # use only taxa that are in common between traits and community
  if (any(str_detect(metrics, "^fd$|^rfd$"))) {
    
    trait_tree <- hclust(trait_distances, method = "average") %>%
      ape::as.phylo()
    
    # Use only taxa that are in common between phylogeny and community
    subsetted_data <- picante::match.phylo.comm(phy = trait_tree, comm = comm_df)
    trait_tree <- subsetted_data[["phy"]]
    comm_df <- subsetted_data[["comm"]]
    
    assert_that(
      isTRUE(
        all.equal(
          sort(colnames(comm_df)),
          sort(trait_tree$tip.label)
        )
      ),
      msg = "Tip names don't match between community and traits"
    )
    
    # Make alternative tree with equal branch lengths
    trait_tree_alt <- trait_tree
    trait_tree_alt$edge.length <- rep(length(trait_tree_alt$edge.length), 1)
    # rescale so total tree length is 1
    trait_tree_alt$edge.length <- trait_tree_alt$edge.length / sum(trait_tree_alt$edge.length)
    # rescale original tree so total length is 1
    trait_tree$edge.length <- trait_tree$edge.length / sum(trait_tree$edge.length)
    
  }
  
  # Make sparse community df
  comm_df_sparse <- phyloregion::dense2sparse(comm_df)
  
  # Calculate biodiversity metrics for random communities
  random_vals <-
    purrr::rerun(
      n_reps, 
      # Use 10,000 iterations (swaps) for each null community
      calc_biodiv_random(comm_df, phy, phy_alt, trait_tree, trait_tree_alt, 10000, metrics = metrics)
    )
  
  # Calculate biodiversity metrics for observed community
  # - set up null vectors first
  ses_pd <- NULL
  ses_pd_alt <- NULL
  ses_rpd <- NULL
  ses_fd <- NULL
  ses_fd_alt <- NULL
  ses_rfd <- NULL
  ses_pe <- NULL
  ses_pe_alt <- NULL
  ses_rpe <- NULL
  
  # - calculate selected metrics
  if ("pd" %in% metrics) {
    pd_obs <- phyloregion::PD(comm_df_sparse, phy)
    ses_pd <- get_ses(random_vals, pd_obs, "pd")}
  
  if ("rpd" %in% metrics) {
    pd_alt_obs <- phyloregion::PD(comm_df_sparse, phy_alt)
    ses_pd_alt <- get_ses(random_vals, pd_alt_obs, "pd_alt")}
  
  if ("rpd" %in% metrics) {
    rpd_obs <- pd_obs / pd_alt_obs
    ses_rpd <- get_ses(random_vals, rpd_obs, "rpd")}
  
  if ("fd" %in% metrics) {
    fd_obs <- phyloregion::PD(comm_df_sparse, trait_tree)
    ses_fd <- get_ses(random_vals, fd_obs, "fd")}
  
  if ("rfd" %in% metrics) {
    fd_alt_obs <- phyloregion::PD(comm_df_sparse, trait_tree_alt)
    ses_fd_alt <- get_ses(random_vals, fd_alt_obs, "fd_alt")}
  
  if ("rfd" %in% metrics) {
    rfd_obs <- fd_obs / fd_alt_obs
    ses_rfd <- get_ses(random_vals, rfd_obs, "rfd")}
  
  if ("pe" %in% metrics) {
    pe_obs <- phyloregion::phylo_endemism(comm_df_sparse, phy, weighted = TRUE)
    ses_pe <- get_ses(random_vals, pe_obs, "pe")}
  
  if ("rpe" %in% metrics) {
    pe_alt_obs <- phyloregion::phylo_endemism(comm_df_sparse, phy_alt, weighted = TRUE)
    ses_pe_alt <- get_ses(random_vals, pe_obs, "pe_alt")}
  
  if ("rpe" %in% metrics) {
    rpe_obs <- pe_obs / pe_alt_obs
    ses_rpe <- get_ses(random_vals, rpe_obs, "rpe")
  }
  
  # Combine results
  bind_cols(
    ses_pd,
    ses_pd_alt,
    ses_rpd,
    ses_fd,
    ses_fd_alt,
    ses_rfd,
    ses_pe,
    ses_pe_alt,
    ses_rpe
  ) %>%
    mutate(site = rownames(comm_df)) %>%
    select(site, everything())
  
}

#' Categorize phylogenetic endemism
#' 
#' see:
#' http://biodiverse-analysis-software.blogspot.com/2014/11/canape-categorical-analysis-of-palaeo.html
#'
#' (here, PE_orig = pe_obs_p, PE_alt = pe_alt_obs_p, and RPE = rpe_obs_p)
#' (don't consider 'super endemism' as it doesn't add much meaning)
#'
#' 1)    If either PE_orig or PE_alt are significantly high then we look for palaeo or neo endemism
#'   a)    If RPE is significantly high then we have palaeo-endemism 
#'         (PE_orig is consistently higher than PE_alt across the random realisations)
#'   b)    Else if RPE is significantly low then we have neo-endemism 
#'         (PE_orig is consistently lower than PE_alt across the random realisations)
#'     c)    Else we have mixed age endemism in which case
#'        i)    If both PE_orig and PE_alt are highly significant (p<0.01) then we 
#'              have super endemism (high in both palaeo and neo)
#'        ii)   Else we have mixed (some mixture of palaeo, neo and non endemic)
#' 2)    Else if neither PE_orig or PE_alt are significantly high then we have a non-endemic cell
#'
#' @param df Input data frame. Must have p-values for pe, pe_alt, and rpe.
#'
#' @return Dataframe with areas of endemism categorized.
#' 
categorize_endemism <- function (df) {
  
  df %>% 
  mutate(
    endem_type = case_when(
      (pe_obs_p >= 0.95 | pe_alt_obs_p >= 0.95) & rpe_obs_p >= 0.975 ~ "paleo",
      (pe_obs_p >= 0.95 | pe_alt_obs_p >= 0.95) & rpe_obs_p <= 0.025 ~ "neo",
      pe_obs_p >= 0.95 | pe_alt_obs_p >= 0.95 ~ "mixed",
      TRUE ~ NA_character_
    )
  )
  
}

# Phyloregion ----

#' Find optimal K for clustering by taxonomy
#'
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names)
#' @param k_max Maximum number of clusters to assign
#'
#' @return a list containing the following as returned from the GMD package (Zhao et al. 2011):
#' - k: optimal number of clusters (bioregions)
#' - totbss: total between-cluster sum-of-square
#' - tss: total sum of squares of the data
#' - ev: explained variance given k
#' 
find_k_taxonomy <- function (comm_df, k_max = 20) {
  
  # Convert input dataframe to sparse
  comm_sparse <- phyloregion::dense2sparse(comm_df)
  
  # Calculate taxonomic beta diversity
  bc <- phyloregion::beta_diss(comm_sparse, index.family = "sorensen")
  
  # Select the best clustering algorithm
  method_scores <- phyloregion::select_linkage(bc[[1]])
  
  # The values used for `method` in phyloregion::phyloregion() are
  # slightly different from the output of phyloregion::select_linkage(),
  # so adjust appropriately
  
  method_names <- c(
    ward.D = "ward.D", ward.D2 = "ward.D2", Single = "single",
    Complete = "complete", UPGMA = "average", WPGMA = "mcquitty",
    WPGMC = "median",  UPGMC = "centroid")
  
  method_select <- method_names[names(method_scores[1])]
  
  # Find the optimal K value
  phyloregion::optimal_phyloregion(bc[[1]], method = method_select, k = k_max)
  
}

#' Find optimal K for clustering by phylogeny
#'
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names)
#' @param k_max Maximum number of clusters to assign
#' @param phy Phylogenetic tree
#'
#' @return a list containing the following as returned from the GMD package (Zhao et al. 2011):
#' - k: optimal number of clusters (bioregions)
#' - totbss: total between-cluster sum-of-square
#' - tss: total sum of squares of the data
#' - ev: explained variance given k
#' 
find_k_phylogeny <- function (comm_df, k_max = 20, phy) {
  
  # Use only taxa that are in common between phylogeny and community
  subsetted_data <- picante::match.phylo.comm(phy = phy, comm = comm_df)
  phy <- subsetted_data[["phy"]]
  comm_df <- subsetted_data[["comm"]]
  
  # Convert input dataframe to sparse
  comm_sparse <- phyloregion::dense2sparse(comm_df)
  
  # Calculate taxonomic beta diversity
  bc <- phyloregion::phylobeta(comm_sparse, phy, index.family = "sorensen")
  
  # Select the best clustering algorithm
  method_scores <- phyloregion::select_linkage(bc[[1]])
  
  # The values used for `method` in phyloregion::phyloregion() are
  # slightly different from the output of phyloregion::select_linkage(),
  # so adjust appropriately
  
  method_names <- c(
    ward.D = "ward.D", ward.D2 = "ward.D2", Single = "single",
    Complete = "complete", UPGMA = "average", WPGMA = "mcquitty",
    WPGMC = "median",  UPGMC = "centroid")
  
  method_select <- method_names[names(method_scores[1])]
  
  # Find the optimal K value
  phyloregion::optimal_phyloregion(bc[[1]], method = method_select, k = k_max)
  
}

#' Cluster regions by taxonomy
#'
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names)
#' @param k Number of clusters to assign
#'
#' @return Dataframe
#' 
cluster_taxonomic_regions <- function (comm_df, k) {
  
  # Convert input dataframe to sparse
  comm_sparse <- phyloregion::dense2sparse(comm_df)
  
  # Calculate taxonomic beta diversity
  bc <- phyloregion::beta_diss(comm_sparse, index.family = "sorensen")
  
  # Select the best clustering algorithm
  method_scores <- phyloregion::select_linkage(bc[[1]])
  
  # The values used for `method` in phyloregion::phyloregion() are
  # slightly different from the output of phyloregion::select_linkage(),
  # so adjust appropriately
  
  method_names <- c(
    ward.D = "ward.D", ward.D2 = "ward.D2", Single = "single",
    Complete = "complete", UPGMA = "average", WPGMA = "mcquitty",
    WPGMC = "median",  UPGMC = "centroid")
  
  method_select <- method_names[names(method_scores[1])]
  
  # Assign taxonomic regions based on clustering
  regions <- phyloregion::phyloregion(bc[[1]], k = k, method = method_select)
  
  regions[["membership"]] %>%
    as_tibble()
}

#' Cluster regions by phylogenetic distance
#'
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names)
#' @param k Number of clusters to assign
#' @param phy Phylogenetic tree
#'
#' @return Dataframe
#' 
cluster_phylo_regions <- function (comm_df, phy, k) {
  
  # Use only taxa that are in common between phylogeny and community
  subsetted_data <- picante::match.phylo.comm(phy = phy, comm = comm_df)
  phy <- subsetted_data[["phy"]]
  comm_df <- subsetted_data[["comm"]]
  
  # Convert input dataframe to sparse
  comm_sparse <- phyloregion::dense2sparse(comm_df)
  
  # Calculate phylogenetic beta diversity
  bc <- phyloregion::phylobeta(comm_sparse, phy, index.family = "sorensen")
  
  # Select the best clustering algorithm
  method_scores <- phyloregion::select_linkage(bc[[1]])
  
  # The values used for `method` in phyloregion::phyloregion() are
  # slightly different from the output of phyloregion::select_linkage(),
  # so adjust appropriately
  
  method_names <- c(
    ward.D = "ward.D", ward.D2 = "ward.D2", Single = "single",
    Complete = "complete", UPGMA = "average", WPGMA = "mcquitty",
    WPGMC = "median",  UPGMC = "centroid")
  
  method_select <- method_names[names(method_scores[1])]
  
  # Assign regions based on clustering
  regions <- phyloregion::phyloregion(bc[[1]], k = k, method = method_select)
  
  regions[["membership"]] %>%
    as_tibble()
}



# Ecostructure ----

#' Count abundance
#'
#' Counts number of times an integer in a vector
#' occurrs out of a given set of integers.
#'
#' @param cells Numeric vector
#' @param all_cells Numeric vector
#'
#' @return Named numeric vector
#'
#' @examples
#' count_abun(c(10, 10, 12), c(10, 12, 14))
count_abun <- function(cells, all_cells) {
  cells <- as.factor(as.factor(cells))
  abun_occur <- tabulate(cells)
  names(abun_occur) <- levels(cells)
  abun <- numeric(length(all_cells))
  names(abun) <- all_cells
  abun[names(abun_occur)] <- abun_occur
  abun
}

#' Make community matrix from species' occurrences
#'
#' @param species_coods Dataframe of species occurrences.
#' No missing values allowed.
#' @param resol Degree of spatial resolution used to construct the presence-absence grid.
#' @param species Name of column to use for species
#' @param lat Name of column to use for latitude
#' @param long Name of column to use for longitude
#' @param crs Character or object of class CRS. PROJ.4 type description of a
#' Coordinate Reference System (map projection) used to construct the
#' presence-absence grid.
#'
#' @return List of two items, "comm_dat" is the community dataframe in
#' sparse format; poly_shp is the distribution map as a simple features dataframe
#'
#' @examples
#' # Make test data set of 10 species with 100 occurrences total.
#' test_data <- data.frame(taxon = letters[1:10],
#'   longitude = runif(100, -2, 2),
#'   decimallatitude = runif(100, -2, 2),
#'   latitude = FALSE)
#' # Presence/absence
#' comm_from_points(test_data)
#' # Abundance
#' comm_from_points(test_data, abun = TRUE)
comm_from_points <- function(species_coods,
                             resol = 1,
                             species = "taxon",
                             lon = "longitude",
                             lat = "latitude",
                             crs = sp::CRS("+proj=longlat +datum=WGS84")) {
  
  species_coods <- as.data.frame(species_coods)
  species_coods <- species_coods[, c(species, lon, lat)]
  names(species_coods) <- c("species", "longitude", "latitude")
  
  checkr::check_data(
    species_coods,
    values = list(
      species = "a",
      longitude = 1,
      latitude = 1
    )
  )
  
  assertr::verify(species_coods, longitude <= 180, success_logical)
  
  assertr::verify(species_coods, latitude <= 90, success_logical)
  
  assertthat::assert_that(!is.factor(species_coods$species))
  
  # Create empty raster
  r <- raster::raster(
    resolution = resol,
    xmn = -180,
    xmx = 180,
    ymn = -90,
    ymx = 90,
    crs = crs)
  
  ### Extract cell IDs from points by species
  # Coordinates must be long, lat for raster
  # (which assumes x, y position in that order)
  
  # Extract cell IDs from points by species
  species_coods_nested <- tidyr::nest(species_coods, data = c(longitude, latitude))
  
  # raster::cellFromXY needs data to be data.frame, not tibble
  species_coods_nested <- dplyr::mutate(species_coods_nested, data = purrr::map(data, as.data.frame))
  
  cells_occur <- purrr::map(
    species_coods_nested$data,
    ~ raster::cellFromXY(xy = ., object = r)
  )
  
  names(cells_occur) <- species_coods_nested$species
  
  # Get vector of non-empty cells
  non_empty_cells <- sort(unique(unlist(cells_occur)))
  
  # Subset raster to non-empty cells
  r_sub <- raster::rasterFromCells(r, non_empty_cells)
  
  # Get its extent
  e <- raster::extent(r_sub)
  
  # Convert to spatial polygons dataframe
  p <- as(e, "SpatialPolygons")
  
  # Make polygons grid at same extent as sample points
  m <- sp::SpatialPolygonsDataFrame(p, data.frame(sp = "x"))
  m <- fishnet(mask = m, res = resol)
  
  # Overlay points onto grid
  dat <- as.data.frame(species_coods)
  dat <- dat[complete.cases(dat), ]
  sp::coordinates(dat) <- ~longitude + latitude
  
  sp::proj4string(dat) <- sp::proj4string(m)
  x <- sp::over(dat, m)
  y <- cbind(as.data.frame(dat), x)
  y <- y[complete.cases(y), ]
  Y <- phyloregion::long2sparse(y)
  tmp <- data.frame(grids = row.names(Y), abundance = rowSums(as.matrix(Y)), 
                    richness = rowSums(as.matrix(Y > 0)))
  z <- sp::merge(m, tmp, by = "grids")
  z <- z[!is.na(z@data$richness), ] %>% sf::st_as_sf()
  
  return(list(comm_dat = Y, poly_shp = z))
  
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
    gather(site, abundance, -species) %>%
    pull(site) %>%
    unique
  
  # Make tibble of all longitudes and latitudes by
  # site name for renaming
  # rows in matrix
  long_lats <-
    all_cells %>%
    transmute(
      site = as.character(site), 
      long_lat = paste(longitude, latitude, sep = "_")) %>%
    unique %>%
    # Only use those actually in the community data so we
    # don't end up with a bunch of NA values after joining.
    filter(site %in% pterido_sites)
  
  comm_pteridos %>%
    gather(site, abundance, -species) %>%
    spread(species, abundance) %>%
    # Check for grid cells in comm data but missing from all_cells data
    verify(site %in% long_lats$site) %>%
    # Check that grid cells in comm data are unique
    assert(is_uniq, site) %>%
    left_join(long_lats, by = "site") %>%
    select(-site) %>%
    column_to_rownames("long_lat") %>%
    as.matrix
  
}

# Slightly tweaked version of ecos_plot_pie that takes sf object as background map
# instead of reading in shape file.
ecos_plot_pie2 <- function (
  omega = NULL, coords = NULL, bgmap_path = NULL, 
  no_map = FALSE, blank_bg = NULL,
  adjust = FALSE,
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
  if (is.null(bgmap_path) & no_map == FALSE) {
    message("reading background map shapefile from inst/extdata/ne_110m_coastline \n            folder")
    GlobalCoast <- sf::st_read(system.file("extdata", "ne_110m_land",
                                           "ne_110m_land.shp", package = "ecostructure"), quiet = T)
  }
  if (!is.null(bgmap_path) & no_map == FALSE) {
    GlobalCoast <- bgmap_path
    glob <- c(xmin = long_lim[1], xmax = long_lim[2], ymin = lat_lim[1],
              ymax = lat_lim[2])
    glob <- sf::st_bbox(glob)
    glob <- structure(glob, crs = sf::st_crs(GlobalCoast))
    GlobalCoast <- suppressWarnings(suppressMessages(sf::st_intersection(GlobalCoast,
                                                                         sf::st_as_sfc(glob))))
  }
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
  if (isTRUE(no_map)) {
    if(is.null(blank_bg)) blank_bg <- "lightblue"
    mapplots::basemap(xlim = long_lim, ylim = lat_lim, bg = blank_bg)
  } else {
    plot(sf::st_geometry(GlobalCoast), axes = T, main = "", reset = F,
         xaxs = "i", yaxs = "i", lwd = coastline_lwd) }
  par(lwd = 0.01)
  invisible(lapply(1:dim(omega)[1], function(r) do.call(mapplots::add.pie,
                                                        append(list(z = as.integer(100 * omega[r, ]), x = coords[r,
                                                                                                                 1], y = coords[r, 2], labels = c("", "", ""), radius = radius,
                                                                    color), pie_control))))
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
  
  # Optionally drop columns and rows that are all zeros
  if (isTRUE(drop_zero)) {
    map_data <- map_data[,colSums(map_data) != 0]
    map_data <- map_data[rowSums(map_data) != 0,]
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

map_theme <- function() {
  theme(
    panel.grid.major = element_line(color = "grey95", size = 0.1),
    panel.background = element_rect(fill = "grey60"),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)
  )
}

map_theme_light <- function () {
  theme(
    panel.grid.major = element_line(color = "white", size = 0.1),
    panel.background = element_rect(fill = "grey90"),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.7)
  )
}

map_theme_pp <- function() {
  theme(
    panel.grid.major = element_line(color = "white", size = 0.1),
    panel.background = element_rect(fill = "grey70"),
    axis.ticks = element_blank()
  )
}

map_theme_pp_light <- function() {
  theme(
    panel.grid.major = element_line(color = "white", size = 0.1),
    panel.background = element_rect(fill = "grey90"),
    axis.ticks = element_blank()
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
#' @param occ_data Occurrence data, with one row per
#' grid cell per taxon, including hybrids.
#' @param div_metric Selected diversity metric to plot.
#' Must one of the column names of div_dat.
#' @param metric_title Character: title to use for legend for
#' diversity metric.
#'
#' @return ggplot object
make_diversity_map <- function (div_data, occ_data, div_metric, metric_title, label) {
  
  div_metric <- sym(div_metric)
  
  ggplot(data = div_data, aes(x = longitude, y = latitude)) +
    geom_tile(aes(fill = !!div_metric)) + 
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor, 
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor, 
               pull(occ_data, latitude) %>% max %>% ceiling)
    ) +
    labs(
      fill = metric_title,
      subtitle = label
    ) +
    # jntools::blank_x_theme() +
    # jntools::blank_y_theme() +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_line(size = 0.1, color = "dark grey"),
      panel.grid.minor = element_blank(),
      plot.subtitle = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 20/.pt),
      legend.text = element_text(size = 16/.pt),
      legend.justification=c(0,1), legend.position=c(0,1)
    ) +
    scale_y_continuous(labels = scales::label_number(suffix = "\u00b0N", accuracy = 1)) +
    scale_x_continuous(labels = scales::label_number(suffix = "\u00b0E", accuracy = 1))
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
#' @param title_color Color to use for background of title
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
                                 title_color = NULL,
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
  
  if(!is.null(title_color)) {
    plot <- plot + theme(
      plot.title = element_textbox_simple(
        fill = title_color,
        size = 13,
        lineheight = 1,
        padding = margin(5.5, 5.5, 5.5, 5.5),
        margin = margin(0, 0, 5.5, 0)
      )
    )
  }
  
  plot
  
}

# Phylogeneic analysis ----

#' Read in Japan rbcL alignment from the downloaded Dryad zip file
#'
#' @param zip_folder Path to downloaded Dryad zip file
#'
#' @return List of class DNAbin; DNA alignment
#' 
read_ja_rbcL_from_zip <- function (zip_folder) {
  
  temp_dir <- tempdir()
  
  unzip(zip_folder, exdir = temp_dir)
  
  unzip(fs::path(temp_dir, "japan_pterido_rbcl_cipres.zip"), exdir = temp_dir)
  
  ape::read.nexus.data(fs::path(temp_dir, "japan_pterido_rbcl_cipres/infile.nex")) %>%
    as.DNAbin
  
}

#' Read in Green List from the downloaded Dryad zip file
#'
#' @param zip_folder Path to downloaded Dryad zip file
#'
#' @return Dataframe
#' 
read_greenlist_from_zip <- function (zip_folder) {
  
  temp_dir <- tempdir()
  
  unzip(zip_folder, exdir = temp_dir)
  
  readxl::read_excel(fs::path(temp_dir, "FernGreenListV1.01E.xls"))
  
}

#' Combine Japanese rbcL sequences with globally sampled genes
#'
#' @param broad_alignment_list List of DNA alignments, including one called "rbcL"
#' @param japan_rbcL List of class DNAbin; rbcL sequences of ferns from Japan
#'
combine_ja_rbcL_with_global <- function (broad_alignment_list, japan_rbcL) {
  
  # Remove Japan names from broad alignment list
  reduced <- map(broad_alignment_list, ~magrittr::extract(., !(rownames(.) %in% rownames(japan_rbcL)),))
  
  # Add Japan rbcL sequences back, align with mafft
  reduced[["rbcL"]] <- c(as.list(reduced[["rbcL"]]), as.list(japan_rbcL)) %>%
    ips::mafft(exec = "/usr/bin/mafft") %>%
    # Trim any column that is >90% gaps
    ips::deleteGaps(gap.max = nrow(.)*0.9)
  
  # Concatenate genes
  concatenate_genes(reduced)
}

#' Concatenate a list of aligned genes
#'
#' @param dna_list List of matrices of class DNAbin
#'
#' @return Matrix of class DNAbin
#'
#' @examples
#' data(woodmouse)
#' gene_1 <- woodmouse[,1:100]
#' gene_2 <- woodmouse[,101:200]
#' woodmouse_genes <- list(gene_1, gene_2)
#' concatenate_genes(woodmouse_genes)
concatenate_genes <- function (dna_list) {
  require(apex)
  
  assertthat::assert_that(is.list(dna_list))
  assertthat::assert_that(all(lapply(dna_list, class) == "DNAbin"), 
                          msg = "All elements of dna_list must be of class DNAbin")
  assertthat::assert_that(all(sapply(dna_list, is.matrix)), 
                          msg = "All elements of dna_list must be matrices")
  
  # Check that there are no duplicate sequence names (species) within a gene
  map_df(dna_list, ~rownames(.) %>% tibble(species = .), .id = "gene") %>%
    assert_rows(col_concat, is_uniq, species, gene, error_fun = assertr::error_stop)
  
  dna_multi <- new("multidna", dna_list) 
  apex::concatenate(dna_multi)
}


#' Format Japan rbcL sequence names
#'
#' @param alignment DNA alignment of rbcL gene for Japanese ferns with names coded
#' as numeric taxon IDs
#' @param taxon_id_map Dataframe of standard taxonomy with
#' columns `taxon_id` and `taxon`
#'
#' @return Renamed DNA alignment with species as names instead of codes
#'
rename_alignment <- function (alignment, taxon_id_map) {
  
  new_names <-
    tibble(name = rownames(alignment)) %>%
    mutate(taxon_id = str_split(name, "_") %>% map_chr(2)) %>%
    left_join(taxon_id_map, by = "taxon_id") %>%
    assert(not_na, taxon) %>%
    assert(is_uniq, taxon) %>%
    pull(taxon)
  
  rownames(alignment) <- new_names
  
  alignment
  
}

#' Concatenate a list of aligned genes
#'
#' @param dna_list List of matrices of class DNAbin
#'
#' @return Matrix of class DNAbin
#'
#' @examples
#' data(woodmouse)
#' gene_1 <- woodmouse[,1:100]
#' gene_2 <- woodmouse[,101:200]
#' woodmouse_genes <- list(gene_1, gene_2)
#' concatenate_genes(woodmouse_genes)
concatenate_genes <- function (dna_list) {
  require(apex)
  
  assertthat::assert_that(is.list(dna_list))
  assertthat::assert_that(all(lapply(dna_list, class) == "DNAbin"), 
                          msg = "All elements of dna_list must be of class DNAbin")
  assertthat::assert_that(all(sapply(dna_list, is.matrix)), 
                          msg = "All elements of dna_list must be matrices")
  
  # Check that there are no duplicate sequence names (species) within a gene
  map_df(dna_list, ~rownames(.) %>% tibble(species = .), .id = "gene") %>%
    assert_rows(col_concat, is_uniq, species, gene, error_fun = assertr::error_stop)
  
  dna_multi <- new("multidna", dna_list) 
  apex::concatenate(dna_multi)
}

# Dating with treePL ----

#' Read in calibration and configure dates for treepl
#'
#' @param date_file_path Path to CSV file with treepl dates.
#' Must include at least the following columns:
#' - clade: name of clade
#' - taxon_1: representative taxon 1
#' - taxon_2: representative taxon 2. The MRCA of the two taxa defines the clade
#' - age: Age to assign to clade (in millions of years)
#' - age_type: 'min', 'max' or 'fixed'
#' 
#' @return Tibble with columns for use in treepl config file.
#'
load_calibration_dates <- function(date_file_path) {
  read_csv(date_file_path) %>%
    janitor::clean_names() %>%
    select(clade, age, age_type, taxon_1, taxon_2) %>%
    assert(is_uniq, clade) %>%
    assert(not_na, clade, age, age_type, taxon_1, taxon_2) %>%
    mutate(mrca = glue("mrca = {clade} {taxon_1} {taxon_2}")) %>%
    mutate(
      min_dates = case_when(
        age_type %in% c("min", "fixed") ~ glue("min = {clade} {age}"),
        TRUE ~ NA_character_),
      max_dates = case_when(
        age_type %in% c("max", "fixed") ~ glue("max = {clade} {age}"),
        TRUE ~ NA_character_))
} 

#' Do an initial treepl run to determine optimal
#' smoothing parameters with random cross-validation.
#' 
#' Requires treepl to be installed and on PATH.
#' 
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phy"; phylogeny
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param write_tree Logical; should the phylogeny be written to working
#' directory before running treePL?
#' @param cvstart Start value for cross-validation
#' @param cvstop Stop value for cross-validation
#' @param cvsimaniter Start the number of cross validation simulated annealing 
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing 
#' iterations, default = 5000
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl_cv <- function (
  phy, alignment, calibration_dates, 
  write_tree = TRUE,
  cvstart = "1000", cvstop = "0.1",
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd, echo) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>% unique
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue::glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to working directory
  phy_name <- deparse(substitute(phy))
  phy_path <- fs::path_ext_set(phy_name, "tre")
  if(isTRUE(write_tree)) {ape::write.tree(phy, fs::path(wd, phy_path))}
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2]
  
  outfile_path <- fs::path_ext_set(paste0(phy_name, "_cv"), "out")
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min_dates)) %>% pull(min_dates),
    calibration_dates %>% filter(!is.na(max_dates)) %>% pull(max_dates),
    glue("cvstart = {cvstart}"),
    glue("cvstop = {cvstop}"),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("cvoutfile = {outfile_path}"),
    glue("seed = {seed}"),
    glue("nthreads = {nthreads}"),
    "randomcv"
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  readr::write_lines(treepl_config, fs::path(wd, "treepl_cv_config"))
  
  # Run treePL
  processx::run("treePL", "treepl_cv_config", wd = wd, echo = echo)
  
  # Return cross-validation results
  read_lines(fs::path(wd, outfile_path))
  
}

#' Do an initial treepl run to determine optimal
#' smoothing parameters with random cross-validation.
#' 
#' Requires treepl to be installed and on PATH.
#' 
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phy"; phylogeny
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param write_tree Logical; should the phylogeny be written to working
#' directory before running treePL? Can be FALSE if it is already there from
#' run_treepl_initial().
#' @param cv_results Output of run_treepl_cv() so the best smoothing
#' parameter can be selected.
#' @param cvsimaniter Start the number of cross validation simulated annealing 
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing 
#' iterations, default = 5000
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl_prime <- function (
  phy, alignment, calibration_dates, 
  cv_results,
  write_tree = FALSE,
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd, echo) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>% unique
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue::glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to working directory
  phy_name <- deparse(substitute(phy))
  phy_path <- fs::path_ext_set(phy_name, "tre")
  if(isTRUE(write_tree)) {ape::write.tree(phy, fs::path(wd, phy_path))}
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2]
  
  # Get best smoothing parameter
  # Select the optimum smoothing value (smallest chisq) from cross-validation
  # The raw output looks like this:
  # chisq: (1000) 6.7037e+30
  # chisq: (100) 3673.45
  # etc.
  best_smooth <-
    tibble(cv_result = cv_results) %>%
    mutate(smooth = str_match(cv_result, "\\((.*)\\)") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    mutate(chisq = str_match(cv_result, "\\) (.*)$") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    arrange(chisq) %>%
    slice(1) %>%
    pull(smooth)
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min_dates)) %>% pull(min_dates),
    calibration_dates %>% filter(!is.na(max_dates)) %>% pull(max_dates),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("seed = {seed}"),
    glue("smooth = {best_smooth}"),
    glue("nthreads = {nthreads}"),
    "prime"
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  readr::write_lines(treepl_config, fs::path(wd, "treepl_prime_config"))
  
  # Run treePL
  results <- processx::run("treePL", "treepl_prime_config", wd = wd, echo = echo)
  
  # Return stdout
  read_lines(results$stdout)
  
}

#' Run treePL
#' 
#' Requires treepl to be installed and on PATH.
#' 
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phy"; phylogeny
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param priming_results Results of running treePL with `prime` option to
#' determine optional parameters. Output of run_treepl_prime().
#' @param cv_results Results of running treePL with cross-validation to
#' determine optimal rate-smoothing parameter. Output of run_treepl_cv().
#' @param write_tree Logical; should the phylogeny be written to working
#' directory before running treePL? Can be FALSE if it is already there from
#' run_treepl_initial().
#' @param cvsimaniter Start the number of cross validation simulated annealing 
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing 
#' iterations, default = 5000
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl <- function (
  phy, alignment, calibration_dates, 
  priming_results,
  cv_results,
  write_tree = FALSE,
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd, echo) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>% unique
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue::glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to working directory
  phy_name <- deparse(substitute(phy))
  phy_path <- fs::path_ext_set(phy_name, "tre")
  if(isTRUE(write_tree)) {ape::write.tree(phy, fs::path(wd, phy_path))}
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2]
  
  # Get best smoothing parameter
  # Select the optimum smoothing value (smallest chisq) from cross-validation
  # The raw output looks like this:
  # chisq: (1000) 6.7037e+30
  # chisq: (100) 3673.45
  # etc.
  best_smooth <-
    tibble(cv_result = cv_results) %>%
    mutate(smooth = str_match(cv_result, "\\((.*)\\)") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    mutate(chisq = str_match(cv_result, "\\) (.*)$") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    arrange(chisq) %>%
    slice(1) %>%
    pull(smooth)
  
  # Set name of output file
  outfile_path <- fs::path_ext_set(paste0(phy_name, "_dated"), "tre")
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    glue("smooth = {best_smooth}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min_dates)) %>% pull(min_dates),
    calibration_dates %>% filter(!is.na(max_dates)) %>% pull(max_dates),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("seed = {seed}"),
    glue("nthreads = {nthreads}"),
    glue("outfile = {outfile_path}"),
    # Include specifications from priming
    priming_results %>% extract(., str_detect(., "opt =")),
    priming_results %>% extract(., str_detect(., "optad =")),
    priming_results %>% extract(., str_detect(., "optcvad ="))
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  readr::write_lines(treepl_config, fs::path(wd, "treepl_config"))
  
  # Run treePL
  processx::run("treePL", "treepl_config", wd = wd, echo = echo)
  
  # Read in tree
  ape::read.tree(fs::path(wd, outfile_path))
  
}

# Manuscript rendering ----

#' Generate a path to save a results file
#' 
#' Only works for figures or tables cited in text, and outputs
#' files to the "results" folder
#'
#' @param result_num Number of result, e.g. "Fig. 2"
#' @param extension Extension to use for file
#'
#' @return String.
#' 
result_file <- function (result_num, extension) {
  
  fs::path(
    here::here("results"),
    result_num %>%
      str_remove_all("\\.") %>% 
      str_replace_all(" ", "_")
  ) %>%
    fs::path_ext_set(extension)
}

