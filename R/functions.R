# Data processing ----

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
#' @return Unzipped data files:
#' - rbcL_clean_sporos.fasta: rbcL sequences of sporophytes from Moorea
#'
unzip_ebihara_2019 <- function (dryad_zip_file, exdir) {
  
  # Unzip only the needed files
  unzip(dryad_zip_file, "FernGreenListV1.01E.xls", exdir = exdir)
  unzip(dryad_zip_file, "ESM1.csv", exdir = exdir)
  unzip(dryad_zip_file, "ESM2.csv", exdir = exdir)
  unzip(dryad_zip_file, "japan_pterido_rbcl_cipres.zip", exdir = exdir)
  unzip(dryad_zip_file, "2_grid_cells_all.csv", exdir = exdir)
  unzip(dryad_zip_file, "ppgi_taxonomy.csv", exdir = exdir)
  
  c(
    fs::path(exdir, "FernGreenListV1.01E.xls"),
    fs::path(exdir, "ESM1.csv"),
    fs::path(exdir, "ESM2.csv"),
    fs::path(exdir, "japan_pterido_rbcl_cipres.zip"),
    fs::path(exdir, "2_grid_cells_all.csv"),
    fs::path(exdir, "ppgi_taxonomy.csv")
  )
  
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
           endemic = Endemism, conservation_status = RL2012,
           hybrid = Hybrid) %>%
    # Simplify taxon names, replace space with underscore
    taxastand::add_parsed_names(scientific_name, taxon) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
    mutate(taxon_id = as.character(taxon_id)) %>%
    # Change hybrid to TRUE/FALSE
    mutate(hybrid = ifelse(is.na(hybrid), FALSE, TRUE)) %>%
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
  
  # Return results as a tibble including the resolution so this
  # can be looped by tar_make() and the results can be selected
  # by resolution
  return(tibble(resol = resol, comm_dat = list(Y), poly_shp = list(z)))
  
}

#' Extract shape dataframe from comm_scaled_list
#'
#' @param data Output of comm_from_points()
#' @param resol_select Target resolution size to use
#'
#' @return Dataframe (sf object)
#' 
shape_from_comm_scaled_list <- function (comm_scaled_list, resol_select) {
  
  comm_scaled_list %>% filter(resol == resol_select) %>% pull(poly_shp) %>% magrittr::extract2(1) %>% sf::st_as_sf()
  
}

#' Extract community dataframe from comm_scaled_list
#'
#' @param data Output of comm_from_points()
#' @param resol_select Target resolution size to use
#'
#' @return Dataframe (sf object)
#' 
comm_from_comm_scaled_list <- function (comm_scaled_list, resol_select) {
  
  comm_scaled_list %>% filter(resol == resol_select) %>% pull(comm_dat) %>% magrittr::extract2(1) %>%
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
#' @param taxon_keep_list Vector of taxon names to keep (others will be dropped)
#'
#' @return Tibble
#'
format_traits <- function(path_to_lucid_traits, taxon_id_map, taxon_keep_list) {
  
  # Read in raw trait data for pteridophytes of Japan.
  # These were originally formatted for lucid dichotomous key software.
  # So they are mostly quantitative traits that have been converted to binary format,
  # or numeric traits. There are a lot of traits. One row per taxon.
  traits_raw <- read_excel(path_to_lucid_traits, skip = 1) %>%
    clean_names() %>% 
    select(-x1, -x222) %>%
    rename(taxon = x2) %>%
    mutate(taxon = str_replace_all(taxon, ":", "_")) %>%
    # Check for NA values
    assert(not_na, everything())
  
  # Separate out into numeric and binary traits
  # (numeric container "number" in name, assume binary otherwise)
  traits_numeric <- select(traits_raw, taxon, contains("number"))
  traits_binary <-  select(traits_raw, -contains("number"))
  
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
    # - indusium: treat as present if either true or false indusium is present
    mutate(
      indusium_present = case_when(
        leaf_sorus_indusium_presence_absence_present == 1 ~ 1,
        leaf_sorus_indusium_presence_absence_present_caducous == 1 ~ 1,
        leaf_sorus_false_indusium == 1 ~ 1,
        TRUE ~ 0
      ))
  
  # Repair `leaf_lamina_texture_thick_herbaceous`, which is for some reason split across two columns
  traits_binary <-
  traits_binary %>%
    rowwise() %>%
    # Combine into a single column
    mutate(leaf_lamina_texture_thick_herbaceous = sum(leaf_lamina_texture_thick_herbaceous_53, leaf_lamina_texture_thick_herbaceous_55, na.rm = TRUE)) %>%
    ungroup() %>%
    # Make sure that worked properly
    assert(in_set(c(0,1)), leaf_lamina_texture_thick_herbaceous) %>%
    select(-leaf_lamina_texture_thick_herbaceous_53, -leaf_lamina_texture_thick_herbaceous_55)

  # Select only putatively functional traits: frond shape, frond texture, margin shape, presence or absence of indusium
  
  # # qualitative traits include:
  # (* indicates trait to use; others discarded)
  # *leaf_lamina_shape_sterile_frond 
  # *leaf_lamina_texture
  # leaf_lamina_color
  # leaf_lamina_vennation
  # leaf_lamina_pseudo_veinlet
  # leaf_lamina_terminal_pinna_distinct
  # leaf_lamina_lateral_pinna_shape_fertile
  # leaf_lamina_lateral_pinna_shape_sterile
  # leaf_lamina_lateral_pinna_stalk_fertile
  # *leaf_lamina_margin
  # leaf_sorus_shape
  # leaf_sorus_indusium_shape
  # leaf_sorus_false_indusium_present
  # leaf_lamina_rhachis_adaxial_side_grooved_present
  # leaf_lamina_terminal_pinna_present
  # leaf_sorus_indusium_present
  
  traits_binary <-
  traits_binary %>% 
    select(
      taxon,
      indusium_present,
      contains("leaf_lamina_shape_sterile_frond"),
      contains("leaf_lamina_texture"),
      contains("leaf_lamina_margin")
    ) %>%
    rename_with(~str_replace_all(.x, "leaf_lamina_shape_sterile_frond", "shape")) %>%
    rename_with(~str_replace_all(.x, "leaf_lamina_texture", "texture")) %>%
    rename_with(~str_replace_all(.x, "leaf_lamina_margin", "margin"))
  
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
  
  # Select only putatively functional traits: frond length + width, stipe length, number of pinna pairs
  # (`leaf_lamina_index_length_width_frond` is the ratio of length to width, but not using)
  traits_numeric_combined <-
    traits_numeric_combined %>%
    select(
      taxon,
      frond_length = leaf_lamina_frond_length_cm,
      frond_width = leaf_lamina_frond_width_cm,
      stipe_length = leaf_lamina_index_lamina_length_stipe_length_frond,
      number_pinna_pairs = leaf_lamina_lateral_pinna_of_pairs_frond_number
    )
  
  # Combine all numeric and categorical traits
  traits_combined <- left_join(traits_numeric_combined, traits_binary, by = "taxon")
  
  # Fix some taxon names (synonyms)
  traits_combined <-
    traits_combined %>%
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
    filter(traits_combined, !(taxon %in% taxon_id_map$taxon)) %>%
    pull(taxon)
  
  # Drop any missing names
  traits_combined <- filter(traits_combined, taxon %in% taxon_id_map$taxon)
  
  msg <- assertthat::validate_that(
    length(missing_taxon_id) == 0,
    msg = glue::glue("The following taxa in the trait data could not be verified in the taxa list and have been dropped: {paste(missing_taxon_id, collapse = ', ')}")
  )
  
  if(is.character(msg)) message(msg)
  
  # Keep only taxa in taxon keep list
  traits_combined <- filter(traits_combined, taxon %in% taxon_keep_list)
  
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
  traits_corr_test <- ggplot2::remove_missing(traits_combined)
  
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
  
  if(length(traits_correlated_to_drop) > 0) message(glue::glue("The following traits have correlation coefficient > 0.6 and will be dropped: {traits_correlated_to_drop}"))
  
  ### Drop zero and low-variance traits from full dataframe ###
  
  # Vector of any static trait (only a single trait state)
  traits_static <- 
    traits_combined %>%
    select(all_of(colnames(traits_binary))) %>% # consider binary traits only
    get_static_traits
  
  # Next make a vector of any trait with less than two occurrences of a given trait state
  traits_low_var <-
    traits_combined %>%
    select(all_of(colnames(traits_binary))) %>% # consider binary traits only
    get_low_var_traits
  
  # Drop correlated and non-varying traits from full dataframe
  traits_to_drop <- c(traits_static, traits_low_var, traits_correlated_to_drop) %>% unique()
  traits_combined <- select(traits_combined, -any_of(traits_to_drop))
  
  # Verify that observed correlations in final data are less than 0.6
  traits_combined %>%
    select(-taxon) %>%
    corrr::correlate() %>%
    pivot_longer(-rowname) %>%
    assert(within_bounds(-0.6, 0.6), value, success_fun = success_logical)
  
  traits_combined
  
}

#' Summarize traits for supplemental info
#'
#' @param traits Formatted trait data
#'
#' @return Tibble
#'
make_trait_summary <- function (traits) {
  
  # Summarize traits
  traits %>%
    select(-taxon) %>%
    colnames() %>%
    tibble(trait = .) %>%
    mutate(trait = str_split(trait, "_") %>% map_chr(1)) %>%
    mutate(
      trait = case_when(
        trait == "frond" ~ "frond_width",
        trait == "number" ~ "number_pinna_pairs",
        trait == "stipe" ~ "stipe_length",
        TRUE ~ trait
      )
    ) %>%
    count(trait) %>%
    mutate(trait_type = case_when(
      trait %in% c("frond_width", "number_pinna_pairs", "stipe_length") ~ "continuous",
      trait == "indusium" ~ "binary",
      TRUE ~ "qualitative"
    )) %>%
    arrange(trait_type, trait)
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


#' Match trait data and tree
#' 
#' Order of species in traits will be rearranged to match the
#' phylogeny.
#'
#' @param traits Dataframe of traits, with 'species' column and
#' additional columns, one for each trait
#' @param phy Phylogeny (list of class "phylo")
#' @param return Type of object to return
#'
#' @return Either a dataframe or a list of class "phylo"; the tree or
#' the traits, pruned so that only species occurring in both datasets
#' are included.
#' @export
#'
#' @examples
match_traits_and_tree <- function (traits, phy, return = c("traits", "tree")) {
  
  assert_that("taxon" %in% colnames(traits))
  
  # Keep only species in phylogeny
  traits <- traits %>%
    filter(taxon %in% phy$tip.label) 
  
  # Trim to only species with trait data
  phy <- drop.tip(phy, setdiff(phy$tip.label, traits$taxon))
  
  # Get traits in same order as tips
  traits <- left_join(
    tibble(taxon = phy$tip.label),
    traits,
    by = "taxon"
  )
  
  # Make sure that worked
  assert_that(isTRUE(all.equal(traits$taxon, phy$tip.label)))
  
  # Return traits or tree
  assert_that(return %in% c("tree", "traits"))
  
  if(return == "tree") { 
    return (phy) 
  } else {
    return (traits)
  }
  
}

#' Analyze phylogenetic signal in a continuous trait of interest
#'
#' @param selected_trait Name of trait to analyze phylogenetic signal
#' @param traits Dataframe including all untransformed traits, with
#' 'taxon' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return List of estimated Blomberg's K and Pagel's lambda and
#' their significance
#' 
analyze_cont_phylosig <- function (selected_trait, traits, phy) {
  
  # Trim data to non-missing trait values and
  # make sure species in same order in tree and traits
  traits_select <- traits %>% select(taxon, all_of(selected_trait)) %>%
    remove_missing(na.rm = TRUE)
  
  traits_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "traits") 
  phy_trim <- match_traits_and_tree(traits = traits_select, phy = phy, "tree") 
  
  # Extract named vector of trait values for phylosig()
  trait_vec <- pull(traits_trim, selected_trait) %>%
    set_names(traits_trim$taxon)
  
  # Run phylosig() on selected trait
  # using Blomberg's K
  k_model <- phytools::phylosig(phy_trim, trait_vec, method = "K", test = TRUE)
  # and Pagel's lambda
  lambda_model <- phytools::phylosig(phy_trim, trait_vec, method = "lambda", test = TRUE)
  
  # get model results
  tibble(
    trait = selected_trait,
    kval = k_model$K,
    k_pval = k_model$P,
    lambda = lambda_model$lambda,
    lambda_pval = lambda_model$P,
    n = length(trait_vec)
  )
  
}

#' Analyze phylogenetic signal in binary traits
#'
#' @param traits Dataframe including all untransformed traits, with
#' 'species' as a column and other traits as other columns.
#' @param phy Phylogeny
#'
#' @return Dataframe of estimated Fritz and Purvis' D values and
#' associated statistics for each binary trait in `traits`
#' 
analyze_binary_phylosig <- function (traits, phy) {
  
  # Duplicated node labels are not allowed, so drop these
  if(any(duplicated(phy$node.label))) phy$node.label <- NULL
  
  # Nest by trait and loop phylo.d()
  # Note: some of the traits include NAs.
  # For phylo.d(), need to have matching tree and trait data with all NAs removed
  # Don't do comparative.data() on all the traits together, or species missing 
  # data for ANY trait will get dropped
  traits %>%
    gather(trait, value, -taxon) %>%
    nest(data = -trait) %>%
    # Remove NAs for individual trait observations
    mutate(data = map(data, ~remove_missing(., na.rm = TRUE))) %>%
    mutate(
      # Construct a comparative data object for each
      # binary trait
      comp_data = map(
        data, 
        ~ caper::comparative.data(
          phy = phy, 
          data = as.data.frame(.), 
          names.col= "taxon", 
          na.omit = TRUE)
      ),
      # Run phylo.d on each binary trait
      phylo_d_out = map(
        comp_data,
        ~caper::phylo.d(data = ., binvar = value) # phylo.d() uses NSE for `binvar`
      ),
      # Extract the useful bits of information from each model fit
      phylo_d_summary = map(
        phylo_d_out,
        ~tibble(
          num_present = pluck(., "StatesTable", 1), 
          num_absent = pluck(., "StatesTable", 2), 
          D = pluck(., "DEstimate"),
          prob_random = pluck(., "Pval1"), 
          prob_brownian = pluck(., "Pval0")
        )
      )
    ) %>%
    select(trait, phylo_d_summary) %>%
    unnest(cols = phylo_d_summary)
  
}

# Community diversity ----

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
calc_biodiv_random <- function (
  comm_df, phy = NULL, phy_alt = NULL, 
  trait_tree = NULL, trait_tree_alt = NULL, 
  null_model = c("frequency", "richness", "independentswap", "trialswap"), 
  n_iterations = 1000, metrics) {
  
  # Make sure selection of metrics is OK
  assert_that(is.character(metrics))
  assert_that(
    length(metrics) > 0,
    msg = "At least one biodiversity metric must be selected")
  assert_that(
    all(metrics %in% c("richness", "pd", "rpd", "fd", "rfd", "pe", "rpe")),
    msg = "Biodiversity metrics may only be selected from 'richness', 'pd', 'rpd', 'fd', 'rfd', 'pe', or 'rpe'"
  )
  
  assert_that(is.character(null_model))
  assert_that(
    null_model %in% c("frequency", "richness", "independentswap", "trialswap"),
    msg = "Null model may only be selected from 'frequency', 'richness', 'independentswap', or 'trialswap'"
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
  random_comm <- picante::randomizeMatrix(comm_df, null.model = null_model, iterations = n_iterations)
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
  richness <- NULL
  
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
  if ("richness" %in% metrics) richness <- vegan::specnumber(random_comm)
  
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
    rpe = rpe,
    richness = richness
  ) %>%
    # Only keep non-NULL results
    purrr::compact()
  
}

#' Count number of times one number is higher than others
#'
#' @param x Number to count
#' @param y Vector of numbers to compare
#' @param na.rm = Logical; should NA values in the comparison
#' vector be removed before making comparison?
#'
#' @return Number of times x is higher than y
#' 
#' @examples {
#' count_higher(4, 1:10)
#' count_higher(4, c(1:10, NaN))
#' }
count_higher <- function (x, y, na.rm = TRUE) {
  
  assertthat::assert_that(assertthat::is.number(x))
  assertthat::assert_that(is.numeric(y))
  
  # remove any NAs before making comparison
  if(isTRUE(na.rm)) y <- y[!is.na(y)]
  
  # if comparison is zero length, return NA
  if(length(y) == 0) return (NaN)
  
  sum((x > y))
}

#' Count number of times one number is lower than others
#'
#' @param x Number to count
#' @param y Vector of numbers to compare
#' @param na.rm = Logical; should NA values in the comparison
#' vector be removed before making comparison?
#'
#' @return Number of times x is lower than y
#' 
#' @examples {
#' count_lower(4, 1:10)
#' count_lower(NaN, 1:10)
#' }
count_lower <- function (x, y, na.rm = TRUE) {
  
  assertthat::assert_that(assertthat::is.number(x))
  assertthat::assert_that(is.numeric(y))
  
  # remove any NAs before making comparison
  if(isTRUE(na.rm)) y <- y[!is.na(y)]
  
  # if comparison is zero length, return NA
  if(length(y) == 0) return (NaN)
  
  sum((x < y))
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
      # Calculate SES
      rand_mean = purrr::map_dbl(random_values, ~mean(., na.rm = TRUE)),
      rand_sd = purrr::map_dbl(random_values, ~sd(., na.rm = TRUE)),
      obs_z = (obs_val - rand_mean) / rand_sd,
      # Count number of times observed value is higher than random values
      obs_c_upper = purrr::map2_dbl(.x = obs_val, .y = random_values, ~count_higher(.x, .y)),
      # Count number of times observed value is lower than random values
      obs_c_lower = purrr::map2_dbl(.x = obs_val, .y = random_values, ~count_lower(.x, .y)),
      # Count the number of non-NA random values used for comparison
      obs_q = purrr::map_dbl(random_values, ~magrittr::extract(., !is.na(.)) %>% length()),
      # Calculate p-value for upper tail
      obs_p_upper = obs_c_upper / obs_q,
      # Calculate p-value for lower tail
      obs_p_lower = obs_c_lower / obs_q
    )
  
  colnames(results) <- paste(metric, colnames(results), sep = "_")
  
  results
}

#' Run randomization analysis for a set of biodiversity metrics
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
#' sample species richness.
#'   
#' @param comm_df Input community matrix in data.frame format (communities as rows,
#' species as columns, with row names and column names) (or a list containing this)
#' @param phy Input phylogeny with total branch length scaled to 1
#' @param null_model Name of null model to use. Must choose from 'frequency', 'richness', 
#' 'independentswap', or 'trialswap' (see picante::randomizeMatrix).
#' @param n_reps Number of random communities to replicate
#' @param n_iterations Number of iterations to use when swapping occurrences to
#' generate each random community
#' @param metrics Character vector names of metrics to calculate (or a list containing this). 
#' Must one or more of 'pd', 'rpd', 'fd', 'rfd', 'pe', or 'rpe'
#' @param dataset_name Name of the dataset
#'
#' @return Tibble. For each of the biodiversity metrics, the observed value (_obs), 
#' mean of the random values (_rand_mean), SD of the random values (_rand_sd), 
#' rank of the observed value vs. the random values (_obs_rank), standard effect size
#' (_obs_z), and p-value (_obs_p) are given.
#' 
run_rand_analysis <- function(comm_df, phy = NULL, trait_distances = NULL, null_model, n_reps, metrics, n_iterations = 10000, dataset_name) {
  
  # When running run_rand_analysis as a loop in {targets}, some inputs will be 
  # length-of-one lists. Unlist these.
  # If comm_df is a list, pull out the dataframe
  if(inherits(comm_df[1], "list")) comm_df <- comm_df[[1]]
  # If metrics is a list, pull out the character vector
  if(inherits(metrics[1], "list")) metrics <- metrics[[1]]
  
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
    phy_alt$edge.length <- rep(x = 1, times = length(phy_alt$edge.length))
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
    trait_tree_alt$edge.length <- rep(x = 1, times = length(trait_tree_alt$edge.length))
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
      calc_biodiv_random(comm_df, phy, phy_alt, trait_tree, trait_tree_alt, null_model = null_model, n_iterations = n_iterations, metrics = metrics)
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
    mutate(site = rownames(comm_df), dataset = dataset_name) %>%
    select(dataset, site, everything())
  
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
        (pe_obs_p_upper >= 0.95 | pe_alt_obs_p_upper >= 0.95) & rpe_obs_p_upper >= 0.975 ~ "paleo",
        (pe_obs_p_upper >= 0.95 | pe_alt_obs_p_upper >= 0.95) & rpe_obs_p_lower >= 0.975~ "neo",
        pe_obs_p_upper >= 0.99 | pe_alt_obs_p_upper >= 0.99 ~ "super",
        pe_obs_p_upper >= 0.95 | pe_alt_obs_p_upper >= 0.95 ~ "mixed",
        TRUE ~ "not significant"
      ),
      endem_type = factor(endem_type, levels = c("paleo", "neo", "not significant", "mixed", "super"))
    )
  
}


#' Categorize results of randomization test
#'
#' @param df Input dataframe
#'
#' @return Dataframe with siginificance of randomization results categorized
categorize_signif <- function (df) {
  
  df %>%
    mutate(
      pd_signif = case_when(
        pd_obs_p_upper > 0.99 ~ "> 0.99",
        pd_obs_p_upper > 0.975 ~ "> 0.975",
        pd_obs_p_lower > 0.99 ~ "< 0.01",
        pd_obs_p_lower > 0.975 ~ "< 0.025",
        TRUE ~ "not significant"
      ),
      pd_signif = factor(pd_signif, levels = c("< 0.01", "< 0.025", "not significant", "> 0.975", "> 0.99"))
    ) %>%
    mutate(
      rpd_signif = case_when(
        rpd_obs_p_upper > 0.99 ~ "> 0.99",
        rpd_obs_p_upper > 0.975 ~ "> 0.975",
        rpd_obs_p_lower > 0.99 ~ "< 0.01",
        rpd_obs_p_lower > 0.975 ~ "< 0.025",
        TRUE ~ "not significant"
      ),
      rpd_signif = factor(rpd_signif, levels = c("< 0.01", "< 0.025", "not significant", "> 0.975", "> 0.99"))
    ) %>%
    mutate(
      fd_signif = case_when(
        fd_obs_p_upper > 0.99 ~ "> 0.99",
        fd_obs_p_upper > 0.975 ~ "> 0.975",
        fd_obs_p_lower > 0.99 ~ "< 0.01",
        fd_obs_p_lower > 0.975 ~ "< 0.025",
        TRUE ~ "not significant"
      ),
      fd_signif = factor(fd_signif, levels = c("< 0.01", "< 0.025", "not significant", "> 0.975", "> 0.99"))
    ) %>%
    mutate(
      rfd_signif = case_when(
        rfd_obs_p_upper > 0.99 ~ "> 0.99",
        rfd_obs_p_upper > 0.975 ~ "> 0.975",
        rfd_obs_p_lower > 0.99 ~ "< 0.01",
        rfd_obs_p_lower > 0.975 ~ "< 0.025",
        TRUE ~ "not significant"
      ),
      rfd_signif = factor(rfd_signif, levels = c("< 0.01", "< 0.025", "not significant", "> 0.975", "> 0.99"))
    )
  
}

# Bioregions ----

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

# Conservation ----

#' Calculate the percent of area within grid cells of high biological diversity
#' that has protected status
#'
#' @param biodiv_ferns_spatial Spatial dataframe including biodiveristy metrics and grid cells
#' @param protected_areas Spatial dataframe including status of protected areas and their shapes
#' @param japan_shp Spatial dataframe including the total land area of Japan
#'
#' @return Dataframe with the percent of protected areas in grid cells with significantly high biodiversity
#' 
calculate_protected_area <- function(biodiv_ferns_spatial, protected_areas, japan_shp) {
  
  # Add percent rank for richness
  biodiv_ferns_spatial <- mutate(biodiv_ferns_spatial, richness_obs_p_upper = dplyr::percent_rank(richness))
  
  # Crop to Japan to spatial div results area
  japan_shp <- sf::st_crop(japan_shp, sf::st_bbox(biodiv_ferns_spatial))
  
  # Set CRS for biodiversity data and protected areas
  japan_crs <- sf::st_crs(japan_shp)
  
  biodiv_ferns_spatial <- sf::st_set_crs(biodiv_ferns_spatial, japan_crs)
  
  protected_areas <- sf::st_set_crs(protected_areas, japan_crs)
  
  # Make dummy variable so sf::as_Spatial() works
  japan_shp <- mutate(japan_shp, admin = "Japan")
  
  # Crop spatial data to only land regions within Japan map
  biodiv_ferns_spatial_cropped <- raster::intersect(sf::as_Spatial(biodiv_ferns_spatial), sf::as_Spatial(japan_shp)) %>% sf::st_as_sf()
  
  biodiv_ferns_spatial_conserv <- sf::st_join(biodiv_ferns_spatial_cropped, protected_areas, left = TRUE)
  
  # Calculate total area of cells with significantly high biodiversity
  signif_cells_total_area <-
    biodiv_ferns_spatial_cropped %>%
    mutate(area = st_area(.) %>% units::set_units(km^2)) %>%
    select(grids, richness_obs_p_upper, pd_obs_p_upper, fd_obs_p_upper, pe_obs_p_upper, area) %>%
    as_tibble() %>%
    select(-geometry) %>%
    pivot_longer(cols = contains("upper"), values_to = "p_value", names_to = "metric") %>%
    filter(p_value > 0.95) %>%
    mutate(metric = str_remove_all(metric, "_obs_p_upper")) %>%
    group_by(metric) %>%
    summarize(
      total_area = sum(area),
      .groups = "drop"
    )
  
  # Calculate area of protected zones within cells with significantly high biodiversity
  biodiv_ferns_spatial_conserv %>%
    select(grids, richness_obs_p_upper, pd_obs_p_upper, fd_obs_p_upper, pe_obs_p_upper, status, area) %>%
    as_tibble() %>%
    select(-geometry) %>%
    pivot_longer(cols = contains("upper"), values_to = "p_value", names_to = "metric") %>%
    filter(p_value > 0.95, status %in% c("medium", "high")) %>%
    mutate(metric = str_remove_all(metric, "_obs_p_upper")) %>%
    group_by(status, metric) %>%
    summarize(
      protected_area = sum(area),
      .groups = "drop"
    ) %>%
    left_join(signif_cells_total_area, by = "metric") %>%
    mutate(percent_protected = protected_area / total_area) %>%
    mutate(percent_protected = as.numeric(percent_protected))
  
}

# Plotting ----

map_theme <- function() {
  theme(
    panel.grid.major = element_line(color = "grey60", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "grey60"),
    axis.title = element_blank(),
    legend.position = "bottom"
  )
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


# Dummy function to track arbitary output from rmarkdown::render()
render_tracked <- function (tracked_output, ...) {
  rmarkdown::render(...)
}

#' Convert latex file to docx
#' 
#' Requires pandoc to be installed and on command line
#'
#' @param latex Path to input latex file.
#' @param docx Path to output docx file.
#' @param template Path to template docx file.
#' @param wd Working directory to run conversion. Should be same as
#' directory containing any files needed to render latex to pdf.
#'
#' @return List including STDOUT of pandoc; externally, the
#' docx file will be rendered in `wd`.
#' 
latex2docx <- function (latex, docx, template = NULL, wd = getwd()) {
  
  assertthat::assert_that(assertthat::is.readable(latex))
  
  assertthat::assert_that(assertthat::is.dir(fs::path_dir(docx)))
  
  latex <- fs::path_abs(latex)
  
  docx <- fs::path_abs(docx)
  
  template <- if (!is.null(template)) {
    glue::glue("--reference-doc={fs::path_abs(template)}")
  } else {
    NULL
  }
  
  processx::run(
    command = "pandoc",
    args = c("-s", latex, template, "-o", docx),
    wd = wd
  )
  
}
