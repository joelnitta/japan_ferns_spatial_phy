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

#' Unzip Fern Tree of Life (FTOL data)
#'
#' @param zip_file Path to the data zip file downloaded from FigShare
#' @param unzip_path Path to directory to put the unzipped
#' contents (will be created if needed).
#' @return Unzipped data files:
#' - ftol_plastid_concat.fasta
#' - ftol_plastid_parts.csv
#'
unzip_ftol <- function (zip_file, exdir) {
  
  # Unzip only the needed files
  unzip(zip_file, "ftol_plastid_concat.fasta", exdir = exdir)
  unzip(zip_file, "ftol_plastid_parts.csv", exdir = exdir)
  
  c(
    fs::path(exdir, "ftol_plastid_concat.fasta"),
    fs::path(exdir, "ftol_plastid_parts.csv")
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

#' Subset a community to only taxa with reproductive mode data
#'
#' @param comm Community matrix with colnames as species
#' @param green_list Dataframe with columns `taxon` and `apomict`
#'
#' @return Subsetted community matrix
#' 
subset_comm_by_repro <- function (comm, repro_data) {
  
  taxa_with_repro <-
    repro_data %>%
    filter(!is.na(apomict)) %>%
    assert(not_na, apomict) %>%
    pull(taxon)
  
  comm[,colnames(comm) %in% taxa_with_repro]
  
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


#' Clean up reproductive mode data for pteridophytes of Japan
#'
#' @param repro_data_raw Tibble. Raw data read-in from Electronic Supp. Mat. 1
#' @param green_list Tibble. List matching taxon ID to species name.
#'
#' @return Tibble. Column `reproductive_mode` is character "unknown", "sexual", 
#' "apomictic", "sex_apo". Other columns are T/F for a particular growth form
#' (evergreen/seasonal green) or reproductive mode.
process_repro_data <- function (repro_data_raw, green_list) {
  
  # Format repro mode data
  repro_data_raw %>%
    clean_names %>% 
    # Reformat reproductive mode according to data codes
    assert(in_set(0,1,2,3), reproductive_mode) %>%
    assert(not_na, reproductive_mode) %>%
    mutate(
      reproductive_mode = case_when(
        reproductive_mode == 0 ~ "unknown",
        reproductive_mode == 1 ~ "sexual", 
        reproductive_mode == 2 ~ "apomictic",
        reproductive_mode == 3 ~ "sex_apo"
      ),
      # Assumes that there are no missing data 
      # (all empty cells mean the trait is not present)
      sexual_diploid = case_when(
        sexual_diploid == 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      sexual_polyploid = case_when(
        sexual_polyploid == 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      evergreen = case_when(
        evergreen == 1 ~ TRUE,
        TRUE ~ FALSE
      ),
      seasonal_green = case_when(
        seasonal_green == 1 ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    select(taxon_id, reproductive_mode, sexual_diploid, sexual_polyploid, evergreen, seasonal_green) %>%
    # Define as apomictic taxon if reproductive mode is apomictic or both
    mutate(
      taxon_id = as.character(taxon_id),
      apomict = case_when(
        reproductive_mode == "apomictic" ~ TRUE,
        reproductive_mode == "sex_apo" ~ TRUE,
        reproductive_mode == "sexual" ~ FALSE,
        # Set as NA if repro mode is unknown
        reproductive_mode == "unknown" ~ NA
      )) %>%
    # Add species names by joining green list
    left_join(
      select(green_list, taxon, taxon_id), by = "taxon_id"
    ) %>%
    # Make sure that worked correctly
    assert(not_na, taxon) %>%
    assert(is_uniq, taxon)
  
}

#' Load WorldClim data cropped to the area including Japan
#'
#' The WorldClim data first needs to be downloaded with this:
#' raster::getData("worldclim", download = TRUE, var = "bio", res = 2.5, path = path)
#'
#' @param path Path to data downloaded with raster::getData()
#'
#' @return Simple features dataframe with four climate variables 
#' 
load_ja_worldclim_data <- function(path) {
  
  # Load the WorldClim dataset at 2.5 minute resolution, crop it to the area around Japan
  raster::getData("worldclim", download = FALSE, var = "bio", res = 2.5, path = path) %>%
    raster::crop(raster::extent(122, 154, 20, 46)) %>%
    # convert to sf
    spex::polygonize() %>%
    # select relevant variables https://www.worldclim.org/data/bioclim.html
    select(
      temp = bio1, # Annual Mean Temperature
      temp_season = bio4, # Temperature Seasonality (standard deviation ×100)
      precip = bio12, # Annual Precipitation
      precip_season = bio15 # Precipitation Seasonality (Coefficient of Variation)
    ) %>%
    # Set CRS to NA (like other spatial dataframes in this analysis)
    st_set_crs(NA)
  
}

#' Join climate data to occurrence data,
#' calculate mean climate values per grid cell
#'
#' @param shape_ferns Simple features data frame; grid of fern occurrence points
#' generated by shape_from_comm_scaled_list()
#' @param ja_climate_data Simple features data frame; climate data in Japan at the
#' 2.5 minute scale
#'
#' @return Tibble with mean values of each climate variable per grid cell
#' 
calc_mean_climate <- function(shape_ferns, ja_climate_data) {
  
  # Check for correctly formatted columns
  assert_that("grids" %in% colnames(shape_ferns))
  assert_that("temp" %in% colnames(ja_climate_data))
  
  shape_ferns %>%
    select(grids) %>%
    st_join(ja_climate_data) %>%
    group_by(grids) %>%
    mutate(across(c(temp, temp_season, precip, precip_season), ~mean(., na.rm = TRUE))) %>%
    ungroup %>%
    select(grids, temp, temp_season, precip, precip_season) %>%
    unique() %>%
    st_drop_geometry()
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
  
  # Set up weighting
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
    # categorical version by name
    mutate(comp_trait = case_when(
      str_detect(trait, "^margin_") ~ "margin",
      str_detect(trait, "^shape_") ~ "shape",
      str_detect(trait, "^texture_") ~ "texture",
      TRUE ~ trait
    )) %>%
    # First weigh by combined trait
    add_count(comp_trait) %>% 
    mutate(weight_by_comp_trait = 1 / n) %>% 
    select(-n) %>%
    # Then weigh by trait type
    add_count(trait_type) %>%
    mutate(weight_by_type = 1 / n) %>%
    select(-n) %>%
    # Then weigh by combination of combined trait and trait type
    mutate(final_weight = weight_by_comp_trait * weight_by_type) %>%
    # Make sure weights have been applied properly: equal within combined trait
    # and equal across trait types (binary vs continuous)
    verify(sum(weight_by_comp_trait) == n_distinct(comp_trait)) %>%
    verify(sum(weight_by_type) == n_distinct(trait_type))
  
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

#' Make a trait dendrogram from a matrix of distances
#'
#' @param trait_distances Object of class "dist"; distances
#' measured from trait data
#'
#' @return List of class "phylo"; trait dendrogram
#' 
make_trait_tree <- function(trait_distances) {
  hclust(trait_distances, method = "average") %>%
    ape::as.phylo()
}

# Community diversity ----

#' Format the output of canaper::cpr_rand_test()
#'
#' @param df Dataframe; Output of canaper::cpr_rand_test()
#'
#' @return Tibble with columns corresponding to site names (`grids`), the observed
#' metric (`*_obs`), z-score of the metric (`*_obs_z`), and percentage of time the
#' observed value was observed greater than random simulations 
#' (`*_obs_p_upper`) or lower than random simulations (`*_obs_p_lower`)
#' 
format_cpr_res <- function(df) {
  df %>% as_tibble(rownames = "grids") %>%
    select(matches("grids|_obs$|_obs_z$|_obs_p_upper$|_obs_p_lower$"))
}

#' Calculate the percentage of apomictic fern species in each community
#'
#' @param comm_ferns Dataframe; fern communities with species as columns and sites as rows
#' @param repro_data Tibble; Reproductive mode data of ferns of Japan, including whether
#' each species in apomictic or not
#'
#' @return tibble; percent of apomictic species per community
calc_perc_apo <- function(comm_ferns, repro_data) {
  comm_ferns %>%
    rownames_to_column("grids") %>%
    as_tibble %>%
    pivot_longer(values_to = "abun", names_to = "taxon", -grids) %>%
    filter(abun > 0) %>%
    left_join(
      select(repro_data, taxon, apomict), by = "taxon"
    ) %>%
    # Drop species that lack reproductive mode data
    filter(!is.na(apomict)) %>%
    group_by(grids) %>%
    summarize(
      richness = n_distinct(taxon),
      n_apo = sum(apomict),
      .groups = "drop"
    ) %>%
    mutate(percent_apo = n_apo / richness) %>%
    # Run tests
    assert(not_na, everything()) %>%
    assert(within_bounds(0, 1), percent_apo) %>%
    select(-richness) # drop richness (won't need after joining other data)
}


#' Classify significance of randomization test
#' 
#' Modified version of canaper::cpr_classify_signif() that
#' does not check the name of the metric
#'
#' @param df Results of canaper::cpr_rand_test()
#' @param metric Selected metric to classify significance.
#' @param one_sided Logical; is hypothesis one-sided?
#' @param upper Logical; should values in the top 5% be counted?
#' If FALSE, values in the bottom 5% are counted. Only applies if
#' `one_sided` is TRUE.
#'
#' @return Dataframe
#' 
classify_signif <- function(df, metric, one_sided = FALSE, upper = FALSE) {
  
  df[[paste0(metric, "_obs_p_lower")]]
  
  if (!isTRUE(one_sided)) {
    signif <- dplyr::case_when(
      df[[paste0(metric, "_obs_p_lower")]] > 0.99 ~ "< 0.01",
      df[[paste0(metric, "_obs_p_lower")]] > 0.975 ~ "< 0.025",
      df[[paste0(metric, "_obs_p_upper")]] > 0.99 ~ "> 0.99",
      df[[paste0(metric, "_obs_p_upper")]] > 0.975 ~ "> 0.975",
      TRUE ~ "not significant"
    ) } else {
      if (isTRUE(upper)) {
        signif <- dplyr::case_when(
          df[[paste0(metric, "_obs_p_upper")]] > 0.99 ~ "> 0.99",
          df[[paste0(metric, "_obs_p_upper")]] > 0.95 ~ "> 0.95",
          TRUE ~ "not significant"
        )
      } else {
        signif <- dplyr::case_when(
          df[[paste0(metric, "_obs_p_lower")]] > 0.99 ~ "< 0.01",
          df[[paste0(metric, "_obs_p_lower")]] > 0.95 ~ "< 0.05",
          TRUE ~ "not significant"
        )
      }
    }
  
  df[[paste0(metric, "_signif")]] <- signif
  
  df
  
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
  biodiv_ferns_spatial <- mutate(biodiv_ferns_spatial, richness_obs_p_upper = dplyr::percent_rank(richness)) %>%
    select(-area) # drop latitudinal area (so we can join with protected area)
  
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
    as.DNAbin %>%
    as.matrix()
  
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
    ips::mafft(x = ., exec = "/usr/bin/mafft", options = "--adjustdirection") %>%
    remove_mafft_r %>%
    # Trim any column that is >90% gaps
    ips::deleteGaps(gap.max = nrow(.)*0.9)
  
  # Concatenate genes
  concatenate_genes(reduced)
}

# Remove sequences that are all missing from an alignment
# - helper function for load_ftol_alignment()
remove_blank_seqs <- function(dna_align) {
  
  # First convert to alignment to character
  aln_char <- as.character(dna_align)
  
  # Get a vector of unique bases for each sequence
  uni_chars <- apply(aln_char, 1, unique)
  
  # Only keep those that are not all missing
  not_missing <- uni_chars != "-"
  
  dna_align[not_missing, ]
}

#' Remove the "_R_" from mafft-generated alignments
#'
#' @param matrix Alignment matrix
#'
#' @return matrix
#' 
remove_mafft_r <- function (matrix) {
  rownames(matrix) <- rownames(matrix) %>% str_remove_all("_R_")
  matrix
}

#' Load list of aligned genes from the Fern Tree of Life (FTOL) project
#'
#' @param ftol_plastid_concat Path to concatenated alignment of all genes
#' @param ftol_plastid_parts Path to CSV file with start and end position
#' of each gene within concatenated alignment
#'
#' @return List of alignments (one per gene)
#' 
load_ftol_alignment <- function (ftol_plastid_concat, ftol_plastid_parts) {
  
  # Load concatenated sequences and start/end positions of each gene
  ftol_plastid_concat_seqs <- ape::read.FASTA(ftol_plastid_concat) %>% as.matrix()
  ftol_plastid_parts_dat <- readr::read_csv(ftol_plastid_parts)
  
  # Split the concatenated sequences into a list of sequences, one per gene
  ftol_plastid_parts_dat %>%
    mutate(subseq = map2(start, end, ~ftol_plastid_concat_seqs %>% magrittr::extract(, .x:.y))) %>%
    pull(subseq) %>%
    set_names(ftol_plastid_parts_dat$gene) %>%
    # Remove empty sequences from each gene
    map(remove_blank_seqs)
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
  
  stderr_path <- Sys.time() %>% 
    str_replace_all(" |:", "_") %>% paste0("treepl_cv_config", ., ".stderr") %>%
    fs::path(wd, .)
  
  stdout_path <- Sys.time() %>% 
    str_replace_all(" |:", "_") %>% paste0("treepl_cv_config", ., ".stdout") %>%
    fs::path(wd, .)
  
  # Run treePL
  processx::run("treePL", "treepl_cv_config", wd = wd, stderr = stderr_path, stdout = stdout_path)
  
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
  
  stderr_path <- Sys.time() %>% 
    str_replace_all(" |:", "_") %>% paste0("treepl_prime_config", ., ".stderr") %>%
    fs::path(wd, .)
  
  # Run treePL
  results <- processx::run("treePL", "treepl_prime_config", wd = wd, stderr = stderr_path)
  
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
  thorough = TRUE, wd) {
  
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
  
  stderr_path <- Sys.time() %>% 
    str_replace_all(" |:", "_") %>% paste0("treepl_", ., ".stderr") %>%
    fs::path(wd, .)
  
  stdout_path <- Sys.time() %>% 
    str_replace_all(" |:", "_") %>% paste0("treepl_", ., ".stdout") %>%
    fs::path(wd, .)
  
  # Run treePL, save stderr and stdout to working dir
  processx::run(
    "treePL", "treepl_config", 
    wd = wd, 
    stderr = stderr_path,
    stdout = stdout_path)
  
  # Read in tree
  ape::read.tree(fs::path(wd, outfile_path))
  
}

# Spatial models ----

#' Filter data for modeling
#'
#' Only rows with no missing values in the selected columns will be retained
#' 
#' @param data Input data
#' @param vars_keep Character vector of column names to keep
#'
#' @return Tibble
filter_data_for_model <- function(data, vars_keep) {
  select(data, all_of(vars_keep)) %>%
    ggplot2::remove_missing()
}

#' Calculate centroids from SF (simple features) data
#'
#' @param sf_data Dataframe of class "sf"
#'
#' @return Tibble with columns "lat" and "long" with latitude and
#' longitude of the centroid of each geometrical feature
#' 
sf_add_centroids <- function(sf_data) {
  
  # Extract centroids
  centroids <-  # Start with Simple feature collection ("sf") dataframe
    sf_data %>%
    # Calculate centroid of each geometry feature
    st_centroid() %>% 
    # Convert centroids to character (e.g., "c(140.9, 45.5)")
    mutate(geom_char = as.character(geometry)) %>% 
    # Drop geometry column
    sf::st_set_geometry(NULL) %>%
    as_tibble() %>%
    # Parse centroids to numeric
    separate(geom_char, c("long", "lat"), sep = ", ") %>%
    mutate(across(c(long, lat), parse_number)) %>%
    # Make sure it worked
    assert(not_na, long, lat) %>%
    assert(within_bounds(-180, 180), long) %>%
    assert(within_bounds(-90, 90), lat) %>%
    select(grids, long, lat)
  
  # Add centroids to original data
  left_join(sf_data, centroids, by = "grids")
}

#' Make a spatial weights list for testing spatial autocorrelation
#'
#' @param data Dataframe with lat and long of centroids of each site ("grids")
#'
#' @return List
#' 
make_dist_list <- function(data) {
  
  # Convert centroids to matrix
  centroids <- data %>%
    # Make sure needed columns are present
    verify(all(c("grids", "lat", "long") %in% colnames(data))) %>%
    # "grids" are the site names
    select(grids, long, lat) %>%
    # set rownames as sites so these carry through to list names in final result
    as_tibble() %>%
    column_to_rownames("grids") %>% 
    as.matrix()
  
  # Make inverted distance matrix (so points far away have high values)
  dist_mat <- 1/as.matrix(dist(centroids))
  
  # Set the diagonal (same site to same site) to zero
  diag(dist_mat) <- 0
  
  # If sites have the exact same lat/long,
  # replace Inf values for these with 0
  dist_mat[is.infinite(dist_mat)] <- 0
  
  # Convert distance matrix to spatial weights list
  spdep::mat2listw(dist_mat, row.names = rownames(dist_mat), style = "M")
  
}

#' Prepare data for calculating Moran's I in a loop
#'
#' @param morans_vars_env Variables to test for spatial autocorrelation in environmental dataset
#' @param biodiv_ferns_cent_env Tibble with metrics of biodiversity for environmental dataset
#' @param dist_list_env List of spatial distances from environmental dataset (spatial weights list)
#' @param morans_vars_repro Variables to test for spatial autocorrelation in repro. dataset
#' @param biodiv_ferns_cent_repro Tibble with metrics of biodiversity for repro. dataset
#' @param dist_list_repro List of spatial distances from repro. dataset (spatial weights list)
#'
#' @return Tibble with columns "vars" (character), "data_type" (character), "data" (list-column),
#' and "dist_list" (list-column)
#' 
prepare_data_for_moran <- function(
  morans_vars_env,
  biodiv_ferns_cent_env,
  dist_list_env,
  morans_vars_repro = "percent_apo",
  biodiv_ferns_cent_repro,
  dist_list_repro
) {
    tibble(
      vars = morans_vars_env,
      data_type = "env",
      data = list(biodiv_ferns_cent_env),
      dist_list = list(dist_list_env)
    )
}

#' Run permutation test for Moran's I statistic
#'
#' @param var_name String; name of variable to test
#' @param biodiv_data Dataframe including the variable to test
#' @param listw Spatial weights list for testing spatial autocorrelation
#' @param nsim Number of simulations to use for calculating Moran's I
#'
#' @return Tibble with Moran's I and p-value
run_moran_mc <- function(var_name, biodiv_data, listw, nsim = 1000) {
  spdep::moran.mc(
    x = biodiv_data[[var_name]], listw = listw, nsim = nsim) %>% 
    broom::tidy() %>%
    rename(morans_i = statistic, pval_i = p.value) %>%
    mutate(var = var_name)
}

#' Run modified t-test accounting for spatial autocorrelation in indepedent variables
#' 
#' @param biodiv_ferns_repro_spatial Spatial dataframe including climate percent apomixis
#' @param mean_climate Dataframe including grid cell IDs and climate variables
#' @param vars_select Selecte variables to include in analysis
#'
#' @return Tibble with modified t-test statistics
#' 
run_mod_ttest_ja <- function(biodiv_ferns_cent_repro, vars_select) {
  
  # Extract long/lat as matrix
  coords <- select(biodiv_ferns_cent_repro, long, lat) %>% as.matrix()
  
  # Make cross table of all unique combinations of independent vars
  t_test_vars <-
    biodiv_ferns_cent_repro %>%
    select(all_of(vars_select)) %>%
    colnames() %>%
    list(var1 = ., var2 = .) %>%
    cross_df() %>%
    # Exclude identical variables
    filter(var1 != var2) %>%
    rowwise() %>%
    # Filter to only unique combinations
    mutate(comb1 = paste(sort(c(var1, var2)), collapse = " ")) %>%
    ungroup() %>%
    select(comb1) %>%
    unique() %>%
    separate(comb1, c("var1", "var2"), sep = " ")
  
  # Helper function for running SpatialPack::modified.ttest() in a loop across variables
  run_mod_ttest <- function(var1, var2, data, coords) {
    res <- SpatialPack::modified.ttest(data[[var1]], data[[var2]], coords)
    tibble(
      var1 = var1,
      var2 = var2,
      p_value = res$p.value,
      corr = res$corr,
      f_stat = res$Fstat,
      dof = res$dof
    )
  }
  
  map2_df(t_test_vars$var1, t_test_vars$var2, ~run_mod_ttest(var1 = .x, var2 = .y, data = biodiv_ferns_cent_repro, coords = coords))
  
}

#' Prepare data for running spaMM in a loop
#' 
#' Includes a quadratic term for temperature only
#'
#' @param resp_var_env Response variables to include in spatial model for environmental dataset
#' @param biodiv_ferns_cent Tibble with metrics of biodiversity
#' @param resp_var_repro Response variables to include in spatial model for repro. dataset
#'
#' @return Tibble with columns "resp_var" (character), "formula" (character), "data" (list-column),
#' "data_type" (character)
#' 
prepare_data_for_spamm <- function(
  resp_var_env = c("fd_obs_z", "pd_obs_z", "rfd_obs_z", "richness", "rpd_obs_z"),
  biodiv_ferns_cent,
  resp_var_repro = c("pd_obs_z", "rpd_obs_z")
) {
  bind_rows(
    tibble(
      resp_var = resp_var_env,
      formula = glue("{resp_var_env} ~ temp + I(temp^2) + precip + precip_season + area + Matern(1|long+lat)")
    ),
    tibble(
      resp_var = resp_var_repro,
      formula = glue("{resp_var_repro} ~ percent_apo + precip + precip_season + area + Matern(1|long+lat)")
    )
  ) %>%
    mutate(data = list(biodiv_ferns_cent))
}

#' Fit a linear mixed model including spatial autocorrelation
#' 
#' Wrapper around spaMM::fitme() so it can be run as a loop in {targets}
#' 
#' A negative binomial model is fit for species richness; gaussian otherwise.
#'
#' @param formula Formula as a character string
#' @param data Data for the model
#' @param resp_var Response variable
#'
#' @return Tibble with columns: "resp_var" (character), "formula" (character),
#' "model" (list)
#' 
run_spamm <- function(formula, data, resp_var) {
  tibble(
    resp_var = resp_var,
    formula = formula,
    # Specify model family according to response variable
    mod_family = case_when(
      resp_var == "richness" ~ "negbin",
      resp_var == "pe_obs_signif" ~ "binomial",
      TRUE ~ "gaussian"
    ),
    model = list(spaMM::fitme(as.formula(formula), data = data, family = mod_family))
  )
}

#' Generate a table of formulas to compare with LRT
#' 
#' Each comparison drops an independent variable from the
#' full model, so the effect of that variable on the model
#' can be assessed.
#' 
#' Comparison will also be made with the null (spatial) model
#' 
#' Helper function for prepare_data_for_lrt()
#'
#' @param resp_var Name of response variable
#' @param indep_vars Vector of independent variables
#'
#' @return Tibble
generate_spatial_comparisons <- function(resp_var, indep_vars) {
  
  full_formula_string <- glue("{resp_var} ~ {paste(indep_vars, collapse = ' + ')} + Matern(1|long+lat)") %>% 
    as.character
  
  tibble(
    resp_var = resp_var,
    indep_var = indep_vars,
    full_formula = full_formula_string) %>%
    # Define string to remove from full formula to create null formula.
    # In case of temperature, remove both temp and quadratic term (temp^2)
    mutate(
      indep_var_rm = case_when(
        indep_var == "temp" ~ "temp + I(temp^2)",
        TRUE ~ indep_var
      ),
      indep_var_rm = fixed(paste(indep_var_rm, "+"))
    ) %>%
    mutate(null_formula = str_remove_all(full_formula_string, indep_var_rm) %>% str_replace_all("  ", " ")) %>%
    select(-indep_var_rm) %>%
    bind_rows(
      tibble(
        resp_var = resp_var,
        indep_var = "null_model",
        full_formula = full_formula_string,
        null_formula = glue("{resp_var} ~ 1 + Matern(1|long+lat)") %>% as.character
      )
    ) %>%
    rename(comparison = indep_var)
}

#' Extract independent variables from a formula string
#' 
#' Also drops Matern random effects
#' 
#' Helper function for prepare_data_for_lrt()
#'
#' @param formula_string Formula as a string
#'
#' @return Character vector of independent variables
#' 
extract_indep_vars <- function(formula_string) {
  formula_string %>%
    str_remove_all("\\+ *Matern\\(.*\\)") %>%
    str_remove_all("^[^~]+~") %>%
    str_remove_all(" ") %>%
    str_split("\\+") %>%
    unlist()
}

#' Make a tibble for comparing log-likelihood values between spatial models
#' 
#' Constructs sets of comparisons between full model and smaller models.
#' 
#' Each comparison drops an independent variable from the
#' full model, so the effect of that variable on the full model
#' can be assessed.
#' 
#' Comparison will also be made with the null (spatial) model.
#'
#' @param spatial_models Tibble of spatial model results with columns
#' 'resp_var', 'formula', and 'log_lik'
#'
#' @return Tibble in wide format with columns 'resp_var', 'comp_group', 'null_formula',
#' and 'full_formula'
#' 
prepare_data_for_lrt <- function(spatial_models, biodiv_ferns_cent) {
  
  # Extract response and independent variables from models
  spatial_models %>%
    mutate(indep_vars = map(formula, extract_indep_vars)) %>%
    select(resp_var, indep_vars) %>%
    # Construct formulas for comparing full vs. null model
    mutate(comp = map2(resp_var, indep_vars, generate_spatial_comparisons)) %>%
    select(-indep_vars, -resp_var) %>%
    unnest(comp) %>%
    # Add data
    mutate(data = list(biodiv_ferns_cent)) %>%
    assert(not_na, everything())
  
}

#' Run a likelihood ratio test on formulas provided as a tibble
#'
#' @param null_formula Null formula as character vector for spaMM::fixedLRT()
#' @param full_formula Full formula as character vector for spaMM::fixedLRT()
#' @param data Data for model
#' @param resp_var Name of response variable
#' @param comparison Name of variable that is being compared (present in full
#' formula but absent in null formula)
#'
#' @return Tibble with columns 'chi2_LR', 'df', 'p_value', 'loglik_null', 
#' 'loglik_full', 'resp_var', and 'comparison'
#' 
run_spamm_lrt <- function(null_formula, full_formula, data, resp_var, comparison) {
  
  # Conduct LRT
  lrt_res <- spaMM::fixedLRT(
    null.formula = as.formula(null_formula), 
    formula = as.formula(full_formula), 
    data = data, 
    method = "ML")
  
  # Extract important statistics from result
  # (chi2, df, p-value, log-likelihoods of null and full model)
  lrt_res %>%
    magrittr::extract2("basicLRT") %>%
    as_tibble() %>%
    mutate(
      loglik_null = logLik(lrt_res$nullfit),
      loglik_full = logLik(lrt_res$fullfit),
      null_formula = null_formula,
      full_formula = full_formula,
      resp_var = resp_var, 
      comparison = comparison)
}

#' Extract beta table (fixed effects) from a model
#'
#' @param model An object of class HLfit, as returned by the fitting functions in spaMM.
#'
#' @return Tibble with fixed effects
#' @examples
#' betas <-
#' spatial_models %>%
#'   mutate(beta_table = map(model, get_beta_table)) %>%
#'   select(resp_var, beta_table) %>%
#'   unnest(beta_table)
#' 
get_beta_table <- function(model) {
  
  # Make version of summary() that won't print to screen
  quiet_summary <- quietly(summary)
  
  model %>%
    # get model summary
    quiet_summary %>%
    # extract the beta table (a matrix)
    pluck("result", "beta_table") %>%
    # convert to tibble
    as.data.frame() %>%
    rownames_to_column("term") %>%
    clean_names() %>%
    as_tibble()
}

#' Extract correlation parameters from a model
#'
#' @param model An object of class HLfit, as returned by the fitting functions in spaMM.
#'
#' @return Tibble with correlation parameters (rhu and nu)
#' 
get_corr_pars <- function(model) {
  get_ranPars(model, which="corrPars") %>%
    pluck(1) %>%
    unlist() %>%
    as.matrix() %>%
    t() %>%
    as_tibble()
}

#' Extract AIC from a model
#'
#' @param model An object of class HLfit, as returned by the fitting functions in spaMM.
#'
#' @return Tibble with AIC values
#' 
get_aic <- function(model) {
  AIC.HLfit.quiet <- quietly(spaMM::AIC.HLfit)
  AIC.HLfit.quiet(model) %>%
    pluck("result") %>%
    as.matrix() %>%
    t() %>%
    as_tibble() %>%
    janitor::clean_names()
}

#' Summarize statistics of spatial models
#' 
#' @param spatial_models Tibble, each with a spatial model.
#' Output of run_spamm()
#'
#' @return Tibble including columns: `resp_var` (response variable),
#' `formula` (text string of formula used for model), `data_type` (dataset used for
#' model; 'env' for environmental, 'repro' for reproductive), `loglik` (log-likelihood),
#' `resid_var_phi` (residual variance), `nu` (correlation parameter of Matern matrix),
#'  `rho` (correlation parameter of Matern matrix)
#' 
get_model_stats <- function(spatial_models) {
  spatial_models %>%
    # Extract log likelihood and residual variance (phi)
    mutate(
      loglik = map_dbl(model, logLik),
      resid_var_phi = map_dbl(model, ~residVar(., "phi") %>% unique)
    ) %>%
    # Extract nu and rho
    mutate(map_df(.$model, get_corr_pars)) %>%
    # Ditch the model
    select(-model)
}

#' Summarize parameters of spatial models
#' 
#' Filters to only environmental models (those based on temperature and precipitation)
#'
#' @param spatial_models Tibble, each with a spatial model.
#' Output of run_spamm()
#' @param lrt_comp_table Tibble, results of likelihood ratio test (LRT)
#' Output of run_spamm_lrt()
#'
#' @return Tibble including columns `resp_var` (response variable), `term` (fixed effect parameter name),
#'  `estimate` (slope of parameter) `cond_se` (conditional standard error), `t_value`, 
#'  `chi2_LR` (chi-squared value of LRT), `df`
#' 
get_model_params <- function(spatial_models) {
  spatial_models %>%
    # Extract fixed effects
    mutate(fixed_effects = map(model, get_beta_table)) %>%
    # Ditch the model
    select(-model) %>%
    unnest(fixed_effects)
}

#' Compare between environmental and reproductive models using cAIC
#'
#' @param spatial_models Tibble, each with a spatial model.
#' Output of run_spamm()
#'
#' @return Tibble including columns `resp_var` (response variable),
#' `model_type` (% apomictic taxa or temperature), `marginal_aic`, `conditional_aic`,
#'  `dispersion_aic`, `effective_df`,  `loglik` (log-likelihood), 
#'  and `resid_var_phi` (residual variance)
compare_aic_env_repro <- function(spatial_models) {
  # Filter to models based on reproductive mode dataset
  spatial_models %>%
    filter(resp_var %in% c("pd_obs_z", "rpd_obs_z")) %>%
    # Extract model type (temp or % apomictic taxa)
    tidyr::extract(formula, "model_type", "~ ([^ ]+) ") %>%
    # Get AIC values and residual variance (phi) for each model
    mutate(
      map_df(model, get_aic),
      loglik = map_dbl(model, logLik),
      resid_var_phi = map_dbl(model, ~residVar(., "phi") %>% unique)
    ) %>%
    # Sort by conditional AIC
    group_by(resp_var) %>%
    arrange(conditional_aic, .by_group = TRUE) %>%
    ungroup() %>%
    # Ditch the model
    select(-model) 
}

#' Run Grubb's test for outliers
#'
#' Returns the results in a tidy dataframe
#'
#' @param x a numeric vector for data values. 
#' @param type Integer value indicating test variant. See `?outliers::grubbs.test`
#' @param ... Other arguments passed to `outliers::grubbs.test()`
#'
#' @return Tibble
grubbs_tidy <- function(x, type = 10, ...) {
  
  outliers::grubbs.test(x = x, type = type, ...) %>%
    broom::tidy() %>%
    dplyr::mutate(stat_name = c("g", "u")) %>%
    tidyr::pivot_wider(values_from = "statistic", names_from = "stat_name") %>%
    janitor::clean_names() %>%
    dplyr::select(g, u, p_value, method, alternative)
  
}


#' Crops an area then calculate its area in sq km
#' 
#' Helper function for calc_area_by_lat()
#'
#' @param x object of class sf or sfc
#' @param ymin minimum y extent of cropping area
#' @param ymax maximum y extent of cropping area
#' @param xmin minimum y extent of cropping area
#' @param xmax maximum y extent of cropping area
#'
#' @return Area in sq km of the cropped area
#' 
area_from_crop <- function(x, ymin, ymax, xmin, xmax) {
  suppressMessages(sf::st_crop(x = x, ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax)) %>%
    sf::st_area() %>%
    units::set_units(km^2) %>%
    as.numeric() %>%
    # sf::st_area returns vec of length zero if no area
    # convert those to an actual zero
    ifelse(length(.) == 0, 0, .)
}

#' Calculate area in latitudinal bands of a geometric shape
#' 
#' Slices the shape into bands of width `lat_cut`, then calculates
#' the area of each band and a rolling mean average of width `lat_window`
#'
#' @param shp Geometric shape; object of class sf or sfc
#' @param lat_cut Width of latitude to slice the shape into
#' @param lat_window Width of rolling window to calulate rolling mean area
#'
#' @return Tibble with columns `ymin`, `ymax`, `area`
calc_area_by_lat <- function(shp, lat_cut = 0.2, lat_window = 1) {
  
  # Get bounding box of the input shape (xmin, xmax, ymin, ymax)
  bounds <- st_bbox(shp) 
  
  # Define latitudes to cut along: every 0.2 degree
  lat_cuts <- seq(
    floor(bounds$ymin), 
    ceiling(bounds$ymax), 
    0.2)
  
  # Prepare tibble for looping
  dat <- tibble(
    ymin = lat_cuts[1:(length(lat_cuts) - 1)],
    ymax = lat_cuts[2:length(lat_cuts)],
    xmin = bounds$xmin %>% floor,
    xmax = bounds$xmax %>% ceiling,
    x = list(shp)
  )
  
  # Determine rolling mean window in units of input data
  # k x lat_cut = lat_window,
  # so lat_window / lat_cut = 1 degree / 0.2 degree = 5
  width_k <- lat_window / lat_cut
  
  # Calculate area in 0.2 degree latitudinal bands,
  # then rolling mean in 1 degree bands
  dat %>%
    mutate(area = pmap_dbl(., area_from_crop)) %>%
    select(ymin, ymax, area) %>%
    # for the missing values on either end, just extend from the first non-missing value
    mutate(area = zoo::rollmean(x = area, k = width_k, fill = c("extend", NA, "extend")))
  
}


#' Add rolling mean area of latidudinal bands to biodiv data
#'
#' @param biodiv_ferns_spatial Spatial dataframe including biodiveristy metrics and grid cells
#' @param lat_area_ja Tibble with area (km2) in latitudinal bands in Japan (ymin and ymax indicate min and
#' max latitude of each band)
#'
#' @return Tibble
#' 
add_roll_area <- function(biodiv_ferns_spatial, lat_area_ja) {
  # interval_inner_join() only works (properly) on integers
  # lat has max 1 decimal, so multiply by 10 to convert to integer
  area_mapped_to_centroids <-
    biodiv_ferns_spatial %>% 
    mutate(ymin = lat*10,
           ymin = as.integer(ymin)) %>%
    mutate(ymax = ymin) %>%
    # fuzzyjoin doesn't like spatial dataframes
    sf::st_set_geometry(NULL) %>%
    select(grids, ymin, ymax) %>%
    fuzzyjoin::interval_inner_join(
      mutate(lat_area_ja, ymin = as.integer(ymin*10), ymax = as.integer(ymax*10)),
      type = "within",
      by = c("ymin", "ymax")) %>%
    # make sure the join worked properly
    verify(nrow(.) == nrow(biodiv_ferns_spatial)) %>%
    select(-matches("ymin|ymax"))
  
  left_join(biodiv_ferns_spatial, area_mapped_to_centroids, by = "grids") %>%
    assert(is_uniq, grids) %>%
    assert(not_na, area)
}

#' Drop one outlier value from biodiv data:
#' single extremely high percent_apo (> 60%) due to small number of species
#'
#' @param data data on biodiversity of ferns in Japan, including
#' columns `percent_apo`, `richness`, others
#'
#' @return data with one row dropped
#' 
drop_apo_outlier <- function (data) {
  data %>%
  verify(max(percent_apo) > 0.6) %>%
  filter(percent_apo != max(percent_apo)) %>%
  verify(max(percent_apo) < 0.6) %>%
  verify(nrow(.) == (nrow(data) - 1))
}


#' Predict fitted values for a spatial model
#'
#' @param spatial_models One row of the spatial_models data frame
#'
#' @return Tibble with four columns: the name of the response variable,
#' the name of the independent variable ussed for fitting,
#' the formula used for the model, and a list-column of fitted values.
#' The fitted values include the independent variable, the predicted response, and 
#' 95% CI lower and upper bounds.
#' 
predict_fit <- function(spatial_models) {
  
  indep_var <- case_when(
    str_detect(spatial_models$formula, "temp") ~ "temp",
    str_detect(spatial_models$formula, "percent_apo") ~ "percent_apo",
    TRUE ~ NA_character_
  )
  
  model <- spatial_models$model[[1]]
  
  fit <- spaMM::pdep_effects(model, focal_var = indep_var)
  
  colnames(fit)[colnames(fit) == "focal_var"] <- indep_var
  colnames(fit)[colnames(fit) == "pointp"] <- model$predictor[[2]] %>% as.character()
  
  tibble(
    resp_var = spatial_models$resp_var,
    indep_var = indep_var,
    formula = spatial_models$formula,
    fit = list(as_tibble(fit))
  )
  
}

#' Add column add selected predictor variable to spatial models
#'
#' @param spatial_models Tibble with spatial models
#'
#' @return `spatial_models` with column `indep_var_select` indicating the name of the
#' independent variable to use for plotting
#' 
add_selected_pred_to_model <- function(spatial_models) {
  
  spatial_models %>%
    mutate(indep_var_select = case_when(
      str_detect(formula, "temp") ~ "temp",
      str_detect(formula, "percent_apo") ~ "percent_apo",
      TRUE ~ NA_character_
    )) %>%
    assert(not_na, indep_var_select)
  
}

# Manuscript rendering ----

#' Rename response variables in data frame
#' 
#' For formatting results tables
#'
#' @param df 
#'
#' @return df with renamed variables
#' 
rename_resp_vars <- function(df) {
  
  lookup <-
    tibble(
      resp_var = c("richness", "fd_obs_z", "pd_obs_z", "rfd_obs_z", "rpd_obs_z"),
      resp_var_print = c("Richness", "SES of FD", "SES of PD", "SES of RFD", "SES of RPD")
    )
  
  df %>%
    left_join(lookup, by = "resp_var") %>%
    select(-resp_var) %>%
    rename(resp_var = resp_var_print)
  
}

#' Rename model type in data frame
#' 
#' For formatting results tables
#'
#' @param df 
#'
#' @return df with renamed model types
#' 
rename_model_type <- function(df) {
  
  lookup <-
    tibble(
      model_type = c("percent_apo", "temp"),
      model_type_print = c("Reproductive", "Environmental")
    )
  
  df %>%
    left_join(lookup, by = "model_type") %>%
    select(-model_type) %>%
    rename(model_type = model_type_print)
  
}

#' Clean a vector of variables for printing
#'
#' @param x Character vector (variables)
#'
#' @return Character vector
#' 
clean_vars <- function(x) {
  x %>%
    gsub("^([^_]+)_obs_z", "SES of \\U\\1", ., perl=TRUE) %>%
    str_replace_all("^percent_apo$", "% apomictic taxa") %>%
    str_replace_all("^temp$", "Temperature") %>%
    str_replace_all("^temp_season$", "Temperature seasonality") %>%
    str_replace_all("^precip_season$", "Precipitation seasonality") %>%
    str_replace_all("^precip$", "Precipitation") %>%
    str_replace_all("^area$", "Area") %>%
    str_replace_all("^richness$", "Richness")
}

#' Make a tibble matching the nodes of a dendrogram to the cluster they belong to
#'
#' @param dendro Dendrogram; list of class "phylo". Dendrogram clustering sites. Part of output of phyloregion::phyloregion()
#' @param tax_clusters Tibble with two columns, "grids" (site names) and "cluster" (ID of cluster for each site)
#' @param cluster_select Value indicating the cluster to select for making the tibble
#'
#' @return Tibble
#' 
map_cluster_to_nodes <- function(dendro, tax_clusters, cluster_select) {
  
  # Get vector of tips matching the selected cluster
  tips <- tax_clusters %>% filter(cluster == cluster_select) %>% pull(grids)
  
  if(length(tips) == 1) {
    # If only one tip, just return that node (tip) number
    nodes <- which(dendro$tip.label == tips)
  } else if (length(tips) > 1) {
    # If multiple tips, get MRCA, then all descendents
    node_mrca <- ape::getMRCA(dendro, tips)
    nodes <- phytools::getDescendants(dendro, node_mrca)
  } else {
    stop("Can't find that cluster")
  }
  
  # Return a tibble mapping nodes to the cluster
  tibble(node = nodes, cluster = cluster_select)
}
