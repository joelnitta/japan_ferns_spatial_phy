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
    mutate_all(~str_replace_all(., "ë", "e"))
  
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
make_pd_highlight_map <- function (div_data, world_map, occ_data) {
  
  # gghighlight only "knows about" data in the most recent layer
  # So make plot with world map on top of highlighted SES of PD,
  # then rearrange layers.
  plot <-
    ggplot(div_data, aes(x = longitude, y = latitude)) +
    geom_tile(aes(fill = ses_pd), color = "black") +
    gghighlight( (pd_obs_p > 0.975 | pd_obs_p < 0.025) & !is.na(ses_pd) ) +
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
      fill = "SES of PD"
    )
  
}
