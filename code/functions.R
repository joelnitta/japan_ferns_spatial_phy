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
  
# Community diversity ----

#' Calculate species richness for all grid cells
#' 
#' The long / lats of the occurrence data
#' are already cell centroids, so it is trivial
#' to compute species richness.
#'
#' @param data Occurrence data, where each row is the 
#' occurrence of a species in a grid cell
#' @return tibble
make_richness_matrix <- function (data) {
  data %>%
  group_by(latitude, longitude) %>%
  summarize(
    richness = n()
  )
}

#' Calculate Faith's PD for a single community in the community matrix
#' 
#' Does not include the root when summing branch lengths. At least two taxa in the 
#' community must be present in the tree
#' 
#' @param single_comm row of a community matrix
#' @param phy phylogenetic tree
#' 
#' @return a number -- Faith's index for that community
calc_faithsPD_single <- function(single_comm, phy) {
  # use ape::drop.tip directly without geiger::treedata()
  phy <- ape::drop.tip(phy, phy$tip.label[!(phy$tip.label %in% names(single_comm[single_comm > 0]))])
  
  # get sum of branches remaining
  sum(phy$edge.length)
}

#' Calculate Faith's PD for all communities in the community matrix
#' 
#' Does not include the root when summing branch lengths. At least two taxa in the 
#' community must be present in the tree.
#' 
#' @param comm community matrix, with species as columns and communities as rows.
#' @param phy phylogenetic tree
#' 
#' @return a numeric vector -- Faith's index for each community in the matrix.
calc_faithsPD <- function (comm, phy) {
  apply(comm, MARGIN = 1, FUN = calc_faithsPD_single, phy = phy)
}

#' Calculates Faith's PD for all communities in the community matrix 
#' after shuffling the tips of the phylogeny
#' 
#' Does not include the root when summing branch lengths. At least two taxa in the 
#' community must be present in the tree. This function is meant to be used 
#' to produce null distributions of communities.
#' 
#' @param comm community matrix, with species as columns and communities as rows.
#' @param phy phylogenetic tree
#' 
#' @return a numeric vector -- Faith's index for each community in the matrix, 
#' after randomly shuffling the tips of the tree.
pd_shuffle <- function (comm, phy) {
  calc_faithsPD(comm, picante::tipShuffle(phy))
}

# ses_pd_shuffle

#' Calculate the standard effect size (SES) of Faith's PD across multiple communities.
#' 
#' Does not include the root when summing branch lengths. At least two taxa in each
#' community must be present in the tree. Null communities are simulated by randomly shuffling 
#' taxon names at the tips of the tree.
#' 
#' @param comm community matrix, with species as columns and communities as rows.
#' @param phy phylogenetic tree
#' @param runs number of null communities to simulate
#' 
#' @return a data frame -- including number of taxa (ntaxa), observed Faith's PD (pd_obs), 
#' mean value of PD across all simulated communities (pd_rand_mean), standard deviation of PD 
#' across all simulated communities (pd_rand_sd), rank of the observed PD relative to simulated 
#' values (pd_obs_rank), standard effect size of PD (ses_pd), and probability of the observed value (pd_obs_p).
ses_pd_shuffle <- function (comm, phy, runs) {
  pd_obs <- calc_faithsPD(comm, phy)
  pd_rand <- replicate(runs, pd_shuffle(comm, phy))
  pd_rand_mean <- apply(pd_rand, 1, mean, na.rm = TRUE)
  pd_rand_sd <- apply(pd_rand, 1, sd, na.rm = TRUE)
  ses_pd <- (pd_obs - pd_rand_mean) / pd_rand_sd
  pd_obs_rank <- apply(cbind(pd_obs, pd_rand), 1, function (x) rank(x)["pd_obs"])
  pd_obs_p <- pd_obs_rank/(runs + 1)
  ntaxa <- apply(comm, 1, function (x) sum(x > 0))
  HUC <- comm$HUC
  tibble(HUC = HUC, 
         ntaxa = ntaxa, 
         pd_obs = pd_obs,
         pd_rand_mean = pd_rand_mean, 
         pd_rand_sd = pd_rand_sd, 
         pd_obs_rank = pd_obs_rank, 
         ses_pd = ses_pd,
         pd_obs_p = pd_obs_p,
         runs = runs)
}

# Plotting ----

#' Make a plot showing species richness on a map of Japan
#'
#' @param richness Species richness in 1km grid cells, 
#' with longitude and latitude
#' @param occ_data Occurrence data used to generate species richness matrix
#' @param world_map Background world map
#'
#' @return ggplot object
make_richness_plot <- function (richness, occ_data, world_map) {
  
  ggplot(world_map, aes(x = long, y = lat)) +
    geom_polygon(aes(group = group), fill = "light grey") +
    geom_tile(data = richness,
              aes(x = longitude, y = latitude, fill = richness)) + 
    coord_quickmap(
      xlim = c(pull(occ_data, longitude) %>% min %>% floor, 
               pull(occ_data, longitude) %>% max %>% ceiling),
      ylim = c(pull(occ_data, latitude) %>% min %>% floor, 
               pull(occ_data, latitude) %>% max %>% ceiling)
    ) +
    scale_fill_viridis_c(na.value="transparent") +
    jntools::blank_x_theme() +
    jntools::blank_y_theme() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    labs(
      fill = "Richness"
    )
  
}