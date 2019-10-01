#### Ecostructure ####

# This submits each ecostructure target as a separate job on a cluster.
# Only run after running `plan.R` and saving `dispersion_fields_matrix_japan`
# as an RDS file in `data/`.

# Ecostructure plan for Japan pteridophytes dataset
ecos_plan <- drake_plan(
  # Geographic motifs
  ecos_fit_no_obs_geo = target(
    run_ecos(
      out_dir = "results", motif = "geo", dataset = "pteridos-no-obs",
      dat = readRDS(file_in("data/dispersion_fields_matrix_japan.RDS")),
      K = K, tol = tol, num_trials = num_trials),
    mem_type = "medium",
    transform = cross(K = !!(2:10), num_trials = 1, tol = 0.1),
  )
) %>%
  # Fill-in memory settings for running on cluster
  specify_resources()
