# Tests

# Run tests to make sure custom functions work as expected
# Returns an error if anything fails, otherwise NULL
run_tests <- function() {
  # extract_indep_vars() ----
  test_that("Extracting independent variables works", {
    expect_equal(
      extract_indep_vars("pe_obs_p_upper ~ temp + precip + precip_season + Matern(1 | long + lat)"),
      c("temp","precip","precip_season")
    )
    expect_equal(
      extract_indep_vars("pe_obs_p_upper ~   temp + I(temp^2) + precip + precip_season+ Matern(1|long + lat)"),
      c("temp","I(temp^2)","precip","precip_season")
    )
  })
}
