# Geocode points of interest in Japan

# This uses a live query of the OpenStreetMap Nominatim API (https://nominatim.org/). 
# In case the API changes in the future, save the data externally, and load for analysis.

# Set vectors of different areas to label
main_islands <- c("Hokkaido", "Shikoku", "Kyushu", "Honshu")
small_islands <- c("Okinawa", "Amamioshima", "Ogasawara", "Yakushima", "Izu Islands")
other <- c("Tokyo")

# Geocode for lat long
japan_points_raw <- tmaptools::geocode_OSM(c(main_islands, small_islands, other))

readr::write_csv(points_raw, "data_raw/japan_points_raw.csv")
