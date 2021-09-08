# README for raw data

This describes the raw data files used in the analysis. These are processed to
clean data (files in `data/`) using the script `R/process_raw_data.R`

forest_area_zip_files: Folder containing 47 zip files, each with SHP files
containing the distribution of forest areas in Japan (one per prefecture)
downloaded from https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-A45.html.
Downloaded using script `R/download_forest_shape.R`

world-clim: Folder containing data downloaded from the WorldClim database at 2.5
minute resolution using the following R command on 2021-04-23.
`raster::getData("worldclim", download = TRUE, var = "bio", res = 2.5, path =
"data_raw/world_clim")`

doi_10.5061_dryad.4362p32__v4.zip: data files from Ebihara and Nitta 2019
downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32

gm-jpn-all_u_2_2.zip: Political map of Japan downloaded from
https://www.gsi.go.jp/kankyochiri/gm_japan_e.html

JP_pterid_excl_hyb200620.xlsx: Raw occurrence points for pteridophytes in Japan
excluding hybrids.

Lucid20210807.xlsx Trait data for pteridophytes in Japan exported from Lucid
(https://www.lucidcentral.org/). For more information on data formatting, see
https://help.lucidcentral.org/lucid/scoring-the-key/

map14-1.zip: Shape and tiff files of Japanese deer distribution in Japan
downloaded from
https://www.biodic.go.jp/biodiversity/activity/policy/map/files/map14-1.zip on
2021-09-01. Provided by Japan Ministry of the Environment. For more information,
see https://www.biodic.go.jp/biodiversity/activity/policy/map/map14/index.html
(in Japanese).

map17.zip: Shape files of protected areas in Japan downloaded from
https://www.biodic.go.jp/biodiversity/activity/policy/map/map17/index.html on
2020-07-15. Provided by Japan Ministry of the Environment.

mesh2.zip: Zipped shp file of second-degree mesh codes for Japan. Downloaded
from
http://gis.biodic.go.jp/BiodicWebGIS/Questionnaires?kind=mesh2&filename=mesh2.zip
on 2020-07-07.
