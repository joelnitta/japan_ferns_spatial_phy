# README for raw data

This describes the raw data files used in the analysis. These are processed to
clean data (files in `data/`) using the script [`R/process_raw_data.R`](../R/process_raw_data.R)
All checksums are 32-byte MD5 hashes generated with tools::md5sum() in R.

- `forest_area_zip_files`: Folder containing 47 zip files, each with SHP files
containing the distribution of forest areas in Japan (one per prefecture)
downloaded from https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-A45.html.
Downloaded using script [`R/download_forest_shape.R`](../R/download_forest_shape.R).
For MD5 checksums, see [forest_area_zip_files_checksums.txt](forest_area_zip_files_checksums.txt)

- `world-clim`: Folder containing data downloaded from the WorldClim database at 2.5
minute resolution using the following R command on 2021-04-23.
`raster::getData("worldclim", download = TRUE, var = "bio", res = 2.5, path =
"data_raw/world_clim")`.
For MD5 checksums, see [world_clim_checksums.txt](world_clim_checksums.txt)

- `doi_10.5061_dryad.4362p32__v4.zip`: data files from Ebihara and Nitta 2019
downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32.
MD5 checksum: 4ae5355c524ab47336d3fe56bf7aa440

- `gm-jpn-all_u_2_2.zip`: Political map of Japan downloaded from
https://www.gsi.go.jp/kankyochiri/gm_japan_e.html.
MD5 checksum: d296fa71a08346c4c0a172cbf235e266

- `JP_pterid_excl_hyb200620.xlsx`: Raw occurrence points for pteridophytes in Japan
excluding hybrids.
MD5 checksum: 8ddc1cd9efd11ae72cfa7106003fea15

- `Lucid20210807.xlsx`: Trait data for pteridophytes in Japan exported from Lucid
(https://www.lucidcentral.org/). For more information on data formatting, see
https://help.lucidcentral.org/lucid/scoring-the-key/.
MD5 checksum: 4d3b891bd5eca15b0e0ce2e47ae7481e

- `map14-1.zip`: Shape and tiff files of Japanese deer distribution in Japan
downloaded from
https://www.biodic.go.jp/biodiversity/activity/policy/map/files/map14-1.zip on
2021-09-01. Provided by Japan Ministry of the Environment. For more information,
see https://www.biodic.go.jp/biodiversity/activity/policy/map/map14/index.html
(in Japanese).
MD5 checksum: 112300bf61659c96605ac25d97653785

- `map17.zip`: Shape files of protected areas in Japan downloaded from
https://www.biodic.go.jp/biodiversity/activity/policy/map/map17/index.html on
2020-07-15. Provided by Japan Ministry of the Environment.
MD5 checksum: 512db536dff2eb825fa9c1ca2c72a157

- `mesh2.zip`: Zipped shp file of second-degree mesh codes for Japan. Downloaded
from
http://gis.biodic.go.jp/BiodicWebGIS/Questionnaires?kind=mesh2&filename=mesh2.zip
on 2020-07-07.
MD5 checksum: 8d429db0b1a36c9132bc171094419d4c
