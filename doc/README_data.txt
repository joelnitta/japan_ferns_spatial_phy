This README.txt file was generated on 09 September, 2021 by Joel Nitta

------------------- GENERAL INFORMATION -----------------

Title of Dataset: Data from: Spatial phylogenetics of Japanese ferns: Patterns,
processes, and implications for conservation

Author Information

Principal Investigator: Joel H. Nitta

Department of Biological Sciences, Graduate School of Science, The University of
Tokyo, 2-11-16 Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan joelnitta@gmail.com

Associate or Co-investigators: Brent D. Mishler, Wataru Iwasaki, Atsushi Ebihara

Geographic location of data collection: Japan

Information about funding sources or sponsorship that supported the collection
of the data: Funding provided in part by the Japan Society for the Promotion of
Science (Kakenhi) Grant Number 16H06279

--------------------------

SHARING/ACCESS INFORMATION

--------------------------

Licenses/restrictions placed on the data, or limitations of reuse: CC0 1.0
Universal (CC0 1.0)

Recommended citation for the data: Nitta JH, Mishler BD, Iwasaki W, Ebihara A
(2021) Data from: Spatial phylogenetics of Japanese ferns: Patterns, processes,
and implications for conservation FIXME: add DOI when available

Citation for and links to publications that cite or use the data: Nitta JH,
Mishler BD, Iwasaki W, Ebihara A (2021) Spatial phylogenetics of Japanese ferns:
Patterns, processes, and implications for conservation FIXME: add journal when
published

Code for analyzing the data is available on github:
https://github.com/joelnitta/japan_ferns_spatial_phy

--------------------

DATA & FILE OVERVIEW

--------------------

File list (filenames, directory structure (for zipped files) and brief
description of all data files):

•	_targets.tar.gz: Tarball (compressed folder) including all workflow results
produced by R targets package

•	japan_climate.gpkg: Climate data in Japan downloaded from WorldClim database

•	japan_deer_range.gpkg: Distribution maps of Japanese deer (Cervus nippon) in
Japan

•	japan_ferns_comm_full.csv: Community matrix (species x sites matrix) of
native, non-hybrid ferns in Japan, full (unfiltered) dataset

•	japan_ferns_inext_results.csv: Results of collection curve analysis on fern
specimens in Japan

•	japan_ferns_occ_summary.csv: Summary of specimen data of native, non-hybrid
ferns in Japan

•	japan_ferns_redundancy_by_res: Taxonomic richness, abundance, and redunancy at
different spatial resolutions for native, non-hybrid ferns of Japan

•	japan_ferns_shape_full.gpkg: Location of grid-cells (sites) for native,
non-hybrid ferns in Japan, full (unfiltered) dataset

•	japan_ferns_traits_lucid.csv: Morphological trait data of native, non-hybrid
ferns in Japan

•	japan_map_points.csv: Latitude and longitude of points used to produce a
labeled map of Japan

•	japan_map.gpkg: Map of Japan

•	japan_protected_areas.gpkg: Protected areas in Japan

•	testo_sundue_2016_calibrations.csv: Fossil calibration points used for
phylogenetic dating analysis

•	results.zip: Selected analysis results. Files include:

–	japan_ferns_biodiv.csv: Biodiversity statistics of native, non-hybrid ferns in
Japan

–	japan_ferns_comm.csv: Community matrix (species x sites matrix) of native,
non-hybrid ferns in Japan used for biodiversity analysis

–	japan_ferns_shape.gpkg: Location of grid-cells (sites) for native, non-hybrid
ferns in Japan used for biodiversity analysis

–	japan_ferns_traits.csv: Trait matrix of native, non-hybrid ferns in Japan used
for functional biodiversity analysis

–	japan_ferns_tree_dated.tre: Maximim-likelihood, ultrametric (dated) phylogeny
of native, non-hybrid ferns in Japan with clades consisting of identical OTUs
collapsed to polytomies

–	japan_ferns_tree_uncollapsed.tre: Maximim-likelihood, ultrametric (dated)
phylogeny of native, non-hybrid ferns in Japan with clades consisting of
identical OTUs not collapsed

–	japan_ferns_tree.tre: Maximim-likelihood phylogeny of native, non-hybrid ferns
in Japan

Additional data used in analysis not included in current data package:

•	doi_10.5061_dryad.4362p32__v4.zip: Data including DNA sequences, taxonomy, and
breeding system data for native Japanese ferns (Ebihara and Nitta 2019)
downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32.
MD5 checksum 4ae5355c524ab47336d3fe56bf7aa440

•	ftol_data_release_v0.0.1.zip: Fern Tree of Life project global DNA alignment
and phylogenetic tree (Nitta et al., in prep.) downloaded from
https://figshare.com/s/87b8ec102303bf9f3212. MD5 checksum
26b95c283eec6529c71b252a47d3a06c FIXME: replace with DOI when available

Checksums are 32-byte MD5 hashes generated with tools::md5sum() in R.

--------------------------

METHODOLOGICAL INFORMATION

--------------------------

Description of methods used for collection/generation of data:

A list of native, non-hybrid fern specimens housed at the herbarium of the
Museum of Science and Nature, Japan was converted to a community data matrix at
four grain sizes (grid-cells spanning Japan, each 0.1, 0.2, 0.3, or 0.4 degrees
per side). The 0.2 degree grain size was selected for further analysis based on
redundancy (ratio of number of specimens to number of taxa per cell).

All taxon names are based on the Green List (http://www.rdplants.org/gl/;
English version available at
https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32).

Traits were measured on each species as described in Ebihara and Nitta (2019).

Phylogenetic analysis was conducted with maximum likelihood in IQ-TREE v1.6.12
(Nguyen et al. 2015) by combining plastid rbcL sequences of each taxon with a
globally sampled data matrix (Nitta et al, in prep). Next, dating analysis was
carried out using treePL v1.0 (Smith and O’Meara 2012) with 26 fossil
calibration points after Testo and Sundue (2016). The dated phylogeny was then
trimmed to include Japanese taxa only.

The community matrix, traits, and phylogeny were used to analyze spatial
patterns of phylogenetic diversity and endemism.

Data files were generated from raw data (not included here) using scripts
available at https://github.com/joelnitta/japan_ferns_spatial_phy, in particular
https://github.com/joelnitta/japan_ferns_spatial_phy/blob/main/R/process_raw_data.R.

For full methods, see Nitta JH, Mishler BD, Iwasaki W, Ebihara A (2021) Spatial
phylogenetics of Japanese ferns: Patterns, processes, and implications for
conservation FIXME: Add journal when published

--------------------------

DATA-SPECIFIC INFORMATION

--------------------------

_targets.tar.gz: Tarball (compressed folder) including all workflow results
produced by R targets package. This is provided to enable inspection of workflow
steps without running the entire workflow from the beginning. To use it, unpack
the tar achive with the command "tar -xzf _targets.tar.gz“. Then, in R,
the”tar_load()" function in the R package “targets” can be used to load any
workflow step (target) defined in _targets.R
(https://github.com/joelnitta/japan_ferns_spatial_phy/blob/main/_targets.R). For
more information on the structure of the _targets folder and how to use it, see
https://github.com/ropensci/targets.

MD5 checksum: b779fa1e2fae29a3ef0b9cbfc9e31c5a

Corresponding commit in repo
(https://github.com/joelnitta/japan_ferns_spatial_phy): e2ba68b1b9cd4b926acddd7586a581854b2aad69

--------------------------

japan_deer_range.gpkg: Distribution maps of Japanese deer (Cervus nippon) in
Japan. Raw data were downloaded from the Japan Ministry of the Environment
(https://www.biodic.go.jp/biodiversity/activity/policy/map/map14/index.html).
The data include range maps based on three types of data: range of Japanese deer
surveyed in 1978, range surveyed in 2003, and estimated range inferred from a
model including snow cover and forest type based on the 2003 survey data. Only
estimated range with < 0.10 movement cost was included.

Number of variables: 2

Number of cases/rows: 3

Coordinate reference system: JGD2000

Variable list:

•	range: Source of range data. ‘1978’ = 1978 survey, ‘2003’ = 2003 survey,
‘estimated’ = range estimated from model based on 2003 survey data

•	geom: Vector describing shape and position of range polygon

Missing data codes: No missing data.

MD5 checksum: 18d0211e46420193bc1b030cbf802a52

--------------------------

japan_climate.gpkg: Climate data in Japan downloaded from WorldClim database
(https://worldclim.org/data/v1.4/formats.html). Raw data were downloaded at the
2.5 minute scale. The intersection of the raw climate data with 0.2 degree
grid-cells (japan_ferns_shape_full.gpkg) were identified, then average value for
each climatic variable was calculated for each grid-cell.

Number of variables: 5

Number of cases/rows: 73836

Coordinate reference system: JGD2000

Variable list:

•	temp: Mean annual temperature, in units of 10 * degrees celsius (BIO1)

•	temp_season: Temperature seasonality (standard deviation ×100) (BIO4)

•	precip: Annual precipitation (mm) (BIO12)

•	precip_season: Precipitation seasonlaity (coefficient of variation) (BIO15)

•	geom: Vector describing shape and position of grid-cell

Missing data codes: No missing data.

MD5 checksum: e37f5382f1c2f2cda6e1d0f29d09d675

--------------------------

japan_ferns_comm_full.csv: Community matrix (species x sites matrix) of native,
non-hybrid ferns in Japan. Sites correspond to 0.2 degree x 0.2 degree
grid-cells covering Japan. For location of grid-cells, see
japan_ferns_shape_full.gpkg. Full dataset including all grid-cells (none
excluded due to low redundancy).

Number of variables: 674

Number of cases/rows: 1386

Variable list:

•	grids: Name of grid-cell.

•	Other columns: Each column is named for a Japanese fern taxon. Values indicate
number of specimens that were observed in each grid-cell for that taxon.

Missing data codes: No missing data.

MD5 checksum: 24410fb49923043a2e72533670768f88

--------------------------

japan_ferns_inext_results.csv: Output from running iNEXT (Hsieh, Ma, and Chao
2016) on specimen data of native, non-hybrid ferns in Japan. iNEXT estimates
diversity (here, species richness) and sampling completeness based on
interpolation and extrapolation from the original data.

Number of variables: 9

Number of cases/rows: 40

Variable list:

•	m: Sample size (number of individuals)

•	method: Method used to estimate diveristy

•	order: Diversity order of Hill number (here 0, equivalent to species richness)

•	qD: Estimated diversity

•	qD.LCL: Lower 95% confidence limit on estimated diversity

•	qD.UCL: Upper 95% confidence limit on estimated diversity

•	SC: Estimated sampling completeness

•	SC.LCL: Lower 95% confidence limit on sampling completeness

•	SC.UCL: Upper 95% confidence limit on sampling completeness

Missing data codes: No missing data.

MD5 checksum: 99e0cd4ecb08aa86a4f31bd3fb5fde2c

--------------------------

japan_ferns_occ_summary.csv: Summary of raw specimen data for native, non-hybrid
ferns of Japan. Raw specimen data consists of one row per specimen, with species
name, latitude, longitude, and collection date (raw specimen data not included
here to protect endangered species). Raw data were filtered to remove duplicate
specimens and specimens that occurred outside of the second-degree mesh
(approximately 10 km sq grid-cell) map of Japan
(http://www.biodic.go.jp/trialSystem/EN/shpddl.html)

Number of variables: 3

Number of cases/rows: 4

Variable list:

•	dataset: Dataset name. ‘occ_point_data_ferns_unfiltered’ indicates data
includes all specimens, ‘occ_point_data_ferns’ indicates data only includes
specimens after filtering

•	variable: Variable. ‘n_taxa’ indicates number of taxa, ‘n_specimens’ indicates
number of specimens

•	value: Value corresponding to each variable

Missing data codes: No missing data.

MD5 checksum: 7dc8f58056b9a13386a4fc2ffc7ae17a

--------------------------

japan_ferns_redundancy_by_res.csv: Results of testing binning raw specimen data
of native, non-hybrid ferns of Japan into grid-cells at different resolutions.
Sampling redundancy is a metric of sampling completeness, calculated as 1 -
(richness/abundance).

Number of variables: 5

Number of cases/rows: 7064

Variable list:

•	res: Degree of resolution (grid-cell side length), from 0.1 to 0.4 degrees

•	grids: Name of grid-cell

•	abundance: Number of specimens occurring in grid-cell

•	richness: Number of taxa occurring in grid-cell

•	redundancy: Redundancy (measure of sampling completeness)

Missing data codes: No missing data.

MD5 checksum: 93676e34db248172e86a780db72b7db6

--------------------------

japan_ferns_shape_full.gpkg: Location of grid-cells (sites) for native,
non-hybrid ferns in Japan in GeoPackage format. Compiled from raw specimen data
by binning specimens into 0.2 degree grid-cells across Japan. Full dataset
including all grid-cells (none excluded due to low redundancy).

Number of variables: 5

Number of cases/rows: 1386

Coordinate reference system: JGD2000

Variable list:

•	grids: Name of grid-cell

•	abundance: Number of specimens occurring in grid-cell

•	richness: Number of taxa occurring in grid-cell

•	redundancy: Redundancy (measure of sampling completeness)

•	geom: Vector describing shape and position of grid-cell

Missing data codes: No missing data.

MD5 checksum: 4527b1965fdbb732f341485c4cbf70b7

--------------------------

japan_ferns_traits_lucid.csv: Morphological trait data of native, non-hybrid
ferns in Japan, originally formatted for Lucid dichotomous key software
(https://www.lucidcentral.org/). For more information on codes used for scoring
of traits, see https://help.lucidcentral.org/lucid/scoring-the-key/

Number of variables: 436

Number of cases/rows: 675

Variable list: ‘taxon’ indicates taxon name. Other columns correspond to traits.
If the column name includes ‘number’ it is numeric, otherwise categorical.

Categorical trait codes as follows:

•	0 = absent

•	1 = common

•	2 = rare

•	3 = uncertain

•	4 = common and misinterpreted

•	5 = rare and misinterpreted

•	6 = not scoped

Numeric trait data formatted as a series of numbers separated by colons (‘:’),
unless missing (0 or 3). For example, 1:9.4:12:14:16. The starting ‘1:’ does not
mean anything. The next series of numbers correspond to the outside minimum,
normal minimum, normal maximum and outside maximum. For example, a plant leaf
may normally be in the range of 4-5 cm, but sometimes as low as 2 cm or as high
as 6 cm. This would be written 2:4:5:6.

Missing data codes: ‘0’ indicates trait is missing; ‘3’ indicates uncertain.

MD5 checksum: f0fdb678e38fa19a650140af9d4c1d06

--------------------------

japan_map_points.csv: Latitude and longitude of points used to produce a labeled
map of Japan. Some points were originally generated using script at
https://github.com/joelnitta/japan_ferns_spatial_phy/blob/main/R/geocode.R

Number of variables: 5

Number of cases/rows: 17

Variable list:

•	query: Name of point of interest

•	lat: Latitude

•	lon: Longitude

•	lat_end: Ending latitude (only for features that are lines, not points)

•	lon_end: Ending longitude (only for features that are lines, not points)

Missing data codes: No data entered between commas in CSV.

MD5 checksum: 17e92adffd7daa40b26494ae680fd4cb

--------------------------

japan_map.gpkg: Map of land area of Japan. Raw data downloaded from Geospatial
Information Authority of Japan under the Creative Commons Attribution License
v4.0 (https://www.gsi.go.jp/kankyochiri/gm_japan_e.html). All polygons were
combined into a single polygon and CRS set to JGD2000 using code available at
https://github.com/joelnitta/japan_ferns_spatial_phy/blob/main/R/process_raw_data.R.

Number of variables: 1

Number of cases/rows: 1

Variable list:

•	geom: Vector describing shape and position of land area polygon

Missing data codes: No missing data.

MD5 checksum: a44b1953457aa44bbf5cd67bbacd0bd0

--------------------------

japan_protected_areas.gpkg: Spatial of protected areas in Japan in GeoPackage
format. Original spatial data downloaded from the Japan Ministry of the
Environment
(https://www.biodic.go.jp/biodiversity/activity/policy/map/map17/index.html) and
Ministry of Land, Infrastructure, Transport and Tourism
(https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-A45.html). Categorization of
protected areas generally followed Kusumoto et al. (2017): either “high” (no
human activities allowed at all) or “medium” status (some economic activities
allowed by permit); other areas not meeting these criteria or not including
plants were not considered. Any overlapping areas within the same protection
level were combined, and any areas overlapping between “medium” and “high” were
considered only “high”. Code used to categorize areas available at
https://github.com/joelnitta/japan_ferns_spatial_phy/blob/main/R/process_raw_data.R

Number of variables: 2

Number of cases/rows: 2

Coordinate reference system: JGD2000

Variable list:

•	status: Protection status (‘high’ or ‘medium’)

•	geom: Vector describing shape and position of protected area polygon

Missing data codes: No missing data.

MD5 checksum: 0c2c7959b0c7cc7a837428ece390e5b8

--------------------------

testo_sundue_2016_calibrations.csv: Fossil calibration points used for
phylogenetic dating analysis with treePL (Smith and O’Meara 2012), based on
Testo and Sundue (2016). MRCA = “most recent common ancestor”.

Number of variables: 9

Number of cases/rows: 27

Variable list:

•	Clade: Name of clade

•	Stem_Crown: Type of node. ‘crown’ (node corresponding to MRCA of extant taxa
only) or ‘stem’ (node corresponding to split between MRCA of extant taxa and its
sister)

•	Fossil: Name of fossil

•	Age: Age of fossil

•	Age_type: Type of age. ‘min’ for minimum age, or ‘fixed’ for fixed age

•	Citation: Citation for the fossil

•	taxon_1: Name of first taxon that is used to identify MRCA

•	taxon_2: Name of second taxon that is used to identify MRCA

•	note: Notes, in particular explaining any differences with Testo and Sundue
(2016)

Missing data codes: No data entered between commas in CSV.

MD5 checksum: 60a5bf932a5b9039b5caf61567d172d0

--------------------------

japan_ferns_biodiv.csv (contained in “results.zip”): Biodiversity statistics of
native, non-hybrid ferns and environmental variables in Japan. Biodiversity
metrics calculated as described in Nitta et al. 2021. Climatic (temperature and
preciptation) variables calculated as described for japan_climate.gpkg. Includes
one row with missing environmental data and one outlier for % apomixis that were
removed prior to spatial modeling analysis in Nitta et al. 2021.

Number of variables: 28

Number of cases/rows: 1241

Variable list:

•	grids: Name of grid-cell

•	lat: Latitude of grid-cell centroid

•	long: Longitude of grid-cell centroid

•	abundance: Number of specimens occurring in grid-cell

•	richness: Number of taxa occurring in grid-cell

•	fd_obs_z: Z-score (standard effect size) of observed functional diversity (FD)

•	pd_obs_z: Z-score of observed phylogenetic diversity (PD)

•	pe_obs_z: Z-score of observed phylogenetic endemisim (PE)

•	rfd_obs_z: Z-score of observed relative FD

•	rpd_obs_z: Z-score of observed relative PD

•	richness_obs_p_upper: Richness percentile

•	fd_obs_p_upper: FD percentile

•	pd_obs_p_upper: PD percentile

•	pe_obs_p_upper: PE percentile

•	pd_signif: Significance of PD compared to null distribution, two-sided test at
alpha 0.05

•	rpd_signif: Significance of RPD compared to null distribution, two-sided test
at alpha 0.05

•	fd_signif: Significance of FD compared to null distribution, two-sided test at
alpha 0.05

•	rfd_signif: Significance of RFD compared to null distribution, two-sided test
at alpha 0.05

•	pe_signif: Significance of PE compared to null distribution, two-sided test at
alpha 0.05

•	taxonomic_cluster: Membership in cluster based on taxonomic distances

•	phylo_cluster: Membership in cluster based on phylogenetic distances

•	endem_type: Endemism type determined by CANAPE

•	lat_area: Area of 1-degree latitudinal band around grid-cell

•	temp: Annual mean temperature, in units of degrees celsius * 10

•	temp_season: Temperature seasonality

•	precip: Annual precipitation, measured in mm

•	precip_season: Precipitation seasonality

•	percent_apo: Percent apomictic taxa

Missing data codes: No data entered between commas in CSV.

MD5 checksum: bdb5c51884fe831049533e4497089ad0

--------------------------

japan_ferns_comm.csv (contained in “results.zip”): Community matrix (species x
sites matrix) of native, non-hybrid ferns in Japan used for biodiversity
analysis. Same as japan_ferns_comm_full.csv, but grid-cells with redundancy <0.1
(indicating under-sampling) excluded.

Number of variables: 673

Number of cases/rows: 1241

Variable list:

•	grids: Name of grid-cell.

•	Other columns: Each column is named for a Japanese fern taxon. Values indicate
number of specimens that were observed in each grid-cell for that taxon.

Missing data codes: No missing data.

MD5 checksum: 7e3320123ddb5edd4c7b7e7ad83ee578

--------------------------

japan_ferns_shape.gpkg (contained in “results.zip”): Location of grid-cells
(sites) for native, non-hybrid ferns in Japan used for biodiversity analysis.
Same as japan_ferns_shape_full.gpkg, but grid-cells with redundancy <0.1
(indicating under-sampling) excluded.

Number of variables: 5

Number of cases/rows: 1241

Coordinate reference system: JGD2000

Variable list:

•	grids: Name of grid-cell

•	abundance: Number of specimens occurring in grid-cell

•	richness: Number of taxa occurring in grid-cell

•	redundancy: Redundancy (measure of sampling completeness)

•	geom: Vector describing shape and position of grid-cell

Missing data codes: No missing data.

MD5 checksum: 59008832d23aafd65f59232d87155e92

--------------------------

japan_ferns_traits.csv (contained in “results.zip”): Trait matrix of native,
non-hybrid ferns in Japan used for functional biodiversity analysis.

Number of variables: 78

Number of cases/rows: 675

Variable list: ‘taxon’ indicates taxon name. Other columns correspond to traits.
‘frond_width’, ‘stipe_length’, ‘number_pinna_pairs’ are numeric; others are
binary (0 or 1)

Missing data codes: Missing data coded as ‘NA’.

MD5 checksum: 58cd40120c19b27bea373bdf51288e70

--------------------------

japan_ferns_tree.tre (contained in “results.zip”): Maximim-likelihood phylogeny
of native, non-hybrid ferns in Japan in newick format inferred with maximum
likelihood IQ-TREE v1.6.12 (Nguyen et al. 2015). Values at nodes are bootstrap
support values calculated using 1000 replicates of ultrafast bootstrap (Nguyen
et al. 2015).

Number of tips: 663

Missing data codes: No missing data.

MD5 checksum: f5905d6c788238b3d1e01a5b75af8c39

--------------------------

japan_ferns_tree_dated_uncollapsed.tre (contained in “results.zip”):
Maximim-likelihood, ultrametric (dated) phylogeny of native, non-hybrid ferns in
Japan in newick format. Tree dated using treePL v1.0 (Smith and O’Meara 2012)
with 26 fossil calibration points after Testo and Sundue (2016). Values at nodes
are bootstrap support values calculated using 1000 replicates of ultrafast
bootstrap (Nguyen et al. 2015). Clades consisting of identical OTUs not
collapsed.

Number of tips: 663

Missing data codes: No missing data.

MD5 checksum: feaf9aa5ff2d55ef80ca830b10d59dcc

--------------------------

japan_ferns_tree_dated.tre (contained in “results.zip”): Maximim-likelihood,
ultrametric (dated) phylogeny of native, non-hybrid ferns in Japan in newick
format. Tree dated using treePL v1.0 (Smith and O’Meara 2012) with 26 fossil
calibration points after Testo and Sundue (2016). Values at nodes are bootstrap
support values calculated using 1000 replicates of ultrafast bootstrap (Nguyen
et al. 2015). Clades consisting of identical sequences have been collapsed to
polytomies with zero branch length between OTUs.

Number of tips: 663

Missing data codes: No missing data.

MD5 checksum: 096e6835cc7160e2968d7598f3bcef20

--------------------------

CHANGE LOG

---

2021-09-09

Generate this README file.

--------------------------

References

Ebihara, Atsushi, and Joel H. Nitta. 2019. “An Update and Reassessment of Fern
and Lycophyte Diversity Data in the Japanese Archipelago.” Journal of Plant
Research 132 (6): 723–38. https://doi.org/10.1007/s10265-019-01137-3.

Hsieh, T. C., K. H. Ma, and Anne Chao. 2016. “iNEXT: An R Package for
Rarefaction and Extrapolation of Species Diversity (Hill Numbers).” Methods in
Ecology and Evolution 7 (12): 1451–6. https://doi.org/10/gdskcx.

Kusumoto, Buntarou, Takayuki Shiono, Masashi Konoshima, Atsushi Yoshimoto,
Takayuki Tanaka, and Yasuhiro Kubota. 2017. “How Well Are Biodiversity Drivers
Reflected in Protected Areas? A Representativeness Assessment of the
Geohistorical Gradients That Shaped Endemic Flora in Japan.” Ecological Research
32 (3): 299–311. https://doi.org/10.1007/s11284-017-1451-6.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh.
2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
Maximum-Likelihood Phylogenies.” Molecular Biology and Evolution 32 (1): 268–74.
https://doi.org/10.1093/molbev/msu300.

Smith, Stephen A, and Brian C. O’Meara. 2012. “treePL: Divergence Time
Estimation Using Penalized Likelihood for Large Phylogenies.” Bioinformatics 28
(20): 2689–90. https://doi.org/10.1093/bioinformatics/bts492.

Testo, Weston L., and Michael A Sundue. 2016. “A 4000-Species Dataset Provides
New Insight into the Evolution of Ferns.” Molecular Phylogenetics and Evolution
105: 200–211. https://doi.org/10.1016/j.ympev.2016.09.003. 
