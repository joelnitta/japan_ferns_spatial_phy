# Data

Data files used in analysis.

2_grid_cells_all.csv: List of 10km grid cells for Japan ("secondary
grid cells"). Columns include: "id", secondary grid cell code; 
"name", name of grid cell in Japanese, "x", longitude, "y", latitude. 
Secondary grid code refers to "Basic Grid Square Codes" defined by the 
statistics Bureau of Japan (http://www.stat.go.jp/english/data/mesh/05.html).
Encoding is Unicode (UTF-8).

ESM1.csv: A list of native fern and lycophyte taxa (species, subspecies and varieties) 
in Japan accepted in this study. Taxon ID refers to that in FernGreenList ver.1.0.1 (http://www.rdplants.org/gl/). Unless otherwise noted, rbcL GenBank accession numbers 
are those used in Ebihara et al. (2010). Asterisks indicate newly generated sequences 
by this study. Voucher information only provided for newly generated sequences. 
Information on reproductive modes and ploidy levels follow those in Ebihara et 
al. (2016, 2017), and only records based on material collected in Japan are used. 
For reproductive mode, irregular meiosis is not considered, 0 = no information, 
1 = sexual, 2 = apomictic and 3 = sexual + apomictic. Encoding is Unicode (UTF-8).

ESM2.csv: A list of fern and lycophyte herbarium specimens from Japan used to
generate the 10 km2 grid cell distribution maps in Ebihara et al. (2016, 2017). 
A single specimen is cited per taxon per grid cell for all native taxa including 
species, subspecies, varieties, and hybrids. Taxon ID refers to that in 
FernGreenList ver. 1.0.1 (http://www.rdplants.org/gl/). Secondary grid code as
in 2_grid_cells_all190705.xlsx. Herbarium acronyms follow Index Herbariorum (http://sweetgum.nybg.org/science/ih/) with additional acronyms defined in the
manuscript. Encoding is Unicode (UTF-8).

FernGreenListV1.01.xls: List of Japanese ferns and lycophytes species including
scientific name, endemic status, conservation status, and other taxonomic data. 
Downloaded from http://www.rdplants.org/gl/FernGreenListV1.01.xls on 2019-07-17.

japan_pterido_rbcl_cipres.zip: Output of Bayesian analysis of Japanese pteridophyte
rbcL data matrix (rbcl_mrbayes.nex) downloaded from the CIPRES cluster. The resulting
tree file is "infile.nex.con.tre" within the zip folder. Tips named with two numbers
separated by an underscore, e.g., 601_1. The first number corresponds to family code,
and the second number to taxon ID in FernGreenList ver. 1.0.1 and ESM1.xlsx.

JpFern_rbcl.txt: DNA alignment in PHYLIP format of Japanese ferns and lycophytes.
Sequence names formatted the same as tip labels in the tree file mentioned above.
This is used by "make_mrbayes_nex.R" in the "code" folder to make rbcl_mrbayes.nex.

ppgi_taxonomy.csv: Spreadsheet of the Pteridophyte Phylogeny I working group 
taxonomic system for pteridophytes at the genus level and above (The Pteridophyte 
Phylogeny Group, 2016. A community-derived classification for extant lycophytes 
and ferns. J Syst Evol 54:563-606). Includes columns for class, order, suborder, 
family, subfamily, and genus. Encoding is Unicode (UTF-8).

rbcl_mrbayes.nex: NEXUS file used for phylogenetic analysis of Japanese fern
and lycophyte taxa with MrBayes. Output by "make_mrbayes_nex.R" in the "code"
folder. This was uploaded to CIPRES for phylogenetic analysis.
