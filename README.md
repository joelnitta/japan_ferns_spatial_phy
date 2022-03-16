# Spatial phylogenetics of Japanese ferns

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "Spatial phylogenetics of Japanese ferns: Patterns, processes, and implications for conservation". https://doi.org/10.1101/2021.08.26.457744 

All code is in [R](https://cran.r-project.org/). The [targets package](https://docs.ropensci.org/targets/index.html) is used to manage the workflow. To run all analyses and generate the manuscript, [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) and run `targets::tar_make()`.

## Data

Data files need to be downloaded from three locations.

1. Dataset on FigShare for this project: https://doi.org/10.6084/m9.figshare.16655263. Cick on the "Download all" icon, download the zipped dataset, then unzip it and put the contents in the `data/` folder in this repo.
2. Dataset on Dryad for Ebihara and Nitta 2019: https://datadryad.org/stash/dataset/doi:10.5061/dryad.4362p32. Download the zipped dataset and put in the `data/` folder directly (without unzipping).
3. Dataset on FigShare for FTOL v0.0.1 (Nitta et al, in prep): https://doi.org/doi:10.6084/m9.figshare.13256801 (LINK NOT LIVE YET). Download the zipped dataset and put in the `data/` folder directly (without unzipping).

For more information about data files, see the READMEs for [raw data](doc/README_data_raw.md) and [processed data](doc/README_data.txt).

## Reproducible analysis with Docker

This project requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/japan_ferns_spatial_phy) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/) and clone this repository.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine):

```
cd /path/to/repo
```

Run `targets::tar_make()`:

```
docker run --rm -v ${PWD}:/tmpdir -w /tmpdir joelnitta/japan_ferns_spatial_phy:latest Rscript -e 'targets::tar_make()'
```

You will see the targets being built by `targets`, and the final manuscript should be compiled at the end as `manuscript.pdf` and `manuscript.docx` in the `results` folder. Other figure and table files will also be compiled.

## Interacting with the code

If you want to interact with the code in the Docker container, you can launch the container in the background using `docker-compose`:

```
docker-compose up -d
```

Navigate to http://localhost:8787/ in your browser of choice (firefox or google chrome recommended). There, you should be able to access an instance of the [RStudio](https://rstudio.com/) IDE, which can be used to inspect and manipulate objects in R. You can click on "Build All" in the "Build" tab to run the workflow. 

When you're done, take down the container:

```
docker-compose down
```

## Targets cache

The [targets package](https://docs.ropensci.org/targets/index.html) manages the workflow and saves all intermediate analysis results to a folder named `_targets`; this is the targets cache.
Normally, you would have to run all of the analyses starting from the original data files to generate all of the analysis results, as described above.
This takes a long time. The longest step is the phylogenetic analysis, which takes about 1 week using 10 cores in parallel.

I have put the targets cache for this project [on github](https://github.com/joelnitta/japan_ferns_spatial_phy_cache) under version control using the [gittargets package](https://github.com/ropensci/gittargets).

So instead of running everything from scratch, you can checkout the exact results matching a specific code version as follows (this assumes we are in the `japan_ferns_spatial_phy` folder and requires git):

1. Clone the targets cache to a folder called `_targets`.

```
git clone https://github.com/joelnitta/japan_ferns_spatial_phy_cache _targets
```

2. Enter the `_targets` directory.

```
cd _targets
```

3. Fetch branches from the remote repo ([each branch corresponds to a selected commit in the code](https://docs.ropensci.org/gittargets/articles/git.html#snapshot-model)).

```
git fetch
```

4. Change to the latest branch (the part of the name after `code=` matches the corresponding commit hash in `japan_ferns_spatial_phy`).

```
git switch code=868d97bc3e205adbf417c74123314f48db87e368
```

5. Move back up to the `japan_ferns_spatial_phy` folder.

```
cd ..
```

You can also change between different snapshots of the targets cache and code using [gittargets](https://github.com/ropensci/gittargets).

When you open the project in R [as described above](#interacting-with-the-code), you can use `targets::tar_load()` to load any target (intermediate workflow step) listed in [`_targets.R`](_targets.R). For more information on how to use the `targets` package, see https://github.com/ropensci/targets.

## Licenses

- Code: [MIT license](LICENSE.md)
- Data: [CC0 1.0 license](https://creativecommons.org/publicdomain/zero/1.0/)
- [Manuscript (preprint)](https://doi.org/10.1101/2021.08.26.457744): [CC BY-NC-ND 4.0 license](https://creativecommons.org/licenses/by-nc-nd/4.0/)
- [Roboto font](https://github.com/google/roboto/): [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)
