# Spatial phylogenetics of Japanese ferns

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "Spatial phylogenetics of Japanese ferns: Patterns, processes, and implications for conservation". https://doi.org/10.1101/2021.08.26.457744 

All code is in [R](https://cran.r-project.org/). The [targets package](https://docs.ropensci.org/targets/index.html) is used to manage the workflow. To run all analyses and generate the manuscript, [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) and run `targets::tar_make()`.

## Data

FIXME: Add a description of how to download data

For more information about data files, see the READMEs for [raw data](doc/README_data_raw.md) and [processed data](README_data.txt).

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

## Using the `targets` cache

The analysis includes some steps that take a long time to run, especially maximum-likelihood phylogenetic analysis (ca. 1 week with 10 cores in parallel). To avoid running the entire workflow from scratch, untar the `_targets.tar.gz` file in the Dryad dataset and place it in the root of this repo as `_targets`:

```
tar -xzf _targets.tar.gz
```

Then, when you open the project in R [as described above](#interacting-with-the-code), you can use `targets::tar_load()` to load any target (intermediate workflow step) listed in [`_targets.R`](_targets.R). For more information on how to use the `targets` package, see https://github.com/ropensci/targets.

## Licenses

- Code: [MIT license](LICENSE.md)
- Data: [CC0 1.0 license](https://creativecommons.org/publicdomain/zero/1.0/)
- [Manuscript (preprint)](https://doi.org/10.1101/2021.08.26.457744): [CC BY-NC-ND 4.0 license](https://creativecommons.org/licenses/by-nc-nd/4.0/)
- [Roboto font](https://github.com/google/roboto/): [Apache 2.0 license](http://www.apache.org/licenses/LICENSE-2.0)

