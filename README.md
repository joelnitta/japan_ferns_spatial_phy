# Spatial phylogenetics of Japanese ferns

Code repostitory to run analyses and generate figures and manuscript for Nitta et al. "Spatial phylogenetics of Japanese ferns: Patterns, processes, and implications for conservation". https://doi.org/10.1101/2021.08.26.457744 

All code is in [R](https://cran.r-project.org/). The [targets package](https://docs.ropensci.org/targets/index.html) is used to manage the workflow. To run all analyses and generate the manuscript, [clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository) and run `targets::tar_make()`.

## Data

FIXME: Add a description of how to download data

## Reproducible analysis with Docker

This project requires various packages to be installed, and may not work properly if package versions have changed. Therefore, a [Docker image is provided](https://hub.docker.com/r/joelnitta/japan_ferns_spatial_phy) to run the code reproducibly.

To use it, first [install docker](https://docs.docker.com/install/) and clone this repository.

Navigate to the cloned repository (where `/path/to/repo` is the path on your machine), and launch the container:

```
cd /path/to/repo
docker-compose up -d
```

Enter the container:

```
docker exec -it japan_ferns_analysis_1 bash
```

Inside the container, run `targets::tar_make()`:

```
Rscript -e "targets::tar_make()"
```

You will see the targets being built by `targets`, and the final manuscript should be compiled at the end as `manuscript.pdf` and `manuscript.docx` in the `results` folder. Other figure and table files will also be compiled.

When it's finished, exit the container and take it down:

```
exit
docker-compose down
```

