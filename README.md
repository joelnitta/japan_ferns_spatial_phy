# Biodiversity of Japanese Ferns and Lycophytes

## Building the analysis environment

`conda` is used to maintain a reproducible analysis environment, and `packrat` is used to maintain R package versions. To create the environment, first install conda, or update to the most recent version:

```
conda update -n base -c defaults conda
```

Build the environment from `environment.yml`:

```
conda env create -f environment.yml
```

Load the environment:

```
source activate japan-ferns-biogeo
```
