## Setting-up a conda environment
> Only for advanced users.

This documentation contains installing instructions for advanced users that opt to setup a conda environment for the transcriptutorial.

### Install miniconda3
> NOTE: only done once (1st time) in your local machine

`miniconda3` is the framework to create and manage `envs`.
Before starting to work with virtual environments (`envs`), you need to install `miniconda3`.

Download the installer of `miniconda3`
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Create a directory where you are going to install your SOFTWARE
```
mkdir /home/$USER/SOFTWARE/
```

Run the installer and define the installation path. See below.
```
bash Miniconda3-latest-Linux-x86_64.sh

# Press <ENTER> to continue an read Miniconda License
# Go down and enter 'yes' to accept Miniconda License
# Importantly! Miniconda3 will be installed into the default location,
# that is your $HOME. However, you could install it in another location
# Thus, specify another location. For instance:
#	/home/<username>/SOFTWARE/miniconda3
# run conda init when asked during the installation
```
If you didnt run `conda init` during the installation, this is a good time to do so. This will make changes into your `~/.bashrc` file needed to setup the conda environment when needed.

Remove the installer
```
rm Miniconda3-latest-Linux-x86_64.sh
```

### Setup of environment

The first time that you install conda, it is installed with a default configuration.
You have to extend this configuration by adding essential channels that are required to install `R`, `CRAN` and `BioConductor` packages.

```bash
conda config --append channels conda-forge
conda config --append channels r
conda config --append channels bioconda
```

Then you can create the conda enviroment from scratch using the command below. Alternatively, you could create it from a recipe (see next chunk).
```bash 
conda create -p envs/transcriptutorial r-base=4.0 python=3.8 bioconductor-CARNIVAL r-cowplot bioconductor-dorothea r-dplyr \ 
	r-ggplot2 r-ggrepel r-gridExtra bioconductor-GSEABase r-hexbin bioconductor-limma r-network bioconductor-OmnipathR \
	r-pheatmap bioconductor-piano r-plyr bioconductor-progeny r-readr r-reshape r-reshape2 r-scales r-tibble r-tidyr \
	r-visNetwork bioconductor-vsn r-rmarkdown
```

[Or create from a recipe]
```bash
conda create --prefix ./envs/transcriptutorial --file ./envs/transcriptutorial.txt
```

Finally, get the latest `devel` version of `OmnipathR`:
```bash
# Terminal
conda activate ./envs/transcriptutorial
$CONDA_PREFIX/bin/R
# R session
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## Last release in Bioconductor
BiocManager::install("OmnipathR", version='devel')
```
### External devtools::packages
In case of external packages to CRAN/Bioconductor, these could be installed
in an interactive R session using `devtools`.

```bash
# Terminal
conda activate ./envs/transcriptutorial
conda install -p ./envs/transcriptutorial r-devtools
$CONDA_PREFIX/bin/R
```

```r
# R session
# Setup for devtools installation
getOption("unzip")
Sys.getenv("TAR")
options(unzip = "/opt/conda/bin/unzip")
Sys.setenv(TAR = "/bin/tar")
# Actual required packages
devtools::install_github("USER/REPOSITORY") 
```
