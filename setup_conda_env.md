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

Then you can create the conda enviroment from scratch using: 
```bash 
conda create -p envs/transcriptutorial r-base=4.0 python=3.8 bioconductor-CARNIVAL r-cowplot r-dplyr bioconductor-viper \ 
	r-ggplot2 r-ggrepel r-gridExtra bioconductor-GSEABase r-hexbin bioconductor-limma r-network bioconductor-OmnipathR \
	r-pheatmap bioconductor-piano r-plyr bioconductor-progeny r-readr r-reshape r-reshape2 r-scales r-tibble r-tidyr \
	r-visNetwork bioconductor-vsn r-ggraph r-tidygraph r-rmarkdown
```

Finally, install manually `dorothea` and the latest `devel` version of `OmnipathR`:
```bash
# Terminal
conda activate ./envs/transcriptutorial
$CONDA_PREFIX/bin/R
# R session
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## dorothea
#NOTE: dorothea is included in bioconda, although installation fails due to wrong indexing to bioconductor URL 
#	with current version of recipe. We have reported the issue. It will be fixed soon.
BiocManager::install("dorothea")

## OmnipathR package
BiocManager::install("OmnipathR", version='devel')
```

> Note: if the installation of any of the packages fails with the conda recipe, just install it manually as the last chunk above.

