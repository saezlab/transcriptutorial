## 1. Clone repository
You will get file structure of the project and main scripts by using:
```
git clone https://github.com/saezlab/transcriptutorial.git
```

## 2. Download data
The data can be downloaded here : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119931
At the bottom of the webpage, download GSE119931_PANC1.FOXA2KO.genes.counts.txt.gz and uncompress it in the data folder.
Do it by using:
```bash
url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119931/suppl/GSE119931_PANC1.FOXA2KO.genes.counts.txt.gz";
destfile="./data/$(basename $url)";

wget -O $destfile $url
gunzip $destfile
```

## 3. Installation of R packages
> NOTE: For advanced users, if you prefer to setup a conda environment instead, visit [setup_conda_env.md](setup_conda_env.md).

First, you need an updated version of `R` (>4.0) (http://www.r-project.org).
Then you can install all required packages following the commands on `R` console:

```R
# Define required packages
cran_pkgs <- c("cowplot", "dplyr", "ggplot2", "ggrepel", "gridExtra", "hexbin", "network", 
"pheatmap", "plyr", "readr", "reshape", "reshape2", "scales", "tibble", "tidyr", 
"ggraph", "tidygraph", "snowfall", "visNetwork", "rmarkdown")

bioC_pkgs <- c("viper", "dorothea", "GSEABase", "limma", "OmnipathR", "piano", "progeny", 
"vsn", "CARNIVAL")

# Install CRAN packages
install.packages(cran_pkgs)

# Install BioConductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioC_pkgs)

# [optional] Install the devel version of omnipathR
BiocManager::install("OmnipathR", version='devel')
```

