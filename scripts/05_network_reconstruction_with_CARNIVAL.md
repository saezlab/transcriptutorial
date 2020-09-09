05\_network\_reconstruction\_with\_CARNIVAL
================
Rosa Hernansaiz-Ballesteros
18/05/2020

### License Info

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

Please check <http://www.gnu.org/licenses/>.

## Introduction

**CARNIVAL** (CAusal Reasoning for Network identification using Integer
VALue programming) is a method for the identification of upstream
regulatory signalling pathways from downstream gene expression (GEX)
[(Liu, Trairatphisan, Gjerga et
al. 2019)](https://doi.org/10.1038/s41540-019-0118-z). The aim of
**CARNIVAL** is to identify a subset of interactions from a prior
knowledge network, a network that represents potential regulated
pathways linking known or potential targets of perturbation, to active
transcription factors derived from GEX data.

In the fifth part of our transcriptomics tutorial series, we demonstrate
how to use **CARNIVAL** based on the transcription factor (TF) activity
derived from transcriptomics data using **DoRothEA** [(Garcia-Alonso et
al. 2019)](https://doi.org/10.1101/gr.240663.118), and the prior
knowledge network obtained from [**Omnipath**](http://omnipathdb.org/),
a literature curated mammalian signaling pathway resource [(Türei et al
2016)](https://www.nature.com/articles/nmeth.4077).

In order to help CARNIVAL to find a solution faster, we can also use
**PROGENy** [(Schubert et
al. 2018)](https://www.nature.com/articles/s41467-017-02391-6) scores
to infer the score of the representative genes for each of the
calculated
pathways.

[**CARNIVAL**](https://www.bioconductor.org/packages/release/bioc/html/CARNIVAL.html),
[**DoRothEA**](http://bioconductor.org/packages/release/data/experiment/html/dorothea.html),
[**PROGENy**](http://bioconductor.org/packages/release/bioc/html/progeny.html)
and
[**OmnipathR**](https://www.bioconductor.org/packages/release/bioc/html/OmnipathR.html)
are available as bioconductor packages. Visit our website for additional
information about each of the tools. <https://saezlab.org>

## Getting Started

We first load the required libraries and the support functions.

``` r
library(progeny)
library(dorothea)
library(CARNIVAL)
library(OmnipathR)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(visNetwork)

library(ggplot2)
library(pheatmap)


## For the volcano plot (related to support functions)
library(ggrepel)

## We also load the support functions
source("assignPROGENyScores.r")
source("generateTFList.r")
source("carnival_visNetwork.r")
```

In addition, we read the results from the previous scripts:

  - Transcription factor activities
    (04\_TranscriptionFactor\_activity\_with\_Dorothea.Rmd)
  - Pathways activity scores (03\_Pathway\_activity\_with\_Progeny.Rmd).

<!-- end list -->

``` r
## We read the normalised counts and the experimental design 
tf_activities <- read_csv("../results/TFActivity_CARNIVALinput.csv")
PathwayActivity <- read_csv("../results/PathwayActivity_CARNIVALinput.csv")
```

## Getting the scaffold network from Omnipath

Before running **CARNIVAL**, we need to create or upload a scaffold
network. This will be *“the map”* that the ILP algorithm will follow to
find the causal network. We use **Omnipath** to obtain the signed and
directed interactions from all the available resources. CARNIVAL
requires this information in a *sif* table (node1, interaction, node2)
format, therefore we use the *consensus* columns of direction
(consensus\_direction) and sign (consensus\_stimulation and
consensus\_inhibition) to extract it.

The query returns 0/1 as logic status of being a stimulation or an
inhibition reaction. Thus, this output is reformulated as 1/-1 to
indicate stimulation or inhibition, respectively. We can keep either the
interactions that are consistent, or both alternatives (e.g. A 1 B; A -1
B). In this example, we keep the consistent ones.

``` r
omniR <- import_Omnipath_Interactions()
```

    ## Warning: 'import_Omnipath_Interactions' is deprecated.
    ## Use 'import_omnipath_interactions' instead.
    ## See help("Deprecated")

``` r
# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
                                (consensus_stimulation == 1 | 
                                 consensus_inhibition == 1
                                 ))
  
# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which( omnipath_sd$consensus_stimulation == 0)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 1)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 0)] = 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
      dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
      unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c('source', 'interaction', 'target')

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)

#save SIF
write_tsv(sif, "../results/omnipath_carnival.tsv")
```

## Transcription Factor and pathway activities for CARNIVAL

We use the supplementary functions *generateTFList.r* and
*assignPROGENyScores.r* to shift the formats of tf\_activities and
PathwayActivity to the one required by CARNIVAL.

``` r
# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$TF
tf_activities_carnival$TF <- NULL
tfList = generateTFList(tf_activities_carnival, top=50, access_idx = 1)

# progeny for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
PathwayActivity_carnival$Pathway <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                            progenyMembers = progenyMembers, 
                                            id = "gene", 
                                            access_idx = 1)
```

## Running CARNIVAL

CARNIVAL has been developed to find the causal link between the
activities of the transcription factors (TFs) and the ‘perturbed’ nodes.
In current version, v1.0.0, we have 3 main inputs that we have to
provide:

  - *measObj*: The TFs’ activities (like the ones we have obtained from
    DoRothEA)
  - *inputObj*: The ‘perturbed’ nodes we want that CARNIVAL connects
    with the activity of TFs. There are 3 ways of using it:

<!-- end list -->

1)  Give the name and sign of the selected nodes;
2)  Give the name only, so the algorithm will select the sign that best
    fit the models,
3)  Give *NULL* as value will create a “Perturbation” node that will try
    both signs for all ‘initial’ nodes of the given network ( *netObj*
    ).

<!-- end list -->

  - *netObj*: The network that will serve as map to connect the TFs’
    activities ( *measObj* ) and the perturbed nodes ( *inputObj* )

Although it is not required, a fourth object called *weightObj* can be
also given. This object gives values ranged from -1 to 1 for a set of
nodes of the network. The aim of *weightObj* is helping the solver to
find optimal solutions faster.

In the present example, we use assign as perturbation nodes all the
“initial” nodes (option 2), and as *weightObj* the PROGENy scores
assigned to the most representative genes of the calculated pathways,

Please, check the [CARNIVAL](https://saezlab.github.io/CARNIVAL/) page
to get some more insight.

``` r
# get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

# run carnival
carnival_result = runCARNIVAL( inputObj= iniciators,
                               measObj = tfList$t, 
                               netObj = sif, 
                               weightObj = progenylist$score, 
                               solverPath = "/Applications/CPLEX_Studio129/cplex/bin/x86-64_osx/cplex", 
                               solver = "cplex",
                               timelimit=7200,
                               mipGAP=0,
                               poolrelGAP=0 )
```

    ## Warning in controlNodeIdentifiers(netObj = netObj): Your network contains identifiers with '-' symbol and they will
    ##             be replaced with '_'

    ## Warning in controlNodeIdentifiers(netObj = netObj): Your network contains identifiers with '/' symbol and they will
    ##             be replaced with '_'

CARNIVAL gives a list of 4 elements:

  - weightedSIF: summary of all interactions found in all models
  - nodesAttributes: summary of all nodes and how many times are they
    found in the different models
  - sifAll: networks of all the models
  - attributesAll: node attributes of all models

We can now visualise the network…

``` r
#transoform to data.frame
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)

carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)

saveRDS(carnival_result,"../results/carnival_result.rds")

# visualization
visNet = carnival_visNet(evis = carnival_result$weightedSIF,
                         nvis = carnival_result$nodesAttributes)
```

    ## Graphical representation of sample

``` r
#visNet
visSave(visNet, file = paste0('carnival_visualization_visNetwork.html'), selfcontained = TRUE)
```

## References

> Liu A., Trairatphisan P., Gjerga E. et al. “From expression footprints
> to causal pathways: contextualizing large signaling networks with
> CARNIVAL”. *npj Systems Biology and Applications*. 2019. DOI:
> [10.1038/s41540-019-0118-z](https://www.nature.com/articles/s41540-019-0118-z)

> Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J.
> “Benchmark and integration of resources for the estimation of human
> transcription factor activities.” *Genome Research*. 2019. DOI:
> [10.1101/gr.240663.118](https://genome.cshlp.org/content/29/8/1363).

> Türei D, Korcsmáros T, & Saez-Rodriguez J, “OmniPath: guidelines and
> gateway for literature-curated signaling pathway resources”. *Nat
> Methods* 2016. DOI:
> [10.1038/nmeth.4077](https://doi.org/10.1038/nmeth.4077).

> Schubert M, Klinger B, Klünemann M, Sieber A, Uhlitz F, Sauer S,
> Garnett MJ, Blüthgen N, Saez-Rodriguez J. “Perturbation-response genes
> reveal signaling footprints in cancer gene expression.” *Nature
> Communications*. 2018. DOI:
> [10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6)

## Session Info Details

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggrepel_0.8.2    pheatmap_1.0.12  ggplot2_3.3.2    visNetwork_2.0.9
    ##  [5] dplyr_1.0.2      tidyr_1.1.1      tibble_3.0.3     readr_1.3.1     
    ##  [9] OmnipathR_1.3.5  jsonlite_1.7.0   igraph_1.2.5     CARNIVAL_1.0.0  
    ## [13] dorothea_1.0.0   progeny_1.10.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2           Biobase_2.48.0       mixtools_1.2.0      
    ##  [4] UniProt.ws_2.28.0    bit64_4.0.2          splines_4.0.2       
    ##  [7] foreach_1.5.0        assertthat_0.2.1     BiocFileCache_1.12.1
    ## [10] stats4_4.0.2         RBGL_1.64.0          blob_1.2.1          
    ## [13] yaml_2.2.1           Category_2.54.0      viper_1.22.0        
    ## [16] pillar_1.4.6         RSQLite_2.2.0        lattice_0.20-41     
    ## [19] glue_1.4.1           digest_0.6.25        RColorBrewer_1.1-2  
    ## [22] colorspace_1.4-1     htmltools_0.5.0      Matrix_1.2-18       
    ## [25] GSEABase_1.50.1      lpSolve_5.6.15       XML_3.99-0.5        
    ## [28] pkgconfig_2.0.3      genefilter_1.70.0    purrr_0.3.4         
    ## [31] xtable_1.8-4         scales_1.1.1         annotate_1.66.0     
    ## [34] generics_0.0.2       IRanges_2.22.2       ellipsis_0.3.1      
    ## [37] withr_2.2.0          BiocGenerics_0.34.0  survival_3.1-12     
    ## [40] magrittr_1.5         crayon_1.3.4         memoise_1.1.0       
    ## [43] evaluate_0.14        doParallel_1.0.15    MASS_7.3-51.6       
    ## [46] segmented_1.2-0      class_7.3-17         graph_1.66.0        
    ## [49] tools_4.0.2          hms_0.5.3            lifecycle_0.2.0     
    ## [52] stringr_1.4.0        bcellViper_1.24.0    S4Vectors_0.26.1    
    ## [55] kernlab_0.9-29       munsell_0.5.0        AnnotationDbi_1.50.3
    ## [58] compiler_4.0.2       e1071_1.7-3          rlang_0.4.7         
    ## [61] grid_4.0.2           RCurl_1.98-1.2       iterators_1.0.12    
    ## [64] htmlwidgets_1.5.1    rappdirs_0.3.1       bitops_1.0-6        
    ## [67] rmarkdown_2.3        gtable_0.3.0         codetools_0.2-16    
    ## [70] curl_4.3             DBI_1.1.0            R6_2.4.1            
    ## [73] gridExtra_2.3        knitr_1.29           bit_4.0.4           
    ## [76] KernSmooth_2.23-17   stringi_1.4.6        parallel_4.0.2      
    ## [79] Rcpp_1.0.5           vctrs_0.3.2          dbplyr_1.4.4        
    ## [82] tidyselect_1.1.0     xfun_0.16
