05_network_reconstruction_with_CARNIVAL
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
al. 2019)](https://doi.org/10.1101/gr.240663.118), and the prior
knowledge network obtained from [**Omnipath**](http://omnipathdb.org/),
a literature curated mammalian signaling pathway resource [(Türei et al
2016)](https://www.nature.com/articles/nmeth.4077).

In order to help CARNIVAL to find a solution faster, we can also use
**PROGENy** [(Schubert et
al. 2018)](https://www.nature.com/articles/s41467-017-02391-6) scores to
infer the score of the representative genes for each of the calculated
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

-   Transcription factor activities
    (04_TranscriptionFactor_activity_with_Dorothea.Rmd)
-   Pathways activity scores (03_Pathway_activity_with_Progeny.Rmd).

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
(consensus_direction) and sign (consensus_stimulation and
consensus_inhibition) to extract it.

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
*assignPROGENyScores.r* to shift the formats of tf_activities and
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

-   *measObj*: The TFs’ activities (like the ones we have obtained from
    DoRothEA)
-   *inputObj*: The ‘perturbed’ nodes we want that CARNIVAL connects
    with the activity of TFs. There are 3 ways of using it:

1)  Give the name and sign of the selected nodes;
2)  Give the name only, so the algorithm will select the sign that best
    fit the models,
3)  Give *NULL* as value will create a “Perturbation” node that will try
    both signs for all ‘initial’ nodes of the given network ( *netObj*
    ).

-   *netObj*: The network that will serve as map to connect the TFs’
    activities ( *measObj* ) and the perturbed nodes ( *inputObj* )

Although it is not required, a fourth object called *weightObj* can be
also given. This object gives values ranged from -1 to 1 for a set of
nodes of the network. The aim of *weightObj* is helping the solver to
find optimal solutions faster.

In the present example, we use assign as perturbation nodes all the
“initial” nodes (option 2), and as *weightObj* the PROGENy scores
assigned to the most representative genes of the calculated pathways,

Please, check the [CARNIVAL](https://saezlab.github.io/CARNIVAL/) page
to get some more insight. For more specific technical details, as how to
install different ILP solvers that CARNIVAL supports, pleasre check
[CARNIVAL
vignette](https://bioconductor.org/packages/release/bioc/vignettes/CARNIVAL/inst/doc/CARNIVAL.html).

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
                               solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex",
                               #solverPath = "/Applications/CPLEX_Studio129/cplex/bin/x86-64_osx/cplex", 
                               solver = "cplex",
                               timelimit=60,
                               mipGAP=0,
                               poolrelGAP=0 )
```

    ## Warning in preprocessPriorKnowledgeNetwork(priorKnowledgeNetwork): self loop(s)
    ## detected and removed from prior knowledge network.

    ## Warning in checkPerturbations(perturbations, nodesPriorKnowledgeNetwork): These
    ## perturbation nodes are provided as NA, their activity will be inferred: MDFI
    ## | HOMER1 | ITPR3 | ASPH | RNF24 | CABP1 | SESTD1 | MX1 | NPHS2 | SNF8 | CNR1 |
    ## PIRT | MAP7 | PACSIN3 | BSPRY | CALB1 | KL | RAB11A | S100A10 | PDZD3 | NIPSNAP1
    ## | NYX | EFHC1 | PTPN13 | PRKCSH | ESPN | ZIC1 | TUB | CDC73 | SYMPK | MYCN |
    ## KLF2 | NPAS2 | ATAD2 | BACE2 | PREX2 | PTPRD | RPL11 | SPNS1 | ATG4B | ADGRV1 |
    ## RALB | PRKACG | PKNOX1 | ATG101 | MIB1 | HRK | UBE2T | KDM2B | IL22 | CXCL16 |
    ## CCL3L1 | ASF1B | NMI | CARD18 | ISG15 | CARD16 | ANO9 | KLF6 | FLT3LG | KITLG |
    ## SERPINB9 | APLF | E2F7 | TLX1 | PLK3 | GTF3C4 | KDM3B | KDM3A | KDM4B | KDM4A |
    ## APIP | CAPNS1 | XRCC5 | KMT2B | PRMT6 | SETD1A | CXCL9 | CCL8 | CCL4 | SETDB2 |
    ## KMT2E | RIOX1 | SETD7 | PRDM2 | DOT1L | ATG4A | VMP1 | SUPT20H | RAB33B | CISD2
    ## | UBE2M | OS9 | SSBP1 | RORB | NR1D2 | CAND1 | GLRX | SIRT6 | NKRF | FRAT2 |
    ## MIR22HG | CDK11A | SMARCB1 | NSMCE2 | TRIM8 | RECQL5 | MARCHF5 | DTX3L | NFKBIZ
    ## | ASH1L | TNFSF14 | CD40LG | TNFSF13B | LTA | LTB | NRG1 | USP22 | NKIRAS1 |
    ## NKIRAS2 | ACAA2 | ING1 | USP17L2 | KMT5C | SIRT7 | EIF2AK1 | CDK5RAP3 | TIMD4 |
    ## CSNK1E | EHMT2 | KAT6A | JARID2 | KDM6A | TRDMT1 | MEN1 | PHF8 | SETD1B | GZMA
    ## | DZIP3 | PAF1 | PHF1 | HELLS | SLC2A4RG | INSIG1 | DUSP22 | FAAP24 | ADRM1 |
    ## MAT1A | NEDD8 | ACTL6A | CEP85 | SENP1 | HDAC10 | SUMO3 | ZNF350 | EP400 | CHD8
    ## | SOX17 | MED14 | CDC14B | BIRC7 | TLR10 | TLR8 | MPG | HIF1AN | SUMO2 | NF1
    ## | RNF135 | MCM10 | FBXW8 | KDM7A | KDM1B | USP9X | GMNN | ING2 | GFI1B | ERCC8
    ## | SHPRH | WSB1 | SUV39H2 | NANOS1 | RPS19BP1 | CTNNBIP1 | AKAP8 | PSME3 | LMO4
    ## | SENP3 | CTSH | CTSK | CTSS | ZFHX3 | CTSB | GNRH1 | CHFR | DENND4A | BRCC3 |
    ## RING1 | HSBP1 | RAD54L | PAX5 | OTUB1 | EEF1E1 | TIPRL | DYNLL2 | BAG4 | TNFSF9
    ## | BAG1 | GADD45B | HUS1 | SEM1 | NAB2 | NAB1 | CD70 | BLVRA | GNL3 | BABAM1 |
    ## UBE2D3 | IL24 | PRKRA | CKS1B | XAF1 | PPP1R13L | PPP1R13B | WRAP53 | PDK2 |
    ## CASP12 | HTATIP2 | FMR1 | CIB1 | RILPL2 | RAC3 | CHP1 | PPM1E | PPM1F | UBR2
    ## | GCKR | CXCL11 | RSF1 | CD80 | CD86 | EMSY | UBR1 | FKBP1A | SMARCAD1 | TAOK2
    ## | NDP | GH1 | TSLP | IL15 | NRG3 | FRZB | IHH | DUSP8 | NRG4 | SHOC2 | IL23A |
    ## CNTF | FBXW11 | CISH | DISP1 | CRIPAK | IL31 | OSM | IL21 | AMH | BMP6 | GDF6 |
    ## PPP4R3A | CDO1 | IL20 | IL26 | IL19 | IL27 | GPC3 | MFAP5 | CTF1 | MTMR4 | HYAL2
    ## | DACT2 | ANKRD6 | PDGFC | FGF4 | FGF7 | FGF9 | VEGFB | EFNA1 | EFNA2 | EFNA5 |
    ## EFNB3 | EFNA3 | EFNA4 | EFNB2 | LTK | ANGPT4 | AGRN | DOK7 | PTK6 | NR0B2 | IL5
    ## | DIXDC1 | TDGF1 | SYNJ2BP | UBE2O | RLIM | RAB23 | SCO2 | THPO | IL7 | DUSP2
    ## | NRG2 | LAMTOR3 | DHH | DUSP9 | WIF1 | BMP8B | DTX4 | IL9 | MBIP | SNIP1 |
    ## WNT10B | DUSP10 | MAP3K13 | EXTL1 | DUSP19 | ELP1 | NUP214 | FOXH1 | CBY1 | IL11
    ## | CNKSR2 | NKD1 | WDR83 | TPBG | TGIF1 | POFUT1 | NEURL1 | MTSS1 | GRK6 | DISP2
    ## | DISP3 | TRIM33 | CNKSR1 | SPOP | MIB2 | GDF1 | BMP15 | NODAL | NUMBL | SOCS2
    ## | CTDSP2 | BCL9 | CTDNEP1 | IFNK | IL36G | IFNL1 | IFNL2 | IFNL3 | MDK | HHAT |
    ## CNTN1 | CCN3 | MFAP2 | PTPRT | INVS | LRRFIP2 | CAPRIN2 | PDGFD | FGF6 | FGF17
    ## | FGF18 | FGF22 | PTN | PTPRB | ANGPTL1 | AATK | RYBP | RFX1 | NR2E3 | NR2C1
    ## | PRKX | GNG4 | PTPA | ALOX5AP | CCL19 | CCL21 | CCL25 | INHBB | INHA | ADM |
    ## ADM2 | CALCB | APELA | NMB | GRP | TNFSF13 | BMP10 | BMP8A | GDF2 | GAST | CCL13
    ## | CCL14 | CCL15 | CCL23 | CCL7 | CCL27 | CCL28 | CCL16 | CCL11 | CCL24 | CCL26
    ## | CCL17 | CCL22 | CCL20 | CCL1 | TNFSF8 | RARRES2 | CLCF1 | CRH | UCN | UCN2 |
    ## UCN3 | YARS1 | CXCL6 | CXCL2 | CXCL3 | CXCL13 | TNFSF15 | POMC | PDYN | PENK |
    ## EDA | EPGN | FGF19 | VEGFC | VEGFD | MT-RNR2 | SAA1 | HEBP1 | FSHB | SPX | GHRL
    ## | GHRH | GIP | GCG | TNFSF18 | GNRH2 | GDF10 | GDF3 | GDF7 | GDF9 | GH2 | NPPA
    ## | NPPB | NPPC | GUCA2A | GUCA2B | GUCA1A | GUCA1B | PSAP | RSPO1 | RSPO2 | RSPO3
    ## | RSPO4 | IFNA10 | IFNA14 | IFNA16 | IFNA17 | IFNA2 | IFNA21 | IFNA4 | IFNA5 |
    ## IFNA6 | IFNA7 | IFNA8 | IFNW1 | IL17C | IL17A | IL17F | IL1RN | IL17B | IL25 |
    ## EBI3 | IL36A | IL36B | IL36RN | KISS1 | CGB8 | LHB | PMCH | AGRP | MLN | OSTN
    ## | NTF4 | TAC1 | TAC3 | NMS | NMU | PNOC | NPB | NPW | NPFF | NPVF | NPS | NTS |
    ## OXT | AVP | TNFSF4 | HCRT | ADCYAP1 | VIP | PROK1 | CSH1 | CSH2 | CSHL1 | PRLH
    ## | PTH2 | QRFP | RLN2 | RLN1 | RLN3 | INSL3 | INSL5 | SCT | CORT | TSHB | TNFSF12
    ## | UTS2 | UTS2B | XCL1 | XCL2 | PYY | PPY | RAPGEF4 | NLGN1 | BCAP31 | ADORA1 |
    ## VDAC2 | GDI1 | ADORA2A | PSD3 | PRKG2 | RASA4 | RAPGEF3 | ARFGEF2 | SIPA1 | RGS4
    ## | SNPH | HOMER2 | HPCA | MYRIP | AP2S1 | RELN | STXBP6 | CNR2 | SRI | DOC2A |
    ## DNAJC5 | AP1B1 | CERK | ENC1 | PTMA | MT3 | GADD45GIP1 | PGAM5 | ERC1 | CSNK2A2
    ## | DHX58 | ZBP1 | CARD17 | PYDC1 | RNF125 | HLA-DMA | HLA-DMB | CASP4 | CARD6
    ## | POP1 | CNOT8 | PYDC2 | POPDC2 | SLC9A3R2 | MAGI3 | CTNND2 | DLGAP2 | TAF12
    ## | PCSK7 | POU2AF1 | SRSF7 | PDK4 | CDKN3 | PTPN21 | DUSP7 | ARHGAP44 | MYLK3
    ## | ARHGAP11B | ADCK5 | LARP7 | GXYLT1 | FILIP1L | GREB1 | ARHGAP29 | ARHGAP15 |
    ## PLEKHG4 | ZC3H12A | BCORL1 | PITRM1 | FGD3 | AMER1 | TOR1AIP1 | IQSEC2 | FFAR4 |
    ## PPM1L | THEM4 | FKBP15 | ARTN | ARHGAP21 | UBAP2 | WLS | MAP3K21 | TTBK1 | FHL5
    ## | ARHGAP40 | TSPO2 | SYDE2 | PDE4DIP | METTL21C | RPGRIP1L | PDPK2P | ITPRIPL1
    ## | NCAPH2 | PLEKHA7 | TTBK2 | PTCRA | STYK1 | ZNF774 | MAST2 | PNCK | ARHGAP11A
    ## | PKN3 | ULK3 | FBXO38 | LRFN4 | SLC30A9 | TRAF7 | APOA5 | TSSK4 | FIP1L1 | NPNT
    ## | DLK2 | BCOR | POGZ | DNMBP | GLDN | IL34 | NBEAL2 | FGD5 | BICDL1 | SRCAP |
    ## RASSF6 | ARHGAP27 | PHLPP2 | SSH2 | ASXL2 | ADAMTS13 | SND1 | ZC3HAV1 | PEX26 |
    ## HAUS6 | MDGA2 | GOLGA7 | ARHGAP22 | IRF2BP2 | SRGAP1 | ARHGAP30 | FGD2 | NRARP |
    ## RNF180 | LDB1 | CDH24 | NLRX1 | CADPS2 | PIM3 | IFTAP | HOOK3 | JAZF1 | LONP2 |
    ## FBXO11 | VRK2 | ERO1B | IRF2BP1 | PLD3 | VRK3 | TIAM2 | ARHGEF19 | TRIM59 | PHF6
    ## | ARHGAP12 | TP53INP2 | FAM20C | RPAP2 | COLGALT2 | SEPTIN12 | DOCK3 | CDK20 |
    ## ABCA7 | NLGN4X | TAGAP | DOCK4 | ARHGEF28 | UNC80 | STK40 | PPM1K | RASGEF1C |
    ## RNF10 | EID2 | BANP | XXYLT1 | COLGALT1 | POGLUT1 | NECAP1 | TSC22D3 | HJURP |
    ## CCNY | BVES | CSNK2A3 | FBXO22 | SYNE1 | NBEA | DNER | NLGN4Y | TRIM58 | PJA1 |
    ## SVIP | TERB2 | WDR48 | GABPB2 | UHMK1 | LNX1 | PLPPR1 | MRAP | MIPOL1 | NEK9 |
    ## DMXL2 | SSH3 | RAPGEF6 | SMCR8 | WHAMM | FNIP1 | UBASH3B | STYX | BRK1 | NUDCD2
    ## | ZFPM2 | APH1B | DYDC1 | ATXN2L | ARAP3 | PHIP | LIPH | STON2 | AUTS2 | ARAP2
    ## | TFAP2B | DDX1 | STARD8 | TM9SF4 | BAP1 | RAPGEF5 | PHF3 | CNOT9 | PXDN | LPAR1
    ## | PRCC | TFG | THRSP | DPF2 | SEMA4D | FGF13 | FGF11 | FGF14 | KCNB2 | WNT2B
    ## | WNT8B | FBXW5 | RNF34 | MUL1 | CDCA5 | ERGIC1 | FAM162A | RAD51AP1 | PHF21A
    ## | DOCK10 | EAF2 | GCC1 | NTNG2 | RAB39B | HDAC11 | ARHGEF26 | HOOK2 | TRIM11 |
    ## PIK3IP1 | MRAP2 | P2RY11 | PDXP | CENPN | ERO1A | INTS4 | PAWR | CFAP53 | EAF1
    ## | ZNF462 | CHAMP1 | ZFP91 | RITA1 | ZNF521 | MEGF10 | RPGRIP1 | NSD1 | TRIM47 |
    ## FGD4 | NSMCE3 | SPATA13 | PTCHD1 | CAMK1G | OXGR1 | ARHGEF17 | PLEKHG4B | NEK1
    ## | CSMD1 | PHF12 | NACC1 | PASK | PIGS | TESK2 | CLCC1 | FIZ1 | RBM15 | GPHA2
    ## | S1PR3 | GPER1 | EYA3 | NEU1 | MMP19 | SCAF11 | NDN | MAP3K6 | NRTN | TBL1Y
    ## | PKP2 | SEMA3C | P2RY13 | FERMT1 | ARHGAP9 | CENPK | NDFIP1 | CRB3 | DUSP26 |
    ## COLEC11 | GPR174 | KLF16 | WNK3 | KDM5D | ACE2 | RBCK1 | NSD3 | ASPSCR1 | UPF3B
    ## | SPRY4 | SETD5 | ZDHHC5 | PHACTR1 | ARHGAP39 | WNT10A | CDH19 | SIL1 | IRF2BPL
    ## | LPAR5 | VPS11 | S1PR5 | PORCN | P2RY12 | CDH23 | ADNP | DPH5 | RANBP17 | PARL
    ## | GHITM | HIPK3 | SAV1 | USP8 | SLA2 | PLEKHF2 | GRPEL1 | CLK4 | PRDM16 | NTN4
    ## | PATZ1 | CDH20 | MOV10 | SDF2L1 | IL1RAPL2 | SSU72 | MMP25 | LTB4R2 | NOP10
    ## | RTN4 | GPR84 | ARHGEF3 | DISC1 | CYSLTR2 | LRRC4B | RNF146 | MDN1 | POT1 |
    ## SDHAF2 | NHP2 | SIRT5 | FBXL12 | GAR1 | DACT1 | FAM13B | CDK12 | PLK2 | NLGN3
    ## | AHSP | MAT2B | C1RL | LRP1B | BFAR | SH3BP4 | GMIP | PIM2 | ARHGAP23 | FNIP2
    ## | PAK5 | STK26 | CTTNBP2NL | ARHGAP20 | PDP2 | ARHGAP28 | RERE | GALP | ASH2L
    ## | DNAJB9 | FAM8A1 | WNT16 | LPAR3 | CLN8 | TJP2 | STK17A | FBXL17 | RIMBP3 |
    ## UBQLN2 | CTNNA3 | LCMT1 | MUTYH | BAZ2B | ZDHHC2 | CDH22 | HOOK1 | STOML2 |
    ## PTBP2 | RALY | MUC12 | BAG5 | MYT1L | PLEKHG1 | GPC6 | HECTD1 | CADPS | AKAP8L
    ## | PCDHA7 | PCDHA6 | PCDHA4 | PCDHA12 | CDC14A | GPR132 | TTF2 | GPR34 | ANKRD26
    ## | SATB2 | CNTN6 | AURKC | CYSLTR1 | BCS1L | DLGAP4 | MAST1 | NTNG1 | KDM2A |
    ## PPP2R2C | THRAP3 | SBDS | TMED5 | BZW2 | SH3BP1 | STARD13 | GRIP1 | WNK2 | RWDD3
    ## | HECTD4 | LAS1L | PCDHA9 | PCDHA8 | PCDHA5 |

    ## Warning in checkMeasurements(measurements, nodesPriorKnowledgeNetwork): These
    ## measurement nodes are not in prior knowledge network and will be ignored: DUX4 |
    ## SIX2 | SNAI2 | SNAPC4 | THAP11 | ZNF263

    ## Warning in checkWeights(weights, nodesPriorKnowledgeNetwork): These nodes are
    ## not in prior knowledge network and will be ignored: P53

CARNIVAL gives a list of 4 elements:

-   weightedSIF: summary of all interactions found in all models
-   nodesAttributes: summary of all nodes and how many times are they
    found in the different models
-   sifAll: networks of all the models
-   attributesAll: node attributes of all models

We can now visualise the network…

``` r
#transoform to data.frame
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)
carnival_result$weightedSIF = carnival_result$weightedSIF[carnival_result$weightedSIF$Weight !=0,]

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

> Liu A., Trairatphisan P., Gjerga E. et al. “From expression footprints
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

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggrepel_0.9.1    pheatmap_1.0.12  ggplot2_3.3.6    visNetwork_2.1.0
    ##  [5] dplyr_1.0.9      tidyr_1.2.0      tibble_3.1.8     readr_2.1.2     
    ##  [9] OmnipathR_3.4.0  CARNIVAL_2.7.2   dorothea_1.8.0   progeny_1.18.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.9         prettyunits_1.1.1  assertthat_0.2.1   digest_0.6.29     
    ##  [5] utf8_1.2.2         R6_2.5.1           cellranger_1.1.0   plyr_1.8.7        
    ##  [9] backports_1.4.1    evaluate_0.16      httr_1.4.4         pillar_1.8.1      
    ## [13] rlang_1.0.4        progress_1.2.2     curl_4.3.2         readxl_1.4.0      
    ## [17] rstudioapi_0.13    checkmate_2.1.0    rmarkdown_2.15     stringr_1.4.1     
    ## [21] bcellViper_1.32.0  htmlwidgets_1.5.4  bit_4.0.4          igraph_1.3.4      
    ## [25] munsell_0.5.0      compiler_4.2.0     xfun_0.32          pkgconfig_2.0.3   
    ## [29] htmltools_0.5.3    tidyselect_1.1.2   gridExtra_2.3      lpSolve_5.6.15    
    ## [33] fansi_1.0.3        crayon_1.5.1       tzdb_0.3.0         withr_2.5.0       
    ## [37] later_1.3.0        rappdirs_0.3.3     grid_4.2.0         jsonlite_1.8.0    
    ## [41] gtable_0.3.0       lifecycle_1.0.1    DBI_1.1.3          magrittr_2.0.3    
    ## [45] scales_1.2.0       vroom_1.5.7        cli_3.3.0          stringi_1.7.8     
    ## [49] reshape2_1.4.4     xml2_1.3.3         logger_0.2.2       ellipsis_0.3.2    
    ## [53] generics_0.1.3     vctrs_0.4.1        RColorBrewer_1.1-3 rjson_0.2.21      
    ## [57] tools_4.2.0        bit64_4.0.5        glue_1.6.2         purrr_0.3.4       
    ## [61] hms_1.1.2          parallel_4.2.0     fastmap_1.1.0      yaml_2.3.5        
    ## [65] colorspace_2.0-3   knitr_1.39
