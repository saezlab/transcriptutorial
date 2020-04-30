#Copyright (C) 2019  Aurelien Dugourd
#Contact : aurelien.dugourd@bioquant.uni-heidelberg.de

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#############################
#######   WHADIZDIS    ######
#############################

#This script is designed to guide the analysis one RNAseq data with the help of footprint based analysis tools
#such as DOROTHEA/Viper, Progeny and CARNIVAL
#Here we use an example dataset to guide the analysis

#Main libraries
library(readr)
library(vsn)
library(limma)
library(viper)

#Support functions also requires
library(ggplot2)
library(reshape)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)
library(hexbin)


working_dir <- "~/Documents/transcriptutorial/" #put whatever is your working directory here
setwd(working_dir)


source("scripts/support_functions.R")

### Preparing dorothea
url <- paste0(
  'http://omnipathdb.org/interactions?',
  'datasets=tfregulons&tfregulons_levels=A,B&genesymbols=1&fields=sources,tfregulons_level'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

##Dorothea/viper
dorothea <- download_omnipath()
dorothea <- dorothea[,c(4,3,6,7)]
dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
dorothea$sign <- ifelse(dorothea$sign == 0, 1, dorothea$sign)
dorothea <- dorothea[,c(1,2,5)]

### Import the raw count dataframe
#downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119931
#download the file : GSE119931_PANC1.FOXA2KO.genes.counts.txt.gz and decompress it in the data folder

GSE119931_PANC1_FOXA2KO_genes_counts <- as.data.frame(
  read_delim("data/GSE119931_PANC1.FOXA2KO.genes.counts.txt", 
                                                   "\t", escape_double = FALSE, trim_ws = TRUE)) 
count_df <- GSE119931_PANC1_FOXA2KO_genes_counts[,c(7:12)]
row.names(count_df) <- GSE119931_PANC1_FOXA2KO_genes_counts$Geneid

#### Pre-processing and normalisation

#First create a dataframe to summarise experimental design called targets

targets <- as.data.frame(matrix(NA,length(names(count_df)),2))
names(targets) <- c("sample","condition")
targets$sample <- names(count_df)
targets$condition <- gsub(".Rep[0-9]$","",targets$sample)

#Make some plots to check what the data looks like after only a log2 transformation

#First we remove rows that contain only 0
count_df <- count_df[rowSums(count_df) > 0,]
#remaining 0 have to be made as NA so that log2 transformation is possible
count_df[count_df == 0] <- NA

#make the plots
magicPlotMaker(df = log2(count_df), outpath = paste(working_dir,"visualisation/raw_log2/", sep = ""), targets = targets)
dev.off()
#from the density plot we define the minimum cutoff of expression. 
#Distributions are weird under log2(count) of around 4

count_df[log2(count_df) < 4 ] <- NA

#remove rows that don't have enough well measured genes in enough samples

count_df <- count_df[rowSums(is.na(count_df[,c(1:3)])) < 2,]
count_df <- count_df[rowSums(is.na(count_df[,c(4:6)])) < 2,]

### VSN normalisation
#now we can normalise the cleaned dataframe using vsn
fit <- vsnMatrix(as.matrix(count_df)) #train vsn parameters

#make sure the mean/sd trend is not going crazy
meanSdPlot(fit)

#if good, normalise data with the trained parameters of vsn
count_df_vsn <- as.data.frame(vsn::predict(fit,as.matrix(count_df)))

#now let's visualise the normalised data
magicPlotMaker(df = count_df_vsn, outpath = paste(working_dir,"visualisation/vsn/", sep = ""), targets = targets) 
dev.off()
#from PCA, we see that conditions are well seprated by 2nd component. So it's ok, we will have some signal.

### Identifier kung-fu (optional)
#since here with have ensembl id but most our ressources are based on either uniprot or gene symbole
#we need to do some identifer kung-fu

#I got this identifer matching dataframe from uniprot
gene_id_mapping_from_uniprot <- as.data.frame(
  read_delim("support/gene_id_mapping_from_uniprot.tab", 
                                           "\t", escape_double = FALSE, trim_ws = TRUE))
gene_id_mapping_from_uniprot <- gene_id_mapping_from_uniprot[!is.na(gene_id_mapping_from_uniprot$`Gene names`),]

#let's make a pseudo dictionary to make the mapping efficient
ensembl_to_symbol <- gsub(" .*","",gene_id_mapping_from_uniprot$`Gene names`)
names(ensembl_to_symbol) <- gene_id_mapping_from_uniprot[,1]

#remove all genes that have no gene symbol from our count dataframe
row.names(count_df_vsn) <- gsub("[.][0-9]*","",row.names(count_df_vsn))
count_df_vsn <- count_df_vsn[row.names(count_df_vsn) %in% names(ensembl_to_symbol),]

#now let's convert ids with the pseudo dictionary
for(i in 1:length(count_df_vsn[,1]))
{
  row.names(count_df_vsn)[i] <- ensembl_to_symbol[row.names(count_df_vsn)[i]]
}

### LIMMA differential analysis
#now let's run a simple differential analysis using a simple wrapper for such situation

#first check the conditions order
unique(targets$condition)

#we want to compare the KO with the WT so we build a comparison list
comparisons <- list(c(2,-1)) #each vector of the list represent the contrasts, here we substract the first condition (-1) to the second one (2)

#now that the comparisons are defined, we can run limma
limmaRes <- runLimma(measurements = count_df_vsn, 
                     targets = targets, 
                     comparisons = comparisons)

#once limma has run, we extract the statistic dataframe summarise the differential analysis
ttop_KOvsWT <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(count_df_vsn[,1]), adjust.method = "fdr"))

##make a qqplot
null_model <- pnorm(rnorm(length(ttop_KOvsWT[,1])))
plot(sort(null_model), sort(ttop_KOvsWT$P.Value) ,xlim = c(1,0), ylim = c(1,0)) #not bad, not great, let's proceed
abline(coef = c(0,1))

### DOROTHEA

#first we need to format the dorothea dataframe we prepared into a proper viper input
viper_regulon <- df_to_viper_regulon(dorothea)

#then we format the input gene level statistic for viper
#here we only have one condition so it will be a named vector
#if you have more than one condition, you need to make a matrix of t_values with gene names as rownames and column as conditions

eset <- ttop_KOvsWT$t
names(eset) <- ttop_KOvsWT$ID

#run viper and you should get a dataframe of TF activities
TF_activity <- as.data.frame(
  viper(eset = eset, regulon = viper_regulon, nes = T, minsize = 5, eset.filter = F)) #most import paramter is eset.filter. With dorthea it sohuld be set to FALSE (see luz paper)

#to interpret the results, we can look at the expression of targets of the most deregulated TFs
volcano_nice(ttop_KOvsWT[ttop_KOvsWT$ID %in% names(viper_regulon$E2F4$tfmode),], #example E2F4
             FCIndex = 2, pValIndex = 5, IDIndex = 1,nlabels = 20, label = T)

### Progeny

#first we import a model matrix
prog_matrix <- as.data.frame(
  read_csv("support/model_NatComm+14_human.csv"))

#then we format the gene level statistic as a dataframe of t_values, here we only have one condition
t_table <- ttop_KOvsWT[,c(1,4)]

#then we run progeny enrichment function
progeny_activities <- runProgenyFast(df = t_table, weight_matrix = prog_matrix, k = 10000, z_scores = T)
names(progeny_activities) <- "NES"

#We can make scatter plot of gene_level statistic against weights for each pathway to hepl interpret the results
scat_plots <- progenyScatter(df = t_table, weight_matrix = prog_matrix, statName = "t_values")

#visualise MAPK responsive genes
plot(scat_plots[[1]]$MAPK)

### NEXT WE DO CARNIVAL

##Import and generate a causal network from OMNIPATH

url <- paste0(
  'http://omnipathdb.org/interactions?',
  'fields=sources,references&genesymbols=1'
)

download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

omnipath <- download_omnipath()
omnipath <- omnipath[omnipath$is_stimulation != 0 | omnipath$is_inhibition != 0,]
omnipath_sif <- omnipath[omnipath$is_stimulation ==1,c(3,6,4)] #First we get the activation
omnipath_sif_2 <- omnipath[omnipath$is_inhibition ==1,c(3,7,4)] #Then we get the inhibtion

names(omnipath_sif) <- c("source","sign","target")
names(omnipath_sif_2) <- c("source","sign","target")

omnipath_sif_2$sign <- omnipath_sif_2$sign * -1

#Then we bind together activations and inhibtion to get the complete network
omnipath_sif <- as.data.frame(rbind(omnipath_sif,omnipath_sif_2)) 

omnipath_sif$source <- gsub("[-+{},;() ]","___",omnipath_sif$source)
omnipath_sif$target <- gsub("[-+{},;() ]","___",omnipath_sif$target)

#We define the end points that we are trying to reach in the network, from the initial perturbation
TF_carni_inputs <- as.data.frame(t(TF_activity))
TF_carni_inputs <- TF_carni_inputs[,c(1:50)]
names(TF_carni_inputs) <- gsub("[-+{},;() ]","___",names(TF_carni_inputs))

#FOXA2 was knocked out so we define it as the initially perturbed node
Perturbation_carni_input <- as.data.frame(matrix(-1,1,1))
names(Perturbation_carni_input) <- "FOXA2" #FOXA2 was knocked out so we define it as the initially perturbed node

result1 = runCARNIVAL(solverPath=solverPath, 
                      netObj = omnipath_sif, 
                      measObj = TF_carni_inputs, 
                      inputObj = Perturbation_carni_input, 
                      weightObj = NULL, 
                      timelimit = 3600, 
                      poolCap = 100,
                      poolReplace = 2, 
                      limitPop = 500,
                      DOTfig = TRUE, 
                      solver = "cplex")
