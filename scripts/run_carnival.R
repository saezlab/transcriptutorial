library(CARNIVAL)
library(readr)
library(dplyr)
library(visNetwork)

source("scripts/carnival_visNet.R")

TF_differential_activities <- as.data.frame(
  read_csv("data/TF_differential_activities.csv"))
TF_differential_activities_input <- TF_differential_activities$treatment.vs.control
names(TF_differential_activities_input) <- TF_differential_activities$...1

omnipath_carnival <- read_delim("support/omnipath_carnival.tsv", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)


iniciator <- 1
names(iniciator) <- "TGFB1"

TF_differential_activities_input <- TF_differential_activities_input[abs(TF_differential_activities_input) > 2]

omnipath_carnival_pruned <- prune_network(omnipath_carnival[,c(1,3,2)], names(iniciator), names(TF_differential_activities_input))[,c(1,3,2)]

# run carnival
carnival_result = runCARNIVAL( inputObj= iniciator, 
                               measObj = TF_differential_activities_input, 
                               netObj = omnipath_carnival_pruned, 
                               solverPath = "/Users/aureliendugourd/Dropbox/cplex/cplex", #depends on your solver path
                               solver = "cplex",
                               timelimit=120,
                               mipGAP=0,
                               poolrelGAP=0,
                               threads = 6) # depends on your computer

#transoform to data.frame
SIF_network <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
SIF_network <- SIF_network[SIF_network$Weight != 0,]

ATTribute_network <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
ATTribute_network <- ATTribute_network[ATTribute_network$ZeroAct != 100,]

carnival_visNet(evis = carnival_result$weightedSIF,
                         nvis = carnival_result$nodesAttributes)



write_csv(SIF_network, file = "results/SIF_network.csv")
write_csv(ATTribute_network, file = "results/ATTribute_network.csv")
