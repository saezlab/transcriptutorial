library(readr)
library(CARNIVAL)

source("assignPROGENyScores.r")
source("generateTFList.r")

working_dir <- "~/Projects/transcriptutorial/scripts" #put whatever is your working directory here
setwd(working_dir)

#files 

tf_activities <- read_csv("../data/TF_act_sample_resolution.csv")
PathwayActivity <- read_csv("../data/progeny_sample_resolution.csv")
sif = read_tsv("../results/omnipath_carnival.tsv")



# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$TF
tf_activities_carnival$TF <- NULL
tfList = generateTFList(tf_activities_carnival, top=50, access_idx = 1:ncol(tf_activities_carnival))

# progeny for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival[,1]
PathwayActivity_carnival[,1] <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                  progenyMembers = progenyMembers, 
                                  id = "gene", 
                                  access_idx = 1:ncol(PathwayActivity_carnival))
# get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

# run carnival for all samples
carnival_result = list()
for ( n in names(tfList) ){
  carnival_result[[n]] = runCARNIVAL( inputObj= iniciators,
                                 measObj = tfList[[n]], 
                                 netObj = sif, 
                                 weightObj = progenylist[[n]], 
                                 solverPath = "/Applications/CPLEX_Studio129/cplex/bin/x86-64_osx/cplex", 
                                 solver = "cplex",
                                 timelimit = 7200,
                                 mipGAP = 0,
                                 poolrelGAP = 0 )
  
  #transoform to data.frame
  carnival_result[[n]]$weightedSIF <- data.frame(carnival_result[[n]]$weightedSIF, stringsAsFactors = F)
  carnival_result[[n]]$weightedSIF$Sign <- as.numeric(carnival_result[[n]]$weightedSIF$Sign)
  carnival_result[[n]]$weightedSIF$Weight <- as.numeric(carnival_result[[n]]$weightedSIF$Weight)
  
  carnival_result[[n]]$nodesAttributes <- data.frame(carnival_result[[n]]$nodesAttributes, stringsAsFactors = F)
  carnival_result[[n]]$nodesAttributes$ZeroAct <- as.numeric(carnival_result[[n]]$nodesAttributes$ZeroAct)
  carnival_result[[n]]$nodesAttributes$UpAct <- as.numeric(carnival_result[[n]]$nodesAttributes$UpAct)
  carnival_result[[n]]$nodesAttributes$DownAct <- as.numeric(carnival_result[[n]]$nodesAttributes$DownAct)
  carnival_result[[n]]$nodesAttributes$AvgAct <- as.numeric(carnival_result[[n]]$nodesAttributes$AvgAct)
 
  saveRDS(carnival_result,"../results/carnival_sample_resolution.rds")   
}
