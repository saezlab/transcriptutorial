library(readr)
library(viper)
library(progeny)

working_dir <- "~/Projects/transcriptutorial/" #put whatever is your working directory here
setwd(working_dir)


source("scripts/support_functions.R")

count_df_vsn <- as.data.frame(read_csv("data/count_df_vsn.csv"))
row.names(count_df_vsn) <- count_df_vsn$gene

count_df_vsn <- count_df_vsn[,-1]
count_df_vsn <- count_df_vsn[complete.cases(count_df_vsn),]

PathwayActivity_counts <- progeny(as.matrix(count_df_vsn), scale=TRUE, 
                                  organism="Human", top = 100, perm = 10000, z_scores = F)
PathwayActivity_counts <- as.data.frame(t(PathwayActivity_counts))

write.csv(PathwayActivity_counts,"~/Projects/transcriptutorial/data/progeny_sample_resolution.csv")
