library(readr)
library(viper)

working_dir <- "~/Documents/transcriptutorial/" #put whatever is your working directory here
setwd(working_dir)


source("scripts/support_functions.R")

count_df_vsn <- as.data.frame(read_csv("Documents/transcriptutorial/data/count_df_vsn.csv"))
row.names(count_df_vsn) <- count_df_vsn$gene

count_df_vsn <- count_df_vsn[,-1]
count_df_vsn <- count_df_vsn[complete.cases(count_df_vsn),]

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

viper_regulon <- df_to_viper_regulon(dorothea)

TF_activity <- as.data.frame(
  viper(eset = count_df_vsn, regulon = viper_regulon, nes = T, minsize = 5, eset.filter = F)) #most import paramter is eset.filter. With dorthea it sohuld be set to FALSE (see luz paper)

TF_activity$TF <- row.names(TF_activity)
TF_activity <- TF_activity[,c(7,1,2,3,4,5,6)]

write_csv(TF_activity,"~/Documents/transcriptutorial/data/TF_act_sample_resolution.csv")
