library(OmnipathR)
library(readr)
library(dplyr)

omniR <- import_Omnipath_Interactions()

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
write_tsv(sif, "support/omnipath_carnival.tsv")
