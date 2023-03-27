# run_Carnival.R:

Identifies potential regulated pathways between perturbed targets and active transcription factors
Uses gene expression data and the CARNIVAL pipeline
Reads a CSV file with differential activities of transcription factors and a TSV file with a subset of signed and directed interactions from the Omnipath interactions database
Infers transcription factors' activities, generates linear constraints to create integer linear programming problems, and applies an ILP solver to identify the sub-network topology with minimized fitting error and model size
Saves resulting weighted SIF and node attribute files as CSV files

# create_omnipath_carnival_pkn.R:

Extracts a subset of signed and directed interactions from the Omnipath interactions database
Ensures interactions are consistent in direction and stimulation/inhibition
Uses OmnipathR, readr, and dplyr libraries
Saves the resulting SIF file as a TSV file in the support folder
The SIF file can be used as input for the CARNIVAL pipeline to identify regulated pathways between perturbed targets and active transcription factors