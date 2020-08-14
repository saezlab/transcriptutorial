#'\code{count_edges_nodes_degree}
#'
#'Calculate the number of edges, nodes and the density of a network based on a sif file.
#'@param sif Data frame of values source, interaction, target
#'@return numeric vector with values of number of edges, nodes and density
#'
#'Rosa Hernansaiz Ballesteros, 2020
count_edges_nodes_degree <- function(sif){
  edg = nrow(sif)
  nod = length(union(sif$Node1, sif$Node2))
  
  # remove the double interactions
  edg_dir = sif %>%
    dplyr::select(Node1, Node2) %>%
    unique() %>% nrow()
  
  if (edg != edg_dir){
    message( paste0("There are ", edg - edg_dir, " extra interactions with different sign") )
  }
  
  den = edg_dir / (nod * (nod - 1))
  
  return(c("edges" = edg, "nodes" = nod, "density" = den))
}

#'\code{degree_count}
#'
#'Calculate the degree distribution  based on a sif file.
#'@param sif Data frame of values source, interaction, target
#'@return data.frame with counts of in/out-degree and per sign
#'
#'Rosa Hernansaiz Ballesteros, 2020
degree_count <- function(sif){
  nod = union(sif$Node1, sif$Node2)
  
  dd = data.frame(do.call(rbind,
                          lapply(nod, function(n){
                            # in-degree
                            in_degree = sif %>% dplyr::filter(Node2 == n) %>% nrow()
                            in_positive = sif %>% dplyr::filter(Node2 == n & Sign == 1)  %>% nrow()
                            in_negative = sif %>% dplyr::filter(Node2 == n & Sign == -1)  %>% nrow()
                            
                            # out-degree
                            out_degree = sif %>% dplyr::filter(Node1 == n) %>% nrow()
                            out_positive = sif %>% dplyr::filter(Node1 == n & Sign == 1)  %>% nrow()
                            out_negative = sif %>% dplyr::filter(Node1 == n & Sign == -1)  %>% nrow()
                            
                            r = c("node" = n, "total_count" = in_degree + out_degree,
                                  "in_count" = in_degree, "in_positive" = in_positive, "in_negative" = in_negative,
                                  "out_count" = out_degree, "out_positive" = out_positive, "out_negative" = out_negative
                            )
                            return(r)
                            
                          })
                          )
                  )
  return(dd)
}

#'\code{getTopology}
#'
#'Creation of an matrix which values are the edge's weight for all samples.
#'Sing of interaction is kept in the edge name. 
#'Positive interactions are edge1+edge2, while negative ones are edge1-edge2.
#'@param scafoldNET Data frame of values source, interaction, target
#'@param networks list of data frames with subsets of interactions from scafoldNET.
#'The list must contain the names of each element
#'Same values as scafoldNET, but it calso contains a column with
#'weights
#'@param inversCARNIVAL Boolean, default False. If true, 'Perturbation' nodes are added to
#'the first leaves.
#'
#'@return Data frame with edges in rows, samples in columns and 
#'which values are the weights
#'@export
#'
#'Rosa Hernansaiz Ballesteros, 2020

getTopology <- function(networks=NULL, scafoldNET=NULL, inversCARNIVAL=F){
  
  # set names in scafoldNET and change - for _ 
  colnames(scafoldNET) = c('source', 'interaction', 'target')
  scafoldNET$source = gsub("-", "_", scafoldNET$source, ignore.case = FALSE, perl = FALSE,
                           fixed = FALSE, useBytes = FALSE)
  scafoldNET$target = gsub("-", "_", scafoldNET$target, ignore.case = FALSE, perl = FALSE,
                           fixed = FALSE, useBytes = FALSE)
  
  # create edges names: edge1SIGNedge2
  edges = unique(paste(scafoldNET$source, 
                       gsub("1","+",gsub("-1","-",scafoldNET$interaction)), 
                       scafoldNET$target, 
                       sep = ''))
  
  # add perturbation nodes
  if(inversCARNIVAL){
    initiator = base::setdiff(scafoldNET$source,scafoldNET$target)
    edges = c(edges, paste("Perturbation", initiator, sep="+"), paste("Perturbation", initiator, sep="-"))
  }
  
  # create matrix
  topologies = data.frame(matrix(data = NA, 
                                 nrow = length(edges), 
                                 ncol = length(networks),
                                 dimnames = list(edges, names(networks))),
                          stringsAsFactors = F)
  
  # set progress bar
  pb = txtProgressBar(min = 0, max = length(names(networks)), initial = 0, style = 3) 
  
  # fill in topologies
  for (i in 1:length(names(networks))){
    #update progress bar
    setTxtProgressBar(pb,i)
    #get samples name
    sample = names(networks)[i]
    
    # set names of the df with weights
    colnames(networks[[sample]]) = c('source', 'interaction', 'target', 'weight')
    
    # get edges names (edge1SIGNedge2) and their weights
    edgeWeight = t(apply(networks[[sample]], 1, function(row){
      edge = paste(row['source'], 
                   gsub("1","+",gsub("-1","-", row['interaction'])), 
                   row['target'], 
                   sep="")
      edge = gsub(" ","", edge)
      weight = as.numeric(row['weight'])
      return(c(edge, weight))
    }))
    # set colnames and type 
    edgeWeight = as.data.frame(edgeWeight, stringsAsFactors=F)
    colnames(edgeWeight) = c('edge', 'weight')
    
    #Fill the matrix with the weight
    for (edge in edges){
      eW = edgeWeight[edgeWeight$edge==edge,]
      # if the edge exists
      if (nrow(eW)!=0){
        topologies[edge,sample] = as.numeric(eW$weight)
      }else{ next }
    }
    
  }
  #remove rows withouth information
  topologies = topologies[apply(topologies, 1, function(row){any(!is.na(row))}),]
  
  return(topologies)
}

#'\code{getCoreInteractions}
#'
#'Comparison of several topologies to extract the shared topology
#'@param topology data frame of interactions "protein1 +/- protein2" per samples (output of getTopology)
#'@param psmpl Numeric [0,100]. 95 by default. Percentage of samples that 
#'share an interaction to consider it as core. 
#'
#'@return data frames which values are source, interaction and target shared among all elements
#'Rosa Hernansaiz Ballesteros, 2020
#'
getCoreInteractions <- function(topology = NULL, psmpl = 95){
  
  #create SIF of shared interactions
  topology[!is.na(topology)] <- 1
  
  interactions_shared = apply(topology, 1, sum, na.rm=T)
  psmpl = psmpl/100 * length(colnames(topology))
  aux = names(interactions_shared)[which(interactions_shared >= psmpl)]
  aux = gsub("+", " + ", aux, fixed = T)
  aux = gsub("-", " - ", aux, fixed = T)
  
  if ( length(aux) != 0 ){
    
    sif = data.frame(do.call(rbind, sapply(aux, strsplit, split = ' ', fixed = T )), stringsAsFactors = F)
    colnames(sif) = c('source', 'interaction', 'target')
    
    sif$interaction[sif$interaction=='+'] = '1'
    sif$interaction[sif$interaction=='-'] = '-1'
    
    writeLines(paste0(nrow(sif), " interactions found in at least ", round(psmpl), " samples out of ", length(colnames(topology))))
    
  }else{
    warning("There are no core interactions found")
    sif = NULL
  }
  
  
  return(sif)
  
}
