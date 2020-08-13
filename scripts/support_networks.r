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
                            
                            r = c("node" = n, "total_degree" = in_degree + out_degree,
                                  "in_degree" = in_degree, "in_positive" = in_positive, "in_negative" = in_negative,
                                  "out_degree" = out_degree, "out_positive" = out_positive, "out_negative" = out_negative
                            )
                            return(r)
                            
                          })
                          )
                  )
  return(dd)
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
  
  counts = data.frame(
          do.call(rbind,
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
  
  
  
  return(counts)
}
