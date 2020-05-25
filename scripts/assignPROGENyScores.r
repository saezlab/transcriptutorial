assignPROGENyScores <- function (progeny = progeny, progenyMembers = progenyMembers, 
          id = "gene", access_idx = 1) 
{
  if (id == "uniprot") {
    idx <- which(names(progenyMembers) == "uniprot")
    progenyMembers <- progenyMembers[[idx]]
  }
  else {
    idx <- which(names(progenyMembers) == "gene")
    progenyMembers <- progenyMembers[[idx]]
  }
  members <- matrix(data = , nrow = 1, ncol = 2)
  pathways <- colnames(progeny)
  ctrl <- intersect(x = access_idx, y = 1:nrow(progeny))
  if (length(ctrl) == 0) {
    stop("The indeces you inserted do not correspond to \n              the number of rows/samples")
  }
  for (ii in 1:length(pathways)) {
    mm <- progenyMembers[[which(names(progenyMembers) == 
                                  pathways[ii])]]
    for (jj in 1:length(mm)) {
      members <- rbind(members, c(pathways[ii], mm[jj]))
    }
  }
  members <- members[-1, ]
  scores <- matrix(data = , nrow = nrow(progeny), ncol = nrow(members))
  colnames(scores) <- members[, 2]
  rownames(scores) <- rownames(progeny)
  members <- unique(members)
  for (i in 1:ncol(scores)) {
    for (j in 1:nrow(scores)) {
      scores[j, i] <- as.numeric(progeny[j, members[which(members[, 
                                                                  2] == colnames(scores)[i]), 1]])
    }
  }
  pxList <- list()
  for (ii in 1:length(access_idx)) {
    pxList[[length(pxList) + 1]] <- as.data.frame(t(as.matrix(scores[access_idx[ii], 
                                                                     ])))
  }
  names(pxList) <- rownames(progeny)[ctrl]
  return(pxList)
}