group <- function(lineMAT)
{
  ## from a matrix of line information given by the function infoLines() or infoLines_inser_del(),
    ## groups sub-sequences that are repetitions of the same sub-sequence Q
  ## we keep groups that represent sub-sequences that are repeated at least 3 times.
  
  ## detects repetitions that happen before or after the focus sub-sequence
  
  ## three outputs :
        # groups : a list that gives for each repeated sub-sequence a list of the "original" form and all its different repeated forms
        # Tocc : a list that gives for each repeated sub-sequence a matrix with all the beginning and ending times of all its occurrences
        # nrep : a matrix that gives for each repeated sub-sequence the number of times it is repeated

  I <- unique(c(as.character(lineMAT$Q),as.character(lineMAT$Qrep))) # the list of all the different sub-sequences that are repeated
  groups <- NULL
  k <- 1
  Tocc <- NULL
  nrep <- NULL
  for (i in I)
  {
    sub <- rbind(lineMAT[lineMAT$Q==i,],
                 lineMAT[lineMAT$Qrep==i,])
    colnames(sub) <- c("lag","BE","END","Q","Qrep","BE","END")
    
    groups[[k]] <- c(as.character(i),
                      as.character(lineMAT$Qrep[lineMAT$Q==i]),
                     as.character(lineMAT$Q[lineMAT$Qrep==i]))
    
    Tocc[[k]] <- sub[,2:3]
    Tocc[[k]] <- rbind(Tocc[[k]],sub[,6:7])
    
    ## remove all replicates in Tocc and groups
    Tocc[[k]] <- Tocc[[k]][!duplicated(Tocc[[k]]),]
    Tocc[[k]] <- Tocc[[k]][order(Tocc[[k]][,1]),] # to order the occurrence times
    groups[[k]] <- unique(groups[[k]])
    
    ## count the number of occurrences
    ## be careful : if two occurrences are overlapping, count only one occurrence of the two
    tocc <- Tocc[[k]]
    n=1
    if(nrow(tocc)>1)
    {
      indice=1
      for(j in 2:nrow(tocc))
      {
       if((tocc[j,1]-tocc[indice,2])> 0)
       {n <- n+1
        indice <- j}
     }
    }
    
    ## keep the group only if the sub-sequence is repeated at least 3 times
    if(n>2)
    {
      nrep <- rbind(nrep,data.frame(Q=i,n=n))
      if(i!=I[length(I)])
      {
        k <- k+1
      }
    }
  }
  
  ## keep the last group only if the sub-sequence is repeated at least 3 times
  if(n<=2)
  {Tocc[[k]] <- NULL
  groups[[k]] <- NULL}
  
  return(list(groups,Tocc,nrep))
}