group_first <- function(lineMAT)
{
  ## from a matrix of line information given by the function infoLines() or infoLines_inser_del(),
     ## groups occurrences that are repetitions of the same sub-sequence Q
  ## this is the first step to count the number of occurrences, it will be followed by a second step to correct for possible
  ## under-estimations of the number of occurrences 
  ## thus, for now, we keep all groups, even those that represent sub-sequences repeated less than three times
  
  ## detects repetitions that happen before or after the "template" sub-sequence
  
  ## three outputs :
  # groups : a list that gives for each repeated sub-sequence a list of the "original" form and all different repeated forms
  # Tocc : a list that gives for each repeated sub-sequence a matrix with all the beginning and ending times of all its occurrences
  # nrep : a matrix that gives for each repeated sub-sequence the number of times it is repeated
   
  # this function should not be used without the second step (using the group() function)

  I <- unique(c(as.character(lineMAT$Q),as.character(lineMAT$Qrep))) # the vector of the different sub-sequences that are repeated
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
    
    ## remove all duplicates in Tocc and groups
    Tocc[[k]] <- Tocc[[k]][!duplicated(Tocc[[k]]),]
    Tocc[[k]] <- Tocc[[k]][order(Tocc[[k]][,1]),] 
    groups[[k]] <- unique(groups[[k]])
    
    ## count the number of occurrences
    ## be careful : if two occurrences are overlapping, count only one occurrence 
    ## (this can happen in special cases where shifted short sub-sequences create longer repeated sub-sequences)
    tocc <- Tocc[[k]]
    if(nrow(tocc)>1)
    {diff <- tocc[2:nrow(tocc),1]-tocc[1:(nrow(tocc)-1),2]
     n=1+length(diff[diff>0])}
    else{n=1}
    
    nrep <- rbind(nrep,data.frame(Q=i,n=n))
      
    k <- k+1
  }
     
  return(list(groups,Tocc,nrep))
}