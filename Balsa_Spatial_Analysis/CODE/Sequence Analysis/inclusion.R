inclusion <- function(Tocc,lineMAT)
{
  # looks for the groups of repetitions that are totally included within another group
    # the criterium for inclusion is : if all the ranks occupied by the repetitions of a sub-sequence are included within
    # the ranks occupied by the repetitions of another sub-sequence, then the first group of sub-sequences is deleted
  # if two groups are contained within each other, choose the one for which the "original" form is the most repeated with no modification
  # if the two groups are still equivalent, choose the one that has the longest "original" form
  # the output is the vector "indices_remove" of group indices that are to be removed
  
  # Tocc is the list of matrices giving, for each repeated sub-sequence, the beginning and ending ranks of its repetitions. It is an output of the function LOOK()
  # lineMAT is the raw information matrix on lines of 1s in the site-lag matrix (see help of the function infoLines for details)

  
  # build a matrix "indicesINC" that gives 
  # in the first column the index of the group that is included within the group of the indice given in the second column
  indicesINC <- NULL
  for(i in 1:length(Tocc))
  {
    for(j in (1:length(Tocc))[-i])
    {
      TEST <- sapply(1:nrow(Tocc[[i]]),function(x) incl_test(Tocc,i,j,x))
    
        if(sum(TEST)==nrow(Tocc[[i]]))
        {
          indicesINC <- rbind(indicesINC,c(i,j))
        }
    }
  }
  
  if(!is.null(indicesINC))
  {
    if(nrow(indicesINC)>1)
    {
      # look if there are groups included within each other
      mutual <- NULL
      for(i in 1:(nrow(indicesINC)-1))
      {
        reverse <- NULL
        for(j in (i+1):nrow(indicesINC))
        {reverse <- rbind(reverse,rev(indicesINC[j,]))}
    
        dup = which(duplicated(rbind(indicesINC[i,],reverse))) 
        if(length(dup)>0)
        {
         mutual <- rbind(mutual,indicesINC[i,]) 
        }
      }
  
      # if yes, investigate if the "original form" of one is more representative than the other
      I <- as.character(unique(lineMAT$Q))
      win <- NULL
      if(!is.null(mutual))
      {
        for(i in 1:nrow(mutual))
        {
          # check if one sub-sequence is repeated more often exactly than the other
          nrep1 <- length(lineMAT$Qrep[lineMAT$Qrep==I[mutual[i,1]]])
          nrep2 <- length(lineMAT$Qrep[lineMAT$Qrep==I[mutual[i,2]]])
          if(nrep1>nrep2)
          {win <- c(win,mutual[i,1])}
          if(nrep1<nrep2)
          {win <- c(win,mutual[i,2])}
     
          # if not, choose the group that has the longest "original" form
          if(nrep1==nrep2)
          {
            Nchar <- sapply(1:2,function(x) nchar(I[mutual[i,x]]))
            if(Nchar[1]==Nchar[2])
            {win <- c(win, sample(mutual[i,],1))}
            else
            {win <- c(win, mutual[i,which.max(Nchar)])}
          }  
        }
      }
  
      ## build the output vector of indices of groups that are to be removed
      indices_remove <- unique(indicesINC[,1])
      if(length(win)>0)
      {
        for(i in 1:length(win))
        {indices_remove <- indices_remove[indices_remove!=win[i]]}
      }
      return(indices_remove)
    }
    else
    {
      indices_remove <- indicesINC[1]
      return(indices_remove)
    }
  }
  else{print("no inclusion")}


}