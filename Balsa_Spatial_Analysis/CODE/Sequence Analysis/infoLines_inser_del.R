infoLines_inser_del <- function(lineMAT,lineMAT2,SEQ,L)
{
  # this function joins lines of lineMAT2 with lines of lineMAT if it detects insertions or deletions
  
  # lineMAT is the matrix of information of lines of repetitions, output of the function infoLines() or infoLines_subMAT()
  # lineMAT2 is either equal to lineMAT, or an output of the function infoLines_inser_del if we authorize for several insertions or deletions
  # L is the minimal length of Q and Qrep (the sub-sequence that is repeated and its repeated form) for a deletion or an addition to be authorized
  
  # the output is a matrix of line information of the same kind of lineMAT and lineMAT2
  # the column "lag" is modified in order to authorize the re-employment of the output matrix to allow for several insertions or deletions
  
  lineMAT3 <- NULL
 
  if(nrow(lineMAT)>0
     &nrow(lineMAT2)>0)
  {
    rownames(lineMAT) <- 1:nrow(lineMAT)
    rownames(lineMAT2) <- 1:nrow(lineMAT2)
  
    indices2 <- 1:nrow(lineMAT2) # to remember the row indices of the lines of lineMAT2 that are not used in the junction process
    indices <- 1:nrow(lineMAT)   # to remember the row indices of the lines of lineMAT that are not used in the junction process
  
    for(i in 1:nrow(lineMAT2))  # scan all the lines of lineMAT2
    {

      #*** for an insertion ***#
    
      sub <- lineMAT[lineMAT$BE==(lineMAT2$END[i]+1)   
                      &lineMAT$lag==(lineMAT2$lag[i]+1),] 
      if(nrow(sub)>0)
      {
        ind <- as.numeric(row.names(sub))
        Q <- paste(lineMAT2$Q[i],paste(SEQ[lineMAT$BE[ind]:lineMAT$END[ind]],collapse=" "),sep=" ")
      
        if(length(unlist(strsplit(Q," "))) >= L) # the junction is done only if the resulting "original" sub-sequence Q is at least of length L
        {
          lineMAT3 <- rbind(lineMAT3,data.frame(lag=(lineMAT2$lag[i]+1),BE=lineMAT2$BE[i],END=lineMAT$END[ind],
                                     Q=Q,Qrep = paste(lineMAT2$Qrep[i]
                                                      ,as.character(SEQ[lineMAT$BE[ind]+lineMAT2$lag[i]])
                                                      ,lineMAT$Qrep[ind],sep=" ")
                                    ,BErep = lineMAT2$BErep[i]
                                    ,ENDrep = lineMAT$ENDrep[ind])) 
     
         # update the lists of indices that were not used in the junction process
         indices <- indices[indices!=ind]
         indices2 <- indices2[indices2!=i]
        }
      }
    
    
      #*** for a deletion ***#
    
      sub <- lineMAT[lineMAT$BE==(lineMAT2$END[i]+2)
                      &lineMAT$lag==(lineMAT2$lag[i]-1),] 
      if(nrow(sub)>0)
      { 
        ind = as.numeric(row.names(sub))
        Qrep = paste(lineMAT2$Qrep[i],lineMAT$Qrep[ind],sep=" ")
        Q = paste(SEQ[lineMAT2$BE[i]:lineMAT$END[ind]],collapse=" ")
        if(length(unlist(strsplit(Qrep," "))) >=L)        # the junction is done only of the resulting repeated form is at least L-site-long
        {
          lineMAT3 <- rbind(lineMAT3,data.frame(lag=(lineMAT2$lag[i]-1),BE=lineMAT2$BE[i],END=lineMAT$END[ind],
                                      Q=Q,Qrep=Qrep
                                      ,BErep = lineMAT2$BErep[i]
                                      ,ENDrep = lineMAT$ENDrep[ind]))
        
          # update the vectors of indices that were not used in the junction process
          indices <- indices[indices!=ind]
          indices2 <- indices2[indices2!=i]
        }
      }
    }
  
    lineMAT3 <- rbind(lineMAT3,lineMAT[indices,],lineMAT2[indices2,])
  }
 
  return(lineMAT3)
}