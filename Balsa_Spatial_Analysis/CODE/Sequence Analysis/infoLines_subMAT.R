infoLines_subMAT <- function(subMAT,t,SEQ)
{
  ## subMAT is a sub-division of the site-lag matrix for a specific occurrence occurring within the rank indices given by t (a vector of length 2)
  ## this function is used to look for all the occurrences of a specific sub-sequence occurring between the rank indices given by t
  ## the output is a matrix containing the information about all series of "1"s contained in subMAT
  
  ## in the output matrix "lineMAT", columns are :
  # lag : the lag at which the repetition occurs
  # BE : the rank indice at which the "original" sub-sequence begins
  # END : the rank indice at which the "original" sub-sequence ends
  # Q : the "original" sub-sequence
  # Qrep : the "repeated" sub-sequence
  # BErep : the rank indice at which the "repeated" sub-sequence begins
  # ENDrep : the rank indice at which the "repeated" sub-sequence ends
   
  
  P <- which(subMAT==1,arr.ind=TRUE) # Which indices in subMAT are "1"s
  P <- P[order(P[,1]),]
  t <- as.numeric(t)
  P[,2] <- P[,2] + t[1] -1  # correct the rank indices
  
  lineMAT <- NULL
  for(lag in unique(P[,1])) # lag stands for the lags at which some repetitions occur
  {
    times <- P[P[,1]==lag,2]   # The rank indices where repetitions occur at this lag
    
    
    if(length(times) > 0)
    {
      if(length(times)>=2)
      {
        diff <- times[2:length(times)]-times[1:(length(times)-1)]
        
        # if all the "1"s are consecutive
        if(max(diff)==1)  
        {
          lineMAT <- rbind(lineMAT,data.frame(lag=lag,BE=times[1],END=times[length(times)]
                                    ,Q=paste(SEQ[times[1]:times[length(times)]],collapse=" ")
                                    ,Qrep=paste(SEQ[(times[1]:times[length(times)])+lag],collapse=" ")
                                    ,BErep=times[1]+lag
                                    ,ENDrep=times[length(times)]+lag))}
        else{
          # if all the "1"s are not consecutive
          cuts <- which(diff>1)  # at which indices of times we need to "cut"
          
          lineMAT <- rbind(lineMAT, data.frame(lag=lag,BE=times[1],END=times[cuts[1]]
                                     ,Q=paste(SEQ[times[1]:times[cuts[1]]],collapse=" ")
                                     ,Qrep=paste(SEQ[(times[1]:times[cuts[1]])+lag],collapse=" ")
                                     ,BErep=times[1]+lag,ENDrep=times[cuts[1]]+lag))
          if(length(cuts)>1)
          {
            for(j in 2:length(cuts))
            {
              lineMAT <- rbind(lineMAT,data.frame(lag=lag,BE=times[cuts[j-1]+1],END=times[cuts[j]]
                                        ,Q=paste(SEQ[times[cuts[j-1]+1]:times[cuts[j]]],collapse=" ")
                                        ,Qrep=paste(SEQ[(times[cuts[j-1]+1]:times[cuts[j]])+lag],collapse=" ")
                                        ,BErep=times[cuts[j-1]+1]+lag,ENDrep=times[cuts[j]]+lag))
            }
          }
          lineMAT <- rbind(lineMAT,data.frame(lag=lag,BE=times[cuts[length(cuts)]+1],END=times[length(times)]
                                    ,Q=paste(SEQ[times[cuts[length(cuts)]+1]:times[length(times)]],collapse=" ")
                                    ,Qrep=paste(SEQ[(times[cuts[length(cuts)]+1]:times[length(times)])+lag],collapse=" ")
                                    ,BErep=times[cuts[length(cuts)]+1]+lag,ENDrep=times[length(times)]+lag))
        }
      }else{
        lineMAT <- rbind(lineMAT,data.frame(lag=lag,BE=times[1],END=times[1]
                                  ,Q=as.character(SEQ[times[1]])
                                  ,Qrep=as.character(SEQ[times[1]+lag])
                                  ,BErep=times[1]+lag,ENDrep=times[1]+lag))
      }
    }
  }
  
  rownames(lineMAT) <- 1:nrow(lineMAT)
  
  return(lineMAT)
}