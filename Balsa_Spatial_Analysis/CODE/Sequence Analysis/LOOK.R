LOOK <- function(SEQ,mutation_corr=TRUE,inserdel=TRUE,L=5)
{
  ## this function performs all the steps to look for repeated sub-sequences in SEQ
 
  # SEQ is the movement sequence, i.e. a sequence of numbers representing the ordered visitation of different sites
  
  # if mutation_corr=TRUE, mutations are authorized given the rule : 10111, 11011, 11101 -> 11111,
    #  with "1" representing an identical site, "0" a different site.
  # if inserdel=TRUE, one insertion or deletion is authorized per sub-sequence that is more than L-site-long
  # the function performs a double check for repeated sub-sequences. This is necessary to avoid the under-estimation of the number of repetitions
  
  # the output repeated sub-sequences are at least 3-site-long and repeated at least 3 times
  
  # the outputs are :
          # groups : a list of vectors giving, for each repeated sub-sequence, all of its different repeated forms
          # Tocc : a list of matrices that give, for each repeated sub-sequence, the time indices of beginnings and endings of its occurrences
          # nrep : a matrix that gives for each repeated sub-sequence (specified in the first column) the number of times it is repeated (in the second column)
          # lineMAT_final : the raw information matrix on lines of 1s in the site-lag matrix (see help of the function infoLines for details)
  
  
  ##----------------------------------------------------------------------------
  ## Construct the similarity matrix
  ##----------------------------------------------------------------------------
  
  MAT <- matrix(0,nrow=length(SEQ),ncol=length(SEQ))
  for(i in 1:(length(SEQ)-1))
  {
    id <- which(SEQ[1:(length(SEQ)-i)]==SEQ[(1+i):length(SEQ)])
    MAT[i,id] <- 1
  }
  
  if(mutation_corr==TRUE)
  {
   # Allow for the replacement of one site by another
   # situations authorized : 10111, 11011, 11101 -> 11111
   MATsimp <- t(sapply(1:length(SEQ),function(x) simp(MAT[x,],x)))
  }else{MATsimp <- MAT}
  
  ##------------------------------------------------------------------------------------
  ## Construct the "line information" matrix "lineMAT"
  ## optionally, authorize for an insertion and/or deletion per repeated sub-sequence
  ##------------------------------------------------------------------------------------
  
  lineMAT <- infoLines(MATsimp,SEQ)
  
  if(inserdel==TRUE)
  {
    # Correct lineMAT for insertions and deletions
    lineMAT2 <- infoLines_inser_del(lineMAT,lineMAT,SEQ,L)
  }else{lineMAT2 <- lineMAT}
  
  ## delete sub-sequences in lineMAT2 that are less than 3-sites-long
  lineMAT2 <- lineMAT2[(lineMAT2$END-lineMAT2$BE)>1,]
  lineMAT2 <- lineMAT2[(lineMAT2$ENDrep-lineMAT2$BErep)>1,]
  
  ##-------------------------------------------------------------------------------------------------------
  ## Group the repetitions by common sub-sequence that is repeated, 
  ## and correct for possible under-estimations of the number of occurrences
  ##-------------------------------------------------------------------------------------------------------
  
  g1 <- group_first(lineMAT2)  # construct the groups from lineMAT2 and keeps all of them, 
                               # even the ones repeated representing sub-sequences repeated less than 3 times
  
  # look in MATsimp if we missed some repetitions
  lineMAT_final <- NULL
  for(ind in 1:length(g1[[1]]))  # for every repeated sub-sequence
  {
    tocc <- g1[[2]][[ind]]       # the matrix of the ranks of occurrences of the indth repeated sub-sequence
    for(j in 1:nrow(tocc))
    {
      t <- as.numeric(tocc[j,])
      subMAT <- MAT[,t[1]:t[2]]
      if(mutation_corr==TRUE)
      {
        subMATsimp <- t(sapply(1:length(SEQ),function(x) simp_subMAT(subMAT[x,])))
      }else{subMATsimp <- subMAT}
        
      if(length(which(subMAT==1))>=3)
      {
        lineMAT_bis <- infoLines_subMAT(subMATsimp,t,SEQ)
        if(inserdel==TRUE)
        {
          lineMAT2_bis <- infoLines_inser_del(lineMAT_bis,lineMAT_bis,SEQ,L)
        }else{lineMAT2_bis <- lineMAT_bis}
          
        ## delete sub-sequences that are less than 3-sites-long
        lineMAT2_bis <- lineMAT2_bis[(lineMAT2_bis$END-lineMAT2_bis$BE)>1,]
        lineMAT2_bis <- lineMAT2_bis[(lineMAT2_bis$ENDrep-lineMAT2_bis$BErep)>1,]
        
        if(nrow(lineMAT2_bis)>0)
        {lineMAT_final <- rbind(lineMAT_final,lineMAT2_bis)}
      }
    }
  }
  g2 <- group(lineMAT_final)  # corrected groups
  
  groups <- g2[[1]]  # the different forms of the sub-sequence that is repeated within each group
  Tocc <- g2[[2]]    # all the beginning and ending times of all the occurrences within each group
  nrep <- g2[[3]]    # summary matrix with for each sub-sequence that is repeated, the number of times it is repeated 
  
  return(list(groups,Tocc,nrep,lineMAT_final))
}