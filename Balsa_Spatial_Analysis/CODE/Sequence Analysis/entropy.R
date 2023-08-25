entropy <- function(SEQ,max3=110,max4=58,max5=29)
{

  # this function computes for any sequence SEQ, its Shannon entropy, and the corrected conditional entropies at orders 1 to 5 following Porta et al. (1998)
  # the output is: 
          # O_corr: the vector of the Shannon entropy followed by the corrected conditional entropies for order 1 to 5
  
  # max4 and max5 are the maximum numbers of different sites visited in SEQ for the conditional entropies of order 4 and 5 to be computed, respectively
  # memory size problems can arise when computing the highest-order conditional transition matrices, for sequences with a high number of different sites used.
       # in this case, an error message is given, you can try to increase the amount of memory used by R using memory.limit()
       # otherwise, max4 and max5 should be set to smaller values
  # if the number of different sites visited is superior to max3, the computation is stopped and returns a series of NA
  
  # In any case, the relevant estimation of the minimum conditional entropy is warranted only if a minimum is reached in O_corr
  # that is, that the conditional entropy does not significantly decrease until the last value (that is, O5, or O4 or O3 if a very high number of different sites are visited)

  
  sites <- unique(SEQ)
  base <- length(sites)

  #-------------------------------
  ## order 0
  #-------------------------------

  prop <- sapply(sites,function(x) length(SEQ[SEQ==x])/length(SEQ))

  O0 <- -sum(sapply(1:base,function(x) prop[x]*log(prop[x],base)))

  #------------------------------
  ## order 1
  #------------------------------

  p1 <- matrix(nrow = base, ncol = base, 0)      # counts the transitions from one site to another
  for (t in 1:(length(SEQ) - 1)) p1[which(sites==SEQ[t]), which(sites==SEQ[t + 1])] <- p1[which(sites==SEQ[t]), which(sites==SEQ[t + 1])] + 1
  
  prop1 <- matrix(nrow = base, ncol = base, 0)   # conditional probabilities
  for (i in 1:base) 
    { if(sum(p1[i,])>0){prop1[i, ] <- p1[i, ] / sum(p1[i, ])   }}

  O1 <- 0
  for(i in 1:base)
  {
    for(j in 1:base)
    {
      if(prop1[i,j]>0)
      {O1 <- O1 - prop[i]*prop1[i,j]*log(prop1[i,j],base)}
    }
  }

  O1_corr <- O1 + (length(which(p1==1))/(length(SEQ)-1))*O0

  rm(prop1)

  #-------------------------------
  ## order 2
  #-------------------------------

  p2 <- matrix(nrow=base*base,ncol=base,0)
  for(t in 1:(length(SEQ)-2))
  {
    row_num <- (which(sites==SEQ[t])-1)*length(sites) + which(sites==SEQ[t+1])
    p2[row_num, which(sites==SEQ[t+2])] <- p2[row_num, which(sites==SEQ[t+2])]+1
  }
  
  prop2 <- matrix(nrow=base*base,ncol=base,0)
  for(i in 1:nrow(p2)) 
  {
    if(sum(p2[i,])>0)
    {prop2[i,] <- p2[i,]/sum(p2[i,])}
  }

  p1_glob <- p1/sum(p1)
  rm(p1)
  
  O2 <- 0
  for(i in 1:base)
  {
    for(j in 1:base)
    {
      row_num <- (i-1)*base + j
      for(k in 1:base)
      {
        if(prop2[row_num,k]>0)
        {O2 <- O2 - p1_glob[i,j]*prop2[row_num,k]*log(prop2[row_num,k],base)}
      }
    }
  }
  O2_corr <- O2 + (length(which(p2==1))/(length(SEQ)-2))*O0
  rm(prop2)
  rm(p1_glob)

  #----------------------------------
  ## order 3
  #----------------------------------

  if(base<=max3)
  {
    p3 <- matrix(nrow=base^3,ncol=base,0)
    for(t in 1:(length(SEQ)-3))
    {
      row_num <- (which(sites==SEQ[t])-1)*(base^2) + (which(sites==SEQ[t+1])-1)*base + which(sites==SEQ[t+2])
      p3[row_num, which(sites==SEQ[t+3])] <- p3[row_num, which(sites==SEQ[t+3])]+1
    }
  
    prop3 <- matrix(nrow=base^3,ncol=base,0)
    for(i in 1:nrow(p3)) 
    {
      if(sum(p3[i,])>0)
      {prop3[i,] <- p3[i,]/sum(p3[i,])}
    }

    p2_glob <- p2/sum(p2)
    rm(p2)

    O3 <- 0
    for(i in 1:base)
    {
      for(j in 1:base)
      {
        row_num_2 <- (i-1)*base + j
        for(k in 1:base)
        {
          row_num <- (i-1)*(base^2) + (j-1)*base + k
          for(l in 1:base)
          {
            if(prop3[row_num,l]>0)
            {O3 <- O3 - p2_glob[row_num_2,k]*prop3[row_num,l]*log(prop3[row_num,l],base)}
          }
        }
      }
    }
    O3_corr <- O3 + (length(which(p3==1))/(length(SEQ)-3))*O0
    rm(prop3)
    rm(p2_glob)


    #----------------------------------
    ## order 4
    #----------------------------------

    if(base <=max4)
    {
      p4 <- matrix(nrow=base^4,ncol=base,0)
      for(t in 1:(length(SEQ)-4))
      {
        row_num <- (which(sites==SEQ[t])-1)*(base^3) + (which(sites==SEQ[t+1])-1)*(base^2) + (which(sites==SEQ[t+2])-1)*base + which(sites==SEQ[t+3])
        p4[row_num, which(sites==SEQ[t+4])] <- p4[row_num, which(sites==SEQ[t+4])]+1
      }
    
      prop4 <- matrix(nrow=base^4,ncol=base,0)
      for(i in 1:nrow(p4)) 
      {
        if(sum(p4[i,])>0)
        {prop4[i,] <- p4[i,]/sum(p4[i,])}
      }

      p3_glob <- p3/sum(p3)
      rm(p3)

      O4 <- 0
      for(i in 1:base)
      {
        for(j in 1:base)
        { 
          for(k in 1:base)
          {
            row_num_3 <- (i-1)*(base^2) + (j-1)*base + k
            for(l in 1:base)
            {
              row_num <- (i-1)*(base^3) + (j-1)*(base^2) + (k-1)*base + l
              for(m in 1:base)
              {
                if(prop4[row_num,m]>0)
                {O4 <- O4 - p3_glob[row_num_3,l]*prop4[row_num,m]*log(prop4[row_num,m],base)}
              }
            }
          }
        }
      }
      O4_corr <- O4 + (length(which(p4==1))/(length(SEQ)-4))*O0
      rm(prop4)
      rm(p3_glob)


      #----------------------------------
      ## ordre 5
      #----------------------------------

      if(base <= max5)
      {
        p5 <- matrix(nrow=base^5,ncol=base,0)
        for(t in 1:(length(SEQ)-5))
        {
          row_num <- (which(sites==SEQ[t])-1)*(base^4) + (which(sites==SEQ[t+1])-1)*(base^3) + (which(sites==SEQ[t+2])-1)*(base^2) + (which(sites==SEQ[t+3])-1)*base + which(sites==SEQ[t+4])
          p5[row_num, which(sites==SEQ[t+5])] <- p5[row_num, which(sites==SEQ[t+5])]+1
        }
    
        prop5 <- matrix(nrow=base^5,ncol=base,0)
        for(i in 1:nrow(p5)) 
        {
          if(sum(p5[i,])>0)
          {prop5[i,] <- p5[i,]/sum(p5[i,])}
        }

        p4_glob <- p4/sum(p4)
        rm(p4)

        O5 <- 0
        for(i in 1:base)
        {
          for(j in 1:base)
          {
            for(k in 1:base)
            {
              for(l in 1:base)
              {
                row_num_4 <- (i-1)*(base^3) + (j-1)*(base^2) + (k-1)*base + l
                for(m in 1:base)
                {
                  row_num <- (i-1)*(base^4) + (j-1)*(base^3) + (k-1)*(base^2) + (l-1)*base + m
                  for(n in 1:base)
                  {
                  if(prop5[row_num,n]>0)
                    {O5 <- O5 - p4_glob[row_num_4,m]*prop5[row_num,n]*log(prop5[row_num,n],base)}
                  }
                }
              }
            }
          }
        }
        O5_corr <- O5 + (length(which(p5==1))/(length(SEQ)-5))*O0
  
        O_corr <- c(O0,O1_corr,O2_corr,O3_corr,O4_corr,O5_corr)
      }
      else
      {
        O_corr <- c(O0,O1_corr,O2_corr,O3_corr,O4_corr)
      }
    }
    else
    {
      O_corr <- c(O0,O1_corr,O2_corr,O3_corr)
    }
  }
  else
  {
    O_corr <- rep(NA, times=6)
  }
    return(O_corr)

}
