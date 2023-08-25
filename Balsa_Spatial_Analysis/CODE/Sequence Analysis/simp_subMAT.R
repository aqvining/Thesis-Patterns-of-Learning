simp_subMAT <- function(line)
{
  ## this function is used to authorize for the replacement of one site by another 
  ## in a row of a subMatrix (used for postdouble-check)
  # line : a row of the sub-site-lag matrix
  ## the function replaces "0"s surrounded by four "1"s
  ## accepted situations : 11011, 10111, 11101 -> 11111
  
  line2 <- line
   
  if(length(line) >= 5)
  {
    if(line[1]==1
       &line[3]==1
       &line[4]==1
       &line[5]==1)
    {line2[2] <- 1}
    
    if(length(line) >= 6)
    {
     if((line[1]==1
         &line[2]==1
         &line[4]==1
         &line[5]==1)|(line[2]==1
                       &line[4]==1
                       &line[5]==1
                       &line[6]==1))
     {line2[3]<-1}
     
     if((line[length(line)]==1
         &line[length(line)-1]==1
         &line[length(line)-3]==1
         &line[length(line)-4]==1)|(line[length(line)-1]==1
                                    &line[length(line)-3]==1
                                    &line[length(line)-4]==1
                                    &line[length(line)-5]==1))
     {line2[length(line)-2] <- 1}
    }
  
  
    if(length(line)>=7)
    {
      for (i in 4:(length(line)-3))
      {
        if((line[i-1]==1
            &line[i-2]==1
            &line[i+1]==1
            &line[i+2]==1)|(line[i-1]==1
                            &line[i+1]==1
                            &line[i+2]==1
                            &line[i+3]==1)|(line[i-3]==1
                                             &line[i-2]==1
                                             &line[i-1]==1
                                             &line[i+1]==1))
        {line2[i] <- 1}
      }
    }

  
    if(line[length(line)-2]==1
       &line[length(line)-3]==1
       &line[length(line)-4]==1
       &line[length(line)]==1)
    {line[length(line)-1] <- 1}
  
    
  }


  return(line2)
}