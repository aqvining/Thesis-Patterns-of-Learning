simp <- function(line,lag)
{
  ## this function is used to authorize for the replacement of one site by another
  # line : a row of the site-lag matrix
  # lag : the corresponding lag (equal to the row number)
  ## the function replaces "0"s surrounded by four "1"s
  ## accepted situations : 11011, 10111, 11101 -> 11111
  
  line2 <- line
   
  if(length(line)-lag > 6)
  {
    if(line[1]==1
       &line[3]==1
       &line[4]==1
       &line[5]==1)
    {line2[2] <- 1}

    if((line[1]==1
        &line[2]==1
        &line[4]==1
        &line[5]==1)|(line[2]==1
                      &line[4]==1
                      &line[5]==1
                      &line[6]==1))
    {line2[3]<-1}
  
  
    if((length(line)-lag-3)>=4)
    {
      for (i in 4:(length(line)-lag-3))
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

  
    if(line[length(line)-lag-2]==1
       &line[length(line)-lag-3]==1
       &line[length(line)-lag-4]==1
       &line[length(line)-lag]==1)
    {line2[length(line)-lag-1] <- 1}
  
    if((line[length(line)-lag]==1
       &line[length(line)-lag-1]==1
       &line[length(line)-lag-3]==1
       &line[length(line)-lag-4]==1)|(line[length(line)-lag-1]==1
                                        &line[length(line)-lag-3]==1
                                        &line[length(line)-lag-4]==1
                                        &line[length(line)-lag-5]==1))
    {line2[length(line)-lag-2] <- 1}
  }


  return(line2)
}