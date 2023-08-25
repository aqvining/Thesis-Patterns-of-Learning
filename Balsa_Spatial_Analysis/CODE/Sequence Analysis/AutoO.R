AutoO <- function(en,thr=0.05)
{
  # this function looks for the most relevant order of dependency
  # under the criteria: O is the most relevant order of dependency
  # if it is the largest value of p leading to a significant decrease in entropy
  # H_p-1 - H_p < thr
  
  # en is the vector of corrected values of conditional entropies
  # thr is the threshold value to consider that a decrease in corrected conditional entropy is significant
  
  
  O1 = which.min(en)

  
  if(O1==1)
  {return(O1-1)}
  else
  {
    O=O1
    for(i in O1:2)
    {
      if(en[i-1]-en[i]<thr)
      {
        O=O-1
      }
      else{break}
    }
    return(O-1)
  }
}