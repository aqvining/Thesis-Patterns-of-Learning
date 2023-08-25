incl_test <- function(Tocc,i,j,k)
{
  # this function tests if the time indices of the kth line of Tocc[[i]]
  # are explained by an occurrence within the jth group
  
  bool <- 0
  for(l in 1:nrow(Tocc[[j]]))
  {
    if(Tocc[[i]][k,1] >= Tocc[[j]][l,1]
       &Tocc[[i]][k,2] <= Tocc[[j]][l,2])
    {bool <- 1}
  }
  return(bool)
}
