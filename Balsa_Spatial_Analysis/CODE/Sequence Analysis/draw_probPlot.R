draw_probPlot <- function(nrep,ind_expl_final)
{
  ## this function draws the output of the function group()
  ## as a graph where each point represents a group of repeated sub-sequences
  ## abscissa : the length of the sub-sequence that is repeated
  ## ordinate : the number of times the sub-sequence is repeated
  
  # nrep : a matrix giving for every sub-sequence that is repeated, the number of times it is repeated
  # ind_expl_final : the indices of the sub-sequences that we retain after deleting the ones for which the repetitions are totally included within another group

  
  Q <- nrep[,1]
  Lengths <- sapply(Q,function(x) length(unlist(strsplit(as.character(x)," "))))
  nr <- nrep[,2]
  
  par(mfrow=c(1,1),mar=c(5,5,4,2)+0.1)
  plot(nr[ind_expl_final]~Lengths[ind_expl_final],col="black",pch=20,cex=3,cex.axis=1.5
        ,ylab="number of repetitions",xlab="length of the repeated sub-sequence",cex.lab=1.5,lwd=3
        ,xlim=c(min(Lengths[ind_expl_final])-1,max(Lengths[ind_expl_final])+1),ylim=c(min(nr[ind_expl_final])-1,max(nr[ind_expl_final])+1))
  
}
  