draw_seq <- function(nr,L,nrep,ind_expl_final,sites_loc,mvt_loc,cex,Plotrange)
{
  ## this function draws sub-sequences of length L repeated nr times in space
  ## this function is used as an exploratory tool to visualize in space the groups highlighted when drawing the graph nr~L with the function draw_probPlot()
  
  # takes as arguments : 
    # nr : the number of times the sub-sequence is repeated
    # L : the length of the sub-sequence that is repeated  
    # nrep : a matrix giving for every sub-sequence that is repeated, the number of times it is repeated
    # ind_expl_final : the indices of the sub-sequences that we retain after deleting the ones for which the repetitions are totally included within another group
    # sites_loc : the site locations
           # If there are only two columns, it is assumed that the identity of the site is given by the row index
           # If there are three columns, the first column contains the identity of the site
    # mvt_loc : the movement locations (optional argument)
    # cex : the "cex" argument for plotting the dots representing the sites
    # Plotrange : a matrix giving the range within which to draw the plots.
           # X range in the first row, Y range in the second row.
           # If it is not given, we will take the range of mvt_loc, and if mvt_loc is not given, the range of sites_loc

  Qrep <- nrep[,1]
  Lengths <- sapply(Qrep,function(x) length(unlist(strsplit(as.character(x)," "))))
  Nr <- nrep[,2]

  # indices of the sub-sequences that are of length L and repeated nr times
  indices <- intersect(which(Nr==nr),which(Lengths==L))
  indices <- intersect(indices,ind_expl_final)
  
  if(length(indices)==1)
  {par(mfrow=c(1,1))}
  if(length(indices)==2)
  {par(mfrow=c(1,2))}
  if(length(indices)>2)
  {par(mfrow=c(2,2))}
  
  if(ncol(sites_loc)==2)
  {
    for(i in indices)
    {
      Q <- as.numeric(strsplit(as.character(Qrep[i])," ")[[1]])
      
      if(!missing(mvt_loc))
      {
        if(missing(Plotrange))
          {plot(mvt_loc,type="l",asp=1,xlab="X",ylab="Y",xlim=range(mvt_loc[,1]),ylim=range(mvt_loc[,2]))}
        else{plot(mvt_loc,type="l",asp=1,xlab="X",ylab="Y",xlim=Plotrange[1,],ylim=Plotrange[2,])}
       points(sites_loc,pch=20,col="red",cex=cex)}
      else{
        if(missing(Plotrange))
          {plot(sites_loc,pch=20,col="red",cex=cex,asp=1,xlab="X",ylab="Y")}
        else{plot(sites_loc,pch=20,col="red",cex=cex,asp=1,xlab="X",ylab="Y",xlim=Plotrange[1,],ylim=Plotrange[2,])}
      }
      
      lines(sites_loc[Q,],type="o",pch=20,col="orange",cex=cex,lwd=3)
      points(sites_loc[Q[1],1],sites_loc[Q[1],2],pch=20,col="green",cex=cex)
      if(Q[1]!=Q[length(Q)])
        {points(sites_loc[Q[length(Q)],1],sites_loc[Q[length(Q)],2],pch=20,col="cyan",cex=cex)}
#       legend("topright",legend=paste(Q,collapse=", "))
    }
  }
  else
  {
    for(i in indices)
    {
      if(class(sites_loc[,1])=="numeric")
      {Q <- as.numeric(strsplit(as.character(Qrep[i])," ")[[1]])}
      else{Q <- strsplit(as.character(Qrep[i])," ")[[1]]}
      ind <- sapply(Q,function(x) which(sites_loc[,1]==x))
      
      if(!missing(mvt_loc))
      {
        if(missing(Plotrange))
          {plot(mvt_loc,type="l",asp=1,xlab="X",ylab="Y",xlim=range(mvt_loc[,1]),ylim=range(mvt_loc[,2]))}
        else{plot(mvt_loc,type="l",asp=1,xlab="X",ylab="Y",xlim=Plotrange[1,],ylim=Plotrange[2,],xlim=Plotrange[1,],ylim=Plotrange[2,])}
      points(sites_loc[,2:3],pch=20,col="red",cex=cex)}
      else{
        if(missing(Plotrange))
          {plot(sites_loc[,2:3],pch=20,col="red",cex=cex,asp=1,xlab="X",ylab="Y")}
        else{plot(sites_loc[,2:3],pch=20,col="red",cex=cex,asp=1,xlab="X",ylab="Y",xlim=Plotrange[1,],ylim=Plotrange[2,])}
      }
      
      lines(sites_loc[ind,2:3],type="o",pch=20,col="orange",cex=cex,lwd=3)
      points(sites_loc[ind[1],2],sites_loc[ind[1],3],pch=20,col="green",cex=cex)
      if(Q[1]!=Q[length(Q)])
      {points(sites_loc[ind[length(ind)],2],sites_loc[ind[length(ind)],3],pch=20,col="cyan",cex=cex)}
#       legend("topright",legend=paste(Q,collapse=", "))
    }
  }
  
  title(main=paste("n = ",as.character(nr),",  L = ",as.character(L),sep=""),outer=TRUE,line=-1.5)

}