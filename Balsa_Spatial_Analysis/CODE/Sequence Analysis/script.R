

## example on the movement sequence corresponding to the last 1,500 time steps of the simulation for the sixth individual
# (represented in black in Fig. 3 in Riotte-Lambert et al. 2016)

SEQ <- read.table("SEQ.txt",sep="\t")
SEQ <- SEQ[,1]
patches= as.matrix(read.table("patches.txt",sep="\t"))

library(compiler)  ## to speed up the computations
enableJIT(3)

source("LOOK.R")
source("simp.R")
source("infoLines.R")
source("infoLines_inser_del.R")
source("group.R")
source("group_first.R")
source("infoLines_subMAT.R")
source("simp_subMAT.R")

source("inclusion.R")
source("incl_test.R")

source("draw_probPlot.R")
source("draw_seq.R")

source("entropy.R")
source("AutoO.R")


##---------------------------------------------------------------------------------------------------
## look for the repeated sub-sequences and make groups
##---------------------------------------------------------------------------------------------------

sortie <- LOOK(SEQ)

Tocc <- sortie[[2]]
lineMAT <- sortie[[4]]
nrep <- sortie[[3]]
groups <- sortie[[1]]

##-----------------------------------------------------------------------------------------
## Determine the indices of the groups that are contained within another group, 
## and that will be removed when plotting the output graph
##------------------------------------------------------------------------------------------

indices_remove <- inclusion(Tocc,lineMAT)
if(class(indices_remove)=="character") {
  ind_expl_final <- 1:length(Tocc)     # the indices of the groups that we retain
  nrep_corr <- nrep
}else{
  ind_expl_final <- (1:length(Tocc))[-indices_remove]  # the indices of the groups that we retain
  nrep_corr <- nrep[ind_expl_final,]
}



##--------------------------------------------------------------------------------------------------
## Draw the output graph
##--------------------------------------------------------------------------------------------------

draw_probPlot(nrep,ind_expl_final)

##-------------------------------------------------------------------------------------------------------------
# Plot specific repeated sub-sequences in space
##-------------------------------------------------------------------------------------------------------------

# to draw a specific plot
draw_seq(nrep,ind_expl_final,nr=4,L=15,sites_loc=patches,cex=1.5)


# to draw all the plots
Lengths <- sapply(nrep_corr[,1],function(x) length(unlist(strsplit(as.character(x)," "))))
for(i in unique(nrep_corr[,2])) {
  for(j in unique(Lengths[which(nrep_corr[,2]==i)])){
    draw_seq(nrep,ind_expl_final,nr=i,L=j,sites_loc=patches,cex=1.5)
  }
}

##-----------------------------------------------------------------------------------------
# calculate entropy and the routine movement index (R)
##--------------------------------------------------------------------------------

e = entropy(SEQ)

# draw the plot of corrected conditional entropies
plot(e,type="o",pch=20,ylim=c(0,1),xaxt="n",yaxs="i",xaxs="i",cex=1.2,xlab="order of dependency",ylab="corrected conditional entropy")
axis(1,at=c(1:6),labels=as.character(0:5))

# The routine movement index
R = 1-min(e)

# The relevant order of dependency O
O = AutoO(e)


