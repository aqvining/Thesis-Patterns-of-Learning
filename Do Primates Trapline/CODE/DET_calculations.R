library(spatstat)
filterout <- function(Ldata){
  for (i in 2:length(Ldata)){
    if(Ldata[i] == Ldata[i-1] ){Ldata[i - 1]= NA}
  }
  Ldata=Ldata[!is.na(Ldata)]
  Ldata
}

getAllTrialDETs <- function(data_full, minl = 3) {
  #purpose: transform a dataframe of sequences location visits into a dataframe of trial by trial DET values
  #input: datafull; a dataframe with the columns "Array" "ID" "Location" and "Trial"
  #output: A dataframe with the same structure as input, but with DET and recursions columns instead of Location and one row per trial 
  results <- data.frame()
  data_full$Location <- as.numeric(data_full$Location)
  for (array in unique(data_full$Array)){
    array_data <- filter(data_full, Array == array)
    for (id in unique(array_data$ID)) {
      print(paste("processing individual", id, "in array", array))
      sequence_data <- filter(array_data, ID ==id)
      if (max(sequence_data$Location) < 1002) {
        warning(paste("only one trial completed, data from", id, "in", array, "omitted"))
      } else {
        trialDETs <- get_trapline_recursions(sequence_data$Location, minl)
        trialDETs$ID <- id
        trialDETs$Array <- array
        results <- rbind(results, trialDETs)
      }
    }
  }
  return(results)
}

get_sequence_DET <- function(sequence_data, minl = 3) {
  #purpose: Get DET calculations from a a dataframe with location sequences
  #input:sequence_data; A dataframe with columns Location (character), ID (Charactor or Factor), Array (Character or Factor), and Trial (numeric)
  #      minl; the minimum length for a sequence match to be considered a repeat (numeric)
  #output: A dataframe with one row per unique Trial, and columns DET, Array, ID, and Trial
  if(length(unique(sequence_data$ID)) != 1 | length(unique(sequence_data$Array)) != 1) stop("This function is only meant to process a sequence from a single individual in a single array")
  trialDETs <- get_trapline_recursions(sequence_data$Location, minl)
  trialDETs$ID <- unique(sequence_data$ID)
  trialDETs$Array <- unique(sequence_data$Array)
  return(trialDETs)
}

get_trapline_recursions <- function(x, minl){
  require(dplyr)
  #purpose: determine the determinism of each trial in a sequence of location visits separated by trials. 
  #input: x; a sequence of locations, inferable as numerics. Trials should be separated in the sequence by the value 1000 + the preceding trial number
  #notes: The code in this function is adapted from Ayers et al. 2015
  
  x = as.numeric(x)
  #Depending on the dataset it may be desirable to filter out consecutive visits 
  #to the same flower. See function below and delete '#' in the line below to use
  x = filterout(Ldata = x)
  
  #-----set up matrix resembling a recurrence plot, where a 1 indicates a repeat 
  #-----visit and 0 indicates the absence of a repeat.
  #if (length(unique(x)) < minl) return(NA)
  #if (length(x) <= 3*minl) return(NA)
  #if (length(unique(x)) == length(x)) return(0)
  det1 = matrix(cbind(rep(x,length(x))),nrow=length(x))
  tdet = t(det1)
  det = ((det1 - tdet) == 0) * 1
  
  #set the main diagonal equal to zero so it won't be included in the calculation
  
  diag(det) = 0
  
  #Use spatstat package to create a 'countour map' of the matrix,
  #which assigns all sets of contiguous 1's a unique number
  yi <- as.im(det)
  ycOut <- connected(yi, background = 0)
  yc <- ycOut$v
  
  #Depending on the dataset it may be desirable to filter out diagonals perpendicular to #the main diagonal. Code is provided for the 'removeperpdiag' function below.
  #Delete "#" from the line below to filter out perpendicular diagonals
  
  #yc = removeperpdiag(yc,minl)
  
  #Note: this code may take several minutes to run for very long sequences
  
  #---- filter out short repeats: a 'trapline' should include more unique resources
  #---- than the minimum cutoff (minl)
  
  #make an alternative DET matrix that contains the resource IDs
  det2 = matrix(rep(x,nrow(det)),nrow=nrow(det),byrow=TRUE)*det
  #make a dataframe with the number of times each resource appears in a diagonal
  listofseq = data.frame(group = yc[1:length(yc)], seq=det2[1:length(det2)])
  #how many unique resources are in the diagonal
  uniquevisits = rowSums((table(listofseq)>0)*1)
  #only count diagonals with at least 'minl' number of unique resources
  longenough = (uniquevisits >= minl)*table(yc)
  #ycOutlongenough <- ycOut
  binaryDET <- data.frame(matrix(0, nrow = sum(det)/2, ncol = 2))
  colnames(binaryDET) <- c("seqRepeat", "Trial")
  recursionNum = 1 #tracks number of recursions that have been found in following for loop (should find nrow(binaryDET) recursions), AQV
  for (i in 1:nrow(ycOut$v)) { #for each visit, AQV
    for (j in 1:i){ #and for each visit that occurred previous to it, AQV
      group = ycOut$v[j,i] #get the unique identifier of the sequence of locations both recursions belong to (if they are not in matching sequences, this will be NA), AQV
      if(!is.na(group)) {
        if (!longenough[[group]] == 0) binaryDET$seqRepeat[recursionNum] <- 1 #if the sequence is longer than minL, this recursion is a repeat, AQV
        binaryDET$Trial[recursionNum] <- sum(i > which(x >=100)) + 1 #counts the number of trial end markers that occurred in sequence prior to current recurrence and adds one to get the current trial number, AQV
        recursionNum <- recursionNum + 1
      }
    }
  }
  #find the numerator:
  #(remember this still includes both the top and bottom halves of the matrix)
  
  #contig = sum(longenough)
  
  #denominator= sum(det)
  
  #This also still includes top and bottom halves of the matrix
  
  
  #------------------- total DET score
  #divide the numerator and denominator in half before calculating DET for just
  #the top half of the matrix
  summary <- binaryDET %>% group_by(Trial) %>% summarise(DET = mean(seqRepeat), Recursions = length(seqRepeat)) #get the number of recursions and the proportion that were repeats for each trial in the sequence 
  return(summary) 
}
