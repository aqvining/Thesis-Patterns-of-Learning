seperateTrials <- function(seqData) {
  #purpose: insert rows between all trials with the location column set to 1000 + trial number. Prevents transitions between trials from influencing DET calculations
  #input: seqData; a dataframe with columns "Location" and "Order". Typically produced by the data_cleaning_DET and/or expand_data functions.
  #output: a dataframe of the same structure as seqData, with additional rows seperating trials.
  seperated <- data.frame()
  if (sum(!complete.cases(seqData$Location)) > 0) {
    warning(paste("incomplete cases removed, rows", paste(which(!complete.cases(seqData)), collapse = ",")))
    seqData <- seqData[complete.cases(seqData),]
  }
  seqData$Location <- as.character(seqData$Location)
  while (nrow(seqData) > 0) {
    trial_data <- filter(seqData, ID == seqData$ID[1], Array == seqData$Array[1], Trial == seqData$Trial[1]) #get data from first trial in data frame
    seqData <- seqData[seqData$ID != seqData$ID[1] | seqData$Array != seqData$Array[1] | seqData$Trial != seqData$Trial[1],] #remove trial data (currently in processing) from seqData
    trial_data <- rbind(trial_data, trial_data[1,]) #add new row to trial data
    trial_data$Location[nrow(trial_data)] <- 1000 + trial_data$Trial[1] #sets location to a unique identifer based on trial number to seperate movements
    trial_data$Order[nrow(trial_data)] <- nrow(trial_data)
    seperated <- rbind(seperated, trial_data)
    print(nrow(seqData))
  }
  return(seperated)
}
