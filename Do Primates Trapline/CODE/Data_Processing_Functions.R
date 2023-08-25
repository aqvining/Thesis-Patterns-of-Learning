require(dplyr)

data_cleaning_DET <- function(raw_data) {
  #script to read, clean, and expand location visit data of primates navigating local arrays for analysis of DET
  #inputs: none, check code to ensure data file and path are correct
  #output: Fully cleaned and expanded vertical dataframe with a row for each location visited across all experiments
  raw_data <- read.csv("./DATA/Visits Raw.csv")
  data_formatted <- dataFormatting(raw_data)
  data_formatted$Sequence <- as.character(data_formatted$Sequence)
  data_cleaned <- expandData(data_formatted)
  return(data_cleaned)
}

dataFormatting <- function(data) {
  data <- data[,1:5]
  colnames(data)[4] = "Sequence"
  data$Trial <- NA
  
  #add trial numbers to data~~~~#
  data$Trial[1] <- 1
  trial <- 1
  for (i in 2:nrow(data)) {
    if (data$ID[i] == data$ID[i-1] & data$Array[i] == data$Array[i-1]) { trial <- trial + 1
    } else trial <- 1
    data$Trial[i] <- trial
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  return(data)
}

expandData <- function(data) {
  #purpose: apply split_trial to each row of the data, tracking the trial number for each individual x Array
  #input: dataframe with cols for Array, Species, ID, Sequence, and Session
  #output: dataframe with 1 row per location visited across trials in Location Column. Trial col becomes numeric, giving the trial number a visit occurred, and Order gives the sequence of visits within a trial
  data$Sequence <- as.character(data$Sequence)
  expanded_data <- data.frame(matrix(nrow = sum(sapply(data$Sequence, FUN = function(X) length(strsplit(X, NULL)[[1]]))), ncol = 7)) #to hold all split data
  colnames(expanded_data) <- c("Array", "Species", "ID", "Location", "Order", "Trial", "Session")
  j <- 1
  #group data frame by array and ID
  for (array in unique(data$Array)) {
    array_data <- filter(data, Array == array)
    for (id in unique(array_data$ID)){
      IND_data <- filter(array_data, ID == id)
      print(paste(array, id))
      #split and order characters from each sequence
      for (i in 1:nrow(IND_data)){
        trial_data <- as.matrix(split_trial(trial_data = IND_data[i,], trial_num = i)) #as.matrix converts all values to char, easier to convert back than deal with dataframe formatting
        expanded_data[j:(j+nrow(trial_data) -1),] <- trial_data
        j <- j + nrow(trial_data)
      }
    }
  }
  expanded_data$Order <- as.numeric(expanded_data$Order)
  expanded_data$Trial <- as.numeric(expanded_data$Trial)
  expanded_data$Session <- as.numeric(expanded_data$Session)
  expanded_data <- expanded_data[complete.cases(expanded_data$Location),] #removal of "S" and "F" by split_trial can leave empty rows at the end
  return(expanded_data)
}

split_trial <- function(trial_data, trial_num) {
  #purpose: split a single row of data w/ all platform visits ordered in a single column into a data frame with a row for each platform visit
  #inputs: trial_data; A 1 row data frame with cols Array, Species, ID, Sequence, Session, and Dist
  #~        trial_num; integer
  #output: A dataframe with nrow = # characters in trial column of input. Added columns for trial and w/in trial order
  visits <- strsplit(as.character(trial_data$Sequence[1]), NULL)[[1]] #create vector with element for each visit
  visits <- visits[!visits == "F" & !visits == "S"] # S and F were used to encode trial starts and ends in some cases but not others. Removes this artefact
  expanded_data <- data.frame(Array = trial_data$Array[1], 
                              Species = trial_data$Species[1], 
                              ID = trial_data$ID[1],
                              Location = visits,
                              Order = 1:length(visits),
                              Trial = trial_num,
                              Session = trial_data$Session,
                              stringsAsFactors = FALSE)
  return(expanded_data)
}

fix_learning_DET <- function(learning_DET) {#Adds species (learing factor) to learning_DET, which was not done during calculations (oops)
  #rename ID to make splitting species from number easier
  learning_DET$ID <- sapply(learning_DET$ID, function(X) gsub("Aye aye", "Aye_aye", X))
  learning_DET$ID <- sapply(learning_DET$ID, function(X) gsub("Mouse Lemur", "Mouse_Lemur", X))
  learning_DET$ID <- sapply(learning_DET$ID, function(X) gsub("Dwarf Lemur", "Dwarf_Lemur", X))
  learning_DET$ID <- sapply(learning_DET$ID, function(X) gsub("Japanese macaque", "Japanese_macaque", X))
  learning_DET$Species <- NA # to be filled in subsequent for loop
  learning_factor <- 1 #default for first row, will update in subsequent for loop
  for (i in 1:nrow(learning_DET)) {
    ID <- strsplit(learning_DET$ID[i], split = " ")[[1]] #creates character vector with species in element 1 and ID number in element 2
    #following conditionals check if simulations have started on a new learning factor, using the transitions from ID 100 to ID 1
    if (i != 1 & ID[2] == 1) {
      if (learning_DET$ID[i - 1] == 100) {
        learning_factor <- learning_factor %% 3 + 1 #adds 1 to learning factor, resetting to 1 after 3
        print(ID)
      }
    }
    learning_DET$Species[i] <- paste(ID[1], " Learning_Factor", c("1", "1.2", "2")[learning_factor], sep = "")
    learning_DET$ID[i] <- ID[2]
  }
  return(learning_DET)
}

fix_data_learning_sim <- function(data) {#modifies Species column to include learning factor AND biological species starting locations are drawn from. Removes the latter element from ID column
  data$ID <- sapply(data$ID, function(X) gsub("Aye aye", "Aye_aye", X))
  data$ID <- sapply(data$ID, function(X) gsub("Mouse Lemur", "Mouse_Lemur", X))
  data$ID <- sapply(data$ID, function(X) gsub("Dwarf Lemur", "Dwarf_Lemur", X))
  data$ID <- sapply(data$ID, function(X) gsub("Japanese macaque", "Japanese_macaque", X))
  data$Species <- as.character(data$Species)
  for(i in 1:nrow(data)){
    ID <- strsplit(data$ID[i], split = " ")[[1]] #gets character vector with species name in first element and individual name in second
    learning_factor <- strsplit(data$Species[i], split = " ")[[1]][2] #gets factor number from full name
    data$Species[i] <- paste(ID[1], " Learning_Factor", learning_factor, sep = "")
    data$ID[i] <- ID[2]
    if (i %% 1000 == 0) print(data$Species[i])
  }
  return(data)
}

giveUniqueID <- function(IDs) {
  #input: A character vector, or coercable as such.
  #output: A character vector
  #Purpose: for each element in the input vector, the output reflects how many changes in value have occurred to it prior in the sequence. Used to give unique IDs to different individuals with the same name
  newID <- 1
  for(i in 1:(length(IDs)-1)){
    ID <- IDs[i]
    IDs[i] <- newID
    if (IDs[i+1] != ID) newID <- newID + 1
  }
  IDs[length(IDs)] <- newID
  return(as.character(IDs))
}

finalClean <- function(data_full) {
  #input: dataframe with column "Species"
  #output: dataframe with same structure as input plus column "Source"
  #purpose: A final standardization of name formats, and a split of simulation parameters from species for easier analysis
  data_full$Species <- gsub("Aye aye", replacement = "Aye_Aye", x = data_full$Species)
  data_full$Species <- gsub("Aye_aye", replacement = "Aye_Aye", x = data_full$Species)
  data_full$Species <- gsub("Mouse Lemur", replacement = "Mouse_Lemur", x = data_full$Species)
  data_full$Species <- gsub("Dwarf Lemur", replacement = "Dwarf_Lemur", x = data_full$Species)
  data_full$Species <- gsub("Japanese macaque", replacement = "Japanese_Macaque", x = data_full$Species)
  data_full$Species <- gsub("Japanese_macaque", replacement = "Japanese_Macaque", x = data_full$Species)
  data_full$Species <- gsub("NA", "None", data_full$Species)
  data_full$Source <- "Experimental"
  for(i in 1:nrow(data_full)){
    rename <- strsplit(data_full$Species[i], split = " ")[[1]]
    data_full$Species[i] <- rename[1]
    if(length(rename) == 2) data_full$Source[i] <- rename[2]
    if(length(rename) > 2) stop("strsplit produced too many elements")
  }
  return(data_full)
}
