setwd("C:/Users/avining/Documents/Manuscripts/Do Primates Trapline")

source("./CODE/seperateTrials.R")
source("./CODE/DET_calculations.R")
source("./CODE/Array_Processing.R")
source("./CODE/Iterative Learning Simulations.R")
source("./CODE/Data_Processing_Functions.R")

data_raw <- read.csv("./DATA/Visits raw.csv")
data_formatted <- dataFormatting(data_raw)
write.csv(data_formatted,"./DATA/Visits_formatted.csv")
#data_formatted <- read.csv("./DATA/Visits_formatted.csv")[,-1]

data_expanded <- expandData(data_formatted) #puts each location in its own row. Allows analysis of DET by each decision (potentially)
data_expanded <- seperateTrials(data_expanded) #sets unique identifier between trials, prevents transitions between trials from affecting DET
write.csv(data_expanded,"./DATA/Visits Cleaned.csv")


DET_results <- getAllTrialDETs(data_expanded)
write.csv(DET_results, "./results/DET_by_trial.csv")
#DET_results <- read.csv("./results/DET_by_trial.csv")[,-1]

data_full <- full_join(data_formatted, DET_results)
data_full <- data_full[complete.cases(data_full$Sequence),] #removes trials from JK in 2017, which were coded redundantly with 2013 trials. This trials should ideal be removed prior to data processing

resource_arrays <- list(Double_Trapezoid, Pentagon, Z) #spatial feature data frames loaded from Array_Processing.R
names(resource_arrays) <- c("DT", "Pentagon", "Zarray")

data_full <- dplyr::mutate(data_full, Distance = mapply(getSequenceDistance, sequence = Sequence, resource_array = resource_arrays[Array])) #calculate total distance traveled for each trial

write.csv(runIterativeLearningSimulations(data_full, resource_arrays, runs = 100), "./Results/learning_simulations.csv") #run iterative learning simulations using the start probabilities of Individuals in empirical data and they arrays each Ind was tested in. Also calculates distance travelled in each simulated trial

data_learning_sim <- read.csv("./Results/learning_simulations.csv")[,-1]
data_learning_sim_expanded <- expandData(data_learning_sim[1:12000,])
data_learning_sim_seperated <- seperateTrials(data_learning_sim_expanded)
for(i in seq(12001, 348001, by = 12000)) data_learning_sim_seperated <- rbind(data_learning_sim_seperated, seperateTrials(expandData(data_learning_sim[i:(i+11999),]))) #unclear why this runs so much faster than doing it all at once, but it does. Splitting by 12000 keeps trials intact

write.csv(data_learning_sim_seperated, "./Results/data_learning_sim_seperated.csv")
data_learning_sim_seperated <- read.csv("./Results/data_learning_sim_seperated.csv")
data_learning_sim_seperated$Location <- as.character(data_learning_sim_seperated$Location)

#~~~~~~~~~~~~~~DET Calculation for simulations~~~~~~~~~~~~~~~#
#The following code calculates the similarity of each trial sequence to all previous trial sequences using the DET code from Ayers et al. 
#As this takes a long time (several days) and I have yet to set up a access to a remote processing server, it is broken into chunks based on when I needed to pause for various reasons
sequenceStarts <- c(0,which(data_learning_sim_seperated$Location == "1120")) + 1 #gets the first row of each Ind by Array sequence. Includes nrow + 1, which is used to handle the end case

#Chunk 1
learning_DET1 <- data.frame()
for (i in 1:length(sequenceStarts)) {
  print(paste("calculating DET for ", data_learning_sim_seperated$ID[sequenceStarts[i]], " in array ", data_learning_sim_seperated$Array[sequenceStarts[i]]))
  learning_DET1 <- rbind(learning_DET1, get_sequence_DET(data_learning_sim_seperated[sequenceStarts[i]:(sequenceStarts[i+1]-1),]))
}
learning_DET1 <- getAllTrialDETs(data_learning_sim_seperated)
write.csv(learning_DET1, "./Results/learning_sim_DET_partial.csv") #stopped and saved partial to clear up memory

#Chunk 3
learning_DET2 <- data.frame()
i #to start where left off, 1072
for (i in 1072:length(sequenceStarts)) {
  print(paste("calculating DET for ", data_learning_sim_seperated$ID[sequenceStarts[i]], " in array ", data_learning_sim_seperated$Array[sequenceStarts[i]]))
  learning_DET2 <- rbind(learning_DET2, get_sequence_DET(data_learning_sim_seperated[sequenceStarts[i]:(sequenceStarts[i+1]-1),]))
}
write.csv(learning_DET2, "./Results/learning_sim_DET_partial2.csv")
#Chunk 3
learning_DET3 <- data.frame()
i #to start where left off, 2336
for (i in 2336:length(sequenceStarts)) {
  print(paste("calculating DET for ", data_learning_sim_seperated$ID[sequenceStarts[i]], " in array ", data_learning_sim_seperated$Array[sequenceStarts[i]]))
  learning_DET3 <- rbind(learning_DET3, get_sequence_DET(data_learning_sim_seperated[sequenceStarts[i]:(sequenceStarts[i+1]-1),]))
}

#Chunk 4
write.csv(learning_DET3, "./Results/learning_sim_DET_partial3.csv")
i #2840
learning_DET4 <- data.frame()
for (i in 2840:length(sequenceStarts)) {
  print(paste("calculating DET for ", data_learning_sim_seperated$ID[sequenceStarts[i]], " in array ", data_learning_sim_seperated$Array[sequenceStarts[i]]))
  learning_DET4 <- rbind(learning_DET4, get_sequence_DET(data_learning_sim_seperated[sequenceStarts[i]:(sequenceStarts[i+1]-1),]))
}
write.csv(learning_DET4, "./Results/learning_sim_DET_partial4.csv")

#Combine
learning_DET1 <- read.csv("./Results/learning_sim_DET_partial.csv")[,-1]
learning_DET2 <- read.csv("./Results/learning_sim_DET_partial2.csv")[,-1]
learning_DET3 <- read.csv("./Results/learning_sim_DET_partial3.csv")[,-1]
learning_DET4 <- read.csv("./Results/learning_sim_DET_partial4.csv")[,-1]
learning_DET <- rbind(learning_DET1,learning_DET2, learning_DET3, learning_DET4)
write.csv(learning_DET, "./Results/learning_sim_DET_full.csv")

#~~~~~~~~~~~~~~~~End DET Calculation Code~~~~~~~~~~~~~~~~#



learning_DET <- fix_learning_DET(learning_DET) #adds species column and reformats ID column, function in Data_Cleaning.R
data_learning_sim <- fix_data_learning_sim(data_learning_sim) #reformats species and ID columns to match learning_DET, function in Data_Cleaning.R


learning_data <- full_join(data_learning_sim, learning_DET)
learning_data$ID <- giveUniqueID(learning_data$ID)

data_full <- rbind(data_full, learning_data)
data_full <- finalClean(data_full)
write.csv(data_full, "./Results/data_full.csv")

###~~~~~Data Processing Round 2, following initial visualizations. Goal is to look at cumulative DET, ie the total proportion of recursions that are in sequence repeats up to and including a given trial, rather than the proportion of recursions within the trial

