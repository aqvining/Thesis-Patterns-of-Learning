source("./CODE/Array_Processing.R")

runIterativeLearningSimulations <- function(emperical_data, resource_arrays, runs = 1000){
  #purpose: Parameterize and run all iterative learning simulations
  #input: emperical_data; data_frame with empirical observations
  #       resource_arrays: named list containing spatial features dataframes
  simulated_data <- data.frame()
  for(array in names(resource_arrays)) {
    for (species in c("NA",unique(as.character(filter(emperical_data, Array == array)$Species)))) {
      for (learning_factor in c(1,1.2,2)) {
        print(paste("running ", species, " in ", array, " with learning factor", learning_factor))
        for (i in 1:runs) {
          if (i %% 50 == 0) print(paste("run ",i, " of ", runs))
          run_data <- iterativeLearningSimulation(resource_arrays[array], learning_factor, trials = 120, starts = getStarts(resource_arrays[array], species, emperical_data))
          run_data$ID <- paste(species, i)
          simulated_data <- rbind(simulated_data, run_data)
        }
      }
    }
  }
  return(simulated_data)
}

iterativeLearningSimulation <- function(resource_array, learning_factor, trials, exclude_backtrack = TRUE, starts = "NA") {
  #purpose: Run a set of trials in which a foraging agent moves between locations in a resource array, minimizing total length traveled by increasing the probability of transitions in a trial that resulted in the shortest path traveled
  #input: resource_arrays; a named list of length one containing a spatial features dataframe with points in the geometry column and a feature NAME
  #       learning_factor; a numeric giving the factor to multiply transition probabilities by upon completion of a shortest-length-traveled trial
  #       trials: the number of trials (full array visits) to simulate
  #       starts: a named vector containing the probability that each element name will be the first location visited. if NULL, starting locations will be selected at random from resource_array$NAME
  #output: a dataframe with colums for Array (character), Species (character, "itLearner_"), Sequence (character), Session (numeric, NA), Trial (numeric), and Distance (Numeric)
  learner1 <- Learner(sequence = "S", transition_matrix = getDistanceWeightedTransitions(resource_array, starts), learning_factor = learning_factor, min_distance_trial = 0)
  data <- data.frame(Array = names(resource_array), Species = paste("itLearner_", learning_factor), Sequence = NA, Session = NA, Trial = 1:trials, Distance = NA)
  for(i in 1:trials) {
    learner1$runTrial(exclude_backtrack)
    data$Sequence[i] <- paste0(learner1$sequence, collapse = "")
    data$Distance[i] <- getSequenceDistance(learner1$sequence, resource_array)
    if(learner1$min_distance_trial == 0) learner1$updateTransitions(data$Distance[i])
    if(data$Distance[i] <= learner1$min_distance_trial) learner1$updateTransitions(data$Distance[i])
  }
  return(data)
}

getStarts <- function(resource_array, species, data) {
  #purpose: get proportion of observed starts at each location in a given array for a given species
  #input: resource_array; a named list of length one contain a spatial features data frame with spatial points in the geometry column and feature NAME
  #       species; character matching level in data$Species
  #       data; a dataframe with column Species (factor), Array (factor), and Sequence (vector of characters seen in resource_array$NAME OR a single numeric/character that can be split into such)
  #output: a single character giving a NAME from resource_array
  if (species == "NA") {start_probs <- rep(1/nrow(resource_array[[1]]), times = nrow(resource_array[[1]]))
    names(start_probs) <- resource_array[[1]]$NAME
    return(start_probs)
  }
  data <- filter(data, Array %in% names(resource_array), Species == species)
  starts <- rep(NA, times = nrow(data))
  for(i in 1:length(starts)) {
    sequence <- data$Sequence[i]
    if (length(sequence) == 1) sequence <- strsplit(as.character(sequence), split = NULL)[[1]] #split character by element
    starts[i] <- sequence[1]
  }
  start_probs <- sapply(resource_array[[1]]$NAME, FUN = function(X, starts) sum(starts == X)/length(starts), starts = starts) #gets proportion of times each location in array was the first in a sequence
  names(start_probs) <- resource_array[[1]]$NAME
  return(start_probs)
}

getDistanceWeightedTransitions <- function(resource_array, starts) {
  #purpose: generate a transition matrix based on distances between resource arrays, with likelihood proportional to 1/d^2
  #input: resource_array; a spatial features data frame with spatial points in the geometry column and NAME feature OR a list of such (only the first element will be used)
  #output: numeric array
  if("list" %in% class(resource_array)) resource_array <- resource_array[[1]]
  distanceMatrix <- getLabelledDistanceMatrix(resource_array) #from Array_Processing.R
  distanceMatrix <- apply(distanceMatrix, MARGIN = c(1,2), FUN = function(X) ifelse(X==0, 0, 1/(X^2))) #transitions proportional to 1/distance^2
  transition_matrix <- t(apply(distanceMatrix, MARGIN = 1, function(X) X/(sum(X)))) #normalize, transposition necessary as apply output is by column
  transition_matrix <- rbind(starts, transition_matrix)
  rownames(transition_matrix)[1] <- "S"
  return(transition_matrix)
}

Learner <- setRefClass("Learner", 
                       field = list(sequence = "character", transition_matrix = "matrix", learning_factor = "numeric", min_distance_trial = "numeric"),
                       method = list(initialize = function(..., sequence = c(), transition_matrix = array(data = NA, dim = 2),learning_factor = 1.2, min_distance_trial = NA) {
                         callSuper(..., sequence = sequence, transition_matrix = transition_matrix, learning_factor = learning_factor, min_distance_trial = min_distance_trial)
                         },
                         runTrial = function(exclude_backtrack = TRUE) { #generates a sequence of visits from transition matrix, ending when all locations have been visited
                           sequence <<- "S" #set starting location
                           while(sum(colnames(transition_matrix) %in% sequence) < ncol(transition_matrix)) { #while not all locations have been visited, simulate more transitions
                             if(exclude_backtrack & length(sequence) >= 2){ previous_location <- sequence[length(sequence)-1]
                                new_transition_matrix <- transition_matrix[,!colnames(transition_matrix) == previous_location]
                                sequence <<- c(sequence, sample(colnames(new_transition_matrix), 1, prob = new_transition_matrix[sequence[length(sequence)],]))
                             } else sequence <<- c(sequence, sample(colnames(transition_matrix), 1, prob = transition_matrix[sequence[length(sequence)],])) #get new location by sampling row of transition matrix determined by last location visited, add new location
                             if(length(sequence) > 50) return()
                           }
                         },
                         updateTransitions = function(sequence_path_length) {
                           min_distance_trial <<- sequence_path_length
                           transitions <- unique(paste0(head(sequence, -1), sequence[-1])) #gets transitions as a vector of ordered pair characters
                           for (transition in transitions) {
                             transition <- strsplit(transition, NULL)[[1]]
                             transition_matrix[transition[1], transition[2]] <<- transition_matrix[transition[1],transition[2]] * learning_factor #increase likelihood of transitions made
                             transition_matrix <<- t(apply(transition_matrix, MARGIN = 1, function(X) X/(sum(X)))) #normalize
                           }
                         })
                       )

#testing
# learner1 <- Learner(sequence = "S", transition_matrix = getDistanceWeightedTransitions(resource_arrays, getStarts(resource_arrays[1], "Vervet", data_full)), learning_factor = 1.2)
# learner1$runTrial()
# learner1
# sequenceDistance <- getSequenceDistance(learner1$sequence, Double_Trapezoid)
# learner1$updateTransitions(sequenceDistance)
# learner1
# 
# testData <- iterativeLearningSimulation(resource_arrays[2], 2, 100, getStarts(resource_arrays[2], "Vervet", data_full))

