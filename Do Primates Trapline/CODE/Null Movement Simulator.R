library(circular)
library(stats)
library(reshape2)
library(dplyr)
library(sf)
library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~Object Set-up~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Environment <- setRefClass("Environment", 
                           field = list(foragers = "list", patches = "list", bounds = "numeric", initForagers = "numeric", initPatches = "numeric"),
                           method = list(initialize = function(..., #initialization function sets default values for fields
                                                               initForagers = 1,
                                                               initPatches = 10,
                                                               bounds = c(50,50),  #bounds can be passed to the creation of Patches and Foragers in Environment to restrict where they are placed. Also restricts movement of foragers in Environment
                                                               foragers = lapply(rep(1, times = initForagers), function(x) Forager(location = c(runif(x, -bounds[1], bounds[1]), runif(x, -bounds[2],bounds[2])))),
                                                               patches = lapply(1:initPatches, function(x) Patch(location = c(runif(1, -bounds[1], bounds[1]), runif(1, -bounds[2],bounds[2])), name = as.character(x)))){ #default creation of two patches
                             callSuper(..., foragers = foragers, patches = patches, bounds = bounds, initForagers = initForagers, initPatches = initPatches)
                            },
                            progress = function(){
                              for (forager in foragers[sample(length(foragers))]) {#operates on each forager in random order
                                forager$setTarget(patches)
                                forager$move(bounds = bounds)
                              }
                            },
                            plotPatches = function(){
                              Patches <- data.frame(x = sapply(patches, function(p) p$location["x"]), y = sapply(patches, function(p) p$location["y"]), names = sapply(patches, function(p) p$name))
                              pPlot <- ggplot(Patches, aes(x = x, y = y)) +
                                geom_point(shape = 8, color = "green", size = 2) +
                                geom_text(aes(label = names), nudge_y = bounds[2]/20, size = 2) +
                                theme_classic() +
                                coord_fixed()
                              return(pPlot)
                            },
                            plotCurrent = function(){
                              pPlot <- plotPatches()
                              fLocations <- data.frame(x = sapply(foragers, function(f) f$location["x"]), y = sapply(foragers, function(f) f$location["y"]))
                              cPlot <- pPlot + geom_point(data = fLocations, aes(x = x, y = y), shape = 13, size = 3, color = "blue")
                              return(cPlot)
                            },
                            plotPaths = function(){
                              cPlot <- plotCurrent()
                              paths <- data.frame(x = numeric(0), y = numeric(0), forager = numeric(0), timestep = numeric(0))
                              for (i in 1:length(foragers)) {
                                foragerPath <- cbind(foragers[[i]]$path, forager = i, timestep = 1:nrow(foragers[[i]]$path))
                                paths <- rbind(paths, foragerPath)
                              }
                              paths$forager <- factor(paths$forager)
                              timesteps <- max(paths$timestep)
                              pathsPlot <- cPlot + geom_path(data = paths, aes(x = x, y = y, color = forager, alpha = timestep/timesteps))
                              return(pathsPlot)
                            })
                          )

ArrayEnvironment <- setRefClass("ArrayEnvironment", fields= list(sequence = "character", array = "character", trials = "numeric"), contains= "Environment",
                                method = list(
                                  initialize = function(..., sequence = character(0), array = "DT", trials = 0) {
                                    callSuper(..., sequence = sequence, array = array, trials = trials)
                                  },
                                  progress = function(){
                                    callSuper()
                                    for (forager in foragers) {
                                      if (length(unique(forager$visitSeq[-c(1:forager$repeatAvoid)])) == length(patches) | nrow(forager$path) >= 10000) { #end conditions for trial
                                        sequence <<- c(sequence, forager$visitSeq[-c(1:forager$repeatAvoid)])
                                        trials <<- trials + 1
                                        forager$location <- c(runif(1, -bounds[1], bounds[1]), runif(1, -bounds[2],bounds[2]))
                                        names(forager$location) <- c("x", "y")
                                        forager$bearing <- as.numeric(rwrappedcauchy(n = 1, mu = circular(0), rho = 0))
                                        forager$visitSeq <- rep("NA", times = forager$repeatAvoid)
                                      }
                                    }
                                  }
                                ))

Patch  <- setRefClass("Patch", 
                      field = list(location ="XY", name = "character"),
                      method = list(initialize =
                                      function(..., location = st_point(runif(2, -50, 50)), name = sample(LETTERS, 1)){
                                        callSuper(...,location = location, name = name)
                                        names(.self$location) <- c("x","y")
                                      })
)


Forager <- setRefClass("Forager",
                      field = list(location="XY", bearing ="numeric", speed = "numeric", turnDev = "numeric", sight = "numeric", path = "data.frame", visitSeq = "character", targeting = "logical", target = "Patch", repeatAvoid = "numeric"),
                      method = list(initialize = function(..., location = st_point(runif(2, -50, 50)), 
                                                          bearing = as.numeric(rwrappedcauchy(n = 1, mu = circular(0), rho = 0)), 
                                                          speed = 1, 
                                                          turnDev = 0.7, 
                                                          sight = 4,
                                                          path = data.frame(x = location[1], y = location[2]),
                                                          repeatAvoid = 2,
                                                          visitSeq = rep("NA", times = repeatAvoid),#NA's prevent errors when checking for recent visits
                                                          targeting = FALSE,
                                                          target = Patch()) {
                                      if ("numeric" %in% class(location) & length(location) == 2) location <<- st_point(location)
                                      callSuper(..., location = location, bearing = bearing, speed = speed, turnDev = turnDev, sight = sight, path = path, visitSeq = visitSeq, targeting = targeting, target = target, repeatAvoid = repeatAvoid)
                                      names(.self$location) <- c("x","y")
                                    },
                                    setTarget = function(patches) { 
                                      if (! class(patches) == "list") stop("the second argument to function setTarget must be of class 'list'")
                                      if (targeting) return()
                                      choices <- getChoices(.self, patches)
                                      if (! class(choices) == "list") return()
                                      target <<- tSelect(choices)
                                      targeting <<- TRUE
                                    },
                                    tSelect = function(choices) {
                                      if (sum(choices$dist == 0) > 0) return(choices$patches[[which(choices$dist == 0)]])     #if on a patch, return that patch as target
                                      choices$prob <- rep(1/length(choices$patches), times = length(choices$patches)) #random choice
                                      return(selector(choices))
                                    },
                                    move = function(bounds = NA){
                                      if (targeting) { #if a target is given, assess proximity and action
                                        if (getDist(location, target$location) == 0) { 
                                          visitSeq <<- c(visitSeq, target$name) #if on the patch, process (add to visitSeq)
                                          targeting <<- FALSE
                                        } else {
                                          bearing <<- unname(atan2(x = target$location["x"] - location["x"], y=target$location["y"] - location["y"])) #if not on the patch, move directly toward it
                                          if (getDist(location, target$location) <= speed) { location <<- target$location #if target is within reach, set location to patch location
                                          } else location <<- location + speed * c(cos(bearing), sin(bearing)) #otherwise, move closer
                                        }
                                      } else { #if not targeting
                                        bearing <<- as.numeric(circular(bearing) + rwrappedcauchy(n = 1, mu = circular(0), rho = turnDev))
                                        location <<- location + speed * c(cos(bearing), sin(bearing))
                                        if (is.numeric(bounds) & length(bounds == 2)){ #check for valid bounds
                                          turnVarIncrease <- 0
                                          while(abs(location["x"]) > bounds[1] | abs(location["y"]) > bounds[2]) {#check if out of bounds
                                            location <<- unlist(path[nrow(path),]) #reset location
                                            if(turnDev-turnVarIncrease > 0.2) turnVarIncrease <- turnVarIncrease + 0.02 #relax directional persistence
                                            bearing <<- as.numeric(circular(bearing) + rwrappedcauchy(n = 1, mu = circular(0), rho = turnDev - turnVarIncrease)) # get new bearing
                                            location <<- location + speed * c(cos(bearing), sin(bearing)) #try moving again
                                          }
                                        }
                                      }
                                      path[nrow(path) + 1,] <<- location
                                    })
                      )

dForager <- setRefClass("distanceForager", fields = list(), contains = "Forager",
                        methods = list(
                          tSelect = function(choices) {
                            if (sum(choices$dist == 0) > 0) return(choices$patches[[which(choices$dist == 0)]])                                         #if on a patch, return that patch as target
                            choices$prob <- choices$dist^3/sum(choices$dist^3) #distance discounted choice
                            return(selector(choices))
                          }))



selector <- function(choices){
  randomizer <- runif(1)
  i = 1
  while(choices$prob[i] <= randomizer) {
    randomizer <- randomizer - choices$prob[i]
    i <- i + 1
  }
  return(choices$patches[[i]])
}
getChoices <- function(forager, patches) {
  recentVisits <- forager$visitSeq[c((length(forager$visitSeq) - (forager$repeatAvoid - 1)):length(forager$visitSeq))]
  patches <- patches[which(! sapply(patches, function(p) p$name) %in% recentVisits)]
  distances <- sapply(patches, function(patch) getDist(forager$location, patch$location))
  if (sum(distances <= forager$sight) == 0) return(NA)                                                                        #if no patches in sight, return no target
  choices <- list("patches" = patches[which(distances <= forager$sight)], dist = distances[which(distances <= forager$sight)])
  return(choices)
}

getDist = function(c1, c2) {
  #input: two sets of coordinates. x and y coordinates must be names
  #output: Euclidean distance between input points
  c1 = unlist(c1)
  c2 = unlist(c2)
  names(c1) = c("x", "y")
  names(c2) = c("x","y")
  if (! is.numeric(c1) || ! is.numeric(c2)) return(NA)
  if (identical(c1, c2)) return(0)
  unname(sqrt(((c1["x"] - c2["x"]) ** 2) + ((c1["y"] - c2["y"]) ** 2))) 		 #pythagorean theorem
}

###~~~~~~Sample Script~~~~~~~~
testForagers <- vector("list", 10)
for (i in 1:length(testForagers)) testForagers[[i]] <- dForager(turnDev = 0.7)
testPatches <- vector("list", 100)
for (i in 1:length(testPatches)) testPatches[[i]] <- Patch(name = as.character(i))
testEnviron <- Environment(foragers = testForagers, patches = testPatches)
for (t in seq(100)) testEnviron$progress()
testEnviron$plotCurrent()
testEnviron$plotPaths()
