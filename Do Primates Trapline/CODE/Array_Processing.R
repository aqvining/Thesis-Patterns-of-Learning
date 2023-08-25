require(sf)
require(plyr)

buildArray <- function(points) {
  #purpose: makes a spatial features collection (sfc) of x,y coordinates from either a list or a matrix of points
  #inputs: sets of coordinate pairs, either in a list or matrix.
  #output: an sfc with spatial points in the geometry column. Point names will be taken from names of list or rownames of matrix
  if (class(points) == "matrix") {
    points_list <- alply(.data = points, .margins = 1, .fun = st_point) #create a list with each element equal to the application of st_point to 1 row of the points matrix
    points_array <- st_sf(geom = st_sfc(points_list), NAME = unname(attr(points_list, "split_labels"))) #turn points_list into a spatial features dataframe
  }
  if (class(points) == "list") {
    points_list <- lapply(points, st_point) #turn each element of points into an st_point object
    if (is.null(names(points_list))) names(points_list) <- as.character(seq(length(points_list))) #name unnamed list
    points_array <- st_sf(geom = points_list, NAME = names(points_list)) #turn list into spatial features dataframe
  }
  return(points_array)
}

getSequenceDistance <- function(sequence, resource_array) {
  #purpose: get total distance between all transitions within an array of resources
  #input: sequence: a vector of characters or numerics representing names of locations in the array OR a single element that can be split into such
  #       resource_Array: a spatial features dataframe with spatial points in the geometry column and the feature NAME
  #output: numeric, total distance between transitions
  distanceMatrix <- getLabelledDistanceMatrix(resource_array)
  distanceTravelled <- 0
  if (length(sequence) == 1) sequence <- strsplit(as.character(sequence), NULL)[[1]]
  if (sequence[1] == "S") sequence <- sequence[-1]
  for (i in 2:length(sequence)) {
    if (sequence[i] == "F") return(distanceTravelled)
    distanceTravelled <- distanceTravelled + distanceMatrix[sequence[i], sequence[i-1]]
  }
  return(distanceTravelled)
}

getLabelledDistanceMatrix <- function(points) {
  #input: a spatial features dataframe with spatial points in the geometries column and feature NAME, a list of length one containing such
  #out: A matrix with distances between each point in the geometry column and axis labelled with NAME
  if("list" %in% class(points)) points <- points[[1]]
  dMatrix <- st_distance(points)
  colnames(dMatrix) <- points$NAME
  rownames(dMatrix) <- points$NAME
  return(dMatrix)
}

#Create Arrays from Traplining Experiments
Double_Trapezoid <- buildArray(matrix(c(1,-7,-3,-11,-11,11,-1,4,8,7,11,-10), byrow = TRUE, nrow = 6))
Pentagon <- buildArray(matrix(c(-5,-6,-8,2,0,8,8,2,5,-6), byrow = TRUE, nrow = 5))
Z <- buildArray(matrix(c(7,-8,-7,-8,-2,-2,2,2,7,8,-7,8), byrow = TRUE, nrow = 6))

#~Testing Script
# DT <- matrix(c(-1,4,8,7,-11,11,1,-7,11,-10,-3,-11), byrow = TRUE, nrow = 6, dimnames = list(as.character(c(4,5,3,1,6,2)))) #named matrix
# buildArray(DT)
# DT <- matrix(c(-1,4,8,7,-11,11,1,-7,11,-10,-3,-11), byrow = TRUE, nrow = 6) #unnamed matrix
# buildArray(DT)
# DT = list(c(-1,4), c(8,7), c(-11,11), c(1,-7), c(11, -10), c(-3,-11)) #unnamed list
# buildArray(DT)
# names(DT) <- as.character(c(4,5,3,1,6,2)) #names list
# buildArray(DT)

