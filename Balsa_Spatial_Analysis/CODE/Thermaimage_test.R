library(Thermimage)

f <- "F:/BCI Ochroma Cam1/Ochroma_19_12_14/FLIR0230.csq" #set file location
x <- frameLocates(f) #doesn't work. Maybe not setup for T 1020 camera?
x<-frameLocates(vidfile = system.file("extdata", "SampleSEQ.seq", package = "Thermimage"))

#Try other functions
convertflirVID(f)
