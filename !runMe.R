library("SpaDES")
usrEmail = if (pemisc::user() %in% c("Tati", "tmichele")) "tati.micheletti@gmail.com" else NULL
googledrive::drive_auth(email = usrEmail)
setPaths(modulePath = dirname(getwd()), cachePath = checkPath(file.path(getwd(), "cache"), create = TRUE),
         outputPath = checkPath(file.path(getwd(), "outputs"), create = TRUE),
         inputPath = checkPath(file.path(getwd(), "data"), create = TRUE))
getPaths() # shows where the 4 relevant paths are

times <- list(start = 1, end = 1)

parameters <- list()

# To determine the tree species to use

modules <- list("fitBirdBiomassModel")
objects <- list(
  "treeSp" = c("Pice_Eng", "Pice_Gla", "Pice_Mar", "Pinu_Con", 
              "Popu_Tre", "Pseu_Men") # REVIEW DEPENDING ON WHAT GREG GIVES ME
)
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)
