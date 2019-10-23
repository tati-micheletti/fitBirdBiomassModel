eBirdDataLoading <- function(pathData,
                             monthsToConsider = 5:8, 
                             yearToConsider = 2005,
                             url, archive, targetFile){
  birdsDT <- reproducible::prepInputs(url = url, 
                                      archive = archive, targetFile = targetFile, 
                                          destinationPath = pathData, fun = "data.table::fread")
  birdsDT <- birdsDT[, c('OBSERVATION COUNT', 'LATITUDE', 'LONGITUDE', 'OBSERVATION DATE')]
  names(birdsDT) <- c("observationCounts", "y", "x", "observationDate")
  birdsDT[, YEAR := substrBoth(strng = observationDate, fromEnd = FALSE, howManyCharacters = 4)]
  birdsDT <- birdsDT[YEAR == yearToConsider, ]
  birdsDT[observationCounts == 'X', observationCounts := 0]
  birdsDT <- birdsDT[, observationCounts := as.numeric(observationCounts)]
  birdsDT[, MONTHS := substrBoth(strng = 
                                       substrBoth(strng = observationDate, 
                                                  fromEnd = FALSE, 
                                                  howManyCharacters = 7), 
                                     fromEnd = TRUE, howManyCharacters = 2)] # Extract Year to match Knn and Months
  birdsDT <- birdsDT[, MONTHS := as.numeric(MONTHS)]
  birdsDT <- birdsDT[MONTHS %in% monthsToConsider,] # Only summer months to avoid winter bias
  birdsDT[, c("observationDate", "YEAR", "MONTHS") := NULL] # remove columns that we don't need
  return(birdsDT)
}