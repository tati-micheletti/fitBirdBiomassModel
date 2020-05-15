defineModule(sim, list(
  name = "fitBirdBiomassModel",
  description = paste0("This module fits a glm using softwood and hardwood biomass",
                       " for bird data from eBirds. It has RBNU as default"),
  keywords = c("forest optimization", "anticipation function"),
  authors = person("Tati", "Micheletti", email = "tati.micheletti@gmail.com", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.6.9002", fitBirdBiomassModel = "1.0.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fitBirdBiomassModel.Rmd"),
  reqdPkgs = list("raster", "data.table", "reproducible", "magrittr", "tati-micheletti/usefun@master"),
  parameters = rbind(
    defineParameter("woodType", "logical", FALSE, min, max, "Should the biomass of hard and softwood be summed (TRUE) or 
                    used as individual species (FALSE)?"),
    defineParameter("monthsToConsider", "numeric", 5:8, min, max, "Which months of the year should we consider for the dataset 
                    (i.e. if we have differences in summer and winter)?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "studyArea", objectClass = "shapefile", 
                 desc = "Should be BC, otherwise need to pass the correct data", 
                 sourceURL = ""),
    expectsInput(objectName = "treeSp", objectClass = "character", 
                 desc = paste0("Species to consider in the LandR format: 'Tree_spe'. The layers for this",
                              "come from the NRCan product largely known as KNN: https://cfs.nrcan.gc.ca/publications?id=35489",
                               "Reference: Mapping attributes of Canada’s forests at moderate resolution through kNN and MODIS ",
                               "imagery. 2014. Beaudoin, A.; Bernier, P.Y.; Guindon, L.; Villemaire, P.; Guo, X.J.; Stinson, G.;",
                               " Bergeron, T.; Magnussen, S.; Hall, R.J. Can. J. For. Res. 44:521-532."), 
                 sourceURL = ""),
    expectsInput(objectName = "dataURL", objectClass = "character",  
                 desc = "url of the bird data in eBird format", 
                 sourceURL = "https://drive.google.com/open?id=1hUvj5PHNDWe1VWZReciBUt4wQHAF3Gw5"),
    expectsInput(objectName = "dataFile", objectClass = "character", 
                 desc = "name of the bird data file with extension: ie. ebd_CA-BC_rebnut_relAug-2019.txt", 
                 sourceURL = NA),
    expectsInput(objectName = "dataArchive", objectClass = "character", 
                 desc = "name of the bird data archive: ie. ebd_CA-BC_rebnut_relAug-2019.zip", 
                 sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "birdModel", objectClass = "list", 
                  desc = "list with model objects of each one of the species"),
    createsOutput(objectName = "covarTable", objectClass = "list", 
                  desc = "List of covariates for each one of the species"),
    createsOutput(objectName = "birdsDT", objectClass = "data.table", 
                  desc = "Bird dataset. Needs to contain counts for each coordinate"),
    createsOutput(objectName = "frstAttStk", objectClass = "RasterStack", 
                  desc = "Knn (Biomass) Layer + LCC05"),
    createsOutput(objectName = "dt", objectClass = "data.table", 
                  desc = "Data table of counts and biomass to fit the models")
    
  )
))

doEvent.fitBirdBiomassModel = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim), "fitBirdBiomassModel", "prepData")
      sim <- scheduleEvent(sim, time(sim), "fitBirdBiomassModel", "createModels")
      sim <- scheduleEvent(sim, time(sim), "fitBirdBiomassModel", "createCovarTables")
      
    },
    prepData = {

      LCC05 <- Cache(prepInputs, url = "https://drive.google.com/open?id=1ziUPnFZMamA5Yi6Hhex9aZKerXLpVxvz",
                                        targetFile = "LCC2005_V1_4a.tif",
                                        studyArea = sim$studyArea, useSAcrs = TRUE,
                                        destinationPath = dataPath(sim), filename2 = "LCC05",
                                        overwrite = TRUE,
                                        method = "ngb", omitArgs = "overwrite")
      names(LCC05) <- "LCC05"
      
      # 2. Get knn data for biomass, crop it to the study area
      biomassLayer <- stack(lapply(X = sim$treeSp, function(sp){
        ras <- Cache(reproducible::prepInputs, url = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                                     "canada-forests-attributes_attributs-",
                                                     "forests-canada/2011-attributes_attributs-2011/",
                                                     "NFI_MODIS250m_2011_kNN_Species_", sp,"_v1.tif"),
                                        destinationPath = dataPath(sim),
                                        studyArea = sim$studyArea, filename2 = paste0("KNN250m_", sp),
                                        rasterToMatch = LCC05, overwrite = TRUE, omitArgs = "overwrite")
        return(ras)
      }))
      names(biomassLayer) <- sim$treeSp
      
      # 3. Get the bird data from eBird's file
      sim$birdsDT <- Cache(eBirdDataLoading, pathData = dataPath(sim),
                           url = sim$dataURL,
                           targetFile = sim$dataFile,
                           archive = sim$dataArchive,
                           monthsToConsider = P(sim)$monthsToConsider, 
                           yearToConsider = 2005, userTags = "eBird")

      # 4. Extract from LCC05 and knn the species' biomasses per species and the LCC cover type
      sim$frstAttStk <- raster::stack(biomassLayer, LCC05)
      valuesStack <- data.table(raster::extract(sim$frstAttStk, 
                                                data.table(x = sim$birdsDT$x, 
                                                           y = sim$birdsDT$y)))

      # 5. Sum the biomass of each type of tree: softwood or hardwood
      if (P(sim)$woodType){
        valuesStack[, softwood := sum(.SD), by = 1:NROW(valuesStack), 
                    .SDcols = sim$softwoodSpecies]
        valuesStack[, hardwood := sum(.SD), by = 1:NROW(valuesStack), 
                    .SDcols = sim$hardwoodSpecies]
      }
      # 6. Make up the table with: observationCounts, biomass softwood, biomass hardwood
      sim$dt <- cbind(valuesStack, observationCounts = sim$birdsDT$observationCounts)
      
    },
    createModels = {

      # 7. Make the glm models:  observationCounts ~ softwood + hardwood biomass
      sim$birdModel <- glm(data = sim$dt,
                           formula = ifelse(P(sim)$woodType,
                                            "observationCounts ~ softwood + hardwood",
                                            paste0("observationCounts ~ ",
                                                         paste(usefun::grepMulti(x = names(sim$dt), patterns = "_"), 
                                                               collapse = " + "))),
                           family = "poisson")
    },
    createCovarTables = {
      # 8. Get the covariates and put it on a table
      sim$covarTable <- sim$birdModel$coefficients
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  if(!suppliedElsewhere("studyArea", sim)){
    sim$studyArea <- usefun::provinceBCRStudyArea(province = "British Columbia", 
                                                  country = "Canada")
  }
  
  if(!suppliedElsewhere("sppEquivCol", sim)){
    sim$sppEquivCol <- "PSP"
  }
  
  if(!suppliedElsewhere("treeSp", sim)){
    sim$treeSp <- c("Pice_Eng", "Pice_Gla", "Pice_Mar", "Pinu_Con", 
                    "Popu_Tre", "Pseu_Men")
    
  }
  
  if(!suppliedElsewhere("sppEquiv", sim)){
    # Equivalency table for tree species
    data("sppEquivalencies_CA", package = "LandR")
    # Make BC spp equivalencies
    sppEquivalencies_CA <- sppEquivalencies_CA[KNN %in% sim$treeSp,]
    sim$sppEquiv <- sppEquivalencies_CA
  }
  
  if (!suppliedElsewhere("softwoodSpecies")){
    sim$softwoodSpecies <- c("Pice_Eng", "Pice_Gla", "Pice_Mar",
                             "Pinu_Con", "Pseu_Men")
    if (!sim$softwoodSpecies %in% sim$treeSp) stop("Softwood species default is not matching treeSp. 
                                                   Please provide 'softwoodSpecies' to match 'treeSp'")
  }
  if (!suppliedElsewhere("hardwoodSpecies")){
    sim$hardwoodSpecies <- c("Popu_Tre")
    if (!sim$hardwoodSpecies %in% sim$treeSp) stop("Hardwood species default is not matching treeSp. 
                                                   Please provide 'softwoodSpecies' to match 'treeSp'")
  }
  if (!suppliedElsewhere("dataURL")){
    sim$dataURL <- extractURL(objectName = "dataURL")
    if (!suppliedElsewhere("dataFile")){
      sim$dataFile <- "ebd_CA-BC_rebnut_relAug-2019.txt"
    } else {
      if (sim$dataFile != "ebd_CA-BC_rebnut_relAug-2019.txt") stop("url was not provided, but dataFile was. Please provide both or none.")
    }
    if (!suppliedElsewhere("dataArchive")){
      sim$dataArchive <- "ebd_CA-BC_rebnut_relAug-2019.zip"
    } else {
      if (sim$dataArchive != "ebd_CA-BC_rebnut_relAug-2019.zip") stop("url was not provided, but dataArchive was. Please provide both or none.")
    }
  }
  
  return(invisible(sim))
}
