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
  #reqdPkgs = list("raster", "data.table", "reproducible", "magrittr", "tati-micheletti/usefun@master"),
  reqdPkgs = list("raster", "data.table", "reproducible", "magrittr", "PredictiveEcology/usefulFuns"),
  parameters = rbind(
    defineParameter("woodType", "logical", FALSE, min, max, "Should the biomass of hard and softwood be summed (TRUE) or 
                    used as individual species (FALSE)?"),
    defineParameter("monthsToConsider", "numeric", 5:8, min, max, "Which months of the year should we consider for the dataset 
                    (i.e. if we have differences in summer and winter)?"),
    defineParameter("yearsToConsider", "numeric", NULL, min, max, 
                    "Which years should we consider for the dataset?"),
    defineParameter("fixKNNageWithFRI", "logical", FALSE, min, max, 
                    "Should ages from FRI be used to correct KNN data?"),
    defineParameter("isBirdDTlonglat", "logical", FALSE, min, max, 
                    paste0("If birdDT IS provided by the user, is it in longlat projection?",
                           "If NOT, then need to be provided in the same projection as rasterTemplate"))
  ),
  inputObjects = bindrows(
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
                 sourceURL = NA),
    expectsInput(objectName = "hardwoodSpecies", objectClass = "character", 
                 desc = "which species of the species list are considered hardwood", 
                 sourceURL = NA),
    expectsInput(objectName = "softwoodSpecies", objectClass = "character", 
                 desc = "which species species list are considered softwood", 
                 sourceURL = NA),
    expectsInput(objectName = "covariateStack", objectClass = "RasterStack", 
                 desc = paste0("If the user wants to provide extra layers, can pass here."), 
                 sourceURL = NA),
    createsOutput(objectName = "birdsDT", objectClass = "data.table", 
                  desc = paste0("Bird dataset. Needs to contain counts for each coordinate.",
                                " Mainly an output if I am using the eBird data. But can be provided."))
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "birdModel", objectClass = "list", 
                  desc = "list with model objects of each one of the species"),
    createsOutput(objectName = "covarTable", objectClass = "list", 
                  desc = "List of covariates for each one of the species"),
    createsOutput(objectName = "birdsDT", objectClass = "data.table", 
                  desc = "Bird dataset. Needs to contain counts for each coordinate"),
    createsOutput(objectName = "frstAttStk", objectClass = "RasterStack", 
                  desc = "Knn (Biomass) Layer + LCC05 + Age"),
    createsOutput(objectName = "dt", objectClass = "data.table", 
                  desc = "Data table of counts and biomass to fit the models"),
    createsOutput(objectName = "birdModel.step", objectClass = "", 
                  desc = "AIC from model fitting")
    
  )
))

doEvent.fitBirdBiomassModel = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (!is.null(sim$covariateStack)){
        covariateStackNames <- names(covariateStack)
        sim$covariateStack <- stack(lapply(sim$covariateStack, function(r){
          ras <- postProces(r, studyArea = sim$studyArea, 
                            rasterToMatch = sim$templateRaster,
                            destinationPath = dataPath(sim),
                            overwrite = TRUE,
                            omitArgs = "overwrite",
                            userTags = c("objectName:r",
                                         "module:fitBirdBiomassModel"))
          return(ras)
        }))
        names(sim$covariateStack) <- covariateStackNames
      }
        
      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim), "fitBirdBiomassModel", "prepData")
      sim <- scheduleEvent(sim, time(sim), "fitBirdBiomassModel", "createModels")
      sim <- scheduleEvent(sim, time(sim), "fitBirdBiomassModel", "createCovarTables")
    },
    prepData = {
      sim$templateRaster[] <- sim$templateRaster[]
      LCC05 <- prepInputs(url = "https://drive.google.com/open?id=1ziUPnFZMamA5Yi6Hhex9aZKerXLpVxvz",
                     targetFile = "LCC2005_V1_4a.tif",
                     studyArea = sim$studyArea, 
                     rasterToMatch = sim$templateRaster,
                     destinationPath = dataPath(sim),
                     method = "ngb", 
                     useGDAL = FALSE,
                     overwrite = TRUE,
                     omitArgs = "overwrite",
                     userTags = c("objectName:LCC05", 
                                  "module:fitBirdBiomassModel"))
      names(LCC05) <- "LCC05"
      
      # 2a. Get NFI kNN data for tree-species-wise proportions of total live above ground biomass (AGB), and crop it to the study area
      treespLayers <- stack(lapply(X = sim$treeSp, 
                                   function(sp){
                                   return(prepInputs(url = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                                             "canada-forests-attributes_attributs-",
                                                             "forests-canada/2011-attributes_attributs-2011/",
                                                             "NFI_MODIS250m_2011_kNN_Species_", 
                                                             sp,
                                                             "_v1.tif"),
                                                destinationPath = dataPath(sim),
                                                useGDAL = FALSE,
                                                filename2 = file.path(dataPath(sim), paste0(sp, "_biomass.tif")),
                                                overwrite = TRUE,
                                                studyArea = sim$studyArea,
                                                rasterToMatch = sim$templateRaster))}))
      names(treespLayers) <- sim$treeSp
      
      # 2b. Get NFI kNN data for merchantable volume, and crop it to the study area
      mvolLayers <- stack(prepInputs(url = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                             "canada-forests-attributes_attributs-",
                                             "forests-canada/2011-attributes_attributs-2011/",
                                             "NFI_MODIS250m_2011_kNN_Structure_Volume_Merch_v1.tif"),
                                destinationPath = dataPath(sim),
                                useGDAL = FALSE,
                                studyArea = sim$studyArea,
                                overwrite = TRUE,
                                rasterToMatch = sim$templateRaster))
      names(mvolLayers) <- c("mvol")
      # 2c. Get both KNN and FRI data for Tree Age, and crop it to the study area. Use FRI data 
      # for where it is available (i.e. 32% of the bird point counts)
      ageLayer <- prepInputs(url = paste0("https://drive.google.com/open?id=19lJM",
                                             "ILt3jee14dkrQZLq1UO_eyRaieUE"),
                                   targetFile = "treeAgeRaster250m_BC.tif",
                                   archive = "treeAgeRaster250m_BC.zip",
                                   destinationPath = dataPath(sim),
                             useGDAL = FALSE,
                                   studyArea = sim$studyArea,
                             overwrite = TRUE,
                                   rasterToMatch = sim$templateRaster)
      names(ageLayer) <- c("ageLayer")
      
      ageLayer_KNN <- prepInputs(url = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                                        "canada-forests-attributes_attributs-forests-canada/2011-",
                                                        "attributes_attributs-2011/NFI_MODIS250m_2011_kNN_Structure",
                                                        "_Stand_Age_v1.tif"),
                                     destinationPath = dataPath(sim),
                                 useGDAL = FALSE,
                                     studyArea = sim$studyArea,
                                     overwrite = TRUE,
                                     rasterToMatch = sim$templateRaster)
      names(ageLayer_KNN) <- c("ageLayer_KNN")

      # 3. Get the bird data from eBird's file
      if (is.null(sim$birdsDT)){
        sim$birdsDT <- Cache(eBirdDataLoading, pathData = dataPath(sim),
                             url = sim$dataURL,
                             targetFile = sim$dataFile,
                             archive = sim$dataArchive,
                             monthsToConsider = P(sim)$monthsToConsider, 
                             yearToConsider = P(sim)$yearsToConsider, 
                             userTags = "eBird")
      }
      if (P(sim)$isBirdDTlonglat){
        # eBird is projected as latlong; we need to fix it to the same projection
        # of the otehr objects in the simulation: we use sim$templateRaster's projection
        projxy <- spTransform(x = SpatialPoints(cbind(sim$birdsDT$x, sim$birdsDT$y), 
                                                proj4string = CRS("+proj=longlat")), 
                              CRSobj = raster::crs(LCC05))
        projxy <- as.data.table(coordinates(projxy))
        names(projxy) <- c("x", "y")
        sim$birdsDT[, c("x", "y") := projxy]
      }

      # 4. Extract from LCC05 and knn the species' biomasses per species and the LCC cover type
      treespLayers <- 0.01 * treespLayers * mvolLayers # convert to merch vol density (m3/ha)
      names(treespLayers) <- sim$treeSp
      browser()
      ageLayer <- postProcess(ageLayer, rasterToMatch = LCC05)
      sim$frstAttStk <- raster::stack(treespLayers, LCC05, ageLayer, ageLayer_KNN, 
                                      sim$covariateStack)
      valuesStack <- data.table(raster::extract(sim$frstAttStk, 
                                                data.table(x = sim$birdsDT$x, 
                                                           y = sim$birdsDT$y)))
      # 6. Make up the table with: observationCounts, biomass softwood, biomass hardwood
      sim$dt <- cbind(valuesStack, sim$birdsDT)
      # Filter points to boreal forest
      sim$dt <- na.omit(sim$dt, col = "ageLayer_KNN")
      # 5. Sum the biomass of each type of tree: softwood or hardwood
      if (P(sim)$woodType){
        sim$dt[, softwood := sum(.SD), by = 1:NROW(sim$dt), 
                    .SDcols = sim$softwoodSpecies]
        sim$dt[, hardwood := sum(.SD), by = 1:NROW(sim$dt), 
                    .SDcols = sim$hardwoodSpecies]
      }
      # Correct age in relation to KNN. 
      sim$dt[, ageLayer_KNN := ageLayer_KNN-(2011-as.numeric(YEAR))] # Correcting the KNN data to match
                                                                     # bird data
      # Negative numbers mean we don't know the real forest age at the time of observation, then, 
      # need to exclude
      sim$dt <- sim$dt[ageLayer_KNN >= 0, ]
      
      if (P(sim)$fixKNNageWithFRI){
        sim$dt[is.na(ageLayer), ageLayer := ageLayer_KNN] # Using FRI age where for we have it 
        sim$dt[, ageLayer_KNN := NULL]
        names(sim$dt)[names(sim$dt) == "ageLayer"] <- "ageLayer_KNN"
      }
      names(sim$dt)[names(sim$dt) == "ageLayer_KNN"] <- "age"
    },
    createModels = {
      if (is.null(sim$formula)){
        if (P(sim)$woodType){
          formula <- "observationCounts ~ softwood + hardwood + age"
        } else {
          #cleanup
          toKeep <- c(sim$treeSp, "observationCounts", "age")
          sim$dt <- sim$dt[, ..toKeep]
          # Remove columns with only 0
          sim$dt
          sim$dt <- sim$dt[, colSums(sim$dt != 0) > 0, with = FALSE]
          formula <- paste0("observationCounts ~ ",
                            paste(names(sim$dt)[names(sim$dt) != "observationCounts"], 
                                  collapse = " + "))
        }
      } else {
        formula <- sim$formula
      }
      # 7. Make the glm models:  observationCounts ~ softwood + hardwood biomass or by species or 
      # own passed model
      sim$birdModel <- glm(data = sim$dt,
                           formula = formula,
                           family = "poisson")
      print(summary(sim$birdModel))

      if (length(sim$birdModel$coefficients)-1 > 1)
        sim$birdModel.step <- MASS::stepAIC(sim$birdModel, trace = TRUE)
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
    sim$studyArea <- Cache(usefulFuns::provinceBCRStudyArea, 
                           province = "British Columbia", 
                           country = "Canada", userTags = c("objectName:studyArea", 
                                                            "module:inputObjects"))
  }
  if(!suppliedElsewhere("templateRaster", sim)){
    sim$templateRaster <- prepInputs(url = paste0("https://drive.google.com/open?id=19lJM",
                                                  "ILt3jee14dkrQZLq1UO_eyRaieUE"),
                                     targetFile = "treeAgeRaster250m_BC.tif",
                                     archive = "treeAgeRaster250m_BC.zip",
                                     destinationPath = Paths$inputPath,
                                         userTags = c("objectName:templateRaster",
                                                      "module:inputObjects")
    )
  }
  if(!suppliedElsewhere("birdsDT", sim)){
    params(sim)$fitBirdBiomassModel$isBirdDTlonglat <- TRUE
    warning("birdsDT was not supplied, so eBird data will be used. 
Therefore, fitBirdBiomassModel will be turned TRUE", immediate = TRUE)
  }
  if(!suppliedElsewhere("treeSp", sim)){
    # from the RIA ws3 model. Choose the species you want by commenting out (#)
    # the ones you don't want. You can also specify other species in other areas, but 
    # don't forget to specify which species are hardwood and which are 
    # softwood if you change the species to another area that not BC.
    sim$treeSp <- c(
                    "Pice_Eng", # engelmann_spruce
                    "Pice_Gla", # white_spruce
                    "Pice_Mar", # black_spruce
                    "Pinu_Con", # lodgepole_pine
                    "Popu_Tre", # aspen
                    "Pseu_Men", # douglas_fir
                    "Abie_Las", # subalpine_fir
                    "Popu_Spp", # poplar
                    "Pinu_Mon", # western_white_pine
                    "Sali_Spp", # willow
                    "Tsug_Mer", # mountain_hemlock
                    "Lari_Lya", # alpine_larch
                    "Pinu_Alb", # whitebark_pine
                    "Pice_Spp", # spruce
                    "Abie_Ama", # amabilis_fir
                    "Popu_Tri", # cottonwood (black?)
                    "Lari_Lar", # tamarack
                    "Tsug_Het", # western_hemlock
                    "Abie_Bal", # balsam_fir
                    "Thuj_Pli", # redcedar
                    "Betu_Pap"  # paper_birch
    )
  }

  if (!suppliedElsewhere("softwoodSpecies")){
    sim$softwoodSpecies <- c(
      "Pice_Eng", # engelmann_spruce
      "Pice_Gla", # white_spruce
      "Pice_Mar", # black_spruce
      "Pinu_Con", # lodgepole_pine
      "Pseu_Men", # douglas_fir
      "Abie_Las", # subalpine_fir
      "Pinu_Mon", # western_white_pine
      "Tsug_Mer", # mountain_hemlock
      "Lari_Lya", # alpine_larch
      "Pinu_Alb", # whitebark_pine
      "Pice_Spp", # spruce
      "Abie_Ama", # amabilis_fir
      "Lari_Lar", # tamarack
      "Tsug_Het", # western_hemlock
      "Abie_Bal", # balsam_fir
      "Thuj_Pli" # redcedar
    )
    if (!all(sim$softwoodSpecies %in% sim$treeSp)) stop("Softwood species default is not matching treeSp.
                                                  Please provide 'softwoodSpecies' to match 'treeSp'")
  }
  if (!suppliedElsewhere("hardwoodSpecies")){
    sim$hardwoodSpecies <- c(
      "Popu_Tre", # aspen
      "Popu_Spp", # poplar
      "Sali_Spp", # willow
      "Popu_Tri", # cottonwood (black?)
      "Betu_Pap"  # paper_birch
    )
    if (!all(sim$hardwoodSpecies %in% sim$treeSp)) stop("Hardwood species default is not matching treeSp.
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
