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
                    (i.e. if we have differences in summer and winter)?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "studyArea", objectClass = "shapefile", 
                 desc = "Should be BC, otherwise need to pass the correct data", 
                 sourceURL = ""),
    expectsInput(objectName = "treeSp", objectClass = "character", 
                 desc = "Species to consider in the LandR format: 'Tree_spe'", 
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

      sim$studyArea <- spTransform(sim$studyArea, CRS('+init=epsg:3005'))
        
      LCC05 <- Cache(reproducible::prepInputs, 
                     url = "https://drive.google.com/open?id=1ziUPnFZMamA5Yi6Hhex9aZKerXLpVxvz",
                     targetFile = "LCC2005_V1_4a.tif",
                     studyArea = sim$studyArea, 
                     #rasterToMatch = sim$templateRaster,
                     useSAcrs = TRUE,
                     destinationPath = dataPath(sim), 
                     filename2 = "LCC05",
                     overwrite = TRUE,
                     method = "ngb", 
                     omitArgs = "overwrite")
      names(LCC05) <- "LCC05"
      #browser()
      #LCC05 <- projectRaster(LCC05, crs=CRS("+init=epsg:3005"))
      
      
      # 2a. Get NFI kNN data for tree-species-wise proportions of total live above ground biomass (AGB), and crop it to the study area
      treespLayers <- stack(lapply(X = sim$treeSp, 
                                   function(sp){
                                   return(Cache(reproducible::prepInputs, 
                                                url = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                                             "canada-forests-attributes_attributs-",
                                                             "forests-canada/2011-attributes_attributs-2011/",
                                                             "NFI_MODIS250m_2011_kNN_Species_", 
                                                             sp,
                                                             "_v1.tif"),
                                                destinationPath = dataPath(sim),
                                                studyArea = sim$studyArea, 
                                                filename2 = paste0("KNN250m_", sp),
                                                rasterToMatch = LCC05,
                                                #rasterToMatch = sim$templateRaster,
                                                overwrite = TRUE, 
                                                omitArgs = "overwrite"))}))
      names(treespLayers) <- sim$treeSp
      
      # 2b. Get NFI kNN data for merchantable volume, and crop it to the study area
      mvolLayers <- stack(Cache(reproducible::prepInputs, 
                                url = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                             "canada-forests-attributes_attributs-",
                                             "forests-canada/2011-attributes_attributs-2011/",
                                             "NFI_MODIS250m_2011_kNN_Structure_Volume_Merch_v1.tif"),
                                destinationPath = dataPath(sim),
                                studyArea = sim$studyArea, 
                                filename2 = "KNN250m_mvol",
                                #rasterToMatch = sim$templateRaster,
                                rasterToMatch = LCC05, 
                                overwrite = TRUE, 
                                omitArgs = "overwrite"))
      names(mvolLayers) <- c("mvol")
       
      # 3. Get the bird data from eBird's file
      sim$birdsDT <- Cache(eBirdDataLoading, pathData = dataPath(sim),
                           url = sim$dataURL,
                           targetFile = sim$dataFile,
                           archive = sim$dataArchive,
                           monthsToConsider = P(sim)$monthsToConsider, 
                           yearToConsider = 2005, userTags = "eBird")
      
      projxy <- spTransform(SpatialPoints(cbind(sim$birdsDT$x, sim$birdsDT$y), 
                                               proj4string=CRS("+proj=longlat")), 
                                 CRS("+init=epsg:3005"))@coords
      
      sim$birdsDT$x <- projxy[, 1]
      sim$birdsDT$y <- projxy[, 2]

      # 4. Extract from LCC05 and knn the species' biomasses per species and the LCC cover type
      treespLayers <- 0.01 * treespLayers * mvolLayers # convert to merch vol density (m3/ha)
      names(treespLayers) <- sim$treeSp
      #browser()
      
      sim$frstAttStk <- raster::stack(treespLayers, LCC05)
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
      #browser()
      # 6. Make up the table with: observationCounts, biomass softwood, biomass hardwood
      sim$dt <- cbind(valuesStack, observationCounts = sim$birdsDT$observationCounts)
      
    },
    createModels = {
      formula <- ifelse(P(sim)$woodType,
                        "observationCounts ~ softwood + hardwood",
                        paste0("observationCounts ~ ",
                               paste(usefulFuns::grepMulti(x = names(sim$dt), patterns = "_"), 
                                     collapse = " + ")))
      #browser()
      # 7. Make the glm models:  observationCounts ~ softwood + hardwood biomass
      sim$birdModel <- glm(data = sim$dt,
                           formula = formula,
                           family = "poisson")
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
    # HACK ##############
    # Not sure where to define sim$studyArea so that is counts as "suppliedElsewhere",
    # so hard-coding it into here to get the thing running for the RIA landbase.
    sim$studyArea <- usefulFuns::provinceBCRStudyArea(province = "British Columbia", country = "Canada")
    #####################
  }
  
  if(!suppliedElsewhere("sppEquivCol", sim)){
    sim$sppEquivCol <- "PSP"
  }
  
  if(!suppliedElsewhere("treeSp", sim)){
    #sim$treeSp <- c("Pice_Eng", "Pice_Gla", "Pice_Mar", "Pinu_Con", "Popu_Tre", "Pseu_Men")
    # from the RIA ws3 model:
    # x  'douglas_fir',
    # x  'white_spruce', 
    # x  'subalpine_fir',
    # x  'poplar',
    # x  'western_white_pine',
    # x  'willow',
    # x  'mountain_hemlock',
    # x  'black_spruce',
    # x  'alpine_larch',
    # x  'whitebark_pine',
    # x  'spruce',
    # x  'amabilis_fir',
    # x  'cottonwood',
    # x  'tamarack',
    # x  'engelmann_spruce',
    # x  'aspen',
    # x  'western_hemlock',
    # x  'balsam_fir',
    # x  'lodgepole_pine',
    # x 'redcedar',
    # x  'paper_birch'
    sim$treeSp <- c(
                    #"Pice_Eng", # engelmann_spruce
                    "Pice_Gla", # while_spruce
                    "Pice_Mar", # black_spruce
                    "Pinu_Con", # lodgepole_pine
                    "Popu_Tre", # aspen
                    #"Pseu_Men", # douglas_fir
                    #"Abie_Las", # subalpine_fir
                    #"Popu_Spp", # poplar
                    #"Pinu_Mon", # western_white_pine
                    #"Sali_Spp", # willow
                    #"Tsug_Mer", # mountain_hemlock
                    #"Lari_Lya", # alpine_larch
                    #"Pinu_Alb", # whitebark_pine
                    "Pice_Spp", # spruce
                    #"Abie_Ama", # amabilis_fir
                    #"Popu_Tri", # cottonwood (black?)
                    #"Lari_Lar", # tamarack
                    #"Tsug_Het", # western_hemlock
                    #"Abie_Bal", # balsam_fir
                    #"Thuj_Pli", # redcedar
                    "Betu_Pap"  # paper_birch
    )
    # sim$treeSp <- c("Pice_Mar", # black_spruce
    #                 "Popu_Spp", # poplar
    #                 "Lari_Lya", # alpine_larch
    #                 "Lari_Lar", # tamarack
    #                 "Pinu_Alb") # whitebark_pine)
  }
  
  if(!suppliedElsewhere("sppEquiv", sim)){
    # Equivalency table for tree species
    data("sppEquivalencies_CA", package = "LandR")
    # Make BC spp equivalencies
    sppEquivalencies_CA <- sppEquivalencies_CA[KNN %in% sim$treeSp,]
    sim$sppEquiv <- sppEquivalencies_CA
  }
  
  if (!suppliedElsewhere("softwoodSpecies")){
    sim$softwoodSpecies <- c("Pice_Mar", "Lari_Lya", "Lari_Lar", "Pinu_Alb")
    #sim$softwoodSpecies <- c("Pice_Eng", "Pice_Gla", "Pice_Mar",
    #                         "Pinu_Con", "Pseu_Men")
    #if (!sim$softwoodSpecies %in% sim$treeSp) stop("Softwood species default is not matching treeSp. 
    #                                               Please provide 'softwoodSpecies' to match 'treeSp'")
  }
  if (!suppliedElsewhere("hardwoodSpecies")){
    sim$hardwoodSpecies <- c("Popu_Spp")
    #sim$hardwoodSpecies <- c("Popu_Tre")
    #if (!sim$hardwoodSpecies %in% sim$treeSp) stop("Hardwood species default is not matching treeSp. 
    #                                               Please provide 'softwoodSpecies' to match 'treeSp'")
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
