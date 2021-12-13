---
title: "fitBirdBiomassModel"
author: "Tati Micheletti"
date: "13 December 2021"
output:
  html_document:
    df_print: paged
---
# Overview

This module creates a model using total biomass of soft and hardwood for Red-breasted nuthatch for the RIA region. 
The covariates can be used in an anticipation function for harvesting.

# Usage

```{r module_usage}

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
objects <- list()
inputs <- list()
outputs <- list()

mySim <- simInit(times = times, params = parameters, modules = modules,
                 objects = objects)

mySimOut <- spades(mySim)
```


# Events

prepData
    1. Get knn data for biomass, crop it to the study area  
    2. Get LCC05 data -- we assume all classifications match 2001  
    3. Get the bird data from eBird's file  
    4. Extract from LCC05 and knn the species' biomasses per species and the LCC cover type  
    5. Sum the biomass of each type of tree: softwood or hardwood  
    6. Make up the table with: observationCounts, biomass softwood, biomass hardwood  
createModels  
    7. Make the glm models:  observationCounts ~ softwood + hardwood biomass  
createCovarTables  
    8. Get the covariates and put it on a table  

# Data dependencies

## Input data

By default, all the data needed is automatically downloaded from GDrive. The data comes from eBirds. The default model fits RBNU data in BC for the year 2005. If you want to provide your own data, you can provide the following objects: `dataURL`, `dataFile` and `dataArchive`. See module's metadata for further information. The syntax is:

```{r url, echo = TRUE, eval = FALSE}
objects <- list(
               "dataURL" = "https://drive.google.com/file/d/1hUvj5PHNDWe1VWZReciBUt4wQHAF3Gw5/view?usp=sharing",
               "dataFile" = "ebd_CA-BC_rebnut_relAug-2019.txt",
               "dataArchive" = "ebd_CA-BC_rebnut_relAug-2019.zip"
)
```


You can also provide which tree species to use with similar syntax:

```{r treeSp, echo = TRUE, eval = FALSE}
objects <- list(
  "treeSp" =  c("Pice_Eng", "Pice_Gla", "Pice_Mar", "Pinu_Con", "Popu_Tre", "Pseu_Men")
)
```


The model can handle both each individual species' biomass, or concatenate soft and hardwood ones. As default, it will use all species provided as covariates in the model. To concatenate and model exclusively soft and hardwood, set the parameter `woodType = TRUE` as:

```{r wood, echo = TRUE, eval = FALSE}
  parameters <- list(
  "fitBirdBiomassModel" = list(
    "woodType" = TRUE 
   )
  )
```


And pass the tree species that compose soft and the ones that compose hardwood. There is a default for all `treeSp` (as described above), `softwood` and `hardwood` as:

```{r softhard, echo = TRUE, eval = FALSE}
  softwoodSpecies <- c("Pice_Eng", "Pice_Gla", "Pice_Mar", "Pinu_Con", "Pseu_Men")
  hardwoodSpecies <- c("Popu_Tre")
```

  However, if you provide a mismatching `treeSp` and `softwoodSpecies` or `hardwoodSpecies` list, it will return an error.
  
  You can also pass another 'studyArea' but if you want to use another species, it is necessary to pass the url, targetFile and archive where the data is in GoogleDrive or another address. Also, the data NEEDS to have the collowing fields, identically written (eBird compatible): 'OBSERVATION COUNT', 'LATITUDE', 'LONGITUDE', 'OBSERVATION DATE'. While LATITUDE and LONGITUDE need to be in latlong format, OBSERVATION DATE needs to be in the format: `YYYY-MM-DD`.

 At last, you can provide different months of data collection to compose the model (i.e. differences among seasons). Default is to summer months `5:8`. 

  
## Output data

The main output is `covarTable`, which contains the covariates of the `glm` model. It can be accessed:

```{r covarTable, echo = TRUE, eval = FALSE}

> summary(sim$birdModel)

Call:
glm(formula = formula, family = "poisson", data = sim$dt)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-4.1769  -0.9952  -0.7171   0.0821  20.1344  

Coefficients:
              Estimate Std. Error z value Pr(>|z|)    
(Intercept)  7.594e-01  5.684e-03  133.60   <2e-16 ***
softwood     2.181e-03  2.114e-05  103.19   <2e-16 ***
hardwood    -1.176e-02  2.697e-04  -43.62   <2e-16 ***
age         -8.235e-04  5.925e-05  -13.90   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 184252  on 63135  degrees of freedom
Residual deviance: 169747  on 63132  degrees of freedom
AIC: 319047

Number of Fisher Scoring iterations: 6

> sim$covarTable
  (Intercept)      softwood      hardwood           age 
 0.7594492894  0.0021811110 -0.0117632087 -0.0008234927
```
