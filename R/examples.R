# example changes!
#' Print a string!
#'
#' @param String character containing a string to
#' print
#' @return A printed string
#' @examples
#' leeR_demo()
leeR_demo <- function(string) {
  print(string)
}

library(dismo) # dependencies: raster, sp
library(plyr)

#' get_species_data
#'
#' @param species_name a length-2 vector containing genus and species name
#' @return a dataframe of presence data for the specified species
#' @examples
#' get_species_data(c("aedes", "aegypti"))
get_species_data <- function(species_name) {
  a = 1
  # name for the dataframe/file containing the presence data for the specific species
  presence_name <- paste(c(species_name, "presence"), collapse="_")
  species_presence_filename <- paste(c(species_name, "presence.rds"), collapse="_")
  
  if (exists(presence_name)) {
    # presence data is already loaded in current session; return that data
    return (get(presence_name))
  } else if (file.exists(species_presence_filename)) {
    # presence data is saved locally; load data and return
    assign(presence_name, readRDS(species_presence_filename), envir=.GlobalEnv)
    return (get(presence_name))
  } else {
    # download presence data from GBIF and save locally
    species_presence <- gbif(genus=species_name[1], species=species_name[2])
    assign(presence_name, species_presence, envir=.GlobalEnv)
    # save data under specific species name
    saveRDS(get(presence_name), species_presence_filename)
    return (get(presence_name))
  }
}


#' get_environment_layers
#'
#' @param resolution the spatial resolution of environmental layers: 10, 5, 2.5
#' @return a list of environment raster stack corresponding to environment conditions in each month
#' @examples
#' get_environment_layers()
get_environment_layers <- function(resolution=10) {
  if (environment_loaded) return
  else {
    tmean <- getData('worldclim', var='tmean', res=resolution)
    tmin <- getData('worldclim', var='tmin', res=resolution)
    tmax <- getData('worldclim', var='tmax', res=resolution)
    prec <- getData('worldclim', var='prec', res=resolution)
    climate_jan <- stack(tmean$tmean1, tmin$tmin1, tmax$tmax1, prec$prec1)
    climate_feb <- stack(tmean$tmean2, tmin$tmin2, tmax$tmax2, prec$prec2)
    climate_mar <- stack(tmean$tmean3, tmin$tmin3, tmax$tmax3, prec$prec3)
    climate_apr <- stack(tmean$tmean4, tmin$tmin4, tmax$tmax4, prec$prec4)
    climate_may <- stack(tmean$tmean5, tmin$tmin5, tmax$tmax5, prec$prec5)
    climate_jun <- stack(tmean$tmean6, tmin$tmin6, tmax$tmax6, prec$prec6)
    climate_jul <- stack(tmean$tmean7, tmin$tmin7, tmax$tmax7, prec$prec7)
    climate_aug <- stack(tmean$tmean8, tmin$tmin8, tmax$tmax8, prec$prec8)
    climate_sep <- stack(tmean$tmean9, tmin$tmin9, tmax$tmax9, prec$prec9)
    climate_oct <- stack(tmean$tmean10, tmin$tmin10, tmax$tmax10, prec$prec10)
    climate_nov <- stack(tmean$tmean11, tmin$tmin11, tmax$tmax11, prec$prec11)
    climate_dec <- stack(tmean$tmean12, tmin$tmin12, tmax$tmax12, prec$prec12)
    
    # save 12-month climates globally
    climates <<- c(climate_jan, climate_feb, climate_mar, climate_apr, 
                   climate_may, climate_jun, climate_jul, climate_aug, 
                   climate_sep, climate_oct, climate_nov, climate_dec)
    environment_loaded <<- TRUE
  }
}


#' get_boundary
#'
#' @param region region name stored in a vector (vector length depends on the granularity level of the region)
#' @return a spatial polygon outlining the region
#' @examples
#' get_boundary(c("USA")) 
#' get_boundary(c("USA","Pennsylvania")) 
#' get_boundary(c("USA","Pennsylvania", "Allegheny"))
get_boundary <- function(region) {
  if (length(region) == 1) {
    target_region <- getData("GADM", country=region, level=0)
  }
  else if (length(region) == 2) {
    nation <- getData("GADM", country=region[1], level=1)
    state <- region[2]
    target_region <- nation[match(toupper(state),toupper(nation$NAME_1)),]
  }
  else if (length(region) == 3) {
    nation <- getData("GADM", country=region[1], level=2)
    state <- region[3]
    target_region <- nation[match(toupper(state),toupper(nation$NAME_2)),]
  }
  return (target_region)
}

#' get_maxent_predict
#'
#' @param month an integer representing the month of the year (from 1 to 12)
#' @param species_presence a dataframe containing longitudes and latitudes of species observation
#' @param climate a raster stack of environmental layers in the corresponding month
#' @param species_name a length-2 vector containing genus and species name
#' @return raw: raw probability raster
#' @return logis: logistic probability raster
#' @examples
#' get_maxent_predict(6, mosquito_presence, climate_jun, c("aedes", "aegypti"))
get_maxent_predict <- function(month, species_presence, climate, species_name) {
  # java setup for maxent modeling
  jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
  # set prediction output filename ### need to include species (and resolution?)
  predict_filename_logis <- paste(paste(species_name, collapse="_"), "prediction_logis", month, ".asc", sep="_")
  predict_filename_raw <- paste(paste(species_name, collapse="_"), "prediction_raw", month, ".asc", sep="_")
  
  # prepare maxent_predict_raw
  if (file.exists(predict_filename_raw) && file.exists(predict_filename_logis)) {
    maxent_predict_raw <- raster(predict_filename_raw)
    maxent_predict_logis <- raster(predict_filename_logis)
  } else {
    model_month <- maxent(climate, species_presence)
    maxent_predict_raw <- predict(model_month, climate, filename=predict_filename_raw, 
                                  overwrite=T, args=c("outputformat=raw"))
    maxent_predict_logis <- predict(model_month, climate, filename=predict_filename_logis, overwrite=T)
  }
  
  projection(maxent_predict_raw) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  projection(maxent_predict_logis) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  return (list(raw=maxent_predict_raw, logis=maxent_predict_logis))
}

#' population_sampling
#'
#' @param maxent_predict_raw character containing a string to
#' print
#' @param maxent_predict_logis
#' @param region
#' @param species
#' @param N
#' @param month_logis_base
#' @param month
#' @return month_samples: sample points drawn from the given month
#' @return month_region_log: logistic probabilities in the given region in given month
#' @return month_N: the number of samples drawn from the given month proportional to the sample size from the reference month
#' @examples
#' leeR_demo()
# helper: population sampling given distribution
population_sampling <- function(maxent_predict_raw, maxent_predict_logis,
                                region=c("ITA"),
                                species, N, month_logis_base, month) {
  # get region's shape (sptialpolygons) and extract corresponding grids from the prediction raster
  target_region <- get_boundary(region)
  region_pred_raw <- mask(crop(maxent_predict_raw, target_region), target_region)
  
  region_pred_logis <- mask(crop(maxent_predict_logis, target_region), target_region)
  region_pred_logis_base <- mask(crop(month_logis_base, target_region), target_region)
  N <- mean(getValues(region_pred_logis), na.rm=T) / mean(getValues(region_pred_logis_base), na.rm=T) * N
  
  # use raw probabilities from maxent model
  #plot(region_pred_raw, main=paste("Raw Probabilities -", region))
  probrast<-region_pred_raw
  
  # normalize the region probability raster by dividing 
  # by the sum of all probabilities within the region:
  probrast<-probrast/sum(getValues(probrast), na.rm=T)
  
  # a function to sample with replacement N points on a raster, with 
  # inclusion probabilities defined by the raster values
  probsel<-function(probrast, N, x_radius, y_radius){
    # extract values (probabilities) from raster
    x <- getValues(probrast)
    # set NA cells in raster to zero
    x[is.na(x)] <- 0
    # sample with replacement N cells
    samp <- sample(nrow(probrast)*ncol(probrast), replace=T, size=N, prob=x)
    # count the number of times each cell is selected
    freq_table <- count(samp)
    samprast <- probrast
    # set value of sampled cell to its count (selected times)
    samprast[freq_table$x] <- freq_table$freq 
    
    # a function to extract centers of selected grid cells from raster
    get_points <- function(i) (rasterToPoints(samprast, fun=function(x){x>=i}))
    
    # extract and combine cell centers 
    # according to the number of times each gets selected
    points <- do.call(rbind, lapply(1:max(freq_table$freq), get_points))
    
    # add within-cell variation to the coordinates selected
    x_vary <- runif(N, min = -x_radius, max = x_radius)
    y_vary <- runif(N, min = -y_radius, max = y_radius)
    points[,"x"] <- points[,"x"] + x_vary
    points[,"y"] <- points[,"y"] + y_vary
    
    # convert to SpatialPoints
    points<-SpatialPoints(points)
    return(points)
  }
  
  x_radius <- res(region_pred_raw)[1]/2
  y_radius <- res(region_pred_raw)[2]/2
  
  # Time: sampling 1000 
  start_time <- Sys.time()
  sample_points_1000<-probsel(probrast, N, x_radius, y_radius)
  end_time <- Sys.time()
  return (list(month_samples=sample_points_1000, month_region_log=region_pred_logis, month_N=N))
}


#' run_simulation
#'
#' @param speciess_name
#' @param month
#' @param N
#' @param region
#' @param resolution
#' @return log
#' @return sample
#' @return range
#' @return N
#' @examples
#' leeR_demo()
run_simulation <- function(species_name, month, N, region, resolution=10) {
  species_presence <- get_species_data(species_name)
  species_month <- species_presence[species_presence$month==month,c("lon", "lat")]
  species_month <- subset(species_month, !is.na(lon) & !is.na(lat))
  get_environment_layers(resolution)
  
  # get environment data by month
  climate_month <- climates[[month]]
  
  # get maxent prediction
  month_predicted_distributions <- get_maxent_predict(month, species_month, climate_month, species_name)
  month_logis_base <- month_predicted_distributions$logis
  
  monthly_logistic_distributions <- stack()
  monthly_samples <- c()
  monthly_N <- c()
  probability_ranges <- list()
  
  for (i in c(1:length(months))) {
    month <- months[i]
    # get species data by month
    species_month <- species_presence[species_presence$month==month,c("lon", "lat")]
    species_month <- subset(species_month, !is.na(lon) & !is.na(lat))
    
    # get environment data by month
    climate_month <- climates[[month]]
    
    # get maxent prediction
    month_predicted_distributions <- get_maxent_predict(month, species_month, climate_month, species_name)
    month_raw <- month_predicted_distributions$raw
    month_logis <- month_predicted_distributions$logis
    
    # sample from predicted distribution
    result <- population_sampling(month_raw, month_logis, region, species_name, N, month_logis_base, month)
    month_samples <- result$month_samples
    month_region_log <- result$month_region_log
    monthly_N <- c(monthly_N, result$month_N)
    monthly_logistic_distributions <- stack(monthly_logistic_distributions, month_region_log)
    monthly_samples <- c(monthly_samples, month_samples)
    probability_ranges[[i]] = range(getValues(month_logis), na.rm=T)
  }
  return  (list(log=monthly_logistic_distributions, sample=monthly_samples, 
                range=probability_ranges, N=monthly_N))
  
}
