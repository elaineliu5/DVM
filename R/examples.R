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

#' @title get_species_data
#' @description The function get_species_data draws species occurrence data from GBIF.
#' get_species_data takes in a species name (a vector containing genus and species name) 
#' and returns a dataframe containing presence data (including longitudes and latitudes)
#' for the specified species drawn from GBIF.
#' @param species_name a vector containing genus and species name
#' @return a dataframe of presence data for the specified species
#' @examples
#' # example: plot mosquito species presence data
#' 
#' mosquito_presence <- get_species_data(c("aedes", "aegypti"))
#' mosquito_location <- mosquito_presence[, c("lon", "lat")]
#' 
#' # install.packages("maptools")
#' library(maptools)
#' data(wrld_simpl)
#' plot(wrld_simpl)
#' points(mosquito_location, col='blue', cex=0.5)
#' title("Aedes Aegypti Observations")
#' 
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


#' @title get_environment_layers
#' @description get_environment_layers draws environmental layers from WorldClim. 
#' It takes in a resolution value and draws monthly environment data from WorldClim.
#' @param resolution the spatial resolution of environmental layers
#' (10, 5, or 2.5 (arc min))
#' @return a list of environment raster stack corresponding to 
#' environment conditions in each month
#' @examples
#' # plot mean temperatures throughout the year
#' 
#' par(mfrow=c(1, 4))
#' environment_data <- get_environment_layers()
#' tmeans <- stack()
#' max_range <- c(0,0)
#' for (i in 1:12) {
#'   tmeans <- stack(tmeans, environment_data[[i]][[1]])
#'   max_range <- range(max_range, range(getValues(environment_data[[i]][[1]]), na.rm=T))
#' }
#' plot(tmeans, breaks=seq(from=min(max_range), to=max(max_range), length.out=100), 
#'      col=rev(terrain.colors(99)), legend=F)
#' plot(tmeans, legend.only=T, 
#'      breaks=seq(from=min(max_range), to=max(max_range), length.out=100),
#'      col=rev(terrain.colors(99)),
#'      axis.args=list(at=seq(max_range[1], max_range[2], length.out=5),
#'                     labels=round(seq(max_range[1], max_range[2], length.out=5),3),
#'                     cex.axis=0.6))
#' 
#' 
#' # plot environment conditions (mean temperature, min temperature, max temperature and precipitation) for June
#' plot(environment_data[[6]])
#' 
get_environment_layers <- function(resolution=10) {
  climates_datafile_name <- "climates.rds"
  if (exists("climates")) return (climates)
  else if (file.exists(climates_datafile_name)) {
    # climate data is saved locally; load data and return
    return (readRDS(climates_datafile_name))
  } else {
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
    climates <- c(climate_jan, climate_feb, climate_mar, climate_apr, 
                  climate_may, climate_jun, climate_jul, climate_aug, 
                  climate_sep, climate_oct, climate_nov, climate_dec)
    saveRDS(climates, climates_datafile_name)
    return (climates)
  }
}

#' @title get_boundary
#' @description get_boundary draws the geographical boundary of a specified region from GADM.
#' It takes in a vector of region name and returns a spatial polygon outlining the region.
#' @param region region name stored in a vector 
#' (vector length depends on the granularity level of the region,
#' from national level down to county level, 
#' i.e. c(“USA”), c(“USA”,“Pennsylvania”), c(“USA
#' ”,“Pennsylvania”, “Allegheny”)), and
#' @return a spatial polygon outlining the region
#' @examples
#' # Plotting the region of Pennsylvania
#' PA <- get_boundary(c("USA", "Pennsylvania"))
#' plot(PA)
#' 
#' # Plotting the environment conditions for Pennsylvania
#' PA_environment_6 <- mask(crop(environment_data[[6]], PA), PA)
#' plot(PA_environment_6)
#' 
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

#' @title get_maxent_predict
#' @description get_maxent_predict uses the Maxent model to predict species distribution in a given month.
#' The model for a species is determined from a set of environmental or climate layers 
#' for a set of grid cells in a landscape, together with a set of sample locations 
#' where the species has been observed. 
#' The computed model is a probability distribution 
#' (indicating the level of suitability for species' living) 
#' over all the grid cells.
#' @param month an integer representing the month of the year (from 1 to 12)
#' @param species_presence a dataframe containing longitudes and latitudes of species observation occurred in the given month
#' @param climate a raster stack of environmental layers in the corresponding month
#' @param species_name a vector containing genus and species name
#' @return raw: raw probability raster
#' @return logis: logistic probability raster
#' @examples
#' # Example: predict mosquito distribution in June
#' month <- 6
#' mosquito_month <- mosquito_presence[mosquito_presence$month==month, c("lon", "lat")]
#' mosquito_month <- subset(mosquito_month, !is.na(lon) & !is.na(lat))
#' mosquito_predictions_6 <- get_maxent_predict(month, mosquito_month, environment_data[[6]], c("aedes", "aegypti"))
#' par(mfrow=c(1,2))
#' plot(mosquito_predictions_6$raw)
#' title("Aedes Aegypti Distribution Raw Predictions")
#' plot(mosquito_predictions_6$logis)
#' title("Aedes Aegypti Distribution Logistic Predictions")
#' 
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

#' @title population_sampling
#' @description Given a specific region and month, 
#' population_sampling function zooms in on the probability distributions 
#' in the region, calculates its average logistic probabilities 
#' and calculate the corresponding sample size for 
#' the given month (based on the ratio of mean logistic probabilities 
#' in the given month and the reference month)
#' @param maxent_predict_raw global raw probabilities raster for the given month
#' @param maxent_predict_logis global logistic probabilities raster for the given month
#' @param region a vector of region names specifying a region of interest 
#' from national levels down to county levels
#' @param species a vector containing genus and species name
#' @param N number of species samples to draw in the reference month
#' @param month_logis_base global logistic probabilities raster for the reference month
#' @param month an interger (1-12) indicating the month of interest
#' @return month_samples: sample points (spatial points) drawn from the given month
#' @return month_region_log: logistic probability raster of the given region in given month
#' @return month_N: the number of samples (integer) drawn from the given month proportional to the sample size from the reference month
#' @examples
#' # Example: draw a sample of 2000 mosquitos from Pennsylvania
#' sampling_results <- population_sampling(mosquito_predictions_6$raw,
#' mosquito_predictions_6$logis,
#' region=c("USA", "Pennsylvania"),
#' species=c("aedes", "aegypti"), 
#' N=2000,
#' month_logis_base=mosquito_predictions_6$logis, 
#' month=6, graph_logis=T, graph_raw=T) 
#' plot(sampling_results$month_region_log)
#' points(sampling_results$month_samples, cex=0.5, pch=16, col="blue")
#' title("Logistic Probabilities and Aedes Aegypti Population Sampling for PA")
# helper: population sampling given distribution
population_sampling <- function(maxent_predict_raw, maxent_predict_logis,
                                region=c("ITA"),
                                species, N, month_logis_base, month, 
                                graph_logis=F, graph_raw=F) {
  # get region's shape (sptialpolygons) and extract corresponding grids from the prediction raster
  target_region <- get_boundary(region)
  region_pred_raw <- mask(crop(maxent_predict_raw, target_region), target_region)
  
  region_pred_logis <- mask(crop(maxent_predict_logis, target_region), target_region)
  if (graph_logis) plot(region_pred_logis, main=paste("Logistic Probabilities -", paste(region, collapse=",")))
  region_pred_logis_base <- mask(crop(month_logis_base, target_region), target_region)
  N <- mean(getValues(region_pred_logis), na.rm=T) / mean(getValues(region_pred_logis_base), na.rm=T) * N
  
  # use raw probabilities from maxent model
  if (graph_raw) plot(region_pred_raw, main=paste("Raw Probabilities -", paste(region, collapse=",")))
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
  
  start_time <- Sys.time()
  sample_points<-probsel(probrast, N, x_radius, y_radius)
  end_time <- Sys.time()
  return (list(month_samples=sample_points, month_region_log=region_pred_logis, month_N=N))
}


#' @title run_simulation
#' @description run_simulation takes in a species name, a month, 
#' the number of samples to draw for that month, 
#' a region and optionally the environment layer resolution. 
#' It performs the simulation for the species' population sampling 
#' in the region across different months.
#' @param species_name a vector containing genus and species name
#' @param month an interger (1-12) indicating the month of interest
#' @param N number of species samples to draw in the reference month
#' @param region a vector of region names specifying a region of interest 
#' from national levels down to county levels
#' @param resolution the spatial resolution of environmental layers
#' (10, 5, or 2.5 (arc min))
#' @return log: a raster stack of logistic probability raster of the given region by month
#' @return sample: a list of sample points (spatial points) drawn from the given region by month
#' @return range: the range of logistic probabilities in the given region across all months
#' useful for setting a common color scale when plotting the probability rasters for all months
#' @return N: number of samples drawn from each month 
#' (proportionate to the reference month sample size 
#' by ratio of mean logistic probabilities)
#' @examples
#' region1 <- c("USA", "Pennsylvania")
#' species_name1 <- c("Ixodes", "ricinus", "Linnaeus")
#' species_name2 <- c("mus", "musculus")
#' result1 <- run_simulation(species_name1, month, N, region1, resolution)
#' result2 <- run_simulation(species_name2, month, N, region1, resolution)
#' 
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

#' @title plot_one_result
#' @description plot_one_result plots the logistic probabilities 
#' and population sample points with a color scale corresponding to 
#' a given probability range
#' @param result a returned value from run_simulation
#' @param max_range the maximum range of probabilities to be used for setting the color scale (breaks/legends)
#' @return None
#' @examples
#' plot_one_result(result1, max_range=range(unlist(result$range)))
#' 
plot_one_result <- function(result, max_range) {
  par(mfrow=c(4,3))
  for (i in 1:nlayers(result$log)) {
    layer <- result$log[[i]]
    plot(layer, breaks=seq(from=min(max_range), to=max(max_range), length.out=100), 
         col=rev(terrain.colors(99)),
         legend=F,
         main=paste(paste(result$species_name, collapse=" "), as.integer(result$N[i]),
                    "samples -",
                    paste(result$region, collapse=", "), months[i]))
    plot(layer, legend.only=T, 
         breaks=seq(from=min(max_range), to=max(max_range), length.out=100), 
         col=rev(terrain.colors(99)),
         axis.args=list(at=seq(max_range[1], max_range[2], length.out=5),
                        labels=round(seq(max_range[1], max_range[2], length.out=5),3),
                        cex.axis=0.6),
         legend.args=list(text='Logistic Probability', side=4, font=2, line=2.5, cex=0.8))
    plot(result$sample[[i]], add=T, pch=16, cex=0.5, col="blue")
  }
}

#' @title plot_results
#' @description plot_results plots the logistic probabilities 
#' and population sample points from multiple species simulations 
#' using a common color scale.
#' @param results a list of returned values from run_simulation
#' @return None
#' @examples 
#' plot_results(list(result1, result2))
#' 
plot_results <- function(results) {
  max_range <- c(0,0)
  for (i in 1:length(results)) {
    result = results[[i]]
    max_range <- range(max_range, range(unlist(result$range)))
  }
  for (i in 1:length(results)) {
    result = results[[i]]
    print(i)
    plot_one_result(result, max_range)
  }
}

