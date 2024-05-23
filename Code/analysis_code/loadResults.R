###############################################################################
# inhale and organize precipitation excursion detection over the Holocene

# Load libraries
library(tidyverse)
library(lipdR)

# Set the directories for input data, output data, and figures
# Directories for input data, output data, and figures
root_dir   <- 'analysis_code/'
output_dir <- file.path(root_dir,'new_output/')
figure_dir <- file.path(output_dir,'figures/')
data_dir <- file.path(output_dir,'data/')
processed_data_dir <- file.path(output_dir,'processed_data/')

filenames <- list.files(path=data_dir,pattern='rds',all.files=FALSE,full.names=TRUE)


#load in hydro data
filenamesHydro <- filenames[str_detect(tolower(filenames),pattern = "_hydro_")]

n_files_hydro <- length(filenamesHydro)
print(paste(' === FOUND FILES HYDRO:',n_files_hydro,'==='))

# Get years of the event detection
data1 <- readRDS(filenamesHydro[1])
data1$event.yr <- as.numeric(data1$event.yr)
event_years <- data1$event.yr
n_years <- length(event_years)


bigDataHydro <-  data1
# Loop through all of the files, loading data
for (file_num in 2:n_files_hydro) {
  cat(paste('Loading hydro file number:',file_num,'\r'))
  data <- readRDS(filenamesHydro[file_num])
  data$event.yr <- as.numeric(data$event.yr)

  bigDataHydro <- bind_rows(bigDataHydro,data)
}

bigDataHydro <- filter(bigDataHydro, !is.na(event_probability))

#Adjust by interpretation direction
hydroData <- googlesheets4::read_sheet(ss = "1eLQRxgUgmrb8OfFc2zOM0445KdaeA2lVAYPn1HNkfVE",sheet = "Hydro")

abHydro <- select(bigDataHydro,paleoData_TSid,contains("positive") | contains("negative")) %>%
  left_join(hydroData,by = "paleoData_TSid")

#align by direction.
isNegativeHydro <- which(startsWith(tolower(abHydro$hydroDirection),"n"))

bigDataHydro$pvalue_positive[isNegativeHydro] <- abHydro$pvalue_negative[isNegativeHydro]
bigDataHydro$null_probability_positive[isNegativeHydro] <- abHydro$null_probability_negative[isNegativeHydro]
bigDataHydro$event_probability_positive[isNegativeHydro] <- abHydro$event_probability_negative[isNegativeHydro]

bigDataHydro$pvalue_negative[isNegativeHydro] <- abHydro$pvalue_positive[isNegativeHydro]
bigDataHydro$null_probability_negative[isNegativeHydro] <- abHydro$null_probability_positive[isNegativeHydro]
bigDataHydro$event_probability_negative[isNegativeHydro] <- abHydro$event_probability_positive[isNegativeHydro]


if(!dir.exists(processed_data_dir)){
  dir.create(processed_data_dir)
}

save(bigDataHydro,file = paste0(processed_data_dir,"/bigDataHydro.RData"))

# repeat for temperature

filenamesTemp <- filenames[str_detect(tolower(filenames),pattern = "_temp_")]

n_files_temp <- length(filenamesTemp)
print(paste(' === FOUND FILES TEMP:',n_files_temp,'==='))

# Get years of the event detection
data1 <- readRDS(filenamesTemp[1])
data1$event.yr <- as.numeric(data1$event.yr)
event_years <- data1$event.yr
n_years <- length(event_years)


bigDataTemp <-  data1
# Loop through all of the files, loading data
for (file_num in 2:n_files_temp) {
  cat(paste('Loading temp file number:',file_num,'\r'))
  data <- readRDS(filenamesTemp[file_num])
  data$event.yr <- as.numeric(data$event.yr)

  bigDataTemp <- bind_rows(bigDataTemp,data)
}

bigDataTemp <- filter(bigDataTemp, !is.na(event_probability))

#temperature
tempData <- googlesheets4::read_sheet(ss = "1eLQRxgUgmrb8OfFc2zOM0445KdaeA2lVAYPn1HNkfVE",sheet = "Temp")

abTemp <- select(bigDataTemp,paleoData_TSid,contains("positive") | contains("negative")) %>%
  left_join(tempData,by = "paleoData_TSid")

#align by direction.
isNegativeTemp <- which(startsWith(tolower(abTemp$tempDirection),"n"))

bigDataTemp$pvalue_positive[isNegativeTemp] <- abTemp$pvalue_negative[isNegativeTemp]
bigDataTemp$null_probability_positive[isNegativeTemp] <- abTemp$null_probability_negative[isNegativeTemp]
bigDataTemp$event_probability_positive[isNegativeTemp] <- abTemp$event_probability_negative[isNegativeTemp]

bigDataTemp$pvalue_negative[isNegativeTemp] <- abTemp$pvalue_positive[isNegativeTemp]
bigDataTemp$null_probability_negative[isNegativeTemp] <- abTemp$null_probability_positive[isNegativeTemp]
bigDataTemp$event_probability_negative[isNegativeTemp] <- abTemp$event_probability_positive[isNegativeTemp]

if(!dir.exists(processed_data_dir)){
  dir.create(processed_data_dir)
}

save(bigDataTemp,file = paste0(processed_data_dir,"/bigDataTemp.RData"))


bigData <- bind_rows(bigDataHydro,bigDataTemp)
save(bigData,file = paste0(processed_data_dir,"/bigData.RData"))


###############################################################################
# gather data for spatial analysis and viz

#organize data for spatial weighting
allData <- select(bigData,lats = geo_latitude,lons = geo_longitude,TSids = paleoData_TSid) %>% distinct()

# run distance grid for temp and precip
wts=allData
grid.resolution = 3
radius.group=2000

allLon <- seq(-179.5, 179.5, grid.resolution)
allLat <- seq(-89.5, 89.5, grid.resolution)
countA <- 0
allDistances <- list()
totalRuns <- length(allLon)*length(seq(-89.5, 89.5, grid.resolution))
dist1 <- rep(NA, nrow(wts))



allDist <- function(lonNow,latNow){
    dist1 <- map2_dbl(wts$lons, wts$lats, \(x,y){ geosphere::distm(c(lonNow, latNow), c(x, y), fun = geosphere::distHaversine)/1000})

    weight1 <- 1 / (dist1 * 10/radius.group)
    Distances <- data.frame("TSid" = wts$TSids,
                            "distance" = dist1,
                            "weight" = weight1)
    Distances <- Distances[Distances$distance < radius.group,]
    return(list(lat=latNow, lon=lonNow, dist=Distances))
  }


future::plan(future::multisession,workers = 20) #this will run in parallel
#future::plan(future::sequential) # change to this if parallelization is causing problems

#get all lat/lon combinations:
allLL <- expand_grid(allLon,allLat)

allDistances <- furrr::future_map2(allLL$allLon,allLL$allLat,allDist,.progress = TRUE)




#save distance grid for future calcs
saveRDS(allDistances, file = paste0(processed_data_dir, "distanceGridAllHAC.rds"))











