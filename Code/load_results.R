###############################################################################
# inhale and organize precipitation excursion detection over the Holocene

# Load libraries
library(tidyverse)
library(lipdR)

# Set the directories for input data, output data, and figures
data_dir   = 'new_output/data/'
figure_dir = 'spatial_figs_final/'
processed_data_dir = 'processed_data_final/'


filenames <- list.files(path=data_dir,pattern='rds',all.files=FALSE,full.names=TRUE)

n_files <- length(filenames)
print(paste(' === FOUND FILES:',n_files,'==='))

# Get years of the event detection
data1 <- readRDS(filenames[1])
event_years <- data1$event.yr
n_years <- length(event_years)


bigData <-  readRDS(filenames[1])
# Loop through all of the files, loading data
for (file_num in 2:n_files) {
  cat(paste('Loading file number:',file_num,'\r'))
  data <- readRDS(filenames[file_num])

  bigData <- bind_rows(bigData,data)
}

bigData <- filter(bigData, !is.na(eventDetectionWithUncertainty))

if(!dir.exists(processed_data_dir)){
  dir.create(processed_data_dir)
}

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


future::plan(future::multisession,workers = 4) #this will run in parallel
#future::plan(future::sequential) # change to this if parallelization is causing problems

#get all lat/lon combinations:
allLL <- expand_grid(allLon,allLat)

allDistances <- furrr::future_map2(allLL$allLon,allLL$allLat,allDist,.progress = TRUE)




#save distance grid for future calcs
saveRDS(allDistances, file = paste0(processed_data_dir, "distanceGridAllHAC.rds"))











