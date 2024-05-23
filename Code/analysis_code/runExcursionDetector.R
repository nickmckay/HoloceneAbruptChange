# Load libraries
library(tidyverse)
library(lipdR)
library(geoChronR)
library(actR)
library(dplyr)

# === SET UP ACTR ===

# Directories for input data, output data, and figures
root_dir   <- 'analysis_code/'
output_dir <- file.path(root_dir,'new_output/')
figure_dir <- file.path(output_dir,'figures/')
data_dir <- file.path(output_dir,'data/')
map_dir <- file.path(figure_dir,"allMaps/")


if(!dir.exists(root_dir)){
  dir.create(root_dir)
}
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
if(!dir.exists(figure_dir)){
  dir.create(figure_dir)
}
if(!dir.exists(data_dir)){
  dir.create(data_dir)
}
if(!dir.exists(map_dir)){
  dir.create(map_dir)
}


# Options for actR
window_mean         <- 400
step                <- 200
n_ens               <- 100
n_ens_null          <- 50
n_points_per_period <- 4

# Sample from the parameter uncertainty
set.seed(1) # Use a seed so that these values are always the same
eventWindows <- rnorm(n_ens,mean = window_mean, sd = 100)
refWindows   <- eventWindows
sdCrit       <- rnorm(n_ens,mean = 2, sd = .25)

# Make a figure of the random numbers
par(mfrow=c(3,1))
hist(eventWindows,main='Length of event windows')
hist(refWindows,main='Length of reference windows')
hist(sdCrit,main='Number of standard deviations to exceed')
dev.print(png,width=800,height=800,paste(output_dir,'figures/parameter_hist.png',sep=''))

future::plan(future::multisession,workers = 20) #this will run in parallel
#future::plan(future::sequential) # change to this if parallelization is causing problems

# Avoid some warnings
options(future.rng.onMisuse = "ignore")

# Select the ages to test for excursions
ages_to_test = seq(600,12000,by=step)

# === LOAD AND PROCESS DATA ===
# Get a list of data to download
googlesheets4::gs4_deauth()
hydroSites <- googlesheets4::read_sheet("1eLQRxgUgmrb8OfFc2zOM0445KdaeA2lVAYPn1HNkfVE",sheet = "Hydro")

#download the corresponding data
hydro <- readLipd(hydroSites)
hydroTs <- as.lipdTsTibble(hydro)

hydroTs <- hydroTs |>
  dplyr::filter(paleoData_TSid %in% hydroSites$paleoData_TSid)

# Avoid some annoying data issues
hydroTs2test <- select(hydroTs,-starts_with("inCompilation"), -ends_with("Ensemble"))


#now prep temperature

tempSites <- googlesheets4::read_sheet("1eLQRxgUgmrb8OfFc2zOM0445KdaeA2lVAYPn1HNkfVE",sheet = "Temp")

#download the corresponding data
temp <- readLipd(tempSites)
tempTs <- as.lipdTsTibble(temp)

tempTs <- tempTs |>
  dplyr::filter(paleoData_TSid %in% tempSites$paleoData_TSid)

# Avoid some annoying data issues
tempTs2test <- select(tempTs,-starts_with("inCompilation"), -ends_with("Ensemble"))

# === RUN ACTR ===

# Start recording results in a file
cat(' === STARTING CALCULATIONS ===',file=paste(output_dir,'logfile_excursion.txt',sep=''),sep='\n')

runRecord <- function(record_num,type_txt,ts2test){
  n_records <- nrow(ts2test)

    print(paste(' === STARTING CALCULATIONS FOR RECORD',record_num,'/',n_records,'==='))

  # Get the record and some metadata
  record <- ts2test[record_num,]
  dataset_txt <- record$dataSetName
  tsid_txt    <- record$paleoData_TSid

  # Create filenames
  output_txt  = paste('result_',type_txt,'_',dataset_txt,'_',tsid_txt,'_window',window_mean,'_step',step,sep='')
  filename_fig  = paste(output_dir,'figures/',output_txt,'.pdf',sep='')
  filename_data = paste(output_dir,'data/',output_txt,'.rds',sep='')

  #check if the file exists
  if(file.exists(filename_data)){
    return(NULL) #then return nothing and move on
  }

  # Put the calculations in a tryCatch statement so that errors don't break the code
  p_value <- tryCatch({

    # Use the excursion detector on the record
    slidingResult <- detectExcursionSlidingWindow(record,
                                                  event.yr.vec      = ages_to_test,
                                                  event.window      = eventWindows,
                                                  ref.window        = refWindows,
                                                  sig.num           = sdCrit,
                                                  min.vals          = n_points_per_period,
                                                  n.consecutive     = 2,
                                                  exc.type          = "either",
                                                  n.ens             = n_ens,
                                                  null.hypothesis.n = n_ens_null,
                                                  simulate.time.uncertainty  = FALSE,
                                                  simulate.paleo.uncertainty = FALSE)



    # Save the output
    saveRDS(slidingResult,file=filename_data)

    # Plot the results
    slide_plot <- actR::plotExcursionSliding(slidingResult)
    ggsave(plot = slide_plot,path = filename_fig)

    # Print the status and update the logfile
    output_txt = paste('Success - Record ',record_num,'/',n_records,': type=',type_txt,', dataSetName=',dataset_txt,', TSid=',tsid_txt,sep='')
    #print(output_txt)
    cat(output_txt,file=paste(output_dir,'logfile_excursion.txt',sep=''),sep='\n',append=TRUE)
  },
  error = function(e) {
    # Print the status and update the logfile
    output_txt = paste('ERROR: FAILED - Record ',record_num,'/',n_records,': type=',type_txt,', dataSetName=',dataset_txt,', TSid=',tsid_txt,sep='')
    #print(output_txt)
    cat(output_txt,file=paste(output_dir,'logfile_excursion.txt',sep=''),sep='\n',append=TRUE)
  })
}


#Process all the hydro data. This will take awhile.
furrr::future_walk(1:nrow(hydroTs2test),runRecord,type_txt = "hydro",ts2test = hydroTs2test,.progress = TRUE)

#now the temperature data. This will take awhile.
furrr::future_walk(1:nrow(tempTs2test),runRecord,type_txt = "temp",ts2test = tempTs2test,.progress = TRUE)

