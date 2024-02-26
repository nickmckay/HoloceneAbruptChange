# Load libraries
library(tidyverse)
library(foreach)
library(lipdR)
library(geoChronR)
library(actR)
library(dplyr)

# === SET UP ACTR ===

# Directories for input data, output data, and figures
root_dir   = 'analysis_code/'
# File here: https://lipdverse.org/HoloceneAbruptChange/0_11_0/HoloceneAbruptChange0_11_0.RData
data_file  = paste(root_dir,'path/to/lipd/data',sep='')
output_dir = paste(root_dir,'new_output/',sep='')

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

future::plan(future::multisession,workers = 5) #this will run in parallel
#future::plan(future::sequential) # change to this if parallelization is causing problems

# Avoid some warnings
options(future.rng.onMisuse = "ignore")

# Select the ages to test for excursions
ages_to_test = seq(600,12000,by=step)


# === LOAD AND PROCESS DATA ===

# Load the data as TS
load(data_file)

# Combine the tts and pts variables
ats <- bind_rows(tts,pts)

# Print all of the interpretation variables
table(ats$interpretation1_variable)

# Make a map
ATS <- as.lipdTs(ats)
plotSummaryTs(ATS)
dev.print(png,width=800,height=800,paste(output_dir,'figures/map_records.png',sep=''))

# Avoid some annoying data issues
ats2test <- select(ats,-starts_with("inCompilation"), -ends_with("Ensemble"))


# === RUN ACTR ===

# Start recording results in a file
cat(' === STARTING CALCULATIONS ===',file=paste(output_dir,'logfile_excursion.txt',sep=''),sep='\n')

# Loop through the records, testing the Holocene for excursions
n_records <- length(ats2test$paleoData_values)
n_tts     <- length(tts$paleoData_values)
record_num = 1
for (record_num in 1:n_records) {
  print(paste(' === STARTING CALCULATIONS FOR RECORD',record_num,'/',n_records,'==='))
  if (record_num <= n_tts) {type_txt <- 'tts'} else {type_txt <- 'pts'}

  # Get the record and some metadata
  record <- ats2test[record_num,]
  dataset_txt <- record$dataSetName
  tsid_txt    <- record$paleoData_TSid

  # Put the calculations in a tryCatch statement so that errors don't break the code
  p_value <- tryCatch({

    # Make a basic plot
    par(mfrow=c(1,1))
    plot(record$age[[1]],record$paleoData_values[[1]],type='b',xlim=c(12000,0))

    # Use the excursion detector on the record
    slidingResult <- actR:::detectExcursionSlidingWindow(record,
                                                         event.yr          = ages_to_test,
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

    # Create filenames
    output_txt  = paste('result_',type_txt,'_',dataset_txt,'_',tsid_txt,'_window',window_mean,'_step',step,sep='')
    filename_fig  = paste(output_dir,'figures/',output_txt,'.png',sep='')
    filename_data = paste(output_dir,'data/',output_txt,'.rds',sep='')

    # Save the output
    saveRDS(slidingResult,file=filename_data)

    # Plot the results
    actR:::plotExcursionSliding(slidingResult)
    dev.print(png,width=1800,height=800,filename_fig)

    # Print the status and update the logfile
    output_txt = paste('Success - Record ',record_num,'/',n_records,': type=',type_txt,', dataSetName=',dataset_txt,', TSid=',tsid_txt,sep='')
    print(output_txt)
    cat(output_txt,file=paste(output_dir,'logfile_excursion.txt',sep=''),sep='\n',append=TRUE)
  },
  error = function(e) {
    # Print the status and update the logfile
    output_txt = paste('ERROR: FAILED - Record ',record_num,'/',n_records,': type=',type_txt,', dataSetName=',dataset_txt,', TSid=',tsid_txt,sep='')
    print(output_txt)
    cat(output_txt,file=paste(output_dir,'logfile_excursion.txt',sep=''),sep='\n',append=TRUE)
  })
}
