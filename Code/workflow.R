#analysis workflow



#load some functions used throughout
source("spatial.R")

#run the excursion detector on the database
source("excursion_detector_all_records.R")

#load in data, create intermediate serializations
source("load_results.R")

#begin analysis
source("analysis.R")

#weighted time detection
source("weightedTimeDetection.R")

