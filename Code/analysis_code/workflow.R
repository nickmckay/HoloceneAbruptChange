#analysis workflow
#run the excursion detector on the database
source("analysis_code/runExcursionDetector.R")

#load in data, create intermediate serializations
source("analysis_code/loadResults.R")

#begin analysis
source("analysis_code/SpatialAnalysisAndMapping.R")

#weighted time detection
source("analysis_code/timeAnalysis.R")

#make the figure with all the postage stamp maps
source("analysis_code/concatenateMaps.R")

#make the figure with all the postage stamp maps
source("analysis_code/createSpaceTimeFigures.R")

#make the figure with all the postage stamp maps
source("analysis_code/CreateTableOfExcursionResults.R")
