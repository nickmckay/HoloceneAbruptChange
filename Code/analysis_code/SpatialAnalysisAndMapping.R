library(actR)
library(purrr)
library(tidyverse)
library(lipdR)

# Set the directories for input data, output data, and figures

root_dir   <- 'analysis_code/'
output_dir <- file.path(root_dir,'new_output/')
figure_dir <- file.path(output_dir,'figures/')
map_dir <- file.path(figure_dir,"allMaps/")
data_dir <- file.path(output_dir,'data/')
processed_data_dir <- file.path(output_dir,'processed_data/')

#load in preprocessed data
if(!exists("bigData")){
  load(paste0(processed_data_dir,"/bigData.RData"))
}

if(!exists("bigDataHydro")){
  load(paste0(processed_data_dir,"/bigDataHydro.RData"))
}

if(!exists("bigDataTemp")){
  load(paste0(processed_data_dir,"/bigDataTemp.RData"))
}

# Ready for analysis
allTimes <- sort(unique(bigData$event.yr))

allDistances <- readRDS(paste0(processed_data_dir, "distanceGridAllHAC.rds"))



# Set some parameters
all.age.to.plot <- seq(600,12000,by = 200)

weight.dist <- 1500

n.colors <- 11
grayCol <-  "#CCCCCC"
dryCol <- RColorBrewer::brewer.pal(n.colors,"BrBG")
dryCol[median(seq_along(dryCol))] <- grayCol
dryCol <- dryCol[rev(seq_along(dryCol))]

hotCol <- RColorBrewer::brewer.pal(n.colors,"RdBu")
hotCol[median(seq_along(hotCol))] <- grayCol

#set up scale transformation
fwd <- \(x) (x-.5)^3
inv <-  \(x) ifelse(x>=0,yes = x^(1/3),no = -1*abs(x)^(1/3)) +.5
cube <- scales::trans_new("cube",transform = fwd, inverse = inv)

makeMaps <- function(age.to.plot,resData,type){
  if(type == "hydro"){
    Col <- dryCol
  }else if(type == "temp"){
    Col <- hotCol
  }else{
    stop("unrecognized type")
  }


pval.grid.weight <- purrr::map(.x = allDistances,
                                           .f = actR::spatialSigTest,
                                           events = filter(resData,event.yr == age.to.plot),
                                           agg.method="robustNull",
                                           use.weights = TRUE,
                                           distance.cutoff = weight.dist,
                                     .progress = FALSE)




plot.name <-glue::glue("{age.to.plot} - {type} - Robust Null - {weight.dist} km")

figout <- plotSignificance(pval.grid.weight,
                           sigTestResults = filter(resData,event.yr == age.to.plot),
                           which.test = "pvalNet",
                           color.breaks = c(0.01,0.05,0.1,0.2),
                           restrict.sites = TRUE,
                           alpha.by.weight = TRUE,
                           cutoff.distance = weight.dist)

    plotOut <- figout$plot +
      scale_fill_stepsn(colours = Col,
                        breaks = c(.05,0.1,.2,.8,.9,.95),
                        trans = cube) +
      ggtitle(plot.name)

    ggsave(
      filename = paste0(map_dir, glue::glue("{plot.name}.pdf")),
      plot = plotOut,
      width = 15, height = 9
    )

}

future::plan(future::multisession,workers = 20) #this will run in parallel
options(future.rng.onMisuse = "ignore")


resDataHydro <- select(bigDataHydro,
                       paleoData_TSid,
                       datasetId,
                       geo_latitude,
                       geo_longitude,
                       event.yr,
                       starts_with("pvalue_"),
                       starts_with("event_probability"),
                       starts_with("null_probability"))

furrr::future_walk(all.age.to.plot,
                   .f = makeMaps,
                   resData = resDataHydro,
                   type = "hydro",
                   .progress = TRUE)

resDataTemp <- select(bigDataTemp,
                      paleoData_TSid,
                      datasetId,
                      geo_latitude,
                      geo_longitude,
                      event.yr,
                      starts_with("pvalue_"),
                      starts_with("event_probability"),
                      starts_with("null_probability"))

furrr::future_walk(all.age.to.plot,
                   .f = makeMaps,
                   resData = resDataTemp,
                   type = "temp",
                   .progress = TRUE)







# Calculate some statistics -----------------------------------------------


#number and fraction of sites with significant excursions at 4.2
thisTimeHydro <- filter(bigDataHydro,event.yr == 4200)
thisTimeTemp <- filter(bigDataTemp,event.yr == 4200)

nSigHydro42 <- sum(thisTimeHydro$pvalue_either < 0.05)
nSigTemp42 <- sum(thisTimeTemp$pvalue_either < 0.05)

fSigHydro42 <- nSigHydro42/sum(!is.na(thisTimeHydro$pvalue_either))
fSigTemp42 <- nSigTemp42/sum(!is.na(thisTimeTemp$pvalue_either))

sitesWithSigExcursions <- bind_rows(thisTimeHydro,thisTimeTemp) %>%
  filter(pvalue_either < 0.05) %>%
  select(dataSetName) %>%
  distinct()

sitesWithSigExcursions <- thisTimeHydro %>%
  filter(pvalue_either < 0.05) %>%
  select(dataSetName) %>%
  distinct()

nSigTemp42 <- sum(thisTimeTemp$pvalue_either < 0.05)

# 4.2 OR 4.0 --------------------------------------------------------------


#make 4.2 OR 4.0 maps
temp40or42 <- filter(bigDataTemp,event.yr %in% c(4200,4000))

hydro40or42 <- filter(bigDataHydro,event.yr %in% c(4200,4000))


pval.grid.weight.temp.40or42 <- map(.x = allDistances,
                                    .f = actR::spatialSigTest,
                                    events = temp40or42,
                                    agg.method="robustNull",
                                    use.weights = TRUE,
                                    distance.cutoff = weight.dist,
                                    .progress = TRUE)

pval.grid.weight.hydro.40or42 <- map(.x = allDistances,
                                     .f = actR::spatialSigTest,
                                     events = hydro40or42,
                                     agg.method="robustNull",
                                     use.weights = TRUE,
                                     distance.cutoff = weight.dist,
                                     .progress = TRUE)

figoutTemp.40or42 <- actR::plotSignificance(pval.grid.weight.temp.40or42,
                                      sigTestResults = temp40or42,
                                      which.test = "pvalNet",
                                      color.breaks = c(0.01,0.05,0.1,0.2),
                                      restrict.sites = TRUE,
                                      alpha.by.weight = TRUE,
                                      cutoff.distance = weight.dist)


plotOut <- figoutTemp.40or42$plot +
  scale_fill_stepsn(colours = hotCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube) +
  ggtitle("Temperature excursions at 4.2 or 4.0 ka")

ggsave(
  filename = paste0(figure_dir, glue::glue("figoutTemp.40or42.pdf")),
  plot = plotOut,
  width = 15, height = 9
)

figoutHydro.40or42 <- actR::plotSignificance(pval.grid.weight.hydro.40or42,
                                       sigTestResults = hydro40or42,
                                       which.test = "pvalNet",
                                       color.breaks = c(0.01,0.05,0.1,0.2),restrict.sites = TRUE,
                                       alpha.by.weight = TRUE,cutoff.distance = weight.dist)


plotOut <- figoutHydro.40or42$plot +
  scale_fill_stepsn(colours = dryCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube) +
  ggtitle("Hydroclimate excursions at 4.2 or 4.0 ka")

ggsave(
  filename = paste0(figure_dir, glue::glue("figoutHydro.40or42.pdf")),
  plot = plotOut,
  width = 15, height = 9
)


