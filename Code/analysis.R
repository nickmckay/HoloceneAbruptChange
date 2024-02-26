library(actR)
library(purrr)
library(tidyverse)
library(lipdR)
setwd('/HoloceneAbruptChange/')
source("spatial.R")
# Set the directories for input data, output data, and figures

data_dir   = '/new_output/data/'
figure_dir = '/HoloceneAbruptChange/spatial_figs_final/'
processed_data_dir = '/HoloceneAbruptChange/processed_data_final/'

#load in preprocessed data
if(!exists("bigData")){
  load(paste0(processed_data_dir,"/bigData.RData"))
}

#remove duplicate datasets
dupsToRemove <- c("CoppermineSaddleback.Nichols.1975",
                  "JR51-GC35.Bendle.2007",
                  "Kaartlamminsuo.Rankama.1988",
                  "AgeroedsMosse.Nilsson.1964","Boone.White.1986",
                  "MarionLake.Mathewes.1973",
                  "FontaineHenry.Lespez.2008","TurovaDacha.EPD",
                  "WilderSeebeimRuhestein.EPD","EtangdelaGruere.EPD",
                  "Basse-Ville.Voeltzel.1987",
                  "Rotsee.Ammann.1111","Lobsigensee.Ammann.1985",
                  "Moossee.Rey.2020",
                  "BuntesMoor.Weirich.1980","Alsopahok.Zatyko.2007",
                  "SaegistalseePollen.Wick.2003",
                  "PortageBog.Anderson.1980","PrazRodet.Shotyk.1997",
                  "LacduMontdOrgeSion.Colombaroli.2013",
                  "HieressurAmby.EPD","LeGrandLemps.Clerc.1988",
                  "Graham.Fuller.1997","High.Fuller.1997",
                  "PathLake.Neil.2014",
                  "HumberPond3.McAndrews.1989",
                  "LagodellAccesa.Colombaroli.2008","SaladaPequena.author.1111",
                  "StewartBog.Jimenez-Moreno.2008",
                  "Ferndale.Albert.1981","HeshangCave.Hu.2007",
                  "CuevadelDiablo.Bernal.2011",
                  "SO189-39KL.Mohtadi.2014")


bigData <- filter(bigData,!dataSetName %in% dupsToRemove)

#filter by autocorrelation
calcDetrendedAutocorrelation <- function(time,paleoData_values,event.yr,...){

  if(NCOL(time) > 1){
    time <- time[,1]
  }

  if(NCOL(paleoData_values) > 1){
    paleoData_values <- paleoData_values[,1]
  }

  hi <- which(between(time,left = event.yr - 750, event.yr + 750) & is.finite(paleoData_values))
  if(length(hi) < 3){
    return(NA)
  }
  time <- time[hi]
  vals <- paleoData_values[hi]

  co <-lm(vals ~ time)
  sv <- predict(co)
  det <- vals - sv

  ar1 <- acf(det,plot = FALSE)$acf[2]

  return(ar1)
}


ar <- pmap_dbl(bigData,calcDetrendedAutocorrelation,.progress = TRUE)

bigData$ar <- ar

bigData <- filter(bigData,ar < 0.98)

#first do the timeseries.
#all datasets either

hydroVars <- c("precipitation","effectivePrecipitation","streamflow")
tempVars <- c("temperature","growingDegreeDays")

#hydroclimate
googlesheets4::gs4_auth("nick.mckay2@gmail.com",cache = ".secret")
hydroData <- googlesheets4::read_sheet(ss = "1nFFTkHMY2ok0Ow1GiWfm4_DSDXstw1PH4tEac6I8svA",sheet = "Hydroclimate")

hi1 <- hydroData$interpretation1_variable %in% hydroVars
hi2 <- hydroData$interpretation2_variable %in% hydroVars
hi3 <- hydroData$interpretation3_variable %in% hydroVars

wi <- matrix(0,nrow = nrow(hydroData))
wi[hi2 & !hi1] <- 2
wi[hi3 & !hi1 &!hi2] <- 3
wi[hi1] <- 1

hydroData$hydroVariable <- NA
hydroData$hydroDirection <- NA

for(i in 1:nrow(hydroData)){
  hydroData$hydroVariable[i] <- hydroData[[paste0("interpretation",wi[i],"_variable")]][i]
  hydroData$hydroDirection[i] <- hydroData[[paste0("interpretation",wi[i],"_direction")]][i]
}

bigDataHydro <- dplyr::filter(bigData,paleoData_TSid %in% hydroData$paleoData_TSid)
abHydro <- select(bigDataHydro,paleoData_TSid,contains("above") | contains("below")) %>%
  left_join(hydroData,by = "paleoData_TSid")

#align by direction.
isNegativeHydro <- which(startsWith(tolower(abHydro$hydroDirection),"n"))

bigDataHydro$pvalue_above[isNegativeHydro] <- abHydro$pvalue_below[isNegativeHydro]
bigDataHydro$nullDetectionWithUncertaintyAbove[isNegativeHydro] <- abHydro$nullDetectionWithUncertaintyBelow[isNegativeHydro]
bigDataHydro$eventDetectionWithUncertaintyAbove[isNegativeHydro] <- abHydro$eventDetectionWithUncertaintyBelow[isNegativeHydro]
bigDataHydro$pvalue_below[isNegativeHydro] <- abHydro$pvalue_above[isNegativeHydro]
bigDataHydro$nullDetectionWithUncertaintyBelow[isNegativeHydro] <- abHydro$nullDetectionWithUncertaintyAbove[isNegativeHydro]
bigDataHydro$eventDetectionWithUncertaintyBelow[isNegativeHydro] <- abHydro$eventDetectionWithUncertaintyAbove[isNegativeHydro]


#temperature
tempData <- googlesheets4::read_sheet(ss = "1nFFTkHMY2ok0Ow1GiWfm4_DSDXstw1PH4tEac6I8svA",sheet = "Temperature")

ti1 <- tempData$interpretation1_variable %in% tempVars
ti2 <- tempData$interpretation2_variable %in% tempVars
ti3 <- tempData$interpretation3_variable %in% tempVars

wi <- matrix(0,nrow = nrow(tempData))
wi[ti2 & !ti1] <- 2
wi[ti3 & !ti1 &!ti2] <- 3
wi[ti1] <- 1

tempData$tempVariable <- NA
tempData$tempDirection <- NA

for(i in 1:nrow(tempData)){
  tempData$tempVariable[i] <- tempData[[paste0("interpretation",wi[i],"_variable")]][i]
  tempData$tempDirection[i] <- tempData[[paste0("interpretation",wi[i],"_direction")]][i]
}

bigDataTemp <- dplyr::filter(bigData,paleoData_TSid %in% tempData$paleoData_TSid)
abTemp <- select(bigDataTemp,paleoData_TSid,contains("above") | contains("below")) %>%
  left_join(tempData,by = "paleoData_TSid")

#align by direction.
isNegativeTemp <- which(startsWith(tolower(abTemp$tempDirection),"n"))

bigDataTemp$pvalue_above[isNegativeTemp] <- abTemp$pvalue_below[isNegativeTemp]
bigDataTemp$nullDetectionWithUncertaintyAbove[isNegativeTemp] <- abTemp$nullDetectionWithUncertaintyBelow[isNegativeTemp]
bigDataTemp$eventDetectionWithUncertaintyAbove[isNegativeTemp] <- abTemp$eventDetectionWithUncertaintyBelow[isNegativeTemp]
bigDataTemp$pvalue_below[isNegativeTemp] <- abTemp$pvalue_above[isNegativeTemp]
bigDataTemp$nullDetectionWithUncertaintyBelow[isNegativeTemp] <- abTemp$nullDetectionWithUncertaintyAbove[isNegativeTemp]
bigDataTemp$eventDetectionWithUncertaintyBelow[isNegativeTemp] <- abTemp$eventDetectionWithUncertaintyAbove[isNegativeTemp]


#ready for analysis

allTimes <- sort(unique(bigData$event.yr))

#allGlobal <- vector(mode = "list",length = length(allTimes))
hydroGlobal <- vector(mode = "list",length = length(allTimes))
tempGlobal <- vector(mode = "list",length = length(allTimes))

for(i in 1:length(allTimes)){
  cat(paste(round(100* i/length(allTimes)),"%\r"))
#thisTime <- filter(bigData,event.yr == allTimes[i])
thisTimeHydro <- filter(bigDataHydro,event.yr == allTimes[i])
thisTimeTemp <- filter(bigDataTemp,event.yr == allTimes[i])


#allGlobal[[i]] <- calculateMultiTestPvalue(thisTime) %>% as.data.frame()
hydroGlobal[[i]] <- calculateMultiTestPvalue(thisTimeHydro) %>% as.data.frame()
tempGlobal[[i]] <- calculateMultiTestPvalue(thisTimeTemp) %>% as.data.frame()


}

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

nSigTemp42 <- sum(thisTimeTemp$pvalue_either < 0.05)

# spatial -----------------------------------------------------------------
allDistances <- readRDS(paste0(processed_data_dir, "distanceGridAllHAC.rds"))

fdr <- function(pvals,qlevel=0.05,useEffN = FALSE){
  n <- length(pvals)
  s <- sort(pvals,index.return = TRUE)
  ar1 <- cor(pvals[-1],pvals[-length(pvals)])
  effN <- n * (1-ar1)/(1+ar1)
  sorted.pvals <- s$x
  sort.index <- s$ix
  if(is.numeric(useEffN)){
    effN <- useEffN
    useEffN <- TRUE
  }

  if(useEffN){
    pass <- (sorted.pvals < qlevel*seq(1/effN,1,length.out = n))
  }else{
    pass <- (sorted.pvals < qlevel*seq_len(n)/n)
  }

  to.return <- rep(FALSE,times = n)
  to.return[sort.index[pass]] <- TRUE
  return(to.return)
}


pullAndAverage <- function(ensNum,nulls,distWeights,areaWeights){
  #get null vals
  av <- map_dbl(nulls,ensNum)

  #calculate weights
  w <- distWeights + areaWeights

  w <- w/sum(w,na.rm = TRUE)

  return(sum(w * av,na.rm = TRUE))
}

gn2 <- function(x,type = "Above"){
  if(is.null(x[[paste0("pval",type)]])){
    return(NA)
  }else{
    if(is.na(x$pvalAbove)){
      return(NA)
    }else{
      return(getNulls(x$resultsUsed[[paste0("nullDetectionWithUncertainty",type)]],
                      nrows = nrow(x$resultsUsed),
                      n.ens = 1000,
                      weights = x$resultsUsed$weight))
    }
  }
}


getSpatialMeanData <- function(pval.grid.weight){
#pull out the relevant data
  pea <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$pvalAbove) == 1,x$pvalAbove,NA))
  peb <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$pvalBelow) == 1,x$pvalBelow,NA))
  pen <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$pvalNet) == 1,x$pvalNet,NA))
# aea <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$allEventAbove) == 1,x$allEventAbove,NA))
# aeb <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$allEventBelow) == 1,x$allEventBelow,NA))
# aen <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$allEventNet) == 1,x$allEventNet,NA))
aweight <- map_dbl(pval.grid.weight,\(x) ifelse(length(x$maxWeight) == 1,x$maxWeight,NA))
area_weight <- cos((map_dbl(pval.grid.weight,"lat") * pi)/180)

good <- which(is.finite(pea))

totalWeights <- aweight[good] + area_weight[good]
totalWeights <- totalWeights/sum(totalWeights)

#fdr
pa05fdr <-fdr(pea[good],0.05)
pb05fdr <-fdr(pea[good],0.05)
pn05fdr <-fdr(pea[good],0.05)


pa05 <- pea[good] < 0.05
pb05 <- peb[good] < 0.05
pn05 <- pen[good] < 0.05


wpea <- sum(pa05 * totalWeights)
wpeb <- sum(pb05 * totalWeights)
wpen <- sum(pn05 * totalWeights)

wpeafdr <- sum(pa05fdr * totalWeights)
wpebfdr <- sum(pb05fdr * totalWeights)
wpenfdr <- sum(pn05fdr * totalWeights)


tots <- data.frame(weightedPassAbove = wpea,
                      weightedPassBelow = wpeb,
                      weightedPassNet = wpen,
                   weightedPassAboveFDR = wpeafdr,
                   weightedPassBelowFDR = wpebfdr,
                   weightedPassNetFDR = wpenfdr,
                   ngridcells = length(good))

return(tots)
}

all.age.to.plot <- seq(600,12000,by = 200)

weight.dist <- 1500
future::plan(future::multisession,workers = 4) #this will run in parallel
options(future.rng.onMisuse = "ignore")


n.colors <- 11
grayCol <-  "#CCCCCC"
dryCol <- RColorBrewer::brewer.pal(n.colors,"BrBG")
dryCol[median(seq_along(dryCol))] <- grayCol
dryCol <- dryCol[rev(seq_along(dryCol))]

hotCol <- RColorBrewer::brewer.pal(n.colors,"RdBu")
hotCol[median(seq_along(hotCol))] <- grayCol

fwd <- \(x) (x-.5)^3
inv <-  \(x) ifelse(x>=0,yes = x^(1/3),no = -1*abs(x)^(1/3)) +.5
cube <- scales::trans_new("cube",transform = fwd, inverse = inv)


#make 4.2 or 4.0 maps
temp40or42 <- filter(bigDataTemp,event.yr %in% c(4200,4000))

hydro40or42 <- filter(bigDataHydro,event.yr %in% c(4200,4000))


pval.grid.weight.temp.40or42 <- furrr::future_map(.x = allDistances,
                                                  .f = spatialSigTest,
                                                  excursion.results = temp40or42,
                                                  agg.method="robustNull",
                                                  use.weights = TRUE,
                                                  distance.cutoff = weight.dist,
                                                  .progress = TRUE)

pval.grid.weight.hydro.40or42 <- furrr::future_map(.x = allDistances,
                                                  .f = spatialSigTest,
                                                  excursion.results = hydro40or42,
                                                  agg.method="robustNull",
                                                  use.weights = TRUE,
                                                  distance.cutoff = weight.dist,
                                                  .progress = TRUE)

figoutTemp.40or42 <- plotSignificance(pval.grid.weight.temp.40or42,
                               sigTestResults = temp40or42,
                               which.test = "pvalNet",
                               color.breaks = c(0.01,0.05,0.1,0.2),restrict.sites = TRUE,
                               alpha.by.weight = TRUE,cutoff.distance = weight.dist)


plotOut <- figoutTemp.40or42$plot +
  scale_fill_stepsn(colours = hotCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube) +
  ggtitle("Temperature excursions at 4.2 or 4.0 ka")

ggsave(
  filename = paste0(figure_dir, glue::glue("figoutTemp.40or42.pdf")),
  plot = plotOut,
  width = 15, height = 9
)

figoutHydro.40or42 <- plotSignificance(pval.grid.weight.hydro.40or42,
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


for(age.to.plot in all.age.to.plot){


pval.grid.weight.temp <- furrr::future_map(.x = allDistances,
                                        .f = spatialSigTest,
                                        excursion.results=filter(bigDataTemp,event.yr == age.to.plot),
                                        agg.method="robustNull",
                                        use.weights = TRUE,
                                      distance.cutoff = weight.dist,
                                        .progress = TRUE)

pval.grid.weight.hydro <- furrr::future_map(.x = allDistances,
                                           .f = spatialSigTest,
                                           excursion.results=filter(bigDataHydro,event.yr == age.to.plot),
                                           agg.method="robustNull",
                                           use.weights = TRUE,
                                           distance.cutoff = weight.dist,

n.colors <- 11
grayCol <-  "#CCCCCC"
dryCol <- RColorBrewer::brewer.pal(n.colors,"BrBG")
dryCol[median(seq_along(dryCol))] <- grayCol
dryCol <- dryCol[rev(seq_along(dryCol))]

hotCol <- RColorBrewer::brewer.pal(n.colors,"RdBu")
hotCol[median(seq_along(hotCol))] <- grayCol

fwd <- \(x) (x-.5)^3
inv <-  \(x) ifelse(x>=0,yes = x^(1/3),no = -1*abs(x)^(1/3)) +.5
cube <- scales::trans_new("cube",transform = fwd, inverse = inv)


pvalOptions <- c("pvalNet")
  for(po in pvalOptions){
    #now with weights
    figoutTemp <- plotSignificance(pval.grid.weight.temp,
                                   sigTestResults = filter(bigDataTemp,event.yr == age.to.plot),
                               which.test = po,
                               color.breaks = c(0.01,0.05,0.1,0.2),restrict.sites = TRUE,
                               alpha.by.weight = TRUE,cutoff.distance = weight.dist)

    if(po == "pvalNet"){
      toplot <- figoutTemp$plot +
        scale_fill_stepsn(colours = hotCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube)

    }else{
      toplot <- figoutTemp$plot
    }

    plot.name <-glue::glue("{age.to.plot} - Temperature - {po} - Robust Null - {weight.dist} km")

    toplot <- toplot + ggtitle(plot.name)


    ggsave(
      filename = paste0(figure_dir, glue::glue("{plot.name}.pdf")),
      plot = toplot,
      width = 15, height = 9
    )


    figoutHydro <- plotSignificance(pval.grid.weight.hydro,
                                   sigTestResults = filter(bigDataHydro,event.yr == age.to.plot),
                                   which.test = po,
                                   color.breaks = c(0.01,0.05,0.1,0.2),restrict.sites = TRUE,
                                   alpha.by.weight = TRUE,cutoff.distance = weight.dist)

    if(po == "pvalNet"){
      toplot <- figoutHydro$plot +
        scale_fill_stepsn(colours = dryCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube)

    }else{
      toplot <- figout$plot
    }

    plot.name <-glue::glue("{age.to.plot} - Hydroclimate - {po} - Robust Null - {weight.dist} km")

    toplot <- toplot + ggtitle(plot.name)


    ggsave(
      filename = paste0(figure_dir, glue::glue("{plot.name}.pdf")),
      plot = toplot,
      width = 15, height = 9
    )

 }
