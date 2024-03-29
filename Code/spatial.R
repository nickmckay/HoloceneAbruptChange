gaspari_cohn <- function(r, rmax) {
  # This function takes two arguments:
  # r: a vector of distances from the point of interest to other points
  # rmax: the maximum distance over which the function is defined

  r <- abs(r)

  # Define constants
  a = 2/3 * rmax
  b = 4/3 * rmax

  # Calculate the function
  out = rep(0, length(r))
  ind1 = which(r <= a)
  out[ind1] = 1 - 7/6 * (r[ind1]/a)^2 + 3/6 * (r[ind1]/a)^3

  ind2 = which(r > a & r <= b)
  out[ind2] = 1/3 * (2 - (r[ind2]/a))^3

  ind3 = which(r > b)
  out[ind3] = 0

  return(out)
}


#' Distance from grid cells to TS objects
#'
#' @param wts a lipd TS object
#' @param grid.resolution grid cell size resolution, eg. 1 deg x 1 deg
#' @param radius.group radius of filter for each grid cell in km
#'
#' @return distances
#' @export
#'
distanceGridPar <- function(wts=NULL,
                            grid.resolution = 10,
                            radius.group=1000){


  allLon <- seq(-179.5, 179.5, grid.resolution)
  allLat <- seq(-89.5, 89.5, grid.resolution)
  countA <- 0
  allDistances <- list()
  totalRuns <- length(allLon)*length(seq(-89.5, 89.5, grid.resolution))

  a1 <- foreach::foreach(lonNow = allLon) %:%
    foreach::foreach(latNow = allLat) %dopar% {
      countA <- countA + 1
      TS1 <- NULL
      dist1 <- apply(wts, 1, function(x) geosphere::distm(c(lonNow, latNow), c(x$geo_longitude, x$geo_latitude), fun = geosphere::distHaversine)/1000)
      weight1 <- 1 / (dist1 * 10/radius.group)
      Distances <- data.frame("TSid" = wts$paleoData_TSid,
                              "distance" = dist1,
                              "weight" = weight1)
      Distances <- Distances[Distances$distance < radius.group,]
      allDistances[[countA]] <- list(lat=latNow, lon=lonNow, dist=Distances)
    }


  b1 <- do.call('c', a1)

  return(b1)
}


#' Distance from grid cells to TS objects
#'
#' @param wts a lipd TS object
#' @param grid.resolution grid cell size resolution, eg. 1 deg x 1 deg
#' @param radius.group radius of filter for each grid cell in km
#'
#' @return distances
#' @export
#'
distanceGrid <- function(wts=NULL,
                         grid.resolution = 1,
                         radius.group=NULL){

  #prep run
  allLon <- seq(-179.5, 179.5, grid.resolution)
  countA <- 0
  allDistances <- list()
  totalRuns <- length(allLon)*length(seq(-89.5, 89.5, grid.resolution))

  for (lonNow in allLon){
    for (latNow in seq(-89.5, 89.5, grid.resolution)){
      countA <- countA + 1
      TS1 <- NULL
      Distances <- data.frame("TSid" = wts$paleoData_TSid,
                              "distance" = apply(wts, 1, function(x) distm(c(lonNow, latNow), c(x$geo_longitude, x$geo_latitude), fun = distHaversine)/1000)
      )
      allDistances[[countA]] <- list(lat=latNow, lon=lonNow, dist=Distances)
    }
    cat(round(countA/totalRuns,2)*100, "%\r")
  }

  allDistances2 <- allDistances[which(unlist(lapply(allDistances, function(x) !is.null(x))))]
  allDistances2 <- allDistances2[which(unlist(lapply(allDistances2, function(x) !is.na(x))))]

  if (!is.null(radius.group)){
    for (i in 1:length(allDistances2)){
      if (!is.null(allDistances2[[i]])){
        allDistances2[[i]]$dist <- allDistances2[[i]]$dist[which(unlist(allDistances2[[i]]$dist$distance < radius.group)),]
      }
    }
  }

  return(allDistances2)
}


#' Data Coverage Plot
#'
#' @param distance.grid as produced by distanceGrid()
#' @param color.breaks vector of breaks for coloring plot
#'
#' @return plot data
#' @export
#'
plotDataDensity <- function(distance.grid=NULL,
                            color.breaks = c(1,2,4,8,16)){

  plotData <- data.frame("TScount"=unlist(lapply(distance.grid, function(x) sum(!is.na(unique(x$dist$TSid))))),
                         "lat"=unlist(lapply(distance.grid, function(x) x$lat)),
                         "lon"=unlist(lapply(distance.grid, function(x) x$lon)))


  world <- map_data("world")


  p1 <- ggplot(data = plotData[plotData$TScount>0,],
               mapping = aes(x=lon,
                             y=lat,
                             fill=TScount)) +
    scale_fill_viridis_b(trans="log", breaks=color.breaks)+
    geom_tile() +
    geom_polygon(data=world, aes(x = long, y = lat, group = group), color="black", fill=alpha("white",0.1), inherit.aes = FALSE) +
    geom_point(inherit.aes = FALSE, data = wts, mapping = aes(y=geo_latitude,x=geo_longitude),color="white")+
    xlab("") +
    ylab("") +
    #title(main = "Data Coverage")+
    theme(legend.title = element_blank(),
          title = element_text("Data Coverage"))+
    coord_map(xlim=c(-180,180), ylim = c(-90,90), projection="mollweide")

  print(p1)

  returns <- list(plot = p1, plot.data=plotData)

  return(returns)
}


#' test for excusions accross many sime series and return error messages
#'
#' @param wts2test a lipd TS object
#' @param nulls numer of null hypothesis tests
#' @param eventWindows length of event window
#' @param refWindows length of reference window
#' @param sdCrit sd threshhold
#'
#' @return excursion test results
#' @export
#'
excursionTestHighRes <- function(wts2test=NULL, nulls=100,
                                 eventWindows=sample(rnorm(100,mean = 400, sd = 100), size = 1),
                                 refWindows=sample(rnorm(100,mean = 600, sd = 100), size = 1),
                                 sdCrit=sample(rnorm(100,mean = 2,sd = .25), size = 1)){



  excursionTrySingle <- function(wts, nulls=100){
    out <- tryCatch(
      {
        testResults <- detectExcursion(wts,
                                       n.ens = 10,
                                       null.hypothesis.n = nulls,
                                       event.yr = 8200,
                                       event.window = eventWindows,
                                       ref.window = refWindows,
                                       sig.num = sdCrit,
                                       min.vals = 4,
                                       simulate.time.uncertainty = FALSE,
                                       simulate.paleo.uncertainty = FALSE)
        testResults$empirical_pvalue
      },
      error=function(cond){
        print(cond)
        return(NA)
      }
    )
  }



  results <- data.frame(matrix(nrow = nrow(wts2test), ncol = 4))
  names(results) <- c("TSid","pval", "lat", "lon")
  for (jj in 1:nrow(wts2test)){
    results[jj,1] <- wts2test[jj,]$paleoData_TSid
    results[jj,2] <- suppressMessages(excursionTrySingle(wts2test[jj,], nulls = 100))
    if (is.na(results[jj,2])){
      results[jj,2] <- NA
    }else{
      if (results[jj,2]==0){
        results[jj,2] <- suppressMessages(excursionTrySingle(wts2test[jj,], nulls = 1000))
        if (results[jj,2]==0){
          results[jj,2] <- .0001
          # results[jj,2] <- suppressMessages(excursionTrySingle(wts2test[jj,], nulls = 10000))
          # if (results[jj,2]==0){
          #   results[jj,2] <- .00001
          #}
        }
      }
    }


    results[jj,3] <- wts2test[jj,]$geo_latitude
    results[jj,4] <- wts2test[jj,]$geo_longitude



    cat("\r\r", round(jj/nrow(wts2test),2)*100, "% ")
  }
  return(results)

}


#' sam as excursionTestHighRes but run in parallel (beta)
#'
#' @param wts2test lipd.ts
#' @param nulls numer of null hypothesis tests
#' @param eventWindows length of event window
#' @param refWindows length of ref window
#' @param sdCrit sd threshold for event
#'
#' @return excursion test results
excursionTestHighResPar <- function(wts2test=NULL, nulls=100,
                                    eventWindows=sample(rnorm(100,mean = 400, sd = 100), size = 1),
                                    refWindows=sample(rnorm(100,mean = 600, sd = 100), size = 1),
                                    sdCrit=sample(rnorm(100,mean = 2,sd = .25), size = 1)){


  resultsDF <- data.frame(matrix(nrow = nrow(wts2test), ncol = 4))
  names(resultsDF) <- c("TSid","pval", "lat", "lon")

  results <- foreach (jj = 1:nrow(wts2test), .errorhandling = 'pass') %dopar%  {
    res1 <- testResults <- actR::detectExcursion(wts2test[jj,],
                                                 n.ens = 10,
                                                 null.hypothesis.n = 100,
                                                 event.yr = 8200,
                                                 event.window = eventWindows,
                                                 ref.window = refWindows,
                                                 sig.num = sdCrit,
                                                 min.vals = 4,
                                                 simulate.time.uncertainty = FALSE,
                                                 simulate.paleo.uncertainty = FALSE)$empirical_pvalue
    if (is.na(res1)){
      NA
    }else if (res1 > 0){
      res1
    }else{
      if (res1==0){
        actR::detectExcursion(wts2test[jj,],
                              n.ens = 10,
                              null.hypothesis.n = 1000,
                              event.yr = 8200,
                              event.window = eventWindows,
                              ref.window = refWindows,
                              sig.num = sdCrit,
                              min.vals = 4,
                              simulate.time.uncertainty = FALSE,
                              simulate.paleo.uncertainty = FALSE)$empirical_pvalue

        if (res1==0){
          actR::detectExcursion(wts2test[jj,],
                                n.ens = 10,
                                null.hypothesis.n = 10000,
                                event.yr = 8200,
                                event.window = eventWindows,
                                ref.window = refWindows,
                                sig.num = sdCrit,
                                min.vals = 4,
                                simulate.time.uncertainty = FALSE,
                                simulate.paleo.uncertainty = FALSE)$empirical_pvalue
          if (res1==0){
            .00001
          }
        }
      }
    }
  }

  resultsDF$TSid <- wts2test$paleoData_TSid
  resultsDF$pval <- unlist(results)
  resultsDF$lat <- wts2test$geo_latitude
  resultsDF$lon <- wts2test$geo_longitude

  return(results)

}

getNulls <- function(allNulls,nrows,n.ens,weights = NA){

  if(all(is.na(weights))){
    weights <- rep(1,times = nrows)
  }

  if(length(allNulls) == 0){
stop("the null is missing")
  }

  weightMat <- rep(weights,times = length(allNulls[[1]])) %>%
    matrix(nrow = nrows, byrow = FALSE)

  nullMat <- matrix(unlist(allNulls),nrow = nrows,byrow = TRUE)

  nulls <- colSums(nullMat * weightMat)

  while(length(nulls) < n.ens){
    #shuffle column orders to increase sample size.
    shuffMat <- nullMat
    for(i in 1:nrow(nullMat)){
      shuffMat[i,] <- nullMat[i,sample.int(ncol(nullMat),replace = FALSE)]
    }
    nulls <- c(nulls,colSums(shuffMat * weightMat))
  }

  return(nulls)
}

calculateMultiTestPvalue <- function(excOut,weights = NA,n.ens = 1000){

  excOut <- dplyr::filter(excOut,!is.na(eventDetectionWithUncertaintyAbove))

  if(nrow(excOut) == 0){
    out <- list()
    out$pvalAbove <- NA
    out$pvalBelow <- NA
    out$pvalEither <- NA
    out$pvalBoth <- NA
    out$pvalNet <- NA
    return(out)
  }


  if(all(is.na(weights))){
    weights <- rep(1,times = nrow(excOut))
  }

  #make weights sum to 1
  weights <- weights/sum(weights)

  #get the sum of the event testing:
  allEventAbove <- sum(excOut$eventDetectionWithUncertaintyAbove * weights,na.rm = FALSE)
  allEventBelow <- sum(excOut$eventDetectionWithUncertaintyBelow * weights,na.rm = FALSE)
  allEventEither <- sum(excOut$eventDetectionWithUncertaintyEither * weights,na.rm = FALSE)
  allEventBoth <- sum(excOut$eventDetectionWithUncertaintyBoth * weights,na.rm = FALSE)
  allEventNet <- allEventAbove - allEventBelow


  #compare the sum of the nulls
  #add in weights
  if(length(excOut$nullDetectionWithUncertaintyAbove) == 0){
    stop("empty null")
  }

  nullsAbove <- getNulls(excOut$nullDetectionWithUncertaintyAbove,nrow(excOut),1000,weights = weights)
  nullsBelow <- getNulls(excOut$nullDetectionWithUncertaintyBelow,nrow(excOut),1000,weights = weights)
  nullsEither <- getNulls(excOut$nullDetectionWithUncertaintyEither,nrow(excOut),1000,weights = weights)
  nullsBoth <- getNulls(excOut$nullDetectionWithUncertaintyBoth,nrow(excOut),1000,weights = weights)
  nullsNet <- nullsAbove - nullsBelow

  out <- list()
  out$pvalAbove <- kdePval(nullsAbove, allEventAbove)$pval
  out$pvalBelow <- kdePval(nullsBelow, allEventBelow)$pval
  out$pvalEither <- kdePval(nullsEither, allEventEither)$pval
  out$pvalBoth <- kdePval(nullsBoth, allEventBoth)$pval
  out$pvalNet <- kdePval(nullsNet,allEventNet)$pval

  out$allEventAbove <- allEventAbove
  out$allEventBelow <- allEventBelow
  out$allEventBoth <- allEventBoth
  out$allEventEither <- allEventEither
  out$allEventNet <- allEventNet

  out$clAbove <- quantile(nullsAbove,probs = c(.975))
  out$clBelow <- quantile(nullsBelow,probs = c(.975))
  #out$clEither <- quantile(nullsEither,probs = c(.9,.95,.99))
  #out$clBoth <- quantile(nullsBoth,probs = c(.9,.95,.99))
  #out$clNet <- quantile(nullsNet,probs = c(.005,.025,.05,.95,.975,.995))


  return(out)
}

#' aggregate excursion tests spatially
#'
#' @param map.grid.cell from distanceGrid()
#' @param excursion.results from excursionTestHighRes()
#' @param agg.method fisher or sidak
#' @param min.pval min p value
#'
#' @return gridded significance results
#' @export
#'
  spatialSigTest <- function(map.grid.cell,
                             excursion.results,
                             agg.method="lancaster",
                             min.pval = 0.001,
                             distance.cutoff = 2000,
                             use.weights = TRUE){

    TSids <- map.grid.cell$dist %>%
      filter(distance < distance.cutoff)  %>%
 #     mutate(weight = 1 / (distance * 10/distance.cutoff)) %>%
      distinct()

    if(nrow(TSids) == 0){
      map.grid.cell$pvalAbove <- NA
      map.grid.cell$pvalBelow <- NA
      map.grid.cell$pvalBoth <- NA
      map.grid.cell$pvalEither <- NA
      return(map.grid.cell)
    }

    #get the data within range for this grid cell
    resultsNow <- filter(excursion.results,paleoData_TSid %in% TSids$TSid) %>%
      rename(TSid = paleoData_TSid) %>%
      left_join(TSids,by = "TSid") %>%
      mutate(weight = gaspari_cohn(distance,distance.cutoff))

    if(agg.method != "robustNull"){
      uniquesWeights <- unique(resultsNow$weight)
      for (ii in uniquesWeights){
        sameWeightsLoc <- which(resultsNow$weight == ii)
        sameWeights <- resultsNow[sameWeightsLoc,]
        if (nrow(sameWeights)>1){
          resultsNow[sameWeightsLoc[1],]$pvalue_either <- aggregation::sidak(sameWeights$pvalue_either)
          resultsNow[sameWeightsLoc[-1],]$pvalue_either <- NA
          resultsNow[sameWeightsLoc[1],]$pvalue_both <- aggregation::sidak(sameWeights$pvalue_both)
          resultsNow[sameWeightsLoc[-1],]$pvalue_both <- NA
          resultsNow[sameWeightsLoc[1],]$pvalue_above <- aggregation::sidak(sameWeights$pvalue_above)
          resultsNow[sameWeightsLoc[-1],]$pvalue_above <- NA
          resultsNow[sameWeightsLoc[1],]$pvalue_below <- aggregation::sidak(sameWeights$pvalue_below)
          resultsNow[sameWeightsLoc[-1],]$pvalue_below <- NA
        }
      }
      resultsNow <- filter(resultsNow,!is.na(pvalue_either))
    }else{
      #only include the data with the lowest p-values for each site
      resultsNow$smallest <- FALSE
      uniquesWeights <- unique(resultsNow$weight)
      for (ii in uniquesWeights){
        sameWeightsLoc <- which(resultsNow$weight == ii)
        sameWeights <- resultsNow[sameWeightsLoc,]
        if (nrow(sameWeights)>1){
          smallest <- which.min(sameWeights$pvalue_either)
          resultsNow$smallest[sameWeightsLoc[smallest]] <- TRUE
        }else{
          resultsNow$smallest[sameWeightsLoc] <- TRUE
        }
      }
      resultsNow <- filter(resultsNow,smallest == TRUE)
    }

    if(min.pval > 0){
      resultsNow <- resultsNow %>%
        rowwise() %>%
        mutate(across(starts_with("pvalue_"),\(x) max(x,min.pval)))
    }

    if(nrow(resultsNow)==0){
      map.grid.cell$pval <- NA
    }else{
      if (agg.method=="fisher"){
        map.grid.cell$pvalAbove <- aggregation::fisher(resultsNow$pvalue_above)
        map.grid.cell$pvalBelow <- aggregation::fisher(resultsNow$pvalue_below)
        map.grid.cell$pvalBoth <- aggregation::fisher(resultsNow$pvalue_both)
        map.grid.cell$pvalEither <- aggregation::fisher(resultsNow$pvalue_either)

      }else if (agg.method=="sidak"){
        map.grid.cell$pvalAbove <- aggregation::sidak(resultsNow$pvalue_above)
        map.grid.cell$pvalBelow <- aggregation::sidak(resultsNow$pvalue_below)
        map.grid.cell$pvalBoth <- aggregation::sidak(resultsNow$pvalue_both)
        map.grid.cell$pvalEither <- aggregation::sidak(resultsNow$pvalue_either)
        }else if (agg.method=="lancaster"){
          map.grid.cell$pvalAbove <- aggregation::lancaster(resultsNow$pvalue_above,weights = resultsNow$weight)
          map.grid.cell$pvalBelow <- aggregation::lancaster(resultsNow$pvalue_below,weights = resultsNow$weight)
          map.grid.cell$pvalBoth <- aggregation::lancaster(resultsNow$pvalue_both,weights = resultsNow$weight)
          map.grid.cell$pvalEither <- aggregation::lancaster(resultsNow$pvalue_either,weights = resultsNow$weight)
      }else if (agg.method=="robustNull"){
        if(use.weights){
          resultsNow$weight <- map_dbl(resultsNow$weight,\(x) max(0.0001,min(x,1)))
          map.grid.cell <- append(map.grid.cell,calculateMultiTestPvalue(resultsNow,weights = resultsNow$weight))
        }else{
          map.grid.cell <- append(map.grid.cell,calculateMultiTestPvalue(resultsNow))
        }
      }else{
        stop("agg.method must be one of: fisher, sidak, lancaster or robustNull")
      }

    }
    map.grid.cell$resultsUsed <- resultsNow
    if(nrow(resultsNow) > 0){
      map.grid.cell$sumWeight <- sum(resultsNow$weight)
      map.grid.cell$maxWeight <- max(resultsNow$weight)
      top2 <- c(sort(resultsNow$weight,decreasing = TRUE),0)[1:2]

      map.grid.cell$sumTop2Weight <- sum(top2)

    }else{
      map.grid.cell$sumWeight <- NA
      map.grid.cell$maxWeight <- NA
      map.grid.cell$sumTop2Weight <- NA

    }

    return(map.grid.cell)
  }



#' plot results of spatialSigTest()
#'
#' @param pval.grid from spatialSigTest()
#' @param sigTestResults from excursionTestHighRes()
#' @param color.breaks
#' @import ggplot2
#'
#' @return plot and plot data
#' @export
#'
plotSignificance <- function(pval.grid=NULL,
                             sigTestResults,
                             color.breaks = c(.001,.01,.05,.1,.2),
                             which.test = "pval",
                             restrict.sites = TRUE,
                             alpha.by.weight = TRUE,
                             cutoff.distance = 1500,
                             projection = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"){


  # plotData <- data.frame("pval"= unlist(lapply(pval.grid, function(x) x$pval)),
  #                        "lat"=unlist(lapply(pval.grid, function(x) x$lat)),
  #                        "lon"=unlist(lapply(pval.grid, function(x) x$lon)))

  getPvals <- function(grid,which.test){
    if(is.null(grid[[which.test]])){
      return(NA)
    }else{
      return(grid[[which.test]])
    }
  }

  pvals <- map_dbl(pval.grid,getPvals,which.test)
  lat <- map_dbl(pval.grid,pluck,"lat")
  lon <- map_dbl(pval.grid,pluck,"lon")
  good <- which(!is.na(pvals))

  weight <- map_dbl(pval.grid[good],pluck,"sumTop2Weight")

  weight[weight < 0] <- 0



  plotData <- data.frame(pval = pvals[good],lat = lat[good],lon = lon[good],weight = weight)

  pvalOptions <- c("pvalAbove","pvalBelow","pvalEither","pvalNet")


  handleNet <- function(above,below){
    if(above > 0.5 & below > 0.5){
      return(0.5)
    }

    out <- ifelse(above > below, 1 - below,  above)

    return(out)
  }

  sigTestResults$pvalPlot <- dplyr::case_when(
    which.test == "pvalAbove" ~  sigTestResults$pvalue_above,
    which.test == "pvalBelow" ~  sigTestResults$pvalue_below,
    which.test == "pvalEither" ~  sigTestResults$pvalue_either,
    which.test == "pvalBoth" ~  sigTestResults$pvalue_both,
    which.test == "pvalNet" ~  purrr::map2_dbl(sigTestResults$pvalue_above, sigTestResults$pvalue_below, handleNet)
  )

  #optionally only plot the most significant
  if(restrict.sites){
    sigTestResults$plot <- FALSE

    udsn <- unique(sigTestResults$datasetId)
    for(u in udsn){
      w <- which(sigTestResults$datasetId == u)
      if(length(w) == 1){
        sigTestResults$plot[w] <- TRUE
      }else{
        if(which.test != "pvalNet"){
          tp <- w[which.min(sigTestResults$pvalPlot[w])]
          sigTestResults$plot[tp] <- TRUE
        }else{
          tp <- w[which.max(abs(sigTestResults$pvalPlot[w]-0.5))]
          sigTestResults$plot[tp] <- TRUE
        }
      }

    }
    sigTestResults <- filter(sigTestResults,plot == TRUE)

  }



  #sort by significance
  if(which.test == "pvalNet"){
  sigTestResults <- sigTestResults %>%
    arrange(abs(pvalPlot - 0.5))
  }else{
    sigTestResults <- sigTestResults %>%
      arrange(desc(pvalPlot))
  }

  cut <- 1.5
  world <- map_data("world") %>%
    filter(abs(long) < 180-cut,
           abs(lat) < 90-cut)

  # world_rob <- project(cbind(world$long,world$lat), proj = projection) %>% as.data.frame()
  # names(world_rob) <- c("x_proj","y_proj")
  # world <- bind_cols(world, world_rob)
  # #crs_goode <- "+proj=igh"

  toPlot <- filter(plotData,
                   !is.na(pval),
                   abs(lat) < 90-cut,
                   abs(lon) < 180-cut)

  #get corners
  # latSpacing <- 3
  # lonSpacing <- 3
  #
  # toPlot$xmin <- toPlot$lon-lonSpacing/2
  # toPlot$xmin[toPlot$xmin < -180] <- -180
  #
  # toPlot$xmax <- toPlot$lon+lonSpacing/2
  # toPlot$xmax[toPlot$xmax > 180] <- 180
  #
  # toPlot$ymin <- toPlot$lat-latSpacing/2
  # toPlot$ymin[toPlot$ymin < -90] <- -90
  #
  # toPlot$ymax <- toPlot$lat+latSpacing/2
  # toPlot$ymax[toPlot$ymax > 90] <- 90
  #
  #
  # toPlot_rob_min <- project(cbind(toPlot$xmin,toPlot$ymin),
  #                       proj = projection) %>% as.data.frame()
  # toPlot_rob_max <- project(cbind(toPlot$xmax,toPlot$ymax),
  #                           proj = projection) %>% as.data.frame()
  #
  # toPlot_rob <- bind_cols(toPlot_rob_min,toPlot_rob_max)
  #
  # names(toPlot_rob) <- c("xmin_proj","ymin_proj","xmax_proj","ymax_proj")
  # toPlot <- bind_cols(toPlot, toPlot_rob)

  p1 <- ggplot(data = toPlot,mapping = aes(x = lon,y = lat, fill = pval))


  if(alpha.by.weight){
    p1 <- p1 + geom_tile(aes(alpha = weight))
  }else{
    p1 <- p1 + geom_tile()
  }

  p1 <- p1 +
    geom_polygon(data=world, aes(x = long, y = lat, group = group), color="gray50", fill=alpha("white",0.1), inherit.aes = FALSE) +
    #geom_point(inherit.aes = FALSE, data = wts2test, mapping = aes(y=geo_latitude,x=geo_longitude),color="white")+
    geom_point(inherit.aes = FALSE, data=sigTestResults, aes(x=geo_longitude, y=geo_latitude, fill=pvalPlot), size = 2, shape=21, color="black")+
    scale_fill_viridis_b(breaks=color.breaks)+
    #geom_point(inherit.aes = FALSE, data=sigTestResults[complete.cases(sigTestResults),], aes(x=lon, y=lat), color="black", size=1.5, shape=8)+
    scale_x_continuous("",breaks = seq(-180,180,by = 60)) +
    scale_y_continuous("",breaks = c(-87.5,seq(-60,60,by = 30),87.5)) +
    #title(main = "Excursion Significance")+
    theme(legend.title = element_blank())+
    cowplot::theme_minimal_grid()+
    #coord_map(xlim=c(-180,180),
              # ylim = c(-90,90),
              # projection="robinson")
    coord_sf(xlim=c(-180,180),
             ylim = c(-90,90),
             crs = projection,
             expand = TRUE,
             datum = sf::st_crs(4326),
             default_crs = sf::st_crs(4326))

  #print(p1)

  returns <- list(plot = p1, plot.data=plotData)

  return(returns)
}











plotSignificanceOrangePeel <- function(pval.grid=NULL,
                                       sigTestResults,
                                       color.breaks = c(.001,.01,.05,.1,.2)){

  #plot projection from: https://wilkelab.org/practicalgg/articles/goode.html

  plotData <- data.frame("pval"=unlist(lapply(pval.grid, function(x) x$pval)),
                         "lat"=unlist(lapply(pval.grid, function(x) x$lat)),
                         "lon"=unlist(lapply(pval.grid, function(x) x$lon)))

  world_sf <- st_as_sf(getMap(resolution = "low"))

  # projection outline in long-lat coordinates
  lats <- c(
    90:-90, # right side down
    -90:0, 0:-90, # third cut bottom
    -90:0, 0:-90, # second cut bottom
    -90:0, 0:-90, # first cut bottom
    -90:90, # left side up
    90:0, 0:90, # cut top
    90 # close
  )
  longs <- c(
    rep(180, 181), # right side down
    rep(c(80.01, 79.99), each = 91), # third cut bottom
    rep(c(-19.99, -20.01), each = 91), # second cut bottom
    rep(c(-99.99, -100.01), each = 91), # first cut bottom
    rep(-180, 181), # left side up
    rep(c(-40.01, -39.99), each = 91), # cut top
    180 # close
  )

  goode_outline <-
    list(cbind(longs, lats)) %>%
    st_polygon() %>%
    st_sfc(
      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    )

  # now we need to work in transformed coordinates, not in long-lat coordinates
  goode_outline <- st_transform(goode_outline, crs = crs_goode)

  # get the bounding box in transformed coordinates and expand by 10%
  xlim <- st_bbox(goode_outline)[c("xmin", "xmax")]*1.1
  ylim <- st_bbox(goode_outline)[c("ymin", "ymax")]*1.1

  # turn into enclosing rectangle
  goode_encl_rect <-
    list(
      cbind(
        c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
      )
    ) %>%
    st_polygon() %>%
    st_sfc(crs = crs_goode)

  # calculate the area outside the earth outline as the difference
  # between the enclosing rectangle and the earth outline
  goode_without <- st_difference(goode_encl_rect, goode_outline)

  d_points <- st_as_sf(x =data.frame(lat=sigTestResults$lat, lon=sigTestResults$lon), coords = c("lon","lat"), crs=4326) %>%
    st_transform(crs=crs_goode)

  # d_points <- data.frame(long = c(-110, -103, -84),
  #                        lat  = c(45, 40, 41)) %>%
  #   st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  #   st_transform(crs = crs_use)

  plotData <- plotData[!is.na(plotData$pval),]
  d_tiles <- st_as_sf(x =data.frame(lat=plotData$lat, lon=plotData$lon), coords = c("lon","lat"), crs=4326) %>%
    st_transform(crs=crs_goode)

  d_tiles <- d_tiles %>%
    mutate(lon = unlist(map(d_tiles$geometry,1)),
           lat = unlist(map(d_tiles$geometry,2)))

  d_tiles$pval <- plotData$pval


  p1 <- ggplot(world_sf) +
    geom_sf(size = 0.5/.pt, fill=alpha("gray",0.1)) +
    geom_sf(data = goode_without, fill = "white", color = "NA") +
    geom_tile(aes(x=lon,y=lat,color= plotData$pval, fill=plotData$pval), data = d_tiles, linewidth=3)+
    geom_sf(aes(fill = sigTestResults$pval), data = st_jitter(d_points, factor = .005),
            size = 4, shape=21, color="red")+

    scale_fill_viridis_b(breaks=color.breaks)+
    scale_color_viridis_b(breaks=color.breaks)+
    coord_sf(crs = crs_goode) +
    theme_minimal_grid() +
    xlab("") +
    ylab("")
  #geom_tile(inherit.aes = FALSE, data=plotData[!is.na(plotData$pval),], mapping = aes(x=lon, y=lat, fill=pval, geometry = geometry), stat = "sf_coordinates") +
  #geom_polygon(data=world, aes(x = long, y = lat, group = group), color="black", fill=alpha("white",0.1), inherit.aes = FALSE) +
  #geom_point(inherit.aes = FALSE, data = wts2test, mapping = aes(y=geo_latitude,x=geo_longitude),color="white")+
  # geom_sf(data = d_points,
  #         color = "red", size = 5)
  #scale_fill_viridis_b(breaks=color.breaks)


  print(p1)

  returns <- list(plot = p1, plot.data=plotData)

  return(returns)
}
