
#estimate spatial weights
#go through each and select one from each site
equalAreaGridWeights <- function(lon,lat,spacing,dggs){
  if(!exists("dggs")){
    dggs   <- dggridR::dgconstruct(spacing=spacing, metric=TRUE,resround = "down")
  }

  df <- data.frame(lon = lon,
                   lat = lat,
                   grid = dggridR::dgGEO_to_SEQNUM(dggs,lon,lat)$seqnum)

  count <- df %>% group_by(grid) %>% summarise(count=n())

  dfb <- dplyr::left_join(df,count, by = "grid") %>%
    mutate(weight = 1/count)

  return(dfb)

}

fdr <- function(pvals,qlevel=0.05,useEffN = FALSE){
  n <- length(pvals)
  s <- sort(pvals,index.return = TRUE)

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


selectAndWeight <- function(thisTime,dist.cut = 500,dggs){
  #pick which TSid for this time
  #first get only the "best" from each dataset
  udsid <- unique(thisTime$datasetId)
  tu <- c()
  for(u in 1:length(udsid)){
    these <- filter(thisTime,datasetId == udsid[u])
    tu[u] <- these$paleoData_TSid[which.min(these$pvalue_either)]
  }




  thisTime <- filter(thisTime,paleoData_TSid %in% tu)
  #now calculate all distances for each remaining

  wdf <- equalAreaGridWeights(lon = thisTime$geo_longitude, lat = thisTime$geo_latitude,spacing = dist.cut,dggs)

  thisTime$weight = wdf$weight

  return(thisTime)
}


calculateBothMultiTestSignificance <- function(events1,events2,weights1 = NA,weights2 = NA,n.ens = 1000){

  #only include events with results
  events1$weights <- weights1

  events1 <- dplyr::filter(events1,!is.na(event_probability))

  weights1 <- events1$weights

  events2$weights <- weights2

  events2 <- dplyr::filter(events2,!is.na(event_probability))

  weights2 <- events2$weights




  if(nrow(events1) == 0 | nrow(events2) == 0){
    out <- list()
    out$pvalAbove <- NA
    out$pvalBelow <- NA
    out$pvalEither <- NA
    out$pvalBoth <- NA
    out$pvalNet <- NA
    return(out)
  }


  if(all(is.na(weights1))){
    weights1 <- rep(1,times = nrow(events1))
  }

  if(all(is.na(weights2))){
    weights2 <- rep(1,times = nrow(events2))
  }

  #make weights sum to 1
  weights1 <- weights1/sum(weights1)
  weights2 <- weights2/sum(weights2)


  #get the sum of the event testing:
  allEventPos1 <- sum(events1$event_probability_positive * weights1,na.rm = FALSE)
  allEventNeg1 <- sum(events1$event_probability_negative * weights1,na.rm = FALSE)
  allEventEither1 <- sum(events1$event_probability_either * weights1,na.rm = FALSE)
  allEventBoth1 <- sum(events1$event_probability_both * weights1,na.rm = FALSE)
  allEventNet1 <- allEventPos1 - allEventNeg1

  allEventPos2 <- sum(events2$event_probability_positive * weights2,na.rm = FALSE)
  allEventNeg2 <- sum(events2$event_probability_negative * weights2,na.rm = FALSE)
  allEventEither2 <- sum(events2$event_probability_either * weights2,na.rm = FALSE)
  allEventBoth2 <- sum(events2$event_probability_both * weights2,na.rm = FALSE)
  allEventNet2 <- allEventPos2 - allEventNeg2


  #compare the sum of the nulls
  #add in weights
  if(length(events1$null_probability) == 0){
    stop("empty null")
  }

  nullsPos1 <- getNulls(events1$null_probability_positive,nrow(events1),n.ens,weights = weights1)
  nullsNeg1 <- getNulls(events1$null_probability_negative,nrow(events1),n.ens,weights = weights1)
  nullsEither1 <- getNulls(events1$null_probability_either,nrow(events1),n.ens,weights = weights1)
  nullsBoth1 <- getNulls(events1$null_probability_both,nrow(events1),n.ens,weights = weights1)
  nullsNet1 <- nullsPos1 - nullsNeg1

  #compare the sum of the nulls
  #add in weights
  if(length(events2$null_probability) == 0){
    stop("empty null")
  }


  nullsPos2 <- getNulls(events2$null_probability_positive,nrow(events2),n.ens,weights = weights2)
  nullsNeg2 <- getNulls(events2$null_probability_negative,nrow(events2),n.ens,weights = weights2)
  nullsEither2 <- getNulls(events2$null_probability_either,nrow(events2),n.ens,weights = weights2)
  nullsBoth2 <- getNulls(events2$null_probability_both,nrow(events2),n.ens,weights = weights2)
  nullsNet2 <- nullsPos2 - nullsNeg2

  #calculate the product of the sums and use that
  allEventPos <- allEventPos1 * allEventPos2
  allEventNeg <- allEventNeg1 * allEventNeg2
  allEventEither <- allEventEither1 * allEventEither2
  allEventBoth <- allEventBoth1 * allEventBoth2
  allEventNet <- allEventPos - allEventNeg

  nullsPos <- apply(cbind(nullsPos1,nullsPos2),1,prod)
  nullsNeg <- apply(cbind(nullsNeg1,nullsNeg2),1,prod)
  nullsEither <- apply(cbind(nullsEither1,nullsEither2),1,prod)
  nullsBoth <- apply(cbind(nullsBoth1,nullsBoth2),1,prod)
  nullsNet <- nullsPos - nullsNeg

  out <- list()
  out$pvalPos <- kdePval(nullsPos, allEventPos)$pval
  out$pvalNeg <- kdePval(nullsNeg, allEventNeg)$pval
  out$pvalEither <- kdePval(nullsEither, allEventEither)$pval
  out$pvalBoth <- kdePval(nullsBoth, allEventBoth)$pval
  out$pvalNet <- kdePval(nullsNet,allEventNet)$pval

  out$allEventPos <- allEventPos
  out$allEventNeg <- allEventNeg
  out$allEventBoth <- allEventBoth
  out$allEventEither <- allEventEither
  out$allEventNet <- allEventNet

  out$clPos <- quantile(nullsPos,probs = c(.975))
  out$clNeg <- quantile(nullsNeg,probs = c(.975))
  out$clEither <- list(clEither = quantile(nullsEither,probs = c(.9,.95,.99)))
  out$clBoth <- list(clBoth = quantile(nullsBoth,probs = c(.9,.95,.99)))
  out$clNet <- list(clNet = quantile(nullsNet,probs = c(.005,.025,.05,.95,.975,.995)))


  return(out)
}

testFdr <- function(df,effN){
  df$fdrPos10 <- fdr(df$pvalPos,useEffN = effN,qlevel = .1)
  df$fdrNeg10 <- fdr(df$pvalNeg,useEffN = effN,qlevel = .1)
  df$fdrPos05 <- fdr(df$pvalPos,useEffN = effN,qlevel = .05)
  df$fdrNeg05 <- fdr(df$pvalNeg,useEffN = effN,qlevel = .05)

  df$fdrSigLevelPos <- "none"
  df$fdrSigLevelPos[df$fdrPos10] <- "0.10"
  df$fdrSigLevelPos[df$fdrPos05] <- "0.05"

  df$fdrSigLevelNeg <- "none"
  df$fdrSigLevelNeg[df$fdrNeg10] <- "0.10"
  df$fdrSigLevelNeg[df$fdrNeg05] <- "0.05"
  return(df)
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


calculateTimeSeries <- function(time,bigDataHydro,bigDataTemp){
  thisHydro <- dplyr::filter(bigDataHydro,event.yr == !!time) |>
    selectAndWeight(dist.cut = weight.distance,dggs)

  hydroGlobal <-  actR::calculateMultiTestSignificance(thisHydro,thisHydro$weight)  |>
    as_tibble()

  thisTemp <- dplyr::filter(bigDataTemp,event.yr == !!time) |>
    selectAndWeight(dist.cut = weight.distance,dggs)

  tempGlobal <-  actR::calculateMultiTestSignificance(thisTemp,thisTemp$weight)  |>
    as_tibble()

  bothGlobal <- calculateBothMultiTestSignificance(events1 = thisHydro,
                                                   events2 = thisTemp,
                                                   weights1 = thisHydro$weight,
                                                   weights2 = thisTemp$weight) |>
    as_tibble()

  return(list(hydroGlobal,tempGlobal,bothGlobal))
}
