#estimate spatial weights
#go through each and select one from each site
source("equalAreaWeights.R")

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

weight.distance <- 1500

dggs   <- dggridR::dgconstruct(spacing=weight.distance, metric=TRUE)


allTimes <- sort(unique(bigData$event.yr))

#allGlobal <- vector(mode = "list",length = length(allTimes))
hydroGlobal <- vector(mode = "list",length = length(allTimes))
tempGlobal <- vector(mode = "list",length = length(allTimes))

allTsidHydro <- allTsidTemp <- c()

for(i in 1:length(allTimes)){
  cat(paste(round(100* i/length(allTimes)),"%\r"))
  thisTimeHydro <- filter(bigDataHydro,event.yr == allTimes[i])

  allTsidHydro <- unique(c(allTsidHydro,thisTimeHydro$paleoData_TSid))

  thisTimeTemp <- filter(bigDataTemp,event.yr == allTimes[i])

  allTsidTemp <- unique(c(allTsidTemp,thisTimeTemp$paleoData_TSid))

  thisTimeHydro <- selectAndWeight(thisTimeHydro,dist.cut = weight.distance,dggs)

  thisTimeTemp <- selectAndWeight(thisTimeTemp,dist.cut = weight.distance,dggs)

  hydroGlobal[[i]] <- calculateMultiTestPvalue(thisTimeHydro,weights = thisTimeHydro$weight) %>% as.data.frame()
  tempGlobal[[i]] <- calculateMultiTestPvalue(thisTimeTemp,weights = thisTimeTemp$weight) %>% as.data.frame()


}

hgdf <- list_rbind(hydroGlobal)
tgdf <- list_rbind(tempGlobal)

hgdf$age <- allTimes
tgdf$age <- allTimes
#

effN <- nrow(tgdf)/2 #account for non-independence of tests


tgdf$fdrAbove10 <- fdr(tgdf$pvalAbove,useEffN = effN,qlevel = .1)
tgdf$fdrBelow10 <- fdr(tgdf$pvalBelow,useEffN = effN,qlevel = .1)
hgdf$fdrAbove10 <- fdr(hgdf$pvalAbove,useEffN = effN,qlevel = .1)
hgdf$fdrBelow10 <- fdr(hgdf$pvalBelow,useEffN = effN,qlevel = .1)

tgdf$fdrAbove05 <- fdr(tgdf$pvalAbove,useEffN = effN,qlevel = .05)
tgdf$fdrBelow05 <- fdr(tgdf$pvalBelow,useEffN = effN,qlevel = .05)
hgdf$fdrAbove05 <- fdr(hgdf$pvalAbove,useEffN = effN,qlevel = .05)
hgdf$fdrBelow05 <- fdr(hgdf$pvalBelow,useEffN = effN,qlevel = .05)

tgdf$fdrSigLevelAbove <- "none"
tgdf$fdrSigLevelAbove[tgdf$fdrAbove10] <- "0.10"
tgdf$fdrSigLevelAbove[tgdf$fdrAbove05] <- "0.05"

tgdf$fdrSigLevelBelow <- "none"
tgdf$fdrSigLevelBelow[tgdf$fdrBelow10] <- "0.10"
tgdf$fdrSigLevelBelow[tgdf$fdrBelow05] <- "0.05"

hgdf$fdrSigLevelAbove <- "none"
hgdf$fdrSigLevelAbove[hgdf$fdrAbove10] <- "0.10"
hgdf$fdrSigLevelAbove[hgdf$fdrAbove05] <- "0.05"

hgdf$fdrSigLevelBelow <- "none"
hgdf$fdrSigLevelBelow[hgdf$fdrBelow10] <- "0.10"
hgdf$fdrSigLevelBelow[hgdf$fdrBelow05] <- "0.05"



patternSpacing <- 0.03
aboveColor <- RColorBrewer::brewer.pal(5,"BrBG")[5]
belowColor <- RColorBrewer::brewer.pal(5,"BrBG")[1]
lineColorAbove <- RColorBrewer::brewer.pal(7,"BrBG")[6]
lineColorBelow <- RColorBrewer::brewer.pal(7,"BrBG")[2]


#make a plot. Any excursion, either direction
hydroDirectional <- ggplot(hgdf) +
  geom_col_pattern(aes(x = age,y = allEventAbove,fill = "Positive",pattern = fdrSigLevelAbove), pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  geom_col_pattern(aes(x = age,y = -allEventBelow,fill = "Negative",pattern = fdrSigLevelBelow),pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "bottom") +
  geom_line(aes(x = age,y = -clBelow), color = lineColorBelow,linetype = 1) +
  geom_line(aes(x = age,y = clAbove), color = lineColorAbove,linetype = 1) +
  geom_label(aes(x = 12300,y = clAbove[58] ),label = "95% cl",color = lineColorAbove) +
  geom_label(aes(x = 12300,y = -clBelow[58]),label = "95% cl",color = lineColorBelow) +
  scale_pattern_discrete("FDR Significance Level",choices = c( "stripe","none"), guide = "legend")+
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
  ylab("Spatially-weighted mean excursion frequency") +
  xlab("Age (BP)") +
  theme_bw()

#
hydroDirectional

ggsave(filename = paste0(figure_dir, "equalArea - hydroDirectional",weight.distance, ".pdf"),
       plot = hydroDirectional,
       width = 10, height = 3)
#
#
#make a plot. Any excursion, either direction

aboveColor <- RColorBrewer::brewer.pal(5,"RdBu")[1]
belowColor <- RColorBrewer::brewer.pal(5,"RdBu")[5]
lineColorAbove <- RColorBrewer::brewer.pal(7,"RdBu")[2]
lineColorBelow <- RColorBrewer::brewer.pal(7,"RdBu")[6]

library(ggpattern)
tempDirectional <- ggplot(tgdf) +
  geom_col_pattern(aes(x = age,y = allEventAbove,fill = "Positive",pattern = fdrSigLevelAbove), pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  geom_col_pattern(aes(x = age,y = -allEventBelow,fill = "Negative",pattern = fdrSigLevelBelow),pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "top") +
  geom_line(aes(x = age,y = -clBelow), color = lineColorBelow,linetype = 1) +
  geom_line(aes(x = age,y = clAbove), color = lineColorAbove,linetype = 1) +
  geom_label(aes(x = 12300,y = clAbove[58] ),label = "95% cl",color = lineColorAbove) +
  geom_label(aes(x = 12300,y = -clBelow[58]),label = "95% cl",color = lineColorBelow) +
  scale_pattern_discrete("FDR Significance Level",choices = c( "crosshatch", "stripe","none"), guide = "legend")+
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
  ylab("Spatially-weighted mean excursion frequency") +
  theme_bw()

tempDirectional

ggsave(filename = paste0(figure_dir, "equalArea - tempDirectional",weight.distance, ".pdf"),
       plot = tempDirectional,
       width = 10, height = 3)



