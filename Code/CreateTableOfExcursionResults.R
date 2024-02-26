#Create table of excursion timing

newBigData <- bind_rows(bigDataHydro,bigDataTemp)

newBigData$midPointTime <- round(rowMeans(select(newBigData,time_start,time_end)),-2)

getSigTimes <- function(time_mid,pval,interpretation1_variable,interpOpts,alpha = 0.05){
  good <- which(pval < alpha & interpretation1_variable %in% interpOpts)
  outStr <- paste(unique(sort(time_mid[good])),collapse = ", ")
  return(outStr)
}

excursionsByDataset <- newBigData %>%
  group_by(dataSetName) %>%
  summarize(significantWarm = getSigTimes(midPointTime,
                                              pvalue_above,
                                              interpretation1_variable,
                                              interpOpts = tempVars),
            significantCold = getSigTimes(midPointTime,
                                              pvalue_below,
                                              interpretation1_variable,
                                              interpOpts = tempVars),
            significantWet = getSigTimes(midPointTime,
                                              pvalue_above,
                                              interpretation1_variable,
                                              interpOpts = hydroVars),
            significantDry = getSigTimes(midPointTime,
                                          pvalue_below,
                                              interpretation1_variable,
                                              interpOpts = hydroVars))




write_csv(excursionsByDataset,"ExcursionsByDataset.csv")

# Fraction significant by time --------------------------------------------


fSignificant <- function(dataSetName,pval,interpretation1_variable,interpOpts,alpha = 0.05){
  good <- which(pval < alpha & interpretation1_variable %in% interpOpts)
  total <- which(interpretation1_variable %in% interpOpts)
  outFrac <- length(unique(dataSetName[good])) / length(unique(dataSetName[total]))
  return(outFrac)
}

fSignificantEither <- function(dataSetName,pval1,pval2,interpretation1_variable,interpOpts,alpha = 0.05){
  good <- which((pval1 < alpha & interpretation1_variable %in% interpOpts) | (pval2 < alpha & interpretation1_variable %in% interpOpts))
  total <- which(interpretation1_variable %in% interpOpts)
  outFrac <- length(unique(dataSetName[good])) / length(unique(dataSetName[total]))
  return(outFrac)
}

excursionsByYear <- newBigData %>%
  group_by(midPointTime) %>%
  summarize(fractionSigWarm = fSignificant(dataSetName,
                                          pvalue_above,
                                          interpretation1_variable,
                                          interpOpts = tempVars),
            fractionSigCold = fSignificant(dataSetName,
                                           pvalue_below,
                                           interpretation1_variable,
                                           interpOpts = tempVars),
            fractionSigWet = fSignificant(dataSetName,
                                           pvalue_above,
                                           interpretation1_variable,
                                           interpOpts = hydroVars),
            fractionSigDry = fSignificant(dataSetName,
                                          pvalue_below,
                                           interpretation1_variable,
                                           interpOpts = hydroVars),
            fractionSigHydro = fSignificantEither(dataSetName,
                                          pvalue_above,
                                          pvalue_below,
                                          interpretation1_variable,
                                          interpOpts = hydroVars),
            fractionSigTemp = fSignificantEither(dataSetName,
                                            pvalue_above,
                                            pvalue_below,
                                            interpretation1_variable,
                                            interpOpts = tempVars)
            )



aboveColorHydro <- RColorBrewer::brewer.pal(5,"BrBG")[5]
belowColorHydro <- RColorBrewer::brewer.pal(5,"BrBG")[1]

aboveColorTemp <- RColorBrewer::brewer.pal(5,"RdBu")[1]
belowColorTemp <- RColorBrewer::brewer.pal(5,"RdBu")[5]



tempExcFraction <- ggplot(excursionsByYear) +
  geom_col(aes(x = midPointTime, y = fractionSigWarm),fill = aboveColorTemp) +
  geom_col(aes(x = midPointTime, y = -fractionSigCold),fill = belowColorTemp) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "top") +
  theme_bw()



hydroExcFraction <- ggplot(excursionsByYear) +
  geom_col(aes(x = midPointTime, y = fractionSigWet),fill =aboveColorHydro) +
  geom_col(aes(x = midPointTime, y = -fractionSigDry),fill = belowColorHydro) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "bottom") +
  theme_bw()

fracBoth <- egg::ggarrange(plots = list(tempExcFraction,hydroExcFraction),nrow = 2)

ggsave(filename = "Supplemental Figure - fraction of significant excursions by time.pdf",
       plot = fracBoth,
       width = 10, height = 6)


apply(excursionsByYear[,-1],2,mean)
apply(excursionsByYear[,-1],2,range)



# get data stats ----------------------------------------------------------
totalDatasetsAnalyzed <- nrow(excursionsByDataset)
hydroDatasets <- unique(bigDataHydro$dataSetName)
nHydroDatasets <- length(hydroDatasets)

tempDatasets <- unique(bigDataTemp$dataSetName)
nTempDatasets <- length(tempDatasets)

bothDatasets <- length(intersect(tempDatasets,hydroDatasets))

seasonTable <- newBigData %>%
  group_by(dataSetName) %>%
  summarize(nSeasonalities = length(unique(interpretation1_seasonality)))

sum(seasonTable$nSeasonalities > 1)

#percent of datasets with at least one excursion
mean(rowSums(excursionsByDataset[,-1] == "") < 4)


# Create S2 table ---------------------------------------------------------


DatasetMetadata <- read_csv("/lipdverse/html/4_2_refs/0_13_0/references.csv")
ExcursionMetadata <- read_csv("/HoloceneAbruptChange/ExcursionsByDataset.csv")

DatasetS1 <- left_join(ExcursionMetadata,DatasetMetadata,by = "dataSetName") %>%
  arrange(desc(Lat)) %>%
  select(everything(), starts_with("significance")) %>%
  write_csv("/HoloceneAbruptChange/DatasetS1.csv")




