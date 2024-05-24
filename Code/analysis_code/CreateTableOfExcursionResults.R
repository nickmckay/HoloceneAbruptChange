#Create table of excursion timing
root_dir   <- 'analysis_code/'
output_dir <- file.path(root_dir,'new_output/')
figure_dir <- file.path(output_dir,'figures/')
map_dir <- file.path(figure_dir,"allMaps/")
data_dir <- file.path(output_dir,'data/')
processed_data_dir <- file.path(output_dir,'processed_data/')

bigDataHydro$midPointTime <- round(rowMeans(select(bigDataHydro,time_start,time_end)),-2)
bigDataTemp$midPointTime <- round(rowMeans(select(bigDataTemp,time_start,time_end)),-2)


getSigTimes <- function(time_mid,pval,alpha = 0.05){
  good <- which(pval < alpha)
  outStr <- paste(unique(sort(time_mid[good])),collapse = ", ")
  return(outStr)
}

excursionsByDatasetHydro <- bigDataHydro %>%
  group_by(dataSetName) %>%
  summarize(significantWet = getSigTimes(midPointTime,
                                              pvalue_positive),
            significantDry = getSigTimes(midPointTime,
                                          pvalue_negative),
            lat = unique(geo_latitude))

excursionsByDatasetTemp <- bigDataTemp %>%
  group_by(dataSetName) %>%
  summarize(significantWarm = getSigTimes(midPointTime,
                                         pvalue_positive),
            significantCold = getSigTimes(midPointTime,
                                         pvalue_negative),
            lat = unique(geo_latitude))

excursionsByDataset <- full_join(excursionsByDatasetTemp,excursionsByDatasetHydro,by = c("dataSetName","lat"))

sum(excursionsByDataset$lat >= 30)/nrow(excursionsByDataset)

write_csv(select(excursionsByDataset,-lat),file.path(processed_data_dir,"ExcursionsByDataset.csv"))

# Fraction significant by time --------------------------------------------


fSignificant <- function(dataSetName,pval,alpha = 0.05){
  good <- which(pval < alpha)
  outFrac <- length(unique(dataSetName[good])) / length(unique(dataSetName))
  return(outFrac)
}

fSignificantEither <- function(dataSetName,pval1,pval2,alpha = 0.05){
  good <- which((pval1 < alpha) | (pval2 < alpha))
  outFrac <- length(unique(dataSetName[good])) / length(unique(dataSetName))
  return(outFrac)
}

excursionsByYearHydro <- bigDataHydro %>%
  group_by(midPointTime) %>%
  summarize(fractionSigWet = fSignificant(dataSetName,
                                           pvalue_positive),
            fractionSigDry = fSignificant(dataSetName,
                                          pvalue_negative),
            fractionSigHydro = fSignificantEither(dataSetName,
                                          pvalue_positive,
                                          pvalue_negative)
            )

excursionsByYearTemp <- bigDataTemp %>%
  group_by(midPointTime) %>%
  summarize(fractionSigWarm = fSignificant(dataSetName,
                                           pvalue_positive),
            fractionSigCold = fSignificant(dataSetName,
                                           pvalue_negative),
            fractionSigTemp = fSignificantEither(dataSetName,
                                                 pvalue_positive,
                                                 pvalue_negative)
  )


aboveColorHydro <- RColorBrewer::brewer.pal(5,"BrBG")[5]
belowColorHydro <- RColorBrewer::brewer.pal(5,"BrBG")[1]

aboveColorTemp <- RColorBrewer::brewer.pal(5,"RdBu")[1]
belowColorTemp <- RColorBrewer::brewer.pal(5,"RdBu")[5]



tempExcFraction <- ggplot(excursionsByYearTemp) +
  geom_col(aes(x = midPointTime, y = fractionSigWarm),fill = aboveColorTemp) +
  geom_col(aes(x = midPointTime, y = -fractionSigCold),fill = belowColorTemp) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "top") +
  theme_bw()



hydroExcFraction <- ggplot(excursionsByYearHydro) +
  geom_col(aes(x = midPointTime, y = fractionSigWet),fill =aboveColorHydro) +
  geom_col(aes(x = midPointTime, y = -fractionSigDry),fill = belowColorHydro) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "bottom") +
  theme_bw()

fracBoth <- egg::ggarrange(plots = list(tempExcFraction,hydroExcFraction),nrow = 2)

ggsave(filename = file.path(figure_dir,"Supplemental Figure - fraction of significant excursions by time.pdf"),
       plot = fracBoth,
       width = 10, height = 6)


apply(excursionsByYearTemp[,-1],2,mean)
apply(excursionsByYearTemp[,-1],2,range)

apply(excursionsByYearHydro[,-1],2,mean)
apply(excursionsByYearHydro[,-1],2,range)


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
mean(rowSums(excursionsByDataset[,-1] == "" | is.na(excursionsByDataset[,-1])) < 4)


# Create S2 table ---------------------------------------------------------


DatasetMetadata <- read_csv("https://lipdverse.org/4_2_refs/0_13_0/references.csv")
ExcursionMetadata <- read_csv(file.path(processed_data_dir,"ExcursionsByDataset.csv"))

  left_join(ExcursionMetadata,DatasetMetadata,by = "dataSetName") %>%
  arrange(desc(Lat)) %>%
  select(!starts_with("significant"),starts_with("significant")) %>%
  write_csv(file.path(processed_data_dir,"DatasetS1.csv"))




