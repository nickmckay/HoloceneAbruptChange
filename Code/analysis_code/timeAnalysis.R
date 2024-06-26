# Set the directories for input data, output data, and figures
library(dggridR)
library(sf)
library(actR)
library(ggpattern)
library(tidyverse)
root_dir   <- 'analysis_code/'
output_dir <- file.path(root_dir,'new_output/')
figure_dir <- file.path(output_dir,'figures/')
map_dir <- file.path(figure_dir,"allMaps/")
data_dir <- file.path(output_dir,'data/')
processed_data_dir <- file.path(output_dir,'processed_data/')

source(file.path(root_dir,"helperFunctions.R"))


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

set.seed(82) #for reproducibility


# Create equal area grid to calculate weights -----------------------------
weight.distance <- 1500
dggs   <- dggridR::dgconstruct(spacing=weight.distance, metric=TRUE)

bigData$event.yr <- as.numeric(bigData$event.yr)
allTimes <- sort(unique(bigData$event.yr))

#use purrr and actr to calculate the fraction and significance of every timestep
globalTimeSeries <- purrr::map(allTimes,
                               calculateTimeSeries,
                               bigDataHydro = bigDataHydro,
                               bigDataTemp = bigDataTemp,
                               .progress = TRUE)

#pull out data
hgdf <- map(globalTimeSeries,pluck,1) |> list_rbind()
tgdf <- map(globalTimeSeries,pluck,2) |> list_rbind()
bgdf <- map(globalTimeSeries,pluck,3) |> list_rbind()

hgdf$age <- tgdf$age <- bgdf$age <-  allTimes

#calculate p-value after adjusting for False Discovery Rate (FDR)
effN <- nrow(hgdf)/2 #account for non-independence of tests

hgdf <- testFdr(hgdf,effN)
tgdf <- testFdr(tgdf,effN)
bgdf <- testFdr(bgdf,effN)



# Create plots ------------------------------------------------------------

patternSpacing <- 0.03
aboveColor <- RColorBrewer::brewer.pal(5,"BrBG")[5]
belowColor <- RColorBrewer::brewer.pal(5,"BrBG")[1]
lineColorPos <- RColorBrewer::brewer.pal(7,"BrBG")[6]
lineColorNeg <- RColorBrewer::brewer.pal(7,"BrBG")[2]

#make a plot. Any excursion, either direction
hydroDirectional <- ggplot(hgdf) +
  geom_col_pattern(aes(x = age,
                       y = allEventPos,
                       fill = "Positive",
                       pattern = fdrSigLevelPos),
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_colour  = 'gray30',
                   pattern_spacing = patternSpacing) +
  geom_col_pattern(aes(x = age,
                       y = -allEventNeg,
                       fill = "Negative",
                       pattern = fdrSigLevelNeg),
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_colour  = 'gray30',
                   pattern_spacing = patternSpacing) +
  geom_line(aes(x = age,y = -clNeg), color = lineColorNeg,linetype = 1) +
  geom_line(aes(x = age,y = clPos), color = lineColorPos,linetype = 1) +
  geom_label(aes(x = 12300,y = clPos[58] ),label = "95% cl",color = lineColorPos) +
  geom_label(aes(x = 12300,y = -clNeg[58]),label = "95% cl",color = lineColorNeg) +
  scale_pattern_discrete("FDR Significance Level",choices = c( "stripe","none"), guide = "legend")+
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "bottom") +
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

PosColor <- RColorBrewer::brewer.pal(5,"RdBu")[1]
belowColor <- RColorBrewer::brewer.pal(5,"RdBu")[5]
lineColorPos <- RColorBrewer::brewer.pal(7,"RdBu")[2]
lineColorNeg <- RColorBrewer::brewer.pal(7,"RdBu")[6]

library(ggpattern)
tempDirectional <- ggplot(tgdf) +
  geom_col_pattern(aes(x = age,y = allEventPos,fill = "Positive",pattern = fdrSigLevelPos), pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  geom_col_pattern(aes(x = age,y = -allEventNeg,fill = "Negative",pattern = fdrSigLevelNeg),pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "top") +
  geom_line(aes(x = age,y = -clNeg), color = lineColorNeg,linetype = 1) +
  geom_line(aes(x = age,y = clPos), color = lineColorPos,linetype = 1) +
  geom_label(aes(x = 12300,y = clPos[58] ),label = "95% cl",color = lineColorPos) +
  geom_label(aes(x = 12300,y = -clNeg[58]),label = "95% cl",color = lineColorNeg) +
  scale_pattern_discrete("FDR Significance Level",choices = c( "crosshatch", "stripe","none"), guide = "legend")+
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(PosColor,belowColor)) +
  ylab("Spatially-weighted mean excursion frequency") +
  theme_bw()

tempDirectional

ggsave(filename = paste0(figure_dir, "equalArea - tempDirectional",weight.distance, ".pdf"),
       plot = tempDirectional,
       width = 10, height = 3)

#plot both

PosColor <- RColorBrewer::brewer.pal(8,"PuOr")[7]
belowColor <- RColorBrewer::brewer.pal(8,"PuOr")[2]
lineColorPos <- RColorBrewer::brewer.pal(8,"PuOr")[6]
lineColorNeg <- RColorBrewer::brewer.pal(8,"PuOr")[3]

library(ggpattern)
bothDirectional <- ggplot(bgdf) +
  geom_col_pattern(aes(x = age,y = allEventPos,fill = "Positive",pattern = fdrSigLevelPos), pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  geom_col_pattern(aes(x = age,y = -allEventNeg,fill = "Negative",pattern = fdrSigLevelNeg),pattern_angle = 45,pattern_density = 0.01,pattern_colour  = 'gray30',pattern_spacing = patternSpacing) +
  scale_x_reverse("Age (BP)",breaks = seq(12000,000,by = -2000),position = "bottom") +
  geom_line(aes(x = age,y = -clNeg), color = lineColorNeg,linetype = 1) +
  geom_line(aes(x = age,y = clPos), color = lineColorPos,linetype = 1) +
  geom_label(aes(x = 12300,y = clPos[58] ),label = "95% cl",color = lineColorPos) +
  geom_label(aes(x = 12300,y = -clNeg[58]),label = "95% cl",color = lineColorNeg) +
  scale_pattern_discrete("FDR Significance Level",choices = c( "crosshatch", "stripe","none"), guide = "legend")+
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(PosColor,belowColor)) +
  scale_y_continuous("Spatially-weighted mean excursion frequency",breaks = seq(-.04,.04,by = .02)) +
  theme_bw()

bothDirectional

ggsave(filename = paste0(figure_dir, "equalArea - temp&hydroDirectional",weight.distance, ".pdf"),
       plot = bothDirectional,
       width = 10, height = 3)


