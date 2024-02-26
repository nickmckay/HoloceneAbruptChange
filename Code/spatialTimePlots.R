#load and plot spatial outputs
library(actR)
library(purrr)
library(tidyverse)
library(lipdR)
source("spatial.R")
# Set the directories for input data, output data, and figures
data_dir   = '/HoloceneAbruptChange/actr_results_all/'
figure_dir = '/HoloceneAbruptChange/spatial_figs/'
processed_data_dir = '/HoloceneAbruptChange/processed_data/'


# First load in the RDS ---------------------------------------------------
allRDS <- list.files(figure_dir,pattern = ".RDS")

allTemp <- allRDS[str_detect(allRDS,"Temp")]
allHydro <- allRDS[str_detect(allRDS,"Hydro")]

#load in all Temp
for(a in allTemp){
  tots.temp <- readRDS(file.path(figure_dir,a))


  if(a == allTemp[1]){
    allTot.temp <- tots.temp
  }else{
    allTot.temp <- bind_rows(allTot.temp,tots.temp)
  }
}

#load in all hydro
for(a in allHydro){
  tots.hydro <- readRDS(file.path(figure_dir,a))


  if(a == allHydro[1]){
    allTot.hydro  <- tots.hydro
  }else{
    allTot.hydro  <- bind_rows(allTot.hydro,tots.hydro)
  }
}

aboveColor <- RColorBrewer::brewer.pal(5,"RdBu")[1]
belowColor <- RColorBrewer::brewer.pal(5,"RdBu")[5]
lineColorAbove <- RColorBrewer::brewer.pal(7,"RdBu")[2]
lineColorBelow <- RColorBrewer::brewer.pal(7,"RdBu")[6]

#old way
# #make a plot. Any excursion, either direction
# tempDirectionalSpatial <-ggplot(allTot.temp) +
#   geom_bar(aes(x = age,y = weightedAbove,fill = "Positive"),stat = "identity") +
#   geom_bar(aes(x = age,y = -weightedBelow,fill = "Negative"),stat = "identity") +
#   #geom_area(aes(x = age,y = weightedNet, fill = "Net"),color = "black") +
#   scale_x_reverse(breaks = seq(12000,2000,by = -2000)) +
#   geom_line(aes(x = age,y = `above_95.`), color = lineColorAbove) +
#   geom_line(aes(x = age,y = -`below_95.`), color = lineColorBelow) +
#   geom_label(aes(x = 1200,y = above_95.[1] + .005),label = "95% cl",color = lineColorAbove) +
#   geom_label(aes(x = 1200,y = -(below_95.[1] + .005)),label = "95% cl",color = lineColorBelow) +
#   scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
#   ylab("Fraction of events detected") +
#   theme_bw() +
#   xlab("Age (yr BP)")

#make a plot. Any excursion, either direction
tempDirectionalSpatial <-ggplot(allTot.temp) +
  geom_bar(aes(x = age,y = weightedPassAbove,fill = "Positive"),stat = "identity") +
  geom_bar(aes(x = age,y = -weightedPassBelow,fill = "Negative"),stat = "identity") +
  #geom_area(aes(x = age,y = weightedNet, fill = "Net"),color = "black") +
  scale_x_reverse(breaks = seq(12000,2000,by = -2000)) +
  geom_line(aes(x = age,y = 0.05), color = lineColorAbove) +
  geom_line(aes(x = age,y = -0.05), color = lineColorBelow) +
  geom_label(aes(x = 1200,y = .055),label = "95% cl",color = lineColorAbove) +
  geom_label(aes(x = 1200,y = -.055),label = "95% cl",color = lineColorBelow) +
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
  ylab("Fraction of events detected") +
  theme_bw() +
  xlab("Age (yr BP)")


tempDirectionalSpatial

ggsave(filename = paste0(figure_dir, "tempDirectional", ".pdf"),
       plot = tempDirectionalSpatial,
       width = 8, height = 6)


#make a plot. Any excursion, either direction
aboveColor <- RColorBrewer::brewer.pal(5,"BrBG")[5]
belowColor <- RColorBrewer::brewer.pal(5,"BrBG")[1]
lineColorAbove <- RColorBrewer::brewer.pal(7,"BrBG")[6]
lineColorBelow <- RColorBrewer::brewer.pal(7,"BrBG")[2]

# hydroDirectionalSpatial <- ggplot(allTot.hydro) +
#   geom_bar(aes(x = age,y = weightedAbove,fill = "Positive"),stat = "identity") +
#   geom_bar(aes(x = age,y = -weightedBelow,fill = "Negative"),stat = "identity") +
#   #geom_area(aes(x = age,y = weightedNet, fill = "Net"),color = "black") +
#   scale_x_reverse(breaks = seq(12000,2000,by = -2000)) +
#   geom_line(aes(x = age,y = `above_95.`), color = lineColorAbove) +
#   geom_line(aes(x = age,y = -`below_95.`), color = lineColorBelow) +
#   geom_label(aes(x = 1200,y = above_95.[1] + .005),label = "95% cl",color = lineColorAbove) +
#   geom_label(aes(x = 1200,y = -(below_95.[1] + .005)),label = "95% cl",color = lineColorBelow) +
#   scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
#   ylab("Fraction of events detected") +
#   theme_bw() +
#   xlab("Age (yr BP)")


hydroDirectionalSpatial <- ggplot(allTot.hydro) +
  geom_bar(aes(x = age,y = weightedPassAbove,fill = "Positive"),stat = "identity") +
  geom_bar(aes(x = age,y = -weightedPassBelow,fill = "Negative"),stat = "identity") +
  #geom_area(aes(x = age,y = weightedNet, fill = "Net"),color = "black") +
  scale_x_reverse(breaks = seq(12000,2000,by = -2000)) +
  geom_line(aes(x = age,y = 0.05), color = lineColorAbove) +
  geom_line(aes(x = age,y = -0.05), color = lineColorBelow) +
  geom_label(aes(x = 1200,y = .055),label = "95% cl",color = lineColorAbove) +
  geom_label(aes(x = 1200,y = -.055),label = "95% cl",color = lineColorBelow) +
  scale_fill_manual("Excursion Direction",limits = c("Positive","Negative"),values = c(aboveColor,belowColor)) +
  ylab("Fraction of events detected") +
  theme_bw() +
  xlab("Age (yr BP)")


hydroDirectionalSpatial

ggsave(filename = paste0(figure_dir, "hydroDirectional", ".pdf"),
       plot = hydroDirectionalSpatial,
       width = 8, height = 6)



