library(tidyverse)

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

#data coverage figures
temp <- select(bigData,
               paleoData_TSid,
               dataSetName,
               geo_latitude,
               geo_longitude,
               archiveType) %>%
  distinct()


temp$isHydro <- (temp$paleoData_TSid %in% bigDataHydro$paleoData_TSid)
temp$isTemp <- (temp$paleoData_TSid %in% bigDataTemp$paleoData_TSid)

toPlot <- temp %>%
  group_by(dataSetName) %>%
  summarize(lat = mean(geo_latitude),
            lon = mean(geo_longitude),
            archiveType = unique(archiveType),
            variable = case_when(any(isHydro) & any(isTemp) ~ "both",
                                 any(isHydro) & !any(isTemp) ~ "hydroclimate",
                                 !any(isHydro) & any(isTemp) ~ "temperature",
                                 !any(isHydro) & !any(isTemp) ~ "none")) %>%
  filter(variable != "none")



uarchs <- unique(temp$archiveType)

#make a map
cut <- 1.5

world <- map_data("world") %>%
  filter(abs(long) < 180-cut,
         abs(lat) < 90-cut)
projection = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

map <- ggplot() +
  geom_polygon(data=world, aes(x = long, y = lat, group = group), color="gray50", fill=alpha("white",0.1), inherit.aes = FALSE) +
  #geom_point(inherit.aes = FALSE, data=toPlot, aes(x=lon, y=lat ,shape = variable), color="black", size = 3)+
  geom_point(inherit.aes = FALSE, data=toPlot, aes(x=lon, y=lat, fill=archiveType, shape = variable), color = "black",size = 2)+
  scale_x_continuous("",breaks = seq(-180,180,by = 60)) +
  scale_y_continuous("",breaks = c(-87.5,seq(-60,60,by = 30),87.5)) +
  scale_fill_brewer(palette = "Paired",limits = uarchs) +
  scale_shape_manual(values = c(24,22,21)) +
  #geom_text(aes(x = -200,y = 0),label = "a)") +
  guides(fill = guide_legend(override.aes=list(shape=21),direction = "horizontal",
                             order = 1),
         shape = guide_legend(direction = "vertical"),)+
  cowplot::theme_minimal_grid()+
  theme(legend.title = element_blank(),
        legend.key.spacing = unit(0.1, "cm"),
        legend.position = "top")+
  coord_sf(xlim=c(-180,180),
           ylim = c(-90,90),
           crs = projection,
           expand = TRUE,
           datum = sf::st_crs(4326),
           default_crs = sf::st_crs(4326))

map

#data coverage figures
temp <- select(bigDataTemp,paleoData_TSid:paleoData_values)

temp <- temp[!duplicated(bigDataTemp$datasetId),]

tempTS <- purrr::transpose(temp)

#make a time availability plot
ta <- geoChronR::plotTimeAvailabilityTs(tempTS,age.range = c(0,12000), age.var = "time", group.var = "archiveType", step = 200)+
  scale_fill_brewer(palette = "Paired",limits = uarchs) +
  theme(legend.position = "none") +
  scale_x_reverse("Age (yr BP)",breaks = seq(12000,0,by = -2000),position = "bottom",expand = c(0.01,0.01)) +
  scale_y_continuous("Count",position = "left",expand = c(0.01,0.01),limits = c(0,800))+
  #geom_text(aes(x = 11000,y = 800),label = "b)") +
  ggtitle("Temperature") +
  theme_bw() +
  theme(legend.position = "none",plot.margin=unit(c(-1,0.25,0.25,0.25), "cm"))

#hydro
hydro <- select(bigDataHydro,paleoData_TSid:paleoData_values)

hydro <- hydro[!duplicated(bigDataHydro$datasetId),]

hydroTS <- purrr::transpose(hydro)

#make a time availability plot
ha <- geoChronR::plotTimeAvailabilityTs(hydroTS,age.range = c(0,12000), age.var = "time", group.var = "archiveType", step = 200)+
  scale_fill_brewer(palette = "Paired",limits = uarchs) +
  theme(legend.position = "none") +
  scale_x_reverse("Age (yr BP)",breaks = seq(12000,0,by = -2000),expand = c(0.01,0.01)) +
  scale_y_continuous("Count",limits = c(0,800),position = "right",expand = c(0.01,0.01))+
  ggtitle("Hydroclimate") +
  #geom_text(aes(x = 11000,y = 800),label = "c)") +
  theme_bw() +
  theme(legend.position = "none",plot.margin=unit(c(-1,0.25,0.25,0.25), "cm"))


#define the layout
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,3,3))

out <- gridExtra::arrangeGrob(grobs = list(map,ta,ha),layout_matrix = lay)

ggsave(out,filename = file.path(figure_dir,"ExcursionDataAbruptChange_spacetime.pdf"),width = 8.5,height = 7)


