#data coverage figures
temp <- select(bigDataTemp,paleoData_TSid:paleoData_values)

temp <- temp[!duplicated(bigDataTemp$datasetId),]

tempTS <- purrr::transpose(temp)

uarchs <- unique(temp$archiveType)

#make a time availability plot
ta <- geoChronR::plotTimeAvailabilityTs(tempTS,age.range = c(0,12000), age.var = "time", group.var = "archiveType", step = 200)+
  scale_fill_brewer(palette = "Paired",limits = uarchs) +
  theme(legend.position = "none") +
  scale_x_reverse("Age (yr BP)",breaks = seq(12000,0,by = -2000))


#define the layout
lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2))

out <- gridExtra::arrangeGrob(grobs = list(tsMap,ta),
                               layout_matrix=lay)

ggsave(out,filename = "Temperature Excursion Data Density.pdf",width = 8,height = 5)




# Now hydroclimate --------------------------------------------------------


hydro <- select(bigDataHydro,paleoData_TSid:paleoData_values) %>%
  distinct()

hydroTS <- purrr::transpose(hydro)

#remove duplicates
dsndup <- which(duplicated(hydro$datasetId))
hydroTS <- hydroTS[-dsndup]

uarchs <- c(unique(temp$archiveType),"Shoreline")

#make a map
tsMap <- geoChronR::mapTs(hydroTS,size = 3) + scale_color_brewer(palette = "Paired",limits = uarchs) +
  ggtitle("Hydroclimate datasets")

#make a time availability plot
taHydro <- geoChronR::plotTimeAvailabilityTs(hydroTS,age.range = c(0,12000), age.var = "time", group.var = "archiveType", step = 200)+
  scale_fill_brewer(palette = "Paired",limits = uarchs) +
  theme(legend.position = "none") +
  scale_x_reverse("Age (yr BP)",breaks = seq(12000,0,by = -2000))


#define the layout
lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(2,2,2),
             c(2,2,2))

out <- gridExtra::arrangeGrob(grobs = list(tsMap,taHydro),
                               layout_matrix=lay)

ggsave(out,filename = "Hydroclimate Excursion Data Density.pdf",width = 8,height = 5)
