#RUn this from insied the equalAreaWeights function

#Get the grid cell boundaries for cells which had quakes
grid          <- dgcellstogrid(dggs,dfb$grid)

#Update the grid cells' properties to include the number of earthquakes
#in each cell
grid          <- merge(grid,dfb,by.x="seqnum",by.y="grid")

#Make adjustments so the output is more visually interesting
grid$count    <- log(grid$count)
cutoff        <- quantile(grid$count,0.9)
#grid          <- grid %>% mutate(count=ifelse(count>cutoff,cutoff,count))

#Get polygons for each country of the world
cut <- 1.5
world <- map_data("world") %>%
  filter(abs(long) < 180-cut,
         abs(lat) < 90-cut)

wrapped_grid = st_wrap_dateline(grid, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

gridMapOut <- ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group), fill=NA, color="black")   +
  geom_sf     (data=wrapped_grid, aes(fill=weight), color=alpha("white", 0.4),alpha = 0.7) +
  geom_point  (data = df,aes(x=lon, y=lat),color = "black",size = 1) +
  scale_fill_viridis_c(direction = -1) +
  theme(legend.title = element_blank())+
  cowplot::theme_minimal_grid()+
  coord_sf(xlim=c(-179,179),
           ylim = c(-89,89),
           crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
           expand = TRUE,
           datum = sf::st_crs(4326),
           default_crs = sf::st_crs(4326)) +
  xlab("") + ylab("")

ggsave("~/Dropbox/HoloceneAbruptChange/EqualAreaGridWeights.pdf",plot = gridMapOut)
