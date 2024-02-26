equalAreaGridWeights <- function(lon,lat,spacing,dggs){
  if(!exists("dggs")){
    dggs   <- dggridR::dgconstruct(spacing=spacing, metric=TRUE)
  }

  df <- data.frame(lon = lon,
                   lat = lat,
                   grid = dggridR::dgGEO_to_SEQNUM(dggs,lon,lat)$seqnum)

  count <- df %>% group_by(grid) %>% summarise(count=n())

  dfb <- dplyr::left_join(df,count, by = "grid") %>%
    mutate(weight = 1/count)

  return(dfb)

}

