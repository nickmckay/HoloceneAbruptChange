all.age.to.plot <- seq(1400,11600,by = 200)
#all.age.to.plot <- c(9200,9400,8200)
weight.dist <- 1500
future::plan(future::multisession,workers = 6) #this will run in parallel
options(future.rng.onMisuse = "ignore")


n.colors <- 11
grayCol <-  "#CCCCCC"
dryCol <- RColorBrewer::brewer.pal(n.colors,"BrBG")
dryCol[median(seq_along(dryCol))] <- grayCol
dryCol <- dryCol[rev(seq_along(dryCol))]

hotCol <- RColorBrewer::brewer.pal(n.colors,"RdBu")
hotCol[median(seq_along(hotCol))] <- grayCol

for(age.to.plot in all.age.to.plot){

#pvalOptions <- c("pvalAbove","pvalBelow","pvalNet")
pvalOptions <- c("pvalNet")
  for(po in pvalOptions){
    #now with weights
    figoutTemp <- plotSignificance(pval.grid.weight.temp,
                                   sigTestResults = filter(bigDataTemp,event.yr == age.to.plot),
                                   which.test = po,
                                   color.breaks = c(0.01,0.05,0.1,0.2),
                                   restrict.sites = TRUE,
                                   alpha.by.weight = TRUE,
                                   cutoff.distance = weight.dist)

    if(po == "pvalNet"){
      fwd <- \(x) (x-.5)^3
      inv <-  \(x) ifelse(x>=0,yes = x^(1/3),no = -1*abs(x)^(1/3)) +.5
      cube <- scales::trans_new("cube",transform = fwd, inverse = inv)



      toplot <- figoutTemp$plot +
       # scale_fill_fermenter(palette = "BrBG",direction = -1,breaks = c(.01,.05,0.1,.2,.8,.9,.95,.99))
      #scale_fill_fermenter(palette = "RdBu",direction = 1,breaks = c(.01,.05,0.1,.2,.8,.9,.95,.99))
      # scale_fill_paletteer_binned("scico::broc",direction = -1,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube)
      #scale_fill_fermenter(palette = "BrBG",direction = -1,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube)
      scale_fill_stepsn(colours = hotCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube)

      #scale_fill_paletteer_c("scico::vik",direction = 1,trans = cube)

    }else{
      toplot <- figoutTemp$plot
    }

    plot.name <-glue::glue("{age.to.plot} - Temperature - {po} - Robust Null - {weight.dist} km")

    toplot <- toplot + ggtitle(plot.name)


    ggsave(
      filename = paste0(figure_dir, glue::glue("{plot.name}.pdf")),
      plot = toplot,
      width = 15, height = 9
    )


    figoutHydro <- plotSignificance(pval.grid.weight.hydro,
                                   sigTestResults = filter(bigDataHydro,event.yr == age.to.plot),
                                   which.test = po,
                                   color.breaks = c(0.01,0.05,0.1,0.2),restrict.sites = TRUE,
                                   alpha.by.weight = TRUE,cutoff.distance = weight.dist)

    if(po == "pvalNet"){
      toplot <- figoutHydro$plot +
        scale_fill_stepsn(colours = dryCol,breaks = c(.05,0.1,.2,.8,.9,.95),trans = cube)

        #scale_fill_fermenter(palette = "BrBG",direction = -1,breaks = c(.01,.05,0.1,.2,.8,.9,.95,.99))
        #scale_fill_fermenter(palette = "RdBu",direction = 1,breaks = c(.01,.05,0.1,.2,.8,.9,.95,.99))
    }else{
      toplot <- figout$plot
    }

    plot.name <-glue::glue("{age.to.plot} - Hydroclimate - {po} - Robust Null - {weight.dist} km")

    toplot <- toplot + ggtitle(plot.name)


    ggsave(
      filename = paste0(figure_dir, glue::glue("{plot.name}.pdf")),
      plot = toplot,
      width = 15, height = 9
    )



  }
}
