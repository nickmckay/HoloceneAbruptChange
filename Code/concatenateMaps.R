library(magick)

am <- list.files("/HoloceneAbruptChange/spatial_figs_final/allMaps/",full.names = TRUE)
#reorder in time

allNumbers <- as.numeric(str_extract(am,"\\d+"))
so <- sort(allNumbers,decreasing = TRUE,index.return = TRUE)

am2 <- am[so$ix]
for(r in 1:6){
  rowvec <- case_when(
    r == 1 ~ (1:18)*2, # temp 1
    r == 3 ~ (1:18)*2 + 36, # temp 2
    r == 5 ~ (1:18)*2 + 72, # temp 3
    r == 2 ~ (1:18)*2 - 1, # hydro 1
    r == 4 ~ (1:18)*2 + 35, # hydro 2
    r == 6 ~ (1:18)*2 + 71 # hydro 3    
    
  )
  
  for(i in rowvec){
    thisImage <- am2[i]
    
    year <- str_c(as.numeric(str_extract(thisImage,"\\d+"))/1000," ka                                        ")
    isTemp <- str_detect(thisImage,pattern = "Temperature")
    
    img <- image_read_pdf(path = thisImage)
    
    scaled <- image_scale(img,"1200")
    
    if(isTemp){
      cropped <- image_crop(scaled, "1010x600+70+20")
      cropped <- image_annotate(image = cropped,year, size = 20, color = "black",boxcolor = "white",location = "+00+00")
      cropped
    }else{
      cropped <- image_crop(scaled, "1010x550+70+100")
    }
    
    if(i == min(rowvec)){
      combined <- cropped
    }else{
      combined <- c(combined,cropped)
    }
  }
  out <- image_append(combined,stack = FALSE)
  if(r == 1){
    stacked <- out
  }else{
    stacked <- c(stacked,out)
  }
  print(r)
}
all <- image_append(stacked,stack = TRUE)

image_write(all,path = "/HoloceneAbruptChange/spatial_figs_final/combined.png",format = "png")
