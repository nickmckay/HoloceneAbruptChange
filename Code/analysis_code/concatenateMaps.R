library(magick)
library(pdftools)

root_dir   <- 'analysis_code/'
output_dir <- file.path(root_dir,'new_output/')
figure_dir <- file.path(output_dir,'figures/')
map_dir <- file.path(figure_dir,"allMaps/")
data_dir <- file.path(output_dir,'data/')
processed_data_dir <- file.path(output_dir,'processed_data/')


am <- list.files(map_dir,pattern = ".pdf",full.names = TRUE)
#reorder in time

allNumbers <- as.numeric(str_extract(am,"\\d+"))
so <- sort(allNumbers,decreasing = TRUE,index.return = TRUE)

am2 <- am[so$ix]
for(r in 1:29){

  rowvec <- c(1:4) + 4*(r-1)

  for(i in rowvec){
    thisImage <- am2[i]

    year <- str_c(as.numeric(str_extract(thisImage,"\\d+"))/1000," ka                                        ")
    isTemp <- str_detect(thisImage,pattern = "temp")

    img <- image_read_pdf(path = thisImage)

    scaled <- image_scale(img,"1200")

    cropped <- image_crop(scaled, "1010x600+70+20")

    if(isTemp){
      cropped <- image_annotate(image = cropped,
                                text = year,
                                size = 20,
                                color = "white",
                                boxcolor = "white",
                                location = "+00+00")
    }else{
      cropped <- image_annotate(image = cropped,
                                text = year,
                                size = 20,
                                color = "black",
                                boxcolor = "white",
                                location = "+00+00")
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

image_write(all,path = file.path(figure_dir,"combinedLong.png"),format = "png")
