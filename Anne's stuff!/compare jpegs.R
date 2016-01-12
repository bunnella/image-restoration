collapse_image <- function(big_array) {
  C <- ncol(big_array)
  R <- nrow(big_array)
  t(big_array[R:1, 1:C, 1])
}

library(bmp)
library(jpeg)

setwd("C:/Users/Anne/Downloads/image-restoration/Anne's stuff!/img")

image_folder <- getwd()

list.files(path = image_folder, pattern = ".*bmp")

filename <- "house.bmp"

for (filename in list.files(path = image_folder, pattern = ".*bmp")) {
  original <- read.bmp(filename)
  
  original <- original/255
  
  writeJPEG(original, target = "test.jpg", quality = 0.3)
  jpeg <- readJPEG("test.jpg")
  
  orig.c <- collapse_image(original)
  jpeg.c <- collapse_image(jpeg)
  
  diff <- orig.c-jpeg.c
  
  #image(orig.c,asp=1, xaxt="n", yaxt="n")
  #image(jpeg.c,asp=1, xaxt="n", yaxt="n")
  jpeg(paste("diffs/", filename, ".jpg", sep=""))
  image(orig.c,asp=1, frame=F, xaxt="n", yaxt="n", main = filename)
  dev.off()
}
