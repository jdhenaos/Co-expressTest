library("GEOquery")

GEO <- function(x,y){
  dir.create(y)
  for(i in test[,1]){
    getGEO(GEO = i, destdir = y)
  }
}

test <- read.table("geos.txt")
GEO(test,"./Alz")
