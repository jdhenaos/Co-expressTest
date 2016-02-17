library("GEOquery")

GEO <- function(x,y){
  if(dir.exists(y)){
    options(wwarn = -1)
  }else{
    dir.create(y)
  }
  for(i in test[,1]){
    getGEO(GEO = i, destdir = y)
  }
}

test <- read.table("geos.txt")
GEO(test,"./Alz")