library("GEOquery")

GEO <- function(x){
  getGEO(GEO = x, destdir = ".")
}

test <- read.table("geos.txt")

for(i in test[,1]){
  GEO(i)
}
