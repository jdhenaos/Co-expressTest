library("GEOquery")

GEO <- function(x,y="."){
  if(dir.exists(y)){
    options(warn = -1)
  }else{
    dir.create(y)
  }
  for(i in x[,1]){
    getGEO(GEO = i, destdir = y)
  }
}

GeneSymbol <- function(GPL){
  
}

#ADGEO <- read.table("Alzheimer_Chips.txt")
#PDGEO <- read.table("Parkinson_Chips.txt")
#MSGEO <- read.table("MultipleSclerosis_Chips.txt")

#GEO(ADGEO,"./Alzheimer_GSE")
#GEO(PDGEO,"./Parkinson_GSE")
#GEO(MSGEO,"./MultipleSclerosis_GSE")

GEO(read.table("geos.txt"),"./Alz")
