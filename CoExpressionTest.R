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

GeneSymbol <- function(GPL, d = "."){
  setwd(d)
  gpl <- getGEO(filename = dir(".")[grep(GPL,dir("."))])
  write.table(table(Table(gpl)$"Gene Symbol"), file = "output")
  setwd("../")
}

ExtractInfo <- function(x,d = "."){
  setwd(d)
  f<-read.csv(gzfile(x),
                     comment.char = "!", 
                     sep = "\t",
                     stringsAsFactors = FALSE)
  setwd("../")
  return(f)
}

DataUnion <- function(d = "."){
  f <- dir(d)[grep("^GSE[0-9]+(_|-GPL570)",dir(d))]
  g <- 0
  for (t in f){
    if(g == 0){
      g <- ExtractInfo(t,d)
      options(warn = -1)
    }else{
      h <- ExtractInfo(t,d)
      g <- merge(g,h)
    }
  }
  return(g)
}

#GEO(read.table("Alzheimer_Chips.txt"),"./Alzheimer_GSE")
#GEO(read.table("Parkinson_Chips.txt"),"./Parkinson_GSE")
#GEO(read.table("MultipleSclerosis_Chips.txt"),"./MultipleSclerosis_GSE")

GEO(read.table("geos.txt"),"./Alz")
GeneSymbol("GPL570","./Alz")
D <- DataUnion("./Alz")
dim(D)

#######################

GSE68527 <- ExtractInfo("GSE68527_series_matrix.txt.gz")
GSE52139 <- ExtractInfo("GSE52139_series_matrix.txt.gz")
GSE16759 <- ExtractInfo("GSE16759-GPL570_series_matrix.txt.gz")


dim(GSE68527)
dim(GSE52139)
dim(GSE16759)

f <- dir("./Alz")[grep("^GSE[0-9]+(_|-GPL570)",dir("./Alz"))]
f
