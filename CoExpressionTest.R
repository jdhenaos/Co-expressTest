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

ExtractInfo <- function(x){
  f<-read.csv(gzfile(x),
                     comment.char = "!", 
                     sep = "\t",
                     stringsAsFactors = FALSE)
  return(f)
}

DataUnion <- function(){
  
}

#GEO(read.table("Alzheimer_Chips.txt"),"./Alzheimer_GSE")
#GEO(read.table("Parkinson_Chips.txt"),"./Parkinson_GSE")
#GEO(read.table("MultipleSclerosis_Chips.txt"),"./MultipleSclerosis_GSE")

GEO(read.table("geos.txt"),"./Alz")
GeneSymbol("GPL570","./Alz")

##################################################################

f <- dir(".")[grep("^GSE[0-9]+(_|-GPL570)",dir("."))]
f

for (t in f){
  g <- 0
  if(g == 0){
    print("Vacio")
  }else{
    print("Lleno")
  }
}

GSE68527 <- ExtractInfo("GSE68527_series_matrix.txt.gz")
GSE52139<-read.csv(gzfile("GSE52139_series_matrix.txt.gz"),
                   comment.char = "!", 
                   sep = "\t",
                   stringsAsFactors = FALSE)


dim(GSE68527)
dim(GSE52139)

class(GSE52139)
class(GSE68527)

A <- data.frame(GSE52139,GSE68527)
dim(A)
tail(A$ID_REF)

B <- merge(GSE52139,GSE68527)
dim(B)
head(B$ID_REF)


