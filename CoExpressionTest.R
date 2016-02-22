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
  f <- dir(".")[grep(GPL,dir("."))]
  f <- f[grep(".soft$",f)]
  gpl <- getGEO(filename = f)
  write.table(table(Table(gpl)$"Gene Symbol"), file = "GeneSymbol.txt")
  setwd("../")
}

ExtractInfo <- function(x,d = "."){
  #setwd(d)
  f<-read.csv(gzfile(x),
                     comment.char = "!", 
                     sep = "\t",
                     stringsAsFactors = FALSE)
  #setwd("../")
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
GeneSymbol("GPL570")
D <- DataUnion("./Alz")

AD <- DataUnion("./Alzheimer_GSE")
PD <- DataUnion("./Parkinson_GSE")
MS <- DataUnion("./MultipleSclerosis_GSE")
dim(PD)
dim(AD)
dim(MS)

#########################################################

GSE20146 <- ExtractInfo("GSE20146_series_matrix.txt.gz")
GSE14711 <- ExtractInfo("GSE14711_series_matrix.txt.gz")
GSE20141 <- ExtractInfo("GSE20141_series_matrix.txt.gz")
GSE20153 <- ExtractInfo("GSE20153_series_matrix.txt.gz")
GSE30792 <- ExtractInfo("GSE30792_series_matrix.txt.gz")
GSE4773 <- ExtractInfo("GSE4773_series_matrix.txt.gz")
GSE49036 <- ExtractInfo("GSE49036_series_matrix.txt.gz")
GSE7621 <- ExtractInfo("GSE7621_series_matrix.txt.gz")
GSE9807 <- ExtractInfo("GSE9807_series_matrix.txt.gz")

GSE4757 <- ExtractInfo("GSE4757_series_matrix.txt.gz")
dim(GSE4757)

dim(GSE20146)
dim(GSE14711)
dim(GSE20141)
dim(GSE20153)
dim(GSE30792)
dim(GSE4773)
dim(GSE49036)
dim(GSE7621)
dim(GSE9807)

t1 <- (GSE49036$ID_REF)
t2 <- (GSE7621$ID_REF)

write(paste(t1,":",t2),file = "compare")
######
f <- dir(".")[grep("^GSE[0-9]+(_|-GPL570)",dir("."))]
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

if(GSE4757[,1][1] == GSE9807[,1][1]){
  test <- cbind(GSE4757[1,],GSE9807[2,])
}