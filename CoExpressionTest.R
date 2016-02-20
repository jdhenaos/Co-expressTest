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

GeneSymbol <- function(GPL,dir="."){

}

#GEO(read.table("Alzheimer_Chips.txt"),"./Alzheimer_GSE")
#GEO(read.table("Parkinson_Chips.txt"),"./Parkinson_GSE")
#GEO(read.table("MultipleSclerosis_Chips.txt"),"./MultipleSclerosis_GSE")

GEO(read.table("geos.txt"),"./Alz")

#########################################################
GSE68527<-read.csv(gzfile("GSE68527_series_matrix.txt.gz"),
                   comment.char = "!", 
                   sep = "\t",
                   stringsAsFactors = FALSE)


setwd("./Alz")
gpl3 <- scan(dir(".")[grep("GPL570",dir("."))], comment.char = "!",
                skip = 1,  sep = "\t", character())

head(as.data.frame(gpl3[grep("[#|^]",gpl3,invert = T)], sep = "t"))

head("./GPL570.soft")

gpl <- getGEO(filename = "GPL570.soft")
colnames(Table(gpl))
table(Table(gpl)$"Gene Symbol")
write(table(Table(gpl)$"Gene Symbol"), file = "output")
paste(Table(gpl)$"Gene Symbol"[1],":",table(Table(gpl)$"Gene Symbol"[1]))


