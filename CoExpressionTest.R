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
  sym <- table(Table(gpl)$"Gene Symbol")
  write.table(table(Table(gpl)$"Gene Symbol"), file = "GeneSymbol.txt")
  setwd("../")
  return(sym)
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
      g <- merge.data.frame(g,h,by.x = "ID_REF",by.y = "ID_REF",all = TRUE)
    }
  }
  ng <- g[,2:dim(g)[2]]
  t <- sapply(ng, as.numeric)
  tt <- as.data.frame(t, row.names = g$ID_REF)
  y <- rowMeans(tt, na.rm = T)
  return(y)
}

#GEO(read.table("Alzheimer_Chips.txt"),"./Alzheimer_GSE")
#GEO(read.table("Parkinson_Chips.txt"),"./Parkinson_GSE")
#GEO(read.table("MultipleSclerosis_Chips.txt"),"./MultipleSclerosis_GSE")

GEO(read.table("geos.txt"),"./Alz")
gene <- GeneSymbol("GPL570")
D <- DataUnion("./Alz")

AD <- DataUnion("./Alzheimer_GSE")
PD <- DataUnion("./Parkinson_GSE")
MS <- DataUnion("./MultipleSclerosis_GSE")
dim(PD)
dim(AD)
dim(MS)

a <- summary(PD)
f <- dir(".")[grep("GPL570",dir("."))]
f <- f[grep(".soft$",f)]
gpl <- getGEO(filename = f)
sym <- Table(gpl)
p <- data.frame(sym$`Gene Symbol`, stringsAsFactors = F)
q <- data.frame(p, names(PD), PD, stringsAsFactors = F)
r <- subset(q, q$PD >= a[5])
s <- data.frame(unique(r$sym..Gene.Symbol.), c(0), stringsAsFactors = F)

write.table(s, file = "vacio.txt")

for(i in as.vector(s[,1])){
  s[grep(i,s$unique.r.sym..Gene.Symbol..),2] <- max(r[grep(i,r$sym..Gene.Symbol.),3])
}
