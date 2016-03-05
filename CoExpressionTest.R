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
  sym <- Table(gpl)
  #write.table(table(Table(gpl)$"Gene Symbol"), file = "GeneSymbol.txt")
  #setwd("../")
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
      g <- merge.data.frame(g,h)
    }
  }
  ng <- g[,2:dim(g)[2]]
  t <- sapply(ng, as.numeric)
  tt <- as.data.frame(t, row.names = g$ID_REF)
  y <- rowMeans(tt)
  return(y)
}

FilterData <- function(fi,gene){
  #a <- summary(fi)
  #p <- data.frame(gene$`Gene Symbol`, stringsAsFactors = F)
  #q <- data.frame(p, names(fi), fi, stringsAsFactors = F)
  #r <- subset(q, q$fi >= a[5])
  #s <- data.frame(unique(q$gene..Gene.Symbol.), c(0), stringsAsFactors = F)
  
  
  
  for(i in as.vector(s[,1])){
    s[grep(paste0("^",i,"$"),s$unique.r.gene..Gene.Symbol..),2] <- 
      max(r[grep(paste0("^",i,"$"),r$gene..Gene.Symbol.),3])
  }
  y <- s[-c(grep("^$",s$unique.r.gene..Gene.Symbol..)),]
  return(y)
}

#GEO(read.table("Alzheimer_Chips.txt"),"./Alzheimer_GSE")
#GEO(read.table("Parkinson_Chips.txt"),"./Parkinson_GSE")
#GEO(read.table("MultipleSclerosis_Chips.txt"),"./MultipleSclerosis_GSE")
gene <- GeneSymbol("GPL570")

AD <- DataUnion("./Alzheimer_GSE")
PD <- DataUnion("./Parkinson_GSE")
MS <- DataUnion("./MultipleSclerosis_GSE")

AD2 <- FilterData(AD,gene)
PD2 <- FilterData(PD,gene)
MS2 <- FilterData(MS,gene)

f <- dir(".")[grep("GPL570",dir("."))]
f <- f[grep(".soft$",f)]
gpl <- getGEO(filename = f)
sym <- Table(gpl)
ta <- data.frame(sym$ID, sym$`Gene Symbol`, stringsAsFactors = F)

da <- data.frame(ta,c(0),stringsAsFactors = F)
for(i in as.vector(s[,1])){
  s[grep(paste0("^",i,"$"),s$unique.r.gene..Gene.Symbol..),2] <- 
    max(r[grep(paste0("^",i,"$"),r$gene..Gene.Symbol.),3])
}
