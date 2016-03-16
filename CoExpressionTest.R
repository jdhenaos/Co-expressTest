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
  ta <- data.frame(sym$ID, sym$`Gene Symbol`, stringsAsFactors = F)
  #write.table(table(Table(gpl)$"Gene Symbol"), file = "GeneSymbol.txt")
  #setwd("../")
  return(ta)
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
  return(tt)
}

FilterData <- function(fi,gene){
  da <- data.frame(gene,c(0),stringsAsFactors = F)
  n <- data.frame(names(fi),fi, stringsAsFactors = F)
  m <- merge.data.frame(n, da, by.x = "names.fi.", by.y = "sym.ID")
  l <- data.frame(m$sym..Gene.Symbol.,m$fi,stringsAsFactors = F)
  k <- l[-c(grep(paste0("^","$"),l[,1])),]
  j <- data.frame(unique(k$m.sym..Gene.Symbol.), c(0), stringsAsFactors = F)
  
  for(i in as.vector(j[,1])){
    j[grep(paste0("^",i,"$"),j[,1]),2] <-
      max(l[grep(paste0("^",i,"$"),l[,1]),2])
  }
  return(j)
}

SummaryFilter <- function(Data,Qu){
  su <- summary(Data[,2])
  
  if(Qu == 1){
    fil <- Data[Data[,2] >= su[2],]
  }else if(Qu == 2){
    fil <- Data[Data[,2] >= su[3],]
  }else if(Qu == 3){
    fil <- Data[Data[,2] >= su[5],]
  }else if(Qu == "Mean"){
    fil <- Data[Data[,2] >= su[4],]
  }else{
    fil <- NULL
    warning("Opcion incorrecta", immediate. = T, noBreaks. = T)
  }
  return(fil)
}

#GEO(read.table("Alzheimer_Chips.txt"),"./Alzheimer_GSE")
#GEO(read.table("Parkinson_Chips.txt"),"./Parkinson_GSE")
#GEO(read.table("MultipleSclerosis_Chips.txt"),"./MultipleSclerosis_GSE")
gene <- GeneSymbol("GPL570")

AD <- DataUnion("./Alzheimer_GSE")
PD <- DataUnion("./Parkinson_GSE")
MS <- DataUnion("./MultipleSclerosis_GSE")

mPD <- rowMeans(PD)

AD2 <- FilterData(rowMeans(AD),gene)
PD2 <- FilterData(rowMeans(PD),gene)
MS2 <- FilterData(rowMeans(MS),gene)

FPD <- SummaryFilter(PD2,"Mean")
FAD <- SummaryFilter(AD2,"Mean")
FMS <- SummaryFilter(MS2,"Mean")

##################################################

y <- names(mPD[grep(PD2[3,2],mPD)])
x <- PD[grep(paste0("^",y,"$"),row.names(PD)),]
ds <- sd(x)
total <- ds/mPD[y]


