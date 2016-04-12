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

FilterData <- function(fi,gene,Median=FALSE){
  mPD <- rowMeans(fi)
  da <- data.frame(gene,c(0),stringsAsFactors = F)
  names(da) <- c("a","b","c")
  n <- data.frame(names(mPD),mPD, stringsAsFactors = F)
  names(n) <- c("a","b")
  m <- merge.data.frame(n, da, by.x = "a", by.y = "a")
  l <- data.frame(m$b.y,m$b.x, row.names = m$a,stringsAsFactors = F)
  k <- l[-c(grep(paste0("^","$"),l[,1])),]
  j <- unique(k[,1])
  g <- data.frame()
  
  i<- k[grep("DDR1",k[,1],fixed = T),]
  h <- i[1,]
  h[1,2] <- median(i[,2])
  
  if(Median){
    for(x in j){
      i<- k[grep(x,k[,1],fixed = T),]
      h <- i[1,]
      h[1,2] <- median(i[,2])
      if(length(g) == 0){
        g <- rbind(h)
      }else{
        g <- rbind(g,h)
      }
    }
  }else{
    for(x in j){
      i<- k[grep(x,k[,1],fixed = T),]
      h <- i[grep(max(k[grep(x,k[,1],fixed = T),2]),i[,2]),]
      if(length(g) == 0){
        g <- rbind(h)
      }else{
        g <- rbind(g,h)
      }
    }
  }
  return(g)
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

CovarFilter <- function(son,gen,x,d=FALSE){
  mPD <- rowMeans(son)
  h <- data.frame()
  for(n in row.names(gen)){
    gen[n,2] <- sd(son[n,])/gen[n,2]
  }
  
  li <- gen[order(gen$m.b.x, decreasing = d),]
  fil <- li[1:as.integer((dim(li)[1]*x)/100),] 
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

FPD <- SummaryFilter(PD2,"Mean")
FAD <- SummaryFilter(AD2,"Mean")
FMS <- SummaryFilter(MS2,"Mean")

CAD <- CovarFilter(AD,AD2,5,TRUE)
CPD <- CovarFilter(PD,PD2,5,TRUE)
CMS <- CovarFilter(MS,MS2,5,TRUE)

CAD2 <- CovarFilter(AD,AD2,10,TRUE)
CPD2 <- CovarFilter(PD,PD2,10,TRUE)
CMS2 <- CovarFilter(MS,MS2,10,TRUE)

##################################################

gene <- GeneSymbol("GPL570")
PD <- DataUnion("./Parkinson_GSE")
PD2 <- FilterData(PD,gene)
PD3 <- FilterData(PD,gene,Median = T)

mPD <- rowMeans(fi)
da <- data.frame(gene,c(0),stringsAsFactors = F)
names(da) <- c("a","b","c")
n <- data.frame(names(mPD),mPD, stringsAsFactors = F)
names(n) <- c("a","b")
m <- merge.data.frame(n, da, by.x = "a", by.y = "a")
l <- data.frame(m$b.y,m$b.x, row.names = m$a,stringsAsFactors = F)
k <- l[-c(grep(paste0("^","$"),l[,1])),]
j <- unique(k[,1])
g <- data.frame()

i<- k[grep("DDR1",k[,1],fixed = T),]
h <- i[1,]
h[1,2] <- median(i[,2])

if(Median){
  for(x in j){
    i<- k[grep(x,k[,1],fixed = T),]
    h <- i[1,]
    h[1,2] <- median(i[,2])
    if(length(g) == 0){
      g <- rbind(h)
    }else{
      g <- rbind(g,h)
    }
  }
}else{
  for(x in j){
    i<- k[grep(x,k[,1],fixed = T),]
    h <- i[grep(max(k[grep(x,k[,1],fixed = T),2]),i[,2]),]
    if(length(g) == 0){
      g <- rbind(h)
    }else{
      g <- rbind(g,h)
    }
  }
}


