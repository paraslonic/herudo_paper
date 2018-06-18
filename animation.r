library("data.table")
library("sp")
library(rgl)
library(scatterplot3d)
## READ -------------------------------------------------------------------

GC = fread("draft1_cov_GC")
colnames(GC) = c("contig","GC")

## genome coverage
proton = fread("proton.cov")
proton[,cov:=V4/V3]
illumina = fread("illumina.cov")
illumina[,cov:=V4/V3]
GCOV=data.table(contig=proton$V1, proton_cov = proton$cov, illumina_cov = illumina$cov)

## rna
rna = fread("rna.cov")
rna[,cov:=V4/V3]
RNA=data.table(contig=rna$V1, rna_reads=rna$V4, rna_cov=rna$cov)

## blast
blast = fread("megan.txt", head = FALSE)
blast[,bacteria:=grepl("Bacteria",blast$V2)]
blast[,eukaryota:=grepl("Eukaryota",blast$V2)]
TAX=data.table(contig=blast$V1,bacteria=blast$bacteria,eukaryota=blast$eukaryota)

# C O M B I N E -----------------------------------------------------------
T = merge(GCOV,TAX,by="contig",all = TRUE)
T = merge(T,RNA,by="contig")
T = merge(T,GC,by="contig")
rownames(T) = T$contig
T[,contig:=NULL]
T[,rna_cov:=NULL]
T.bac = T[,"bacteria",with=FALSE]
T.euk = T[,"eukaryota",with=FALSE]
T$illumina_cov = T$illumina_cov*median(T$proton_cov/T$illumina_cov)
hist(T$proton_cov/T$illumina_cov, breaks = 100000, xlim = c(0,10))

# PLOT -------------------------------------------------------------------
colors = rep(rgb(0.1,0.1,0.1,0.1),nrow(T))
colors[T$eukaryota]=rgb(156/255,196/255,217/255,1)
colors[T$bacteria]=rgb(178/255,223/255,138/255,0.9)

gc = T$GC
cov = log(T$illumina_cov+T$proton_cov)
plot(gc,cov,col = colors, cex = 0.4, pch = 19, axes = FALSE, xlim = c(8, 80),
     ylab = "WGS reads coverage, log", xlab = "GC", cex.lab = 1.2)
plot3d(gc, cov, (T$rna_reads), col = colors)
plot3d(gc, cov, log(1+T$rna_reads), col = colors,size = 1, 
       ylab = "WGS reads coverage, log", xlab = "GC", zlab="RNA read counts")

for (i in 1:180) {
  view3d(userMatrix=rotationMatrix(2*pi * (i)/180, 1, -1, -1))
  rgl.snapshot(filename=paste("animation/frame-",
                              sprintf("%03d", i), ".png", sep=""))
}
rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',i,'plot.png',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',i,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', i,'plot.png', sep=''))
  }
}

for(i in seq(1,360,2)){
  png(rename(i))
  scatterplot3d(gc, cov, log(1+T$rna_reads), angle = i, pch =19, color = colors, cex.symbols =  0.3,
                ylab = "WGS reads coverage, log", xlab = "GC", zlab="RNA read counts")
dev.off()
  }

