library("data.table")
library("car")
library("sp")

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

# PLOT -------------------------------------------------------------------
colors = rep(rgb(0.1,0.1,0.1,0.1),nrow(T))
colors[T$eukaryota]=rgb(156/255,196/255,217/255,1)
colors[T$bacteria]=rgb(178/255,223/255,138/255,0.9)

gc = T$GC
cov = log(T$illumina_cov+T$proton_cov)
pdf("bining.pdf")
plot(gc,cov,col = colors, cex = 0.4, pch = 19, axes = FALSE, xlim = c(8, 80),
     ylab = "WGS reads coverage, log", xlab = "GC", cex.lab = 1.2)

axis(side=2, at=c(-5:5), col = "gray60")
axis(side=1, at=seq(10,80,10), col = "gray60")
legend("topleft",c("Bacteria","Eukaryota"), pch=19, col=c(rgb(178/255,223/255,138/255), rgb(156/255,196/255,217/255)), bty = "n", cex = 1.2)


gc.host = gc[T$rna_reads>10]
cov.host = cov[T$rna_reads>10]
el = dataEllipse(gc.host, cov.host, levels=c(0.99999,0.99999),center.cex = 0.1, 
            col = rgb(31/255,120/255,180/255),  lwd = 3, grid=F, plot.points = FALSE)
el 
dev.off()
# SELECT HOST -------------------------------------------------------------

host = point.in.polygon(gc,cov,el[[1]][,1],el[[1]][,2])==1
host[T$rna_reads > 10] = TRUE
host[which(T$bacteria==TRUE)] = FALSE
host[which(T$eukaryota==TRUE)] = TRUE
write.table(rownames(T)[host], "host_contigs",quote = FALSE, row.names=FALSE, col.names=FALSE)
write.table(rownames(T)[!host], "not_host_contigs",quote = FALSE, row.names=FALSE, col.names=FALSE)
bacteria <- rownames(T)[T$bacteria]
write.table(bacteria[complete.cases(bacteria)], "bacteria_contigs",quote = FALSE, row.names=FALSE, col.names=FALSE)




