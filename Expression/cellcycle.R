#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/cellcycle genes/"
dir.create(dirp)
setwd(dirp)

rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)

library(RColorBrewer)
library(gplots)

imp<-c("Mki67","Mxi1","Uvrag","Rb1","E2f2","E2f3","E2f6","E2f7")
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

png(filename="cellcycle_small.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

png(filename="cellcycle_big.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

imp_genes_ordered<-imp_genes[order(imp_genes$S1_zscore),]
png(filename="cellcycle_small_ordered.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes_ordered[,which(names(imp_genes_ordered) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes_ordered$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

png(filename="cellcycle_big_ordered.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes_ordered[,which(names(imp_genes_ordered) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes_ordered$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

##########################
a<-c("Ccnd1", "Ccne1", "Cdk4", "Cdkn1b") #G1 Phase & G1/S Transition
b<-c("Mcm2", "Mcm4") #S Phase & DNA Replication
c<-c("Ccnb1", "Cdk5rap1", "Cks1b") #G2 Phase & G2/M Transition
d<-c("Ccnf", "Mre11a", "Rad51") #M Phase


imp<-c(a,b,c,d)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
lab<-c("G1 Phase & G1/S Transition ", "S Phase & DNA Replication ",  "G2 Phase & G2/M Transition ", "M Phase ")

d<-c(rep(1, length(a)), rep(2, length(b)) , rep(3, length(c)), rep(4, length(d)))
n<-length(lab)
col1<-colorRampPalette(brewer.pal(11,"Paired"))
rc<-col1(n)[d]

png(filename="cellcycle_genes_big.png", units="in", width=3, height=3 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1.3,cex.axis=0.5,cexCol=1.3,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          margins = c(2,20),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.6            # line width
)
dev.off()


png(filename="cellcycle_genes_small.png", units="in", width=1.5, height=1.5 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1.3,cex.axis=0.5,cexCol=1.3,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          margins = c(2,20),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.6            # line width
)
dev.off()


png(filename="cell_cycle_zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1, xlim=c(-2.5,2.5))
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cycle_zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cycle_zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cycle_zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()




a<-c("Atm", "Brca1", "Brca2", "Cdk2", "Cdkn1a", "Cdkn1b", "Chek1", "Gadd45a", "Rad9a", "Trp53")#Cell Cycle Checkpoint & Cell Cycle Arrest
b<-c("Rb1", "Uvrag","Foxm1","Mxi1","Cul1") #cell cylce inhibition
c<-c("Pcna","Myc","Mki67")#Cyle progression genes
d<-c("E2f2","E2f3") #DNA Synthesis
e<-c("E2f6","E2f7")#Inhibitor

imp<-c(a,b,c,d,e)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
lab<-c("Cell Cycle Checkpoints ", "Cell Cylce Inhibitors ",  "Proliferative genes ", "DNA Synthesis ","DNA Synthesis Inhibitors ")

d<-c(rep(1, length(a)), rep(2, length(b)) , rep(3, length(c)), rep(4, length(d)),rep(5, length(e)))
n<-length(lab)
col1<-colorRampPalette(brewer.pal(11,"Paired"))
rc<-col1(n)[d]

png(filename="cellcyclereg_genes_big.png", units="in", width=3, height=4 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1.3,cex.axis=0.5,cexCol=1.3,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          margins = c(2,22),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.6            # line width
)
dev.off()


png(filename="cellcycle_genesreg_small.png", units="in", width=1.5, height=2 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1.3,cex.axis=0.5,cexCol=1.3,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          margins = c(2,20),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.6            # line width
)
dev.off()

a<-c("Atm", "Brca1", "Brca2", "Cdk2", "Cdkn1a", "Cdkn1b", "Chek1", "Gadd45a", "Rad9a", "Trp53")#Cell Cycle Checkpoint & Cell Cycle Arrest
b<-c("Rb1", "Uvrag","Foxm1","Mxi1","Cul1") #cell cylce inhibition
c<-c("Pcna","Myc","Mki67")#Cyle progression genes
d<-c("E2f2","E2f3") #DNA Synthesis
e<-c("E2f6","E2f7")#Inhibitor

imp<-c(b,e)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

png(filename="cell_cyclereg_down_zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1, xlim=c(-2.5,2.5))
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cyclereg_down_zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cyclereg_down_zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cyclereg_down_zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()



a<-c("Atm", "Brca1", "Brca2", "Cdk2", "Cdkn1a", "Cdkn1b", "Chek1", "Gadd45a", "Rad9a", "Trp53")#Cell Cycle Checkpoint & Cell Cycle Arrest
b<-c("Rb1", "Uvrag","Foxm1","Mxi1","Cul1") #cell cylce inhibition
c<-c("Pcna","Myc","Mki67")#Cyle progression genes
d<-c("E2f2","E2f3") #DNA Synthesis
e<-c("E2f6","E2f7")#Inhibitor

imp<-c(c,d)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

png(filename="cell_cyclereg_up_zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
plot(d, col="green4", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1, xlim=c(-2.5,2.5))
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
lines(d,col="red",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cyclereg_up_zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
plot(d, col="green4", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
lines(d,col="red",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cyclereg_up_zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
plot(d, col="green4", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
lines(d,col="red",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="cell_cyclereg_up_zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
plot(d, col="green4", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1,xlim=c(-2.5,2.5))
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
lines(d,col="red",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()




