#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/autophagy and other pathways/"
dir.create(dirp)
setwd(dirp)

rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)

library(RColorBrewer)
library(gplots)


##########AUTOPHAGY MACHINERY#####


a<-c("Atg12", "Atg16l1", "Atg4a", "Atg4b", "Atg4c", "Atg4d", "Atg5", "Atg9a", "Atg9b", "Becn1", "Gabarap", "Gabarapl1", "Gabarapl2", 
     "Irgm1", "Map1lc3a", "Map1lc3b", "Rgs19", "Ulk1", "Wipi1")
b<-c("Atg4a", "Atg4b", "Atg4c", "Atg4d", "Gabarap")
c<-c( "Atg10", "Atg16l1", "Atg16l2", "Atg3", "Atg4a", "Atg4b", "Atg4c", "Atg4d", "Atg7", "Atg9a", "Gabarap", "Gabarapl2", "Rab24")
d<-c("Dram1", "Gabarap", "Lamp1", "Npc1")
e<-c("Atg3", "Atg7", "Hdac6")
f<-c("Atg4a", "Atg4b", "Atg4c", "Atg4d")

imp<-c(a,b,c,d,e,f)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
lab<-c("Autophagic Vacuole Formation","Vacuole Targeting","Protein Transport","Autophagosome-Lysosome Linkage","Ubiquitination","Proteases")

d<-c(rep(1, length(a)), rep(2, length(b)) , rep(3, length(c)), rep(4, length(d)),rep(5, length(e)),rep(6, length(f)))
n<-length(lab)
col1<-colorRampPalette(brewer.pal(11,"Paired"))
rc<-col1(n)[d]

png(filename="autophagy_genes_big.png", units="in", width=3, height=6 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1.3,cex.axis=0.5,cexCol=1.3,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.4,key.par=list(mar=c(1,1,5,0.5),cex=0.6, font=2),
          margins = c(2,25),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.6            # line width
)
dev.off()


png(filename="autophagy_small.png", units="in", width=1.5, height=3 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1.3,cex.axis=0.5,cexCol=1.3,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.4,key.par=list(mar=c(1,1,5,0.5),cex=0.6, font=2),
          margins = c(2,25),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.6            # line width
)
dev.off()


png(filename="autophagy_zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="autophagy_zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="autophagy_zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="autophagy_zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()



########## #############

imp<-c("Atr","Atm","Parp1", "Pten", "Foxo6", "Foxo3", "Akt3","Mtor","Atg3", "Cdc25b", "Gabarapl2", "Pik3c3")
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

png(filename="Pathways.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.2,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)), margins=c(2,8))

dev.off()

png(filename="Pathways_big.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.2,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)), margins=c(2,8))

dev.off()
