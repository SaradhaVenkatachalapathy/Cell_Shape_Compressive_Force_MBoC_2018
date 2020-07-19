#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/quiescence genes/"
dir.create(dirp)
setwd(dirp)

rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)
#important genes

imp<-c( "Tob1","Tob2","Btg2","Arid5a",#T cell quiencence
        "Foxo3", "Klf2","Mir29a","Mir221","Mirlet7a-1","Sod2" )

library(gplots)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
png(filename="Quiescence.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

png(filename="Quiescence_big.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

imp_genes_ordered<-imp_genes[order(imp_genes$S1_zscore),]
png(filename="Quiescence_ordered.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes_ordered[,which(names(imp_genes_ordered) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes_ordered$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

png(filename="Quiescence_big_ordered.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes_ordered[,which(names(imp_genes_ordered) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes_ordered$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()


#mir29a targets
imp<-c("Arrdc4","Blmh","Cdk6","Col1a1","Col3a1","Col5a2","Fbn1","Lamc1","Ppic","Serpinh1","Sparc","Fstl1")
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

library(gplots)
png(filename="mir29a_target_small.png", units="in", width=1, height=1.5 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=1, key.par=list(mar=c(2,1,1,1)),lhei=c(0.5,4), lwid=c(1.5,4))

dev.off()

png(filename="mir29a_target_big.png", units="in", width=2, height=3 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=1, key.par=list(mar=c(2,1,1,1)),lhei=c(0.5,4), lwid=c(1.5,4))

dev.off()


png(filename="mir29a_target_zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="mir29a_target_zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="mir29a_target_zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="mir29a_target_zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

########################################################################
#Molecular regulation of stem cell quiescence
cell_cycle<-c("Anln", "Birc5", "Ccna2", "Ccnb1", "Ccne2", "Sgol1","Ccnd3","Pdk1")
DNA<-c("Mcm4","Pcna","Rrm2","Top2a")
mito<-c("Cycs","Mtch2","Slc25a5")
trasn<-c("Foxo3","Ezh1","Prdm5","Ptov1","Zfp30","Zbtb20","Phf1","Ctdsp1","Thra","Tef")
RNA<-c("Ddx39","Dicer1")
imp<-c(cell_cycle,DNA,mito,trasn,RNA)

imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

library(gplots)
png(filename="stem_cell_quie.png", units="in", width=2, height=3 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.6, key.par=list(mar=c(2,1,1,1)))

dev.off()

png(filename="stem_cell_quie_zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="stem_cell_quie_zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="stem_cell_quie_zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="stem_cell_quie_zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()


########################################################################
#A New Description of Cellular Quiescence, Hilary A. Coller, PLosOne, 2006

a<-c("Ugcg", "Gdf5", "Igfbp3")
b<-c("Hmgcl", "Aldh2", "Jarid2", "Ppp2r5a")
c<-c("Pask")
d<-c("Slc20a2")
e<-c("Rgs2","Cdk10","Map3k5","Itpka", "Tmem115", "Nr4a1","Nr4a2")
f<-c("Cul2", "Cenpa", "Lmnb1", "Lbr", "Fer","Dhx9","Fasn","Prps1","Umps","Fntb","Mrpl12","Kpna3","Pcna","Mettl1","Spfq","Sfrs2","Sf3a1","Ddx30","Bop1","Inpp5a")
h<-c("Arg2","Cul1","Prkcm", "Kiaa0092","Pcyt2","Prkar1b","Hist2h4", "Map2k5","Fgf4", "Ccne1", "Crebl1", "Adora2b")
i<-c("Np", "Mrpl19", "Mphosph6")
j<-"Etfdh"
k<-c("Asl", "Slc1a4", "Umpk", "Lpin1", "Mvd", "Idi1", "Scap", "Fgf5", "Slc12a4", "Snrpc")


imp<-c("Ugcg", "Gdf5", "Igfbp3", #up in mitogen withdrawal1
       "Hmgcl", "Aldh2", "Jarid2", "Ppp2r5a", #up in contact inhibition 2
       "Pask", # up loss of adhesion3
       "Slc20a2",# up loss of adhesion and mitogen withdrawal4
       "Rgs2","Cdk10","Map3k5","Itpka", "Tmem115", "Nr4a1","Nr4a2", #down in mitogen withdrawal5
       "Cul2", "Cenpa", "Lmnb1", "Lbr", "Fer","Dhx9","Fasn","Prps1","Umps","Fntb",  #down in contact inhibition 
       "Mrpl12","Kpna3","Pcna","Mettl1","Sfpq","Srsf2","Sf3a1","Dhx30","Bop1","Inpp5a", #down in contact inhibition 6
       "Arg2","Cul1","Prkd1", "Cep57","Pcyt2","Prkar1b","Hist2h4", "Map2k5","Fgf4", "Ccne1", "Atf6b", "Adora2b", #down in adhesion 8
       "Pnp", "Mrpl19", "Mphosph6", #down in mitogen and contact inhitiobion 9
       "Etfdh", #down in mitogen and adhesion 10
       "Asl", "Slc1a4", "Cmpk1", "Lpin1", "Mvd", "Idi1", "Scap", "Fgf5", "Slc12a4", "Snrpc"#down in adhesion and contact inhitiobion 11
)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]


library(RColorBrewer)
lab<-c("up in mitogen withdrawal", "up in contact inhibition","up in loss of adhesion", "up loss of adhesion and mitogen withdrawal",
       "down in mitogen withdrawal", "down in contact inhibition","down in adhesion loss",  "down in mitogen withdrawal and contact inhibition", 
       "down in mitogen withdrawal and adhesion loss","down in adhesion loss and contact inhibition" )

d<-c(rep(1, length(a)), rep(2, length(b)) , rep(3, length(c)), rep(4, length(d)), rep(5, length(e)), rep(6, length(f)),
     rep(7, length(h)), rep(8, length(i)),rep(9, length(j)),rep(10, length(k)))
n<-length(lab)
col1<-colorRampPalette(brewer.pal(10,"Paired"))
rc<-col1(n)[d]


png(filename="quiescent_genes_Colleretal.png", units="in", width=3, height=3.5 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1,cex.axis=0.5,cexCol=1,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          margins = c(2,25),col=redgreen(100), RowSideColors = rc)

par(lend = 1)           # square line ends for the color legend
legend("bottomright",      # location of the legend on the heatmap plot
       legend = lab, # category labels
       col = col1(n),  # color key
       lty= 1,             # line style
       lwd = 4,cex=0.5            # line width
)
dev.off()




imp<-c("Ugcg", "Gdf5", "Igfbp3", #up in mitogen withdrawal1
       "Hmgcl", "Aldh2", "Jarid2", "Ppp2r5a", #up in contact inhibition 2
       "Pask", # up loss of adhesion3
       "Slc20a2"# up loss of adhesion and mitogen withdrawal4
)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]


png(filename="quiescent_genes_Colleretal_UP.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1,cex.axis=0.5,cexCol=1,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=1,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          col=redgreen(100))
dev.off()

png(filename="Colleretal_quie_zscore_density_small_circle_UP.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="Colleretal_quie_zscore_density_small_rect_UP.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()


imp<-c("Rgs2","Cdk10","Map3k5","Itpka", "Tmem115", "Nr4a1","Nr4a2", #down in mitogen withdrawal5
       "Cul2", "Cenpa", "Lmnb1", "Lbr", "Fer","Dhx9","Fasn","Prps1","Umps","Fntb",  #down in contact inhibition 
       "Mrpl12","Kpna3","Pcna","Mettl1","Sfpq","Srsf2","Sf3a1","Dhx30","Bop1","Inpp5a", #down in contact inhibition 6
       "Arg2","Cul1","Prkd1", "Cep57","Pcyt2","Prkar1b","Hist2h4", "Map2k5","Fgf4", "Ccne1", "Atf6b", "Adora2b", #down in adhesion 8
       "Pnp", "Mrpl19", "Mphosph6", #down in mitogen and contact inhitiobion 9
       "Etfdh", #down in mitogen and adhesion 10
       "Asl", "Slc1a4", "Cmpk1", "Lpin1", "Mvd", "Idi1", "Scap", "Fgf5", "Slc12a4", "Snrpc"#down in adhesion and contact inhitiobion 11
)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]


png(filename="quiescent_genes_Colleretal_down.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), 
          scale="none",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
          cexRow=1,cex.axis=0.5,cexCol=1,srtCol =0,ylab = "",
          labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),
          key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", 
          keysize=1,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
          col=redgreen(100))
dev.off()

png(filename="Colleretal_quie_zscore_density_small_circle_DOWN.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="Colleretal_quie_zscore_density_small_rect_DOWN.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((imp_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((imp_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

