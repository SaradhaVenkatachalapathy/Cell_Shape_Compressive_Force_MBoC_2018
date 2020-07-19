#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/types genes/"
dir.create(dirp)
setwd(dirp)

rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)
#important genes

imp<-c( "Pik3c3","Rad51b","Slc6a9","Asap2","Myo5c","Ttc9",
        "Necab3", "Foxn3","Dlx3","Map4k3" )

library(gplots)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
png(filename="geometry.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

imp<-c( "Ccdc80" ,"Klc3", "Apex1", 
        "Rhof","Ucp3",
        "Galnt11","Ucn2","Kdm4b",
        "Mmp3","Ttc9")

library(gplots)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
png(filename="geometry_load.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()

imp<-c( "S100a2" ,"Pik3c3", "Rad51b", 
        "Nadsyn1","Ttc9","Myo5c",
        "Rangap1","Rtfdc1","Agfg2",
        "G0s2","Hoxb4","Klc3",
        "Actl7b")

library(gplots)
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]
png(filename="geometry_load2.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
             trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
             ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
             keysize=0.9, key.par=list(mar=c(2,1,2,1)))

dev.off()



x<-heatmap.2(as.matrix(log2(imp_genes[,which(names(imp_genes) %in% c("S1", "S2", "S3","S4"))])), col=redgreen(75),
              scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1.3,cex.axis=0.5,
              trace='none',dendrogram="none",cexCol=1.3,srtCol =0,key.title=" ",key.xlab = "",
              ylab = "",labRow = imp_genes$SYMBOL,labCol=c("R"," RL"," C","  CL"),density.info="none", 
              keysize=0.9, key.par=list(mar=c(2,1,2,1)))
