TF_target_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/TF_target_genes.csv", stringsAsFactors = FALSE)
library("gplots")
#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part3/"
dir.create(dirp)
setwd(dirp)
rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)

# take a list of mkl genes
#mkl was list was obtained from Esnault et.al Direct and Near

mkl_gene<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% TF_target_genes[,2])
png(filename="mkl_Zscore_heatmap_Direct_near.png", units="in", width=2, height=8.5 , pointsize=5, res=2200)
par(font.axis=2,font.lab=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.15,
             trace='none',dendrogram="none",cexCol=0.7,srtCol =90,key.title=" ",
             ylab = "",labRow = mkl_gene$SYMBOL,labCol=c("S1","S2","S3","S4"),density.info="none")

dev.off()
png(filename="MKL_Zscore_summed_barplot_Direct_near.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("S1","S2","S3","S4"),ylim=c(x_down,x_up))
box()
abline(h=0)
dev.off()

#mkl was list was obtained from Esnault et.al Direct
mkl_gene<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% TF_target_genes[,3])
png(filename="mkl_Zscore_heatmap_Direct.png", units="in", width=2, height=8.5 , pointsize=5, res=2200)
par(font.axis=2,font.lab=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.3,
             trace='none',dendrogram="none",cexCol=0.7,srtCol =90,key.title=" ",
             ylab = "",labRow = mkl_gene$SYMBOL,labCol=c("S1","S2","S3","S4"),density.info="none")

dev.off()
png(filename="MKL_Zscore_summed_barplot_Direct.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("S1","S2","S3","S4"),ylim=c(x_down,x_up))
box()
abline(h=0)
dev.off()

#mkl was list was obtained from Esnault et.al Stringent
mkl_gene<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% TF_target_genes[,4])
png(filename="mkl_Zscore_heatmap_Stringent.png", units="in", width=2, height=8 , pointsize=5, res=2200)
par(font.axis=2,font.lab=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.2,
             trace='none',dendrogram="none",cexCol=0.7,srtCol =90,key.title=" ",
             ylab = "",labRow = mkl_gene$SYMBOL,labCol=c("S1","S2","S3","S4"),density.info="none")
dev.off()
png(filename="MKL_Zscore_summed_barplot_Stringent.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("S1","S2","S3","S4"),ylim=c(x_down,x_up))
box()
abline(h=0)
dev.off()

#mkl was list was obtained from Esnault et.al Stringent Direct
mkl_gene<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% TF_target_genes[,5])
png(filename="mkl_Zscore_heatmap_Stringent_Direct.png", units="in", width=2, height=7.5 , pointsize=5, res=2200)
par(font.axis=2,font.lab=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.3,
             trace='none',dendrogram="none",cexCol=0.7,srtCol =90,key.title=" ",
             ylab = "",labRow = mkl_gene$SYMBOL,labCol=c("S1","S2","S3","S4"),density.info="none")
dev.off()

png(filename="MKL_Zscore_summed_barplot_Stringent_Direct.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("S1","S2","S3","S4"),ylim=c(x_down,x_up))
box()
abline(h=0)
dev.off()

#mkl was list was obtained from Selvam et.al RT_PCR

mkl_gene<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% TF_target_genes[,6])

png(filename="mkl_Zscore_heatmap_Stringent_Direct_RT_PCR.png", units="in", width=2, height=2 , pointsize=5, res=2200)
par(font.axis=2,font.lab=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.8,
             trace='none',dendrogram="none",cexCol=0.7,srtCol =90,key.title=" ",
             ylab = "",labRow = mkl_gene$SYMBOL,labCol=c("S1","S2","S3","S4"),density.info="none")

dev.off()

png(filename="MKL_Zscore_summed_barplot_Stringent_DirectRT_PCR.png", units="in", width=1.5, height=1.5 , pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("S1","S2","S3","S4"),ylim=c(x_down,x_up))
box()
abline(h=0)
dev.off()


#mkl was list was obtained from combining all the data
temp<-unique(c(TF_target_genes[,2],TF_target_genes[,3],TF_target_genes[,4],TF_target_genes[,5],TF_target_genes[,6]))
mkl_gene<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% temp)

png(filename="mkl_Zscore_heatmap_all.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.15,
             trace='none',dendrogram="none",cexCol=1,srtCol =0,key.title="zscore",
             ylab = "",labRow =" ",labCol=c("R","RL","C","CL"),density.info="none",
             key.par=list(mar=c(1,1,4,1), cex.axis=0.8), cex.axis=0.8, key.xlab = "",keysize=1.1,font=2)

dev.off()

png(filename="MKL_Zscore_summed_barplot_all.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("R","RL","C","CL"),ylim=c(x_down,x_up),cex.names = 1.2, cex.lab=1.5, cex.axis=1.2)
box()
abline(h=0)
dev.off()

png(filename="mkl_Zscore_heatmap_all_small.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
mkl_gene<-mkl_gene[order(mkl_gene$S1_zscore),]
x<-heatmap.2(as.matrix(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redgreen(75),
             scale="none",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.15,
             trace='none',dendrogram="none",cexCol=1,srtCol =0,key.title="zscore",
             ylab = "",labRow =" ",labCol=c("R","RL","C","CL"),density.info="none",
             key.par=list(mar=c(1,1,4,1), cex.axis=0.8), cex.axis=0.8, key.xlab = "",keysize=1.1,font=2)

dev.off()

png(filename="MKL_Zscore_summed_barplot_all_small.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-colSums(mkl_gene[,which(names(mkl_gene) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))],na.rm=T)
x_down<-min(x)-5
x_up<-max(x)+5
barplot(x,las=1,ylab="Summed Zscore",names =c("R","RL","C","CL"),ylim=c(x_down,x_up),cex.names = 1.2, cex.lab=1.5, cex.axis=1.2)
box()
abline(h=0)
dev.off()

write.csv(mkl_gene, file="MKL_gene_transcription.csv")