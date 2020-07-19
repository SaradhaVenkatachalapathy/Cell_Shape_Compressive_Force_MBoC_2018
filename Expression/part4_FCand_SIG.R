#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part4_FC_SIG/"
dir.create(dirp)
setwd(dirp)

rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)

#Perform paired wilcox test on the genes normalised with TPKM
#Wilcoxon tests (as opposed to t-tests, e.g.) are non-parametric. They do not assume anything about the distributions. When the data are very non-normal, t-tests may not be appropriate and one alternative is Wilcoxon. In addition, Wilcoxon tests entire distributions whereas t-tests are tests of the means. It is also known as Mann-Whitney test. Significant genes are identified with p.value < 0.1

x<-c("SYMBOL","S1_B1_TPKM","S2_B1_TPKM","S3_B1_TPKM","S4_B1_TPKM","S1_B2_TPKM","S2_B2_TPKM",
     "S3_B2_TPKM","S4_B2_TPKM","S1_B3_TPKM","S2_B3_TPKM","S3_B3_TPKM","S4_B3_TPKM")

#Wilcox tests
d<-as.matrix(rnaseq_mapped_genes[, which(names(rnaseq_mapped_genes) %in% x)])
row.names(d)<-d[,1]
{
  wil_S1S2<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(2,6,10)]),as.numeric(row[c(3,7,11)]))[c("p.value","statistic")]))
  wil_S1S2<-as.data.frame(t(wil_S1S2))
  wil_S1S2[,3]<-as.character(rownames(wil_S1S2))
  S1S2_num<-nrow(subset(wil_S1S2,wil_S1S2[,1]<0.1))
  S1S2_gene<-subset(wil_S1S2,wil_S1S2[,1]<0.1)[,3]
  S1S2<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S1S2_gene)
  FC_S1S2<-subset(S1S2,(S1S2$FC_S2_S1)>=4 |(S1S2$FC_S2_S1)<=0.25)
  S1S2_num<-nrow(FC_S1S2)
  S1S2_gene<-FC_S1S2$SYMBOL
  
  wil_S1S3<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(2,6,10)]),as.numeric(row[c(4,8,12)]))[c("p.value","statistic")]))
  wil_S1S3<-as.data.frame(t(wil_S1S3))
  wil_S1S3[,3]<-as.character(rownames(wil_S1S3))
  S1S3_num<-nrow(subset(wil_S1S3,wil_S1S3[,1]<0.1))
  S1S3_gene<-subset(wil_S1S3,wil_S1S3[,1]<0.1)[,3]
  S1S3<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S1S3_gene)
  FC_S1S3<-subset(S1S3,(S1S3$FC_S3_S1)>=4 |(S1S3$FC_S3_S1)<=0.25)
  S1S3_num<-nrow(FC_S1S3)
  S1S3_gene<-FC_S1S3$SYMBOL
  
  wil_S1S4<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(2,6,10)]),as.numeric(row[c(5,9,13)]))[c("p.value","statistic")]))
  wil_S1S4<-as.data.frame(t(wil_S1S4))
  wil_S1S4[,3]<-as.character(rownames(wil_S1S4))
  S1S4_num<-nrow(subset(wil_S1S4,wil_S1S4[,1]<0.1))
  S1S4_gene<-subset(wil_S1S4,wil_S1S4[,1]<0.1)[,3]
  S1S4<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S1S4_gene)
  FC_S1S4<-subset(S1S4,(S1S4$FC_S4_S1)>=4 |(S1S4$FC_S4_S1)<=0.25)
  S1S4_num<-nrow(FC_S1S4)
  S1S4_gene<-FC_S1S4$SYMBOL
  
  wil_S2S3<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(3,7,11)]),as.numeric(row[c(4,8,12)]))[c("p.value","statistic")]))
  wil_S2S3<-as.data.frame(t(wil_S2S3))
  wil_S2S3[,3]<-as.character(rownames(wil_S2S3))
  S2S3_num<-nrow(subset(wil_S2S3,wil_S2S3[,1]<0.1))
  S2S3_gene<-subset(wil_S2S3,wil_S2S3[,1]<0.1)[,3]
  S2S3<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S2S3_gene)
  FC_S2S3<-subset(S2S3,(S2S3$FC_S3_S2)>=4 |(S2S3$FC_S3_S2)<=0.25)
  S2S3_num<-nrow(FC_S2S3)
  S2S3_gene<-FC_S2S3$SYMBOL
  
  wil_S2S4<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(3,7,11)]),as.numeric(row[c(5,9,13)]))[c("p.value","statistic")]))
  wil_S2S4<-as.data.frame(t(wil_S2S4))
  wil_S2S4[,3]<-as.character(rownames(wil_S2S4))
  S2S4_num<-nrow(subset(wil_S2S4,wil_S2S4[,1]<0.1))
  S2S4_gene<-subset(wil_S2S4,wil_S2S4[,1]<0.1)[,3]
  S2S4<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S2S4_gene)
  FC_S2S4<-subset(S2S4,(S2S4$FC_S4_S2)>=4 |(S2S4$FC_S4_S2)<=0.25)
  S2S4_num<-nrow(FC_S2S4)
  S2S4_gene<-FC_S2S4$SYMBOL
  
  wil_S3S4<-apply(d, 1, function(row) unlist(wilcox.test(as.numeric(row[c(4,8,12)]),as.numeric(row[c(5,9,13)]))[c("p.value","statistic")]))
  wil_S3S4<-as.data.frame(t(wil_S3S4))
  wil_S3S4[,3]<-as.character(rownames(wil_S3S4))
  S3S4_num<-nrow(subset(wil_S3S4,wil_S3S4[,1]<0.1))
  S3S4_gene<-subset(wil_S3S4,wil_S3S4[,1]<0.1)[,3]
  S3S4<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S3S4_gene)
  FC_S3S4<-subset(S3S4,(S3S4$FC_S4_S3)>=4 |(S3S4$FC_S4_S3)<=0.25)
  S3S4_num<-nrow(FC_S3S4)
  S3S4_gene<-FC_S3S4$SYMBOL
  
}

#number of DE genes
{
  Wilcox_number = matrix( c(NA,S1S2_num,S1S3_num,S1S4_num,
                            S1S2_num,NA,S2S3_num,S2S4_num,
                            S1S3_num,S2S3_num,NA,S3S4_num,
                            S1S4_num,S2S4_num,S3S4_num,NA), 
                          nrow=4, ncol=4) 
  Wilcox_number<-as.data.frame(Wilcox_number)
  colnames(Wilcox_number)<-c("S1","S2","S3","S4")
  rownames(Wilcox_number)<-c("S1","S2","S3","S4")
  as.matrix(Wilcox_number)
  library(gplots)
  library(RColorBrewer)
  
  col1<-colorRampPalette(brewer.pal(8,"Spectral"), bias=1)
  png(filename="Number_DE_genes_Wilcox_number.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(Wilcox_number), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",labRow = c("Rect","Rect\nLoad","Circ","Circ\nLoad"),labCol=c("Rect","Rect\n Load","Circ","Circ\n Load"),density.info="none"
            ,key.par=list(mar=c(4,1,4,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2)
  dev.off()
  
  png(filename="Number_DE_genes_Wilcox_number_numered.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(Wilcox_number), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",labRow = c("Rect","Rect\nLoad","Circ","Circ\nLoad"),labCol=c("Rect","Rect\n Load","Circ","Circ\n Load"),density.info="none"
            ,key.par=list(mar=c(4,1,4,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2,cellnote=Wilcox_number, notecex=1.0,notecol="black",
            na.color=par("bg"))
  dev.off()
  
}

#save the genelists
{
  S1S2<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S1S2_gene)
  S1S3<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S1S3_gene)
  S1S4<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S1S4_gene)
  S2S3<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S2S3_gene)
  S2S4<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S2S4_gene)
  S3S4<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% S3S4_gene)
  significant_gene_data<-rbind(S1S2,S1S3)
  significant_gene_data<-rbind(significant_gene_data,S1S4)
  significant_gene_data<-rbind(significant_gene_data,S2S3)
  significant_gene_data<-rbind(significant_gene_data,S3S4)
  write.csv(significant_gene_data,file="Significant_genes_data.csv")
  
  genelist <- list(S1S2$SYMBOL,S1S3$SYMBOL, S1S4$SYMBOL,S2S3$SYMBOL, S2S4$SYMBOL, S3S4$SYMBOL )
  ngene <- sapply(genelist, length)
  maxg <- seq_len(max(ngene))
  significant_gene <- as.data.frame(sapply(genelist, "[", i = maxg))
  names(significant_gene)<-c("S1S2","S1S3","S1S4","S2S3","S2S4","S3S4")
  write.csv(significant_gene, file="significant_genes.csv")
}

#Fold change box plot
{
  cd<-c("RL/R","C/R","CL/R","C/RL","CL/RL","CL/C")
  x<-which(names(rnaseq_mapped_genes) %in% c("FC_S2_S1","FC_S3_S1","FC_S4_S1","FC_S3_S2","FC_S4_S2","FC_S4_S3"))
  png(filename="Foldchange_transcript_S1S2.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  boxplot(log(S1S2[,x], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18,names=cd)
  dev.off()
  png(filename="Foldchange_transcript_S1S3.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  boxplot(log(S1S3[,x], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18,names=cd)
  dev.off()
  png(filename="Foldchange_transcript_S1_S4.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  boxplot(log(S1S4[,x], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18,names=cd)
  dev.off()
  png(filename="Foldchange_transcript_S2S3.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  boxplot(log(S2S3[,x], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18,names=cd)
  dev.off()
  png(filename="Foldchange_transcript_S2S4.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  boxplot(log(S2S4[,x], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18,names=cd)
  dev.off()
  png(filename="Foldchange_transcript_S3S4.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  boxplot(log(S3S4[,x], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18,names=cd)
  dev.off()
  
  
}

# Fold change hist plot
{
  png(filename="Foldchange_transcript_S1S2_hist.png", units="in", width=1.5, height=1.5, pointsize=5, res=1200)
  par(font.axis=2,font.lab=2,font=2)
  hist(log(S1S2$FC_S2_S1, base=2), col="gray", las=1, cex.axis=1,xlab="FC(RL/R) sig_genes",main="", breaks=40)
  box()
  dev.off()
  
  png(filename="Foldchange_transcript_S1S3_hist.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  hist(log(S1S3$FC_S3_S1, base=2), col="gray", las=1, cex.axis=0.8,xlab="FC(C/R) sig_genes",main="", breaks=30)
  box()
  dev.off()
  
  png(filename="Foldchange_transcript_S1S4_hist.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  hist(log(S1S4$FC_S4_S1, base=2), col="gray", las=1, cex.axis=0.8,xlab="FC(CL/R) sig_genes",main="", breaks=30)
  box()
  dev.off()
  
  png(filename="Foldchange_transcript_S3S2_hist.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  hist(log(S2S3$FC_S3_S2, base=2), col="gray", las=1, cex.axis=0.8,xlab="FC(C/RL) sig_genes",main="", breaks=30)
  box()
  dev.off()
  
  
  png(filename="Foldchange_transcript_S4S2_hist.png", units="in", width=3, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  hist(log(S2S4$FC_S4_S2, base=2), col="gray", las=1, cex.axis=0.8,xlab="FC(CL/RL) sig_genes",main="", breaks=30)
  box()
  dev.off()
  
  
  png(filename="Foldchange_transcript_S4S3_hist.png", units="in", width=1.5, height=1.5, pointsize=5, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  hist(log(S3S4$FC_S4_S3, base=2), col="gray", las=1, cex.axis=1,xlab="FC(CL/C) sig_genes",main="", breaks=30)
  box()
  dev.off()
}

#Sig plot expression patterns 
{
  png(filename="Sig_S1S2.png", units="in", width=6, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, mfrow=c(1,3))
  plot(log2(S1S2$FC_S4_S3)~log2(S1S2$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="CL/C",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S1S2$FC_S3_S1)~log2(S1S2$FC_S4_S2), pch=18, las=1, xlab="CL/RL",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S1S2$FC_S3_S1)~log2(S1S2$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="Sig_S1S3.png", units="in", width=6, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, mfrow=c(1,3))
  plot(log2(S1S3$FC_S4_S3)~log2(S1S3$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="CL/C",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S1S3$FC_S3_S1)~log2(S1S3$FC_S4_S2), pch=18, las=1, xlab="CL/RL",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S1S3$FC_S3_S1)~log2(S1S3$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="Sig_S1S4.png", units="in", width=6, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, mfrow=c(1,3))
  plot(log2(S1S4$FC_S4_S3)~log2(S1S4$FC_S2_S1), pch=18, las=1,xlab="RL/R",ylab="CL/C",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S1S4$FC_S3_S1)~log2(S1S4$FC_S4_S2), pch=18, las=1, xlab="CL/RL",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S1S4$FC_S3_S1)~log2(S1S4$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="Sig_S2S3.png", units="in", width=6, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, mfrow=c(1,3))
  plot(log2(S2S3$FC_S4_S3)~log2(S2S3$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="CL/C",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S2S3$FC_S3_S1)~log2(S2S3$FC_S4_S2), pch=18, las=1, xlab="CL/RL",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S2S3$FC_S3_S1)~log2(S2S3$FC_S2_S1), pch=18, las=1,xlab="RL/R",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="Sig_S2S4.png", units="in", width=6, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, mfrow=c(1,3))
  plot(log2(S2S4$FC_S4_S3)~log2(S2S4$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="CL/C",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S2S4$FC_S3_S1)~log2(S2S4$FC_S4_S2), pch=18, las=1, xlab="CL/RL",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S2S4$FC_S3_S1)~log2(S2S4$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="Sig_S3S4.png", units="in", width=6, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, mfrow=c(1,3))
  plot(log2(S3S4$FC_S4_S3)~log2(S3S4$FC_S2_S1), pch=18, las=1, xlab="RL/R",ylab="CL/C",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S3S4$FC_S3_S1)~log2(S3S4$FC_S4_S2), pch=18, las=1,xlab="CL/RL",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  plot(log2(S3S4$FC_S3_S1)~log2(S3S4$FC_S2_S1), pch=18, las=1,xlab="RL/R",ylab="C/R",ylim=c(-5,5),xlim=c(-5,5))
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
}

#Sig genes heatmaps
{
  library(gplots)
  x<-which(names(S1S2)%in%c("S1","S2","S3","S4"))
  x<-which(names(S1S2)%in%c("S1","S2","S3","S4"))
  S3S41<-S3S4[is.finite(S3S4$FC_S4_S3),]
  S3S41<-S3S41[order(S3S41$S4_zscore, decreasing = T), ]
  n<-length(levels(as.factor(S3S41$GENEBIOTYPE)))
  col1<-colorRampPalette(brewer.pal(11,"Paired"))
  dd<-as.numeric(as.factor(S3S41$GENEBIOTYPE))
  rc<-col1(n)[dd]
  png(filename="S3S4_Significant_genes.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  heatmap.2(as.matrix(S3S41[,x]), scale="row",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
            cexRow=0.35,cex.axis=0.5,cexCol=0.8,srtCol =0,ylab = "",
            labRow = S3S41$SYMBOL,labCol=c("R"," RL"," C","  CL"),
            key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
            margins = c(2,14),col=redgreen(100), RowSideColors = rc)
  
  par(lend = 1)           # square line ends for the color legend
  legend("bottomright",      # location of the legend on the heatmap plot
         legend = levels(as.factor(S3S41$GENEBIOTYPE)), # category labels
         col = col1(n),  # color key
         lty= 1,             # line style
         lwd = 4,cex=0.3            # line width
  )
  dev.off()
  
  
  x<-which(names(S1S2)%in%c("S1","S2","S3","S4"))
  S1S21<-S1S2[is.finite(S1S2$FC_S2_S1),]
  S1S21<-S1S21[order(S1S21$S2_zscore, decreasing = T), ]
  n<-length(levels(as.factor(S1S21$GENEBIOTYPE)))
  col1<-colorRampPalette(brewer.pal(11,"Paired"))
  dd<-as.numeric(as.factor(S1S21$GENEBIOTYPE))
  rc<-col1(n)[dd]
  png(filename="S1S2_Significant_genes.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  heatmap.2(as.matrix(S1S21[,x]), scale="row",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
            cexRow=0.35,cex.axis=0.5,cexCol=0.8,srtCol =0,ylab = "",
            labRow = S1S21$SYMBOL,labCol=c("R"," RL"," C","  CL"),
            key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
            margins = c(2,14),col=redgreen(100), RowSideColors = rc)
  
  par(lend = 1)           # square line ends for the color legend
  legend("bottomright",      # location of the legend on the heatmap plot
         legend = levels(as.factor(S1S21$GENEBIOTYPE)), # category labels
         col = col1(n),  # color key
         lty= 1,             # line style
         lwd = 4,cex=0.3            # line width
  )
  dev.off()
  
  
  x<-which(names(S1S2)%in%c("S1","S2","S3","S4"))
  S1S31<-S1S3[is.finite(S1S3$FC_S3_S1),]
  S1S31<-S1S31[order(S1S31$S3_zscore, decreasing = T), ]
  n<-length(levels(as.factor(S1S31$GENEBIOTYPE)))
  col1<-colorRampPalette(brewer.pal(11,"Paired"))
  dd<-as.numeric(as.factor(S1S31$GENEBIOTYPE))
  rc<-col1(n)[dd]
  png(filename="S1S3_Significant_genes.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  heatmap.2(as.matrix(S1S31[,x]), scale="row",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
            cexRow=0.35,cex.axis=0.5,cexCol=0.8,srtCol =0,ylab = "",
            labRow = S1S31$SYMBOL,labCol=c("R"," RL"," C","  CL"),
            key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
            margins = c(2,14),col=redgreen(100), RowSideColors = rc)
  
  par(lend = 1)           # square line ends for the color legend
  legend("bottomright",      # location of the legend on the heatmap plot
         legend = levels(as.factor(S1S31$GENEBIOTYPE)), # category labels
         col = col1(n),  # color key
         lty= 1,             # line style
         lwd = 4,cex=0.3            # line width
  )
  dev.off()
  
  x<-which(names(S1S4)%in%c("S1","S2","S3","S4"))
  S1S41<-S1S4[is.finite(S1S4$FC_S4_S1),]
  S1S41<-S1S41[order(S1S41$S4_zscore, decreasing = T), ]
  n<-length(levels(as.factor(S1S41$GENEBIOTYPE)))
  col1<-colorRampPalette(brewer.pal(11,"Paired"))
  dd<-as.numeric(as.factor(S1S41$GENEBIOTYPE))
  rc<-col1(n)[dd]
  png(filename="S1S4_Significant_genes.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  heatmap.2(as.matrix(S1S41[,x]), scale="row",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
            cexRow=0.35,cex.axis=0.5,cexCol=0.8,srtCol =0,ylab = "",
            labRow = S1S41$SYMBOL,labCol=c("R"," RL"," C","  CL"),
            key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
            margins = c(2,14),col=redgreen(100), RowSideColors = rc)
  
  par(lend = 1)           # square line ends for the color legend
  legend("bottomright",      # location of the legend on the heatmap plot
         legend = levels(as.factor(S1S41$GENEBIOTYPE)), # category labels
         col = col1(n),  # color key
         lty= 1,             # line style
         lwd = 4,cex=0.3            # line width
  )
  dev.off()
  
  x<-which(names(S2S3)%in%c("S1","S2","S3","S4"))
  S2S31<-S2S3[is.finite(S2S3$FC_S3_S2),]
  S2S31<-S2S31[order(S2S31$S3_zscore, decreasing = T), ]
  n<-length(levels(as.factor(S2S31$GENEBIOTYPE)))
  col1<-colorRampPalette(brewer.pal(11,"Paired"))
  dd<-as.numeric(as.factor(S2S31$GENEBIOTYPE))
  rc<-col1(n)[dd]
  png(filename="S2S3_Significant_genes.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  heatmap.2(as.matrix(S2S31[,x]), scale="row",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
            cexRow=0.35,cex.axis=0.5,cexCol=0.8,srtCol =0,ylab = "",
            labRow = S2S31$SYMBOL,labCol=c("R"," RL"," C","  CL"),
            key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
            margins = c(2,14),col=redgreen(100), RowSideColors = rc)
  
  par(lend = 1)           # square line ends for the color legend
  legend("bottomright",      # location of the legend on the heatmap plot
         legend = levels(as.factor(S2S31$GENEBIOTYPE)), # category labels
         col = col1(n),  # color key
         lty= 1,             # line style
         lwd = 4,cex=0.3            # line width
  )
  dev.off()
  
  
  x<-which(names(S2S4)%in%c("S1","S2","S3","S4"))
  S2S41<-S2S4[is.finite(S2S4$FC_S4_S2),]
  S2S41<-S2S41[order(S2S41$S4_zscore, decreasing = T), ]
  n<-length(levels(as.factor(S2S41$GENEBIOTYPE)))
  col1<-colorRampPalette(brewer.pal(11,"Paired"))
  dd<-as.numeric(as.factor(S2S41$GENEBIOTYPE))
  rc<-col1(n)[dd]
  png(filename="S2S4_Significant_genes.png", units="in", width=2, height=3.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  heatmap.2(as.matrix(S2S41[,x]), scale="row",Colv = NA, Rowv = NA,dendrogram="none",trace='none',
            cexRow=0.35,cex.axis=0.5,cexCol=0.8,srtCol =0,ylab = "",
            labRow = S2S41$SYMBOL,labCol=c("R"," RL"," C","  CL"),
            key=T, symkey=T,key.title=" ",key.xlab = "", density.info="none", keysize=0.5,key.par=list(mar=c(1,1,5,0.5),cex=0.4, font=2),
            margins = c(2,14),col=redgreen(100), RowSideColors = rc)
  
  par(lend = 1)           # square line ends for the color legend
  legend("bottomright",      # location of the legend on the heatmap plot
         legend = levels(as.factor(S2S41$GENEBIOTYPE)), # category labels
         col = col1(n),  # color key
         lty= 1,             # line style
         lwd = 4,cex=0.3            # line width
  )
  dev.off()
  
  
}



#enrichment function
library(gskb)
data(mm_pathway)
data(mm_GO)
data(mm_TF)
data(mm_miRNA)
data(mm_metabolic)
data(mm_location)
data(mm_other)

enrichment<-function(expression_data,de_geneid,datasets,type_name, sampledirection){
  test<-subset(expression_data, expression_data$SYMBOL %in% de_geneid)
  test$SYMBOL<-toupper(test$SYMBOL)
  x<-which(names(test) %in% c("SYMBOL","S1","S2","S3","S4"))
  test<-test[,x]
  test<-subset(test,!duplicated(test$SYMBOL))
  
  a_c<-length(de_geneid)
  b_d<-nrow(expression_data)-length(de_geneid)
  datasets<-datasets
  n<-length(datasets)
  dataset_test<-as.data.frame(matrix(nrow=n,ncol=14))
  colnames(dataset_test)<-c(type_name,"pvalue","fdr","overlap_percent","overlap_genes","odds_ratio","no_de_genes_in_path",
                            "no_de_genes","no_genes_path","estimate_oddsratio","hypergeo_pvalue","Expected","Enrichment",
                            "Type")
  
  if(nrow(test)>0){
    for(i in 1:n){
      temp<-na.omit(subset(test,test$SYMBOL %in% (datasets[[i]][3:length(datasets[[i]])])))
      dataset_test[i,1]<-datasets[[i]][1]
      a <- nrow(temp)
      b <- (length(datasets[[i]])-2)-a
      c <- a_c - a
      d <- b_d - b
      dataset_test$odds_ratio[i]=((a/b)/(c/d))
      dataset_test$no_de_genes_in_path[i]=a
      dataset_test$no_de_genes[i]=a_c
      dataset_test$no_genes_path[i]=b
      dataset_test$overlap_percent[i]=round((a/(a+b) )* 100, digit=2)
      dataset_test$overlap_genes[i]= paste(temp$SYMBOL, collapse=",")
      fisherRes<-fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2),alternative="greater")
      dataset_test$pvalue[i]=fisherRes$p.value
      dataset_test$estimate_oddsratio[i]=fisherRes$estimate
      dataset_test$hypergeo_pvalue[i]=phyper(a,(a+b),(c+d),(a+c),lower.tail=T,log.p = FALSE)
      dataset_test$Expected[i]=((a+b)/(a+b+c+d))*(a+c)
      dataset_test$Enrichment[i]=a/ dataset_test$Expected[i]
      dataset_test$fdr[i]=0
    }
  }
  dataset_test<-subset(dataset_test,dataset_test$no_de_genes_in_path>0)
  
  if(nrow(dataset_test)>0){
    # Multiple-test correction
    dataset_test$fdr <- p.adjust(dataset_test$pvalue, method="BH")
    dataset_test <- dataset_test[order(dataset_test$pvalue, decreasing=F),]
    dataset_test$Type<-sampledirection
    
  }
  
  return(dataset_test)
  
}  
#pathways
a<-mm_pathway
c<-"pathway"
f1<-"pathway_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}

#GO
a<-mm_GO
c<-"GO"
f1<-"GO_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}

#TF
a<-mm_TF
c<-"TF"
f1<-"TF_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}

##miRNA
a<-mm_miRNA
c<-"miRNA"
f1<-"miRNA_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}

##metabolic
a<-mm_metabolic
c<-"metabolic"
f1<-"metabolic_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}

#location
a<-mm_location
c<-"location"
f1<-"location_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}

#other
a<-mm_other
c<-"other"
f1<-"other_types.csv"
{
  S1S2_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1>1)$SYMBOL, a,c,"S1S2up")
  S1S2_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S2,S1S2$FC_S2_S1<1)$SYMBOL, a,c,"S1S2down")
  S1S2_path<-rbind(S1S2_path_up,S1S2_path_down)
  
  S3S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3>1)$SYMBOL, a,c,"S3S4up")
  S3S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S3S4,S3S4$FC_S4_S3<1)$SYMBOL, a,c,"S3S4down")
  S3S4_path<-rbind(S3S4_path_up,S3S4_path_down)
  
  S1S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1>1)$SYMBOL, a,c,"S1S3up")
  S1S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S3,S1S3$FC_S3_S1<1)$SYMBOL, a,c,"S1S3down")
  S1S3_path<-rbind(S1S3_path_up,S1S3_path_down)
  
  S1S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1>1)$SYMBOL,  a,c,"S1S4up")
  S1S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S1S4,S1S4$FC_S4_S1<1)$SYMBOL,  a,c,"S1S4down")
  S1S4_path<-rbind(S1S4_path_up,S1S4_path_down)
  
  S2S3_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2>1)$SYMBOL,  a,c,"S2S3up")
  S2S3_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S3,S2S3$FC_S3_S2<1)$SYMBOL,  a,c,"S2S3down")
  S2S3_path<-rbind(S2S3_path_up,S2S3_path_down)
  
  S2S4_path_up<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2>1)$SYMBOL,a,c,"S2S4up")
  S2S4_path_down<-enrichment(rnaseq_mapped_genes,subset(S2S4,S2S4$FC_S4_S2<1)$SYMBOL, a,c,"S2S4down")
  S2S4_path<-rbind(S2S4_path_up,S2S4_path_down)
  
  pathway_types<-rbind(S1S2_path,S1S3_path)
  pathway_types<-rbind(pathway_types,S1S4_path)
  pathway_types<-rbind(pathway_types,S2S3_path)
  pathway_types<-rbind(pathway_types,S2S4_path)
  pathway_types<-rbind(pathway_types,S3S4_path)
  write.csv(pathway_types,file=f1)
}



dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part4_FC_SIG/expression_patterns_fc_and_sig"
dir.create(dirp)
setwd(dirp)

#Expression patterns gene sets
{
  diff_in_circ_and_rect<-unique(c(subset(S1S2$SYMBOL, (S1S2$SYMBOL %in% S3S4$SYMBOL)),subset(S3S4$SYMBOL, (S3S4$SYMBOL %in% S1S2$SYMBOL))))
  diff_in_rect_not_circ<-subset(S1S2$SYMBOL, !(S1S2$SYMBOL %in% S3S4$SYMBOL))
  diff_in_circ_not_rect<-subset(S3S4$SYMBOL, !(S3S4$SYMBOL %in% S1S2$SYMBOL))
  
  diff_in_before_and_after<-unique(c(subset(S2S4$SYMBOL, (S2S4$SYMBOL %in% S1S3$SYMBOL)),subset(S1S3$SYMBOL, (S1S3$SYMBOL %in% S2S4$SYMBOL)) ))
  diff_in_after_not_before<-subset(S2S4$SYMBOL, !(S2S4$SYMBOL %in% S1S3$SYMBOL))
  diff_in_before_not_after<-subset(S1S3$SYMBOL, !(S1S3$SYMBOL %in% S2S4$SYMBOL))
  
  diff_in_circ_and_rectload<-unique(c(subset(S1S2$SYMBOL, !(S1S2$SYMBOL %in% S1S3$SYMBOL)), subset(S1S3$SYMBOL, !(S1S3$SYMBOL %in% S1S2$SYMBOL))))
  diff_in_rectload_not_circ<-subset(S1S2$SYMBOL, !(S1S2$SYMBOL %in% S1S3$SYMBOL))
  diff_in_circ_not_rectload<-subset(S1S3$SYMBOL, !(S1S3$SYMBOL %in% S1S2$SYMBOL))
  
  genelist <- list(diff_in_circ_and_rect,diff_in_rect_not_circ, diff_in_circ_not_rect
                   ,diff_in_before_and_after, diff_in_after_not_before, diff_in_before_not_after,
                   diff_in_circ_and_rectload,diff_in_rectload_not_circ,diff_in_circ_not_rectload)
  ngene <- sapply(genelist, length)
  maxg <- seq_len(max(ngene))
  significant_gene <- as.data.frame(sapply(genelist, "[", i = maxg))
  names(significant_gene)<-c("diff_in_circ_and_rect","diff_in_rect_not_circ", "diff_in_circ_not_rect"
                             ,"diff_in_before_and_after", "diff_in_after_not_before", "diff_in_before_not_after",
                             "diff_in_circ_and_rectload","diff_in_rectload_not_circ","diff_in_circ_not_rectload")
  write.csv(significant_gene, file="significant_genes_types_symbols.csv")
  
  diff_in_circ_and_rect<-unique(c(subset(S1S2$SYMBOL, (S1S2$SYMBOL %in% S3S4$SYMBOL)),subset(S3S4$SYMBOL, (S3S4$SYMBOL %in% S1S2$SYMBOL))))
  diff_in_rect_not_circ<-subset(S1S2$SYMBOL, !(S1S2$SYMBOL %in% S3S4$SYMBOL))
  diff_in_circ_not_rect<-subset(S3S4$SYMBOL, !(S3S4$SYMBOL %in% S1S2$SYMBOL))
  
  diff_in_before_and_after<-unique(c(subset(S2S4$SYMBOL, (S2S4$SYMBOL %in% S1S3$SYMBOL)), subset(S1S3$SYMBOL, (S1S3$SYMBOL %in% S2S4$SYMBOL))))
  diff_in_after_not_before<-subset(S2S4$SYMBOL, !(S2S4$SYMBOL %in% S1S3$SYMBOL))
  diff_in_before_not_after<-subset(S1S3$SYMBOL, !(S1S3$SYMBOL %in% S2S4$SYMBOL))
  
  diff_in_circ_and_rectload<-unique(c(subset(S1S2$SYMBOL, (S1S2$SYMBOL %in% S1S3$SYMBOL)), subset(S1S3$SYMBOL, (S1S3$SYMBOL %in% S1S2$SYMBOL)) ))
  diff_in_rectload_not_circ<-subset(S1S2$SYMBOL, !(S1S2$SYMBOL %in% S1S3$SYMBOL))
  diff_in_circ_not_rectload<-subset(S1S3$SYMBOL, !(S1S3$SYMBOL %in% S1S2$SYMBOL))
  
  genelist <- list(diff_in_circ_and_rect,diff_in_rect_not_circ, diff_in_circ_not_rect
                   ,diff_in_before_and_after, diff_in_after_not_before, diff_in_before_not_after,
                   diff_in_circ_and_rectload,diff_in_rectload_not_circ,diff_in_circ_not_rectload)
  ngene <- sapply(genelist, length)
  maxg <- seq_len(max(ngene))
  significant_gene <- as.data.frame(sapply(genelist, "[", i = maxg))
  names(significant_gene)<-c("diff_in_circ_and_rect","diff_in_rect_not_circ", "diff_in_circ_not_rect"
                             ,"diff_in_before_and_after", "diff_in_after_not_before", "diff_in_before_not_after",
                             "diff_in_circ_and_rectload","diff_in_rectload_not_circ","diff_in_circ_not_rectload")
  write.csv(significant_gene, file="significant_genes_types_geneid.csv")
  
  diff_in_circ_and_rect<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_and_rect)
  diff_in_rect_not_circ<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_rect_not_circ)
  diff_in_circ_not_rect<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_not_rect)
  diff_in_before_and_after<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_before_and_after)
  diff_in_after_not_before<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_after_not_before)
  diff_in_before_not_after<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_before_not_after)
  diff_in_circ_and_rectload<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_and_rectload)
  diff_in_rectload_not_circ<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_rectload_not_circ)
  diff_in_circ_not_rectload<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_not_rectload)
  
  diff_in_circ_not_rect$Type<-"diff_in_circ_not_rec"
  diff_in_rect_not_circ$Type<-"diff_rec_not_circ"
  diff_in_circ_and_rect$Type<-"diff_in_circ_and_rect"
  diff_in_before_and_after$Type<-"diff_before_not_after"
  diff_in_after_not_before$Type<-"diff_after_not_before"
  diff_in_before_not_after$Type<-"diff_before_and_after"
  diff_in_circ_not_rectload$Type<-"diff_in_circ_not_rectload"
  diff_in_rectload_not_circ$Type<-"diff_rectload_not_circ"
  diff_in_circ_and_rectload$Type<-"diff_rectload_and_circ"
  
  different_types_sig_genes<-rbind(diff_in_circ_not_rect,diff_in_rect_not_circ)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_circ_and_rect)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_before_and_after)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_after_not_before)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_before_not_after)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_circ_not_rectload)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_rectload_not_circ)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_circ_and_rectload)
  write.csv(different_types_sig_genes,"different_types_sig_genes_details.csv")
  
}

#VennDiagram
{
  library(VennDiagram)
  png(filename="Venn_unique_to_geomtry.png", units="in", width=2, height=2 , pointsize=4, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  b<-nrow(diff_in_circ_not_rect)+nrow(diff_in_circ_and_rect)
  a<-nrow(diff_in_rect_not_circ)+nrow(diff_in_circ_and_rect)
  c<-nrow(diff_in_circ_and_rect)
  draw.pairwise.venn(a, b, c, category = c("Changing in Rectangle after load", "Changing in Circle after load"), lty = rep("blank",2), 
                     fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
  dev.off()
  
  png(filename="Venn_unique_to_geomtry_load.png", units="in", width=2, height=2 , pointsize=4, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  a<-nrow(diff_in_after_not_before)+nrow(diff_in_before_and_after)
  b<-nrow(diff_in_before_not_after)+nrow(diff_in_before_and_after)
  c<-nrow(diff_in_before_and_after)
  draw.pairwise.venn(a, b, c, category = c("Different after load", "Different before load"), lty = rep("blank",2), 
                     fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2))
  dev.off()
  
  png(filename="Venn_unique_to_shavevsload.png", units="in", width=2, height=2 , pointsize=4, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  a<-nrow(diff_in_circ_not_rectload)+nrow(diff_in_circ_and_rectload)
  b<-nrow(diff_in_rectload_not_circ)+nrow(diff_in_circ_and_rectload)
  c<-nrow(diff_in_circ_and_rectload)
  draw.pairwise.venn(a, b, c, category = c("Different with shape", "Different with load(rect)"), lty = rep("blank",2), 
                     fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(-10,10), cat.dist = rep(0.025, 2))
  dev.off()
}

#number_hetmaps
{
  number<-matrix(nrow=4,ncol=3)
  colnames(number)<-c("up","down","total")
  rownames(number)<-c("diff_in_circ_not_rect","diff_in_rect_not_circ","different_in_circ_and_rec_circ","different_in_circ_and_rec_rec")
  
  number[1,1]<-nrow(subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1))
  number[1,2]<-nrow(subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1))
  number[1,3]<-nrow(diff_in_circ_not_rect)
  number[2,1]<-nrow(subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1))
  number[2,2]<-nrow(subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1))
  number[2,3]<-nrow(diff_in_rect_not_circ)
  number[3,1]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S4_S3>1))
  number[3,2]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S4_S3<1))
  number[3,3]<-nrow(diff_in_circ_and_rect)
  number[4,1]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1))
  number[4,2]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1))
  number[4,3]<-nrow(diff_in_circ_and_rect)
  circ_rect<-number
  
  col1<-colorRampPalette(brewer.pal(8,"RdPu"), bias=1)
  png(filename="circ_rect.png", units="in", width=1.3, height=1 , pointsize=3, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(circ_rect), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",density.info="none",cellnote=circ_rect, notecex=1.0,notecol="black",
            na.color=par("bg") ,key.par=list(mar=c(4,1,2,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2, margins=c(4,15))
  dev.off()
  
  
  
  
  
}

{
  number<-matrix(nrow=4,ncol=3)
  colnames(number)<-c("up","down","total")
  rownames(number)<-c("diff_in_after_not_before","diff_in_before_not_after","diff_in_before_and_after_before","diff_in_before_and_after_after")
  
  number[1,1]<-nrow(subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1))
  number[1,2]<-nrow(subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1))
  number[1,3]<-nrow(diff_in_after_not_before)
  number[2,1]<-nrow(subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1))
  number[2,2]<-nrow(subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1))
  number[2,3]<-nrow(diff_in_before_not_after)
  number[3,1]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1))
  number[3,2]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1))
  number[3,3]<-nrow(diff_in_before_and_after)
  number[4,1]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S4_S2>1))
  number[4,2]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S4_S2<1))
  number[4,3]<-nrow(diff_in_before_and_after)
  after_before<-number
  
  col1<-colorRampPalette(brewer.pal(8,"RdPu"), bias=1)
  png(filename="after_before.png", units="in", width=1.3, height=1 , pointsize=3, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(after_before), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",density.info="none",cellnote=circ_rect, notecex=1.0,notecol="black",
            na.color=par("bg") ,key.par=list(mar=c(4,1,2,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2, margins=c(4,15))
  dev.off()
}

{
  number<-matrix(nrow=4,ncol=3)
  colnames(number)<-c("up","down","total")
  rownames(number)<-c("diff_in_circ_not_rectload","diff_in_rectload_not_circ","diff_in_circ_and_rectload_circ","diff_in_circ_and_rectload_rectl")
  
  number[1,1]<-nrow(subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1))
  number[1,2]<-nrow(subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1))
  number[1,3]<-nrow(diff_in_circ_not_rectload)
  number[2,1]<-nrow(subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1))
  number[2,2]<-nrow(subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1))
  number[2,3]<-nrow(diff_in_rectload_not_circ)
  number[3,1]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1))
  number[3,2]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1))
  number[3,3]<-nrow(diff_in_circ_and_rectload)
  number[4,1]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S2_S1>1))
  number[4,2]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S2_S1<1))
  number[4,3]<-nrow(diff_in_circ_and_rectload)
  load_circ<-number
  
  col1<-colorRampPalette(brewer.pal(8,"RdPu"), bias=1)
  png(filename="load_circ.png", units="in", width=1.3, height=1 , pointsize=3, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(load_circ), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",density.info="none",cellnote=circ_rect, notecex=1.0,notecol="black",
            na.color=par("bg") ,key.par=list(mar=c(4,1,2,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2, margins=c(4,15))
  dev.off()
}

#Expressionpattern plots
{ 
  png(filename="Circ_rect_load_sig.png", units="in", width=2, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  plot(log2(diff_in_circ_not_rect$FC_S4_S3)~log2(diff_in_circ_not_rect$FC_S2_S1), pch=18, las=1, xlab="Rect+load/Rect",
       ylab="Circ+load/Circ",ylim=c(-5,5),xlim=c(-5,5), col=adjustcolor("purple", alpha.f = 0.4), cex=0.7)
  points(log2(diff_in_rect_not_circ$FC_S4_S3)~log2(diff_in_rect_not_circ$FC_S2_S1),pch=18, col=adjustcolor("green4", alpha.f = 0.4), cex=0.7)
  points(log2(diff_in_circ_and_rect$FC_S4_S3)~log2(diff_in_circ_and_rect$FC_S2_S1),pch=18, col=adjustcolor("red", alpha.f = 0.4), cex=0.7)
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="before_after_sig.png", units="in", width=2, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  plot(log2(diff_in_after_not_before$FC_S3_S1)~log2(diff_in_after_not_before$FC_S4_S2), pch=18, las=1, xlab="Circ+load/Rect+load",
       ylab="Rect/Circ",ylim=c(-5,5),xlim=c(-5,5), col=adjustcolor("purple", alpha.f = 0.4), cex=0.7)
  points(log2(diff_in_before_not_after$FC_S3_S1)~log2(diff_in_before_not_after$FC_S4_S2),pch=18, col=adjustcolor("green4", alpha.f = 0.4), cex=0.7)
  points(log2(diff_in_before_and_after$FC_S3_S1)~log2(diff_in_before_and_after$FC_S4_S2),pch=18, col=adjustcolor("red", alpha.f = 0.4), cex=0.7)
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="geo_and_load_sig.png", units="in", width=2, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  plot(log2(diff_in_circ_not_rectload$FC_S3_S1)~log2(diff_in_circ_not_rectload$FC_S2_S1), pch=18, las=1, xlab="Rect+load/Rect",
       ylab="Rect/Circ",ylim=c(-5,5),xlim=c(-5,5), col=adjustcolor("purple", alpha.f = 0.4), cex=0.7)
  points(log2(diff_in_rectload_not_circ$FC_S3_S1)~log2(diff_in_rectload_not_circ$FC_S2_S1),pch=18, col=adjustcolor("green4", alpha.f = 0.4), cex=0.7)
  points(log2(diff_in_circ_and_rectload$FC_S3_S1)~log2(diff_in_circ_and_rectload$FC_S2_S1),pch=18, col=adjustcolor("red", alpha.f = 0.4), cex=0.7)
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
}

#pathways
a<-mm_pathway
c<-"pathway"
f1<-"pathway_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}

#TF
a<-mm_TF
c<-"TF"
f1<-"TF_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}

#GO
a<-mm_GO
c<-"GO"
f1<-"GO_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}

#miRNA
a<-mm_miRNA
c<-"miRNA"
f1<-"miRNA_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


#metabolic
a<-mm_metabolic
c<-"metabolic"
f1<-"metabolic_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


#location
a<-mm_location
c<-"location"
f1<-"location_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


#other
a<-mm_other
c<-"other"
f1<-"other_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part4_FC_SIG/expression_patterns_fc_and_sig_with_strict_FC"
dir.create(dirp)
setwd(dirp)

#Expression patterns gene sets
{
  diff_in_circ_not_rect<-subset(S3S4, (log2(S3S4$FC_S4_S3) >=2 | log2(S3S4$FC_S4_S3) <=-2) & S3S4$FC_S4_S3!=0 &
                                  (log2(S3S4$FC_S2_S1) < 2 & log2(S3S4$FC_S2_S1) >-2))
  diff_in_circ_not_rect<-subset(diff_in_circ_not_rect$SYMBOL, !(diff_in_circ_not_rect$SYMBOL %in% S1S2$SYMBOL))
  diff_in_rect_not_circ<-subset(S1S2, (log2(S1S2$FC_S2_S1) >=2 | log2(S1S2$FC_S2_S1) <=-2) & S1S2$FC_S2_S1!=0 &
                                  (log2(S1S2$FC_S4_S3) < 2 & log2(S1S2$FC_S4_S3) >-2))
  diff_in_rect_not_circ<-subset(diff_in_rect_not_circ$SYMBOL, !(diff_in_rect_not_circ$SYMBOL %in% S3S4$SYMBOL))
  diff_in_circ_and_rect_circ<-(subset(S3S4, (log2(S3S4$FC_S4_S3) >=2 | log2(S3S4$FC_S4_S3) <=-2) & S3S4$FC_S4_S3!=0 &
                                        (log2(S3S4$FC_S2_S1) >=2 | log2(S3S4$FC_S2_S1) <=-2)))
  diff_in_circ_and_rect_rect<-(subset(S1S2, (log2(S1S2$FC_S4_S3) >=2 | log2(S1S2$FC_S4_S3) <=-2) & S1S2$FC_S2_S1!=0 &
                                        (log2(S1S2$FC_S2_S1) >=2 | log2(S1S2$FC_S2_S1) <=-2)))
  diff_in_circ_and_rect<-unique(c(subset(diff_in_circ_and_rect_rect$SYMBOL, (diff_in_circ_and_rect_rect$SYMBOL %in% S3S4$SYMBOL)),
                                  subset(diff_in_circ_and_rect_circ$SYMBOL, (diff_in_circ_and_rect_circ$SYMBOL %in% S1S2$SYMBOL))))
  
  diff_in_before_not_after<-subset(S1S3, (log2(S1S3$FC_S3_S1) >=2 | log2(S1S3$FC_S3_S1) <=-2) & S1S3$FC_S3_S1 !=0 &
                                     (log2(S1S3$FC_S4_S2) < 2 & log2(S1S3$FC_S4_S2) >-2))
  diff_in_before_not_after<-subset(diff_in_before_not_after$SYMBOL, !(diff_in_before_not_after$SYMBOL %in% S2S4$SYMBOL))
  diff_in_after_not_before<-subset(S2S4, (log2(S2S4$FC_S4_S2) >=2 | log2(S2S4$FC_S4_S2) <=-2) & S2S4$FC_S4_S2 !=0 &
                                     (log2(S2S4$FC_S3_S1) < 2 & log2(S2S4$FC_S3_S1) >-2))
  diff_in_after_not_before<-subset(diff_in_after_not_before$SYMBOL, !(diff_in_after_not_before$SYMBOL %in% S1S3$SYMBOL))
  diff_in_before_and_after_before<-(subset(S1S3, (log2(S1S3$FC_S3_S1) >=2 | log2(S1S3$FC_S3_S1) <=-2) & S1S3$FC_S3_S1 !=0 &
                                             (log2(S1S3$FC_S4_S2) >=2 | log2(S1S3$FC_S4_S2) <=-2)))
  diff_in_before_and_after_after<-(subset(S2S4, (log2(S2S4$FC_S3_S1) >=2 | log2(S2S4$FC_S3_S1) <=-2) & S2S4$FC_S4_S2 !=0 &
                                            (log2(S2S4$FC_S4_S2) >=2 | log2(S2S4$FC_S4_S2) <=-2)))
  diff_in_before_and_after<-unique(c(subset(diff_in_before_and_after_after$SYMBOL, (diff_in_before_and_after_after$SYMBOL %in% S1S3$SYMBOL)),
                                     subset(diff_in_before_and_after_before$SYMBOL, (diff_in_before_and_after_before$SYMBOL %in% S2S4$SYMBOL))))

  diff_in_circ_not_rectload<-subset(S1S3, (log2(S1S3$FC_S3_S1) >=2 | log2(S1S3$FC_S3_S1) <=-2) & S1S3$FC_S3_S1 !=0 &
                                      (log2(S1S3$FC_S2_S1) < 2 & log2(S1S3$FC_S2_S1) >-2))
  diff_in_circ_not_rectload<-subset(diff_in_circ_not_rectload$SYMBOL, !(diff_in_circ_not_rectload$SYMBOL %in% S2S4$SYMBOL))
  diff_in_rectload_not_circ<-subset(S1S2, (log2(S1S2$FC_S2_S1) >=2 | log2(S1S2$FC_S2_S1) <= -2) & S1S2$FC_S2_S1 !=0 &
                                      (log2(S1S2$FC_S3_S1) < 2 & log2(S1S2$FC_S3_S1) > -2))
  diff_in_rectload_not_circ<-subset(diff_in_rectload_not_circ$SYMBOL, !(diff_in_rectload_not_circ$SYMBOL %in% S1S3$SYMBOL))
  diff_in_circ_and_rectload_circ<-(subset(S1S3, (log2(S1S3$FC_S3_S1) >=2 | log2(S1S3$FC_S3_S1) <=-2) & S1S3$FC_S3_S1 !=0 &
                                            (log2(S1S3$FC_S2_S1) >=2 | log2(S1S3$FC_S2_S1) <=-2)))
  diff_in_circ_and_rectload_rect<-(subset(S2S4, (log2(S2S4$FC_S3_S1) >=2 | log2(S2S4$FC_S3_S1) <=-2) & S2S4$FC_S2_S1 !=0 &
                                            (log2(S2S4$FC_S2_S1) >=2 | log2(S2S4$FC_S2_S1) <=-2)))
  diff_in_circ_and_rectload<-unique(c(subset(diff_in_circ_and_rectload_circ$SYMBOL, (diff_in_circ_and_rectload_circ$SYMBOL %in% S1S3$SYMBOL)),
                                      subset(diff_in_circ_and_rectload_rect$SYMBOL, (diff_in_circ_and_rectload_rect$SYMBOL %in% S2S4$SYMBOL))))
  
  
  genelist <- list(diff_in_circ_and_rect,diff_in_rect_not_circ, diff_in_circ_not_rect
                   ,diff_in_before_and_after, diff_in_after_not_before, diff_in_before_not_after,
                   diff_in_circ_and_rectload,diff_in_rectload_not_circ,diff_in_circ_not_rectload)
  ngene <- sapply(genelist, length)
  maxg <- seq_len(max(ngene))
  significant_gene <- as.data.frame(sapply(genelist, "[", i = maxg))
  names(significant_gene)<-c("diff_in_circ_and_rect","diff_in_rect_not_circ", "diff_in_circ_not_rect"
                             ,"diff_in_before_and_after", "diff_in_after_not_before", "diff_in_before_not_after",
                             "diff_in_circ_and_rectload","diff_in_rectload_not_circ","diff_in_circ_not_rectload")
  write.csv(significant_gene, file="significant_genes_types_symbols.csv")
  
  
  
  diff_in_circ_and_rect<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_and_rect)
  diff_in_rect_not_circ<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_rect_not_circ)
  diff_in_circ_not_rect<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_not_rect)
  diff_in_before_and_after<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_before_and_after)
  diff_in_after_not_before<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_after_not_before)
  diff_in_before_not_after<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_before_not_after)
  diff_in_circ_and_rectload<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_and_rectload)
  diff_in_rectload_not_circ<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_rectload_not_circ)
  diff_in_circ_not_rectload<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% diff_in_circ_not_rectload)
  
  diff_in_circ_not_rect$Type<-"diff_in_circ_not_rec"
  diff_in_rect_not_circ$Type<-"diff_rec_not_circ"
  diff_in_circ_and_rect$Type<-"diff_in_circ_and_rect"
  diff_in_before_and_after$Type<-"diff_before_and_after"
  diff_in_after_not_before$Type<-"diff_after_not_before"
  diff_in_before_not_after$Type<-"diff_before_not_after"
  diff_in_circ_not_rectload$Type<-"diff_in_circ_not_rectload"
  diff_in_rectload_not_circ$Type<-"diff_rectload_not_circ"
  diff_in_circ_and_rectload$Type<-"diff_rectload_and_circ"
  
  different_types_sig_genes<-rbind(diff_in_circ_not_rect,diff_in_rect_not_circ)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_circ_and_rect)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_before_and_after)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_after_not_before)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_before_not_after)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_circ_not_rectload)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_rectload_not_circ)
  different_types_sig_genes<-rbind(different_types_sig_genes,diff_in_circ_and_rectload)
  write.csv(different_types_sig_genes,"different_types_sig_genes_details.csv")
  
}

#VennDiagram
{
  library(VennDiagram)
  png(filename="Venn_unique_to_geomtry.png", units="in", width=1.5, height=1.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  b<-nrow(diff_in_circ_not_rect)+nrow(diff_in_circ_and_rect)
  a<-nrow(diff_in_rect_not_circ)+nrow(diff_in_circ_and_rect)
  c<-nrow(diff_in_circ_and_rect)
  draw.pairwise.venn(a, b, c, category = c("Rectangle Specific", "Circle Specific"), lty = rep(1,2), 
                     fill = c("green4", "purple"), alpha = c(0.5, 0.5), cat.pos = c(0,0), cat.dist = rep(0.025, 2),
                     col=c("green4", "purple"), lwd=0.5, ext.text=F, label.col=c("green4", "red","purple"))
  dev.off()
  
  png(filename="Venn_unique_to_geomtry_load.png", units="in", width=1.5, height=1.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  a<-nrow(diff_in_after_not_before)+nrow(diff_in_before_and_after)
  b<-nrow(diff_in_before_not_after)+nrow(diff_in_before_and_after)
  c<-nrow(diff_in_before_and_after)
  draw.pairwise.venn(a, b, c, category = c("Different with load", "Similar with load"), lty = rep(1,2), 
                     fill = c("purple", "green4"), alpha = c(0.5, 0.5), cat.pos = c(0,0), cat.dist = rep(0.025, 2),
                     col=c("purple", "green4"), lwd=0.5, ext.text=F, label.col=c("purple", "red","green4"))
  dev.off()
  
  png(filename="Venn_unique_to_shavevsload.png", units="in", width=1.5, height=1.5 , pointsize=6, res=1200)
  par(font.axis=2,font.lab=2, font=2)
  a<-nrow(diff_in_circ_not_rectload)+nrow(diff_in_circ_and_rectload)
  b<-nrow(diff_in_rectload_not_circ)+nrow(diff_in_circ_and_rectload)
  c<-nrow(diff_in_circ_and_rectload)
  draw.pairwise.venn(a, b, c, category = c("Different with shape", "Different with load"), lty = rep("blank",2), 
                     fill = c("purple", "green4"), alpha = c(0.5, 0.5), cat.pos = c(-20,20), cat.dist = rep(0.025, 2),
                     col=c("purple", "green4"), lwd=0.5, ext.text=F, label.col=c("purple", "red","green4"))
  dev.off()
}

#number_hetmaps
{
  number<-matrix(nrow=4,ncol=3)
  colnames(number)<-c("up","down","total")
  rownames(number)<-c("diff_in_circ_not_rect","diff_in_rect_not_circ","different_in_circ_and_rec_circ","different_in_circ_and_rec_rec")
  
  number[1,1]<-nrow(subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1))
  number[1,2]<-nrow(subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1))
  number[1,3]<-nrow(diff_in_circ_not_rect)
  number[2,1]<-nrow(subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1))
  number[2,2]<-nrow(subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1))
  number[2,3]<-nrow(diff_in_rect_not_circ)
  number[3,1]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S4_S3>1))
  number[3,2]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S4_S3<1))
  number[3,3]<-nrow(diff_in_circ_and_rect)
  number[4,1]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1))
  number[4,2]<-nrow(subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1))
  number[4,3]<-nrow(diff_in_circ_and_rect)
  circ_rect<-number
  
  col1<-colorRampPalette(brewer.pal(8,"RdPu"), bias=1)
  png(filename="circ_rect.png", units="in", width=1.3, height=1 , pointsize=3, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(circ_rect), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",density.info="none",cellnote=circ_rect, notecex=1.0,notecol="black",
            na.color=par("bg") ,key.par=list(mar=c(4,1,2,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2, margins=c(4,15))
  dev.off()
  
  
  
  
  
}

{
  number<-matrix(nrow=4,ncol=3)
  colnames(number)<-c("up","down","total")
  rownames(number)<-c("diff_in_after_not_before","diff_in_before_not_after","diff_in_before_and_after_before","diff_in_before_and_after_after")
  
  number[1,1]<-nrow(subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1))
  number[1,2]<-nrow(subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1))
  number[1,3]<-nrow(diff_in_after_not_before)
  number[2,1]<-nrow(subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1))
  number[2,2]<-nrow(subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1))
  number[2,3]<-nrow(diff_in_before_not_after)
  number[3,1]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1))
  number[3,2]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1))
  number[3,3]<-nrow(diff_in_before_and_after)
  number[4,1]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S4_S2>1))
  number[4,2]<-nrow(subset(diff_in_before_and_after,diff_in_before_and_after$FC_S4_S2<1))
  number[4,3]<-nrow(diff_in_before_and_after)
  after_before<-number
  
  col1<-colorRampPalette(brewer.pal(8,"RdPu"), bias=1)
  png(filename="after_before.png", units="in", width=1.3, height=1 , pointsize=3, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(after_before), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",density.info="none",cellnote=after_before, notecex=1.0,notecol="black",
            na.color=par("bg") ,key.par=list(mar=c(4,1,2,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2, margins=c(4,15))
  dev.off()
}

{
  number<-matrix(nrow=4,ncol=3)
  colnames(number)<-c("up","down","total")
  rownames(number)<-c("diff_in_circ_not_rectload","diff_in_rectload_not_circ","diff_in_circ_and_rectload_circ","diff_in_circ_and_rectload_rectl")
  
  number[1,1]<-nrow(subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1))
  number[1,2]<-nrow(subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1))
  number[1,3]<-nrow(diff_in_circ_not_rectload)
  number[2,1]<-nrow(subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1))
  number[2,2]<-nrow(subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1))
  number[2,3]<-nrow(diff_in_rectload_not_circ)
  number[3,1]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1))
  number[3,2]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1))
  number[3,3]<-nrow(diff_in_circ_and_rectload)
  number[4,1]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S2_S1>1))
  number[4,2]<-nrow(subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S2_S1<1))
  number[4,3]<-nrow(diff_in_circ_and_rectload)
  load_circ<-number
  
  col1<-colorRampPalette(brewer.pal(8,"RdPu"), bias=1)
  png(filename="load_circ.png", units="in", width=1.3, height=1 , pointsize=3, res=1200)
  par(font.lab=2, font=2, font.axis=2)
  heatmap.2(as.matrix(load_circ), col=col1(160), scale="none", 
            key=T, symkey=F, Rowv=F,Colv=FALSE,cexRow=1,
            trace='none',dendrogram="none",cexCol=0.9,srtCol =0,key.title=" ",key.xlab="Number of genes",cex.axis=0.7,
            ylab = "",density.info="none",cellnote=load_circ, notecex=1.0,notecol="black",
            na.color=par("bg") ,key.par=list(mar=c(4,1,2,1)),sepcolor="black",sepwidth=c(0.01,0.01),colsep=0:6,
            rowsep=0:6, font=2, margins=c(4,15))
  dev.off()
}

#Expressionpattern plots
{ 
  png(filename="Circ_rect_load_sig.png", units="in", width=2, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  plot(log2(diff_in_circ_not_rect$FC_S4_S3)~log2(diff_in_circ_not_rect$FC_S2_S1), pch=18, las=1, xlab="Rect+load/Rect",
       ylab="Circ+load/Circ",ylim=c(-5,5),xlim=c(-5,5), col=adjustcolor("purple", alpha.f = 0.8), cex=0.7)
  points(log2(diff_in_rect_not_circ$FC_S4_S3)~log2(diff_in_rect_not_circ$FC_S2_S1),pch=18, col=adjustcolor("green4", alpha.f = 0.8), cex=0.7)
  points(log2(diff_in_circ_and_rect$FC_S4_S3)~log2(diff_in_circ_and_rect$FC_S2_S1),pch=18, col=adjustcolor("red", alpha.f = 1), cex=0.7)
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="before_after_sig.png", units="in", width=2, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  plot(log2(diff_in_after_not_before$FC_S3_S1)~log2(diff_in_after_not_before$FC_S4_S2), pch=18, las=1, xlab="Circ+load/Rect+load",
       ylab="Circ/Rect",ylim=c(-5,5),xlim=c(-5,5), col=adjustcolor("purple", alpha.f = 0.8), cex=0.7)
  points(log2(diff_in_before_not_after$FC_S3_S1)~log2(diff_in_before_not_after$FC_S4_S2),pch=18, col=adjustcolor("green4", alpha.f = 0.8), cex=0.7)
  points(log2(diff_in_before_and_after$FC_S3_S1)~log2(diff_in_before_and_after$FC_S4_S2),pch=18, col=adjustcolor("red", alpha.f = 1), cex=0.7)
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
  png(filename="geo_and_load_sig.png", units="in", width=2, height=2, pointsize=6, res=1200)
  par(font.axis=2,font.lab=2)
  plot(log2(diff_in_circ_not_rectload$FC_S3_S1)~log2(diff_in_circ_not_rectload$FC_S2_S1), pch=18, las=1, xlab="Rect+load/Rect",
       ylab="Circ/Rect",ylim=c(-5,5),xlim=c(-5,5), col=adjustcolor("purple", alpha.f = 0.8), cex=0.7)
  points(log2(diff_in_rectload_not_circ$FC_S3_S1)~log2(diff_in_rectload_not_circ$FC_S2_S1),pch=18, col=adjustcolor("green4", alpha.f = 0.8), cex=0.7)
  points(log2(diff_in_circ_and_rectload$FC_S3_S1)~log2(diff_in_circ_and_rectload$FC_S2_S1),pch=18, col=adjustcolor("red", alpha.f = 1), cex=0.7)
  abline(v=2)
  abline(v=-2)
  abline(h=2)
  abline(h=-2)
  dev.off()
  
}

#pathways
a<-mm_pathway
c<-"pathway"
f1<-"pathway_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}

#TF
a<-mm_TF
c<-"TF"
f1<-"TF_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}

#GO
a<-mm_GO
c<-"GO"
f1<-"GO_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}

#miRNA
a<-mm_miRNA
c<-"miRNA"
f1<-"miRNA_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


#metabolic
a<-mm_metabolic
c<-"metabolic"
f1<-"metabolic_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


#location
a<-mm_location
c<-"location"
f1<-"location_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


#other
a<-mm_other
c<-"other"
f1<-"other_types.csv"
{
  
  diff_in_circ_and_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1>1 &
                                                                         diff_in_circ_and_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_and_rect_up")
  diff_in_circ_and_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rect,diff_in_circ_and_rect$FC_S2_S1<1 &
                                                                           diff_in_circ_and_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_and_rect_down")
  diff_in_circ_and_rect_path<-(diff_in_circ_and_rect_path_down)
  diff_in_rect_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rect_not_circ_up")
  diff_in_rect_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rect_not_circ,diff_in_rect_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rect_not_circ_down")
  diff_in_rect_not_circ_path<-rbind(diff_in_rect_not_circ_path_down,diff_in_rect_not_circ_path_up)
  diff_in_circ_not_rect_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3>1)$SYMBOL,a,c,"diff_in_circ_not_rect_up")
  diff_in_circ_not_rect_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rect,diff_in_circ_not_rect$FC_S4_S3<1)$SYMBOL,a,c,"diff_in_circ_not_rect_down")
  diff_in_circ_not_rect_path<-rbind(diff_in_circ_not_rect_path_down,diff_in_circ_not_rect_path_up)
  
  
  diff_in_before_and_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1>1 &
                                                                            diff_in_before_and_after$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_before_and_after_up")
  diff_in_before_and_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_and_after,diff_in_before_and_after$FC_S3_S1<1 &
                                                                              diff_in_before_and_after$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_before_and_after_down")
  diff_in_before_and_after_path<-rbind(diff_in_before_and_after_path_down,diff_in_before_and_after_path_up)
  diff_in_after_not_before_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2>1)$SYMBOL,a,c,"diff_in_after_not_before_up")
  diff_in_after_not_before_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_after_not_before,diff_in_after_not_before$FC_S4_S2<1)$SYMBOL,a,c,"diff_in_after_not_before_down")
  diff_in_after_not_before_path<-rbind(diff_in_after_not_before_path_down,diff_in_after_not_before_path_up)
  diff_in_before_not_after_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_before_not_after_up")
  diff_in_before_not_after_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_before_not_after,diff_in_before_not_after$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_before_not_after_down")
  diff_in_before_not_after_path<-rbind(diff_in_before_not_after_path_down,diff_in_before_not_after_path_up)
  
  
  diff_in_circ_and_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1>1 &
                                                                             diff_in_circ_and_rectload$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_circ_and_rectloadup")
  diff_in_circ_and_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_and_rectload,diff_in_circ_and_rectload$FC_S3_S1<1 &
                                                                               diff_in_circ_and_rectload$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_circ_and_rectload_down")
  diff_in_circ_and_rectload_path<-rbind(diff_in_circ_and_rectload_path_down,diff_in_circ_and_rectload_path_up)
  diff_in_circ_not_rectload_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1>1)$SYMBOL,a,c,"diff_in_circ_not_rectload_up")
  diff_in_circ_not_rectload_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_circ_not_rectload,diff_in_circ_not_rectload$FC_S3_S1<1)$SYMBOL,a,c,"diff_in_circ_not_rectload_down")
  diff_in_circ_not_rectload_path<-rbind(diff_in_circ_not_rectload_path_down,diff_in_circ_not_rectload_path_up)
  diff_in_rectload_not_circ_path_up<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1>1)$SYMBOL,a,c,"diff_in_rectload_not_circ_up")
  diff_in_rectload_not_circ_path_down<-enrichment(rnaseq_mapped_genes,subset(diff_in_rectload_not_circ,diff_in_rectload_not_circ$FC_S2_S1<1)$SYMBOL,a,c,"diff_in_rectload_not_circ_down")
  diff_in_rectload_not_circ_path<-rbind(diff_in_rectload_not_circ_path_down,diff_in_rectload_not_circ_path_up)
  
  pathway_types<-rbind(diff_in_circ_and_rect_path,diff_in_rect_not_circ_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rect_path)
  pathway_types<-rbind(pathway_types,diff_in_before_and_after_path)
  pathway_types<-rbind(pathway_types,diff_in_after_not_before_path)
  pathway_types<-rbind(pathway_types,diff_in_before_not_after_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_and_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_circ_not_rectload_path)
  pathway_types<-rbind(pathway_types,diff_in_rectload_not_circ_path)
  write.csv(pathway_types,file=f1)
}


