rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/gene_Expression_levels_ranked.csv", stringsAsFactors = FALSE)
#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part5_path_TF/"
dir.create(dirp)
setwd(dirp)
library(gskb)
library(gplots)
cols<-c("SYMBOL","S1","S2","S3","S4")
rnaseq_1<-rnaseq_mapped_genes[,cols]
rnaseq_1$SYMBOL<-toupper(rnaseq_1$SYMBOL)


{
  dirp1<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part5_path_TF/pathwayheat/"
  dir.create(dirp1)
  setwd(dirp1)
  
  data("mm_pathway")
  
  pathway_mean<-as.data.frame(matrix(nrow=892, ncol=5))
  colnames(pathway_mean)<-c("pathway","Rect","Rect_Load","Circ","Circ_Load")
  for(i in 1:892){
    temp<-subset(rnaseq_1,rnaseq_1$SYMBOL%in% (mm_pathway[[i]][3:length(mm_pathway[[i]])]))
    pathway_mean[i,1]<-mm_pathway[[i]][1]
    pathway_mean[i,2:5]<-apply(temp[,2:5], 2, FUN = median)
    pathway_mean[i,6]<-mm_pathway[[i]][2]
    
    if(nrow(temp)>2){
      png(filename=paste(mm_pathway[[i]][1],".png", sep=""), units="in", width=3, height=5 , pointsize=5, res=1200)
      par(font.axis=2,font.lab=2)
      x<-heatmap.2(as.matrix(temp[,2:5]), col=redgreen(75),
                   scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
                   trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="Row ZScore",cex.axis=0.7,
                   ylab = "",labRow = temp[,1],labCol=c("R","RL","C","CL"),density.info="none",
                   keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(3,8))
      dev.off()
    }
   
  }
  
  pathway_mean$RL_R<-log2(pathway_mean$Rect_Load/pathway_mean$Rect)
  pathway_mean$CL_C<-log2(pathway_mean$Circ_Load/pathway_mean$Circ)
  pathway_mean$S_R<-log2(pathway_mean$Circ/pathway_mean$Rect)
  pathway_mean$SL_RL<-log2(pathway_mean$Circ_Load/pathway_mean$Rect_Load)
  
  setwd(dirp)
  pathway_mean<-pathway_mean[order(pathway_mean$RL_R, decreasing=T),]
  png(filename="pathway_activity_increasing_Rect.png", units="in", width=5, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(pathway_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="Row ZScore",cex.axis=0.7,
               ylab = "",labRow = pathway_mean[1:30,1],labCol=c("R","RL","C","CL"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(4,50))
  dev.off()
  a<-pathway_mean[1:30,1]
  pathway_mean<-pathway_mean[order(pathway_mean$RL_R, decreasing=F),]
  png(filename="pathway_activity_decreasing_Rect.png", units="in", width=5, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(pathway_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="Row ZScore",cex.axis=0.7,
               ylab = "",labRow = pathway_mean[1:30,1],labCol=c("R","RL","C","CL"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(4,50))
  dev.off()
  b<-pathway_mean[1:30,1]
  pathway_mean<-pathway_mean[order(pathway_mean$CL_C, decreasing=T),]
  png(filename="pathway_activity_increasing_Circ.png", units="in", width=5, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(pathway_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="ZScore",cex.axis=0.7,
               ylab = "",labRow = pathway_mean[1:30,1],labCol=c("R","RL","C","CL"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(4,50))
  dev.off()
  c<-pathway_mean[1:30,1]
  
  pathway_mean<-pathway_mean[order(pathway_mean$CL_C, decreasing=F),]
  png(filename="pathway_activity_decreasing_Circ.png", units="in", width=5, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(pathway_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="ZScore",cex.axis=0.7,
               ylab = "",labRow = pathway_mean[1:30,1],labCol=c("R","RL","C","CL"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(4,50))
  dev.off()
  d<-pathway_mean[1:30,1]
  
  genelist <- list(a,b,c,d )
  ngene <- sapply(genelist, length)
  maxg <- seq_len(max(ngene))
  top30_pathway <- as.data.frame(sapply(genelist, "[", i = maxg))
  names(top30_pathway)<-c("incr_R","decr_R","incr_C","decr_C")
  write.csv(top30_pathway,file="Top30 pathway.csv")
  
  write.csv(pathway_mean,file="pathway_mean.csv")
}

{
  dirp1<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part5_path_TF/TFheat/"
  dir.create(dirp1)
  setwd(dirp1)
  
  data("mm_TF")
  
  TF_mean<-as.data.frame(matrix(nrow=372, ncol=5))
  colnames(TF_mean)<-c("pathway","Rect","Rect_Load","Circ","Circ_Load")
  for(i in 1:372){
    temp<-subset(rnaseq_1,rnaseq_1$SYMBOL%in% (mm_TF[[i]][3:length(mm_TF[[i]])]))
    TF_mean[i,1]<-mm_TF[[i]][1]
    TF_mean[i,2:5]<-apply(temp[,2:5], 2, FUN = median)
    TF_mean[i,6]<-mm_TF[[i]][2]
    if(nrow(temp)>4 ){
      png(filename=paste(mm_TF[[i]][1],".png", sep=""), units="in", width=3, height=5 , pointsize=5, res=1200)
      par(font.axis=2,font.lab=2)
      x<-heatmap.2(as.matrix(temp[,2:5]), col=redgreen(75),
                   scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=0.7,
                   trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="Row ZScore",cex.axis=0.7,
                   ylab = "",labRow = temp[,1],labCol=c("R","RL","C","CL"),density.info="none",
                   keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(3,8))
      dev.off()
    }
    
  }
  
  TF_mean$RL_R<-log2(TF_mean$Rect_Load/TF_mean$Rect)
  TF_mean$CL_C<-log2(TF_mean$Circ_Load/TF_mean$Circ)
  TF_mean$S_R<-log2(TF_mean$Circ/TF_mean$Rect)
  TF_mean$SL_RL<-log2(TF_mean$Circ_Load/TF_mean$Rect_Load)
  
  setwd(dirp)
  TF_mean<-TF_mean[order(TF_mean$RL_R, decreasing=T),]
  png(filename="TF_activity_increasing_Rect.png", units="in", width=4, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(TF_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="Row ZScore",cex.axis=0.7,
               ylab = "",labRow = TF_mean[1:30,1],labCol=c("Rect","Rect_Load","Circ","Circ_Load"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(6,20))
  dev.off()
  a<-TF_mean[1:30,1]
  TF_mean<-TF_mean[order(TF_mean$RL_R, decreasing=F),]
  png(filename="TF_activity_decreasing_Rect.png", units="in", width=4, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(TF_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="Row ZScore",cex.axis=0.7,
               ylab = "",labRow = TF_mean[1:30,1],labCol=c("Rect","Rect_Load","Circ","Circ_Load"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(6,20))
  dev.off()
  b<-TF_mean[1:30,1]
  TF_mean<-TF_mean[order(TF_mean$CL_C, decreasing=T),]
  png(filename="TF_activity_increasing_Circ.png", units="in", width=4, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(TF_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="ZScore",cex.axis=0.7,
               ylab = "",labRow = TF_mean[1:30,1],labCol=c("Rect","Rect_Load","Circ","Circ_Load"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(6,20))
  dev.off()
  c<-TF_mean[1:30,1]
  
  TF_mean<-TF_mean[order(TF_mean$CL_C, decreasing=F),]
  png(filename="TF_activity_decreasing_Circ.png", units="in", width=4, height=5 , pointsize=5, res=1200)
  par(font.axis=2,font.lab=2)
  x<-heatmap.2(as.matrix(TF_mean[1:30,2:5]), col=redgreen(75),
               scale="row",key=T, symkey=T, Rowv=F,Colv=FALSE,cexRow=1,
               trace='none',dendrogram="none",cexCol=1.2,srtCol =90,key.title=" ",key.xlab="ZScore",cex.axis=0.7,
               ylab = "",labRow = TF_mean[1:30,1],labCol=c("Rect","Rect_Load","Circ","Circ_Load"),density.info="none",
               keysize=0.6, key.par=list(mar=c(4,1,4,1)), margins=c(6,20))
  dev.off()
  d<-TF_mean[1:30,1]
  
  genelist <- list(a,b,c,d )
  ngene <- sapply(genelist, length)
  maxg <- seq_len(max(ngene))
  top30_pathway <- as.data.frame(sapply(genelist, "[", i = maxg))
  names(top30_pathway)<-c("incr_R","decr_R","incr_C","decr_C")
  write.csv(top30_pathway,file="Top30 pathway.csv")
  
  write.csv(TF_mean,file="TF_mean.csv")
}