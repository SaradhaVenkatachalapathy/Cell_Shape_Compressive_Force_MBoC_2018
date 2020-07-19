imp<-c("Atr","Atm","Parp1", "Pten", 
       "Bad","Casp7","Casp8","Mcl1",
       "Akt3","Mtor","Foxo6", "Foxo3", 
       "Lmnb1","Ctnnb1",
       "Atg3", "Cdc25b", "Gabarapl2", "Pik3c3")
imp_genes1<-subset(rnaseq_mapped_genes,rnaseq_mapped_genes$SYMBOL %in% imp)
imp_genes<-imp_genes1[match(imp, imp_genes1$SYMBOL),]

png(filename="Pathways_new.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2)
x<-heatmap.2(as.matrix(imp_genes[,which(names(imp_genes) %in% c("S1_zscore", "S2_zscore", "S3_zscore","S4_zscore"))]), col=redblue(75),
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

