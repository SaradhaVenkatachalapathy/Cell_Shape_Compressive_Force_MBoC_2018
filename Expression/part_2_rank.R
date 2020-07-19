#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part2/"
dir.create(dirp)
setwd(dirp)

rnaseq_mapped_genes <- read.csv("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/allcombined_gene_level_expression_mapped_to_gene_name.csv", stringsAsFactors = FALSE)

# function to return the rank of a vector
calcRank <- function(x){
  sorted <- x[order(x)]
  ranks <- sapply(x, function(x) which(sorted==x))
  return(rank(x))
}

# caculate the relative differential expression of each gene by getting the ranks
rnaseq_mapped_genes$Rank_S1<-0
rnaseq_mapped_genes$Rank_S2<-0
rnaseq_mapped_genes$Rank_S3<-0
rnaseq_mapped_genes$Rank_S4<-0
for( i in 1:nrow(rnaseq_mapped_genes)){
  rnaseq_mapped_genes[i,grepl( "Rank" , names( rnaseq_mapped_genes ) )]<-calcRank(rnaseq_mapped_genes[i, which(names(rnaseq_mapped_genes) %in% c("S1", "S2", "S3","S4"))])
}

# generating the string to represent the trend in gene expression
rnaseq_mapped_genes$raw_order<-NA
rnaseq_mapped_genes$order<-NA

for( i in 1:nrow(rnaseq_mapped_genes)){
  rnaseq_mapped_genes[i,]$raw_order<-paste(rnaseq_mapped_genes[i,]$Rank_S1,rnaseq_mapped_genes[i,]$Rank_S2,rnaseq_mapped_genes[i,]$Rank_S3,rnaseq_mapped_genes[i,]$Rank_S4,sep="-")
  t<-rnaseq_mapped_genes[i,which(names(rnaseq_mapped_genes) %in% c("Rank_S1", "Rank_S2", "Rank_S3","Rank_S4"))]
  names<-c("S1","S2","S3","S4")
  t1<-names[order(t)]
  rnaseq_mapped_genes[i,]$order<-paste(t1[1],t1[2],t1[3],t1[4],sep="<")
}

# plot the global trends of expression

x<-as.data.frame(table(rnaseq_mapped_genes$raw_order))
x<-(subset(x,!grepl( ".5-" , x[,1] )))
x<-x[order(x[,2]),]
png(filename="global_rank_frequencies.png", units="in", width=2, height=2 , pointsize=3, res=300)
par(mar=c(4,15,3,2))
barplot(x[,2], main = "",horiz = T,las=1,names=x[,1],xlab="number of genes")
box()
dev.off()

#plot the frequency of rank as a stacked barblot
counts<-as.data.frame(table(rnaseq_mapped_genes$Rank_S1))
counts[,2]<-counts[,2]/nrow(rnaseq_mapped_genes)
counts[,3]<-as.data.frame(table(rnaseq_mapped_genes$Rank_S2))[,2]/nrow(rnaseq_mapped_genes)
counts[,4]<-as.data.frame(table(rnaseq_mapped_genes$Rank_S3))[,2]/nrow(rnaseq_mapped_genes)
counts[,5]<-as.data.frame(table(rnaseq_mapped_genes$Rank_S4))[,2]/nrow(rnaseq_mapped_genes)
colnames(counts)[c(2:5)]<-c("S1","S2","S3","S4")
counts<-subset(counts,counts[,1]%in% c(1,2,3,4))

library(gplots)
png(filename="global_rank_frequencies_stacked_barplot.png", units="in", width=2, height=2 , pointsize=4, res=300)
barplot(as.matrix(counts[,2:5]), las=1, xlim=c(0,6),col=redgreen(5),ylab="Fraction of the total number of genes")
legend("bottomright",legend = c("Rank 4","Rank3","Rank2","Rank1"), fill = redgreen(5))
box()
dev.off()

library(gplots)
png(filename="global_rank_frequencies_stacked_barplot_named_small.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2, font=2, font.lab=2,mar = c(4,4,2,2))
barplot(as.matrix(counts[,2:5]), las=1, xlim=c(0,7),col=redgreen(5),ylab="Fraction of genes", 
        names=c("R", "RL","C", "CL"), cex.names = 1.2, cex.lab=1.5, cex.axis=1.2)
legend("bottomright",legend = c("Rank 4","Rank3","Rank2","Rank1"), fill = redgreen(5), cex=0.88,border="white")
box()
dev.off()

png(filename="global_rank_frequencies_stacked_barplot_named_big.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2, font=2, font.lab=2,mar = c(4,4,2,2))
barplot(as.matrix(counts[,2:5]), las=1, xlim=c(0,7),col=redgreen(5),ylab="Fraction of genes", 
        names=c("R", "RL","C", "CL"), cex.names = 1.2, cex.lab=1.5, cex.axis=1)
legend("bottomright",legend = c("Rank 4","Rank3","Rank2","Rank1"), fill = redgreen(5), cex=0.88)
box()
dev.off()


png(filename="zscores_gene_box_small.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2, font=2,mar = c(4,7,2,2))
boxplot(rnaseq_mapped_genes[,grepl( "zscore" , names( rnaseq_mapped_genes ) )],col="gray", las=1,ylab="Relative gene expression\nZScore", lty=1, pch=18, 
        names=c("R","RL","C","CL"),cex.names = 1.2, cex.lab=1.5, cex.axis=1.2, lwd=0.8)
dev.off()

png(filename="zscores_gene_box_big.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2, font=2,mar = c(4,7,2,2))
boxplot(rnaseq_mapped_genes[,grepl( "zscore" , names( rnaseq_mapped_genes ) )],col="gray", las=1,ylab="Relative gene expression\nZScore", lty=1, pch=18, 
        names=c("R","RL","C","CL"),cex.names = 1.2, cex.lab=1.5, cex.axis=1.2)
dev.off()

png(filename="zscore_density_big_circle.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((rnaseq_mapped_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((rnaseq_mapped_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="zscore_density_small_circle.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((rnaseq_mapped_genes$S4_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((rnaseq_mapped_genes$S3_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("C","CL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="zscore_density_big_rect.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((rnaseq_mapped_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((rnaseq_mapped_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

png(filename="zscore_density_small_rect.png", units="in", width=1, height=1 , pointsize=3, res=1200)
par(font.axis=2,font.lab=2,mar = c(6,4,2,2)+0.2, font=2)
d <- density((rnaseq_mapped_genes$S2_zscore), na.rm=T) # returns the density data 
plot(d, col="red", las=1, lwd=2, main="", xlab="Relative gene expression(ZScore)", cex.lab=1.3, cex.axis=1)
d <- density((rnaseq_mapped_genes$S1_zscore), na.rm=T) # returns the density data 
lines(d,col="green4",lwd=2)
legend("topleft",c("R","RL"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green4","red"), cex=0.8) 
dev.off()

setwd(dird)
write.csv(rnaseq_mapped_genes, file="gene_Expression_levels_ranked.csv")
setwd(dirp)

