#create folder to store plots
dird<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/"
dirp<-"/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/rna_seq_analysis/plots_part1/"
dir.create(dirp)
setwd(dirp)

# Read in all the RNA seq read counts files (output of htseq)
S1ai<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B1T1.tabular", header=FALSE)
S1aii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B1T2.tabular", header=FALSE)
S1aiii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B1T3.tabular", header=FALSE)
S1aiv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B1T4.tabular", header=FALSE)
colnames(S1ai)<-c("TXID","S1_B1_L1")
colnames(S1aii)<-c("TXID","S1_B1_L2")
colnames(S1aiii)<-c("TXID","S1_B1_L3")
colnames(S1aiv)<-c("TXID","S1_B1_L4")

S2ai<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B1T1.tabular", header=FALSE)
S2aii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B1T2.tabular", header=FALSE)
S2aiii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B1T3.tabular", header=FALSE)
S2aiv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B1T4.tabular", header=FALSE)
colnames(S2ai)<-c("TXID","S2_B1_L1")
colnames(S2aii)<-c("TXID","S2_B1_L2")
colnames(S2aiii)<-c("TXID","S2_B1_L3")
colnames(S2aiv)<-c("TXID","S2_B1_L4")

S3ai<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B1T1.tabular", header=FALSE)
S3aii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B1T2.tabular", header=FALSE)
S3aiii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B1T3.tabular", header=FALSE)
S3aiv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B1T4.tabular", header=FALSE)
colnames(S3ai)<-c("TXID","S3_B1_L1")
colnames(S3aii)<-c("TXID","S3_B1_L2")
colnames(S3aiii)<-c("TXID","S3_B1_L3")
colnames(S3aiv)<-c("TXID","S3_B1_L4")

S4ai<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B1T1.tabular", header=FALSE)
S4aii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B1T2.tabular", header=FALSE)
S4aiii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B1T3.tabular", header=FALSE)
S4aiv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B1T4.tabular", header=FALSE)
colnames(S4ai)<-c("TXID","S4_B1_L1")
colnames(S4aii)<-c("TXID","S4_B1_L2")
colnames(S4aiii)<-c("TXID","S4_B1_L3")
colnames(S4aiv)<-c("TXID","S4_B1_L4")

S1bi<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B2T1.tabular", header=FALSE)
S1bii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B2T2.tabular", header=FALSE)
S1biii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B2T3.tabular", header=FALSE)
S1biv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B2T4.tabular", header=FALSE)
colnames(S1bi)<-c("TXID","S1_B2_L1")
colnames(S1bii)<-c("TXID","S1_B2_L2")
colnames(S1biii)<-c("TXID","S1_B2_L3")
colnames(S1biv)<-c("TXID","S1_B2_L4")

S2bi<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B2T1.tabular", header=FALSE)
S2bii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B2T2.tabular", header=FALSE)
S2biii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B2T3.tabular", header=FALSE)
S2biv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B2T4.tabular", header=FALSE)
colnames(S2bi)<-c("TXID","S2_B2_L1")
colnames(S2bii)<-c("TXID","S2_B2_L2")
colnames(S2biii)<-c("TXID","S2_B2_L3")
colnames(S2biv)<-c("TXID","S2_B2_L4")

S3bi<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B2T1.tabular", header=FALSE)
S3bii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B2T2.tabular", header=FALSE)
S3biii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B2T3.tabular", header=FALSE)
S3biv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B2T4.tabular", header=FALSE)
colnames(S3bi)<-c("TXID","S3_B2_L1")
colnames(S3bii)<-c("TXID","S3_B2_L2")
colnames(S3biii)<-c("TXID","S3_B2_L3")
colnames(S3biv)<-c("TXID","S3_B2_L4")

S4bi<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B2T1.tabular", header=FALSE)
S4bii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B2T2.tabular", header=FALSE)
S4biii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B2T3.tabular", header=FALSE)
S4biv<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B2T4.tabular", header=FALSE)
colnames(S4bi)<-c("TXID","S4_B2_L1")
colnames(S4bii)<-c("TXID","S4_B2_L2")
colnames(S4biii)<-c("TXID","S4_B2_L3")
colnames(S4biv)<-c("TXID","S4_B2_L4")

S1ci<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B3T1.tabular", header=FALSE)
S1cii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B3T2.tabular", header=FALSE)
S1ciii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B3T3.tabular", header=FALSE)
S1civ<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S1B3T4.tabular", header=FALSE)
colnames(S1ci)<-c("TXID","S1_B3_L1")
colnames(S1cii)<-c("TXID","S1_B3_L2")
colnames(S1ciii)<-c("TXID","S1_B3_L3")
colnames(S1civ)<-c("TXID","S1_B3_L4")

S2ci<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B3T1.tabular", header=FALSE)
S2cii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B3T2.tabular", header=FALSE)
S2ciii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B3T3.tabular", header=FALSE)
S2civ<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S2B3T4.tabular", header=FALSE)
colnames(S2ci)<-c("TXID","S2_B3_L1")
colnames(S2cii)<-c("TXID","S2_B3_L2")
colnames(S2ciii)<-c("TXID","S2_B3_L3")
colnames(S2civ)<-c("TXID","S2_B3_L4")

S3ci<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B3T1.tabular", header=FALSE)
S3cii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B3T2.tabular", header=FALSE)
S3ciii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B3T3.tabular", header=FALSE)
S3civ<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S3B3T4.tabular", header=FALSE)
colnames(S3ci)<-c("TXID","S3_B3_L1")
colnames(S3cii)<-c("TXID","S3_B3_L2")
colnames(S3ciii)<-c("TXID","S3_B3_L3")
colnames(S3civ)<-c("TXID","S3_B3_L4")

S4ci<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B3T1.tabular", header=FALSE)
S4cii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B3T2.tabular", header=FALSE)
S4ciii<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B3T3.tabular", header=FALSE)
S4civ<- read.delim("/Users/saradhavenkatachalapathy/Desktop/kathik rna seq/htseq_tophat/S4B3T4.tabular", header=FALSE)
colnames(S4ci)<-c("TXID","S4_B3_L1")
colnames(S4cii)<-c("TXID","S4_B3_L2")
colnames(S4ciii)<-c("TXID","S4_B3_L3")
colnames(S4civ)<-c("TXID","S4_B3_L4")

# merge all the counts into one large dataframe by the transcript id. 

rnaseq <- merge(S1ai,S1aii,by="TXID")
rnaseq <- merge(rnaseq,S1aiii,by="TXID")
rnaseq <- merge(rnaseq,S1aiv,by="TXID")
rnaseq <- merge(rnaseq,S2ai,by="TXID")
rnaseq <- merge(rnaseq,S2aii,by="TXID")
rnaseq <- merge(rnaseq,S2aiii,by="TXID")
rnaseq <- merge(rnaseq,S2aiv,by="TXID")
rnaseq <- merge(rnaseq,S3ai,by="TXID")
rnaseq <- merge(rnaseq,S3aii,by="TXID")
rnaseq <- merge(rnaseq,S3aiii,by="TXID")
rnaseq <- merge(rnaseq,S3aiv,by="TXID")
rnaseq <- merge(rnaseq,S4ai,by="TXID")
rnaseq <- merge(rnaseq,S4aii,by="TXID")
rnaseq <- merge(rnaseq,S4aiii,by="TXID")
rnaseq <- merge(rnaseq,S4aiv,by="TXID")
rnaseq <- merge(rnaseq,S1bi,by="TXID")
rnaseq <- merge(rnaseq,S1bii,by="TXID")
rnaseq <- merge(rnaseq,S1biii,by="TXID")
rnaseq <- merge(rnaseq,S1biv,by="TXID")
rnaseq <- merge(rnaseq,S2bi,by="TXID")
rnaseq <- merge(rnaseq,S2bii,by="TXID")
rnaseq <- merge(rnaseq,S2biii,by="TXID")
rnaseq <- merge(rnaseq,S2biv,by="TXID")
rnaseq <- merge(rnaseq,S3bi,by="TXID")
rnaseq <- merge(rnaseq,S3bii,by="TXID")
rnaseq <- merge(rnaseq,S3biii,by="TXID")
rnaseq <- merge(rnaseq,S3biv,by="TXID")
rnaseq <- merge(rnaseq,S4bi,by="TXID")
rnaseq <- merge(rnaseq,S4bii,by="TXID")
rnaseq <- merge(rnaseq,S4biii,by="TXID")
rnaseq <- merge(rnaseq,S4biv,by="TXID")
rnaseq <- merge(rnaseq,S1ci,by="TXID")
rnaseq <- merge(rnaseq,S1cii,by="TXID")
rnaseq <- merge(rnaseq,S1ciii,by="TXID")
rnaseq <- merge(rnaseq,S1civ,by="TXID")
rnaseq <- merge(rnaseq,S2ci,by="TXID")
rnaseq <- merge(rnaseq,S2cii,by="TXID")
rnaseq <- merge(rnaseq,S2ciii,by="TXID")
rnaseq <- merge(rnaseq,S2civ,by="TXID")
rnaseq <- merge(rnaseq,S3ci,by="TXID")
rnaseq <- merge(rnaseq,S3cii,by="TXID")
rnaseq <- merge(rnaseq,S3ciii,by="TXID")
rnaseq <- merge(rnaseq,S3civ,by="TXID")
rnaseq <- merge(rnaseq,S4ci,by="TXID")
rnaseq <- merge(rnaseq,S4cii,by="TXID")
rnaseq <- merge(rnaseq,S4ciii,by="TXID")
rnaseq <- merge(rnaseq,S4civ,by="TXID")
# sum up the technical replicates

rnaseq$S1_B1 <- rowSums(rnaseq[,c("S1_B1_L1", "S1_B1_L2","S1_B1_L3", "S1_B1_L4")], na.rm=T)
rnaseq$S2_B1 <- rowSums(rnaseq[,c("S2_B1_L1", "S2_B1_L2","S2_B1_L3", "S2_B1_L4")], na.rm=T)
rnaseq$S3_B1 <- rowSums(rnaseq[,c("S3_B1_L1", "S3_B1_L2","S3_B1_L3", "S3_B1_L4")], na.rm=T)
rnaseq$S4_B1 <- rowSums(rnaseq[,c("S4_B1_L1", "S4_B1_L2","S4_B1_L3", "S4_B1_L4")], na.rm=T)

rnaseq$S1_B2 <- rowSums(rnaseq[,c("S1_B2_L1", "S1_B2_L2","S1_B2_L3", "S1_B2_L4")], na.rm=T)
rnaseq$S2_B2 <- rowSums(rnaseq[,c("S2_B2_L1", "S2_B2_L2","S2_B2_L3", "S2_B2_L4")], na.rm=T)
rnaseq$S3_B2 <- rowSums(rnaseq[,c("S3_B2_L1", "S3_B2_L2","S3_B2_L3", "S3_B2_L4")], na.rm=T)
rnaseq$S4_B2 <- rowSums(rnaseq[,c("S4_B2_L1", "S4_B2_L2","S4_B2_L3", "S4_B2_L4")], na.rm=T)

rnaseq$S1_B3 <- rowSums(rnaseq[,c("S1_B3_L1", "S1_B3_L2","S1_B3_L3", "S1_B3_L4")], na.rm=T)
rnaseq$S2_B3 <- rowSums(rnaseq[,c("S2_B3_L1", "S2_B3_L2","S2_B3_L3", "S2_B3_L4")], na.rm=T)
rnaseq$S3_B3 <- rowSums(rnaseq[,c("S3_B3_L1", "S3_B3_L2","S3_B3_L3", "S3_B3_L4")], na.rm=T)
rnaseq$S4_B3 <- rowSums(rnaseq[,c("S4_B3_L1", "S4_B3_L2","S4_B3_L3", "S4_B3_L4")], na.rm=T)

setwd(dirp)
#source("https://bioconductor.org/biocLite.R")
#biocLite("EnsDb.Mmusculus.v79", dependencies=TRUE)
library(EnsDb.Mmusculus.v79)
keys <- keys(EnsDb.Mmusculus.v79)
columns(EnsDb.Mmusculus.v79)
anno<-select(EnsDb.Mmusculus.v79, keys=keys, columns=c("GENEID","SYMBOL","GENENAME","ENTREZID",
                                                       "GENESEQEND","GENESEQSTART","GENEBIOTYPE",
                                                       "TXID","TXSEQEND","TXSEQSTART","TXBIOTYPE",
                                                       "SEQNAME"),keytype="GENEID")
anno$txlength<-anno$TXSEQEND-anno$TXSEQSTART
rnaseq <- merge(rnaseq,anno,by="TXID")
png(filename="histogram_lengths.png", units="in", width=2, height=2, pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(log(rnaseq$txlength, base=2), col="gray", las=1, cex.axis=0.8,xlab="length of transcript",main="")
t<-quantile(rnaseq$txlength)
text(5,12000,paste("Median=",t[3],sep=""), cex=1)
box()
dev.off()


# normalising the data 
rnaseq$S1_B1_TPKM<-rnaseq$S1_B1/(rnaseq$txlength/1000)
rnaseq$S2_B1_TPKM<-rnaseq$S2_B1/(rnaseq$txlength/1000)
rnaseq$S3_B1_TPKM<-rnaseq$S3_B1/(rnaseq$txlength/1000)
rnaseq$S4_B1_TPKM<-rnaseq$S4_B1/(rnaseq$txlength/1000)
rnaseq$S1_B2_TPKM<-rnaseq$S1_B2/(rnaseq$txlength/1000)
rnaseq$S2_B2_TPKM<-rnaseq$S2_B2/(rnaseq$txlength/1000)
rnaseq$S3_B2_TPKM<-rnaseq$S3_B2/(rnaseq$txlength/1000)
rnaseq$S4_B2_TPKM<-rnaseq$S4_B2/(rnaseq$txlength/1000)
rnaseq$S1_B3_TPKM<-rnaseq$S1_B3/(rnaseq$txlength/1000)
rnaseq$S2_B3_TPKM<-rnaseq$S2_B3/(rnaseq$txlength/1000)
rnaseq$S3_B3_TPKM<-rnaseq$S3_B3/(rnaseq$txlength/1000)
rnaseq$S4_B3_TPKM<-rnaseq$S4_B3/(rnaseq$txlength/1000)

total_reads<-colSums(rnaseq[,grepl( "TPKM" , names( rnaseq ) )])/10^6

rnaseq$S1_B1_TPKM<-rnaseq$S1_B1_TPKM/total_reads[1]
rnaseq$S2_B1_TPKM<-rnaseq$S2_B1_TPKM/total_reads[2]
rnaseq$S3_B1_TPKM<-rnaseq$S3_B1_TPKM/total_reads[3]
rnaseq$S4_B1_TPKM<-rnaseq$S4_B1_TPKM/total_reads[4]
rnaseq$S1_B2_TPKM<-rnaseq$S1_B2_TPKM/total_reads[5]
rnaseq$S2_B2_TPKM<-rnaseq$S2_B2_TPKM/total_reads[6]
rnaseq$S3_B2_TPKM<-rnaseq$S3_B2_TPKM/total_reads[7]
rnaseq$S4_B2_TPKM<-rnaseq$S4_B2_TPKM/total_reads[8]
rnaseq$S1_B3_TPKM<-rnaseq$S1_B3_TPKM/total_reads[9]
rnaseq$S2_B3_TPKM<-rnaseq$S2_B3_TPKM/total_reads[10]
rnaseq$S3_B3_TPKM<-rnaseq$S3_B3_TPKM/total_reads[11]
rnaseq$S4_B3_TPKM<-rnaseq$S4_B3_TPKM/total_reads[12]

png(filename="totalcountsbeforeandafternormalisation.png", units="in", width=4, height=2 , pointsize=5, res=1200)
par(mfrow=c(1,2),font.axis=2,font.lab=2, mar=c(5,6,4,4))
barplot((total_reads*10^6), las=1, horiz = T, xlab="total reads", cex.axis=0.7,cex.names=0.8, main="before")
box()
barplot(colSums(rnaseq[,grepl( "TPKM" , names( rnaseq ) )]),las=1, horiz = T,xlab="total reads", cex.axis=0.7,cex.names=0.8, main="after")
box()
dev.off()

# mean bio_replicates
rnaseq$S1<-rowMeans(subset(rnaseq, select = c(S1_B1_TPKM,S1_B2_TPKM,S1_B3_TPKM)), na.rm = TRUE)
rnaseq$S2<-rowMeans(subset(rnaseq, select = c(S2_B1_TPKM,S2_B2_TPKM,S2_B3_TPKM)), na.rm = TRUE)
rnaseq$S3<-rowMeans(subset(rnaseq, select = c(S3_B1_TPKM,S3_B2_TPKM,S3_B3_TPKM)), na.rm = TRUE)
rnaseq$S4<-rowMeans(subset(rnaseq, select = c(S4_B1_TPKM,S4_B2_TPKM,S4_B3_TPKM)), na.rm = TRUE)
rnaseq$Pop_Mean<-apply(rnaseq[, which(names(rnaseq) %in% c("S1", "S2", "S3","S4"))],1,mean) 
rnaseq$Pop_Std_Dev<-apply(rnaseq[, which(names(rnaseq) %in% c("S1", "S2", "S3","S4"))],1,sd) 

rnaseq$S1_zscore<-(rnaseq$S1-rnaseq$Pop_Mean)/rnaseq$Pop_Std_Dev
rnaseq$S2_zscore<-(rnaseq$S2-rnaseq$Pop_Mean)/rnaseq$Pop_Std_Dev
rnaseq$S3_zscore<-(rnaseq$S3-rnaseq$Pop_Mean)/rnaseq$Pop_Std_Dev
rnaseq$S4_zscore<-(rnaseq$S4-rnaseq$Pop_Mean)/rnaseq$Pop_Std_Dev

# Fold changes

rnaseq$FC_S2_S1<-rnaseq$S2/rnaseq$S1#S2/S1
rnaseq$FC_S3_S1<-rnaseq$S3/rnaseq$S1#S3/S1
rnaseq$FC_S4_S1<-rnaseq$S4/rnaseq$S1#S4/S1
rnaseq$FC_S3_S2<-rnaseq$S3/rnaseq$S2#S3/S2
rnaseq$FC_S4_S2<-rnaseq$S4/rnaseq$S2#S4/S2
rnaseq$FC_S4_S3<-rnaseq$S4/rnaseq$S3#S4/S3


png(filename="Foldchange_transcript.png", units="in", width=3, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
boxplot(log(rnaseq[,grepl( "FC" , names( rnaseq ) )], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18)
dev.off()
png(filename="zscores_transcript.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
boxplot(rnaseq[,grepl( "zscore" , names( rnaseq ) )],col="gray", las=1,ylab="Zscores", lty=1, pch=18, names=c("S1","S2","S3","S4"))
dev.off()
png(filename="histogram_pop_mean_transcript_levels.png", units="in", width=2, height=2, pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(log(rnaseq$Pop_Mean, base=2), col="gray", las=1, cex.axis=0.8,xlab="log2 mean transcript expression",main="")
t<-quantile(rnaseq$Pop_Mean, na.rm=T)
text(6,6000,paste("Median=",t[3],sep=""), cex=1)
box()
dev.off()
png(filename="histogram_pop_SD_transcript_levels.png", units="in", width=2, height=2, pointsize=5, res=1200)
par(font.axis=2,font.lab=2)
hist(log(rnaseq$Pop_Std_Dev, base=2), col="gray", las=1, cex.axis=0.8,xlab="log2 mean SD transcript expression",main="")
t<-quantile(rnaseq[,36])
text(6,6000,paste("Median=",t[3],sep=""), cex=1)
box()
dev.off()

setwd(dird)
write.csv(rnaseq, file="allcombined_transcript_level_expression.csv")
setwd(dirp)

# correlations

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr");on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(10^x, 10^y), digits=4)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 18)
}


# Create the plots
png(filename="S1_B1.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B1_L1","S1_B1_L2","S1_B1_L3","S1_B1_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S2_B1.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S2_B1_L1","S2_B1_L2","S2_B1_L3","S2_B1_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S3_B1.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S3_B1_L1","S3_B1_L2","S3_B1_L3","S3_B1_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S4_B1.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S4_B1_L1","S4_B1_L2","S4_B1_L3","S4_B1_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()

png(filename="S1_B2_i.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B2_L1","S1_B2_L2","S1_B2_L3","S1_B2_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S2_B2.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S2_B2_L1","S2_B2_L2","S2_B2_L3","S2_B2_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S3_B2.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S3_B2_L1","S3_B2_L2","S3_B2_L3","S3_B2_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S4_B2.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S4_B2_L1","S4_B2_L2","S4_B2_L3","S4_B2_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()

png(filename="S1_B3.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B3_L1","S1_B3_L2","S1_B3_L3","S1_B3_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S2_B3.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S2_B3_L1","S2_B3_L2","S2_B3_L3","S2_B3_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S3_B3.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S3_B3_L1","S3_B3_L2","S3_B3_L3","S3_B3_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S4_B3.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S4_B3_L1","S4_B3_L2","S4_B3_L3","S4_B3_L4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()


png(filename="B1.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B1","S2_B1","S3_B1","S4_B1")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="B2.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B2","S2_B2","S3_B2","S4_B2")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="B3.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B3","S2_B3","S3_B3","S4_B3")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()


png(filename="S1.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B1","S1_B2","S2_B3")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S2.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S2_B1","S2_B2","S2_B3")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S3.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S3_B1","S3_B2","S3_B3")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S4.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S4_B1","S4_B2","S4_B3")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()

png(filename="B1_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B1_TPKM","S2_B1_TPKM","S3_B1_TPKM","S4_B1_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="B2_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B2_TPKM","S2_B2_TPKM","S3_B2_TPKM","S4_B2_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="B3_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B3_TPKM","S2_B3_TPKM","S3_B3_TPKM","S4_B3_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S1_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1_B1_TPKM","S1_B2_TPKM","S1_B3_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S2_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S2_B1_TPKM","S2_B2_TPKM","S2_B3_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S3_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S3_B1_TPKM","S3_B2_TPKM","S3_B3_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="S4_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S4_B1_TPKM","S4_B2_TPKM","S4_B3_TPKM")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="combined_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
plot(log2(rnaseq[,grepl( "TPKM" , names( rnaseq ) )]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()
png(filename="ALL_norm.png", units="in", width=2, height=2 , pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
t<-c("S1","S2","S3","S4")
plot(log2(rnaseq[,which(names(rnaseq) %in% t)]), lower.panel = panel.cor,upper.panel = upper.panel)
dev.off()

library(corrplot)
t<-cor(rnaseq[,grepl( "_L" , names( rnaseq ) )])
png(filename="corrplot_all_samples_before_norm.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
corrplot(t,is.cor=F, tl.col="black",cl.pos="b", cl.ratio=0.4, method="color")
dev.off()
x<-c("S1_B1","S2_B1","S3_B1","S4_B1","S1_B2","S2_B2","S3_B2","S4_B2","S1_B3","S2_B3","S3_B3","S4_B3")
t<-cor(rnaseq[, which(names(rnaseq) %in% x)])
png(filename="corrplot_samples_before_norm.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
corrplot(t, is.cor=F,tl.col="black",cl.pos="b", cl.ratio=0.4, method="color")
dev.off()
t<-cor(rnaseq[,grepl( "_TPKM" , names( rnaseq ) )])
png(filename="corrplot_samples_after_norm.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
corrplot(t, is.cor=F, tl.col="black",cl.pos="b", cl.ratio=0.4, method="color")
dev.off()
x<-c("S1","S2","S3","S4")
t<-cor(rnaseq[, which(names(rnaseq) %in% x)])
png(filename="corrplot_after_norm.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
corrplot(t, is.cor=F, tl.col="black",cl.pos="b", cl.ratio=0.4, method="color")
dev.off()

# log transform tpkm
exon <- log2(rnaseq[,grepl( "_TPKM" , names( rnaseq ) )])
is.na(exon) <- do.call(cbind,lapply(exon, is.infinite)) 
exon.pca <- prcomp(na.omit(exon),scale= T) 
t<-as.data.frame(exon.pca$rotation)
png(filename="pca_tpkm_log_scaled.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plot(t[,2]~t[,1],pch=18,col=c("purple","blue","red","green4"), cex.axis=1,las=1, cex=1.2, 
     ylab="PC2",xlab="PC1")
legend(-0.28845761,0.44247002,c("S1","S2","S3","S4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("purple","blue","red","green4")) 
dev.off()
png(filename="pca_tpkm_log_scaleed_allPC.png", units="in", width=4, height=4 , pointsize=5, res=1200)
plot(t,pch=18,col=c("purple","blue","red","green4"), cex.axis=1,las=1, cex=1)
dev.off()

library(corrplot)
png(filename="Sample_distances_log_tpkm.png", units="in", width=4, height=4 , pointsize=5, res=1200)
par(mfrow=c(2,2))
t<-as.matrix(dist(t(exon), upper=T),method = "euclidean")
corrplot(t, is.corr=F, tl.col="black",method="color")
plot(hclust(dist(t(exon), method="euclidean")))
t<-as.matrix(dist(t(exon), upper=F,method = "maximum"))
corrplot(t, is.corr=F, tl.col="black",method="color")
plot(hclust(dist(t(exon), method="maximum")))
dev.off()

# log transform unnormalised
x<-c("S1_B1","S2_B1","S3_B1","S4_B1","S1_B2","S2_B2","S3_B2","S4_B2","S1_B3","S2_B3","S3_B3","S4_B3")
exon <- log2(rnaseq[, which(names(rnaseq) %in% x)])
is.na(exon) <- do.call(cbind,lapply(exon, is.infinite)) 
exon.pca <- prcomp(na.omit(exon),scale= F) 
t<-as.data.frame(exon.pca$rotation)
png(filename="pca_unnorm_log_scaled.png", units="in", width=2, height=2 , pointsize=5, res=1200)
plot(t[,2]~t[,1],pch=18,col=c("purple","blue","red","green4"), cex.axis=1,las=1, cex=1.2, 
     ylab="PC2",xlab="PC1")
legend(-0.2875761,0.08247002,c("S1","S2","S3","S4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("purple","blue","red","green4")) 
dev.off()
png(filename="pca_unnorm_log_scaleed_allPC.png", units="in", width=4, height=4 , pointsize=5, res=1200)
plot(t,pch=18,col=c("purple","blue","red","green4"), cex.axis=1,las=1, cex=1)
dev.off()
library(corrplot)
png(filename="Sample_distances_log_unnorm.png", units="in", width=4, height=4 , pointsize=5, res=1200)
par(mfrow=c(2,2))
t<-as.matrix(dist(t(exon), upper=T),method = "euclidean")
corrplot(t, is.corr=F, tl.col="black",method="color")
plot(hclust(dist(t(exon), method="euclidean")))
t<-as.matrix(dist(t(exon), upper=F,method = "maximum"))
corrplot(t, is.corr=F, tl.col="black",method="color")
plot(hclust(dist(t(exon), method="maximum")))
dev.off()


temp<-as.data.frame(rowMeans(subset(rnaseq, select = c(S1_B1,S1_B2,S1_B3)), na.rm = TRUE))
temp[,2]<-rowMeans(subset(rnaseq, select = c(S2_B1,S2_B2,S2_B3)), na.rm = TRUE)
temp[,3]<-rowMeans(subset(rnaseq, select = c(S3_B1,S3_B2,S3_B3)), na.rm = TRUE)
temp[,4]<-rowMeans(subset(rnaseq, select = c(S4_B1,S4_B2,S4_B3)), na.rm = TRUE)
png(filename="rawcount_histogram.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
par(mfrow=c(2,2))
hist(log(temp[,1], base=2), col="gray", las=1, cex.axis=0.8,xlab="log2counts",main="S1",xlim=c(-5,20))
box()
hist(log(temp[,2], base=2), col="gray", las=1, cex.axis=0.8,xlab="log2counts",main="S2",xlim=c(-5,20))
box()
hist(log(temp[,3], base=2), col="gray", las=1, cex.axis=0.8,xlab="log2counts",main="S3",xlim=c(-5,20))
box()
hist(log(temp[,4], base=2), col="gray", las=1, cex.axis=0.8,xlab="log2counts",main="S4",xlim=c(-5,20))
box()
dev.off()

png(filename="rawcount_density.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
d <- density(log(temp[,4], base=2)) # returns the density data 
plot(d, col="green4", las=1, lwd=1, main="")
d <- density(log(temp[,3], base=2)) # returns the density data 
lines(d,col="red",lwd=1)
d <- density(log(temp[,2], base=2)) # returns the density data 
lines(d,col="blue",lwd=1)
d <- density(log(temp[,1], base=2)) # returns the density data 
lines(d,col="purple",lwd=1)
legend(9.335815,0.04086374,c("S1","S2","S3","S4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("purple","blue","red","green4")) 
dev.off()

png(filename="rawcount_vs_tpkm.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
par(mfrow=c(2,2))
plot(log(temp[,1], base=2)~log(rnaseq$S1, base=2), las=1, pch=18, main="S1", xlab="TPKM (log2)",ylab="rawcounts (log2)")
usr <- par( "usr" )
a=round(cor(temp[,1],rnaseq$S1), digits=3)
text( usr[ 1 ], usr[ 4 ], a,adj = c( -1, 2 ))
plot(log(temp[,2], base=2)~log(rnaseq$S2, base=2), las=1, pch=18, main="S2",xlab="TPKM (log2)",ylab="rawcounts (log2)")
usr <- par( "usr" )
a=round(cor(temp[,2],rnaseq$S2), digits=3)
text( usr[ 1 ], usr[ 4 ], a,adj = c( -1, 2 ))
plot(log(temp[,3], base=2)~log(rnaseq$S3, base=2), las=1, pch=18, main="S3",xlab="TPKM (log2)",ylab="rawcounts (log2)")
usr <- par( "usr" )
a=round(cor(temp[,3],rnaseq$S3), digits=3)
text( usr[ 1 ], usr[ 4 ], a,adj = c( -1, 2 ))
plot(log(temp[,4], base=2)~log(rnaseq$S4, base=2), las=1, pch=18, main="S4",xlab="TPKM (log2)",ylab="rawcounts (log2)")
usr <- par( "usr" )
a=round(cor(temp[,1],rnaseq$S4), digits=3)
text( usr[ 1 ], usr[ 4 ], a,adj = c( -1, 2 ))
dev.off()
png(filename="tpkm_histogram.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
par(mfrow=c(2,2))
hist(log(rnaseq$S1, base=2), col="gray", las=1, cex.axis=0.8,xlab="TPKM",main="S1",xlim=c(-20,20))
box()
hist(log(rnaseq$S2, base=2), col="gray", las=1, cex.axis=0.8,xlab="TPKM",main="S2",xlim=c(-20,20))
box()
hist(log(rnaseq$S3, base=2), col="gray", las=1, cex.axis=0.8,xlab="TPKM",main="S3",xlim=c(-20,20))
box()
hist(log(rnaseq$S4, base=2), col="gray", las=1, cex.axis=0.8,xlab="TPKM",main="S4",xlim=c(-20,20))
box()
dev.off()

png(filename="tpkm_density.png", units="in", width=2, height=2, pointsize=4, res=1200)
par(font.axis=2,font.lab=2)
d <- density(log(rnaseq$S4, base=2)) # returns the density data 
plot(d, col="green4", las=1, lwd=1, main="")
d <- density(log(rnaseq$S3, base=2)) # returns the density data 
lines(d,col="red",lwd=1)
d <- density(log(rnaseq$S2, base=2)) # returns the density data 
lines(d,col="blue",lwd=1)
d <- density(log(rnaseq$S1, base=2)) # returns the density data 
lines(d,col="purple",lwd=1)
legend(3.003672,0.03176209,c("S1","S2","S3","S4"), lty=c(1,1), lwd=c(2.5,2.5),col=c("purple","blue","red","green4")) 
dev.off()

x<-c("GENEID","S1_B1_TPKM","S2_B1_TPKM","S3_B1_TPKM","S4_B1_TPKM","S1_B2_TPKM","S2_B2_TPKM",
     "S3_B2_TPKM","S4_B2_TPKM","S1_B3_TPKM","S2_B3_TPKM","S3_B3_TPKM","S4_B3_TPKM")

rnaseq_mapped_genes<-aggregate(. ~ GENEID, data = rnaseq[, which(names(rnaseq) %in% x)],sum,na.rm=T)
# mean bio_replicates
rnaseq_mapped_genes$S1<-rowMeans(subset(rnaseq_mapped_genes, select = c(S1_B1_TPKM,S1_B2_TPKM,S1_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$S2<-rowMeans(subset(rnaseq_mapped_genes, select = c(S2_B1_TPKM,S2_B2_TPKM,S2_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$S3<-rowMeans(subset(rnaseq_mapped_genes, select = c(S3_B1_TPKM,S3_B2_TPKM,S3_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$S4<-rowMeans(subset(rnaseq_mapped_genes, select = c(S4_B1_TPKM,S4_B2_TPKM,S4_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$Pop_Mean<-apply(rnaseq_mapped_genes[, which(names(rnaseq_mapped_genes) %in% c("S1", "S2", "S3","S4"))],1,mean) 
rnaseq_mapped_genes$Pop_Std_Dev<-apply(rnaseq_mapped_genes[, which(names(rnaseq_mapped_genes) %in% c("S1", "S2", "S3","S4"))],1,sd) 

rnaseq_mapped_genes$S1_zscore<-(rnaseq_mapped_genes$S1-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev
rnaseq_mapped_genes$S2_zscore<-(rnaseq_mapped_genes$S2-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev
rnaseq_mapped_genes$S3_zscore<-(rnaseq_mapped_genes$S3-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev
rnaseq_mapped_genes$S4_zscore<-(rnaseq_mapped_genes$S4-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev

# Fold changes
rnaseq_mapped_genes$FC_S2_S1<-rnaseq_mapped_genes$S2/rnaseq_mapped_genes$S1#S2/S1
rnaseq_mapped_genes$FC_S3_S1<-rnaseq_mapped_genes$S3/rnaseq_mapped_genes$S1#S3/S1
rnaseq_mapped_genes$FC_S4_S1<-rnaseq_mapped_genes$S4/rnaseq_mapped_genes$S1#S4/S1
rnaseq_mapped_genes$FC_S3_S2<-rnaseq_mapped_genes$S3/rnaseq_mapped_genes$S2#S3/S2
rnaseq_mapped_genes$FC_S4_S2<-rnaseq_mapped_genes$S4/rnaseq_mapped_genes$S2#S4/S2
rnaseq_mapped_genes$FC_S4_S3<-rnaseq_mapped_genes$S4/rnaseq_mapped_genes$S3#S4/S3


rnaseq_mapped_genes[,1] <- as.character(rnaseq_mapped_genes[,1])
png(filename="Foldchange_gene.png", units="in", width=3, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
boxplot(log(rnaseq_mapped_genes[,grepl( "FC" , names( rnaseq_mapped_genes ) )], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18)
dev.off()
png(filename="zscores_gene.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
boxplot(rnaseq_mapped_genes[,grepl( "zscore" , names( rnaseq_mapped_genes ) )],col="gray", las=1,ylab="Zscores", lty=1, pch=18, names=c("S1","S2","S3","S4"))
dev.off()

library(EnsDb.Mmusculus.v79)
keys <- keys(EnsDb.Mmusculus.v79)
columns(EnsDb.Mmusculus.v79)
anno<-select(EnsDb.Mmusculus.v79, keys=keys, columns=c("GENEID","SYMBOL","ENTREZID","GENEBIOTYPE",
                                                       "SEQNAME"),keytype="GENEID")
rnaseq_mapped_genes <- merge(rnaseq_mapped_genes,anno,by="GENEID")

temp<-as.data.frame(table(rnaseq$GENEID))
png(filename="histogram_transcriptspergene.png", units="in", width=2, height=1.5, pointsize=3, res=1200)
par(font.axis=2,font.lab=2)
par(mfrow=c(1,2))
hist(log(temp[,2], base=2), col="gray", las=1, cex.axis=0.8,xlab="transcriptspergene(log2)",main="")
t<-quantile(temp[,2])
text(3,20041.06,paste("Median=",t[3],sep=""), cex=1)
box()
boxplot(temp[,2], col="gray", las=1, cex.axis=0.8, pch=18,lty=1, ylab="transcriptspergene")
dev.off()

setwd(dird)
write.csv(rnaseq_mapped_genes, file="allcombined_gene_level_expression_mapped_to_gene_id.csv")
setwd(dirp)

x<-c("SYMBOL","S1_B1_TPKM","S2_B1_TPKM","S3_B1_TPKM","S4_B1_TPKM","S1_B2_TPKM","S2_B2_TPKM",
     "S3_B2_TPKM","S4_B2_TPKM","S1_B3_TPKM","S2_B3_TPKM","S3_B3_TPKM","S4_B3_TPKM")

rnaseq_mapped_genes<-aggregate(. ~ SYMBOL, data = rnaseq_mapped_genes[, which(names(rnaseq_mapped_genes) %in% x)],sum,na.rm=T)
# mean bio_replicates
rnaseq_mapped_genes$S1<-rowMeans(subset(rnaseq_mapped_genes, select = c(S1_B1_TPKM,S1_B2_TPKM,S1_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$S2<-rowMeans(subset(rnaseq_mapped_genes, select = c(S2_B1_TPKM,S2_B2_TPKM,S2_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$S3<-rowMeans(subset(rnaseq_mapped_genes, select = c(S3_B1_TPKM,S3_B2_TPKM,S3_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$S4<-rowMeans(subset(rnaseq_mapped_genes, select = c(S4_B1_TPKM,S4_B2_TPKM,S4_B3_TPKM)), na.rm = TRUE)
rnaseq_mapped_genes$Pop_Mean<-apply(rnaseq_mapped_genes[, which(names(rnaseq_mapped_genes) %in% c("S1", "S2", "S3","S4"))],1,mean) 
rnaseq_mapped_genes$Pop_Std_Dev<-apply(rnaseq_mapped_genes[, which(names(rnaseq_mapped_genes) %in% c("S1", "S2", "S3","S4"))],1,sd) 

rnaseq_mapped_genes$S1_zscore<-(rnaseq_mapped_genes$S1-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev
rnaseq_mapped_genes$S2_zscore<-(rnaseq_mapped_genes$S2-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev
rnaseq_mapped_genes$S3_zscore<-(rnaseq_mapped_genes$S3-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev
rnaseq_mapped_genes$S4_zscore<-(rnaseq_mapped_genes$S4-rnaseq_mapped_genes$Pop_Mean)/rnaseq_mapped_genes$Pop_Std_Dev

# Fold changes
rnaseq_mapped_genes$FC_S2_S1<-rnaseq_mapped_genes$S2/rnaseq_mapped_genes$S1#S2/S1
rnaseq_mapped_genes$FC_S3_S1<-rnaseq_mapped_genes$S3/rnaseq_mapped_genes$S1#S3/S1
rnaseq_mapped_genes$FC_S4_S1<-rnaseq_mapped_genes$S4/rnaseq_mapped_genes$S1#S4/S1
rnaseq_mapped_genes$FC_S3_S2<-rnaseq_mapped_genes$S3/rnaseq_mapped_genes$S2#S3/S2
rnaseq_mapped_genes$FC_S4_S2<-rnaseq_mapped_genes$S4/rnaseq_mapped_genes$S2#S4/S2
rnaseq_mapped_genes$FC_S4_S3<-rnaseq_mapped_genes$S4/rnaseq_mapped_genes$S3#S4/S3


rnaseq_mapped_genes[,1] <- as.character(rnaseq_mapped_genes[,1])
png(filename="Foldchange_genename.png", units="in", width=3, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
boxplot(log(rnaseq_mapped_genes[,grepl( "FC" , names( rnaseq_mapped_genes ) )], base=2),col="gray", las=2,ylab="Foldchange (log2 scale)", lty=1, pch=18)
dev.off()
png(filename="zscores_genename.png", units="in", width=2, height=2, pointsize=6, res=1200)
par(font.axis=2,font.lab=2)
boxplot(rnaseq_mapped_genes[,grepl( "zscore" , names( rnaseq_mapped_genes ) )],col="gray", las=1,ylab="Zscores", lty=1, pch=18, names=c("S1","S2","S3","S4"))
dev.off()

library(EnsDb.Mmusculus.v79)
keys <- keys(EnsDb.Mmusculus.v79)
columns(EnsDb.Mmusculus.v79)
anno<-select(EnsDb.Mmusculus.v79, keys=keys, columns=c("SYMBOL","ENTREZID","GENEBIOTYPE",
                                                       "SEQNAME"),keytype="GENEID")
rnaseq_mapped_genes <- merge(rnaseq_mapped_genes,anno,by="SYMBOL")

temp<-as.data.frame(table(rnaseq$GENEID))
png(filename="histogram_transcriptspergenename.png", units="in", width=2, height=1.5, pointsize=3, res=1200)
par(font.axis=2,font.lab=2)
par(mfrow=c(1,2))
hist(log(temp[,2], base=2), col="gray", las=1, cex.axis=0.8,xlab="transcriptspergene(log2)",main="")
t<-quantile(temp[,2])
text(3,20041.06,paste("Median=",t[3],sep=""), cex=1)
box()
boxplot(temp[,2], col="gray", las=1, cex.axis=0.8, pch=18,lty=1, ylab="transcriptspergene")
dev.off()

setwd(dird)
write.csv(rnaseq_mapped_genes, file="allcombined_gene_level_expression_mapped_to_gene_name.csv")
setwd(dirp)