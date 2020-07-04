load(file = "./eSet.NT_FPKM.RData")
sampleinfo <- read.csv("./sampleinfo_allsample.csv", 
                       stringsAsFactors = F )
rownames(sampleinfo) <-gsub("-", ".", sampleinfo$ID)
sampleinfo <- sampleinfo[match(colnames(eSet.NT), rownames(sampleinfo) ),]
Tumor.idx <- which(!grepl("11", colnames(eSet.NT)))
eSet <- eSet.NT[,Tumor.idx]
sampleinfo <- sampleinfo[Tumor.idx,]

sampleinfo$gleason_level <- sampleinfo$gleason_score
sampleinfo$gleason_level <- gsub(6, "low",sampleinfo$gleason_level)
sampleinfo$gleason_level <- gsub(7, "intermediate",sampleinfo$gleason_level)
sampleinfo$gleason_level <- gsub(8, "high",sampleinfo$gleason_level)
sampleinfo$gleason_level <- gsub(9, "high",sampleinfo$gleason_level)
sampleinfo$gleason_level <- gsub(10, "high",sampleinfo$gleason_level)


#####CC clustering#####
library(ConsensusClusterPlus)
eSet_cluster = sweep(eSet,1, apply(eSet,1,median,na.rm=T)) #median center dataset

title= c("./Cluster_FPKM/")
eSet_cluster <- as.matrix(eSet_cluster)
results = ConsensusClusterPlus(eSet_cluster, maxK=6,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",
                               seed=71279,plot="pdf")
results[[3]][["consensusMatrix"]][1:5,1:5]
results[[4]][["consensusTree"]]
results[[4]][["consensusClass"]][1:5]
consensus.Mat <- results[[3]][["consensusMatrix"]]
# k = 3 for CC
save(consensus.Mat ,  file = "./Cluster_FPKM/02ConsensusCluster.RData")

icl = calcICL(results,title=title,plot="pdf")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]

#clustering k = 3 is ideal subgroup for downstream analysis
subtype <- results[[3]][["consensusClass"]]
subtype <- paste("Subtype", subtype , sep = "")
sampleinfo$subtype <- subtype

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Dark2")

library(limma)
pdf("./plots/02PCA_cluster_FPKM.pdf", width = 5, height = 5)
plotMDS(eSet, top=1000, gene.selection="common", 
        col=pal.subtype[factor(sampleinfo$subtype)], 
        pch = 19, cex = 0.8,
        main = "Clustering by m6A regulators")
legend("topright", legend=levels(factor(sampleinfo$subtype)), text.col=pal.subtype,
       bg="white", cex=1)

plotMDS(eSet, top=1000, gene.selection="common", 
        col=pal.grade[factor(sampleinfo$gleason_level)], 
        pch = 19, cex = 0.8,
        main = "Clustering by m6A regulators")
legend("topright", legend=levels(factor(sampleinfo$gleason_level)), text.col=pal.grade,
       bg="white", cex=1)

dev.off()

###PCA score#####
m6A.reg.ov <- c("HNRNPA2B1", "FMR1", "METTL14", "KIAA1429","YTHDF1")
ttest.df <- read.csv("./ttest_NT.csv", header = T, stringsAsFactors = F)
m6A.reg.NTsig <- ttest.df[ttest.df$pval< 0.05,]$gene

library(limma)
pdf("./plots/02PCA_m6Ascore.ov_FPKM.pdf", width = 5, height = 5)
plotMDS(eSet[m6A.reg.ov,], top=1000, gene.selection="common", 
        col=pal.subtype[factor(sampleinfo$subtype)], 
        pch = 19, cex = 0.8,
        main = "Clustering by m6A regulators")
legend("topright", legend=levels(factor(sampleinfo$subtype)), text.col=pal.subtype,
       bg="white", cex=1)
dev.off()

pca.out = prcomp(t(eSet[m6A.reg.ov,]), scale=TRUE)
pca.out$x
pca.score <- pca.out$x[,1]+pca.out$x[,2]
pca.score 
sampleinfo$m6Ascore <- pca.score

summary(aov(m6Ascore~ subtype, data = sampleinfo))
TukeyHSD(aov(m6Ascore~ subtype, data = sampleinfo))

cor.test(sampleinfo$psa_value, sampleinfo$m6Ascore, 
         method = c("pearson"), use="complete.obs")
cor.test(sampleinfo$gleason_score, sampleinfo$m6Ascore, 
         method = c("pearson"), use="complete.obs")

save(eSet, sampleinfo, file ="./eSet_regulator_FPKM.RData")

##
library(RColorBrewer)
pal.grade <- rev(brewer.pal(6, "Accent"))
pal.subtype <- brewer.pal(8, "Dark2")
pal.gleason <- brewer.pal(8, "Set1")
pdf("./plots/02Boxplot_m6Ascore_subtypes.pdf", width= 5, height = 5)
boxplot(m6Ascore~subtype, data=sampleinfo, notch=TRUE,
        col=pal.subtype,
        ylim = c(-8, 15),
        main="m6A score between subtypes", 
        xlab="ANOVA, p = 3.47e-06")
dev.off()

plot(sampleinfo$m6Ascore, sampleinfo$gleason_score)

pdf("./plots/02StackBarplot_subtype.pdf", width = 7, height = 6)
barplot(table(sampleinfo$gleason_score,sampleinfo$subtype ), 
        col = pal.gleason, 
        legend = rownames(table(sampleinfo$gleason_score,sampleinfo$subtype )), 
        main = "Gleason Score between Subtypes", 
        ylab = "Frequency")

dev.off()

######CancerSubtype SigClust######
library(CancerSubtypes)
#data(GeneExp)
#data(miRNAExp)
#data(time)
#data(status)
#GBM=list(GeneExp=GeneExp,miRNAExp=miRNAExp)
#result=ExecuteSNF(GBM, clusterNum=3, K=20, alpha=0.5, t=20)
#group=result$group
#sigclust1=sigclustTest(miRNAExp,group, nsim=500, nrep=1, icovest=3)
#sigclust2=sigclustTest(miRNAExp,group, nsim=1000, nrep=1, icovest=1)

group <- sampleinfo$subtype
pdf("./plots/02sigclust_FPKM.pdf", width = 5, height = 5)
sigclust1 <- sigclustTest(eSet,group, nsim=500, nrep=1, icovest=3)
sigclust2 <- sigclustTest(eSet_cluster,group, nsim=500, nrep=1, icovest=3)
dev.off()
#save(group, eSet, eSet_cluster, file = "./Subtype_cluster/02CancerSubtype_sigclust.RData")



######Silhouette width######
library(CancerSubtypes)
library(RColorBrewer)

group.s <- as.numeric(substr(group, start = 8L, stop = 8L))[order(group)]
consensus.Mat.s <- consensus.Mat[order(group),order(group)]
sil=silhouette_SimilarityMatrix(group.s, consensus.Mat.s)
pdf("./plots/02SilhouettePlot_FPKM.pdf", width = 5, height = 4)
plot(sil , col = pal.subtype[1:3])
dev.off()


#####heatmap #####
library(pheatmap)
pal.subtype <- pal.subtype[1:3]
eSet <- eSet[all.regulators,]
annotation_col = data.frame(
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)
sampleinfo <- sampleinfo[order(sampleinfo$subtype),]
eSet <- eSet[,order(sampleinfo$subtype)]
sampleinfo$ID <- gsub("-", ".", sampleinfo$ID)
eSet <- eSet[,match(sampleinfo$ID, colnames(eSet))]

names(pal.subtype) <- levels(as.factor(sampleinfo$subtype))
breaksList = seq(-3, 3, by = 0.01)

pdf("./plots/02Heatmap_CC_FPKM.pdf", width = 6, height = 4.5)
pheatmap(eSet, 
         scale = "row",
         #kmeans_k = 3,
         #clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         #clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = T,
         annotation_legend = T,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList,
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Time = c("white", "firebrick"),
           # subtype = c(subtype1 = "#A6CEE3", subtype2 = "#1F78B4",
           #            subtype3 = "#B2DF8A", subtype4 = "#33A02C"),
           Subtype = pal.subtype
         ),
         main = "m6A regulators"
)
dev.off()


####correlation ######
load(file ="./eSet_regulator_FPKM.RData")

writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

idx.sub1 <- which(sampleinfo$subtype == "Subtype1")
idx.sub2 <- which(sampleinfo$subtype == "Subtype2")
idx.sub3 <- which(sampleinfo$subtype == "Subtype3")
eSet1.cor <- as.data.frame(t(eSet[,idx.sub1]))
eSet2.cor <- as.data.frame(t(eSet[,idx.sub2]))
eSet3.cor <- as.data.frame(t(eSet[,idx.sub3]))
eSet1.cor <- eSet1.cor[,match(all.regulators, colnames(eSet1.cor))]
eSet2.cor <- eSet2.cor[,match(all.regulators, colnames(eSet2.cor))]
eSet3.cor <- eSet3.cor[,match(all.regulators, colnames(eSet3.cor))]

library(pheatmap)
#col.cor<- colorRampPalette(c("blue", "white", "red"))(20)
breaksList = seq(-1, 1, by = 0.1)

pdf("./plots/02Heatmap_cor_FPKM.pdf", width = 5, height = 5)
cormat<-signif(cor(eSet1.cor),2)
cormat
pheatmap(cormat,  scale = "none", 
         show_rownames = T, show_colnames = T,
         fontsize_row = 7, fontsize_col = 7,
         cluster_cols = F,cluster_rows  = F,
         color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),  breaks = breaksList)


cormat<-signif(cor(eSet2.cor),2)
cormat
pheatmap(cormat,  scale = "none", 
         show_rownames = T, show_colnames = T,
         fontsize_row = 7, fontsize_col = 7,
         cluster_cols = F,cluster_rows  = F,
         color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),  breaks = breaksList)

cormat<-signif(cor(eSet3.cor),2)
cormat
pheatmap(cormat,  scale = "none", 
         show_rownames = T, show_colnames = T,
         fontsize_row = 7, fontsize_col = 7,
         cluster_cols = F,cluster_rows  = F,
         color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),  breaks = breaksList)
dev.off()



#### clinical sample informatin #####
library(RColorBrewer)
pal.grade <- rev(brewer.pal(6, "Accent"))
pal.gleason <- brewer.pal(6,"Set1")
pal.subtype <- brewer.pal(8, "Dark2")
which(sampleinfo$psa_value > 10 )
pdf("./plots/02Boxplot_subtypes_Gleason.pdf", width= 5, height = 5)
boxplot(gleason_score~subtype, data=sampleinfo, notch=F,
        col=pal.subtype,
        main="Gleason Score between subtypes", 
        xlab="")
dev.off()
