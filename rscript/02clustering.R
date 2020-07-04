load(file = "./eSet.NT.RData")
sampleinfo <- read.csv("/datasets/work/hb-diab-cfdna/work/Data/level3/TCGA_mRNA/PRAD/sampleinfo_allsample.csv", 
                       stringsAsFactors = F )
eSet <- eSet.NT[,53:551]
eSet <- eSet[,colnames(eSet) %in% sampleinfo$ID]
sampleinfo <- sampleinfo[sampleinfo$ID  %in% colnames(eSet),]
sampleinfo <- sampleinfo[match( colnames(eSet), sampleinfo$ID),]


library(ConsensusClusterPlus)

eSet_cluster = sweep(eSet,1, apply(eSet,1,median,na.rm=T)) #median center dataset

title= c("./Cluster/")
eSet_cluster <- as.matrix(eSet_cluster)
results = ConsensusClusterPlus(eSet_cluster, maxK=6,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",
                               seed=71279,plot="pdf")
results[[2]][["consensusMatrix"]][1:5,1:5]
results[[4]][["consensusTree"]]
results[[4]][["consensusClass"]][1:5]
consensus.Mat <- results[[2]][["consensusMatrix"]]
save(consensus.Mat ,  file = "./Cluster/02ConsensusCluster.RData")

icl = calcICL(results,title=title,plot="pdf")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]

#clustering k = 4 is ideal subgroup for downstream analysis
subtype <- results[[2]][["consensusClass"]]
subtype <- paste("Subtype", subtype , sep = "")
sampleinfo$subtype <- subtype

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Dark2")
pal.batch <- brewer.pal(8, "")
pal.subtype <- brewer.pal(8, "Dark2")

library(limma)
pdf("./plots/02PCA_cluster.pdf", width = 5, height = 5)
plotMDS(eSet, top=1000, gene.selection="common", 
        col=pal.subtype[factor(sampleinfo$subtype)], 
        pch = 19, cex = 0.8,
        main = "Clustering by m6A regulators")
legend("topleft", legend=levels(factor(sampleinfo$subtype)), text.col=pal.subtype,
       bg="white", cex=1)
dev.off()

save(eSet, sampleinfo, file ="./eSet_regulator.RData")


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
pdf("./plots/02sigclust.pdf", width = 5, height = 5)
sigclust1 <- sigclustTest(eSet,group, nsim=500, nrep=1, icovest=3)
sigclust2 <- sigclustTest(eSet_cluster,group, nsim=500, nrep=1, icovest=3)
dev.off()
#save(group, eSet, eSet_cluster, file = "./Subtype_cluster/02CancerSubtype_sigclust.RData")



######Silhouette width######
library(CancerSubtypes)
library(RColorBrewer)

group.s <- as.numeric(substr(group[order(group)], start = 8L, stop = 8L))
consensus.Mat.s <- consensus.Mat[order(group),order(group)]
sil=silhouette_SimilarityMatrix(group.s, consensus.Mat.s)
pdf("./plots/02SilhouettePlot.pdf", width = 5, height = 4)
plot(sil , col = pal.subtype[1:2])
dev.off()


#####heatmap #####
library(pheatmap)
pal.subtype <- pal.subtype[1:2]

annotation_col = data.frame(
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet)

names(pal.subtype) <- levels(as.factor(sampleinfo$subtype))
breaksList = seq(-3, 3, by = 0.01)

pdf("./plots/02Heatmap_cluster.pdf", width = 6, height = 6)
pheatmap(eSet, 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = T,
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
load(file ="./eSet_regulator.RData")

writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

idx.sub1 <- which(sampleinfo$subtype == "Subtype1")
idx.sub2 <- which(sampleinfo$subtype == "Subtype2")
eSet1.cor <- as.data.frame(t(eSet[,idx.sub1]))
eSet2.cor <- as.data.frame(t(eSet[,idx.sub2]))
eSet1.cor <- eSet1.cor[,match(all.regulators, colnames(eSet1.cor))]
eSet2.cor <- eSet2.cor[,match(all.regulators, colnames(eSet2.cor))]

library(pheatmap)
#col.cor<- colorRampPalette(c("blue", "white", "red"))(20)
breaksList = seq(-1, 1, by = 0.1)

pdf("./plots/02Heatmap_cor.pdf", width = 5, height = 5)
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
dev.off()
