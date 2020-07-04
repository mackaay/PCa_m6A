rm(list = ls())

library(limma)
library(RColorBrewer)
pal.grade <- rev(brewer.pal(6, "Accent"))
pal.subtype <- brewer.pal(8, "Dark2")

load(file ="./eSet_regulator_FPKM.RData")
FPKM <- read.csv("./eSet.csv", header = T, stringsAsFactors = F)
FPKM <- FPKM[!duplicated(FPKM$id),]
rownames(FPKM) <- FPKM$id
FPKM <- FPKM[,2:551]
Tumor.idx <- which(!grepl("11", colnames(FPKM)))
FPKM <- FPKM[,Tumor.idx]
eSet <- log2(FPKM+1)

#remove metatastasis sample 
eSet <- eSet[,-377]
sampleinfo <- sampleinfo[-377,]

save(eSet, sampleinfo, file = "./eSet_TumorALL.RData")

####DEG#####
cellType <- factor(sampleinfo$subtype)
table(cellType)

design <- model.matrix(~0+cellType, data=sampleinfo)
colnames(design) <- c(levels(cellType))

fit <- lmFit(eSet, design)

contMatrix <- makeContrasts(Subtype1-Subtype2,
                            Subtype1-Subtype3,
                            Subtype2-Subtype3,
                            levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

DEG1 <- topTable(fit2, num=Inf, coef=1)
DEG1 <- DEG1[which(abs(DEG1$logFC) >2 & DEG1$adj.P.Val <0.05),]
DEG2 <- topTable(fit2, num=Inf, coef=2)
DEG2 <- DEG2[which(abs(DEG2$logFC) >2 & DEG2$adj.P.Val <0.05),]
DEG3 <- topTable(fit2, num=Inf, coef=3)
DEG3 <- DEG3[which(abs(DEG3$logFC) >2 & DEG3$adj.P.Val <0.05),]

DEG.name <- unique(c(rownames(DEG1), rownames(DEG2), rownames(DEG3)))

write.csv(DEG.name, "./DEG_subtype.csv")

library(pheatmap)

annotation_col = data.frame(
  GleasonRisk = sampleinfo$gleason_level,
  #Batch = sampleinfo$batch,
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(eSet[DEG.name,])
breaksList = seq(-3, 3, by = 0.01)

pdf("./plots/03Heatmap_DEG_FPKM.pdf", width = 7, height = 5)
pheatmap(eSet[DEG.name,order(sampleinfo$subtype)], 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(Subtype1 = pal.subtype[1], Subtype2=pal.subtype[2], 
                       Subtype3 = pal.subtype[3]),
           GleasonRisk = c(low = pal.grade[1], intermediate =pal.grade[2], high=pal.grade[3])),
         main = "DEG of PRAD"
)

pheatmap(eSet[DEG.name,], 
         scale = "row",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = F,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = T,
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,
         annotation_legend = T,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Batch = c(GSE85133 =pal.batch[1], GSE136661 =pal.batch[2]),
           Subtype = c(Subtype1 = pal.subtype[1], Subtype2=pal.subtype[2], 
                       Subtype3 = pal.subtype[3]),
           GleasonRisk = c(low = pal.grade[1], intermediate =pal.grade[2], high=pal.grade[3])),
         main = "DEG of PRAD"
)
dev.off()



#####correlation of m6A score and genes #####
cor.df <- as.data.frame(t(eSet))
cor.df$m6Ascore <- sampleinfo$m6Ascore
cor <- c()
pvalue <- c()
for (i in 1:(ncol(cor.df)-1)) {
  cor.test <- cor.test(cor.df[,i], cor.df$m6Ascore)
  cor <- c(cor, cor.test$estimate)
  pvalue <- c(pvalue, cor.test$p.value)
}
cor.df.p <- data.frame(cor= cor, pvalue= pvalue)
rownames(cor.df.p) <- rownames(eSet)
cor.df.p[which(cor.df.p$pvalue < 0.01 & abs(cor.df.p$cor) >0.7),]

writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)
ttest.df <- read.csv("./ttest_NT.csv", header = T, stringsAsFactors = F)
m6A.reg.NTsig <- ttest.df[ttest.df$pval< 0.05,]$gene
cor.df.p[all.regulators,]

mydata <- cor.df[,all.regulators]
mydata$m6Ascore <- sampleinfo$m6Ascore
cormat<-signif(cor(mydata),2)
cormat

pdf("./plots/03Heatmap_m6Ascore_corr.pdf", width = 6, height = 6)
col<- colorRampPalette(c("blue", "white", "red"))
breaksList = seq(-1, 1, by = 0.01)
pheatmap(cormat, 
         color = col(length(breaksList)), 
         breaks = breaksList)
dev.off()



#####m6A with Gleason and Stemness index######
load(file ="./eSet_regulator_FPKM.RData")
library(limma)
library(RColorBrewer)
pal.grade <- rev(brewer.pal(6, "Accent"))
pal.subtype <- brewer.pal(8, "Dark2")
pal.gleason <- brewer.pal(5, "Set1")

#stemness score
library(grid)
library(ggplot2)
pdf("./plots/03Scatter_cor_m6Ascore&Stemness.pdf", width = 4, height = 4)
m6Ascore <- sampleinfo$m6Ascore
StemIndex <- sampleinfo$StemIndex
grob3 = grobTree(textGrob(paste("Pearson Correlation:", 
                                round(cor(m6Ascore, StemIndex), 4) ), 
                          x =0, y = 0.1, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))
grob4 = grobTree(textGrob(paste("p = ", cor.test(m6Ascore, StemIndex)$p.value ), 
                          x =0, y = 0.05, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))

ggplot(sampleinfo, aes(x=m6Ascore, y=StemIndex)) + 
  geom_point() + 
  ggtitle("m6Ascore vs Gleason score") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "m6Ascore", 
                     #limits = c(-2, 4), 
                     #breaks = seq(-4, 4, 1)
                     ) + 
  scale_y_continuous(name = "Stemness Index", 
                     #limits = c(-4, 4), 
                     #breaks = seq(-4, 4, 1)
                     ) + 
  annotation_custom(grob3) + 
  annotation_custom(grob4) + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_blank(), 
        axis.line = element_line(color="black"), 
        axis.line.x = element_line(color="black"))
dev.off()



##Gleason
library(grid)
library(ggplot2)
pdf("./plots/03Scatter_cor_m6Ascore&Gleason.pdf", width = 4, height = 4)
m6Ascore <- sampleinfo$m6Ascore
gleason_score <- sampleinfo$gleason_score
grob3 = grobTree(textGrob(paste("Pearson Correlation:", 
                                round(cor(m6Ascore, gleason_score), 4) ), 
                          x =0, y = 0.1, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))
grob4 = grobTree(textGrob(paste("p = ", cor.test(m6Ascore, gleason_score)$p.value ), 
                          x =0, y = 0.05, hjust = 0, 
                          gp = gpar(col = "red", fontsize = 10, fontface = "bold")))

ggplot(sampleinfo, aes(x=m6Ascore, y=gleason_score)) + 
  geom_point() + 
  ggtitle("m6Ascore vs Gleason score") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "m6Ascore", 
                     #limits = c(-2, 4), 
                     #breaks = seq(-4, 4, 1)
  ) + 
  scale_y_continuous(name = "Gleason Score", 
                     #limits = c(-4, 4), 
                     #breaks = seq(-4, 4, 1)
  ) + 
  annotation_custom(grob3) + 
  annotation_custom(grob4) + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_blank(), 
        axis.line = element_line(color="black"), 
        axis.line.x = element_line(color="black"))
dev.off()

pdf("./plots/03Boxplot_m6A&Gleason.pdf", width = 5, height = 4)
boxplot(m6Ascore ~ gleason_score, data = sampleinfo, 
        col = pal.gleason, notch = F,
        ylim = c(-10, 16), 
        main = "m6A score between Gleason Score"
        )
dev.off()

summary(aov(sampleinfo$m6Ascore ~ sampleinfo$gleason_score))
#p = 8.07e-10
TukeyHSD(aov(sampleinfo$m6Ascore ~ as.factor(sampleinfo$gleason_score)))
