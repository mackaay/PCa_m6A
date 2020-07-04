FPKM <- read.csv("./eSet.csv", header = T, stringsAsFactors = F)
FPKM <- FPKM[!duplicated(FPKM$id),]
rownames(FPKM) <- FPKM$id
FPKM <- FPKM[,2:551]

writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

count.reg <- FPKM[rownames(FPKM) %in% all.regulators,]
colnames(count.reg) <- substr(colnames(count.reg), start = 1L, stop = 15L)

eSet <- log2(count.reg+1)
grep("11", colnames(eSet))

sampleinfo <- read.csv("./sampleinfo_allsample.csv", 
                       stringsAsFactors = F )
sampleinfo.NT <- data.frame(id = colnames(eSet), 
                            CellType = c(rep("Tumor", 550)), 
                            stringsAsFactors = F)
sampleinfo.NT$CellType[grep("11", colnames(eSet))] <- "Normal"


library(limma)
library(RColorBrewer)
library(scales)
pal.normal <- brewer.pal(8, "Set1")

pdf("./plots/01PCA_normal_FPKM.pdf", height = 5, width = 5)
plotMDS(eSet, top=1000, gene.selection="common", pch = 19,
        col=alpha(pal.normal[factor(sampleinfo.NT$CellType)], 0.5))
legend("topright", legend=levels(factor(sampleinfo.NT$CellType)), text.col=pal.normal,
       bg="white", cex=1)
dev.off()

#######t test between tumor and normal#######
NT.idx <- grep("11", colnames(eSet))
Tumor.idx <- which(!grepl("11", colnames(eSet)))
eSet.NT <- eSet
pval.NT <- c()
for (i in 1:nrow(eSet.NT)) {
  tmp <- round(t.test(eSet.NT[i, NT.idx], eSet.NT[i,Tumor.idx])$p.value, 3)
  pval.NT <- c(pval.NT, tmp)
}
ttest.df <- data.frame(gene = rownames(eSet.NT), pval = pval.NT)
write.csv(ttest.df, file = "./ttest_NT.csv", row.names = F)


######violin plot ######
eSet.NT <- eSet.NT[match(all.regulators, rownames(eSet.NT)),]
violinplot <- as.data.frame(t(eSet.NT))
violinplot$CellType <- sampleinfo.NT$CellType
CellType <- as.character(rep(violinplot$CellType, 21))
tmp <- c()
for (i in 1:21) {
  tmp1 <- violinplot[,i]
  tmp <- c(tmp, tmp1)
}

temp <- c()
for (i in 1:21) {
  temp1 <- rep(colnames(violinplot)[i],nrow(violinplot))
  temp <- c(temp, temp1)
}

violinplot <- data.frame(expression = tmp , 
                         gene = temp, 
                         CellType = CellType, 
                         stringsAsFactors = F)
violinplot$gene <- factor(violinplot$gene, levels = all.regulators)

pdf("./plots/01Boxplot_NT_FPKM.pdf", width = 10, height = 5)
boxplot(expression ~ CellType+gene, data = violinplot,
        #at = c(1:3, 5:7), 
        col = pal.normal[1:2],
        #names = c("", "A", "", "", "B", ""), 
        xaxs = FALSE, las=2,
        ylim = c(0,17),
        cex.lab=0.01, cex.axis=1,
        main= "m6A Regulators Expression in TCGA")
legend("bottomleft", fill = pal.normal[1:2], 
       legend = c("Normal", "Tumor"), cex = 0.8,
       horiz = F)
dev.off()


save(eSet.NT, sampleinfo.NT, file = "./eSet.NT_FPKM.RData")



library(ggplot2)
library(tidyverse)
pdf("./plots/01Violinplot_normal_FPKM.pdf", width = 10, height = 7)
violinplot %>%
  ggplot(aes(x= CellType, y = expression, fill = CellType))+
  geom_boxplot(outlier.size = 0)+
  #geom_jitter(aes(color = response),width = 0.3, alpha = 0.7)+
  xlab("type")+
  scale_fill_manual(values=pal.normal)+
  scale_color_manual(values=pal.normal)+
  facet_wrap(~gene)+
  #ylim(c(0,12))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("")
dev.off()

