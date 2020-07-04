writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

count <- read.delim("/datasets/work/hb-diab-cfdna/work/Data/level3/TCGA_mRNA/PRAD/mRNA.txt", 
                    sep = "\t", stringsAsFactors = F)
count <- count[!duplicated(count$id),]
rownames(count) <- count$id
count <- count[,2:552]
gene_length <- read.csv("/datasets/work/hb-diab-cfdna/work/Data/level3/TCGA_mRNA/gene_length.csv",
                        stringsAsFactors = F)
gene_length <- gene_length[gene_length$Geneid %in% rownames(count),]
count <- count[rownames(count) %in% gene_length$Geneid,]
gene_length <- gene_length[match(rownames(count), gene_length$Geneid),]
gene_kilobase <- as.numeric(gene_length$Length)/1000

RPK <- count/gene_kilobase
coverage <- colSums(RPK[,1:ncol(RPK)])/1000000
RPK <- t(RPK)
TPM <- t(RPK/coverage)
TPM <- round(TPM, 4)
colSums(TPM)
TPM <- as.data.frame(TPM)

save(TPM , file = "./TPM.RData")

count.reg <- TPM[rownames(TPM) %in% all.regulators,]
colnames(count.reg) <- substr(colnames(count.reg), start = 1L, stop = 15L)

eSet <- log2(count.reg+1)

sampleinfo <- read.csv("/datasets/work/hb-diab-cfdna/work/Data/level3/TCGA_mRNA/PRAD/sampleinfo_allsample.csv", 
                        stringsAsFactors = F )
sampleinfo.NT <- data.frame(id = colnames(eSet), 
                            CellType = c(rep("Normal", 52), rep("Tumor", 499)))

library(limma)
library(RColorBrewer)
library(scales)
pal.normal <- brewer.pal(8, "Set1")

pdf("./plots/01PCA_normal.pdf", height = 5, width = 5)
plotMDS(eSet, top=1000, gene.selection="common", pch = 19,
        col=alpha(pal.normal[factor(sampleinfo.NT$CellType)], 0.5))
legend("bottomright", legend=levels(factor(sampleinfo.NT$CellType)), text.col=pal.normal,
       bg="white", cex=1)
dev.off()
eSet.NT <- eSet
save(eSet.NT, sampleinfo.NT, file = "./eSet.NT.RData")

#######t test between tumor and normal#######
pval.NT <- c()
for (i in 1:nrow(eSet.NT)) {
  tmp <- round(t.test(eSet.NT[i, 1:52], eSet.NT[i,53:551])$p.value, 3)
  pval.NT <- c(pval.NT, tmp)
}
ttest.df <- data.frame(gene = rownames(eSet.NT), pval = pval.NT)


######violin plot ######
violinplot <- as.data.frame(t(eSet))
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
                         CellType = CellType)

library(ggplot2)
library(tidyverse)
pdf("./plots/01Violinplot_normal.pdf", width = 10, height = 7)
violinplot %>%
  ggplot(aes(x= CellType, y = expression, fill = CellType))+
  geom_boxplot(outlier.size = 0)+
  #geom_jitter(aes(color = response),width = 0.3, alpha = 0.7)+
  xlab("type")+
  scale_fill_manual(values=pal.normal)+
  scale_color_manual(values=pal.normal)+
  facet_wrap(~gene)+
  ylim(c(0,12))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("")
dev.off()




