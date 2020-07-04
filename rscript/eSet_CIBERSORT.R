load(file = "./TPM.RData")
load(file ="./eSet_regulator.RData")

colnames(TPM) <- substr(colnames(TPM), start = 1L, stop = 15L)
TPM <- TPM[,match(sampleinfo$ID, colnames(TPM))]

write.csv(TPM, file = "./TPM.csv")
write.csv(TPM[,1:235], file = "./TPM_subtype1.csv")
write.csv(TPM[,235:495], file = "./TPM_subtype2.csv")

library(RColorBrewer)
pal.subtype <- brewer.pal(8, "Dark2")
library(limma)
pdf("./plots/01PCA_tumor.pdf", width = 5, height = 5)
plotMDS(log2(TPM+1), top=1000, gene.selection="common", 
        col=pal.subtype[factor(sampleinfo$subtype)], 
        pch = 19, cex = 0.8,
        main = "")
legend("topleft", legend=levels(factor(sampleinfo$subtype)), text.col=pal.subtype,
       bg="white", cex=1)
dev.off()
