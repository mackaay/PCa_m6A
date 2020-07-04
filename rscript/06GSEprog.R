writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

ttest.df <- read.csv("./ttest_NT.csv", header = T, stringsAsFactors = F)
m6A.reg.NTsig <- ttest.df[ttest.df$pval< 0.05,]$gene


######GSE147493####
library(GEOquery)
GSEset <- getGEO("GSE147493", destdir = "./GSE147493_Garraway/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pheno_GSE147493 <- pData(phenoData(GSEset[[1]]))[,c(2,1, 40,41)]
pheno_GSE147493

eSet_GSE147493 <- read.csv("./GSE147493_Garraway/GSE147493_TMM_normalized_data.csv", 
                           stringsAsFactors = F)
eSet_GSE147493 <- eSet_GSE147493[!duplicated(eSet_GSE147493$Symbol),]
rownames(eSet_GSE147493) <- eSet_GSE147493$Symbol
eSet_GSE147493 <- eSet_GSE147493[,-1]

eSet_GSE147493.m6A <- eSet_GSE147493[all.regulators,]
eSet_GSE147493.m6A <- eSet_GSE147493.m6A[-16,]# remove NAs

library(limma)
library(RColorBrewer)
library(scales)
pal.prog <- brewer.pal(8, "Accent")

pdf("./plots/06PCA_GSE147493_Garraway.pdf", height = 5, width = 5)
plotMDS(eSet_GSE147493.m6A, top=1000, gene.selection="common", pch = 19,
        col=alpha(pal.prog[factor(pheno_GSE147493$`metastasis:ch1`)], 0.5))
legend("topright", legend=levels(factor(pheno_GSE147493$`metastasis:ch1`)), text.col=pal.prog,
       bg="white", cex=1)
dev.off()

#ttest
pval <- c()
for (i in 1:nrow(eSet_GSE147493.m6A)) {
  tmp <- round(t.test(eSet_GSE147493.m6A[i, 1:37], eSet_GSE147493.m6A[i,38:99])$p.value, 3)
  pval <- c(pval, tmp)
}
ttest.df <- data.frame(gene = rownames(eSet_GSE147493.m6A), pval = pval)

write.csv(ttest.df, file = "./ttest_GSE147493.csv")
#HNRNPA2B1, FMR1, METTL14 0.052, KIAA1429 0.054,YTHDF1 0.076
meta.regulators <- c("HNRNPA2B1", "FMR1", "METTL14", "KIAA1429","YTHDF1", 
                     "ALKBH5", "HNRNPC")
m6A.reg.ov <- intersect(meta.regulators, m6A.reg.NTsig)

pdf("./plots/06PCA_GSE147493_Garraway_m6A.ov.reg.pdf", height = 5, width = 5)
plotMDS(eSet_GSE147493.m6A[m6A.reg.ov,], top=1000, gene.selection="common", pch = 19,
        col=alpha(pal.prog[factor(pheno_GSE147493$`metastasis:ch1`)], 0.5))
legend("topright", legend=levels(factor(pheno_GSE147493$`metastasis:ch1`)), text.col=pal.prog,
       bg="white", cex=1)
dev.off()


#boxplot
eSet_GSE147493.m6A <- eSet_GSE147493.m6A[match(all.regulators, rownames(eSet_GSE147493.m6A)),]
violinplot <- as.data.frame(t(eSet_GSE147493.m6A))
violinplot$CellType <- pheno_GSE147493$`metastasis:ch1`
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

pdf("./plots/06Boxplot_GSE147493.pdf", width = 10, height = 5)

boxplot(expression ~ CellType+gene, data = violinplot,
        #at = c(1:3, 5:7), 
        col = pal.prog[7:8],
        #names = c("", "A", "", "", "B", ""), 
        xaxs = FALSE, las=2,
        ylim = c(2,19),
        beside=TRUE,
        cex.lab=0.01, cex.axis=1,
        main= "m6A Regulators Expression in GSE147493")
legend("bottomleft", fill = pal.prog[7:8], 
       legend = c("Metastasis", "Non-Metastasis"), cex = 0.8,
       horiz = F)
dev.off()


###m6A score 
pca.out = prcomp(t(eSet_GSE147493.m6A[m6A.reg.ov,]), scale=TRUE)
pca.out$x
pca.score <- pca.out$x[,1]
pca.score 
pheno_GSE147493$m6Ascore <- pca.score

t.test(m6Ascore~ `metastasis:ch1`, data = pheno_GSE147493)

pdf("./plots/06Boxplot_GSE145493_m6Ascore.reg.ov.pdf", width = 4, height = 4)
boxplot(m6Ascore ~ `metastasis:ch1`, data = pheno_GSE147493,
        col = pal.prog[7:8],
        #names = c("", "A", "", "", "B", ""), 
        xaxs = FALSE, las=1,
        ylab = "m6A score",
        ylim = c(-6,6),
        xlab = "p = 0.02512",
        cex.lab=1.2, 
        cex.axis=1.2,
        main= "m6A Score (m6A.reg.ov) in GSE147493")
dev.off()

#ROC of m6A score 
pheno_GSE147493
library(pROC)
grep("SNPH", rownames(eSet_GSE147493))
risk <- data.frame(score= pheno_GSE147493$m6Ascore, stringsAsFactors = F)
risk$SNPH <- t(eSet_GSE147493["SNPH",])
risk$recurrence <- ifelse(pheno_GSE147493$`metastasis:ch1` == "Metastasis", 1, 0)
risk$recurrence <- as.numeric(risk$recurrence)


pdf("./plots/06ROC_m6Ascore_metas_GSE147493.pdf", width = 4, height = 4)
modelroc <- roc(risk$recurrence, risk$score)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     xlim = c(1,0),
     main = "ROC for prediction of recurrence by m6A score")

modelroc <- roc(risk$recurrence, risk$SNPH)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     xlim = c(1,0),
     main = "ROC for prediction of recurrence by SNPH")
dev.off()



######GSE116918####
library(GEOquery)
GSEset <- getGEO("GSE116918", destdir = "./GSE116918_Jain/", getGPL = F)
show(GSEset)
eSet_GSE116918 <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pheno_GSE116918 <- pData(phenoData(GSEset[[1]]))[,c(2,1, 36:44)]
pheno_GSE116918
pheno_GSE116918$`met event (1=yes, 0=no):ch1` <- ifelse(pheno_GSE116918$`met event (1=yes, 0=no):ch1` == 1, "Metastasis", "Non metastasis")

rownames(eSet_GSE116918)[5000:5300]
gpl <- getGEO('GPL25318', destdir="./GSE116918_Jain/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(1,7)]) ## you need to check this , which column do you need
probe2symbol=Table(gpl)[,c(1,7)]
probe2symbol <- probe2symbol[match(rownames(eSet_GSE116918), probe2symbol$ID),]
rownames(eSet_GSE116918) <-  probe2symbol$`Gene Symbol`
eSet_GSE116918 <- aggregate(x = eSet_GSE116918,by = list(rownames(eSet_GSE116918)), FUN = median)
eSet_GSE116918 <- eSet_GSE116918[!grepl("---", eSet_GSE116918$Group.1),]#remove ---
eSet_GSE116918 <- eSet_GSE116918[!grepl("///", eSet_GSE116918$Group.1),]#remove ///
rownames(eSet_GSE116918) <- eSet_GSE116918$Group.1
eSet_GSE116918 <- eSet_GSE116918[,-1]

eSet_GSE116918.m6A <- eSet_GSE116918[which(rownames(eSet_GSE116918) %in% all.regulators),]
eSet_GSE116918.m6A <- aggregate(x = eSet_GSE116918.m6A,by = list(rownames(eSet_GSE116918.m6A)), FUN = median)
rownames(eSet_GSE116918.m6A) <- eSet_GSE116918.m6A$Group.1
eSet_GSE116918.m6A <- eSet_GSE116918.m6A[,-1]

ttest.df <- as.data.frame(t(eSet_GSE116918.m6A) )
ttest.df$type <- pheno_GSE116918$`met event (1=yes, 0=no):ch1`
#select the overlap m6A regulators
pval <- c()
for (i in 1:(ncol(ttest.df)-1)) {
  tmp <- round(t.test(ttest.df[,i]~ ttest.df$type, na.rm =T)$p.value, 4)
  pval <- c(pval, tmp)
}
ttest.df <- data.frame(gene = colnames(ttest.df)[1:(ncol(ttest.df)-1)], pval = pval)
write.csv(ttest.df, file = "./ttest_GSE116918.csv")

###m6A score 
pca.out = prcomp(t(eSet_GSE116918.m6A[m6A.reg.ov,]), scale=TRUE)
pca.out$x
pca.score <- pca.out$x[,1]
pca.score 
pheno_GSE116918$m6Ascore <- pca.score
t.test(m6Ascore~ `met event (1=yes, 0=no):ch1`, data = pheno_GSE116918)

pdf("./plots/06Boxplot_GSE116918_m6Ascore.reg.ov.pdf", width = 4, height = 4)
boxplot(m6Ascore ~ `met event (1=yes, 0=no):ch1`, data = pheno_GSE116918,
        col = pal.prog[7:8],
        #names = c("", "A", "", "", "B", ""), 
        xaxs = FALSE, las=1,
        ylab = "m6A score",
        ylim = c(-6,6),
        xlab = "p = 0.005646",
        cex.lab=0.01, cex.axis=1,
        main= "m6A Score (m6A.reg.ov) in GSE147493")
dev.off()

#ROC of m6A score 
pheno_GSE116918
library(pROC)
grep("SNPH", rownames(eSet_GSE116918))
risk <- data.frame(score= pheno_GSE116918$m6Ascore, stringsAsFactors = F)
eSet_GSE116918 <- as.data.frame(eSet_GSE116918)
risk$SNPH <- eSet_GSE116918["SNPH",]
risk$recurrence <- ifelse(pheno_GSE116918$`metastasis:ch1` == "Metastasis", 1, 0)
risk$recurrence <- as.numeric(risk$recurrence)
modelroc <- roc(risk$recurrence, risk$score)

pdf("./plots/06ROC_m6Ascore_metas_GSE116918.pdf", width = 4, height = 4)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     xlim = c(1,0),
     main = "ROC for prediction of recurrence")

modelroc <- roc(risk$recurrence, risk$SNPH)
plot(modelroc, print.auc=TRUE, auc.polygon=F, grid=c(0.1, 0.2),
     #grid.col=c("green", "red"),
     max.auc.polygon=F,
     auc.polygon.col="skyblue", 
     print.thres=F,
     xlim = c(1,0),
     main = "ROC for prediction of recurrence")
dev.off()



######GSE46691####
GSEset <- getGEO("GSE46691", destdir = "./GSE46691_Erho/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pheno_GSE46691 <- pData(phenoData(GSEset[[1]]))[,c(2,1, 36,37)]
pheno_GSE46691

library(ff)
library(data.table)
eSet_GSE46691 <- fread("./GSE46691_Erho/GSE46691_quantile_normalized.txt.gz")
eSet_GSE46691 <- as.data.frame(eSet_GSE46691)
rownames(eSet_GSE46691) <- eSet_GSE46691$ID_REF
eSet_GSE46691 <- eSet_GSE46691[,-1]
colnames(eSet_GSE46691) <- pheno_GSE46691$geo_accession


gpl <- getGEO('GPL5188', destdir="./GSE46691_Erho/")
colnames(Table(gpl)) ## [1] 49395 24
head(Table(gpl)[,c(10,13)]) ## you need to check this , which column do you need
probe2symbol=Table(gpl)[,c(10,13)]
probes <- rownames(eSet_GSE46691)

library("biomaRt")
ensembl <- useMart("ensembl")
head(listDatasets(ensembl))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl))
probe2symbol <- getBM(attributes = c("affy_huex_1_0_st_v2", "hgnc_symbol"), 
      filters = "affy_huex_1_0_st_v2", 
      values = probes, 
      mart = ensembl)



