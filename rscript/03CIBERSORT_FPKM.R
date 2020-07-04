load(file ="./eSet_regulator_FPKM.RData")
cibersort.all <- read.csv("./CIBERSORT.Output_PRADall.csv", header = T,
                          stringsAsFactors = F)
rownames(cibersort.all) <- gsub("-", ".", cibersort.all$Input.Sample)
cibersort.all <- cibersort.all[match(rownames(sampleinfo),rownames(cibersort.all)),]
colnames(cibersort.all)
cibersort.all <- cibersort.all[,2:23]
cibersort.all$subtype <- sampleinfo$subtype

#t test
pval <- c()
for (i in 1:22) {
  tmp <- aov(cibersort.all[,i]~subtype, data = cibersort.all)
  tmp <- round(summary(tmp)[[1]][["Pr(>F)"]][1], 3)
  pval <- c(pval, tmp)
}
aov.df <- data.frame(gene = colnames(cibersort.all)[1:22], pval = pval, 
                     stringsAsFactors = F)


#violin plot
tmp <- c()
for (i in 1:22) {
  tmp1 <- cibersort.all[,i]
  tmp <- c(tmp, tmp1)
}
temp <- c()
for (i in 1:22) {
  temp1 <- rep(colnames(cibersort.all)[i],nrow(cibersort.all))
  temp <- c(temp, temp1)
}
subtype <- rep(cibersort.all$subtype, 22)

violin.df <- data.frame(Proportion =tmp, 
                        CellType = temp, 
                        Subtype = subtype, 
                        stringsAsFactors = F)

library(RColorBrewer)
pal.grade <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Dark2")

pdf("./plots/03Boxplot_subtype_immune_FPKM.pdf", width = 10, height = 5)
boxplot(Proportion ~ Subtype+CellType, data = violin.df,
        #at = c(1:3, 5:7), 
        col = pal.subtype[1:3],
        #names = c("", "A", "", "", "B", ""), 
        xaxs = FALSE, las=2,
        ylim = c(-0.02,0.6),
        cex.lab=0.01, cex.axis=0.5,
        main= "m6A Regulators Expression in TCGA")
legend("topleft", fill = pal.subtype[1:3], 
       legend = c("Subtype1", "Subtype2", "Subtype3"), cex = 0.8,
       horiz = F)
dev.off()


######PD-L1 ######
load( file = "./eSet_TumorALL.RData")
PDL1 <- eSet["CD247",]
PDL1 <- as.data.frame(t(PDL1))
PDL1$subtype <- sampleinfo$subtype
summary(aov(CD247~subtype, data = PDL1))

