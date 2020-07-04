writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

ttest.df <- read.csv("./ttest_NT.csv", header = T, stringsAsFactors = F)
m6A.reg.NTsig <- ttest.df[ttest.df$pval< 0.05,]$gene

meta.regulators <- c("HNRNPA2B1", "FMR1", "METTL14", "KIAA1429","YTHDF1", 
                     "ALKBH5", "HNRNPC")
m6A.reg.ov <- intersect(meta.regulators, m6A.reg.NTsig)

library(maftools)
load(file = "./mutation/maf.PRAD.RData")

plotmafSummary(maf = plotmafSummary(maf = maf.PRAD, rmOutlier = TRUE, 
                                    addStat = 'median', 
                                    dashboard = TRUE, titvRaw = FALSE), 
               rmOutlier = TRUE, addStat = 'median', 
               dashboard = TRUE, titvRaw = FALSE)
#oncoplot for top ten mutated genes.
oncoplot(maf = maf.PRAD, top = 10)
pdf("./plots/08Oncoplot_regulator.pdf", width = 6, height = 6)
oncoplot(maf = maf.PRAD, genes = all.regulators)
dev.off()
