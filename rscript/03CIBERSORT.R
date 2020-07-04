load(file ="./eSet_regulator.RData")
cibersort.all <- read.csv("./CIBERSORT.Output_PRADall.csv", header = T,
                          stringsAsFactors = F)
rownames(cibersort.all) <- gsub("-", ".", cibersort.all$Input.Sample)
cibersort.all <- cibersort.all[match(sampleinfo$ID,rownames(cibersort.all)),]
colnames(cibersort.all)
cibersort.all <- cibersort.all[,2:23]
cibersort.all$subtype <- sampleinfo$subtype

#t test
pval <- c()
for (i in 1:22) {
  tmp <- round(t.test(cibersort.all[1:235,i], cibersort.all[236:495,i])$p.value, 3)
  pval <- c(pval, tmp)
}
ttest.df <- data.frame(gene = colnames(cibersort.all)[1:22], pval = pval)


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
                        Subtype = subtype)

library(ggplot2)
library(tidyverse)
library(RColorBrewer)
pal.subtype <- brewer.pal(8, "Dark2")
pdf("./plots/03Violinplot_CIBERSORT.pdf", width = 10, height = 7)
violin.df %>%
  ggplot(aes(x= Subtype, y = Proportion, fill = Subtype))+
  geom_boxplot(outlier.size = 0)+
  #geom_jitter(aes(color = response),width = 0.3, alpha = 0.7)+
  xlab("type")+
  scale_fill_manual(values=pal.subtype)+
  scale_color_manual(values=pal.subtype)+
  facet_wrap(~CellType)+
  #ylim(c(0,1))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("")
dev.off()




