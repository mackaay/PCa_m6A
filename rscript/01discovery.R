rm(list = ls())
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

writers <- c("METTL3", "METTL14", "RBM15", "RBM15B", "WTAP", 
             "KIAA1429", "CBLL1", "ZC3H13")
erasers <- c("ALKBH5", "FTO") 
readers <- c("YTHDC1", "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", 
             "IGF2BP1", "HNRNPA2B1", "HNRNPC", "FMR1", "LRPPRC", "ELAVL1")
all.regulators <- c(writers, erasers, readers)

#select all regulator ann450k
ann450k.regulator <- ann450k[grep("METTL3|METTL14|RBM15|RBM15B|WTAP|KIAA1429|CBLL1|ZC3H13|ALKBH5|FTO|YTHDC1|YTHDC2|YTHDF1|YTHDF2|YTHDF3|IGF2BP1|HNRNPA2B1|HNRNPC|FMR1|LRPPRC|ELAVL1", ann450k$UCSC_RefGene_Name),]
#213 cg probes
#select TSS1500 
ann450k.regulator <- ann450k.regulator[grep("TSS1500|TSS200", ann450k.regulator$UCSC_RefGene_Group),]
#95 remain
table(ann450k.regulator$UCSC_RefGene_Name)

write.csv(ann450k.regulator, file = "./ann450k.regulator.csv", sep = ",")




####GSE127985 dataset#####
load("./GSE127985/GSE127985.RData")
ann450k.reg <- read.csv("./ann450k.regulator.csv", stringsAsFactors = F)
ann450k.reg <- ann450k.reg[order(ann450k.reg$genesymbol),]
bVal.GSE127985 <- as.data.frame(bVal.GSE127985)
bVal.reg <- bVal.GSE127985[which(rownames(bVal.GSE127985) %in% ann450k.reg$X),]

library(pheatmap)
library(RColorBrewer)
col.outcome <- brewer.pal(8, "Dark2")
annotation_col = data.frame(
  Outcome = sampleinfo.GSE127985$`prognosis:ch1`
)
rownames(annotation_col) = colnames(bVal.reg)
names(col.outcome) <- levels(as.factor(sampleinfo.GSE127985$`prognosis:ch1`))

pheatmap(bVal.reg, 
         scale = "none",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = T,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = T,
         cluster_rows = F,
         annotation_legend = T,
         color = colorRampPalette(c("darkblue",  "yellow"))(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Time = c("white", "firebrick"),
           # subtype = c(subtype1 = "#A6CEE3", subtype2 = "#1F78B4",
           #            subtype3 = "#B2DF8A", subtype4 = "#33A02C"),
           Outcome = col.outcome
         ),
         main = "m6A regulators"
)
