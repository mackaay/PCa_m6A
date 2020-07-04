library(ff)
library(data.table)
#read DNA methylation data
bVal.PRAD <- fread("./DNAmethy_PRAD.csv.gz")
bVal.PRAD <- as.data.frame(bVal.PRAD)
rownames(bVal.PRAD) <- bVal.PRAD$V1
bVal.PRAD <- bVal.PRAD[,-1]
colnames(bVal.PRAD) <- gsub("CGA", "TCGA",colnames(bVal.PRAD) )
colnames(bVal.PRAD) <- gsub("TTCGA", "TCGA",colnames(bVal.PRAD) )
colnames(bVal.PRAD) <- substr(colnames(bVal.PRAD), start = 1L, stop = 15L)
#read annotation 450K dataset
ann450.regulator <- read.csv(file = "./ann450k.regulator.csv", header = T, stringsAsFactors = F)
ann450.regulator <- ann450.regulator[order(ann450.regulator$genesymbol),]
#retain m6A score info
meta.regulators <- c("HNRNPA2B1", "FMR1", "METTL14", "KIAA1429","YTHDF1", 
                     "ALKBH5", "HNRNPC")
m6A.reg.ov <- intersect(meta.regulators, m6A.reg.NTsig)
ann450.regulator <- ann450.regulator[ann450.regulator$genesymbol %in% m6A.reg.ov,]

bVal.PRAD.reg <- bVal.PRAD[match(ann450.regulator$X, rownames(bVal.PRAD)),]
colnames(bVal.PRAD.reg) <- gsub("-",".", colnames(bVal.PRAD.reg))

load(file ="./eSet_regulator_FPKM.RData")
bVal.PRAD.reg <- bVal.PRAD.reg[,match(colnames(eSet), colnames(bVal.PRAD.reg) )]
bVal.PRAD.reg <- bVal.PRAD.reg[complete.cases(bVal.PRAD.reg),]

bVal.PRAD.reg <- bVal.PRAD.reg[,order(sampleinfo$subtype)]
sampleinfo <- sampleinfo[order(sampleinfo$subtype),]

library(pheatmap)
library(RColorBrewer)
pal.grade <- brewer.pal(8, "Accent")
pal.subtype <- brewer.pal(8, "Dark2")

annotation_col = data.frame(
  Subtype = sampleinfo$subtype
)
rownames(annotation_col) = colnames(bVal.PRAD.reg)
names(pal.subtype) <- levels(as.factor(sampleinfo$subtype))

pdf("./plots/07Heatmap_DNAme_m6Ascore.pdf", width = 10, height = 10)
pheatmap(bVal.PRAD.reg, 
         scale = "none",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = T, show_colnames = T,
         fontsize_row = 10, fontsize_col = 1,
         cluster_cols = F,
         cluster_rows = F,
         annotation_legend = T,
         color = colorRampPalette(c("darkblue",  "yellow"))(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors  = list(
           #Time = c("white", "firebrick"),
           # subtype = c(subtype1 = "#A6CEE3", subtype2 = "#1F78B4",
           #            subtype3 = "#B2DF8A", subtype4 = "#33A02C"),
           Subtype = pal.subtype[1:3]
         ),
         main = "m6A regulators"
)
dev.off()

##all DNA methylation
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k.promoter <- ann450k[grep("TSS1500|TSS200", ann450k$UCSC_RefGene_Group),]

bVal.PRAD.promoter <- bVal.PRAD[match(rownames(ann450k.promoter), rownames(bVal.PRAD)),]
colnames(bVal.PRAD.promoter) <- gsub("-",".", colnames(bVal.PRAD.promoter))
bVal.PRAD.promoter <- bVal.PRAD.promoter[,match(colnames(eSet), colnames(bVal.PRAD.promoter) )]
mVal.PRAD.promoter <- log2(bVal.PRAD.promoter/(1- bVal.PRAD.promoter))

library(limma)
library(scales)
pdf("./plots/04PCA_DNApromoter_FPKM.pdf", height = 5, width = 5)
plotMDS(mVal.PRAD.promoter, top=1000, gene.selection="common", pch = 19,
        col=alpha(pal.subtype[factor(sampleinfo$subtype)], 0.5))
legend("topright", legend=levels(factor(sampleinfo$subtype)), text.col=pal.subtype,
       bg="white", cex=1)
dev.off()
