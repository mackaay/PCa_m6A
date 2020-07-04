######GSE46691####
library(GEOquery)
GSEset <- getGEO("GSE46691", destdir = "./GSE46691_Erho/", getGPL = F)
show(GSEset)
pData(phenoData(GSEset[[1]]))[,c(2,1)]
pData(phenoData(GSEset[[1]]))
colnames(pData(phenoData(GSEset[[1]])))
pheno_GSE46691 <- pData(phenoData(GSEset[[1]]))[,c(2,1, 36,37)]

library(ff)
library(data.table)
eSet_GSE46691 <- fread("./GSE46691_Erho/GSE46691_quantile_normalized.txt.gz")
eSet_GSE46691 <- as.data.frame(eSet_GSE46691)
rownames(eSet_GSE46691) <- eSet_GSE46691$ID_REF
eSet_GSE46691 <- eSet_GSE46691[,-1]
colnames(eSet_GSE46691) <- pheno_GSE46691$geo_accession

probes <- rownames(eSet_GSE46691)

library("biomaRt")
ensembl <- useMart("ensembl")
head(listDatasets(ensembl))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
print("-------------start getBM-------------------")
probe2symbol <- getBM(attributes = c("affy_huex_1_0_st_v2", "hgnc_symbol"), 
                      filters = "affy_huex_1_0_st_v2", 
                      values = probes[1:700000], 
                      mart = ensembl)
save(probe2symbol, "./probe2symbol.RData")