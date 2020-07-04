library(GEOquery)
GSEset <- getGEO("GSE127985", destdir = "./GSE127985/", getGPL = F)
show(GSEset)
eSet <- exprs(GSEset[[1]])
pData(phenoData(GSEset[[1]]))[,c(2,1)]
colnames(pData(phenoData(GSEset[[1]])))
sampleinfo <- pData(phenoData(GSEset[[1]]))[,c(2,1,35,38)]
sampleinfo$title
bVal.GSE127985 <- eSet
sampleinfo.GSE127985 <- sampleinfo
save(bVal.GSE127985, sampleinfo.GSE127985, 
     file = "./GSE127985/GSE127985.RData")
