rm(list = ls())

deps <- c("gelnet","dplyr","gdata","DT")
for(pkg in deps)  if (!pkg %in% installed.packages()) install.packages(pkg, dependencies = TRUE)

setwd("./mRNAsi/")
library(gelnet)
library(dplyr)
library(biomaRt)

load(file ="./eSet_regulator_FPKM.RData")
FPKM <- read.csv("./eSet.csv", header = T, stringsAsFactors = F)
FPKM <- FPKM[!duplicated(FPKM$id),]
rownames(FPKM) <- FPKM$id
FPKM <- FPKM[,2:551]
Tumor.idx <- which(!grepl("11", colnames(FPKM)))
FPKM <- FPKM[,Tumor.idx]
eSet <- log2(FPKM+1)
write.table(eSet, file = "./mRNAsi/eSet.txt", sep = "\t", row.names = T)

###mRNAsi

##training
genes2hugo <- function( v, srcType = "ensembl_gene_id" )
{
  ## Retrieve the EMSEMBL -> HUGO mapping
  ensembl <- biomaRt::useMart( "ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl" )
  ID <- biomaRt::getBM( attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl )
  
  ## Make sure there was at least one mapping
  if( nrow(ID) < 1 ) stop( "No IDs mapped successfully" )
  
  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  stopifnot( all( ID[,1] %in% v ) )
  
  ID
}

main.train <- function( fnOut = "pcbc-stemsig.tsv", fnGenes = NULL )
{
  ## Load RNAseq data
  #  synRNA <- synGet( "syn2701943", downloadLocation = "/data/PCBC" )
  #  X <- read.delim( synRNA@filePath ) %>%
  X <- read.delim( "./rnaseq_norm.tsv" ) %>%
    tibble::column_to_rownames( "tracking_id" ) %>%
    as.matrix()
  
  ## Retrieve metadata
  #synMeta <- synTableQuery( "SELECT UID, Diffname_short FROM syn3156503" )
  Y <- read.csv("./PCBC_metadata.csv") %>%
    dplyr::select(UID, Diffname_short) %>%
    mutate( UID = gsub("-", ".", UID) ) %>%
    tibble::column_to_rownames( "UID" )
  
  ## Retrieve the labels from the metadata
  y <- Y[colnames(X),]
  names(y) <- colnames(X)
  
  ## Fix the missing labels by hand
  y["SC11.014BEB.133.5.6.11"] <- "EB"
  y["SC12.039ECTO.420.436.92.16"] <- "ECTO"
  
  ## Drop the splice form ID from the gene names
  v <- strsplit( rownames(X), "\\." ) %>% lapply( "[[", 1 ) %>% unlist()
  rownames(X) <- v
  
  ## Map Ensembl IDs to HUGO
  V <- genes2hugo( rownames(X) )
  X <- X[V[,1],]
  rownames(X) <- V[,2]
  
  ## Reduce the gene set to the provided list (if applicable)
  if( is.null( fnGenes ) == FALSE )
  {
    vGenes <- read.delim( fnGenes, header=FALSE ) %>% as.matrix() %>% drop()
    VE <- genes2hugo( vGenes, "entrezgene" )
    X <- X[intersect( rownames(X), VE[,2] ),]
  }
  
  ## Mean-center the data
  m <- apply( X, 1, mean )
  X <- X - m
  
  ## Identify stem cell samples
  j <- which( y == "SC" )
  X.tr <- X[,j]
  X.bk <- X[,-j]
  
  ## Train a one-class model
  mm <- gelnet( t(X.tr), NULL, 0, 1 )
  
  ## Store the signature to a file
  write.table(mm$w, file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
  
  ## Perform leave-one-out cross-validation
  auc <- c()
  for( i in 1:ncol(X.tr) )
  {
    ## Train a model on non-left-out data
    X1 <- X.tr[,-i]
    m1 <- gelnet( t(X1), NULL, 0, 1 )
    
    ## Score the left-out sample against the background
    s.bk <- apply( X.bk, 2, function(z) {cor( m1$w, z, method="sp" )} )
    s1 <- cor( m1$w, X.tr[,i], method="sp" )
    
    ## AUC = P( left-out sample is scored above the background )
    auc[i] <- sum( s1 > s.bk ) / length(s.bk)
    cat( "Current AUC: ", auc[i], "\n" )
    cat( "Average AUC: ", mean(auc), "\n" )
  }
  
  return(auc)
}

main.train()



###prediction
library(dplyr)
#mRNAsi
main.predict <- function( fnSig = "pcbc-stemsig.tsv", fnOut = "mRNA_StemScore.tsv" )
{
  ## Load the signature
  w <- read.delim( fnSig, header=FALSE, row.names=1 ) %>% as.matrix() %>% drop()
  
  ## Reduces HUGO|POSITION gene IDs to just HUGO
  f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )
  
  s <- "eSet.txt"
  X <- read.delim( s, as.is=TRUE, check.names=FALSE ) 
  X <- X[rownames(X) %in% names(w), ]
  #%>%  ## Read the raw values
    #filter( !grepl( "\\?", gene_symbol ) ) %>%      ## Drop genes with no mapping to HUGO
    #mutate( gene_symbol = f( gene_symbol ) ) %>%        ## Clip gene ids to HUGO
    #filter( gene_symbol %in% names(w) )         ## Reduce to the signature's gene set
  
  ## SLC35E2 has multiple entries with the same HUGO id
  ## Keep the first entry only
  j <- grep( "SLC35E2", rownames(X) )
  if( length(j) > 1 )
    X <- X[-j[-1],]
  
  ## Convert to a matrix
  X$symbol <- rownames(X) 
  rownames(X) <- NULL
  X <- X %>% tibble::column_to_rownames( "symbol" ) %>% as.matrix()
  
  ## Reduce the signature to the common set of genes
  stopifnot( all( rownames(X) %in% names(w) ) )
  w <- w[ rownames(X) ]
  
  ####### Score via Spearman correlation
  s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )
  
  ## Scale the scores to be between 0 and 1
  s <- s - min(s)
  s <- s / max(s)
  
  write.table(cbind(s), file = fnOut, sep = "\t", quote = FALSE, col.names = FALSE)
}

main.predict("pcbc-stemsig.tsv", "mRNA_StemScore.tsv")



####add mRNAsi into sampleinfo file
setwd("../")
library(RColorBrewer)
pal.grade <- rev(brewer.pal(6, "Accent"))
pal.subtype <- brewer.pal(8, "Dark2")
load(file ="./eSet_regulator_FPKM.RData")

StemScore <- read.table("./mRNAsi/mRNA_StemScore.tsv")
sampleinfo$StemIndex <- StemScore$V2
save(eSet, sampleinfo, 
     file = "./eSet_regulator_FPKM.RData")

pdf("./plots/05Boxplot_subtype_StemIndex.pdf", width = 5, height = 5)
boxplot(sampleinfo$StemIndex ~ sampleinfo$subtype, col = pal.subtype, 
        ylab = "Stem Index", 
        xlab = "ANOVA, p = 0.00332", 
        main = "Stem Index between Subtypes", 
        ylim = c(0, 1.2))
summary(aov(sampleinfo$StemIndex ~ sampleinfo$subtype))
TukeyHSD(aov(sampleinfo$StemIndex ~ sampleinfo$subtype))
dev.off()

violin.df <- data.frame(Grade = sampleinfo$grade, 
                        Subtype = sampleinfo$subtype, 
                        StemIndex = sampleinfo$mRNAsi, 
                        Gender = sampleinfo$gender,
                        Age_number = sampleinfo$age,
                        stringsAsFactors = F)

library(ggplot2)
pdf("./plots/05Violinplot_mRNAsi.pdf", width = 5, height = 5)
summary(aov(StemIndex ~ Grade, data = violin.df))
TukeyHSD(aov(StemIndex ~ Grade, data = violin.df))

ggplot(data = violin.df, aes(x = Grade, y = StemIndex))+
  geom_violin(aes(fill = factor(Grade)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Dark2")+
  labs(y = "Stem Index", 
       x = "ANOVA: p = 0.00195 ", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 

summary(aov(StemIndex ~ Subtype, data = violin.df))
TukeyHSD(aov(StemIndex ~ Subtype, data = violin.df))

ggplot(data = violin.df, aes(x = Subtype, y = StemIndex))+
  geom_violin(aes(fill = factor(Subtype)))+
  geom_boxplot(width = 0.3, outlier.size = -1)+
  geom_jitter(height = 0.002, width = 0.2, cex = 1.2)+
  scale_fill_brewer(palette="Set1")+
  labs(y = "Stem Index", 
       x = "ANOVA: p = 6.79e-13", 
       title = "")+ 
  #ylim(0,1)+
  #theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none") 
dev.off()
