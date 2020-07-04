library(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf)


dir.file <- "./mutation/gdac.broadinstitute.org_PRAD.Mutation_Packager_Calls.Level_3.2016012800.0.0/"
filename <- list.files(dir.file, pattern = "*maf*")
filename <- filename[!filename == "TCGA-EJ-A7NG-01.maf.txt.16.txt"]
maf.PRAD <- c()
for (i in filename) {
  tmp <- read.maf(paste(dir.file,i, sep = ""))
  maf.PRAD <- merge_mafs(maf = c(maf.PRAD, tmp))
}


save(maf.PRAD, file = "./mutation/maf.PRAD.RData")

