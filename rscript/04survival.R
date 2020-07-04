library(RColorBrewer)
pal.grade <- rev(brewer.pal(6, "Accent"))
pal.subtype <- brewer.pal(8, "Dark2")
library(survival)
load(file ="./eSet_regulator_FPKM.RData")

sampleinfo$OS.status <- sampleinfo$vital_status 
sampleinfo$OS.months <- sampleinfo$days_to_last_followup

my.surv <- Surv(sampleinfo$OS.months,sampleinfo$OS.status== "dead")
kmfit1 <- survfit(my.surv~subtype,data = sampleinfo) 
survdiff(my.surv~subtype,data = sampleinfo)

pdf("./plots/05KMplot_OS_FPKM.pdf", width = 4, height = 3)
plot(kmfit1,
     col = pal.subtype, 
     main = "Overall Survival (TCGA)", 
     xlab = "Days", 
     ylab = "Cumulative Survival Rate")
legend("bottomleft", legend = c("Subtype1", "Subtype2", "Subtype3"), 
       col = pal.subtype, lty = c(1,1), cex = 0.75)
dev.off()



###RFS
my.surv <- Surv(sampleinfo$days_to_first_biochemical_recurrence,sampleinfo$biochemical_recurrence== "yes")
kmfit1 <- survfit(my.surv~subtype,data = sampleinfo) 
survdiff(my.surv~subtype,data = sampleinfo)
coxph(my.surv~subtype,data = sampleinfo)

pdf("./plots/05KMplot_RFS_FPKM.pdf", width = 4, height = 3)
plot(kmfit1,
     col = pal.subtype, 
     main = "Recurrent Free Survival (TCGA)", 
     xlab = "Days", 
     ylab = "Cumulative Survival Rate")
legend("bottomleft", legend = c("Subtype1", "Subtype2", "Subtype3"), 
       col = pal.subtype, lty = c(1,1), cex = 0.5)
dev.off()


sampleinfo$psa.recur <- ifelse(sampleinfo$psa_value <10, "no", "yes")

my.surv <- Surv(sampleinfo$days_to_psa,sampleinfo$psa.recur== "yes")
kmfit1 <- survfit(my.surv~subtype,data = sampleinfo) 

pdf("./plots/05KMplot_PSAsurvival_FPKM.pdf", width = 4, height = 3)
plot(kmfit1,pal.subtype, 
     main = "Recurrent Free Survival (TCGA)", 
     xlab = "Days", 
     ylab = "Cumulative Survival Rate")
legend("bottomleft", legend = c("Subtype1", "Subtype2", "Subtype3"), 
       col = pal.subtype, lty = c(1,1), cex = 0.75)
dev.off()




##table 
table(sampleinfo[sampleinfo$subtype=="Subtype3",]$gleason_level)


#


my.surv <- Surv(sampleinfo$days_to_first_biochemical_recurrence,sampleinfo$biochemical_recurrence== "yes")
coxph(my.surv~m6Ascore,data = sampleinfo)
