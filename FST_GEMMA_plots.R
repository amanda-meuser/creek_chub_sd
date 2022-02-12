
#load color palette
library(RColorBrewer)
colors <- brewer.pal(12, "Paired")

#---------------------------------------------------------------------------------
#load and assess data, with the 16020, 3052, and 1132-marker data sets, respectively
#---------------------------------------------------------------------------------

fst1 <- (read_tsv("./sd_3052.weir.fst"))
fst2 <- (read_tsv("./sd_16020.weir.fst"))
fst3 <- (read_tsv("./sd_1132.weir.fst"))

gemma1 <- (read_tsv("./data_lmm_3052.assoc.txt"))
gemma2 <- (read_tsv("./data_lmm_16020.assoc.txt"))
gemma3 <- (read_tsv("./data_lmm_1132.assoc.txt"))

#---------------------------------------------------------------------------------
#merge the Fst and GWAS dataframes then plot them, to look for association
#---------------------------------------------------------------------------------

pdf("fst_gemma.pdf", width=5, height=9)
par(mfrow=c(3,1))

merged1 <- merge(fst1, gemma1, by.x="CHROM", by.y="chr", all=F)
plot(merged1$WEIR_AND_COCKERHAM_FST, -log10(merged1$p_wald), 
     col = colors[2], 
     ylim = c(0, 4), 
     xlim = c(-0.15, 0.3), 
     xlab = "Fst Value", 
     ylab = "-log10(p-value)", 
     main = "3052 marker data set")

merged2 <- merge(fst2, gemma2, by.x="CHROM", by.y="chr", all=F)
plot(merged2$WEIR_AND_COCKERHAM_FST, -log10(merged2$p_wald), 
     col = colors[8], 
     ylim = c(0, 4), 
     xlim = c(-0.15, 0.3), 
     xlab = "Fst Value", 
     ylab = "-log10(p-value)", 
     main = "16020 marker data set")


merged3 <- merge(fst3, gemma3, by.x="CHROM", by.y="chr", all=F)
plot(merged3$WEIR_AND_COCKERHAM_FST, -log10(merged3$p_wald), 
     col = colors[4], 
     ylim = c(0, 4), 
     xlim = c(-0.15, 0.3), 
     xlab = "Fst Value", 
     ylab = "-log10(p-value)", 
     main = "1132 marker data set")

dev.off()

#---------------------------------------------------------------------------------
#Perform SLR on GWAS/Fst association
#---------------------------------------------------------------------------------

p_value1 <- (-log10(merged1$p_wald))
gwas.fst.1.lm <- lm(p_value1 ~ merged1$WEIR_AND_COCKERHAM_FST, data = merged1)
summary(gwas.fst.1.lm)

p_value2 <- (-log10(merged2$p_wald))
gwas.fst.2.lm <- lm(p_value2 ~ merged2$WEIR_AND_COCKERHAM_FST, data = merged2)
summary(gwas.fst.2.lm)

p_value3 <- (-log10(merged3$p_wald))
gwas.fst.3.lm <- lm(p_value3 ~ merged3$WEIR_AND_COCKERHAM_FST, data = merged3)
summary(gwas.fst.3.lm)
