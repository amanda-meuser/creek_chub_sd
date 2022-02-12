#this script produces both plots which have scaffold names included on the x-axis. these were removed for clarity in the journal submission, using PowerPoint, as the text was illegible and order of the scaffolds is irrelevant for the assembly used by these two data sets. 

library(RColorBrewer)
library(readr)

#create colour palette
colors <- brewer.pal(12, "Paired")

#import files
data <- (read_tsv("./data_lmm_3052.assoc.txt"))
data2 <- (read_tsv("./data_lmm_16020.assoc.txt"))

pdf("manhattan_plots.pdf", width=10, height=5)
par(mfrow=c(2,1))

data$chr = gsub(";.*", "", data$chr)
ann <- rep(1, length(data$p_wald))
ann[with(data, chr=="Sc1KuOr_12343")]<-2
ann<-factor(ann,levels= 1:2, labels=c(""," "))
manhattan.plot(data$chr, data$ps, data$p_wald, sig.level = 5e-6, annotate = ann, col = colors[2])


data2$chr = gsub(";.*", "", data2$chr)
ann2 <- rep(1, length(data2$p_wald))
ann2[with(data, chr=="Sc1KuOr_12343")]<-2
ann2<-factor(ann2,levels= 1:2, labels=c(""," "))
manhattan.plot(data2$chr, data2$ps, data2$p_wald, sig.level = 5e-6, annotate= ann2, col = colors[8])

dev.off()
