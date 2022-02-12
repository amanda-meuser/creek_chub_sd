#to create the .csv file of chromosome lengths (so that the script can put them in order), run the following line on your assembly file (.fasta or .fa):
#awk '/^>/ {if (seqlen){print seqlen}; printf $0"\t";seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' FILENAME.fasta > chromosome_lengths.csv


library(qqman)
library(RColorBrewer)
library(readr)

colors <- brewer.pal(12, "Paired")

##to subset to just the top 50 scaffolds 
##use chromosome_lengths csv file created from assembly fasta file
scaff <- read_csv("./chromosome_lengths.csv")
scaffolds <- scaff$scaffold

##import GEMMA output
data <- read_tsv("./data_lmm.assoc_1132.txt")

##start with first 2, then iterate over full data set
 data1 <- subset(data, data$chr == "Sc8LyPk_34;HRSCAF=465")
data2 <- subset(data, data$chr == "Sc8LyPk_143;HRSCAF=732")
 new <- rbind(data1,data2)
 i=3
for(i in 1:25){
	s = as.character(scaffolds[i])
	new_subset <- subset(data, data$chr == s)
	new <- rbind(new,new_subset)

}
##create dataset to become "chromosome" with all remaining bits not in 'new'
xtra = data
xtra = xtra[!(xtra$chr %in% new$chr),]
xtra$chr = gsub(".*_", "", xtra$chr)
xtra$chr = gsub(";.*", "", xtra$chr)
colnames(xtra) <- c("CHR", "SNP", "BP", "n_miss", "allele1", "allele0", "af", "beta", "se", "log1_H1", "l_remle", "P")
xtra$CHR <- as.numeric(xtra$CHR)
##label as chromosome U
xtra$CHR <- NA
xtra$CHR[is.na(xtra$CHR)]<-0

##prep for QQMan
new$chr = gsub(".*_", "", new$chr)
new$chr = gsub(";.*", "", new$chr)
colnames(new) <- c("CHR", "SNP", "BP", "n_miss", "allele1", "allele0", "af", "beta", "se", "log1_H1", "l_remle", "P")
new$CHR <- as.numeric(new$CHR)

##rename 1-25
count = 1	
c = new$CHR
for(i in 1:length(unique(c))){
  this.c <- unique(c)[i]
  new$CHR[new$CHR == this.c] <- count
  count = count + 1 
}

##add 'extra chromosome' back into new
new <- rbind(new,xtra)

manhattan(new, ylim = c(0, 6.5), col = colors[3:4], genomewideline = -log10(5e-6), suggestiveline = F)
