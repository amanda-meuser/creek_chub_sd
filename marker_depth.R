#histogram of depth per individual, across all markers, for each vcf file

#load color palette
library(RColorBrewer)
colors <- brewer.pal(12, "Paired")

# read in depth data
depth0.2 <- read.csv("./CCCS_cccsref0.9test_miss0.2_mac3_Q20_knownsexinds_thin90_maf001.idepth", sep="\t")
depth0.4 <- read.csv("./CCCS_cccsref0.9test_miss0.4_mac3_Q30_DP3_ind95_maf001_1snpcontig_sexinds.recode.idepth", sep="\t")
depth_dacealign <- read.csv("./depth_per_ind_aligned_1132.idepth", sep="\t")

#look at mean and median for each data set
mean(depth0.4$MEAN_DEPTH)
median(depth0.4$MEAN_DEPTH)

mean(depth0.2$MEAN_DEPTH)
median(depth0.2$MEAN_DEPTH)

mean(depth_dacealign$MEAN_DEPTH)
median(depth_dacealign$MEAN_DEPTH)

##make the histograms, with lines at mean and median
pdf("depth_histograms.pdf", width=5, height=9)
par(mfrow=c(3,1))

hist(depth0.4$MEAN_DEPTH,
     col = colors[2:1],
     breaks = seq(0,120,8),
     ylim = c(0, 25),
     xlim = c(0, 120),
     ylab = "Frequency" ,
     xlab = "Depth",
     main = "3052-marker data set")

abline(v = 34, col = "black", lwd = 3, lty = 1)
abline(v = 31, col = "black", lwd = 3, lty = 2)


hist(depth0.2$MEAN_DEPTH,
     col = colors[8:7],
     breaks = seq(0,50,3),
     ylim = c(0, 20),
     xlim = c(0, 50),
     ylab = "Frequency" ,
     xlab = "Depth",
     main = "16,020-marker data set")

abline(v = 11, col = "black", lwd = 3, lty = 1)
abline(v = 8, col = "black", lwd = 3, lty = 2)


hist(depth_dacealign$MEAN_DEPTH,
     col = colors[3:4],
     breaks = seq(0,120,8),
     ylim = c(0, 30),
     xlim = c(0, 120),
     ylab = "Frequency" ,
     xlab = "Depth",
     main = "1132-marker data set")

abline(v = 25, col = "black", lwd = 3, lty = 1)
abline(v = 19, col = "black", lwd = 3, lty = 2)

dev.off()
