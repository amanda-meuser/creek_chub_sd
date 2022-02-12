
#load packages
library(RColorBrewer)
library(EnvStats)

#create colorblind-friendly palette
colors <- brewer.pal(12, "Paired")

#---------------------------------------------------------------------------------
#load and assess data, with the 16020, 3052, and 1132-marker data sets, respectively
#---------------------------------------------------------------------------------

#pull out Fst values from .fst files, after importing to R
file1 <- (read_tsv("./fst_males_females.weir.fst"))
fst1 <- file1$WEIR_AND_COCKERHAM_FST
file2 <- (read_tsv("./fst_males_females2.weir.fst"))
fst2 <- file2$WEIR_AND_COCKERHAM_FST
file3 <- (read_tsv("./fst_aligned.weir.fst"))
fst3 <- file3$WEIR_AND_COCKERHAM_FST

#remove NaNs
fst1 <- na.omit(fst1)
fst2 <- na.omit(fst2)
fst3 <- na.omit(fst3)

#find the 90th and 99th quantiles
quantile(fst1, c(0.9, 0.99))
quantile(fst2, c(0.9, 0.99))
quantile(fst3, c(0.9, 0.99))

#---------------------------------------------------------------------------------
#make a histogram of the values
#---------------------------------------------------------------------------------

pdf("fst_hist.pdf", width=5, height=9)
par(mfrow=c(3,1))

hist(fst1, 
     xlab = "Fst Value", 
     col = colors[1:2], 
     xlim = c(-0.2, 0.3), 
     ylim = c(0, 1200),
     breaks = 20, 
     main = "Histogram of Fst for 3052 marker data set")


hist(fst2, 
     xlab = "Fst Value", 
     col = colors[7:8], 
     xlim = c(-0.2, 0.4), 
     ylim = c(0, 6000),
     breaks = 20, 
     main = "Histogram of Fst for 16020 marker data set")


hist(fst3, 
     main = "",
     xlab = "Fst Value", 
     col = colors[3:4], 
     xlim = c(-0.2, 0.4), 
     ylim = c(0, 500),
     breaks = 20)

dev.off()

#---------------------------------------------------------------------------------
#create a scatter plot, in order
#---------------------------------------------------------------------------------

pdf("scatter_plots.pdf", width=5, height=9)
par(mfrow=c(3,1))

plot(sort(fst1, decreasing = FALSE), 
     abline(h = 0, col = "black", lty=2),
     pch = 20, 
     col = colors[2], 
     ylim = c(-0.15, 0.4), 
     xlab = "Marker Number", 
     ylab = "Fst Value", 
     main = "3052 marker data set")


plot(sort(fst2, decreasing = FALSE), 
     abline(h = 0, col = "black", lty=2),
     pch = 20,
     col = colors[8], 
     ylim = c(-0.15, 0.4), 
     xlab = "Marker Number", 
     ylab = "Fst Value", 
     main = "16020 marker data set")

plot(sort(fst3, decreasing = FALSE), 
     abline(h = 0, col = "black", lty=1),
     abline(h = 0.1336070, col = "black", lwd = 1, lty = 2),
     pch = 20,
     col = colors[4], 
     ylim = c(-0.15, 0.4), 
     xlab = "Marker Number", 
     ylab = "Fst Value", 
     main = "1132 marker data set")

dev.off()

#---------------------------------------------------------------------------------
#Assessing significance of the higest values as outliers 
#---------------------------------------------------------------------------------

rosnerTest(fst1, k=100)
rosnerTest(fst2, k=100)
rosnerTest(fst3, k=100)
