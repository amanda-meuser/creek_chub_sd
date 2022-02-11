# Load packages
library(RColorBrewer)
library(VennDiagram)

# Create colour palette
colors <- brewer.pal(12, "Paired")

# Import text file containing columns of significant markers, per analysis, per data set
fishy <- read.delim("./FST_DAPC_outliers.txt", header = TRUE, sep = "\t")

#----------------------------------------------------------------
# 3052-marker data set
#----------------------------------------------------------------

# Calculate the number of entries for FST, DAPC, and overlap
(venn3052 <- calculate.overlap(fishy[5:6]))

# Remove empty list entries
(venn3052 <- lapply(venn3052, function(x) x[!x %in% ""]))

#----------------------------------------------------------------
# 16020-marker data set
#----------------------------------------------------------------

# Calculate the number of entries for FST, DAPC, and overlap
(venn16020 <- calculate.overlap(fishy[3:4]))

# Remove empty list entries
(venn16020 <- lapply(venn16020, function(x) x[!x %in% ""]))

#----------------------------------------------------------------
# 1132-marker data set
#----------------------------------------------------------------

# Calculate the number of entries for FST, DAPC, and overlap
(venn1132 <- calculate.overlap(fishy[1:2]))

# Remove empty list entries
venn1132 <- lapply(venn1132, function(x) x[!x %in% ""])

#----------------------------------------------------------------
# Create Venn Diagrams
#----------------------------------------------------------------

# 3052
plot3052 <- draw.pairwise.venn(area1 = length(venn3052[[1]]), area2 = length(venn3052[[2]]), cross.area = length(venn3052[[3]]), category = c("FST Outliers", "DAPC Outliers"), lty = rep("blank", 2), fill = colors[1:2], alpha = rep(0.75, 2), ext.text = F, cat.pos = c(180,180), cat.dist = rep(0.025, 2), margin = 0.01)


# 16020
plot16020 <- draw.pairwise.venn(area1 = length(venn16020[[1]]), area2 = length(venn16020[[2]]), cross.area = length(venn16020[[3]]), category = c("FST Outliers", "DAPC Outliers"), lty = rep("blank", 2), fill = colors[7:8], alpha = rep(0.75, 2), cat.pos = c(180,180), cat.dist = rep(0.025, 2), margin = 0.01)


# 1132
plot1132 <- draw.pairwise.venn(area1 = length(venn1132[[1]]), area2 = length(venn1132[[2]]), cross.area = length(venn1132[[3]]), category = c("FST Outliers", "DAPC Outliers"), lty = rep("blank", 2), fill = colors[3:4], alpha = rep(0.75, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2), margin = 0.01)

# Create PDF with plots
pdf("venndiagrams_histograms.pdf", width=5, height=12)

pushViewport(plotViewport(layout=grid.layout(3, 1)))
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=1))
grid.draw(plot3052)
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=2))
grid.draw(plot16020)
popViewport()
pushViewport(plotViewport(layout.pos.col=1, layout.pos.row=3))
grid.draw(plot1132)

dev.off()
