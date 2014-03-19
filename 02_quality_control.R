# Task 2: Quality control
# created by Omar AlOmeir
# date: March 18, 2014

library(RColorBrewer)
library(lattice)
library(preprocessCore)

dat <- read.table("data/GSE1710-data.tsv")
des <- readRDS("data/GSE1710-design.rds")

dat <- na.exclude(dat)
mat <- as.matrix(dat)

# heatmap of the forst 100 rows
bluesFun <- colorRampPalette(brewer.pal(n = 9, "Blues"))
heatmap(mat[1:100,], Rowv = NA, Colv = NA, col = bluesFun(256))

# smaples correlation heatmap
sColors <- cm.colors(nrow(cor(mat)))
heatmap(cor(mat), Rowv = NA, Colv = NA, scale = "none", col = bluesFun(256), ColSideColors = sColors, RowSideColors = sColors, margins = c(10, 10))

# boxplot of means of correlation to check for outliers
boxplot(apply(cor(mat), 2, mean))
which.min(apply(cor(mat), 2, mean))

# GSM29606 has the lowest correlation of all, and can be removed as an outlier
outlier <- subset(des, samples == "GSM29606")
olGroup <- subset(des, group == outlier$group & sex == outlier$sex)
splom(dat[olGroup$samples])

# removing the outlier
roDes <- subset(des, samples != "GSM29606")
roDat <- dat[roDes$samples]
heatmap(cor(as.matrix(roDat)), Rowv = NA, Colv = NA, scale = "none", col = bluesFun(256))
boxplot(apply(cor(as.matrix(roDat)), 2, mean))

# normalizing the data
nMat <- normalize.quantiles(as.matrix(roDat))
colnames(nMat) <- colnames(as.matrix(roDat))
rownames(nMat) <- rownames(as.matrix(roDat))
nDat <- as.data.frame(nMat)

# saving the normalized data without the outlier
write.table(nDat, "data/GSE1710-normalized-data.tsv")
saveRDS(roDes, "data/GSE1710-outlier-removed-design.rds")
