# Task 2: Quality control
# created by Omar AlOmeir
# date: March 18, 2014

library(RColorBrewer)
library(lattice)
library(preprocessCore)
library(impute)

# function to extract genes
prepareData <- function(g) {
  pDat0 <- data.frame()
  for (i in 1:length(g)) {
    pDat <- data.frame(des, gExp = as.vector(t(as.matrix(dat[g[i], ]))), 
                       gene = g[i])
    pDat0 <- rbind(pDat0, pDat)
  }
  pDat0
}

dat <- read.table("data/GSE1710-data.tsv")
des <- readRDS("data/GSE1710-design.rds")
str(dat, max.level=0)
str(des)

# Plot two random genes
stripplot(gExp ~ group | gene, prepareData(sample(rownames(dat),2)), group = sex, jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), asp =1, grid = TRUE)

# heatmap of the first 100 rows before dealing with missing values
mat <- as.matrix(dat)
bluesFun <- colorRampPalette(brewer.pal(n = 9, "Blues"))
heatmap(mat[1:100,], Rowv = NA, Colv = NA, col = bluesFun(256))

# dealing with missing values

# remove rows that contain more than 8 (25%) missing values
dat <- dat[rowSums(is.na(dat))<=8,]
str(dat, max.level=0)

# use knn to impute other NAs
mat <- as.matrix(dat)
mat <- impute.knn(mat)
dat <- as.data.frame(mat$data)
str(dat, max.level=0)

# heatmap of the first 100 rows
mat <- as.matrix(dat)
bluesFun <- colorRampPalette(brewer.pal(n = 9, "Blues"))
heatmap(mat[1:100,], Rowv = NA, Colv = NA, col = bluesFun(256))

# density plot of groups
gDat <- data.frame(des, gExp = as.vector(t(as.matrix(dat[c(1:100), ]))))
densityplot(~ gExp, gDat, groups = group, auto.key = TRUE, n=400, type = "p", plot.points = FALSE, grid = TRUE)

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
which.min(apply(cor(as.matrix(roDat)), 2, mean))

# a boxplot of data before normalization
boxplot(roDat, pch = ".", ylab = "Expression level", las = 3)

# normalizing the data
nMat <- normalize.quantiles(as.matrix(roDat))
colnames(nMat) <- colnames(as.matrix(roDat))
rownames(nMat) <- rownames(as.matrix(roDat))
nDat <- as.data.frame(nMat)

# a boxplot of data before normalization
boxplot(nDat, pch = ".", ylab = "Expression level", las = 3)

# saving the normalized data without the outlier
write.table(nDat, "data/GSE1710-normalized-data.tsv")
saveRDS(roDes, "data/GSE1710-outlier-removed-design.rds")