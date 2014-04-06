# Task 5: PCA and clustering
# created by Yiming Zhang
# date: April 1, 2014

#load the package 
library(RColorBrewer)
library(cluster)
library(pvclust)
library(xtable)
library(limma)
library(plyr)
library(ggplot2)
library(car)
library(scatterplot3d)

#read the design matrix and normalized data and have a simple check
dat <- read.table("data/GSE1710-normalized-data.tsv")
des <- readRDS("data/GSE1710-outlier-removed-design.rds")

#separate age into three groups
des$age<- as.numeric( levels ( des$age))[ as.integer (des$age)]
young <- which(des$age<=30)
middleage <- which(des$age>30&des$age<50)
old <- which(des$age>=50)
des$age[young] <- "young"
des$age[middleage] <- "middleage"
des$age[old] <- "old"
des$age<- as.factor(des$age)
des$age <- recode(des$age, "",levels = c('young', 'middleage', 'old'))

#check the data
str(dat, max.level = 0)
str(des)



##simple way to do PCA
pcs <- prcomp(dat, center = F, scale = F) 
plot(pcs) 
prinComp <- cbind(des, pcs$rotation) 

#check how the first few PCs relate to covariates
prinComp <- cbind(des, pcs$rotation)
plot(prinComp[ ,c("group", "sex", "age", "PC1", "PC2", "PC3")],
     pch = 19, cex = 0.8) 

#Get the PC's percentage
svar<- (pcs$sdev)^2
(PC_percent<- svar/sum(svar))


#look at the principle components in experimental group
Group <- des$group
ggplot(prinComp, aes(x=PC1, y=PC2,colour=Group))+geom_point()


#Plot the 3-D scatterplot for PC1,PC2,PC4
# create column indicating point color
# code is modified from code at http://www.r-bloggers.com/getting-fancy-with-3-d-scatterplots/
prinComp$pcolor[Group=="NC"] <- "red"
prinComp$pcolor[Group=="CD"] <- "blue"
prinComp$pcolor[Group=="UC"] <- "green"
with(prinComp, {
  s3d <- scatterplot3d(PC2, PC4, PC1,        # x y and z axis
                       color=pcolor, pch=19,        # circle color indicates no. of cylinders
                       type="h", lty.hplot=2,       # lines to the horizontal plane
                       xlab="PC2 - 3.24%",
                       ylab="PC4 - 0.91%",
                       zlab="PC1 - 91.07%")
  s3d.coords <- s3d$xyz.convert(PC2, PC4, PC1)

  # add the legend
  legend("left", inset=.05,      # location and inset
         bty="n", cex=.5,              # suppress legend box, shrink text 50%
         title="Groups",
         c("NC", "UC", "CD"), fill=c("red", "blue", "green"))
})


#then in sex group
Sex <- des$sex
ggplot(prinComp, aes(x=PC1, y=PC2,colour=Sex))+geom_point()
#PC2 found correlated with sex

#At last in age group
Age <- des$age
ggplot(prinComp, aes(x=PC1, y=PC2,colour=Age))+geom_point()
ggplot(prinComp, aes(x=PC1, y=PC3,colour=Age))+geom_point()
ggplot(prinComp, aes(x=PC1, y=PC4,colour=Age))+geom_point()
ggplot(prinComp, aes(x=PC1, y=PC5,colour=Age))+geom_point()
# No correlated PC is found



##Clustering 
#rescale the rows in data to make later things easier

sdat <- t(scale(t(dat)))
str(sdat, max.level = 0, give.attr = FALSE)

round(data.frame(avgBefore = rowMeans(head(dat)),
                 avgAfter = rowMeans(head(sdat)),
                 varBefore = apply(head(dat), 1, var),
                 varAfter = apply(head(sdat), 1, var)), 2)

# there is nothing useful in sample culstering analysis, so I ignore it 

# gene clustering
#Load the top hits data to do gene clustering
hits <- read.table("data/groups-hits.tsv")
topDat <- subset(dat, rownames(dat) %in% rownames(hits))
#Hierarchical clustering
geneC.dis <- dist(topDat, method = "euclidean")
geneC.hc.a <- hclust(geneC.dis, method = "ward")
plot(geneC.hc.a, labels = FALSE, main = "Hierarchical with Ward Linkage", xlab = "")
#This plot is ugly though...

#Partitioning 
set.seed(1234)
k <- 5
kmeans.genes <- kmeans(topDat, centers = k)

# choose which cluster we want(I have checked all five clusters, but only cluster 4 gave me interesting answer)
clusterNum <- 4
plot(kmeans.genes$centers[clusterNum, ], ylim = c(-0.1, 5), type = "n", xlab = "Samples", 
     ylab = "Relative expression")
# Plot the expression of all the genes in the selected cluster in grey.
matlines(y = t(topDat[kmeans.genes$cluster == clusterNum, ]), col = "grey")
# Add the cluster center. 
points(kmeans.genes$centers[clusterNum, ], type = "l")
points(kmeans.genes$centers[clusterNum, ], col = des$group, pch = 20)
# add the legend
legend("topleft", inset=.05,      # location and inset
       bty="n", cex=.5,              # suppress legend box, shrink text 50%
       title="Groups",
       c("NC", "UC", "CD"), fill=c("black","red","green"))


