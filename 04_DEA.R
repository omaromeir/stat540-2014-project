# Task 4: Differential Expression Analysis
# created by Omar AlOmeir and Abrar Wafa
# date: March 20, 2014

library(limma)
library(lattice)
library(car)

dat <- read.table("../stat540-project/data/GSE1710-normalized-data.tsv")
des <- readRDS("../stat540-project/data/GSE1710-outlier-removed-design.rds")

# sperate age into groups using Yiming's code 
des$age <- as.numeric(levels(des$age))[as.integer(des$age)]
young <- which(des$age <= 30)
middleage <- which(des$age > 30 & des$age < 50)
old <- which(des$age >= 50)
des$age[young] <- "young"
des$age[middleage] <- "middleage"
des$age[old] <- "old"
des$age <- as.factor(des$age)
# add one step to Yiming's code to reorder the factor levels
des$age <- recode(des$age, "",levels = c('young', 'middleage', 'old'))

str(des)

# create function to extract genes
prepareData <- function(g) {
  pDat0 <- data.frame()
  for (i in 1:length(g)) {
    pDat <- data.frame(des, gExp = as.vector(t(as.matrix(dat[g[i], ]))), 
                       gene = g[i])
    pDat0 <- rbind(pDat0, pDat)
  }
  pDat0
}


##############################################################################################################################################

# Set up design matrix to get 'reference + treatment effects' fit to test for group only:
# i.e. testing: normal vs. CD, normal vs. UC, and UC vs. DC
desMat <- model.matrix(~group, des)
# Fitting two-way ANOVA to all probesets at once:
fit <- lmFit(dat, desMat)
ebFit <- eBayes(fit)
colnames(ebFit)
# Use topTable to get the hits:
hits <- topTable(ebFit, coef = grep("group", colnames(coef(ebFit))), number = nrow(dat), p.value=1e-2)
str(hits) # we're getting 1283 hits
write.table(hits, "data/groups-hits.tsv")
# explore hits
stripplot(gExp ~ group | gene, prepareData(head(rownames(hits), 6)), jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), grid = TRUE)
# explore non-hits
allprobes <- topTable(ebFit, coef = grep("group", colnames(coef(ebFit))), number = nrow(dat))
stripplot(gExp ~ group | gene, prepareData(tail(rownames(allprobes ), 6)), jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), grid = TRUE)


# Do it using 'cell means' instead of 'reference + treatment effects' to check our results:

# Set up design matrix to get 'cell means' fit to test for group only:
cMat <- model.matrix(~group+0, des)
cFit <- lmFit(dat, cMat)
colnames(cFit)
contMatrix <- makeContrasts(NCvsDC = groupCD - groupNC, NCvsUC = groupUC - groupNC, CDvsUC = groupUC - groupCD, levels = cMat)
fitCont <- contrasts.fit(cFit, contMatrix)
ebFitCont <- eBayes(fitCont)
contHits <- topTable(ebFitCont, number = nrow(dat),  p.value=1e-2)
str(contHits) # we're getting 1283 hits (same as previous method)
# explore hits
stripplot(gExp ~ group | gene, prepareData(head(rownames(contHits), 6)), jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), grid = TRUE)
# explore non-hits
allProbes <- topTable(ebFitCont, number = nrow(dat))
stripplot(gExp ~ group | gene, prepareData(tail(rownames(allProbes), 6)), jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), grid = TRUE)

##############################################################################################################################################


# using 'reference + treatment effects' fit to test for group and age:
aMat <- model.matrix(~group*age, des)
afit <- lmFit(dat, aMat)
aebFit <- eBayes(afit)
colnames(aebFit)
# Test for any effect of group and/or age:
aghits <- topTable(aebFit, coef = c("groupCD", "groupUC", "agemiddleage", "ageold", "groupCD:agemiddleage", "groupUC:agemiddleage", "groupCD:ageold", "groupUC:ageold"), number = nrow(dat), p.value = 0.1)
str(aghits) # from here we are getting 305 observations 
# Test for any effect of group:
ghits <- topTable(aebFit, coef = grep("group", colnames(coef(aebFit))),  number = nrow(dat), p.value = 0.1)
str(ghits) # from here we are getting 921 observations 
write.table(ghits, "data/groupVsage-hits.tsv") # saved the resutls upon Yiming's request
# Test for any effect of age:
ahits <- topTable(aebFit, coef = grep("age", colnames(coef(aebFit))),  number = nrow(dat), p.value = 0.1)
str(ahits) # from here we are getting only 4 observations at the same threshold which means there is no significant effect for age


# using 'reference + treatment effects' fit to test for group and sex:
sMat <- model.matrix(~group*sex, des)
sfit <- lmFit(dat, sMat)
sebFit <- eBayes(sfit)
colnames(sebFit)
# Test for any effect of group and/or sex:
sghits <- topTable(sebFit, coef = c("groupCD", "groupUC", "sexmale", "groupCD:sexmale", "groupUC:sexmale"), number = nrow(dat), p.value = 0.1)
str(sghits) # here we're getting 5096 observations
# Test for any effect of group:
gHits <- topTable(sebFit, coef = grep("group", colnames(coef(sebFit))),  number = nrow(dat), p.value = 0.1)
str(gHits) # here we're getting 5359 observations
# Test for any effect of sex:
shits <- topTable(sebFit, coef = grep("sex", colnames(coef(sebFit))),  number = nrow(dat), p.value = 0.1)
str(shits) # here we are getting 219 observations

