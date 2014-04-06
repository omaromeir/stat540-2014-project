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

# Use topTable to get the hits of all groups:
hits <- topTable(ebFit, coef = grep("group", colnames(coef(ebFit))), number = nrow(dat), p.value=1e-2)
(hits)
str(hits) # we're getting 1283 hits
write.table(hits, "data/groups-hits.tsv") # saving the resutls upon Yiming's request

# Test for any effect of CDgroup:
CDhits <- topTable(ebFit, coef = grep("CD", colnames(coef(ebFit))), number = nrow(dat), p.value=1e-2)
str(CDhits) # we have 1453 hits
write.table(CDhits, "data/CD-hits.tsv")

# Test for any effect of UCgroup:
UChits <- topTable(ebFit, coef = grep("UC", colnames(coef(ebFit))), number = nrow(dat), p.value=1e-2)
str(UChits) # we have 527 hits
write.table(UChits, "data/UC-hits.tsv")

# explore hits of all groups 
stripplot(gExp ~ group | gene, prepareData(c("28I19", "83C07", "33F10", "18B09")), jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), grid = TRUE, asp =1)
# explore non-hits of all groups
allprobes <- topTable(ebFit, coef = grep("group", colnames(coef(ebFit))), number = nrow(dat))
stripplot(gExp ~ group | gene, prepareData(tail(rownames(allprobes ), 4)), jitter.data = TRUE, auto.key = TRUE, type = c("p", "a"), grid = TRUE, asp=1)


##############################################################################################################################################

# Use 'cell means' to test for differentially expressed genes between CD and UC:
# Set up design matrix to get 'cell means' fit:
cMat <- model.matrix(~group+0, des)
cFit <- lmFit(dat, cMat)
contMatrix <- makeContrasts(CDvsUC = groupUC - groupCD,  levels = cMat)
fitCont <- contrasts.fit(cFit, contMatrix)
ebFitCont <- eBayes(fitCont)
contHits <- topTable(ebFitCont, number = nrow(dat),  p.value=1e-2)
str(contHits) # we're getting no hits

##############################################################################################################################################


# using 'reference + treatment effects' fit to test for group and age:
aMat <- model.matrix(~group*age, des)
afit <- lmFit(dat, aMat)
aebFit <- eBayes(afit)
colnames(aebFit)
# Test for any effect of group and/or age:
aghits <- topTable(aebFit, coef = c("groupCD", "groupUC", "agemiddleage", "ageold", "groupCD:agemiddleage", "groupUC:agemiddleage", "groupCD:ageold", "groupUC:ageold"), number = nrow(dat), p.value = 0.1)
str(aghits) # from here we are getting 305 observations, setting the threshold to 0.05, we get only 19 observations 
# Test for any effect of group:
ghits <- topTable(aebFit, coef = grep("group", colnames(coef(aebFit))),  number = nrow(dat), p.value = 0.1)
str(ghits) # from here we are getting 921 observations .. after setting threshold to 0.05, we get 183 observations
write.table(ghits, "data/groupVsage-hits.tsv") # saving the resutls upon Yiming's request
# Test for any effect of age:
ahits <- topTable(aebFit, coef = grep("age", colnames(coef(aebFit))),  number = nrow(dat), p.value = 0.1)
str(ahits) # from here we are getting only 4 observations at the same threshold which means there is no significant effect for age


##############################################################################################################################################


# using 'reference + treatment effects' fit to test for group and sex:
sMat <- model.matrix(~group*sex, des)
sfit <- lmFit(dat, sMat)
sebFit <- eBayes(sfit)
colnames(sebFit)
# Test for any effect of group and/or sex:
sghits <- topTable(sebFit, coef = c("groupCD", "groupUC", "sexmale", "groupCD:sexmale", "groupUC:sexmale"), number = nrow(dat), p.value = 1e-2)
str(sghits) # here we're getting 244 observations
# Test for any effect of group:
gHits <- topTable(sebFit, coef = grep("group", colnames(coef(sebFit))),  number = nrow(dat), p.value = 1e-2)
str(gHits) # here we're getting 417 observations
# Test for any effect of sex:
shits <- topTable(sebFit, coef = grep("sex", colnames(coef(sebFit))),  number = nrow(dat), p.value = 1e-2)
str(shits) # here we are getting 1 observation

