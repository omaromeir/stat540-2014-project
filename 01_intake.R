# Task 1: Data intake
# created by Abrar Wafa
# date: March 16, 2014

library(GEOquery)
library(car)

# load the GEO file
GSE <- getGEO("GSE1710", GSEMatrix=FALSE, destdir="data/")

# creating data file
gsms <- GSE@gsms
mat <- do.call('cbind',lapply(gsms,function(x) {
  tab <- Table(x)
  return(as.numeric(tab$VALUE))
}))
rownames(mat) <- Table(gsms[[1]])$ID_REF
dat <- as.data.frame(mat, stringsAsFactors=FALSE) 

# save data to a file
write.table(dat, "data/GSE1710-data.tsv")

# creating design file
group <- c(rep("NC", each = 11), rep("CD", each=10), rep("UC", each = 10))
sex <- c("male", "female", "male", "male", "female", "female", "female", "female", "female", "male", "female", "male", "female", "male", "female", "female", "male", "male", "male", "male", "female", "male", "male", "male", "male", "male", "female", "male", "male", "female", "male" )
age <- c("44", "47", "27", "55", "46", "57", "84", "28", "76", "72", "41", "20", "38", "19", "32", "24", "66", "25", "44", "50", "18", "29", "38", "27", "31", "29", "38", "33", "51", "29", "43")
des <- data.frame(samples=names(gsms), group = group, sex = sex, age = age)
des$group <- recode(des$group, "'NC'='NC'",levels = c('NC', 'CD', 'UC'))

# save design to a file
saveRDS(des, file = "data/GSE1710-design.rds")

# sanity checks
all(Table(gsms[[2]])$ID_REF == rownames(dat))
identical(as.numeric(Table(gsms$GSM29595)$VALUE), dat$GSM29595)
all(colnames(dat) == names(gsms))
str(des)
summary(des$group)
