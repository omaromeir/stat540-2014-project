### 01 Data intake
### Created by: Abrar Wafa
### Date: March 16, 2014


This is an R Markdown document for loading, cleaning, and sorting the data.


```r
library(GEOquery)
library(car)
```


Load the GEO file:

```r
GSE <- getGEO("GSE1710", GSEMatrix = FALSE, destdir = "../data/")
```

```
## Parsing....
```


Creating data file:

```r
gsms <- GSE@gsms
mat <- do.call("cbind", lapply(gsms, function(x) {
    tab <- Table(x)
    return(as.numeric(tab$VALUE))
}))
rownames(mat) <- Table(gsms[[1]])$ID_REF
dat <- as.data.frame(mat, stringsAsFactors = FALSE)
```


Save data to a file:

```r
write.table(dat, "../data/GSE1710-data.tsv")
```


Creating design file:

```r
group <- c(rep("NC", each = 11), rep("CD", each = 10), rep("UC", each = 10))
sex <- c("male", "female", "male", "male", "female", "female", "female", "female", 
    "female", "male", "female", "male", "female", "male", "female", "female", 
    "male", "male", "male", "male", "female", "male", "male", "male", "male", 
    "male", "female", "male", "male", "female", "male")
age <- c("44", "47", "27", "55", "46", "57", "84", "28", "76", "72", "41", "20", 
    "38", "19", "32", "24", "66", "25", "44", "50", "18", "29", "38", "27", 
    "31", "29", "38", "33", "51", "29", "43")
des <- data.frame(samples = names(gsms), group = group, sex = sex, age = age)
des$group <- recode(des$group, "'NC'='NC'", levels = c("NC", "CD", "UC"))
```


Save design to a file:

```r
saveRDS(des, file = "../data/GSE1710-design.rds")
```


Sanity checks:

```r
all(Table(gsms[[2]])$ID_REF == rownames(dat))
```

```
## [1] TRUE
```

```r
identical(as.numeric(Table(gsms$GSM29595)$VALUE), dat$GSM29595)
```

```
## [1] TRUE
```

```r
all(colnames(dat) == names(gsms))
```

```
## [1] TRUE
```

```r
str(des)
```

```
## 'data.frame':	31 obs. of  4 variables:
##  $ samples: Factor w/ 31 levels "GSM29595","GSM29596",..: 1 2 3 4 5 6 7 8 9 10 ...
##  $ group  : Factor w/ 3 levels "NC","CD","UC": 1 1 1 1 1 1 1 1 1 1 ...
##  $ sex    : Factor w/ 2 levels "female","male": 2 1 2 2 1 1 1 1 1 2 ...
##  $ age    : Factor w/ 25 levels "18","19","20",..: 15 17 6 20 16 21 25 7 24 23 ...
```

```r
summary(des$group)
```

```
## NC CD UC 
## 11 10 10
```



