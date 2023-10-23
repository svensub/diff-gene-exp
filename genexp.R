################################################################################
#                FINDING A DIFFERENTIAL IN GENE EXPRESSION COUNTS              #
################################################################################

# 1. set working directory
setwd("path")

# 2. install packages and related libraries
mypackages <- c("tidyverse", "gplots", "ggplot2", "BiocManager", "fdrtool", "xtable")
for (p in mypackages){
  if(!require(p, character.only = TRUE)){
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

BiocManager :: install ("limma", force = TRUE)
BiocManager :: install ("affy", force = TRUE)

library(affy)
library(limma)

# 3. load data

DATA <- readRDS("DATA.rds")
COUNTS <- readRDS("COUNTS.rds")

# 4. find and replace missing values with 0

sum(is.na(COUNTS)) # number of missing values is 884, one can assume since very 
                   # small proportion of the data is missing within a large 
                   # matrix, that the NA values may be due to low luminosity

COUNTS[is.na(COUNTS)] <- 0 

# 5. create the expression set object 

COUNTS <- t(COUNTS) # to create the object, must transpose the count matrix 
                    # using t(count)

eset <- new("ExpressionSet", exprs = COUNTS)
     # creating the expression set object
phenoData(eset) <- new("AnnotatedDataFrame", data = DATA) 
                # to store & manipulate phenotypic data along with its metadata 
eset # summary of object
pData(eset) # viewing phenotypic data
exprs(eset)[1:10, ] # viewing expression values

# 6. fit linear model (with only status)

mm <- model.matrix(~eset$Status) # design matrix specifying disease status 
                                 # (+ intercept) as the only coefficients
linmod <- lmFit(eset, design=mm) # fitting linear model for each gene

# 7. perform eBayes

eb <- eBayes(linmod)

# 8. use topTable to adjust for p-values using Bonferroni

tt <- topTable(eb, number=2500, coef=2, adjust = "bonferroni")
tt[1:20, ]

# assuming alpha = 0.05, the bonferroni threshold would c = alpha/m = 0.05/2500
# = 2e-05

sum(tt$P.Value<=2e-05) # no of p values below this threshold is 132

# 9. to repeat above steps (6:8) with covariates to observe differences

mm_cov <- model.matrix(~ Status+Age+Sex, data = eset)
linmod_cov <- lmFit(eset, design=mm_cov)
eb_cov <- eBayes(linmod_cov)
tt_cov <- topTable(eb_cov, number=2500, coef=2, adjust = "bonferroni")
sum(tt_cov$adj.P.Val<=0.05) # no of p values below this threshold is 130

sum(tt_cov$adj.P.Val<=0.001)

# 10. another correction method: FDR

tt_BH <- topTable(eb_cov, number=2500, coef=2, adjust="BH")
sum(tt_BH$adj.P.Val<=0.05)

# 11. creates figures and tables

## TABLE 1: 20 highest differentially expressed genes with p-values

sort(tt_cov$adj.P.Val, decreasing = FALSE)
xtable(tt_cov[1:10,])


## FIGURE 1: histogram of p-values/adj. p-values

hist(tt$P.Value,
     main = "",
     xlab = "p-values",
     cex.lab = 1.3,
     cex.axis = 1.2)
hist(tt_cov$P.Value,
     main = "",
     xlab = "p-values",
     cex.lab = 1.3,
     cex.axis = 1.2)
hist(tt_cov$adj.P.Val,
     main = "",
     xlab = "Adjusted p-values",
     cex.lab = 1.3,
     cex.axis = 1.2)
hist(tt_BH$adj.P.Val,
     main = "",
     xlab = "Adjusted p-values",
     cex.lab = 1.3,
     cex.axis = 1.2)

## FIGURE 2: QQ-plot

plot(-log(seq(0, 1, length.out=2500)), -log(sort(tt$P.Value)),
     main = "",
     xlab = "Expected -log10(p-value)",
     ylab = "Observed -log10(p-value)",
     cex.lab = 1.3,
     cex.axis = 1.2)
abline(0, 1)

plot(-log(seq(0, 1, length.out=2500)), -log(sort(tt_cov$P.Value)),
     main = "",
     xlab = "Expected -log10(p-value)",
     ylab = "Observed -log10(p-value)",
     cex.lab = 1.3,
     cex.axis = 1.2)
abline(0, 1)

## FIGURE 3: heat map

selected <- p.adjust(eb_cov$p.value[, 2])<0.05
esetSel <- eset[selected, ]
coolmap(exprs(esetSel))


