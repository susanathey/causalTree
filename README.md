# causalTree Introduction

The _causalTree_ function builds a regression model and returns an _rpart_ object, which is the object derived from _rpart_ package, implemneting many ideas in the CART (Classification and Regression Trees), written by Breiman, Friedman, Olshen and Stone. Like _rpart_, _causalTree_ builds a binary regression tree model in two stages, but focuses on estimating heterogeneous causal effect.


To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools) 
install_github("susanathey/causalTree")
```

Example usage:

```R
library(causalTree)
tree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1, treatment = simulation.1$treatment,
                   split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                   xval = 5, cp = 0, minsize = 20, propensity = 0.5)
                  
opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]

opfit <- prune(tree, opcp)

rpart.plot(opfit)

```

For More details, please check out briefintro.pdf.

#### References
Susan Athey, Guido Imbens. <b>Recursive Partitioning for Heterogeneous Causal Effects.</b> [<a href="http://arxiv.org/abs/1504.01132">link</a>]
